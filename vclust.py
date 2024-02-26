#!/usr/bin/env python3
"""Compute Average Nucleotide Identity (ANI) and cluster virus genome sequences.

https://github.com/refresh-bio/vclust-dev
"""

from pathlib import Path
import argparse
import csv
import multiprocessing
import shutil
import subprocess
import sys
import typing
import uuid

__version__ = '0.5'

DEFAULT_THREAD_COUNT = min(multiprocessing.cpu_count(), 64)

VCLUST_DIR = Path(__file__).resolve().parent

# Default paths to third-party binaries
BIN_DIR = VCLUST_DIR / 'bin'
BIN_KMERDB = BIN_DIR / 'kmer-db'
BIN_LZANI = BIN_DIR / 'lz-ani'
BIN_RAPIDCLUSTER = BIN_DIR / 'rapid-cluster'
BIN_FASTASPLIT = BIN_DIR / 'multi-fasta-split'

# lz-ani output columns
ALIGN_FIELDS = [
    'idx1', 'idx2', 'id1', 'id2',  'tani', 'gani', 'ani', 'cov',
    'num_alns', 'len_ratio', 'len1', 'len2',  'nt_match', 'nt_mismatch', 
]
# vclust align output formats
ALIGN_OUTFMT = {
    'lite': ALIGN_FIELDS[:2] + ALIGN_FIELDS[4:9],
    'standard': ALIGN_FIELDS[:10],
    'complete': ALIGN_FIELDS[:],
}

def get_parser() -> argparse.ArgumentParser:
    """Return an argument parser."""
    parser = argparse.ArgumentParser(
        description=f'%(prog)s v.{__version__}: calculate ANI and cluster '
                     'virus (meta)genome sequences',
        add_help=False,
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'v.{__version__}',
        help="Display the tool's version and exit"
    )
    parser.add_argument(
        '-h', '--help',
        action='help',
        help='Show this help message and exit'
    )

    subparsers = parser.add_subparsers(dest='command')

    # Prefilter parser    
    prefilter_parser = subparsers.add_parser(
        'prefilter',
        help='Prefilter genome pairs for alignment',
        formatter_class=fmt,
        add_help=False,
    )
    prefilter_parser._optionals.title = "Options"
    prefilter_parser.add_argument(
        '-i', '--in',
        metavar='<file>',
        dest='input',
        help='Input FASTA file or directory with FASTA files',
        required=True
    )
    prefilter_parser.add_argument(
        '-o', '--out',
        metavar='<file>',
        dest='output',
        help='Output filename',
        required=True
    )
    prefilter_parser.add_argument(
        '-k', '--k',
        metavar="<int>",
        type=int,
        default=18,
        help="Size of k-mer for kmer-db [%(default)s]"
    )
    prefilter_parser.add_argument(
        '-t', '--threads',
        metavar="<int>",
        dest="num_threads",
        type=int,
        default=DEFAULT_THREAD_COUNT,
        help='Number of threads (all by default) [%(default)s]'
    )
    prefilter_parser.add_argument(
        '--min-kmers',
        metavar="<int>",
        type=int,
        default=2,
        help='Filter genome pairs based on minimum number of shared k-mers '
             '[%(default)s]'
    )
    prefilter_parser.add_argument(
        '--min-ident',
        metavar="<float>",
        type=float,
        default=0.1,
        help='Filter genome pairs based on minimum sequence identity of '
        'the shorter sequence (0-1) [%(default)s]'
    )
    prefilter_parser.add_argument(
        '--batch-size',
        metavar="<int>",
        type=int,
        default=0,
        help='Split a multifasta file into smaller files of n FASTA sequences. '
        'This option reduces memory at the expense of speed. By default, no '
        'batch [%(default)s]'
    )
    prefilter_parser.add_argument(
        '--keep_temp',
        action="store_true",
        help='Keep temporary kmer-db files [%(default)s]'
    )
    prefilter_parser.add_argument(
        '--bin',
        metavar='<file>',
        dest="bin_kmerdb",
        default=f'{BIN_KMERDB.relative_to(VCLUST_DIR)}',
        help='Path to the kmer-db binary [%(default)s]'
    )
    prefilter_parser.add_argument(
        '--bin-fasta',
        metavar='<file>',
        dest="bin_fastasplit",
        default=f'{BIN_FASTASPLIT.relative_to(VCLUST_DIR)}',
        help='Path to the multi-fasta-split binary [%(default)s]'
    )
    prefilter_parser.add_argument(
        '-v', '--verbose',
        action="store_true",
        help="Show kmer-db output"
    )
    prefilter_parser.add_argument(
        '-h', '--help',
        action='help',
        help='Show this help message and exit'
    )

    # Align parser
    align_parser = subparsers.add_parser(
        'align',
        help='Align genome sequences and calculate ANI metrics',
        formatter_class=fmt,
        add_help=False,
    )
    align_parser._optionals.title = "Options"
    align_parser.add_argument(
        '-i', '--in',
        metavar='<file>',
        dest='input',
        help='Input FASTA file or directory with FASTA files',
        required=True
    )
    align_parser.add_argument(
        '-o', '--out',
        metavar='<file>',
        dest='output',
        help='Output filename',
        required=True
    )
    align_parser.add_argument(
        '--outfmt',
        metavar='<str>',
        choices=ALIGN_OUTFMT.keys(),
        dest='outfmt',
        default='standard',
        help='Output format [%(default)s]\n'
        f'choices: {",".join(ALIGN_OUTFMT.keys())}'
    )
    align_parser.add_argument(
        '--out-ani',
        dest='ani',
        metavar='<float>',
        type=float,
        default=0,
        help='Min. ANI to output (0-1) [%(default)s]'
    )
    align_parser.add_argument(
        '--out-tani',
        dest='tani',
        metavar='<float>',
        type=float,
        default=0,
        help='Min. tANI to output (0-1) [%(default)s]'
    )
    align_parser.add_argument(
        '--out-gani',
        dest='gani',
        metavar='<float>',
        type=float,
        default=0,
        help='Min. gANI to output (0-1) [%(default)s]'
    )
    align_parser.add_argument(
        '--out-cov',
        dest='cov',
        metavar='<float>',
        type=float,
        default=0,
        help='Min. coverage to output (0-1) [%(default)s]'
    )
    align_parser.add_argument(
        '--filter',
        metavar='<file>',
        dest="filter_path",
        help='Path to filter file (output of prefilter)'
    )
    align_parser.add_argument(
        '--filter-threshold',
        metavar='<float>',
        dest='filter_threshold',
        type=float,
        default=0,
        help='Align genome pairs above the threshold (0-1) [%(default)s]'
    )
    align_parser.add_argument(
        '-t', '--threads',
        metavar='<int>',
        dest='num_threads',
        type=int,
        default=DEFAULT_THREAD_COUNT,
        help='Number of threads (all by default) [%(default)s]'
    )
    align_parser.add_argument(
        '--bin',
        metavar='<file>',
        dest='bin_lzani',
        default=f'{BIN_LZANI.relative_to(VCLUST_DIR)}',
        help='Path to the lz-ani binary [%(default)s]'
    )
    align_parser.add_argument(
        '--mal',
        metavar='<int>',
        type=int,
        default=11,
        help='Min. anchor length [%(default)s]'
    )       
    align_parser.add_argument(
        '--msl',
        metavar='<int>',
        type=int,
        default=7,
        help='Min. seed length [%(default)s]'
    )
    align_parser.add_argument(
        '--mrd',
        metavar='<int>',
        type=int,
        default=128,
        help='Max. dist. between approx. matches in reference [%(default)s]'
    )
    align_parser.add_argument(
        '--mqd',
        metavar='<int>',
        type=int,
        default=32,
        help='Max. dist. between approx. matches in query [%(default)s]'
    )
    align_parser.add_argument(
        '--reg',
        metavar='<int>',
        type=int,
        default=80,
        help='Min. considered region length [%(default)s]'
    )
    align_parser.add_argument(
        '--aw',
        metavar='<int>',
        type=int,
        default=16,
        help='Approx. window length [%(default)s]'
    )
    align_parser.add_argument(
        '--am',
        metavar='<int>',
        type=int,
        default=6,
        help='Max. no. of mismatches in approx. window [%(default)s]'
    )
    align_parser.add_argument(
        '--ar',
        metavar='<int>',
        type=int,
        default=2,
        help='Min. length of run ending approx. extension [%(default)s]'
    )
    align_parser.add_argument(
        '-v', '--verbose',
        action="store_true",
        help="Show lz-ani standard output"
    )
    align_parser.add_argument(
        '-h', '--help',
        action='help',
        help='Show this help message and exit'
    )

    # Cluster parser
    cluster_parser = subparsers.add_parser(
        'cluster',
        help='Cluster genomes based on ANI thresholds',
        formatter_class=fmt,
        add_help=False,
    )
    cluster_parser._optionals.title = "Options"
    cluster_parser.add_argument(
        '-i', '--in',
        metavar='<file>',
        dest='input',
        help='Input file with ANI metrics (tsv)',
        required=True
    )
    cluster_parser.add_argument(
        '-o', '--out',
        metavar='<file>',
        dest='output',
        help='Output filename',
        required=True
    )
    cluster_parser.add_argument(
        '--ids',
        metavar='<file>',
        dest='ids_path',
        help='Input file with sequence identifiers (tsv)',
        required=True
    )
    cluster_parser.add_argument(
        '-r', '--out-repr',
        action='store_true',
        dest='representatives',
        help='Output a representative genome for each cluster instead of '
        'numerical cluster identifiers. The representative genome is selected '
        'as the one with the longest sequence. [%(default)s]'
    )
    choices = [
        'single', 'complete', 'uclust', 'cd-hit',
        'set-cover', 'connected-component', 'leiden'
    ]
    cluster_parser.add_argument(
        '--algorithm',
        metavar='<str>',
        dest="algorithm",
        choices=choices,
        default='single',
        help='Clustering algorithm [%(default)s]\n'
        '* single: Single-linkage\n'
        '* complete: Complete-linkage\n'
        '* uclust: UCLUST\n'
        '* cd-hit: Greedy incremental\n'
        '* set-cover: Greedy set-cover (MMseqs2)\n'
        '* connected-component: Connected component (BLASTclust)\n'
        '* leiden: the Leiden algorithm'
    )
    choices = ['tani','gani','ani']
    cluster_parser.add_argument(
        '--metric',
        metavar='<str>',
        dest='metric',
        choices=choices,
        default='tani',
        help='Similarity metric for clustering [%(default)s]\n'
        f'choices: {",".join(choices)}'
    )
    cluster_parser.add_argument(
        '--tani',
        metavar='<float>',
        dest='tani',
        type=float,
        default=0,
        help='Min. total ANI (0-1) [%(default)s]\n'
    )
    cluster_parser.add_argument(
        '--gani',
        metavar='<float>',
        dest='gani',
        type=float,
        default=0,
        help='Min. global ANI (0-1) [%(default)s]\n'
    )
    cluster_parser.add_argument(
        '--ani',
        metavar='<float>',
        dest='ani',
        type=float,
        default=0,
        help='Min. ANI (0-1) [%(default)s]\n'
    )
    cluster_parser.add_argument(
        '--cov',
        metavar='<float>',
        dest='cov',
        type=float,
        default=0,
        help='Min. coverage (0-1) [%(default)s]\n'
    )
    cluster_parser.add_argument(
        '--num_alns',
        metavar='<int>',
        dest='num_alns',
        type=int,
        default=0,
        help='Max. number of local alignments between two genomes; 0 means all '
        'genome pairs are allowed. [%(default)s]'
    )
    cluster_parser.add_argument(
        '--bin',
        metavar='<file>',
        dest="bin_rapidcluster",
        default=f'{BIN_RAPIDCLUSTER.relative_to(VCLUST_DIR)}',
        help='Path to the rapid-cluster binary [%(default)s]'
    )
    cluster_parser.add_argument(
        '-h', '--help',
        action='help',
        help='Show this help message and exit'
    )

    # Display help message if the script is executed without any arguments.
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    return parser


def get_uuid() -> str:
    """Return unique identifier"""
    return f'vclust-{str(uuid.uuid4().hex)[:16]}'


def verify_binary(bin_path: typing.Union[str, Path]) -> Path:
    """Verify if a binary file exists and is executable.

    Args:
        bin_path (str or Path): The path to the executable binary file.

    Returns:
        Path to the binary file.

    Raises:
        SystemExit: If the binary file is not found.

    """    
    bin_path = Path(bin_path) if isinstance(bin_path, str) else bin_path
    try:
        subprocess.run(
            [f'{bin_path}'],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        exit(f'error: executable not found: {bin_path.resolve()}')
    return bin_path


def validate_args(args, parser):
    args.input = Path(args.input)
    if not args.input.exists():
        parser.error(f'input does not exist: {args.input}')
    args.output = Path(args.output)
    return args


def validate_args_fasta_input(args, parser):
    args.is_multifasta = True
    args.fasta_paths = [args.input]
    if args.input.is_dir():
        args.is_multifasta = False
        args.fasta_paths = sorted(
            f for f in args.input.iterdir() if f.is_file()
        )
    if not args.is_multifasta and (n := len(args.fasta_paths)) < 2:
        parser.error(f'too few fasta files ({n}) in {args.input}')
    return args


def validate_args_prefilter(args, parser):
    if args.min_ident < 0 or args.min_ident > 1:
        parser.error('--min-ident must be between 0 and 1, inclusive.')
    if args.k < 15 or args.k > 30:
        parser.error(f'--k must be in range 15-30')       
    if args.batch_size and args.input.is_dir():
        parser.error('--batch-size only handles a multi-fasta file'
            ', not a directory.')
    return args


def validate_args_align(args, parser):
    if args.filter_path:
        args.filter_path = Path(args.filter_path)
        if not args.filter_path.exists():
            parser.error(f'file does not exist: {args.filter_path}')

    args_dict = vars(args)
    options = ['tani', 'gani', 'ani', 'cov']
    # Verify if option values are within the appropriate range.
    for name in options:
        value = args_dict[name]
        if value < 0 or value > 1:
            parser.error(f'--out-{name} must be between 0 and 1, inclusive.')
    return args


def validate_args_cluster(args, parser):
    # Check if the ids file exists.
    args.ids_path = Path(args.ids_path)
    if not args.ids_path.exists():
        parser.error(f'file does not exist: {args.ids_path}')

    # Validate the ani.tsv file.
    with open(args.input) as fh:
        header = fh.readline().split()
        if 'idx1' not in header and 'idx2' not in header:
            parser.error(f'missing columns `idx1` and `idx2` in {args.input}')

        args_dict = vars(args)
        args.metric_threshold = args_dict.get(args.metric, 0)
        # Check if the threshold for the metric is above 0.
        if not args.metric_threshold:
            parser.error(f'{args.metric} threshold must be above 0. '
                f'Specify the option: --{args.metric}')
        
        options = ['tani', 'gani', 'ani', 'cov', 'num_alns']
        # Verify if option values are within the appropriate range.
        for name in options[:-1]:
            value = args_dict[name]
            if value < 0 or value > 1:
                parser.error(f'--{name} must be between 0 and 1, inclusive.')

        # Verify presence of columns for options with thresholds greater than 0.
        for name in options:
            value = args_dict[name]
            if value != 0 and name not in header:
                parser.error(f'missing column `{name}` in {args.input}')
    return args


def validate_process(func):
    """Decorator to validate the result of a subprocess execution."""
    def wrapper(*args, **kwargs):
        process = func(*args, **kwargs)
        if process.returncode:
            print(f'error while running: {" ".join(process.args)}')
            print(f'error message: {process.stderr}')
            exit(1)
        return process
    return wrapper


@validate_process
def run_fastasplit(input_fasta, out_dir, n, verbose, bin_path = BIN_FASTASPLIT):
    # Run kmer-db build.
    cmd = [
        f"{bin_path}", 
        "-n", f"{n}",
        "--out-prefix",
        f"{out_dir}/part",
        f'{input_fasta}',
    ]
    process = subprocess.run(
        cmd,  
        stdout=None if verbose else subprocess.DEVNULL, 
        stderr=subprocess.PIPE,
        text=True,
    )
    return process    


@validate_process
def run_kmerdb_build(
        input_paths: Path,
        txt_path: typing.Union[Path, str],
        db_path: typing.Union[Path, str],
        is_multisample_fasta: bool,
        kmer_size: int,
        num_threads: int,
        verbose: bool = False,
        bin_path: typing.Union[Path, str] = BIN_KMERDB
    ) -> subprocess.CompletedProcess:
    """Runs kmer-db build to create a database from FASTA file(s).

    Args:
        input_fasta (Path):
            Path to the input FASTA file or directory with input FASTA files.
        outfile_txt (Path or str):
            Path to the output text file that will list the input FASTA files.
        outfile_db (Path or str):
            Path to the output kmer-db database file.
        kmer_size (int):
            k-mer size.
        num_threads (int):
            Number of threads to use in kmer-db.
        verbose (bool):
            Whether to display verbose output.
        BIN_KMERDB (Path or str):
            Path to the kmer-db executable.

    Returns:
        subprocess.CompletedProcess: Completed process information.

    """
    # Create a text file listing input FASTA files.
    with open(txt_path, 'w') as oh:
        for f in input_paths:
            oh.write(f'{f}\n')
        
    # Run kmer-db build.
    cmd = [
        f"{bin_path}", 
        "build",
        "-k", f"{kmer_size}",    
        "-t", f"{num_threads}",
        f'{txt_path}',
        f'{db_path}',
    ]
    if is_multisample_fasta:
        cmd.insert(2, '-multisample-fasta')
    process = subprocess.run(
        cmd,  
        stdout=None if verbose else subprocess.DEVNULL, 
        stderr=subprocess.PIPE,
        text=True,
    )
    return process


@validate_process
def run_kmerdb_all2all(
        db_paths,
        db_list_path,
        outfile_all2all: typing.Union[Path, str],
        kmer_count: int,
        num_threads: int,
        verbose: bool = False,
        bin_path: typing.Union[Path, str] = BIN_KMERDB
    ) -> subprocess.CompletedProcess:
    """Runs kmer-db all2all to find shared k-mers between genomic sequences.

    Args:
        infile_db (Path or str):
            Path to the input kmer-db database file.
        outfile_all2all (Path or str):
            Path to the output all2all file.
        kmer_count (int):
            Minimum number of shared k-mers to report in all2all output.
        num_threads (int):
            Number of threads to use in kmer-db.
        verbose (bool):
            Whether to display verbose output.
        bin_path (Path or str):
            Path to the kmer-db executable.

    Returns:
        subprocess.CompletedProcess: Completed process information.

    """
    with open(db_list_path, 'w') as oh:
        for db_path in db_paths:
            oh.write(f'{db_path}\n')

    cmd = [
        f"{bin_path}", 
        'all2all-parts' if len(db_paths) > 1 else 'all2all',
        '-sparse',
        '-above', f'{kmer_count}', 
        "-t", f"{num_threads}",
        f'{db_list_path}' if len(db_paths) > 1 else f'{db_paths[0]}',
        f'{outfile_all2all}',
    ]
    process = subprocess.run(
        cmd,  
        stdout=None if verbose else subprocess.DEVNULL, 
        stderr=subprocess.PIPE,
        text=True,
    )
    return process


@validate_process
def run_kmerdb_distance(
        infile_all2all: typing.Union[Path, str],
        min_value: float,
        num_threads: int,
        verbose: bool = False,
        bin_path: typing.Union[Path, str] = BIN_KMERDB
    ) -> subprocess.CompletedProcess:
    """Runs kmer-db all2all to find shared k-mers between genomic sequences.

    Args:
        infile_db (Path or str):
            Path to the input kmer-db database file.
        outfile_all2all (Path or str):
            Path to the output all2all file.
        kmer_count (int):
            Minimum number of shared k-mers to report in all2all output.
        num_threads (int):
            Number of threads to use in kmer-db.
        verbose (bool):
            Whether to display verbose output.
        bin_path (Path or str):
            Path to the kmer-db executable.

    Returns:
        subprocess.CompletedProcess: Completed process information.

    """
    cmd = [
        f"{bin_path}", 
        "distance",
        "ani-shorter",
        "-sparse",
        '-above', f'{min_value}',
        "-t", f"{num_threads}",
        f'{infile_all2all}',
    ]
    process = subprocess.run(
        cmd,  
        stdout=None if verbose else subprocess.DEVNULL, 
        stderr=subprocess.PIPE,
        text=True,
    )
    return process


@validate_process
def run_lzani(
        input_paths,
        txt_path,
        output_path,
        out_format,
        out_tani,
        out_gani,
        out_ani,
        out_cov,
        filter_file,
        filter_threshold,
        mal,
        msl,
        mrd,
        mqd,
        reg,
        aw,
        am,
        ar,
        num_threads: int,
        verbose: bool = False,
        bin_path: typing.Union[Path, str] = BIN_LZANI
    ) -> subprocess.CompletedProcess:
    """Runs lz-ani to align genomic sequences.

    Args:
        input_fasta (Path):
            Path to the input FASTA file or directory with input FASTA files.
        infile_txt (Path or str):
            Path to the input text file listing FASTA files.
    Returns:
        subprocess.CompletedProcess: Completed process information.

    """
    # Create a text file listing input FASTA files.
    with open(txt_path, 'w') as oh:
        for f in input_paths:
            oh.write(f'{f}\n')

    cmd = [
        f'{bin_path}',
        'all2all',
        '--in-txt',
        f'{txt_path}',
        '-o',
        f'{output_path}',
        '-t',
        f'{num_threads}',
        '--verbose', f'{int(verbose) + 1}',
        '--mal', f'{mal}',
        '--msl', f'{msl}',
        '--mrd', f'{mrd}',
        '--mqd', f'{mqd}',
        '--reg', f'{reg}',
        '--aw', f'{aw}',
        '--am', f'{am}',
        '--ar', f'{ar}',
        '--multisample-fasta',
        'true' if len(input_paths) == 1 else 'false',
        '--out-type', 'tsv',
        '--out-format',
        ','.join(out_format),
    ]
    if filter_file:
        cmd.extend(['--flt-kmerdb', f'{filter_file}', f'{filter_threshold}'])

    cols = [
        ('tani', out_tani), ('gani', out_gani),
        ('ani', out_ani), ('cov', out_cov)
    ]
    for name, value in cols:
        if value > 0:
            cmd.extend(['--out-filter', f'{name}', f'{value}'])

    process = subprocess.run(
        cmd,  
        stdout=None if verbose else subprocess.DEVNULL, 
        stderr=None if verbose else subprocess.PIPE,
        text=True,
    )
    return process


@validate_process
def run_rapidcluster(
        input_path,
        ids_path,
        output_path,
        algorithm,
        metric,
        metric_threshold,
        tani,
        gani,
        ani,
        cov,
        num_alns,
        is_representatives,
        verbose: bool = False,
        bin_path=BIN_RAPIDCLUSTER,
    ):
    cmd = [
        f'{bin_path}',
        '--objects-file',
        f'{ids_path}',
        '--algo',
        f'{algorithm}',
        f'--id-cols',
        'idx1', 'idx2',
        '--distance-col',
        f'{metric}',
        '--similarity',
        '--numeric-ids',
    ]
    cols = [('tani', tani), ('gani', gani), ('ani', ani), ('cov', cov)]
    for name, value in cols:
        if value > 0:
            cmd.extend(['--min', f'{name}', f'{value}'])
    if num_alns > 0:
        cmd.extend(['--max', 'num_alns', f'{num_alns}'])
    if is_representatives:
        cmd.append('--out-representatives')    
    cmd.extend([f'{input_path}', f'{output_path}'])
    process = subprocess.run(
        cmd,  
        stdout=None if verbose else subprocess.DEVNULL, 
        stderr=None if verbose else subprocess.PIPE,
        text=True,
    )
    return process



class CustomHelpFormatter(argparse.HelpFormatter):
    """Custom help message formatter for argparse."""

    def _format_action_invocation(self, action):
        # Allows options with arguments to be formatted as
        # "arg1, arg2 metavar" instead of the "arg1 metavar, arg2 metavar".
        # https://stackoverflow.com/a/31124505
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ', '.join(action.option_strings) + ' ' + args_string

    def _split_lines(self, text, width):
        # Allows inserting new lines on argparse help text.
        # https://stackoverflow.com/a/56865996
        r = []
        for t in text.splitlines():
            r.extend(argparse.HelpFormatter._split_lines(self, t, width))
        return r

fmt = lambda prog: CustomHelpFormatter(prog)


def main():
    parser = get_parser()
    args = parser.parse_args()
    args = validate_args(args, parser)

    # Prefilter
    if args.command == 'prefilter':
        
        args.bin_kmerdb = verify_binary(args.bin_kmerdb)
        args = validate_args_prefilter(args, parser)
        args = validate_args_fasta_input(args, parser)

        out_dir = Path(args.output).parent / get_uuid()
        out_dir.mkdir(parents=True, exist_ok=True)

        batches = []
        # Input is a directory of fasta files.
        if not args.is_multifasta:
            batches.append(args.fasta_paths)
        else:
            # Split multi-fasta file.
            if args.batch_size:
                args.bin_fastasplit = verify_binary(args.bin_fastasplit)
                p = run_fastasplit(
                    input_fasta=args.input, 
                    out_dir=out_dir,
                    n=args.batch_size,
                    verbose=args.verbose,
                    bin_path=args.bin_fastasplit,
                )
                for f in out_dir.glob('part_*'):
                    batches.append([f])
                batches.sort()
            # Do not split multi-fasta file.
            else:
                batches.append([args.input])

        num_batches = len(batches)
        db_paths = []
        for i, batch in enumerate(batches):
            batch_id = f'part_{i:05d}' if num_batches > 1 else 'whole'
            txt_path = out_dir / f'{batch_id}.txt'
            db_path = out_dir / f'{batch_id}.kdb'

            # Run kmer-db build.
            p = run_kmerdb_build(
                input_paths=batch, 
                txt_path=txt_path,
                db_path=db_path,
                is_multisample_fasta=args.is_multifasta,
                kmer_size=args.k,
                num_threads=args.num_threads,
                verbose=args.verbose,
                bin_path=args.bin_kmerdb,
            )
            db_paths.append(db_path)

            # Regardless of verbosity, always delete the partial FASTA file 
            # after building the corresponding partial k-mer database.
            if num_batches > 1:
                batch[0].unlink()

        # Run kmer-db all2all.
        db_list_path = out_dir / 'db_list.txt'
        all2all_path = out_dir / 'all2all.txt'

        p = run_kmerdb_all2all(
            db_paths=db_paths,
            db_list_path=db_list_path,
            outfile_all2all=all2all_path,
            kmer_count=args.min_kmers,
            num_threads=args.num_threads,
            verbose=args.verbose,
            bin_path=args.bin_kmerdb,
        )

        # Run kmer-db distance.
        p = run_kmerdb_distance(
            infile_all2all=all2all_path,
            min_value=args.min_ident,
            num_threads=args.num_threads,
            verbose=args.verbose,
            bin_path=args.bin_kmerdb,
        )

        out_ani = out_dir / 'all2all.txt.ani-shorter'
        out_ani.rename(args.output)

        if not args.keep_temp:
            if out_dir.exists():
                shutil.rmtree(out_dir)

    # Align
    elif args.command == 'align':
        args.bin_lzani = verify_binary(args.bin_lzani)
        args = validate_args_align(args, parser)
        args = validate_args_fasta_input(args, parser)

        out_dir = Path(args.output).parent / get_uuid()
        out_dir.mkdir(parents=True, exist_ok=True)
        txt_path = out_dir / 'ids.txt'

        # Run lz-ani.
        p = run_lzani(
            input_paths=args.fasta_paths,
            txt_path=txt_path,
            output_path=args.output,
            out_format=ALIGN_OUTFMT[args.outfmt],
            out_tani=args.tani,
            out_gani=args.gani,
            out_ani=args.ani,
            out_cov=args.cov,
            filter_file=args.filter_path,
            filter_threshold=args.filter_threshold,
            mal=args.mal,
            msl=args.msl,
            mrd=args.mrd,
            mqd=args.mqd,
            reg=args.reg,
            aw=args.aw,
            am=args.am,
            ar=args.ar,
            num_threads=args.num_threads,
            verbose=args.verbose,
            bin_path=args.bin_lzani,
        )

        if out_dir.exists():
            shutil.rmtree(out_dir)

    # Cluster
    elif args.command == 'cluster':

        args.bin_rapidcluster = verify_binary(args.bin_rapidcluster)
        args = validate_args_cluster(args, parser)

        # Run rapid-cluster
        p = run_rapidcluster(
            input_path=args.input,
            ids_path=args.ids_path,
            output_path=args.output,
            algorithm=args.algorithm,
            metric=args.metric,
            metric_threshold=args.metric_threshold,
            tani=args.tani,
            gani=args.gani,
            ani=args.ani,
            cov=args.cov,
            num_alns=args.num_alns,
            is_representatives=args.representatives,
            verbose=1,
            bin_path=args.bin_rapidcluster,
        )


if __name__ == '__main__':
    main()