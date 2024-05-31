#!/usr/bin/env python3
"""Compute Average Nucleotide Identity (ANI) and cluster virus genome sequences.

https://github.com/refresh-bio/vclust-dev
"""

import argparse
import logging
import multiprocessing
import pathlib
import shutil
import subprocess
import sys
import typing
import uuid

__version__ = '1.0'

DEFAULT_THREAD_COUNT = min(multiprocessing.cpu_count(), 64)

VCLUST_DIR = pathlib.Path(__file__).resolve().parent

# Default paths to third-party binaries
BIN_DIR = VCLUST_DIR / 'bin'
BIN_KMERDB = BIN_DIR / 'kmer-db'
BIN_LZANI = BIN_DIR / 'lz-ani'
BIN_CLUSTY = BIN_DIR / 'clusty'
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

    fmt = lambda prog: CustomHelpFormatter(prog, max_help_position=32)

    def input_path_type(value):
        path = pathlib.Path(value)
        if not path.exists():
            msg = f'input does not exist: {value}'
            raise argparse.ArgumentTypeError(msg)
        return path

    def ranged_float_type(value):
        f = float(value)
        if f < 0 or f > 1:
            raise argparse.ArgumentTypeError('must be between 0 and 1')
        return f

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

    prefilter_optional = prefilter_parser._action_groups.pop()
    prefilter_required = prefilter_parser.add_argument_group('required arguments')
    prefilter_parser._action_groups.append(prefilter_optional)

    prefilter_required.add_argument(
        '-i', '--in',
        metavar='<file>',
        type=input_path_type,
        dest='input_path',
        help='Input FASTA file or directory with FASTA files',
        required=True
    )
    prefilter_required.add_argument(
        '-o', '--out',
        metavar='<file>',
        type=pathlib.Path,
        dest='output_path',
        help='Output filename',
        required=True
    )
    prefilter_parser.add_argument(
        '-k', '--k',
        metavar="<int>",
        type=int,
        default=25,
        choices=range(15, 31),
        help="Size of k-mer for Kmer-db [%(default)s]"
    )
    prefilter_parser.add_argument(
        '--min-kmers',
        metavar="<int>",
        type=int,
        default=10,
        help='Filter genome pairs based on minimum number of shared k-mers '
             '[%(default)s]'
    )
    prefilter_parser.add_argument(
        '--min-ident',
        metavar="<float>",
        type=ranged_float_type,
        default=0.7,
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
        help='Keep temporary Kmer-db files [%(default)s]'
    )
    prefilter_parser.add_argument(
        '--bin',
        metavar='<file>',
        type=pathlib.Path,
        dest="bin_kmerdb",
        default=f'{BIN_KMERDB}',
        help='Path to the Kmer-db binary [%(default)s]'
    )
    prefilter_parser.add_argument(
        '--bin-fasta',
        metavar='<file>',
        type=pathlib.Path,
        dest="bin_fastasplit",
        default=f'{BIN_FASTASPLIT}',
        help='Path to the multi-fasta-split binary [%(default)s]'
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
        '-v', '--verbose',
        action="store_true",
        help="Show Kmer-db progress"
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
    align_optional = align_parser._action_groups.pop()
    align_required = align_parser.add_argument_group('required arguments')
    align_parser._action_groups.append(align_optional)
    align_required.add_argument(
        '-i', '--in',
        metavar='<file>',
        type=input_path_type,
        dest='input_path',
        help='Input FASTA file or directory with FASTA files',
        required=True
    )
    align_required.add_argument(
        '-o', '--out',
        metavar='<file>',
        type=pathlib.Path,
        dest='output_path',
        help='Output filename',
        required=True
    )
    align_parser.add_argument(
        '--filter',
        metavar='<file>',
        type=input_path_type,
        dest="filter_path",
        help='Path to filter file (output of prefilter)'
    )
    align_parser.add_argument(
        '--filter-threshold',
        metavar='<float>',
        dest='filter_threshold',
        type=ranged_float_type,
        default=0,
        help='Align genome pairs above the threshold (0-1) [%(default)s]'
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
        type=ranged_float_type,
        default=0,
        help='Min. ANI to output (0-1) [%(default)s]'
    )
    align_parser.add_argument(
        '--out-tani',
        dest='tani',
        metavar='<float>',
        type=ranged_float_type,
        default=0,
        help='Min. tANI to output (0-1) [%(default)s]'
    )
    align_parser.add_argument(
        '--out-gani',
        dest='gani',
        metavar='<float>',
        type=ranged_float_type,
        default=0,
        help='Min. gANI to output (0-1) [%(default)s]'
    )
    align_parser.add_argument(
        '--out-cov',
        dest='cov',
        metavar='<float>',
        type=ranged_float_type,
        default=0,
        help='Min. coverage (aligned fraction) to output (0-1) [%(default)s]'
    )
    align_parser.add_argument(
        '--bin',
        metavar='<file>',
        type=pathlib.Path,
        dest='bin_lzani',
        default=f'{BIN_LZANI}',
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
        default=40,
        help='Max. dist. between approx. matches in reference [%(default)s]'
    )
    align_parser.add_argument(
        '--mqd',
        metavar='<int>',
        type=int,
        default=40,
        help='Max. dist. between approx. matches in query [%(default)s]'
    )
    align_parser.add_argument(
        '--reg',
        metavar='<int>',
        type=int,
        default=35,
        help='Min. considered region length [%(default)s]'
    )
    align_parser.add_argument(
        '--aw',
        metavar='<int>',
        type=int,
        default=15,
        help='Approx. window length [%(default)s]'
    )
    align_parser.add_argument(
        '--am',
        metavar='<int>',
        type=int,
        default=7,
        help='Max. no. of mismatches in approx. window [%(default)s]'
    )
    align_parser.add_argument(
        '--ar',
        metavar='<int>',
        type=int,
        default=3,
        help='Min. length of run ending approx. extension [%(default)s]'
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
        '-v', '--verbose',
        action="store_true",
        help="Show LZ-ANI progress"
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
    cluster_optional = cluster_parser._action_groups.pop()
    cluster_required = cluster_parser.add_argument_group('required arguments')
    cluster_parser._action_groups.append(cluster_optional)

    cluster_required.add_argument(
        '-i', '--in',
        metavar='<file>',
        type=input_path_type,
        dest='input_path',
        help='Input file with ANI metrics (tsv)',
        required=True
    )
    cluster_required.add_argument(
        '-o', '--out',
        metavar='<file>',
        type=pathlib.Path,
        dest='output_path',
        help='Output filename',
        required=True
    )
    cluster_required.add_argument(
        '--ids',
        metavar='<file>',
        type=input_path_type,
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
    choices = ['single', 'complete', 'uclust', 'cd-hit', 'set-cover', 'leiden']
    cluster_parser.add_argument(
        '--algorithm',
        metavar='<str>',
        dest="algorithm",
        choices=choices,
        default='single',
        help='Clustering algorithm [%(default)s]\n'
        '* single: Single-linkage (connected component)\n'
        '* complete: Complete-linkage\n'
        '* uclust: UCLUST\n'
        '* cd-hit: Greedy incremental\n'
        '* set-cover: Greedy set-cover (MMseqs2)\n'
        '* leiden: Leiden algorithm'
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
        type=ranged_float_type,
        default=0,
        help='Min. total ANI (0-1) [%(default)s]\n'
    )
    cluster_parser.add_argument(
        '--gani',
        metavar='<float>',
        dest='gani',
        type=ranged_float_type,
        default=0,
        help='Min. global ANI (0-1) [%(default)s]\n'
    )
    cluster_parser.add_argument(
        '--ani',
        metavar='<float>',
        dest='ani',
        type=ranged_float_type,
        default=0,
        help='Min. ANI (0-1) [%(default)s]\n'
    )
    cluster_parser.add_argument(
        '--cov',
        metavar='<float>',
        dest='cov',
        type=ranged_float_type,
        default=0,
        help='Min. coverage/aligned fraction (0-1) [%(default)s]\n'
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
        '--leiden-resolution',
        metavar='<float>',
        type=ranged_float_type,
        default=0.7,
        help='Resolution parameter for the Leiden algorithm [%(default)s]\n'
    )
    cluster_parser.add_argument(
        '--bin',
        metavar='<file>',
        type=pathlib.Path,
        dest="BIN_CLUSTY",
        default=f'{BIN_CLUSTY}',
        help='Path to the Clusty binary [%(default)s]'
    )
    cluster_parser.add_argument(
        '-v', '--verbose',
        action="store_true",
        help="Show Clusty progress"
    )
    cluster_parser.add_argument(
        '-h', '--help',
        action='help',
        help='Show this help message and exit'
    )

    # Show help message if the script is run without any arguments.
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    # Show subparser help message if the script is run without any arguments.
    subparsers = [
        ('prefilter', prefilter_parser),
        ('align', align_parser),
        ('cluster', cluster_parser),
    ]
    for name, subparser in subparsers:
        if sys.argv[-1] == name:
            subparser.print_help()
            parser.exit()

    return parser


def create_logger(name: str, log_level: int = logging.INFO) -> logging.Logger:
    """Returns a logger to log events.

    Args:
        name:
            Name of the logger.
        log_level:
            The numeric level of the logging event (one of DEBUG, INFO etc.).

    """
    logger = logging.getLogger(name)
    logger.setLevel(log_level)

    # Set log format to handlers
    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')

    # Create stream logger handler
    sh = logging.StreamHandler()
    sh.setLevel(log_level)
    sh.setFormatter(formatter)
    logger.addHandler(sh)

    return logger


def get_uuid() -> str:
    """Returns a unique string identifier."""
    return f'vclust-{str(uuid.uuid4().hex)[:16]}'


def validate_binary(bin_path: pathlib.Path) -> pathlib.Path:
    """Verifies if a binary file exists and is executable.

    Args:
        bin_path (Path): Path to the executable binary file.

    Returns:
        Path to the binary file.

    Raises:
        SystemExit: If the binary file is not found or not executable.

    """
    if not bin_path.exists():
        exit(f'error: executable not found: {bin_path.resolve()}')
    if not shutil.which(bin_path):
        exit(f'error: file not executable: {bin_path.resolve()}')
    return bin_path


def validate_args_fasta_input(args, parser) -> argparse.Namespace:
    """Validates the arguments for FASTA input."""
    args.is_multifasta = True
    args.fasta_paths = [args.input_path]

    if args.input_path.is_dir():
        args.is_multifasta = False
        args.fasta_paths = sorted(
            f for f in args.input_path.iterdir() if f.is_file()
        )

    if not args.is_multifasta and len(args.fasta_paths) < 2:
        parser.error(f'Too few fasta files found in {args.input_path}. '
                     f'Expected at least 2, but found {len(args.fasta_paths)}.')

    return args


def validate_args_prefilter(args, parser) -> argparse.Namespace:
    """Validates the arguments for the prefilter command."""
    if args.batch_size and args.input_path.is_dir():
        parser.error('--batch-size only handles a multi-fasta file'
            ', not a directory.')
    return args


def validate_args_cluster(args, parser) -> argparse.Namespace:
    """Validates the arguments for the cluster command."""
    # Check the metric and its threshold.
    args_dict = vars(args)
    args.metric_threshold = args_dict.get(args.metric, 0)
    if not args.metric_threshold:
        parser.error(f'{args.metric} threshold must be above 0. '
            f'Specify the option: --{args.metric}')

    # Check if the input TSV file has the required columns.
    with open(args.input_path) as fh:
        header = fh.readline().split()
        if 'idx1' not in header and 'idx2' not in header:
            parser.error(
                f'missing columns `idx1` and `idx2` in {args.input_path}')
        options = ['tani', 'gani', 'ani', 'cov', 'num_alns']
        for name in options:
            value = args_dict[name]
            if value != 0 and name not in header:
                parser.error(f'missing column `{name}` in {args.input_path}')
    return args


def run(
        cmd: typing.List[str],
        verbose: bool,
        logger: logging.Logger
    ) -> subprocess.CompletedProcess:
    """Runs a given command as a subprocess and handle logging.

    Returns:
        subprocess.CompletedProcess: Completed process information.

    """
    logger.info(f'Running: {" ".join(cmd)}')
    process = subprocess.run(
        cmd,  
        stdout=subprocess.DEVNULL, 
        stderr=None if verbose else subprocess.PIPE,
        text=True,
    )
    if process.returncode:
        logger.error(f'While running: {" ".join(process.args)}')
        logger.error(f'Error message: {process.stderr}')
        exit(1)
    logger.info(f'Done')
    return process


def cmd_fastasplit(
        input_fasta: pathlib.Path,
        out_dir: pathlib.Path,
        n: int,
        verbose: bool,
        bin_path = BIN_FASTASPLIT
    ) -> typing.List[str]:
    """Constructs the command line for multi-fasta-split. 
    
    Args:
        input_fasta (Path):
            Path to the input FASTA file.
        out_dir (Path):
            Path to the output directory.
        n (int):
            Number of sequences per output FASTA file.
        verbose (bool):
            Whether to display verbose output.
        bin_path (Path):
            Path to the multi-fasta-split executable.
        
    Returns:
        list: The constructed command as a list of strings.
        
    """
    cmd = [
        f'{bin_path}', 
        '-n', f'{n}',
        f'--verbosity',
        f'{int(verbose)}',
        '--out-prefix',
        f'{out_dir}/part',
        f'{input_fasta}',
    ]
    return cmd 


def cmd_kmerdb_build(
        input_paths: pathlib.Path,
        txt_path: pathlib.Path,
        db_path: pathlib.Path,
        is_multisample_fasta: bool,
        kmer_size: int,
        num_threads: int,
        verbose: bool = False,
        bin_path: pathlib.Path = BIN_KMERDB
    ) -> typing.List[str]:
    """Constructs the command line for Kmer-db build.

    Args:
        input_fasta (Path):
            Path to the input FASTA file or directory with input FASTA files.
        outfile_txt (Path):
            Path to the output text file that will list the input FASTA files.
        outfile_db (Path):
            Path to the output kmer-db database file.
        kmer_size (int):
            k-mer size.
        num_threads (int):
            Number of threads to use in kmer-db.
        verbose (bool):
            Whether to display verbose output.
        bin_path (Path):
            Path to the kmer-db executable.

    Returns:
        list: The constructed command as a list of strings.

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
    return cmd


def cmd_kmerdb_all2all(
        db_paths: typing.List[pathlib.Path],
        db_list_path: pathlib.Path,
        outfile_all2all: pathlib.Path,
        kmer_count: int,
        num_threads: int,
        verbose: bool,
        bin_path: pathlib.Path = BIN_KMERDB
    ) -> typing.List[str]:
    """Constructs the command line for Kmer-db all2all.

    Args:
        db_paths (list[Path]):
            List of paths to the input kmer-db database files.
        db_list_path (Path):
            Path to the output text file listing the kmer-db database files.
        outfile_all2all (Path):
            Path to the output all2all file.
        kmer_count (int):
            Minimum number of shared k-mers to report in all2all output.
        num_threads (int):
            Number of threads to use in kmer-db.
        verbose (bool):
            Whether to display verbose output.
        bin_path (Path):
            Path to the kmer-db executable.
        
    Returns:
        list: The constructed command as a list of strings.

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
    return cmd


def cmd_kmerdb_distance(
        infile_all2all: pathlib.Path,
        min_value: float,
        num_threads: int,
        verbose: bool = False,
        bin_path: pathlib.Path = BIN_KMERDB
    ) -> typing.List[str]:
    """Constructs the command line for Kmer-db distance.

    Args:
        infile_all2all (Path):
            Path to the input all2all file.
        min_value (float):
            Minimum value to filter all2all output.
        num_threads (int):
            Number of threads to use in kmer-db.
        verbose (bool):
            Whether to display verbose output.
        bin_path (Path):
            Path to the kmer-db executable.

    Returns:
        list: The constructed command as a list of strings.

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
    return cmd


def cmd_lzani(
        input_paths: typing.List[pathlib.Path],
        txt_path: pathlib.Path,
        output_path: pathlib.Path,
        out_format: typing.List[str],
        out_tani: float,
        out_gani: float,
        out_ani: float,
        out_cov: float,
        filter_file: pathlib.Path,
        filter_threshold: float,
        mal: int,
        msl: int,
        mrd: int,
        mqd: int,
        reg: int,
        aw: int,
        am: int,
        ar: int,
        num_threads: int,
        verbose: bool,
        bin_path: pathlib.Path = BIN_LZANI
    ) -> typing.List[str]:
    """Constructs the command line for LZ-ANI.

    Args:
        input_paths (List[Path]):
            List of paths to the input FASTA files.
        txt_path (Path):
            Path to the output text file listing the input FASTA files.
        output_path (Path):
            Path to the output ANI file.
        out_format (List[str]):
            List of LZ-ANI column names.
        out_tani (float):
            Minimum tANI to output.
        out_gani (float):
            Minimum gANI to output.
        out_ani (float):
            Minimum ANI to output.
        out_cov (float):
            Minimum coverage (aligned fraction) to output.
        filter_file (Path):
            Path to the filter file (prefilter's output).
        filter_threshold (float):
            Filter threshold.
        mal (int):
            Minimum anchor length.
        msl (int):
            Minimum seed length.
        mrd (int):
            Maximum distance between approximate matches in reference.
        mqd (int):
            Maximum distance between approximate matches in query.
        reg (int):
            Minimum considered region length.
        aw (int):
            Approximate window length.
        am (int):
            Maximum number of mismatches in approximate window.
        ar (int):
            Minimum length of run ending approximate extension.
        num_threads (int):
            Number of threads to use in lz-ani.
        verbose (bool):
            Whether to display verbose output.
        bin_path (Path):
            Path to the lz-ani executable.

    Returns:
        list: The constructed command as a list of strings.

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

    if verbose: 
        cmd.extend(['--verbose', '2'])

    return cmd


def cmd_clusty(
        input_path: pathlib.Path,
        ids_path: pathlib.Path,
        output_path: pathlib.Path,
        algorithm: str,
        metric: str,
        tani: float,
        gani: float,
        ani: float,
        cov: float,
        num_alns: int,
        is_representatives: bool,
        leiden_resolution: float,
        bin_path=BIN_CLUSTY,
    ) -> typing.List[str]:
    """Constructs the command line for Clusty.

    Args:
        input_path (Path):
            Path to the input ANI file.
        ids_path (Path):
            Path to the input file with sequence identifiers.
        output_path (Path):
            Path to the output file.
        algorithm (str):
            Clustering algorithm.
        metric (str):
            Similarity metric for clustering.
        metric_threshold (float):
            Similarity threshold.
        tani (float):
            Minimum tANI.
        gani (float):
            Minimum gANI.
        ani (float):
            Minimum ANI.
        cov (float):
            Minimum coverage (aligned fraction).
        num_alns (int):
            Maximum number of local alignments between two genomes.
        is_representatives (bool):
            Whether to output a representative genome for each cluster.
        leiden_resolution (float):
            Resolution parameter for the Leiden algorithm.
        verbose (bool):
            Whether to display verbose output.
        bin_path (Path):
            Path to the clusty executable.

    Returns:
        list: The constructed command as a list of strings.
    
    """
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
    if algorithm == 'leiden':
        cmd.extend(['--leiden-resolution', f'{leiden_resolution}'])

    cmd.extend([f'{input_path}', f'{output_path}'])
    return cmd


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


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Initialize logger
    logger = create_logger(
        name='vclust',
        log_level=logging.INFO if args.verbose else logging.ERROR,
    )

    # Prefilter
    if args.command == 'prefilter':
        args.bin_kmerdb = validate_binary(args.bin_kmerdb)
        args = validate_args_prefilter(args, parser)
        args = validate_args_fasta_input(args, parser)

        out_dir = args.output_path.parent / get_uuid()
        out_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f'Creating a temp directory: {out_dir}')

        batches = []
        # Input is a directory of fasta files.
        if not args.is_multifasta:
            batches.append(args.fasta_paths)
        else:
            # Split multi-fasta file.
            if args.batch_size:
                args.bin_fastasplit = validate_binary(args.bin_fastasplit)
                cmd = cmd_fastasplit(
                    input_fasta=args.input_path, 
                    out_dir=out_dir,
                    n=args.batch_size,
                    verbose=args.verbose,
                    bin_path=args.bin_fastasplit,
                )
                p = run(cmd, args.verbose, logger)
                for f in out_dir.glob('part_*'):
                    batches.append([f])
                batches.sort()
            # Do not split multi-fasta file.
            else:
                batches.append([args.input_path])

        num_batches = len(batches)
        db_paths = []
        for i, batch in enumerate(batches):
            logger.info(f'Processing batch: {i+1} / {num_batches}')
            batch_id = f'part_{i:05d}' if num_batches > 1 else 'whole'
            txt_path = out_dir / f'{batch_id}.txt'
            db_path = out_dir / f'{batch_id}.kdb'

            # kmer-db build.
            cmd = cmd_kmerdb_build(
                input_paths=batch, 
                txt_path=txt_path,
                db_path=db_path,
                is_multisample_fasta=args.is_multifasta,
                kmer_size=args.k,
                num_threads=args.num_threads,
                verbose=args.verbose,
                bin_path=args.bin_kmerdb,
            )
            p = run(cmd, args.verbose, logger)
            db_paths.append(db_path)

            # Regardless of verbosity, always delete the partial FASTA file 
            # after building the corresponding partial k-mer database.
            if num_batches > 1:
                batch[0].unlink()

        # Run kmer-db all2all.
        db_list_path = out_dir / 'db_list.txt'
        all2all_path = out_dir / 'all2all.txt'

        cmd = cmd_kmerdb_all2all(
            db_paths=db_paths,
            db_list_path=db_list_path,
            outfile_all2all=all2all_path,
            kmer_count=args.min_kmers,
            num_threads=args.num_threads,
            verbose=args.verbose,
            bin_path=args.bin_kmerdb,
        )
        p = run(cmd, args.verbose, logger)

        cmd = cmd_kmerdb_distance(
            infile_all2all=all2all_path,
            min_value=args.min_ident,
            num_threads=args.num_threads,
            verbose=args.verbose,
            bin_path=args.bin_kmerdb,
        )
        p = run(cmd, args.verbose, logger)

        out_ani = out_dir / 'all2all.txt.ani-shorter'
        out_ani.rename(args.output_path)

        if not args.keep_temp:
            if out_dir.exists():
                logger.info(f'Removing directory: {out_dir}')
                shutil.rmtree(out_dir)

    # Align
    elif args.command == 'align':
        args.bin_lzani = validate_binary(args.bin_lzani)
        args = validate_args_fasta_input(args, parser)

        out_dir = args.output_path.parent / get_uuid()
        out_dir.mkdir(parents=True, exist_ok=True)
        txt_path = out_dir / 'ids.txt'

        logger.info(f'Creating temporary directory: {out_dir}')

        # Run lz-ani.
        cmd = cmd_lzani(
            input_paths=args.fasta_paths,
            txt_path=txt_path,
            output_path=args.output_path,
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
        p = run(cmd, args.verbose, logger)

        if out_dir.exists():
            logger.info(f'Removing directory: {out_dir}')
            shutil.rmtree(out_dir)

    # Cluster
    elif args.command == 'cluster':
        args.BIN_CLUSTY = validate_binary(args.BIN_CLUSTY)
        args = validate_args_cluster(args, parser)

        cmd = cmd_clusty(
            input_path=args.input_path,
            ids_path=args.ids_path,
            output_path=args.output_path,
            algorithm=args.algorithm,
            metric=args.metric,
            tani=args.tani,
            gani=args.gani,
            ani=args.ani,
            cov=args.cov,
            num_alns=args.num_alns,
            is_representatives=args.representatives,
            leiden_resolution=args.leiden_resolution,
            bin_path=args.BIN_CLUSTY,
        )
        p = run(cmd, args.verbose, logger)

if __name__ == '__main__':
    main()