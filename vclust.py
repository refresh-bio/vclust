#!/usr/bin/env python3
"""Compute Average Nucleotide Identity (ANI) and cluster virus genome sequences.

https://github.com/refresh-bio/vclust
"""

import argparse
import logging
import multiprocessing
import os
import pathlib
import platform
import shutil
import subprocess
import sys
import tempfile
import typing

__version__ = '1.3.1'

CITATION = '''Zielezinski A, Gudys A et al. (2025) 
Ultrafast and accurate sequence alignment and clustering of viral genomes. 
Nat Methods. doi: https://doi.org/10.1038/s41592-025-02701-7'''

DEFAULT_THREAD_COUNT = min(multiprocessing.cpu_count(), 64)

VCLUST_SCRIPT_PATH = pathlib.Path(__file__)
VCLUST_DIR = VCLUST_SCRIPT_PATH.resolve().parent

# Default paths to third-party binaries
BIN_DIR = VCLUST_DIR / 'bin'
BIN_KMERDB = BIN_DIR / 'kmer-db'
BIN_LZANI = BIN_DIR / 'lz-ani'
BIN_CLUSTY = BIN_DIR / 'clusty'
BIN_MFASTA = BIN_DIR / 'mfasta-tool'

# LZ-ANI output columns
ALIGN_FIELDS = [
    'qidx', 'ridx', 'query', 'reference',  'tani', 'gani', 'ani', 'qcov',
    'rcov', 'num_alns', 'len_ratio', 'qlen', 'rlen',  'nt_match', 'nt_mismatch', 
]
# Vclust align output formats
ALIGN_OUTFMT = {
    'lite': ALIGN_FIELDS[:2] + ALIGN_FIELDS[4:11],
    'standard': ALIGN_FIELDS[:11],
    'complete': ALIGN_FIELDS[:],
}

def get_parser() -> argparse.ArgumentParser:
    """Return an argument parser."""

    # Constants and shared utilities
    formatter_class = lambda prog: CustomHelpFormatter(
        prog, max_help_position=32, width=100
    )

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
    
    def add_verbosity_argument(parser):
        parser.add_argument(
            '-v',
            metavar='<int>',
            dest='verbosity_level',
            type=int,
            choices=[0, 1, 2],
            default=1,
            help='Verbosity level [%(default)s]:\n'
            '0: Errors only\n1: Info\n2: Debug'
        )

    def add_help_argument(parser):
        parser.add_argument(
            '-h', '--help',
            action='help',
            help='Show this help message and exit'
        )

    def add_threads_argument(parser):
        parser.add_argument(
            '-t', '--threads',
            metavar="<int>",
            dest="num_threads",
            type=int,
            default=DEFAULT_THREAD_COUNT,
            help='Number of threads [%(default)s]'
        )

    # Create the main parser
    parser = argparse.ArgumentParser(
        description=f'''%(prog)s v{__version__}: calculate ANI and cluster 
        virus (meta)genome sequences''',
        formatter_class=formatter_class,
        add_help=False,
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'v{__version__}',
        help="Display the tool's version and exit"
    )
    add_help_argument(parser)

    subparsers = parser.add_subparsers(dest='command')

    # Deduplicate parser
    deduplicate_epilog = [
        'Examples:',
        'vclust deduplicate -i refseq.fna genbank.fna -o nr.fna --add-prefixes',
        'vclust deduplicate -i refseq.fna genbank.fna -o nr.fna.gz --gzip-out '
        '--add-prefixes',
    ]
    deduplicate_parser = subparsers.add_parser(
        'deduplicate',
        help='Deduplicate and merge genome sequences from multiple FASTA files',
        formatter_class=formatter_class,
        add_help=False,
        epilog='\n'.join(deduplicate_epilog)
    )
    deduplicate_optional = deduplicate_parser._action_groups.pop()
    deduplicate_required = deduplicate_parser.add_argument_group(
        'required arguments'
    )
    deduplicate_parser._action_groups.append(deduplicate_optional)

    deduplicate_required.add_argument(
        '-i', '--in',
        metavar='<file>',
        type=input_path_type,
        dest='input_path',
        nargs='+',
        help='Space-separated input FASTA files (gzipped or uncompressed)',
    )
    deduplicate_required.add_argument(
        '-o', '--out',
        metavar='<file>',
        type=pathlib.Path,
        dest='output_path',
        help='Output FASTA file with unique, nonredundant sequences, removing '
        'duplicates and their reverse complements',
        required=True
    )
    deduplicate_parser.add_argument(
        '--add-prefixes',
        metavar='<str>',
        nargs='*',
        default=False,
        help='Add prefixes to sequence IDs. Without any values, prefixes are '
        'set to input file names. Provide space-separated prefixes to override.'
        ' Default: prefixes are disabled [%(default)s]'
    )
    deduplicate_parser.add_argument(
        '--gzip-output',
        action="store_true",
        help='Compress the output FASTA file with gzip [%(default)s]'
    )
    deduplicate_parser.add_argument(
        '--gzip-level',
        metavar='<int>',
        type=int,
        default=4,
        help='Compression level for gzip (1-9) [%(default)s]'
    )
    add_threads_argument(deduplicate_parser)
    add_verbosity_argument(deduplicate_parser)
    add_help_argument(deduplicate_parser)

    # Prefilter parser    
    prefilter_parser = subparsers.add_parser(
        'prefilter',
        help='Prefilter genome pairs for alignment',
        formatter_class=formatter_class,
        add_help=False,
    )

    prefilter_optional = prefilter_parser._action_groups.pop()
    prefilter_required = prefilter_parser.add_argument_group(
        'required arguments'
    )
    prefilter_parser._action_groups.append(prefilter_optional)

    prefilter_required.add_argument(
        '-i', '--in',
        metavar='<file>',
        type=input_path_type,
        dest='input_path',
        help='Input FASTA file or directory of files (gzipped or uncompressed)',
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
        default=20,
        help='Minimum number of shared k-mers between two genomes [%(default)s]'
    )
    prefilter_parser.add_argument(
        '--min-ident',
        metavar="<float>",
        type=ranged_float_type,
        default=0.7,
        help='Minimum sequence identity (0-1) between two genomes. Calculated '
        'based on the shorter sequence [%(default)s]'
    )
    prefilter_parser.add_argument(
        '--batch-size',
        metavar="<int>",
        type=int,
        default=0,
        help='Process a multifasta file in smaller batches of <int> FASTA '
        'sequences. This option reduces memory at the expense of speed. '
        'By default, no batch [%(default)s]'
    )
    prefilter_parser.add_argument(
        '--kmers-fraction',
        metavar="<float>",
        type=ranged_float_type,
        default=1.0,
        help='Fraction of k-mers to analyze in each genome (0-1). A lower '
        'value reduces RAM usage and speeds up processing. By default, all '
        'k-mers [%(default)s]'
    )
    prefilter_parser.add_argument(
        '--max-seqs',
        metavar="<int>",
        type=int,
        default=0,
        help='Maximum number of sequences allowed to pass the prefilter per '
        'query. Only the sequences with the highest identity to the query are '
        'reported. This option reduces RAM usage and speeds up processing '
        '(affects sensitivity). By default, all sequences that pass the '
        'prefilter are reported [%(default)s]'
    )
    add_threads_argument(prefilter_parser)
    add_verbosity_argument(prefilter_parser)
    add_help_argument(prefilter_parser)

    # Align parser
    align_parser = subparsers.add_parser(
        'align',
        help='Align genome sequence pairs and calculate ANI measures',
        formatter_class=formatter_class,
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
        help='Input FASTA file or directory of files (gzipped or uncompressed)',
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
        '--out-aln',
        metavar='<file>',
        type=pathlib.Path,
        dest='aln_path',
        help='Write alignments to the tsv <file>',
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
        '--out-qcov',
        dest='qcov',
        metavar='<float>',
        type=ranged_float_type,
        default=0,
        help='Min. query coverage (aligned fraction) to output (0-1) '
        '[%(default)s]'
    )
    align_parser.add_argument(
        '--out-rcov',
        dest='rcov',
        metavar='<float>',
        type=ranged_float_type,
        default=0,
        help='Min. reference coverage (aligned fraction) to output (0-1) '
        '[%(default)s]'
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
    add_threads_argument(align_parser)
    add_verbosity_argument(align_parser)
    add_help_argument(align_parser)

    # Cluster parser
    cluster_parser = subparsers.add_parser(
        'cluster',
        help='Cluster genomes based on ANI thresholds',
        formatter_class=formatter_class,
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
        help='Min. total ANI (0-1) [%(default)s]'
    )
    cluster_parser.add_argument(
        '--gani',
        metavar='<float>',
        dest='gani',
        type=ranged_float_type,
        default=0,
        help='Min. global ANI (0-1) [%(default)s]'
    )
    cluster_parser.add_argument(
        '--ani',
        metavar='<float>',
        dest='ani',
        type=ranged_float_type,
        default=0,
        help='Min. ANI (0-1) [%(default)s]'
    )
    cluster_parser.add_argument(
        '--qcov',
        metavar='<float>',
        dest='qcov',
        type=ranged_float_type,
        default=0,
        help='Min. query coverage/aligned fraction (0-1) [%(default)s]'
    )
    cluster_parser.add_argument(
        '--rcov',
        metavar='<float>',
        dest='rcov',
        type=ranged_float_type,
        default=0,
        help='Min. reference coverage/aligned fraction (0-1) [%(default)s]'
    )
    cluster_parser.add_argument(
        '--len_ratio',
        metavar='<float>',
        dest='len_ratio',
        type=ranged_float_type,
        default=0,
        help='Min. length ratio between shorter and longer sequence (0-1) '
        '[%(default)s]'
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
        help='Resolution parameter for the Leiden algorithm [%(default)s]'
    )
    cluster_parser.add_argument(
        '--leiden-beta',
        metavar='<float>',
        type=ranged_float_type,
        default=0.01,
        help='Beta parameter for the Leiden algorithm [%(default)s]'
    )
    cluster_parser.add_argument(
        '--leiden-iterations',
        metavar='<int>',
        type=int,
        default=2,
        help='Number of iterations for the Leiden algorithm [%(default)s]'
    )
    add_verbosity_argument(cluster_parser)
    add_help_argument(cluster_parser)

    # Info parser
    info_parser = subparsers.add_parser(
        'info',
        help='Show information about the tool and its dependencies',
        formatter_class=formatter_class,
        add_help=False,
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
        ('deduplicate', deduplicate_parser),
    ]
    for name, subparser in subparsers:
        if sys.argv[-1] == name:
            subparser.print_help()
            parser.exit()

    return parser


def create_logger(name: str, verbosity_level: int) -> logging.Logger:
    """Creates and configures a logger with the specified verbosity level.

    Args:
        name (str): Name of the logger.
        verbosity_level (int): Verbosity level:
            - 0: ERROR (default)
            - 1: INFO
            - 2: DEBUG

    Returns:
        logging.Logger: Configured logger instance.
    """
    # Map verbosity levels to logging levels
    verbosity_to_log_level = {
        0: logging.ERROR,
        1: logging.INFO,
        2: logging.DEBUG
    }
    log_level = verbosity_to_log_level.get(verbosity_level, logging.ERROR)

    # Create logger
    logger = logging.getLogger(name)
    logger.setLevel(log_level)

    # Prevent duplicate handlers in case the logger is reconfigured
    if not logger.handlers:
        # Create a stream handler with a custom formatter
        handler = logging.StreamHandler()
        handler.setLevel(log_level)
        handler.setFormatter(CustomLoggerFormatter())
        logger.addHandler(handler)

    return logger


def _validate_binary(bin_path: pathlib.Path) -> pathlib.Path:
    """Validates the presence and executability of a binary file.

    This function checks if the provided path points to an existing binary file
    and if it is executable. It also attempts to run the binary to ensure it
    operates without errors.

    Args:
        bin_path:
            The path to the executable binary file.

    Returns:
        pathlib.Path: The resolved path to the binary file.

    Raises:
        RuntimeError: If the binary file does not exist, is not executable, 
                      or if running the binary encounters an error.
    """
    bin_path = bin_path.resolve()

    if not bin_path.exists():
        raise RuntimeError(f'File not found: {bin_path}')
    
    if not bin_path.is_file() or not os.access(bin_path, os.X_OK):
        raise RuntimeError(f'Binary file not executable: {bin_path}')
    
    try:
        subprocess.run(
            [str(bin_path)],  
            stdout=subprocess.DEVNULL, 
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f'Running {bin_path} failed: {e.stderr}')
    except OSError as e:
        raise RuntimeError(f'OSError in {bin_path}: {e}')
    except Exception as e:
        raise RuntimeError(f'Unexpected error in {bin_path}: {e}')
    return bin_path


def validate_binary(bin_path: pathlib.Path) -> pathlib.Path:
    try:
        return _validate_binary(bin_path)
    except RuntimeError as e:
        sys.exit(f'error: {e}')


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
                     f'Expected at least 2, found {len(args.fasta_paths)}.')

    return args


def validate_args_deduplicate(args, parser) -> argparse.Namespace:
    """Validates the arguments for the deduplicate command."""
    # Validate dataset prefixes
    if args.add_prefixes:
        if len(args.add_prefixes) != len(args.input_path):
            parser.error(
                'Number of prefixes must match the number of input files.')
        if any(',' in prefix for prefix in args.add_prefixes):
            parser.error('Prefixes cannot contain commas.')
    elif args.add_prefixes == []:
        args.add_prefixes = [f'{f.stem.split(".")[0]}|' for f in args.input_path]
    
    if not (1 <= args.gzip_level <= 9):
        parser.error('Compression level must be between 1 and 9.')

    # Ensure the output file has a .gz extension if gzip compression is enabled.
    if args.gzip_output and not args.output_path.suffix == '.gz':
        args.output_path = pathlib.Path(f'{args.output_path}.gz')

    # Create the output file name for storing information on duplicates.
    args.output_duplicates_path = pathlib.Path(
        f'{args.output_path}.duplicates.txt')

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
        if 'qidx' not in header and 'ridx' not in header:
            parser.error(
                f'missing columns `qidx` and `ridx` in {args.input_path}')
        cols = ['tani', 'gani', 'ani', 'qcov', 'rcov', 'len_ratio', 'num_alns']
        for name in cols:
            value = args_dict[name]
            if value != 0 and name not in header:
                parser.error(f'missing column `{name}` in {args.input_path}')
    return args


def run(
        cmd: typing.List[str],
        verbosity_level: int,
        logger: logging.Logger
    ) -> subprocess.CompletedProcess:
    """Executes a given command as a subprocess and handles logging.

    This function runs the specified command, logs the execution details,
    and manages errors. If verbose mode is enabled, the command's standard
    error output is not suppressed. Otherwise, the standard error is piped
    and logged in case of failure.

    Args:
        cmd:
            The command to run as a list of strings.
        verbose: 
            Flag indicating whether to run the command in verbose mode.
        logger:
            The logger instance for logging information and errors.

    Returns:
        subprocess.CompletedProcess: The completed process information.

    Raises:
        SystemExit: If the command fails to execute or an error occurs.
    """
    logger.info(f'Running: {" ".join(cmd)}')
    try:
        process = subprocess.run(
            cmd,  
            stdout=None if verbosity_level else subprocess.DEVNULL, 
            stderr=None if verbosity_level else subprocess.PIPE,
            text=True,
            check=True
        )
    except subprocess.CalledProcessError as e:
        logger.error(f'Process {" ".join(cmd)} failed with message: {e.stderr}')
        sys.exit(1)
    except OSError as e:
        logger.error(f'OSError: {" ".join(cmd)} failed with message: {e}')
        sys.exit(1)
    except Exception as e:
        logger.error(f'Unexpected: {" ".join(cmd)} failed with message: {e}')
        sys.exit(1)  
    logger.info(f'Completed')
    return process


def cmd_mfasta_deduplicate(
        input_paths: typing.List[pathlib.Path],
        output_path: pathlib.Path,
        output_duplicates_path: pathlib.Path,
        add_prefixes: typing.List[pathlib.Path],
        gzip_output: bool,
        gzip_level: int,
        num_threads: int,
        verbosity_level: int,
        bin_path: pathlib.Path = BIN_MFASTA
    ) -> typing.List[str]:
    """Builds the command for mfasta-tool to deduplicate FASTA files.
    
    Args:
        input_paths (List[Path]):
            List of paths to the input FASTA files.
        output_path (Path):
            Path to the output FASTA file.
        output_duplicates_path (Path):
            Path to the output file listing duplicate sequences.
        add_prefixes (List[Path]):
            List of prefixes to add to sequence IDs.
        gzip_output (bool):
            Whether to compress the output FASTA file.
        gzip_level (int):
            Compression level for gzip (1-9).
        num_threads (int):
            Number of threads to use in mfasta-tool.
        bin_path (Path):
            Path to the mfasta-tool executable.
    
    Returns:
        list: The constructed command as a list of strings.

    """
    cmd = [
        f'{bin_path}',
        'mrds',
        '-i',
        ",".join([str(f) for f in input_paths]),
        '-o',
        f'{output_path}',
        '--out-duplicates',
        f'{output_duplicates_path}',
        '--remove-duplicates',
        '--mark-duplicates-orientation',
        '--rev-comp-as-equivalent',
        '-t',
        f'{num_threads}',
    ]
    if verbosity_level:
        cmd.extend(['--verbosity', '1'])
    if add_prefixes:
        cmd.extend(['--in-prefixes', ",".join(add_prefixes)])
    if gzip_output:
        cmd.extend(['--gzipped-output', '--gzip-level', f'{gzip_level}'])
    return cmd


def cmd_mfasta_split(
        input_fasta: pathlib.Path,
        out_dir: pathlib.Path,
        n: int,
        num_threads: int,
        verbosity_level: int,
        bin_path = BIN_MFASTA
    ) -> typing.List[str]:
    """Builds the command for mfasta-tool to split the input FASTA file. 
    
    Args:
        input_fasta (Path):
            Path to the input FASTA file.
        out_dir (Path):
            Path to the output directory.
        n (int):
            Number of sequences per output FASTA file.
        num_threads (int):
            Number of threads
        bin_path (Path):
            Path to the mfasta-tool executable.
        
    Returns:
        list: The constructed command as a list of strings.
        
    """
    cmd = [
        f'{bin_path}',
        'mrds',
        '-n',
        f'{n}',
        '--out-prefix',
        f'{out_dir}/part',
        '--out-suffix',
        'fna.gz',
        '--gzipped-output',
        '-i',
        f'{input_fasta}',
        '-t',
        f'{num_threads}',
    ]
    if verbosity_level:
        cmd.extend(['--verbosity', '1'])
    return cmd 


def cmd_kmerdb_build(
        input_paths: pathlib.Path,
        txt_path: pathlib.Path,
        db_path: pathlib.Path,
        is_multisample_fasta: bool,
        kmer_size: int,
        kmers_fraction: float,
        num_threads: int,
        bin_path: pathlib.Path = BIN_KMERDB
    ) -> typing.List[str]:
    """Builds the command for Kmer-db build.

    Args:
        input_fasta (Path):
            Path to the input FASTA file or directory with input FASTA files.
        outfile_txt (Path):
            Path to the output text file that will list the input FASTA files.
        outfile_db (Path):
            Path to the output kmer-db database file.
        kmer_size (int):
            k-mer size.
        kmers_fraction (float):
            Fraction of k-mers to analyze for each genome (0-1).
        num_threads (int):
            Number of threads to use in kmer-db.
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
        "-f", f"{kmers_fraction}",
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
        min_kmers: int,
        min_ident: float,
        max_seqs: int,
        num_threads: int,
        bin_path: pathlib.Path = BIN_KMERDB
    ) -> typing.List[str]:
    """Builds the command for Kmer-db all2all.

    Args:
        db_paths (list[Path]):
            List of paths to the input kmer-db database files.
        db_list_path (Path):
            Path to the output text file listing the kmer-db database files.
        outfile_all2all (Path):
            Path to the output all2all file.
        min_kmers (int):
            Minimum number of shared k-mers to report in all2all output.
        min_ident (float):
            Minimum sequence identity of the shorter sequence.
        max_seqs (int):
            Maximum number of sequences allowed to pass the prefilter per query.
        num_threads (int):
            Number of threads to use in kmer-db.
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
        'all2all-parts' if len(db_paths) > 1 else 'all2all-sp',
        '-sparse',
        '-min', f'num-kmers:{min_kmers}',
        '-min', f'ani-shorter:{min_ident}',
        "-t", f"{num_threads}",
        f'{db_list_path}' if len(db_paths) > 1 else f'{db_paths[0]}',
        f'{outfile_all2all}',
    ]
    if max_seqs > 0:
        cmd[5:5] = ['-sample-rows', f'ani-shorter:{max_seqs}']
    return cmd


def cmd_kmerdb_distance(
        infile_all2all: pathlib.Path,
        outfile_distance: pathlib.Path,
        min_ident: float,
        num_threads: int,
        bin_path: pathlib.Path = BIN_KMERDB
    ) -> typing.List[str]:
    """Builds the command for Kmer-db distance.

    Args:
        infile_all2all (Path):
            Path to the input all2all file.
        outfile_distance (Path):
            Path to the output distance (file) file.
        min_ident (float):
            Minimum sequence identity to report in output.
        num_threads (int):
            Number of threads to use in kmer-db.
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
        '-min', f'{min_ident}',
        "-t", f"{num_threads}",
        f'{infile_all2all}',
        f'{outfile_distance}',
    ]
    return cmd


def cmd_lzani(
        input_paths: typing.List[pathlib.Path],
        txt_path: pathlib.Path,
        output_path: pathlib.Path,
        out_format: typing.List[str],
        out_aln_path: pathlib.Path,
        out_tani: float,
        out_gani: float,
        out_ani: float,
        out_qcov: float,
        out_rcov: float,
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
        verbosity_level: bool,
        bin_path: pathlib.Path = BIN_LZANI
    ) -> typing.List[str]:
    """Builds the command for LZ-ANI.

    Args:
        input_paths (List[Path]):
            List of paths to the input FASTA files.
        txt_path (Path):
            Path to the output text file listing the input FASTA files.
        output_path (Path):
            Path to the output ANI file.
        out_format (List[str]):
            List of LZ-ANI column names.
        out_aln_path (Path):
            Path to the output alignment file.
        out_tani (float):
            Minimum tANI to output.
        out_gani (float):
            Minimum gANI to output.
        out_ani (float):
            Minimum ANI to output.
        out_qcov (float):
            Minimum query coverage (aligned fraction) to output.
        out_rcov (float):
            Minimum reference coverage (aligned fraction) to output.
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
    if out_aln_path:
        cmd.extend(['--out-alignment', f'{out_aln_path}'])

    cols = [
        ('tani', out_tani), ('gani', out_gani), ('ani', out_ani), 
        ('qcov', out_qcov), ('rcov', out_rcov)
    ]
    for name, value in cols:
        if value > 0:
            cmd.extend(['--out-filter', f'{name}', f'{value}'])

    if verbosity_level: 
        cmd.extend(['--verbose', f'{verbosity_level + 1}'])

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
        qcov: float,
        rcov: float,
        num_alns: int,
        len_ratio: float,
        is_representatives: bool,
        leiden_resolution: float,
        leiden_beta: float,
        leiden_iterations: int,
        bin_path=BIN_CLUSTY,
    ) -> typing.List[str]:
    """Builds the command for Clusty.

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
        qcov (float):
            Minimum query coverage (aligned fraction).
        rcov (float):
            Minimum reference coverage (aligned fraction).
        num_alns (int):
            Maximum number of local alignments between two genomes.
        len_ratio (float):
            Ratio between genome sequence lengths (shorter / longer)
        is_representatives (bool):
            Whether to output a representative genome for each cluster.
        leiden_resolution (float):
            Resolution parameter for the Leiden algorithm.
        leiden_beta (float):
            Beta parameter for the Leiden algorithm.
        leiden_iterations (int):
            Number of iterations for the Leiden algorithm.
        bin_path (Path):
            Path to the Clusty executable.

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
        'qidx', 'ridx',
        '--distance-col',
        f'{metric}',
        '--similarity',
        '--numeric-ids',
    ]
    cols = [('tani', tani), ('gani', gani), ('ani', ani), ('qcov', qcov), 
    ('rcov', rcov), ('len_ratio', len_ratio)]
    for name, value in cols:
        if value > 0:
            cmd.extend(['--min', f'{name}', f'{value}'])
    if num_alns > 0:
        cmd.extend(['--max', 'num_alns', f'{num_alns}'])
    if is_representatives:
        cmd.append('--out-representatives')
    if algorithm == 'leiden':
        leiden_options = [
            '--leiden-resolution', f'{leiden_resolution}',
            '--leiden-beta', f'{leiden_beta}',
            '--leiden-iterations', f'{leiden_iterations}',
        ]
        cmd.extend(leiden_options)

    cmd.extend([f'{input_path}', f'{output_path}'])
    return cmd


def handle_info(args, parser, logger) -> None:
    """    
    Displays the Vclust version, installation paths, and binary dependencies. 
    Checks for the presence and executable status of required binaries.

    Exits with a non-zero status if any dependencies are missing or 
    not executable.

    Returns:
        None
    
    Raises:
        SystemExit: If any binary dependencies are missing or not executable.
    
    """
    binaries = {
        'Kmer-db': BIN_KMERDB,
        'LZ-ANI': BIN_LZANI,
        'Clusty': BIN_CLUSTY,
        'mfasta': BIN_MFASTA,
    }

    vclustpy_size_kb = f'{VCLUST_SCRIPT_PATH.stat().st_size / 1024:.1f} KB'
    output_lines = [
        f'Vclust version {__version__} (Python {platform.python_version()})',
        '',
        'Citation:',
        CITATION,
        '',
        'Installed at:',
        f'   {VCLUST_SCRIPT_PATH.resolve()} ({vclustpy_size_kb})',
        f'   {BIN_DIR.resolve()}',
        '',
        'Binary dependencies:',
    ]

    errors = []  # List to collect any errors encountered during binary checks.

    # Check each binary's presence and version.
    for name, path in binaries.items():
        try:
            _validate_binary(path)
            version = subprocess.run(
                [str(path), '-version' if name == 'Kmer-db' else '--version'],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            ).stderr.strip()
            size_mb = f'{path.stat().st_size / 1024 / 1024:.1f} MB'
            output_lines.append(f'   {name:<20} v{version:<7} ({size_mb})')
        except Exception as e:
            output_lines.append(f'   {name:<20} [error]')
            errors.append((name, e))

    # Append the status summary based on any encountered errors.
    output_lines.append('')
    if errors:
        output_lines.append(f'\033[31mStatus: error\033[0m')
        output_lines.extend(f"   - {name}: {error}" for name, error in errors)
    else:
        output_lines.append(f'\033[32;1mStatus: ready\033[0m')

    # Output the complete information.
    print('\n'.join(output_lines))

    if errors:
        sys.exit(1)


def handle_deduplicate(args, parser, logger) -> None:
    """
    Deduplicates input FASTA files by validating inputs, constructing commands, 
    and executing the deduplication process.

    Args:
        args (argparse.Namespace):
            Parsed command-line arguments.
        parser (argparse.ArgumentParser):
            Argument parser for validation.
        logger (logging.Logger):
            Logger instance for logging messages.

    """
    validate_binary(BIN_MFASTA)
    args = validate_args_deduplicate(args, parser)
    cmd = cmd_mfasta_deduplicate(
        input_paths=args.input_path,
        output_path=args.output_path,
        output_duplicates_path=args.output_duplicates_path,
        add_prefixes=args.add_prefixes,
        gzip_output=args.gzip_output,
        gzip_level=args.gzip_level,
        verbosity_level=args.verbosity_level,
        num_threads=args.num_threads,
    )
    run(cmd, args.verbosity_level, logger)


def handle_prefilter(args, parser, logger) -> None:
    """
    Preprocesses FASTA inputs by building Kmer-db databases and computing 
    pairwise distances between sequences.

    Args:
        args (argparse.Namespace):
            Parsed command-line arguments.
        parser (argparse.ArgumentParser):
            Argument parser for validation.
        logger (logging.Logger):
            Logger instance for logging messages.

    """
    validate_binary(BIN_KMERDB)
    args = validate_args_prefilter(args, parser)
    args = validate_args_fasta_input(args, parser)

    out_dir = args.output_path.parent

    with tempfile.TemporaryDirectory(dir=out_dir) as temp_dir:
        out_dir = pathlib.Path(temp_dir)
        logger.info(f'Temp directory created: {out_dir}')
        batches = []
        # Input is a directory of fasta files.
        if not args.is_multifasta:
            batches.append(args.fasta_paths)
        else:
            # Split multi-fasta file.
            if args.batch_size:
                validate_binary(BIN_MFASTA)
                cmd = cmd_mfasta_split(
                    input_fasta=args.input_path, 
                    out_dir=out_dir,
                    n=args.batch_size,
                    verbosity_level=args.verbosity_level,
                    num_threads=args.num_threads,
                )
                run(cmd, args.verbosity_level, logger)
                batches.extend(sorted([f] for f in out_dir.glob('part_*')))
            # Do not split multi-fasta file.
            else:
                batches.append([args.input_path])

        db_paths = []
        num_batches = len(batches)
        for i, batch in enumerate(batches):
            logger.info(f'Processing batch: {i+1}/{num_batches}')
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
                kmers_fraction=args.kmers_fraction,
                num_threads=args.num_threads,
            )
            run(cmd, args.verbosity_level, logger)
            db_paths.append(db_path)

            # Delete the partial FASTA file after building the corresponding
            # partial Kmer-db database.
            if num_batches > 1:
                batch[0].unlink()

        # Run kmer-db all2all.
        db_list_path = out_dir / 'db_list.txt'
        all2all_path = out_dir / 'all2all.txt'
        cmd = cmd_kmerdb_all2all(
            db_paths=db_paths,
            db_list_path=db_list_path,
            outfile_all2all=all2all_path,
            min_kmers=args.min_kmers,
            min_ident=args.min_ident,
            max_seqs=args.max_seqs,
            num_threads=args.num_threads,
        )
        run(cmd, args.verbosity_level, logger)

        # Run kmer-db distance.
        cmd = cmd_kmerdb_distance(
            infile_all2all=all2all_path,
            outfile_distance=args.output_path,
            min_ident=args.min_ident,
            num_threads=args.num_threads,
        )
        run(cmd, args.verbosity_level, logger)    


def handle_align(args, parser, logger) -> None:
    """Aligns input sequences using LZ-ANI and outputs alignment metrics.

    Args:
        args (argparse.Namespace):
            Parsed command-line arguments.
        parser (argparse.ArgumentParser):
            Argument parser for validation.
        logger (logging.Logger):
            Logger instance for logging messages.
    
    """
    validate_binary(BIN_LZANI)
    args = validate_args_fasta_input(args, parser)

    out_dir = args.output_path.parent

    with tempfile.TemporaryDirectory(dir=out_dir) as temp_dir:
        out_dir = pathlib.Path(temp_dir)
        logger.info(f'Temp directory created: {out_dir}')
        txt_path = out_dir / 'ids.txt'

        # Run lz-ani.
        cmd = cmd_lzani(
            input_paths=args.fasta_paths,
            txt_path=txt_path,
            output_path=args.output_path,
            out_format=ALIGN_OUTFMT[args.outfmt],
            out_aln_path=args.aln_path,
            out_tani=args.tani,
            out_gani=args.gani,
            out_ani=args.ani,
            out_qcov=args.qcov,
            out_rcov=args.rcov,
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
            verbosity_level=args.verbosity_level,
        )
        run(cmd, args.verbosity_level, logger)


def handle_cluster(args, parser, logger):
    """Clusters sequences based on input metrics and parameters using the 
    specified clustering algorithm.

    Args:
        args (argparse.Namespace):
            Parsed command-line arguments.
        parser (argparse.ArgumentParser):
            Argument parser for validation.
        logger (logging.Logger):
            Logger instance for logging messages.
    """

    validate_binary(BIN_CLUSTY)
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
        qcov=args.qcov,
        rcov=args.rcov,
        num_alns=args.num_alns,
        len_ratio=args.len_ratio,
        is_representatives=args.representatives,
        leiden_resolution=args.leiden_resolution,
        leiden_beta=args.leiden_beta,
        leiden_iterations=args.leiden_iterations,
    )
    run(cmd, args.verbosity_level, logger)


class CustomLoggerFormatter(logging.Formatter):
    """Custom logging formatter with log messages with different colors"""

    FMT = '{asctime} [{levelname:^7}] {message}'
    FORMATS = {
        logging.DEBUG: FMT,
        logging.INFO: f"\33[36m{FMT}\33[0m",
        logging.WARNING: f"\33[33m{FMT}\33[0m",
        logging.ERROR: f"\33[31m{FMT}\33[0m",
    }

    def format(self, record):
        log_fmt = self.FORMATS[record.levelno]
        formatter = logging.Formatter(log_fmt, style="{")
        return formatter.format(record)


class CustomHelpFormatter(argparse.RawDescriptionHelpFormatter):
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
    logger = create_logger('Vclust', getattr(args, 'verbosity_level', 0))
    
    # Command handlers
    command_map = {
        'info': handle_info,
        'deduplicate': handle_deduplicate,
        'prefilter': handle_prefilter,
        'align': handle_align,
        'cluster': handle_cluster,
    }

    # Dispatch the appropriate command
    if args.command in command_map:
        command_map[args.command](args, parser, logger)


if __name__ == '__main__':
    main()
