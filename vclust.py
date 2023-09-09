#!/usr/bin/env python3
"""Compute Average Nucleotide Identity (ANI) and cluster virus genome sequences.

https://github.com/refresh-bio/vclust-dev
"""

import argparse
import csv
import multiprocessing
import subprocess
from pathlib import Path
import sys
import uuid
import typing

__version__ = '0.1'

DEFAULT_THREAD_COUNT = min(multiprocessing.cpu_count(), 32)

SCRIPT_DIR = Path(__file__).resolve().parent
KMERDB_EXEC = SCRIPT_DIR / 'bin' / 'kmer-db'
LZANI_EXEC = SCRIPT_DIR / 'bin' / 'lz-ani'


def get_parser() -> argparse.ArgumentParser:
    """Returns an argument parser."""
    parser = argparse.ArgumentParser(
        description=f'%(prog)s v.{__version__}: calculate ANI and cluster '
                     'virus (meta)genome sequences'
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'v.{__version__}',
        help="display the tool's version and exit"
    )

    subparsers = parser.add_subparsers(dest='command')

    # Align subparser
    align_parser = subparsers.add_parser(
        'align',
        help='prefilter and align virus genomic sequences',
        formatter_class=fmt
    )
    align_parser.add_argument(
        'input',
        help='path to input fasta file or directory with input fasta files'
    )
    align_parser.add_argument(
        'output',
        help='path to output alignment file'
    )
    align_parser.add_argument(
        '-k', '--kmer_size',
        metavar="INT",
        dest="kmer_size",
        type=int,
        default=18,
        help="k-mer size to be used for prefiltering in kmer-db [%(default)s]"
    )
    align_parser.add_argument(
        '-c', '--kmer_count',
        metavar="INT",
        dest="kmer_count",
        type=int,
        default=3,
        help='minimum number of shared k-mers required to start alignment '
             '[%(default)s]'
    )
    align_parser.add_argument(
        '-t', '--threads',
        metavar="INT",
        dest="num_threads",
        type=int,
        default=DEFAULT_THREAD_COUNT,
        help='number of threads to use (all by default) [%(default)s]'
    )
    align_parser.add_argument(
        '-v', '--verbose',
        action="store_true",
        help="show progress"
    )
    align_parser.add_argument(
        '--kmerdb_exec',
        metavar='FILE',
        dest="kmerdb_exec",
        default=f'{KMERDB_EXEC.relative_to(SCRIPT_DIR)}',
        help='path to kmer-db executable [%(default)s]'
    )
    align_parser.add_argument(
        '--lzani_exec',
        metavar='FILE',
        dest="lzani_exec",
        default=f'{LZANI_EXEC.relative_to(SCRIPT_DIR)}',
        help='path to lz-ani executable [%(default)s]'
    )

    # Calculate ANI subparser
    calcani_parser = subparsers.add_parser(
        'calcani',
        help='calculate ANI from sequence alignments',
        formatter_class=fmt,
    )
    calcani_parser.add_argument(
        'input',
        help='path to the input file with sequence alignments'
    )
    calcani_parser.add_argument(
        'output',
        help='path to the output file with ANI values'
    )

    # Cluster subparser
    cluster_parser = subparsers.add_parser(
        'cluster',
        help='cluster virus genomes based on ANI values',
        formatter_class=fmt
    )
    cluster_parser.add_argument(
        'input',
        help='path to the input file with ANI values'
    )
    cluster_parser.add_argument(
        'output',
        help='path to the output file with virus clusters'
    )
    cluster_parser.add_argument(
        '-m', '--method',
        metavar="INT",
        dest="method",
        type=int,
        default=1,
        help='select clustering method (default: %(default)s)\n'
             '1: single-linkage clustering\n'
             '2: complete-linkage clustering\n'
             '3: UCLUST-like clustering\n'
             '4: set-cover (greedy) mmseqs clustering\n'
             '5: connected component (BLASTclust)'
    )
    cluster_parser.add_argument(
        '-t', '--threads',
        metavar="INT",
        dest="num_threads",
        type=int,
        default=DEFAULT_THREAD_COUNT,
        help='number of threads to use (default: %(default)s)'
    )

    # Display help message if the script is executed without any arguments.
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    return parser


def validate_executable(exec_path: typing.Union[str, Path]) -> Path:
    """Checks if an executable file exists and is valid.

    Args:
        exec_path (str or Path): The path to the executable.

    Returns:
        Path to the executable.

    Raises:
        SystemExit: If the executable file is not found.

    """    
    exec_path = Path(exec_path) if isinstance(exec_path, str) else exec_path
    try:
        subprocess.run(
            [f'{exec_path}'],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        exit(f'error: executable not found: {exec_path.resolve()}')
    return exec_path


def validate_process(process: subprocess.CompletedProcess):
    """Validates the result of a subprocess execution.

    This function checks the return code of a completed subprocess and, 
    if it's non-zero, prints an error message and exits the script with
    an exit code of 1.

    Args:
        process (subprocess.CompletedProcess): 
            The result of a completed subprocess.

    Returns:
        subprocess.CompletedProcess: The same completed process information.

    Raises:
        SystemExit: If the subprocess return code is non-zero.

    """    
    if process.returncode:
        print(f'error while running: {" ".join(process.args)}')
        print(f'error message: {process.stderr}')
        exit(1)
    return process


def run_kmerdb_build(
        input_fasta: Path,
        outfile_txt: typing.Union[Path, str],
        outfile_db: typing.Union[Path, str],
        kmer_size: int,
        num_threads: int,
        verbose: bool = False,
        kmerdb_exec: typing.Union[Path, str] = KMERDB_EXEC
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
        kmerdb_exec (Path or str):
            Path to the kmer-db executable.

    Returns:
        subprocess.CompletedProcess: Completed process information.

    """
    # Create a text file listing input FASTA files.
    with open(outfile_txt, 'w') as oh:
        if input_fasta.is_dir():
            is_multifasta = False
            for f in input_fasta.iterdir():
                oh.write(f'{f}\n')
        else:
            is_multifasta = True
            oh.write(f'{input_fasta}')

    # Run kmer-db build.
    cmd = [
        f"{kmerdb_exec}", 
        "build",
        "-k", f"{kmer_size}",    
        "-t", f"{num_threads}",
        f'{outfile_txt}',
        f'{outfile_db}',
    ]
    if is_multifasta:
        cmd.insert(1, '-multisample-fasta')
    process = subprocess.run(
        cmd,  
        stdout=None if verbose else subprocess.DEVNULL, 
        stderr=subprocess.PIPE,
        text=True,
    )
    return process


def run_kmerdb_all2all(
        infile_db: typing.Union[Path, str],
        outfile_all2all: typing.Union[Path, str],
        kmer_count: int,
        num_threads: int,
        verbose: bool = False,
        kmerdb_exec: typing.Union[Path, str] = KMERDB_EXEC
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
        kmerdb_exec (Path or str):
            Path to the kmer-db executable.

    Returns:
        subprocess.CompletedProcess: Completed process information.

    """
    cmd = [
        f"{kmerdb_exec}", 
        "all2all", 
        "-sparse",
        '-above', f'{kmer_count}',
        "-t", f"{num_threads}",
        f'{infile_db}',
        f'{outfile_all2all}',
    ]
    process = subprocess.run(
        cmd,  
        stdout=None if verbose else subprocess.DEVNULL, 
        stderr=subprocess.PIPE,
        text=True,
    )
    return process


def run_lzani(
        input_fasta: Path,
        infile_txt: typing.Union[Path, str],
        infile_all2all: typing.Union[Path, str],
        outfile_aln: typing.Union[Path, str],
        kmer_count: int,
        num_threads: int,
        verbose: bool = False,
        lzani_exec: typing.Union[Path, str] = LZANI_EXEC
    ) -> subprocess.CompletedProcess:
    """Runs lz-ani to align genomic sequences.

    Args:
        input_fasta (Path):
            Path to the input FASTA file or directory with input FASTA files.
        infile_txt (Path or str):
            Path to the input text file listing FASTA files.
        infile_all2all (Path or str):
            Path to the input all2all file.
        outfile_aln (Path or str):
            Path to the output alignment file.
        kmer_count (int):
            Minimum number of shared k-mers between two sequences
            to start alignment.
        num_threads (int):
            Number of threads to use in lz-ani.
        verbose (Path or str):
            Whether to display verbose output.
        lzani_exec (Path or str):
            Path to the lz-ani executable.

    Returns:
        subprocess.CompletedProcess: Completed process information.

    """
    cmd = [
        f'{lzani_exec}',
        'all2all',
        '-t',
        f'{num_threads}',
        '--verbose', f'{int(verbose) + 1}',
        '-bs',
        '-cd', '128',
        '-reg', '80',
        '-mml', '5',
        '-mdl', '11',
        '-mlrim', '32',
        '-aw', '16',
        '-am', '6',
        '-ar', '2',

        '-out', f'{outfile_aln}',
        '-filter', f'{infile_all2all}', f'{kmer_count}'
    ]
    if input_fasta.is_dir():
        cmd.extend(['--in-file-names', f'{infile_txt}'])
    else:
        cmd.extend(['--one-file-name', f'{input_fasta}'])
    process = subprocess.run(
        cmd,  
        stdout=None if verbose else subprocess.DEVNULL, 
        stderr=None if verbose else subprocess.PIPE,
        text=True,
    )
    return process


def align(
        input_fasta: Path,
        outfile_aln: typing.Union[Path, str],
        kmer_size: int,
        kmer_count: int,
        num_threads: int,
        verbose: bool = False,
        kmerdb_exec: typing.Union[Path, str] = KMERDB_EXEC,
        lzani_exec: typing.Union[Path, str] = LZANI_EXEC
    ) -> None:
    """Prefilters FASTA sequences using kmer-db and aligns them using lz-ani.

    Args:
        input_fasta (Path):
            Path to the input FASTA file or directory with input FASTA files.
        outfile_aln (Path or str):
            Path to the output lz-ani alignment file.
        kmer_size (int):
            k-mer size used in kmer-db
        kmer_count (int):
            Minimum number of shared k-mers to report by kmer-db and required
            to start alignment in lz-ani.
        num_threads (int):
            Number of threads to use in kmer-db and lz-ani.
        verbose (bool):
            Whether to display verbose output.
        kmerdb_exec (Path or str):
            Path to the kmer-db executable.
        lzani_exec (Path or str):
            Path to the lz-ani executable.

    Returns:
        str: unique random identifier (token)

    """
    uuid_token = str(uuid.uuid4().hex)
    out_dir = Path(outfile_aln).parent
    txt_path = out_dir / f'{uuid_token}.input'
    db_path = out_dir / f'{uuid_token}.db'
    all2all_path = out_dir / f'{uuid_token}.all2all'

    # Run kmer-db build.
    p = validate_process(run_kmerdb_build(
        input_fasta=input_fasta, 
        outfile_txt=txt_path,
        outfile_db=db_path,
        kmer_size=kmer_size,
        num_threads=num_threads,
        verbose=verbose,
        kmerdb_exec=kmerdb_exec,
    ))
    # Run kmer-db all2all.
    p = validate_process(run_kmerdb_all2all(
        infile_db=db_path,
        outfile_all2all=all2all_path,
        kmer_count=kmer_count,
        num_threads=num_threads,
        verbose=verbose,
        kmerdb_exec=kmerdb_exec,
    ))
    # Run lz-ani.
    p = validate_process(run_lzani(
        input_fasta=input_fasta,
        infile_txt=txt_path,
        infile_all2all=all2all_path,
        outfile_aln=outfile_aln,
        kmer_count=kmer_count,
        num_threads=num_threads,
        verbose=verbose,
        lzani_exec=lzani_exec,
    ))
    # Remove temporary kmer-db related files.
    for file_to_remove in [txt_path, db_path, all2all_path]:
        file_to_remove.unlink()

    return uuid_token


def calcani(
        infile_aln: typing.Union[Path, str],
        outfile_ani: typing.Union[Path, str]
    ) -> None:
    """Calculates ANI values from lz-ani alignments and saves in csv file. 

    Args:
        infile_aln (Path or str): input lz-ani alignment file.
        outfile_ani (Path or str): output csv file.

    Returns:
        None

    """
    with open(infile_aln) as fh, open(outfile_ani, 'w') as oh:
        csv_out = csv.writer(oh)
        csv_out.writerow([
            'id1', 'id2', 'tANI', 'gANI1', 'cov1', 
            'gANI2', 'cov2', 'length_ratio'
        ])
        
        # Skip lines until the section '[input_sequences]' is found.
        for line in fh:
            if line.strip() == '[input_sequences]':
                break

        # Read genome sequence ids and the genome lengths.
        idx2id = {}
        lengths = {}
        for idx, line in enumerate(fh):
            if line.rstrip() == '[lz_similarities]':
                break
            cols = line.split()
            id = cols[1]
            length = int(cols[2])
            idx2id[idx] = id
            lengths[idx] = length

        # Loop over lz-ani alignment lines, calculate ANI, and save to file.
        for line in fh:
            cols = line.split()
            idx1 = int(cols[0])
            idx2 = int(cols[1])
            id1 = idx2id[idx1]
            id2 = idx2id[idx2]
            length1 = lengths[idx1]
            length2 = lengths[idx2]

            identities1 = int(cols[2])
            mismatches1 = int(cols[3])
            segments1 = int(cols[4])
            identities2 = int(cols[5])
            mismatches2 = int(cols[6])
            segments2 = int(cols[7])

            # Calculate tANI, gANI1, cov1, gANI2, cov2, length_ratio
            tani = (identities1 + identities2) / (length1 + length2)
            gani1 = identities1 / length1
            cov1 = (identities1 + mismatches1) / length1
            gani2 = identities2 / length2
            cov2 = (identities2 + mismatches2) / length2
            length_ratio = min(length1, length2) / max(length1, length2)

            csv_out.writerow([
                id1,
                id2,
                f'{tani:.5f}',
                f'{gani1:.5f}',
                f'{cov1:.5f}',
                f'{gani2:.5f}',
                f'{cov2:.5f}',
                f'{length_ratio:.5f}'
            ])


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
    args.input = Path(args.input)
    args.output = Path(args.output)
    if not args.input.exists():
        parser.error(f'Input does not exist: {args.input.resolve()}')
    if args.command == 'align':
        args.kmerdb_exec = validate_executable(args.kmerdb_exec)
        args.lzani_exec = validate_executable(args.lzani_exec)
        align(
            input_fasta=args.input,
            outfile_aln=args.output,
            kmer_size=args.kmer_size,
            kmer_count=args.kmer_count,
            num_threads=args.num_threads,
            verbose=args.verbose,
            kmerdb_exec=args.kmerdb_exec,
            lzani_exec=args.lzani_exec
        )
    elif args.command == 'calcani':
        calcani(infile_aln=args.input, outfile_ani=args.output)
    elif args.command == 'cluster':
        pass
    return args


if __name__ == '__main__':
    main()