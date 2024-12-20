#!/usr/bin/env python3

from pathlib import Path
import shutil
import subprocess

import vclust
import pytest

VCLUST = Path('vclust.py')

DATA_DIR = Path('example')
TMP_DIR = DATA_DIR / 'tmp'
FASTA_DIR = DATA_DIR / 'fna'
FASTA_FILE = DATA_DIR / 'multifasta.fna'
FASTAGZ_FILE = DATA_DIR / 'multifasta.fna.gz' 
ANI_FILE = DATA_DIR / 'output'/ 'ani.tsv'
ALN_FILE = DATA_DIR / 'output' / 'ani.aln.tsv'
IDS_FILE = DATA_DIR / 'output' / 'ani.ids.tsv'
FLTR_FILE = DATA_DIR / 'output' / 'fltr.txt'
DATASET_DIR = DATA_DIR / 'datasets'
DATASET_FILES = [
    DATASET_DIR / 'refseq.fna',
    DATASET_DIR / 'genbank.fna',
    DATASET_DIR / 'other.fna',
]


@pytest.fixture
def test_dir():
    """Returns a testing data"""
    #print('setup')
    if TMP_DIR.exists():
        shutil.rmtree(TMP_DIR)
    TMP_DIR.mkdir(parents=True)
    yield TMP_DIR
    #print('teardown')
    shutil.rmtree(TMP_DIR)


@pytest.mark.parametrize('subcommand',[
    'deduplicate', 'prefilter', 'align', 'cluster', 'info',
])
def test_subcommands(subcommand):
    cmd = [
        f'{VCLUST.resolve()}',
        f'{subcommand}',
    ]
    p = subprocess.run(cmd,
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True)
    assert p.returncode == 0
    assert not p.stderr
    assert p.stdout


@pytest.mark.parametrize('input,params,error_msg',[
    (['missing_file1.fna', 'missing_file2.fna'], [], 'does not exist'),
    (DATASET_FILES, ['--add-prefixes', 'refseq', 'genbank'], 'error:'),
    (DATASET_FILES, ['--gzip-level', '0'], 'between 1 and 9'),
])
def test_parser_error_deduplicate(test_dir, input, params, error_msg):
    out_file = test_dir.joinpath('nr.fna')
    cmd = [
        f'{VCLUST.resolve()}',
        'deduplicate',
        '-i',
        *[str(f) for f in input],
        '-o',
        f'{out_file}',
        '-v',
        '0',
    ]
    cmd.extend(params)
    p = subprocess.run(cmd,
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True)
    assert p.returncode == 2
    assert error_msg in p.stderr
    assert not p.stdout


@pytest.mark.parametrize('input,params,error_msg',[
    (FASTA_DIR, ['--batch-size', '4'], 'error: --batch-size'),
    (FASTA_DIR, ['--min-ident', '95'], 'between 0 and 1'),
    (FASTA_DIR, ['--kmers-fraction', '10'], 'between 0 and 1'),
    (FASTA_DIR, ['--k', '2'], 'invalid choice'),
    (Path('missing_file.fna'), [], 'does not exist'),
])
def test_parser_error_prefilter(test_dir, input, params, error_msg):
    out_file = test_dir.joinpath('filter.txt')
    cmd = [
        f'{VCLUST.resolve()}',
        'prefilter',
        '-i',
        f'{input}',
        '-o',
        f'{out_file}',
    ]
    cmd.extend(params)
    p = subprocess.run(cmd,
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True)
    assert p.returncode == 2
    assert error_msg in p.stderr
    assert not p.stdout


@pytest.mark.parametrize('input,params,error_msg',[
    (FASTA_FILE, ['--out-tani', '40'], 'between 0 and 1'),
    (Path('missing_file.fna'), [], 'does not exist'),
])
def test_parser_error_align(test_dir, input, params, error_msg):
    out_file = test_dir.joinpath('ani.tsv')
    cmd = [
        f'{VCLUST.resolve()}',
        'align',
        '-i',
        f'{input}',
        '-o',
        f'{out_file}'
    ]
    cmd.extend(params)
    p = subprocess.run(cmd,
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True)
    assert p.returncode == 2
    assert error_msg in p.stderr
    assert not p.stdout


@pytest.mark.parametrize('params,error_msg',[
    (['--metric', 'tani'], 'error: tani threshold'),
    (['--metric', 'ani', '--ani', '95'], 'between 0 and 1'),
])
def test_parser_error_cluster(test_dir, params, error_msg):
    out_file = test_dir.joinpath('clusters.tsv')
    cmd = [
        f'{VCLUST.resolve()}',
        'cluster',
        '-i',
        f'{ANI_FILE}',
        '-o',
        f'{out_file}',
        '--ids',
        f'{IDS_FILE}'
    ]
    cmd.extend(params)
    p = subprocess.run(cmd,
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True)
    assert p.returncode == 2
    assert error_msg in p.stderr
    assert not p.stdout


def test_deduplicate_default(test_dir):
    out_file = test_dir.joinpath('nr.fna')
    out_duplicates_file = test_dir.joinpath('nr.fna.duplicates.txt')
    ref_ids = ['NC_002486.1', 'NC_005091.2', 'NC_010807.1', 'NC_025457.1', 
    'KJ473423.1', 'MN428048.1', 'MK937595.1', 'Mushuvirus', 'Mushuvirus_copy']
    cmd = [
        f'{VCLUST.resolve()}',
        'deduplicate',
        '-i',
        *[str(f) for f in DATASET_FILES],
        '-o',
        f'{out_file}',
        '-v',
        '0',
    ]
    p = subprocess.run(cmd,
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True)
    assert p.returncode == 0
    assert not p.stderr
    assert out_file.exists()
    assert out_duplicates_file.exists()
    assert out_file.stat().st_size
    assert out_duplicates_file.stat().st_size
    seq_ids = []
    with open(out_file) as fh:
        for line in fh:
            if line.startswith('>'):
                seq_ids.append(line.split()[0].lstrip('>'))
    assert len(seq_ids) == len(ref_ids)
    assert seq_ids == ref_ids


def test_deduplicate_default(test_dir):
    out_file = test_dir.joinpath('nr.fna')
    out_duplicates_file = test_dir.joinpath('nr.fna.duplicates.txt')
    ref_ids = ['NC_002486.1', 'NC_005091.2', 'NC_010807.1', 'NC_025457.1', 
               'MN428048.1', 'MK937595.1', 'Mushuvirus']
    ref_duplicates = {
        'Mushuvirus -Mushuvirus_copy',
        'NC_025457.1 -KJ473423.1',
        'NC_010807.1 -EU547803.1 -NC_010807.1_duplicate',
        'NC_005091.2 -AY357582.2 -AY357582.2_duplicate',
        'MN428048.1 +MN428048.1_revcomp',
        'NC_002486.1 -AB044554.1',
    }
    cmd = [
        f'{VCLUST.resolve()}',
        'deduplicate',
        '-i',
        *[str(f) for f in DATASET_FILES],
        '-o',
        f'{out_file}',
        '-v',
        '0',
    ]
    p = subprocess.run(cmd,
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True)
    assert p.returncode == 0
    assert not p.stderr
    assert out_file.exists()
    assert out_duplicates_file.exists()
    assert out_file.stat().st_size
    assert out_duplicates_file.stat().st_size
    seq_ids = []
    with open(out_file) as fh:
        for line in fh:
            if line.startswith('>'):
                seq_ids.append(line.split()[0].lstrip('>'))
    assert len(seq_ids) == len(ref_ids)
    assert seq_ids == ref_ids
    duplicates = set()
    with open(out_duplicates_file) as fh:
        for line in fh:
            line = line.strip()
            if line:
                duplicates.add(line)
    assert duplicates == ref_duplicates


@pytest.mark.parametrize('input,params',[
    (DATASET_FILES, ['--add-prefixes']),
    (DATASET_FILES, ['--add-prefixes', 'refseq|', 'genbank|', 'other|']),
])
def test_deduplicate_add_prefixes(test_dir, input, params):
    out_file = test_dir.joinpath('nr.fna')
    out_duplicates_file = test_dir.joinpath('nr.fna.duplicates.txt')
    cmd = [
        f'{VCLUST.resolve()}',
        'deduplicate',
        '-i',
        *[str(f) for f in DATASET_FILES],
        '-o',
        f'{out_file}',
        '-v',
        '0',
    ]
    cmd.extend(params)
    p = subprocess.run(cmd,
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True)
    assert p.returncode == 0
    assert not p.stderr
    assert out_file.exists()
    assert out_duplicates_file.exists()
    assert out_file.stat().st_size
    assert out_duplicates_file.stat().st_size
    prefixes = []
    with open(out_file) as fh:
        for line in fh:
            if line.startswith('>'):
                prefixes.append(line.split()[0].split('|')[0].lstrip('>'))
    assert len(prefixes) == 7
    assert set(prefixes) == {'refseq', 'genbank', 'other'}
 

@pytest.mark.parametrize('input,params,out_file',[
    (DATASET_FILES, ['--gzip-output'], 'nr.fna.gz'),
    (DATASET_FILES, ['--gzip-output'], 'nr.fna'),
])
def test_deduplicate_gzip(test_dir, input, params, out_file):
    out_file = test_dir.joinpath(out_file)
    ref_file = test_dir.joinpath('nr.fna.gz')
    ref_duplicates_file = test_dir.joinpath('nr.fna.gz.duplicates.txt')
    cmd = [
        f'{VCLUST.resolve()}',
        'deduplicate',
        '-i',
        *[str(f) for f in DATASET_FILES],
        '-o',
        f'{out_file}',
        '-v',
        '0',
    ]
    cmd.extend(params)
    p = subprocess.run(cmd,
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True)
    assert p.returncode == 0
    assert not p.stderr
    assert ref_file.exists()
    assert ref_file.stat().st_size
    assert ref_duplicates_file.exists()
    assert ref_duplicates_file.stat().st_size


def test_deduplicate_verbose(test_dir):
    out_file = test_dir.joinpath('nr.fna')
    out_duplicates_file = test_dir.joinpath('nr.fna.duplicates.txt')
    cmd = [
        f'{VCLUST.resolve()}',
        'deduplicate',
        '-i',
        *[str(f) for f in DATASET_FILES],
        '-o',
        f'{out_file}',
    ]
    p = subprocess.run(cmd,
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True)
    assert p.returncode == 0
    assert all(word in p.stderr for word in ['Running', 'Completed', 'INFO'])
    assert out_file.exists()
    assert out_duplicates_file.exists()
    assert out_file.stat().st_size
    assert out_duplicates_file.stat().st_size


@pytest.mark.parametrize('input,params',[
    (FASTA_DIR, []),
    (FASTA_FILE, []),
    (FASTA_FILE, ['--batch-size', '4']),
    (FASTAGZ_FILE, []),
])
def test_prefilter_default(test_dir, input, params):
    out_file = test_dir.joinpath('filter.txt')
    cmd = [
        f'{VCLUST.resolve()}',
        'prefilter',
        '-i',
        f'{input}',
        '-o',
        f'{out_file}',
        '-v',
        '0',
    ]
    cmd.extend(params)
    p = subprocess.run(cmd,
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True)
    assert p.returncode == 0
    assert not p.stderr
    assert out_file.exists()
    assert out_file.stat().st_size
    
    with open(out_file) as fh:
        vids = fh.readline().strip().rstrip(',').split(',')[1:]
        idx2vid = {i: vid.rstrip('.fna') for i, vid in enumerate(vids, start=1)}
        results = {}
        for line in fh:
            cols = line.rstrip().rstrip(',').split(',')
            vid1 = cols[0].rstrip('.fna')
            fields = cols[1:]
            for field in fields:
                value = field.split(':')
                idx = int(value[0])
                ani = float(value[1])
                vid2 = idx2vid[idx]
                results[(vid1, vid2)] = ani
                results[(vid2, vid1)] = ani
        assert results[('NC_010807.alt1', 'NC_010807')] == 0.99848
        assert results[('NC_010807.alt2', 'NC_010807.alt3')] == 0.992238
        assert results[('NC_025457', 'NC_025457.alt1')] == 0.990832
        assert results[('NC_010807.alt1', 'NC_010807.alt3')] == 0.996723
        assert results[('NC_025457.alt2', 'NC_025457.alt1')] == 0.94527
        assert results[('NC_002486', 'NC_002486.alt')] == 0.999979
        assert len(results) == 26


@pytest.mark.parametrize('input,params',[
    (FASTA_FILE, ['--kmers-fraction', '0.2']),
    (FASTA_FILE, ['--max-seqs', '2']),
    (FASTA_FILE, ['-k', '20']),
])
def test_prefilter_params(test_dir, input, params):
    out_file = test_dir.joinpath('filter.txt')
    cmd = [
        f'{VCLUST.resolve()}',
        'prefilter',
        '-i',
        f'{input}',
        '-o',
        f'{out_file}',
        '-v',
        '0',
    ]
    cmd.extend(params)
    p = subprocess.run(cmd,
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True)
    assert p.returncode == 0
    assert not p.stderr
    assert out_file.exists()
    assert out_file.stat().st_size


def test_prefilter_verbose(test_dir):
    out_file = test_dir.joinpath('filter.txt')
    cmd = [
        f'{VCLUST.resolve()}',
        'prefilter',
        '-i',
        f'{FASTA_FILE}',
        '-o',
        f'{out_file}',
    ]
    p = subprocess.run(cmd,
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True)
    assert p.returncode == 0
    assert all(word in p.stderr for word in ['Running', 'Completed', 'INFO'])
    assert out_file.exists()
    assert out_file.stat().st_size


@pytest.mark.parametrize('input,params',[
    (FASTA_DIR, []),
    (FASTA_FILE, []),
    (FASTAGZ_FILE, []),
])
def test_align_default(test_dir, input, params):
    out_file = test_dir.joinpath('ani.tsv')
    cmd = [
        f'{VCLUST.resolve()}',
        'align',
        '-i',
        f'{input}',
        '-o',
        f'{out_file}'
    ]
    p = subprocess.run(cmd)
    assert p.returncode == 0
    assert p.stderr == None
    assert out_file.exists()
    assert out_file.stat().st_size
    ref_pairs = {
        ('NC_010807', 'NC_010807.alt1'): 0.99753,
        ('NC_010807', 'NC_010807.alt2'): 0.98985,
        ('NC_010807', 'NC_010807.alt3'): 0.98384,
        ('NC_005091', 'NC_005091.alt1'): 0.97161,
        ('NC_005091', 'NC_005091.alt2'): 0.96707,
        ('NC_025457', 'NC_025457.alt1'): 0.80607,
        ('NC_025457', 'NC_025457.alt2'): 0.75921,
        ('NC_002486', 'NC_002486.alt'): 1.00000,
    }
    pairs = {}
    with open(out_file) as fh:
        next(fh)
        for line in fh:
            cols = line.split()
            id1 = cols[2].rstrip('.fna')
            id2 = cols[3].rstrip('.fna')
            tani = float(cols[4])
            pairs[(id1, id2)] = tani
    for ref_pair, ref_tani in ref_pairs.items():
        tani = pairs[ref_pair]
        assert abs(tani - ref_tani) < 0.007
    

@pytest.mark.parametrize('outfmt,ref_cols',[
    ('standard', vclust.ALIGN_OUTFMT['standard']),
    ('lite', vclust.ALIGN_OUTFMT['lite']),
    ('complete', vclust.ALIGN_OUTFMT['complete']),
])
def test_align_outfmt(test_dir, outfmt, ref_cols):
    out_file = test_dir.joinpath('ani.tsv')
    cmd = [
        f'{VCLUST.resolve()}',
        'align',
        '-i',
        f'{FASTA_FILE}',
        '-o',
        f'{out_file}',
        '--outfmt',
        f'{outfmt}'
    ]
    p = subprocess.run(cmd)
    with open(out_file) as fh:
        cols = fh.readline().split()
        assert cols == ref_cols


@pytest.mark.parametrize('input,params',[
    (FASTA_DIR, []),
    (FASTA_FILE, []),
])
def test_align_alignments(test_dir, input, params):
    out_file = test_dir.joinpath('ani.tsv')
    out_aln_file = test_dir.joinpath('ani.aln.tsv')
    cmd = [
        f'{VCLUST.resolve()}',
        'align',
        '-i',
        f'{input}',
        '-o',
        f'{out_file}',
        '--out-aln',
        f'{out_aln_file}',

    ]
    p = subprocess.run(cmd)
    assert p.returncode == 0
    assert p.stderr == None
    assert out_aln_file.exists()
    assert out_aln_file.stat().st_size
    with open(out_aln_file) as fh:
        header = fh.readline().split()
        assert len(header) == 10
        assert fh.readlines()


def test_align_verbose(test_dir):
    out_file = test_dir.joinpath('ani.tsv')
    cmd = [
        f'{VCLUST.resolve()}',
        'align',
        '-i',
        f'{FASTA_FILE}',
        '-o',
        f'{out_file}',
    ]
    p = subprocess.run(cmd,
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True)
    assert p.returncode == 0
    assert all(word in p.stderr for word in ['Running', 'Completed', 'INFO'])
    assert out_file.exists()


@pytest.mark.parametrize('input,params',[
    (FASTA_DIR, []),
    (FASTA_FILE, []),
    (FASTA_FILE, ['--batch-size', '4']),
])
def test_workflow_prefilter_align(test_dir, input, params):
    filter_file = test_dir.joinpath('filter.txt')
    ani_file = test_dir.joinpath('ani.tsv')
    cmd = [
        f'{VCLUST.resolve()}',
        'prefilter',
        '-i',
        f'{input}',
        '-o',
        f'{filter_file}'
    ]
    cmd.extend(params)
    p = subprocess.run(cmd)
    assert p.returncode == 0
    assert p.stderr == None
    assert filter_file.exists()
    assert filter_file.stat().st_size

    cmd = [
        f'{VCLUST.resolve()}',
        'align',
        '-i',
        f'{input}',
        '-o',
        f'{ani_file}',
        '--filter',
        f'{filter_file}'
    ]
    p = subprocess.run(cmd)
    assert p.returncode == 0
    assert p.stderr == None
    assert ani_file.exists()
    assert ani_file.stat().st_size


@pytest.mark.parametrize('algorithm',[
    'single',
    'complete',
    'uclust',
    'cd-hit',
    'set-cover',
])
def test_cluster_algorithm(test_dir, algorithm):
    out_file = test_dir / 'clusters.tsv'
    cmd = [
        f'{VCLUST.resolve()}',
        'cluster',
        '-i',
        f'{ANI_FILE}',
        '--ids',
        f'{IDS_FILE}',
        '-o',
        f'{out_file}',
        '--algorithm',
        f'{algorithm}',
        '--metric',
        'tani',
        '--tani',
        '0.95',
        '-v',
        '0',
    ]
    p = subprocess.run(cmd, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True)
    assert p.returncode == 0
    assert not p.stderr
    assert out_file.exists()
    assert out_file.stat().st_size


@pytest.mark.parametrize('filtering_measure',[
    'tani',
    'gani',
    'ani',
    'qcov',
    'rcov',
])
def test_cluster_filtering(test_dir, filtering_measure):
    out_file = test_dir / 'clusters.tsv'
    cmd = [
        f'{VCLUST.resolve()}',
        'cluster',
        '-i',
        f'{ANI_FILE}',
        '--ids',
        f'{IDS_FILE}',
        '-o',
        f'{out_file}',
        '--algorithm',
        'single',
        '--metric',
        'tani',
        '--tani',
        '0.95',
        f'--{filtering_measure}',
        '0.85',
    ]
    p = subprocess.run(cmd, 
        stdout=subprocess.DEVNULL, 
        stderr=subprocess.DEVNULL)
    assert p.returncode == 0
    assert p.stderr == None
    assert out_file.exists()
    assert out_file.stat().st_size


def test_cluster_verbose(test_dir):
    out_file = test_dir / 'clusters.tsv'
    cmd = [
        f'{VCLUST.resolve()}',
        'cluster',
        '-i',
        f'{ANI_FILE}',
        '--ids',
        f'{IDS_FILE}',
        '-o',
        f'{out_file}',
        '--algorithm',
        'single',
        '--metric',
        'tani',
        '--tani',
        '0.95',
    ]
    p = subprocess.run(cmd, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True)
    assert p.returncode == 0
    assert all(word in p.stderr for word in ['Running', 'Completed', 'INFO'])
    assert out_file.exists()
    assert out_file.stat().st_size


@pytest.mark.parametrize('params',[
    ([]),
    (['--leiden-resolution', '0.8', '--leiden-iterations', '3']),
    (['--leiden-resolution', '0.8', '--leiden-beta', '0.001']),
])
def test_cluster_algorithm_leiden(test_dir, params):
    out_file = test_dir / 'clusters.tsv'
    cmd = [
        f'{VCLUST.resolve()}',
        'cluster',
        '-i',
        f'{ANI_FILE}',
        '--ids',
        f'{IDS_FILE}',
        '-o',
        f'{out_file}',
        '--algorithm',
        'leiden',
        '--metric',
        'tani',
        '--tani',
        '0.95',
    ]
    cmd.extend(params)
    p = subprocess.run(cmd, 
        stdout=subprocess.DEVNULL, 
        stderr=subprocess.DEVNULL)
    assert p.returncode == 0
    assert p.stderr == None
    assert out_file.exists()
    assert out_file.stat().st_size