#!/usr/bin/env python3

from pathlib import Path
import shutil
import subprocess

import vclust
import pytest

VCLUST = Path('vclust.py')

DATA_DIR = Path('example')
FASTA_DIR = DATA_DIR / 'fna'
FASTA_FILE = DATA_DIR / 'multifasta.fna'
ANI_FILE = DATA_DIR / 'output'/ 'ani.tsv'
IDS_FILE = DATA_DIR / 'output' / 'ani.ids.tsv'
FLTR_FILE = DATA_DIR / 'output' / 'fltr.txt'

TMP_DIR = DATA_DIR / 'tmp'


@pytest.fixture
def test_dir():
    '''Returns a testing data'''
    #print('setup')
    if TMP_DIR.exists():
        shutil.rmtree(TMP_DIR)
    TMP_DIR.mkdir(parents=True)
    yield TMP_DIR
    #print('teardown')
    shutil.rmtree(TMP_DIR)


@pytest.mark.parametrize('input,params,error_msg',[
    (FASTA_DIR, ['--batch-size', '4'], 'error: --batch-size'),
    (FASTA_DIR, ['--min-ident', '95'], 'error: --min-ident'),
    (FASTA_DIR, ['--k', '2'], 'error: --k'),
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


@pytest.mark.parametrize('params,error_msg',[
    (['--out-tani', '40'], 'error: --out-tani must be between'),
])
def test_parser_error_align(test_dir, params, error_msg):
    out_file = test_dir.joinpath('ani.tsv')
    cmd = [
        f'{VCLUST.resolve()}',
        'align',
        '-i',
        f'{FASTA_FILE}',
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
    (['--metric', 'ani', '--ani', '95'], '--ani must be between 0 and 1'),
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


@pytest.mark.parametrize('input,params',[
    (FASTA_DIR, []),
    (FASTA_FILE, []),
    (FASTA_FILE, ['--batch-size', '4']),
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


@pytest.mark.parametrize('input,params',[
    (FASTA_DIR, []),
    (FASTA_FILE, []),
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



def test_cluster_default(test_dir):
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
        stdout=subprocess.DEVNULL, 
        stderr=subprocess.DEVNULL)
    assert p.returncode == 0
    assert p.stderr == None
    assert out_file.exists()
    assert out_file.stat().st_size
