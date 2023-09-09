#!/usr/bin/env python3

from pathlib import Path
import subprocess
import re
import unittest

class TestScript(unittest.TestCase):

    def setUp(self):
        self.dir = Path('test')
        self.input_fasta_dir = self.dir / 'fna'
        self.input_fasta_file = self.dir / 'multifasta.fna'
        self.aln_path = self.dir / 'alignment.txt'
        self.ani_path = self.dir / 'ani.csv'

    def test_align1(self):
        """Default parameters, check output file"""
        p = subprocess.run([
            f'./vclust.py',
            'align',
            '--threads', '1',
            f'{self.input_fasta_dir}',
            f'{self.aln_path}'
        ])
        self.assertEqual(p.returncode, 0)
        self.assertEqual(p.stderr, None)
        self.assertTrue(self.aln_path.exists())
        with open(self.aln_path) as fh:
            lines = fh.readlines()
        self.assertEqual(len(lines), 42)
        for i, line in enumerate(lines):
            if line.rstrip() == '[no_input_sequences]':
                break
        self.assertEqual(lines[i+1].rstrip(), '12')

    def test_align2(self):
        """Lower k-mer count, check output file"""
        p = subprocess.run([
            f'./vclust.py',
            'align',
            '--threads', '1',
            '--kmer_count', '0', 
            f'{self.input_fasta_dir}',
            f'{self.aln_path}'
        ])
        self.assertEqual(p.returncode, 0)
        self.assertEqual(p.stderr, None)
        self.assertTrue(self.aln_path.exists())
        with open(self.aln_path) as fh:
            lines = fh.readlines()
        self.assertGreater(len(lines), 42)


    def _test_calcani(self, input_fna):
        """Private"""
        p = subprocess.run([
            f'./vclust.py',
            'align',
            '--threads', '1',
            f'{input_fna}',
            f'{self.aln_path}'
        ])
        p = subprocess.run([
            f'./vclust.py',
            'calcani',
            f'{self.aln_path}',
            f'{self.ani_path}'
        ])
        self.assertEqual(p.returncode, 0)
        self.assertEqual(p.stderr, None)
        self.assertTrue(self.ani_path.exists())
        ref = {
            ('01.fna', '02.fna'): '1.00000,1.00000,1.00000,1.00000,1.00000,1.00000',
            ('02.fna', '03.fna'): '1.00000,1.00000,1.00000,1.00000,1.00000,1.00000',
            ('01.fna', '03.fna'): '1.00000,1.00000,1.00000,1.00000,1.00000,1.00000',
            ('04.fna', '05.fna'): '0.99100,0.99100,1.00000,0.99100,1.00000,1.00000',
            ('06.fna', '07.fna'): '0.98100,0.98100,1.00000,0.98100,1.00000,1.00000',
            ('08.fna', '09.fna'): '0.95108,0.95108,1.00000,0.95108,0.99998,1.00000',
            ('11.fna', '12.fna'): '0.66667,0.50000,0.50000,1.00000,1.00000,0.50000',
        }
        results = {}
        with open(self.ani_path) as fh:
            fh.readline()
            for line in fh:
                cols = line.rstrip().split(',')
                id1, id2 = sorted(cols[:2])
                id1 = f'{id1}.fna' if not id1.endswith('.fna') else id1
                id2 = f'{id2}.fna' if not id2.endswith('.fna') else id2
                results[(id1, id2)] = ",".join(cols[2:])

        for id_pair in ref:
            self.assertTrue(results[id_pair], ref[id_pair])

    def test_calcani_fastadir(self):
        """Fasta dir"""
        self._test_calcani(self.input_fasta_dir)

    def test_calcani_fastafile(self):
        """Single multifasta file"""
        self._test_calcani(self.input_fasta_file)

    def tearDown(self):
        for file_to_remove in [self.aln_path, self.ani_path]:
            if file_to_remove.exists():
                file_to_remove.unlink()


if __name__ == '__main__':
    unittest.main()
