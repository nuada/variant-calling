#!/usr/bin/python

import pyhgvs as hgvs
import pyhgvs.utils as hgvs_utils
from pygr.seqdb import SequenceFileDB
from itertools import imap
import gzip
import sys

# TODO test with ensembl transcripts!

__GENOME__ = None
__TRANSCRIPTS__ = None

def read_snpeff_transcripts(bin_file):
    # Load transcrips info as map: snpeff.id ->
    # 'chrom'
    # 'start' transcription start
    # 'end' transcription end
    # 'id' transcript name
    # 'strand'
    # 'cds_start' coding region start
    # 'cds_end' coding region end
    # 'gene_name'
    # 'exons' list of tuples (ex_start, ex_end)
    chrom = dict()
    genes = dict()
    transcripts = dict()
    for line in bin_file:
        row = line.strip().split('\t')
        if row[0] == 'CHROMOSOME':
            # id -> name
            chrom[int(row[1])] = row[5]
        elif row[0] == 'GENE':
            genes[int(row[1])] = {
                'gene_name': row[8],
                'chrom': int(row[2]),
                'strand': {'-1':'-', '1':'+'}[row[6]] }
        elif row[0] == 'TRANSCRIPT':
            if not int(row[1]) in transcripts:
                transcripts[int(row[1])] = dict()
            transcripts[int(row[1])].update({
                'name': row[5],
                'id': '.'.join(row[5].split('.')[:2]),
                'start': int(row[3]),
                'end': int(row[4])+1,
                'gene': int(row[2]) })
            # By default assume cds_{start,end} to be same as end
            if not 'cds_start' in transcripts[int(row[1])]:
                transcripts[int(row[1])]['cds_start'] = int(row[4])+1
            if not 'cds_end' in transcripts[int(row[1])]:
                transcripts[int(row[1])]['cds_end'] = int(row[4])+1
        elif row[0] == 'EXON':
            if not int(row[2]) in transcripts:
                transcripts[int(row[2])] = dict()
            if not 'exons' in transcripts[int(row[2])]:
                transcripts[int(row[2])]['exons'] = list()
            # (ex.no, ex.start, ex.end)
            transcripts[int(row[2])]['exons'].append((int(row[8]), int(row[3]), int(row[4])+1))
        elif row[0] == 'CDS':
            if not int(row[2]) in transcripts:
                transcripts[int(row[2])] = dict()
            transcripts[int(row[2])]['cds_start'] = min(transcripts[int(row[2])].get('cds_start', sys.maxint), int(row[3]))
            transcripts[int(row[2])]['cds_end'] = max(transcripts[int(row[2])].get('cds_end', 0), int(row[4])+1)

    db = dict()
    for t in transcripts.itervalues():
        t.update(genes[t['gene']])
        t['chrom'] = chrom[t['chrom']]
        t['exons'] = map(lambda x: x[1:], sorted(t['exons']))
        db[t['name']] = hgvs_utils.make_transcript(t)
    return db

def initialize(hg_fasta, snpeff_predictor_bin):
    ''' Load required databases:
         * human genome reference
         * refGene.txt'''
    global __GENOME__, __TRANSCRIPTS__

    __GENOME__ = SequenceFileDB(hg_fasta)
    with gzip.open(snpeff_predictor_bin) as f:
        __TRANSCRIPTS__ = read_snpeff_transcripts(f)

def get_transcript(name):
    ''' Get Transcript object for transcript name,
        search both transcript and transcript with version names '''
    return __TRANSCRIPTS__.get(name)

def mk_hgvs(chrom, zero_based_start, transcript_name, ref, alt, use_gene=True):
    ''' Return HGVS descrption of VCF record '''
    t = get_transcript(transcript_name)
    # Do not generate HGVS names without transcripts
    if not t:
        return ''
    return hgvs.format_hgvs_name(chrom, zero_based_start+1, ref, alt, __GENOME__, t, use_gene=use_gene)

if __name__ == '__main__':
    initialize('/resources/hg19/ucsc.hg19.fasta', '/resources/snpEff/3.6/hg19/snpEffectPredictor.bin')
    print mk_hgvs('chr14', 105177343, 'NM_001031714.3', 'C', 'T')
    sys.exit()
    import unittest

    class TestMkHGVS(unittest.TestCase):
        def setUp(self):
#            initialize('/resources/hg19/ucsc.hg19.fasta', '/resources/snpEff/3.6/hg19/snpEffectPredictor.bin')
            initialize('/resources/hg19/ucsc.hg19.fasta', '/resources/snpEff/3.6/GRCh37.75/snpEffectPredictor.bin')
            
        def test_simple_snps(self):
            self.assertEqual(mk_hgvs('chr14', 105177343, 'NM_001031714.3', 'C', 'T', use_gene=False), 'NM_001031714.3:c.2310C>T')

        def test_snps(self):
            snps = list()
            with open('mk_hgvs.test.snps') as f:
                for line in f:
                    row = line.strip().split('\t')
                    row[1] = int(row[1])
                    if '/' in row[5]:
                        row[5] = row[5].split('/')[1]
                    snps.append(row)

            for snp in snps:
                res = mk_hgvs(*snp[:5], use_gene=False)
                if res:
                    try:
                        self.assertEqual(res, ':'.join((snp[2], snp[5])))
                    except:
                        print(snp, res)

    unittest.main()
