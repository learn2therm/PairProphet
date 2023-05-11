import re
import pandas as pd
import numpy as np

from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastpCommandline
from Bio.Blast import NCBIXML

def make_blast_df(data):

    original_cols = data.columns

    data.rename(columns = {original_cols[0]: 'query', original_cols[1]: 'subject'}, inplace = True)

    n_unique_1 = np.unique(data['query'])
    n_unique_2 = np.unique(data['subject'])

    qid_dict = dict(zip(n_unique_1, range(len(n_unique_1))))
    sid_dict = dict(zip(n_unique_2, range(len(n_unique_2))))

    data['qid'] = [qid_dict[i] for i in data.iloc[:,0]]
    data['sid'] = [sid_dict[i] for i in data.iloc[:,1]]

    data['pid'] = data.index

    m = [
    'local_gap_compressed_percent_id',
    'scaled_local_query_percent_id',
    'scaled_local_symmetric_percent_id',
    'query_align_len',
    'query_align_cov',
    'subject_align_len',
    'subject_align_cov',
    'bit_score'
    ]

    final = pd.DataFrame()

    db = open('db.fasta','w')

    for row in data.iterrows():

        db = open('db.fasta','w')
        db.write(f">{row[1]['sid']}\n{row[1]['subject']}\n")
        db.close()

        query = open('query.fasta','w')
        query.write(f">{row[1]['qid']}\n{row[1]['query']}\n")
        query.close()

        cmd = NcbimakeblastdbCommandline(dbtype='prot', input_file='db.fasta', parse_seqids=True)()

        NcbiblastpCommandline(query='query.fasta', db='db.fasta', outfmt=5, out='test_out', word_size=3,
                             gapopen = 11, gapextend = 1)()

        result_handle = open('test_out')

        blast_records = NCBIXML.parse(result_handle)

        test = BlastMet(blast_records)

        outdf = test.compute_metric('ids')[['query_id','subject_id']]

        for metric in m:
            outdf[metric] = test.compute_metric(metric)[metric]

        final = pd.concat([final, outdf], ignore_index = True)

    final = final.astype({"query_id": int, "subject_id": int})

    blast_df = data.merge(final, left_on = ['qid','sid'], right_on = ['query_id','subject_id'])

    return blast_df


class BlastMet:
    """Handles computation of metrics for each alignment in a blast record.

    The HSP with the largest average sequence coverage is used for local metrics.

    Parameters
    ----------
    blast_record : the record containing all hits for a query
    """
    def __init__(self, blast_record):
        recs = []
        for rec in blast_record:
            recs.append(rec)           # This is clunky
        self.record = recs[0]
        self.qid = self.record.query.split(' ')[0]

    def id_hsp_best_cov(self, alignment):
        """Determine HSP with the most average coverage of both sequences.

        Returns
        -------
        Index of HSP with max average seq coverage
        Max average coverage
        """
        scores = []

        for hsp in alignment.hsps:
            scores.append(
                ((hsp.query_end +1 - hsp.query_start)/self.record.query_length +
                 (hsp.sbjct_end +1 - hsp.sbjct_start)/alignment.length)/2)
        return np.argmax(scores), max(scores)

    @staticmethod
    def raw_gap_compressed_percent_id(n_matches, n_gaps, n_columns, n_compressed_gaps):
        """Percent matches in sequence, including but compressing gaps.

        Parameters
        ----------
        n_matches : int, number of matches in match columns
        n_gaps : number of gaps in match columns
        n_columns : total number of alignment match columns
        n_compressed_gaps : number of compressed gaps in match columns
        """
        return n_matches / (n_columns - n_gaps + n_compressed_gaps)

    def compute_metric(self,metric_name: str):
        """Compute the metric with specified name for each alignment"""
        if not hasattr(self, metric_name):
            raise ValueError(f"No metric found with name : {metric_name}")
        else:
            metric = getattr(self, metric_name)

        outputs = []
        for alignment in self.record.alignments:
            hsp_id, _ = self.id_hsp_best_cov(alignment)
            hsp = alignment.hsps[hsp_id]
            outputs.append((self.qid, alignment.hit_id.split('|')[-1], metric(alignment, hsp)))
        return pd.DataFrame(data=outputs, columns=['query_id', 'subject_id', metric_name])

    def ids(self, alignment, hsp):
        return 0

    def subject_align_cov(self, alignment, hsp):
        """Fraction of AA on query string taken up by alignment"""
        return (hsp.sbjct_end +1 - hsp.sbjct_start)/alignment.length

    def subject_align_len(self, alignment, hsp):
        """Length of AA on query string taken up by alignment"""
        return int(hsp.sbjct_end +1 - hsp.sbjct_start)

    def query_align_cov(self, alignment, hsp):
        """Fraction of AA on query string taken up by alignment"""
        return (hsp.query_end +1 - hsp.query_start)/self.record.query_length

    def query_align_len(self, alignment, hsp):
        """Length of AA on query string taken up by alignment"""
        return int(hsp.query_end +1 - hsp.query_start)

    def bit_score(self, alignment, hsp):
        return hsp.score

    def local_gap_compressed_percent_id(self, alignment, hsp):
        """Percent matches in match sequence, including but compressing gaps.

        The largest local HSP score is used
        """
        n_matches = hsp.identities
        n_gaps = hsp.gaps
        n_columns = len(hsp.query)
        n_compressed_gaps = len(re.findall('-+', hsp.query))+len(re.findall('-+', hsp.sbjct))
        return self.raw_gap_compressed_percent_id(n_matches, n_gaps, n_columns, n_compressed_gaps)

    def scaled_local_query_percent_id(self, alignment, hsp):
        """Percent matches in query sequence based on best HSP."""
        return hsp.identities/self.record.query_length

    def scaled_local_symmetric_percent_id(self, alignment, hsp):
        """Percent matches compared to average seq length of query and subject based on best HSP"""
        return 2*hsp.identities/(self.record.query_length + alignment.length)
