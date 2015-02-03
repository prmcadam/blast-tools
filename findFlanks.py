from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
import sys

def runBlastn(q_seq, s_seq, out_file):
	blastn_cline = NcbiblastnCommandline(query = q_seq, db = s_seq, evalue=0.001, outfmt=5, out = out_file)
	stdout,stderr=blastn_cline()
	return
	
def parseBlastResults(result_handle):
	blast_record = NCBIXML.read(result_handle)
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			print hsp
	return
	
	
runBlastn(sys.argv[1], sys.argv[2], sys.argv[3])

parseBlastResults(sys.argv[3])