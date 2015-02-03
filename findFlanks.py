from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
import sys

def runBlastn(q_seq, s_seq, out_file):
	blastn_cline = NcbiblastnCommandline(query = q_seq, subject = s_seq, evalue=0.001, outfmt=5, out = out_file)
	stdout,stderr=blastn_cline()
	return
	
def parseBlastResults(result_handle):
	with open(result_handle,'r') as f:
		blast_records = NCBIXML.parse(f)
		for alignment in blast_record.alignments:
			print alignment
			for hsp in alignment.hsps:
				print hsp
	return
	
	
runBlastn(sys.argv[1], sys.argv[2], sys.argv[3])


parseBlastResults(sys.argv[3])
