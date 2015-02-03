from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
import sys



def runBlastn(q_seq, s_seq, out_file):
	blastn_cline = NcbiblastnCommandline(query = q_seq, subject = s_seq, evalue=0.001, outfmt=5, out = out_file)
	stdout,stderr=blastn_cline()
	return
	
def parseBlastResults(blast_record):
#	with open(result_handle,'r') as f:
#		blast_records = NCBIXML.parse(f)
	print blast_record
	for alignment in blast_record.alignments:
		print alignment
	return
	
	





###########################

if __name__ == '__main__':
	runBlastn(sys.argv[1], sys.argv[2], sys.argv[3])
	with open(sys.argv[3], 'r') as f:
		blast_records = NCBIXML.parse(f)
		for blast_record in blast_records:
		parseBlastResults(blast_record)
	
	