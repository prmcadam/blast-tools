from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
import sys



def runBlastn(q_seq, s_seq, out_file):
	blastn_cline = NcbiblastnCommandline(query = q_seq, subject = s_seq, evalue=0.001, outfmt=5, out = out_file)
	stdout,stderr=blastn_cline()
	return

def parseBlastResults(blast_record, results_dict):
	for alignment in blast_record.alignments:
		results_dict[str(blast_record.query)] = []
		for hsp in alignment.hsps:
			results_dict[blast_record.query].append(str(alignment.hit_def))
			results_dict[blast_record.query].append(hsp.sbjct_start)
                        results_dict[blast_record.query].append(hsp.sbjct_end)
	return results_dict

def calculateRegionOfInterest(results_dict):
	coordinates = []
	contig = set([])
	for x in results_dict:
		contig.add(results_dict[x][0])
		coordinates.append(results_dict[x][1])
                coordinates.append(results_dict[x][2])
	coordinates.sort()
	coordinates.remove(max(coordinates))
        coordinates.remove(min(coordinates))
	return contig, coordinates

def returnMiddleSequence(in_handle, contig, coordinates):
	start, end = coordinates[0], coordinates[1]
	with open(in_handle,'r') as f:
		records = SeqIO.parse(f, 'fasta')
		for record in records:
			if record.id == contig:
				print record.seq[start:end]



def oneContig(contigs):
	if len(contigs) == 1:
		return True
	else:
		return False

###########################

if __name__ == '__main__':
	results_dict = {}
	runBlastn(sys.argv[1], sys.argv[2], sys.argv[3])
	with open(sys.argv[3], 'r') as f:
		blast_records = NCBIXML.parse(f)
		for blast_record in blast_records:
			results_dict = parseBlastResults(blast_record, results_dict)

	contigs, coordinates = calculateRegionOfInterest(results_dict)	

	if oneContig(contigs)==True:
		returnMiddleSequence(sys.argv[2], contigs.pop(), coordinates)
