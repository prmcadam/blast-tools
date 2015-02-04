from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
from optparse import OptionParser
from Bio.SeqRecord import SeqRecord

def opts():
	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	parser.add_option('-f', '--flanking', dest='flanks', action='store', help='fasta file of 2 flanking sequences')
	parser.add_option('-c', '--contigs', dest='contigs', action='store', help='fasta file of contigs')
	parser.add_option('-b', '--blast_output', dest='blast_output', action='store', help='name for blast output xml')
	parser.add_option('-p', '--regions', dest='regions', action='store', help='fasta file of regions of interest')
	return parser.parse_args()

def runBlastn(q_seq, s_seq, out_file):
	blastn_cline = NcbiblastnCommandline(query = q_seq, subject = s_seq, evalue=0.001, outfmt=5, out = out_file, max_target_seqs = 1)
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
	start, end = coordinates[0]-1, coordinates[1]
	with open(in_handle,'r') as f:
		records = SeqIO.parse(f, 'fasta')
		for record in records:
			if record.id == contig:
				new_seq = record.seq[start:end]
				new_id = record.id+'|'+str(start-1)+':'+str(end)
	new_record = SeqRecord(new_seq, id=new_id, description='')
	return(new_record)

def oneContig(contigs):
	if len(contigs) == 1:
		return True
	else:
		return False

def bestHit():
	return


###########################

if __name__ == '__main__':
	results_dict = {}
	regions_dict = {}

	(options,args) = opts()
	flanking=options.flanks.strip()
	contigs_file=options.contigs.strip()
	blast_xml=options.blast_output.strip()+'.1.xml'
        regions_xml=options.blast_output.strip()+'.2.xml'
	results=options.blast_output.strip()+'.results'
	potential_regions=options.regions.strip()
	
	runBlastn(flanking, contigs_file, blast_xml)
	with open(blast_xml, 'r') as f:
		blast_records = NCBIXML.parse(f)
		for blast_record in blast_records:
			results_dict = parseBlastResults(blast_record, results_dict)

	contigs, coordinates = calculateRegionOfInterest(results_dict)	

	if len(contigs)==1:
		with open(results+'.fa','w') as temp_handle:
			record=returnMiddleSequence(contigs_file, contigs.pop(), coordinates)
			SeqIO.write(record,temp_handle,'fasta')
	elif len(contigs) > 1:
		records = []
		with open(results+'.fa','w') as temp_handle:
			for contig in contigs:
				record=returnMiddleSequence(contigs_file, contig, coordinates)
				SeqIO.write(records,temp_handle,'fasta')

		runBlastn(results+'.fa', potential_regions, regions_xml)
		with open(regions_xml, 'r') as f:
			blast_records = NCBIXML.parse(f)
			for blast_record in blast_records:
				regions_dict = parseBlastResults(blast_record, regions_dict)


	for i in regions_dict:
		result_string = results+','+i+','+','.join([str(x) for x in regions_dict[i]])+'\n'
		with open(results+'.csv','w') as f:
			f.write(result_string)
