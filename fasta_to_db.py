#!/home/eng/mablelab/miniconda3/envs/biopython/bin/python3
from Bio import SeqIO
from Bio import Blast
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
import re
import csv
import logging

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

Entrez.api_key = "a736c6d5906199057ec0dfec7a8b8d519708"

class InfoOnlyFilter(logging.Filter):
	def filter(self, record):
		return record.levelno <= logging.INFO

fh = logging.FileHandler('log_troubleshoot.txt', mode='w')
fh.setLevel(logging.WARNING)
formatter = logging.Formatter('%(levelname)s - %(asctime)s - %(message)s')
fh.setFormatter(formatter)
logger.addHandler(fh)

sh = logging.StreamHandler()
sh.setLevel(logging.DEBUG)
sh.addFilter(InfoOnlyFilter())
sh.setFormatter(formatter)
logger.addHandler(sh)

Entrez.email="thomasharrissnell@gmail.com"

# === INPUTS ===
fasta_file = "Tanzania_hostrefeseq_28-7-25_forTom_no_gaps.fasta"  # your input FASTA file
output_tsv = "taxonomy_re_run.tsv"          # output TSV file
default_taxonomy = "k__Animalia;p__Chordata;c__Unknown;o__Unknown;f__Unknown;g__Unknown;s__Unknown"

# === Optional: Mapping of sequence IDs to known taxonomies ===
# Format: {sequence_id: taxonomy_string}
# Fill this dict manually if you know the taxonomies
id_to_taxonomy = {}

def fetch_taxonomy(accession):
	handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="xml")
	record = Entrez.read(handle)
	handle.close()

	try:
		taxid = record[0]['GBSeq_taxonomy']
		organism = record[0]['GBSeq_organism']
		print (organism)
		print (taxid)
		org_parts = organism.split()
		out = taxid + '; ' + org_parts[1]
		return (out)
	except:
		logger.error("Could not extract taxonomy for {}\n".format(accession))

for record in SeqIO.parse(fasta_file, "fasta"):
	seq_id = record.id
	logger.info("Finding taxonomy for species {}...\n".format(seq_id))
	result_stream = NCBIWWW.qblast("blastn", "nt", record.seq, format_type='XML', hitlist_size = 1, megablast=True)
	blast_record = NCBIXML.read(result_stream)
	if blast_record.alignments:
		hit = blast_record.alignments[0]
		title = hit.title
		match = re.search(r'(?:gb|tpg|ref|emb)\|(\w+)\.', title)
		if match:
			accession = match.group(1)
			tax = fetch_taxonomy(accession)
			print (tax)
		else:
			logger.error("Could not extract accession number for {}\n".format(title))
	else:
		logger.error("No hits found in DB for {}\n".format(seq_id))
	try:
		tax_list=tax.split(';')
#		pos=[3,7,10,12,13,14]
		pos=[3,7,10,-2,-1,-0]
		tax_list= ['Animalia'] + [tax_list[i-1].strip() for i in pos]
		logger.info("Taxonomy in appropiate format: {}".format(tax_list))
		default_tax_list=default_taxonomy.split(';')
		logger.info("Taxonomy format is: {}\n".format(default_tax_list))
	except Exception as e:
		if tax is None:
			logger.error("Could not fetch tax for {}\n".format(seq_id))
			logger.info("Could not fetch tax for {}\n".format(seq_id))
			continue
		else:
			tax_list=tax.split(';')
			logger.info("Different organism taxonomy structure in EntrezDB\nOrganism: {}\n Tax: {}\n".format(seq_id, tax_list))
			logger.error("This organism doesnt have the same taxonomy format:\n{}\n{}\nError Type: {}\n".format(seq_id, tax_list, e))
	new_tax=[]
	for old, new in zip(default_tax_list, tax_list):
		prefix=old[0:2]
		new_name=prefix + new
		new_tax.append(new_name)
	new_tax=(';').join(new_tax)
	id_to_taxonomy[seq_id] = new_tax
	result_stream.close()

with open(output_tsv, "w", newline="") as tsvfile:
	writer = csv.writer(tsvfile, delimiter='\t')
	writer.writerow(["Feature ID", "taxonomy"])  # QIIME 2 header
	for seq_id, taxonomy in id_to_taxonomy.items():
		writer.writerow([seq_id, taxonomy])
logger.info(f"✅ Taxonomy file written to: {output_tsv}")
