#!/home/tom/snap/miniconda3/envs/amplicon/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import re
import copy
from collections import defaultdict

#wanted to compare between the results of old vs revised primers
#look into pseudogene hits, save data in csv file to be formatted in excel for results

parser = argparse.ArgumentParser(description='Not sure yet', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--c', default=False, help='load in appropiate consensussequence tsv file')
args = parser.parse_args()
blast_tsv = args.c

col_names=('amplicon sorter ID', 'Feature ID', 'perident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'Description')
data=pd.read_csv(blast_tsv, sep='\t', names=col_names)
taxonomy=pd.read_csv('taxonomy_re_run.tsv', sep='\t', header=0, index_col=0)
expected=pd.read_csv('expected.tsv', sep='\t', header=0, index_col=0)
merged_data=data.merge(taxonomy, how='left', on='Feature ID')
colours = {'Elephant': '#1f77b4','Buffalo': '#ff7f0e','Human': '#2ca02c','Goat': '#9467bd','Spotted Hyaena': '#e377c2','Giraffe': '#7f7f7f','Warthog': '#bcbd22','Antelope': '#aec7e8','Hippopotomus': '#ffbb78', 'Cattle': '#8c564b','Mouse': '#c49c94', 'Panther' : '#d62728', 'Gazelle' : '#17becf'}
normal_names={'Loxodonta africana' : 'Elephant', 'Syncerus caffer' : 'Buffalo', 'Homo sapiens' : 'Human', 'Capra hircus' : 'Goat', 'Crocuta crocuta' : 'Spotted Hyaena', 'Giraffa giraffa' : 'Giraffe', 'Phacochoerus africanus' : 'Warthog', 'Tragelaphus scriptus' : 'Antelope', 'Hippotragus equinus' : 'Antelope', 'Taurotragus derbianus' : 'Antelope', 'Panthera leo' : 'Lion', 'Tragelaphus oryx' : 'Antelope', 'Gazellax subgutturosa' : 'Antelope', 'Bos indicus' : 'Cattle', 'Gazella subgutturosa' : 'Antelope'}
barcodes=['barcode01','barcode02','barcode03','barcode04','barcode05','barcode06','barcode07','barcode08']
barcode_to_pool = {'barcode01' : 'Pool3', 'barcode02' : 'Pools 4 & 5', 'barcode03' : 'Pools 4 & 5', 'barcode04' : 'Pool6', 'barcode05' : 'Pool3', 'barcode06' : 'Pools 4 & 5', 'barcode07' : 'Pools 4 & 5', 'barcode08' : 'Pool6'}

barcode_out_cols = ['pool', 'barcode', 'species identified', 'latin name', 'best hit', 'number of species above 90%', 'cytb:pseudogene counts', 'read number']

species_out_cols = ['Species', 'No. of Associated Sequences', 'Total Read Count', 'Pseudogene Read Count', 'No. of Ambigious Sequences', 'Total Ambigious Read Count', 'Pseudogene Ambigious Read Count', 'Average BLAST Score']

barcode_out = 'barcode_info.csv'
with open (barcode_out, 'w') as out_fh:
	out_fh.write(f"{','.join(barcode_out_cols)}\n")

species_out = 'species_info.csv'
with open (species_out, 'w') as out_fh:
	out_fh.write(f"{','.join(species_out_cols)}\n")

class Barcode:
	def __init__(self, pool, barcode, species_id, latin_name, best_hit_score, num_of_species, cytb_pseudo_ratio, read_number):
		self._pool = pool
		self._barcode = barcode
		self._species_id = species_id
		self._latin_name = latin_name
		self._best_hit_score = float(best_hit_score)
		self._num_of_species = int(num_of_species)
		self._cytb_pseudo_ratio = cytb_pseudo_ratio
		self._read_number = int(read_number)
	@property
	def pool(self):
		return self._pool
	@property
	def barcode(self):
		return self._barcode
	@property
	def species_id(self):
		return self._species_id
	@property
	def latin_name(self):
		return self._latin_name
	@property
	def best_hit_score(self):
		return self._best_hit_score
	@property
	def num_of_species(self):
		return self._num_of_species
	@property
	def cytb_pseudo_ratio(self):
		return self._cytb_pseudo_ratio
	@property
	def read_number(self):
		return self._read_number
	def write_out_csv(self):
		with open (barcode_out, 'a') as out_fh:
			out_fh.write(f"{self.pool},{self.barcode},{self.species_id},{self.latin_name},{self.best_hit_score},{self.num_of_species},{self.cytb_pseudo_ratio},{self.read_number}\n")
class Species:
	def __init__(self, species, num_of_barcodes_associated, read_count, read_count_positive_pseudo_hits, num_of_ambigious_barcodes, ambigious_read_count, read_count_ambigious_pseudo_hits, avg_BLAST):
		self._species = species
		self._num_of_barcodes_associated = num_of_barcodes_associated
		self._read_count = read_count
		self._read_count_positive_pseudo_hits = read_count_positive_pseudo_hits
		self._num_of_ambigious_barcodes = num_of_ambigious_barcodes
		self._ambigious_read_count = ambigious_read_count
		self._read_count_ambigious_pseudo_hits = read_count_ambigious_pseudo_hits
		self._avg_BLAST = avg_BLAST
	@property
	def species(self):
		return self._species
	@property
	def num_of_barcodes_associated(self):
		return self._num_of_barcodes_associated
	@property
	def read_count(self):
		return self._read_count
	@property
	def read_count_positive_pseudo_hits(self):
		return self._read_count_positive_pseudo_hits
	@property
	def num_of_ambigious_barcodes(self):
		return self._num_of_ambigious_barcodes
	@property
	def ambigious_read_count(self):
		return self._ambigious_read_count
	@property
	def read_count_ambigious_pseudo_hits(self):
		return self._read_count_ambigious_pseudo_hits
	@property
	def avg_BLAST(self):
		return self._avg_BLAST
	def write_out_csv(self):
		with open (species_out, 'a') as out_fh:
			out_fh.write(f"{self.species},{self.num_of_barcodes_associated},{self.read_count},{self.read_count_positive_pseudo_hits},{self.num_of_ambigious_barcodes},{self.ambigious_read_count},{self.read_count_ambigious_pseudo_hits},{self.avg_BLAST}\n")
def combine_duplicates(names, values):
	combined = defaultdict(float)
	for name, value in zip(names, values):
		combined[name] += value
	return list(combined.keys()), list(combined.values())

def file_string(input):
	temp = input.split('_Q15_')
	temp = temp[0].split('local_blast_out/')
	return (temp[1])
def read_num(amplicon):
	amplicon=amplicon.split('(')
	amplicon=amplicon[-1]
	return (int(amplicon[:-1]))
def add_to_dict(d, key, value):
	d[key] = d.get(key, 0) + value
def clean_species(taxon_list):
	taxon_list_ind=taxon_list[0]
	species=taxon_list_ind.split(';')
	species=species[-2:]
	new_species=[]
	for item in species:
		new_item=item[2:]
		new_species.append(new_item)
	species_fin=' '.join(new_species)
	return (species_fin)
def convert_name(latin):
	new = normal_names[latin]
	return (new)
def pool_select(pool_num):
	pool_bar={3 : 'barcode01,barcode05', 4 : 'barcode02,barcode03,barcode06,barcode07', 5 : 'barcode02,barcode03,barcode06,barcode07', 6 : 'barcode04,barcode08'}
	pool_bar = pool_bar[pool_num].split(',')
	frames = []
	for barcode in pool_bar:
		barcode_data = merged_data[merged_data['amplicon sorter ID'].str.contains(r'{}'.format(barcode))]
		frames.append(barcode_data)
	df = pd.concat(frames)
	return (df)

pool3_df = pool_select(3)
pool4n5_df = pool_select(4)
pool6_df = pool_select(6)

def pseudo_count(id_list):
	tot = len(id_list)
	pseudo_count = len(re.findall(r'pseudogene', ','.join(id_list)))
	norm_count = tot - pseudo_count
	out = ('{}:{}'.format(norm_count, pseudo_count))
	return (out)

def get_species_read_counts(pool_df):
	species_metrics = defaultdict(lambda: {
		'read_count': 0,
		'barcode_count_set': set(),
		'pseudogene_reads': 0,
		'bad_data_reads': 0,
		'bad_data_pseudogene_reads': 0,
		'bad_data_barcode_set': set(),
		'blast_scores' : []
	})

	amplicon_id_list=pool_df['amplicon sorter ID'].tolist()
	uniq_amp_id=list(dict.fromkeys(amplicon_id_list))

	out_data={}
	bad_data={}
	print (uniq_amp_id)
	for amp in uniq_amp_id:
		read_count = read_num(amp)
		save_read_count = read_count
		data_sub=pool_df[(pool_df['amplicon sorter ID'] == amp) & (pool_df['perident'] >= 90)]
		temp_df = pool_df[pool_df['amplicon sorter ID'] == amp]
		species_scores = defaultdict(list)
		for _, row in data_sub.iterrows():
			sp = convert_name(clean_species([row['Taxon']]))
			species_scores[sp].append(row['perident'])
		top5_scores = temp_df['perident'].tolist()
		save_hit = temp_df['perident'].max()
		temp = re.search(r'barcode0\d', amp)
		save_pool = barcode_to_pool[temp.group(0)]
		save_barcode = (temp.group(0))
		taxon_list=data_sub['Taxon'].tolist()
		id_list=list(set(data_sub['Description'].tolist()))
		hits = (len(list(set(taxon_list))))
		save_hits = hits

		if hits == 1:
			species_fin=clean_species(taxon_list)
			save_species = convert_name(species_fin)
			save_latin = species_fin
			top_scores = sorted(species_scores[save_species], reverse=True)[:5]
			if not top_scores:
				top_scores = sorted(temp_df['perident'].tolist(), reverse=True)[:5]
			add_to_dict(out_data, species_fin, read_count)
			species_metrics[save_species]['read_count'] += read_count
			species_metrics[save_species]['barcode_count_set'].add(amp)
			species_metrics[save_species]['blast_scores'].extend(top_scores)
			x = re.search('pseudogene', ','.join(id_list).lower())
			if x:
				save_norm_count = pseudo_count(id_list)
				species_metrics[save_species]['pseudogene_reads'] += read_count
			else:
				save_norm_count = pseudo_count(id_list)
		elif hits > 1:
			save_species = 'multiple species'
			save_latin = 'N/A'
			species_fin=clean_species(taxon_list)
			add_to_dict(bad_data, species_fin, read_count)
			for sp in list(set(taxon_list)):
				sp = convert_name(clean_species([sp]))
				top_scores = sorted(species_scores[sp], reverse=True)[:5]
				if not top_scores:
					top_scores = sorted(temp_df['perident'].tolist(), reverse=True)[:5]
				species_metrics[sp]['bad_data_reads'] += read_count
				species_metrics[sp]['bad_data_barcode_set'].add(amp)
				species_metrics[sp]['blast_scores'].extend(top_scores)
				x = re.search('pseudogene', ','.join(id_list).lower())
				if x:
					species_metrics[sp]['bad_data_pseudogene_reads'] += read_count
			if 'pseudogene' in ','.join(id_list).lower():
				save_norm_count = pseudo_count(id_list)
			else:
				save_norm_count = pseudo_count(id_list)
		else:
			save_species = 'no species'
			save_latin = 'N/A'
			save_norm_count = pseudo_count(id_list)
			continue
		collective = Barcode(save_pool, save_barcode, save_species, save_latin, save_hit, save_hits, save_norm_count, save_read_count)
		collective.write_out_csv()
	final_metrics = {}
	for sp, data in species_metrics.items():
		final_metrics[sp] = {
			'read_count': data['read_count'],
			'barcode_count': len(data['barcode_count_set']),
			'pseudogene_reads': data['pseudogene_reads'],
			'bad_data_reads': data['bad_data_reads'],
			'bad_data_pseudogene_reads': data['bad_data_pseudogene_reads'],
			'bad_data_barcode_count': len(data['bad_data_barcode_set']),
			'avg_top5_blast_score': np.mean(data['blast_scores'][:5]) if data['blast_scores'] else 0
		}
	return (out_data, bad_data, final_metrics)

pool3_readcount_all = get_species_read_counts(pool3_df)
pool3_readcounts = pool3_readcount_all[0]
pool3_readcounts_bad = pool3_readcount_all[1]
pool3_species = pool3_readcount_all[2]
pool4n5_readcounts_all = get_species_read_counts(pool4n5_df)
pool4n5_readcounts = pool4n5_readcounts_all[0]
pool4n5_readcounts_bad = pool4n5_readcounts_all[1]
pool4n5_species = pool4n5_readcounts_all[2]
pool6_readcounts_all = get_species_read_counts(pool6_df)
pool6_readcounts = pool6_readcounts_all[0]
pool6_readcounts_bad = pool6_readcounts_all[1]
pool6_species = pool6_readcounts_all[2]

all_dicts = [pool3_species, pool4n5_species, pool6_species]

# Setup output container
combined = defaultdict(lambda: {
	'read_count': 0,
	'barcode_count': 0,
	'pseudogene_reads': 0,
	'bad_data_reads': 0,
	'bad_data_pseudogene_reads': 0,
	'bad_data_barcode_count': 0,
	'blast_score_sum': 0.0,
	'blast_score_weight': 0
})

for d in all_dicts:
	for species, metrics in d.items():
		combined[species]['read_count'] += metrics['read_count']
		combined[species]['barcode_count'] += metrics['barcode_count']
		combined[species]['pseudogene_reads'] += metrics['pseudogene_reads']
		combined[species]['bad_data_reads'] += metrics['bad_data_reads']
		combined[species]['bad_data_pseudogene_reads'] += metrics['bad_data_pseudogene_reads']
		combined[species]['bad_data_barcode_count'] += metrics['bad_data_barcode_count']
		# Weight by barcode count (or any other weight you prefer)
		weight = metrics['barcode_count']
		combined[species]['blast_score_sum'] += metrics['avg_top5_blast_score'] * weight
		combined[species]['blast_score_weight'] += weight

# Final pass to compute weighted average for blast scores
for species in combined:
	total_weight = combined[species]['blast_score_weight']
	if total_weight > 0:
	        avg_score = combined[species]['blast_score_sum'] / total_weight
	else:
		avg_score = 0.0
	combined[species]['avg_top5_blast_score'] = round(avg_score, 4)

	del combined[species]['blast_score_sum']
	del combined[species]['blast_score_weight']

# Convert back to regular dict if needed
final_combined_dict = dict(combined)
for sp, metrics in final_combined_dict.items():
	collective = Species(sp, metrics['barcode_count'], metrics['read_count'], metrics['pseudogene_reads'], metrics['bad_data_barcode_count'], metrics['bad_data_reads'], metrics['bad_data_pseudogene_reads'], metrics['avg_top5_blast_score'])
	collective.write_out_csv()

species_dict={'Loxodonta africana' : 0, 'Syncerus caffer' : 0, 'Homo sapiens' : 0, 'Capra hircus' : 0, 'Crocuta crocuta' : 0, 'Giraffa giraffa' : 0, 'Phacochoerus africanus' : 0, 'Tragelaphus scriptus' : 0, 'Hippotragus equinus' : 0, 'Taurotragus derbianus' : 0}
names2=['Warthog','Buffalo','Spotted Hyaena','Giraffe','Elephant','Cattle','Antelope','Human','Goat','Mouse']
nice_dict = {normal_names[i]: 0 for i in species_dict if species_dict[i] == 0}
nice_dict2 = {name: 0 for name in names2}

def calc_pct(dict, template):
	total = np.sum(list(dict.values()))
	result = copy.deepcopy(template)
	for species, reads in dict.items():
		absolute = ((1/total) * reads)
		result[species] = float(absolute)
	return (result)

pool3_pcts = calc_pct(pool3_readcounts, species_dict)
pool4n5_pcts = calc_pct(pool4n5_readcounts, species_dict)
pool6_pcts = calc_pct(pool6_readcounts, species_dict)

def nice_names(in_data, normal_namez):
	out_name = []
	out_pct = []
	for name, pct in sorted(in_data.items()):
		x = normal_namez[name]
		out_name.append(x)
		out_pct.append(pct)
	out_names, out_vals = combine_duplicates(out_name, out_pct)
	return ([out_names, out_vals])

pool3_readcounts = nice_names(pool3_readcounts, normal_names)
pool4n5_readcounts = nice_names(pool4n5_readcounts, normal_names)
pool6_readcounts = nice_names(pool6_readcounts, normal_names)

pool3_readcounts_bad = nice_names(pool3_readcounts_bad, normal_names)
pool4n5_readcounts_bad = nice_names(pool4n5_readcounts_bad, normal_names)
pool6_readcounts_bad = nice_names(pool6_readcounts_bad, normal_names)



pool4n5 = nice_names(pool4n5_pcts, normal_names)
names = pool4n5[0]
pool4n5 = pool4n5[1]
pool3 = nice_names(pool3_pcts, normal_names)
pool3 = pool3[1]
pool6 = nice_names(pool6_pcts, normal_names)
pool6 = pool6[1]

def combined(x1,x2,x3,names):
	out = [[] for _ in range(len(names))]
	for index in range(len(names)):
		out[index].append(x1[index])
		out[index].append(x2[index])
		out[index].append(x3[index])
	return (out)

def expected_data(num):
	out = {}
	if num == 4 or num ==5:
		pool = expected.loc[[4, 5]]
	else:
		pool = expected.loc[[num]]
	item = pool['Dominant_host'].tolist()
	item = ','.join(item)
	list = item.split(',')
	for val in list:
		name, count = val.rstrip(')').split('(')
		add_to_dict(out, name, int(count))
	return (out)

temp2 = expected_data(4)
temp = expected_data(3)
temp3 = expected_data(6)
expected1 = calc_pct(temp, nice_dict2)
expected2 = calc_pct(temp2, nice_dict2)
expected3 = calc_pct(temp3, nice_dict2)

y = combined(pool3, pool4n5, pool6, names)

new_y = []
new_x = ['Pool3', 'Pools 4 & 5', 'Pool6']

for list in y:
	new_list = np.array(list)
	new_y.append(new_list)
stack_data = np.array(new_y).T
bottom = np.zeros(len(new_x))
for i, name in enumerate(names):
	plt.bar(new_x, stack_data[:, i], bottom=bottom, color=colours[name], label=name)
	bottom += stack_data[:, i]
SIZE_DEFAULT = 10

plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.size"] = SIZE_DEFAULT
plt.rcParams["axes.titlesize"] = SIZE_DEFAULT
plt.rcParams["axes.labelsize"] = SIZE_DEFAULT
plt.rcParams["xtick.labelsize"] = SIZE_DEFAULT
plt.rcParams["ytick.labelsize"] = SIZE_DEFAULT
plt.xlabel("Pool")
plt.ylabel("Fraction")
plt.title("Amplicon Species Composition")
plt.legend(title="Species", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("./actual_amplicon_barplot.png", dpi=300, bbox_inches="tight")


names2=['Warthog','Buffalo','Spotted Hyaena','Giraffe','Elephant','Cattle','Antelope','Human','Goat','Mouse']

expected_combined = combined(
	[expected1[name] for name in names2],
	[expected2[name] for name in names2],
	[expected3[name] for name in names2],
	names2
)
expected_stack_data = np.array(expected_combined).T

# Plot expected stacked bar chart
plt.figure()
bottom = np.zeros(len(new_x))
for i, name in enumerate(names2):
	plt.bar(new_x, expected_stack_data[:, i], bottom=bottom, color=colours[name], label=name)
	bottom += expected_stack_data[:, i]
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.size"] = SIZE_DEFAULT
plt.rcParams["axes.titlesize"] = SIZE_DEFAULT
plt.rcParams["axes.labelsize"] = SIZE_DEFAULT
plt.rcParams["xtick.labelsize"] = SIZE_DEFAULT
plt.rcParams["ytick.labelsize"] = SIZE_DEFAULT
plt.xlabel("Pool")
plt.ylabel("Fraction")
plt.title("Expected Species Composition")
plt.legend(title="Species", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("./expected_barplot.png", dpi=300, bbox_inches="tight")
