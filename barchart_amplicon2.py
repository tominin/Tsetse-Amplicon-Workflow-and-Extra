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
colours = {'Elephant': '#1f77b4','Buffalo': '#ff7f0e','Human': '#2ca02c','Goat': '#9467bd','Spotted Hyaena': '#e377c2','Giraffe': '#7f7f7f','Warthog': '#bcbd22','Antelope': '#aec7e8','Hippopotomus': '#ffbb78', 'Cattle': '#8c564b','Mouse': '#c49c94', 'Gazelle' : '#17becf', 'Lion' : '#d62728'}
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
def pool_select(barcode):
	barcode_data = merged_data[merged_data['amplicon sorter ID'].str.contains(r'{}'.format(barcode))]
	return (barcode_data)

pool3_barcode01 = pool_select('barcode01')
pool3_barcode05 = pool_select('barcode05')

pool4n5_barcode02 = pool_select('barcode02')
pool4n5_barcode03 = pool_select('barcode03')
pool4n5_barcode06 = pool_select('barcode06')
pool4n5_barcode07 = pool_select('barcode07')

pool6_barcode04 = pool_select('barcode04')
pool6_barcode08 = pool_select('barcode08')

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

barcode01_readcount_all = get_species_read_counts(pool3_barcode01)
barcode01_readcounts = barcode01_readcount_all[0]
barcode01_bad = barcode01_readcount_all[1]
barcode01_species = barcode01_readcount_all[2]

barcode05_readcount_all = get_species_read_counts(pool3_barcode05)
barcode05_readcounts = barcode05_readcount_all[0]
barcode05_bad = barcode05_readcount_all[1]
barcode05_species = barcode05_readcount_all[2]

barcode02_readcount_all = get_species_read_counts(pool4n5_barcode02)
barcode02_readcounts = barcode02_readcount_all[0]
barcode02_bad = barcode02_readcount_all[1]
barcode02_species = barcode02_readcount_all[2]

barcode03_readcount_all = get_species_read_counts(pool4n5_barcode03)
barcode03_readcounts = barcode03_readcount_all[0]
barcode03_bad = barcode03_readcount_all[1]
barcode03_species = barcode03_readcount_all[2]

barcode06_readcount_all = get_species_read_counts(pool4n5_barcode06)
barcode06_readcounts = barcode06_readcount_all[0]
barcode06_bad = barcode06_readcount_all[1]
barcode06_species = barcode06_readcount_all[2]

barcode07_readcount_all = get_species_read_counts(pool4n5_barcode07)
barcode07_readcounts = barcode07_readcount_all[0]
barcode07_bad = barcode07_readcount_all[1]
barcode07_species = barcode07_readcount_all[2]

barcode04_readcount_all = get_species_read_counts(pool6_barcode04)
barcode04_readcounts = barcode04_readcount_all[0]
barcode04_bad = barcode04_readcount_all[1]
barcode04_species = barcode04_readcount_all[2]

barcode08_readcount_all = get_species_read_counts(pool6_barcode08)
barcode08_readcounts = barcode08_readcount_all[0]
barcode08_bad = barcode08_readcount_all[1]
barcode08_species = barcode08_readcount_all[2]

all_dicts = [barcode01_species, barcode02_species, barcode03_species, barcode04_species, barcode05_species, barcode06_species,barcode07_species ,barcode08_species]

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

species_dict={'Loxodonta africana' : 0, 'Syncerus caffer' : 0, 'Homo sapiens' : 0, 'Capra hircus' : 0, 'Crocuta crocuta' : 0, 'Giraffa giraffa' : 0, 'Phacochoerus africanus' : 0, 'Tragelaphus scriptus' : 0, 'Hippotragus equinus' : 0, 'Taurotragus derbianus' : 0, 'Panthera leo' : 0, 'Tragelaphus oryx' : 0, 'Gazellax subgutturosa' : 0, 'Bos indicus' : 0, 'Gazella subgutturosa' : 0}
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

barcode01_pcts = calc_pct(barcode01_readcounts, species_dict)
barcode02_pcts = calc_pct(barcode02_readcounts, species_dict)
barcode03_pcts = calc_pct(barcode03_readcounts, species_dict)
barcode04_pcts = calc_pct(barcode04_readcounts, species_dict)
barcode05_pcts = calc_pct(barcode05_readcounts, species_dict)
barcode06_pcts = calc_pct(barcode06_readcounts, species_dict)
barcode07_pcts = calc_pct(barcode07_readcounts, species_dict)
barcode08_pcts = calc_pct(barcode08_readcounts, species_dict)

def nice_names(in_data, normal_namez):
	out_name = []
	out_pct = []
	for name, pct in sorted(in_data.items()):
		x = normal_namez[name]
		out_name.append(x)
		out_pct.append(pct)
	out_names, out_vals = combine_duplicates(out_name, out_pct)
	return ([out_names, out_vals])

barcode01_readcounts = nice_names(barcode01_readcounts, normal_names)
barcode02_readcounts = nice_names(barcode02_readcounts, normal_names)
barcode03_readcounts = nice_names(barcode03_readcounts, normal_names)
barcode04_readcounts = nice_names(barcode04_readcounts, normal_names)
barcode05_readcounts = nice_names(barcode05_readcounts, normal_names)
barcode06_readcounts = nice_names(barcode06_readcounts, normal_names)
barcode07_readcounts = nice_names(barcode07_readcounts, normal_names)
barcode08_readcounts = nice_names(barcode08_readcounts, normal_names)

barcode01 = nice_names(barcode01_pcts, normal_names)
names01 = barcode01[0]
barcode01 = barcode01[1]
barcode02 = nice_names(barcode02_pcts, normal_names)
names02 = barcode02[0]
barcode02 = barcode02[1]
barcode03 = nice_names(barcode03_pcts, normal_names)
names03 = barcode03[0]
barcode03 = barcode03[1]
barcode04 = nice_names(barcode04_pcts, normal_names)
names04 = barcode04[0]
barcode04 = barcode04[1]
barcode05 = nice_names(barcode05_pcts, normal_names)
names05 = barcode05[0]
barcode05 = barcode05[1]
barcode06 = nice_names(barcode06_pcts, normal_names)
names06 = barcode06[0]
barcode06 = barcode06[1]
barcode07 = nice_names(barcode07_pcts, normal_names)
names07 = barcode07[0]
barcode07 = barcode07[1]
barcode08 = nice_names(barcode08_pcts, normal_names)
names08 = barcode08[0]
barcode08 = barcode08[1]
names = list(set(names01 + names02 + names03 + names04 + names05 + names06 + names07 + names08))


def combined(*barcodes, names):
	out = []
	for name in names:
		row = []
		for bc_names, bc_values in barcodes:
			mapping = dict(zip(bc_names, bc_values))
			row.append(mapping.get(name, 0))  # 0 if species not present
		out.append(row)
	return out


y = combined(
	(names01, barcode01),
	(names02, barcode02),
	(names03, barcode03),
	(names04, barcode04),
	(names05, barcode05),
	(names06, barcode06),
	(names07, barcode07),
	(names08, barcode08),
	names=names
)

SIZE_DEFAULT = 10

plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.size"] = SIZE_DEFAULT
plt.rcParams["axes.titlesize"] = SIZE_DEFAULT
plt.rcParams["axes.labelsize"] = SIZE_DEFAULT
plt.rcParams["xtick.labelsize"] = SIZE_DEFAULT
plt.rcParams["ytick.labelsize"] = SIZE_DEFAULT


new_y = [np.array(lst) for lst in y]
stack_data = np.array(new_y).T

#new_x = ['Barcode 1', 'Barcode 2', 'Barcode 3', 'Barcode 4', 'Barcode 5', 'Barcode 6', 'Barcode 7', 'Barcode 8']
new_x = ['Old', ' Old', 'Old ', 'Old', 'Revised', ' Revised', 'Revised ', 'Revised']
groups = [
	([0, 4], "Pool3"),                      # indices 0,4
	([1, 2, 5, 6], "Pools 4&5"),              # indices 1,2,5,6
	([3, 7], "Pool6")                       # indices 3,7
]

fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
for ax, (indices, title) in zip(axes, groups):
	print (ax)
	print (indices, title)
	bottom = np.zeros(len(indices))
	print (bottom)
	for i, name in enumerate(names):
		ax.bar([new_x[j] for j in indices],
		stack_data[indices, i],
		bottom=bottom,
		color=colours[name],
		label=name)
		bottom += stack_data[indices, i]
		ax.set_xlabel(title)
		ax.set_ylabel("Fraction")

# Only one legend for all
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, title="Species",
	bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig("./actual_amplicon_barplot_subplots.png", dpi=300, bbox_inches="tight")

