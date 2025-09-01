#!/bin/bash
data="/shared2/Thomas_Harris_Snell/ONT_Amplicon/Tsetse_pool3-6_2Aug2023/20230802_1326_X1_FAX05124_bdbff985"
out_pass="/shared2/Thomas_Harris_Snell/ONT_Amplicon/amplicon_sorter_histogram/pass"
amplicon_sorter="/shared2/Thomas_Harris_Snell/amplicon_sorter.py"
cd "${data}/fastq_pass"
for dir in *; do
	cd "${dir}/Q15"
	mkdir -p "${out_pass}/${dir}"
	out_temp="${out_pass}/${dir}"
	result=""
	for txtfile in *.txt; do
		if [[ -f "$txtfile" ]]; then
			result=$(awk '{ print $1*100}' "$txtfile")
			echo "$result"
		fi
	done
	for fastq in *Q15.fastq; do
		if [[ -f "$fastq" ]]; then
			echo "$result"
			python3 "${amplicon_sorter}" -i "$fastq" -ra -maxr "$result" -np 15 -o "${out_temp}" -ho  #ho optional
		fi
	done
	cd ../..
done
