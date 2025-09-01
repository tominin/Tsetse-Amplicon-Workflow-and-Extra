#!/bin/bash
data="/shared2/Thomas_Harris_Snell/Tanzania2024_Amplicon/barcodes"
Q15="/shared2/Thomas_Harris_Snell/Tanzania2024_Amplicon/Q15"
cd ${data}"
for dir in *; do
	cd ${dir}
	mkdir -p Q15
	for file in *.fastq.gz; do
		base_name="${file%.fastq.gz}"
		gunzip -c "$file" | NanoFilt -q 15 | gzip  > ./Q15/"$base_name"_QC.fastq.gz
		echo "Processed $file"
	done
	cd Q15
	zcat *.fastq.gz > "$dir"_Q15.fastq
	rm *.gz
	rm *Q12.fastq
	cat "$dir"_Q15.fastq | grep -c '^@' > "$base_name"_read_num.txt
	cd ../..
done
