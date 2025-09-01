#!/bin/bash
dir='/shared2/Thomas_Harris_Snell/ONT_Amplicon/amplicon_sorter_out/copy_pass'
db='/shared2/Thomas_Harris_Snell/AmpliconZilla/qiime_db/Tanzania_unaligned_BLCA.fasta'
out_dir='/shared2/Thomas_Harris_Snell/AmpliconZilla/local_blast_out/new_new'
cd ${dir}
for dir in *; do
	cd "${dir}/${dir}_Q15"
	for file in *_Q15_consensussequences.fasta; do
		out="${out_dir}/$(basename "$file" .fasta)_blast.tsv"
		blastn -query "$file" -db "$db" -out "$out" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" -num_threads 8 -max_target_seqs 5 -max_hsps 1
	done
	cd ../..
done
