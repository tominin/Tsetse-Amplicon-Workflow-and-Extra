#!/bin/bash
vcf="/shared5/Thomas_Harris_Snell/Tsetse_variant_MQ20_SNPonly.vcf.gz"
pca_out="/shared5/Thomas_Harris_Snell/pca"
linkage_out="$pca_out/Tsetse_SNPs"
pruned_out="${linkage_out}_pruned"
plink --vcf "$vcf" --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out "$linkage_out"
plink --vcf "$vcf" --extract "${linkage_out}.prune.in" --set-missing-var-ids @:# --allow-extra-chr --make-bed --out "$pruned_out"
plink --bfile "$pruned_out" --set-missing-var-ids @:# --allow-extra-chr --pca --out "$pca_out/Plink_PCA"
