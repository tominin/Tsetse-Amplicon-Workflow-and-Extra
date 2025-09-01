#!/bin/bash
bed="/shared5/Thomas_Harris_Snell/admixture/Tsetse_variant_bedfile.bed"
anc="/shared3/Anubhab/Tsetse_flies/single_fly_Tsetse_4-1-2023/01.RawData/Glo_pallidipes.fna"
for K in {1..10}
do
	admixture --cv ${bed} $K | tee log${K}.out
done
