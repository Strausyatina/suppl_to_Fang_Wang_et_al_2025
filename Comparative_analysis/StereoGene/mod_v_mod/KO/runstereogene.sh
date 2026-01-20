#!/bin/bash

# Model pairs
pairs=(
"KO_H3K27ac.model,KO_H3K27me3.model"
"KO_H3K27me3.model,KO_H3K4me1.model"
"KO_H3K4me1.model,KO_H3K4me2.model"
"KO_H3K27ac.model,KO_H3K4me1.model"
"KO_H3K36me3.model,KO_H3K4me1.model"
"KO_H3K4me1.model,KO_H3K4me3.model"
)

# StereoGene fixed parameters
stereogenes=<Path to StereoGene>
profPath=Profiles
trackPath=exampledata
chrom=<Path to genome.chrom.sizes>
wSize=100K
report=report
plotType=pdf
writeDistr=DETAIL

# Parameter combinations
bins=(1000)
kernels=(1500)

for idx in ${!bins[@]}; do
    bin=${bins[$idx]}
    KernelSigma=${kernels[$idx]}

    echo "==========================================="
    echo " Running parameter set $((idx+1)): bin=$bin, KernelSigma=$KernelSigma"
    echo "==========================================="

    for pair in "${pairs[@]}"; do
        IFS=',' read -r m1 m2 <<< "$pair"
        base1=$(basename $m1 .model)
        base2=$(basename $m2 .model)
        resPath="results_Avg_${base1}_${base2}_bin${bin}_sigma${KernelSigma}"

        echo "Running StereoGene for $m1 and $m2 with bin=$bin and KernelSigma=$KernelSigma"

            $stereogenes \
            profPath=$profPath \
            trackPath=$trackPath \
            chrom=$chrom \
            resPath=$resPath \
            wSize=$wSize \
            bin=$bin \
            KernelSigma=${KernelSigma} \
            report=$report \
            plotType=$plotType \
            writeDistr=$writeDistr \
            $m1 $m2
    done
done

