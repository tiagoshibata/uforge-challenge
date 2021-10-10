#!/bin/bash
for f in aligned_to_reference/*.paf ; do
    sort -k6,6 -k8,8n $f | paftools.js call -f S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa - > ${f/paf/vcf}
    yes | bgzip ${f/paf/vcf}
    bcftools index ${f/paf/vcf.gz}
    bcftools consensus -f S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa ${f/paf/vcf.gz} -o ${f/paf/fa}
done
