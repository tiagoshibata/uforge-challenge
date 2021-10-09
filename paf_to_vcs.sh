#!/bin/bash
for f in aligned_to_reference/*.paf ; do
    sort -k6,6 -k8,8n $f | paftools.js call -f S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa - > ${f/paf/vcs}
done
