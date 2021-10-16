#!/bin/bash
cd "$(dirname "$0")"
echo Unzipping draft_genomas_leveduras.zip...
unzip -d aligned_to_reference draft_genomas_leveduras.zip

echo Decompressing S288C_reference_genome_Current_Release.tgz
tar xf S288C_reference_genome_Current_Release.tgz
