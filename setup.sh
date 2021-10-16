#!/bin/bash
set -e

cd "$(dirname "$0")"
echo Unzipping draft_genomas_leveduras.zip...
unzip -d draft draft_genomas_leveduras.zip

echo Decompressing S288C_reference_genome_Current_Release.tgz
tar xf S288C_reference_genome_Current_Release.tgz
gunzip S288C_reference_genome_R64-3-1_20210421/*.gz

pypy3 -m venv venv
. venv/bin/activate
pip install -r requirements.txt
