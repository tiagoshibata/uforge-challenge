Prerequisites: pypy3, minimap2, muscle

* On first run, run `. ./setup.sh` to decompress genomes and create a pypy3 virtual environment
* Run `. ./venv/bin/activate` to activate the virtual environment
* Run `./finisher.py` to assemble drafts to .sam/.paf files in `./aligned_to_reference`
* Run `./paf_to_vcf.sh` to make a consensus among overlapping regions and create .vcf and a FASTA file for each variation
* Run `./large_variants.py` to print large variations between each country and the reference genome, and generate images in the images/ directory
