Prerequisites: pypy3, minimap2, muscle

* On first run, run `. ./setup.sh` to decompress genomes and create a pypy3 virtual environment
* Run `. ./venv/bin/activate` to activate the virtual environment
* Run `./finisher.py` to assemble drafts to .sam/.paf files in `./aligned_to_reference`
* Run `./paf_to_vcf.sh` to make a consensus among overlapping regions and create .vcf and a FASTA file for each variation
* Run `./large_variants.py` to print large variations between each country and the reference genome, and generate images in the images/ directory

Possible phylogenetic tree:

```
Chromosome ref|NC_001134|, position in reference 29638, indel of length 5964 removed in relation to reference - ['Africa_do_Sul', 'Alemanha', 'Australia', 'Brasil', 'Colombia', 'EUA', 'Espanha', 'Holanda', 'Inglaterra', 'Irlanda', 'Nomina', 'Republica_Tcheca', 'Senegal']
Chromosome ref|NC_001134|, position in reference 221031, indel of length 5921 removed in relation to reference - ['Africa_do_Sul', 'Alemanha', 'Argentina', 'Australia', 'Belgica', 'Brasil', 'Colombia', 'EUA', 'Holanda', 'Inglaterra', 'Irlanda', 'Japao', 'Nomina', 'Republica_Tcheca', 'Senegal']
Chromosome ref|NC_001136|, position in reference 871815, indel of length 5964 removed in relation to reference - ['Alemanha', 'Argentina', 'Belgica', 'Brasil', 'EUA', 'Inglaterra', 'Irlanda', 'Japao', 'Nomina', 'Republica_Tcheca', 'Romenia']
Chromosome ref|NC_001136|, position in reference 878297, indel of length 5923 removed in relation to reference - ['Alemanha', 'Argentina', 'Belgica', 'Brasil', 'EUA', 'Inglaterra', 'Irlanda', 'Japao', 'Nomina', 'Republica_Tcheca', 'Romenia']
Chromosome ref|NC_001136|, position in reference 1095760, indel of length 5931 removed in relation to reference - ['Africa_do_Sul', 'Argentina', 'Belgica', 'Brasil', 'Colombia', 'Estonia', 'Inglaterra', 'Irlanda', 'Japao', 'Nomina', 'Romenia']
Chromosome ref|NC_001137|, position in reference 492433, indel of length 6077 removed in relation to reference - ['Australia', 'Espanha', 'Nomina']
Chromosome ref|NC_001142|, position in reference 472457, indel of length 11515 removed in relation to reference - ['Colombia', 'Inglaterra', 'Irlanda', 'Mexico', 'Nomina', 'Romenia', 'Senegal']
Chromosome ref|NC_001144|, position in reference 593139, indel of length 5923 removed in relation to reference - ['Alemanha', 'Australia', 'Brasil', 'EUA', 'Holanda', 'Inglaterra', 'Irlanda', 'Japao', 'Nomina', 'Republica_Tcheca']
Chromosome ref|NC_001144|, position in reference 650818, indel of length 5925 removed in relation to reference - ['Africa_do_Sul', 'Argentina', 'Australia', 'Belgica', 'Colombia', 'EUA', 'Espanha', 'Holanda', 'Inglaterra', 'Irlanda', 'Japao', 'Nomina', 'Republica_Tcheca', 'Romenia', 'Senegal']
Chromosome ref|NC_001144|, position in reference 941186, indel of length 5964 removed in relation to reference - ['Africa_do_Sul', 'Alemanha', 'Argentina', 'Australia', 'Belgica', 'Brasil', 'Colombia', 'EUA', 'Espanha', 'Estonia', 'Holanda', 'Inglaterra', 'Irlanda', 'Japao', 'Mexico', 'Nomina', 'Republica_Tcheca', 'Romenia', 'Russia', 'Senegal']
Chromosome ref|NC_001147|, position in reference 594814, indel of length 5919 removed in relation to reference - ['Alemanha', 'Argentina', 'Belgica', 'Brasil', 'Colombia', 'EUA', 'Inglaterra', 'Irlanda', 'Mexico', 'Nomina', 'Republica_Tcheca', 'Romenia']
Chromosome ref|NC_001148|, position in reference 804639, indel of length 5923 removed in relation to reference - ['Africa_do_Sul', 'Alemanha', 'Australia', 'Brasil', 'Colombia', 'EUA', 'Espanha', 'Holanda', 'Inglaterra', 'Irlanda', 'Japao', 'Mexico', 'Nomina', 'Republica_Tcheca']
Chromosome ref|NC_001134|, position in reference 259571, indel of length 5922 removed in relation to reference - ['Alemanha', 'Australia', 'Brasil', 'Colombia', 'EUA', 'Espanha', 'Estonia', 'Holanda', 'Inglaterra', 'Irlanda', 'Mexico', 'Romenia', 'Senegal']
Chromosome ref|NC_001136|, position in reference 1206698, indel of length 5923 removed in relation to reference - ['Africa_do_Sul', 'Alemanha', 'Argentina', 'Belgica', 'Brasil', 'Colombia', 'EUA', 'Espanha', 'Estonia', 'Holanda', 'Inglaterra', 'Irlanda', 'Japao', 'Romenia', 'Russia']
Chromosome ref|NC_001137|, position in reference 492690, indel of length 5731 removed in relation to reference - ['Alemanha', 'Irlanda']
Chromosome ref|NC_001139|, position in reference 535754, indel of length 5931 removed in relation to reference - ['Africa_do_Sul', 'Alemanha', 'Australia', 'Colombia', 'EUA', 'Espanha', 'Estonia', 'Holanda', 'Inglaterra', 'Irlanda', 'Japao', 'Mexico', 'Senegal']
Chromosome ref|NC_001139|, position in reference 707189, indel of length 5356 removed in relation to reference - ['Australia', 'Belgica', 'Brasil', 'Colombia', 'Espanha', 'Inglaterra', 'Irlanda', 'Japao', 'Mexico']
Chromosome ref|NC_001139|, position in reference 811440, indel of length 11883 removed in relation to reference - ['Alemanha', 'Argentina', 'Australia', 'Belgica', 'Colombia', 'EUA', 'Espanha', 'Holanda', 'Irlanda', 'Romenia']
Chromosome ref|NC_001144|, position in reference 215073, indel of length 5931 removed in relation to reference - ['Africa_do_Sul', 'Alemanha', 'Argentina', 'Australia', 'Belgica', 'Brasil', 'Colombia', 'EUA', 'Espanha', 'Estonia', 'Holanda', 'Inglaterra', 'Irlanda', 'Japao', 'Mexico', 'Republica_Tcheca', 'Romenia', 'Senegal']
Chromosome ref|NC_001146|, position in reference 519157, indel of length 5905 removed in relation to reference - ['Alemanha', 'Argentina', 'Australia', 'Belgica', 'Brasil', 'Colombia', 'EUA', 'Espanha', 'Inglaterra', 'Irlanda', 'Japao', 'Mexico', 'Republica_Tcheca', 'Romenia', 'Senegal']
Chromosome ref|NC_001147|, position in reference 703767, indel of length 6341 removed in relation to reference - ['Alemanha', 'Australia', 'Belgica', 'EUA', 'Espanha', 'Holanda', 'Inglaterra', 'Irlanda', 'Japao', 'Republica_Tcheca', 'Romenia', 'Senegal']
Chromosome ref|NC_001148|, position in reference 436479, indel of length 6650 removed in relation to reference - ['Alemanha', 'Belgica', 'Brasil', 'Colombia', 'Estonia', 'Holanda', 'Inglaterra', 'Irlanda', 'Japao']
Chromosome ref|NC_001136|, position in reference 645234, indel of length 6269 removed in relation to reference - ['Africa_do_Sul', 'Mexico', 'Republica_Tcheca', 'Senegal']
Chromosome ref|NC_001138|, position in reference 137907, indel of length 5964 removed in relation to reference - ['Mexico']
Chromosome ref|NC_001139|, position in reference 561837, indel of length 5923 removed in relation to reference - ['Alemanha', 'Australia', 'Brasil', 'Colombia', 'EUA', 'Espanha', 'Mexico', 'Republica_Tcheca']
Chromosome ref|NC_001139|, position in reference 568735, indel of length 5965 removed in relation to reference - ['Alemanha', 'Australia', 'Brasil', 'Colombia', 'EUA', 'Espanha', 'Mexico', 'Republica_Tcheca']
Chromosome ref|NC_001145|, position in reference 356998, indel of length 5919 removed in relation to reference - ['Africa_do_Sul', 'Alemanha', 'Brasil', 'Colombia', 'EUA', 'Espanha', 'Holanda', 'Inglaterra', 'Mexico', 'Senegal']
Chromosome ref|NC_001145|, position in reference 372690, indel of length 5931 removed in relation to reference - ['Africa_do_Sul', 'Alemanha', 'Australia', 'Brasil', 'Colombia', 'EUA', 'Espanha', 'Holanda', 'Inglaterra', 'Mexico', 'Republica_Tcheca', 'Senegal']
Chromosome ref|NC_001147|, position in reference 117697, indel of length 5931 removed in relation to reference - ['Alemanha', 'Australia', 'Brasil', 'Inglaterra', 'Mexico', 'Republica_Tcheca', 'Senegal']
Chromosome ref|NC_001148|, position in reference 436884, indel of length 6228 removed in relation to reference - ['Australia', 'EUA', 'Espanha', 'Mexico', 'Republica_Tcheca', 'Senegal']
Chromosome ref|NC_001136|, position in reference 645496, indel of length 5923 removed in relation to reference - ['Alemanha', 'Colombia', 'Inglaterra']
Chromosome ref|NC_001141|, position in reference 205214, indel of length 5433 removed in relation to reference - ['Alemanha', 'Australia', 'Brasil', 'Colombia', 'EUA', 'Espanha', 'Republica_Tcheca']
Chromosome ref|NC_001147|, position in reference 970281, indel of length 5964 removed in relation to reference - ['Colombia', 'EUA']
Chromosome ref|NC_001148|, position in reference 56445, indel of length 5929 removed in relation to reference - ['Colombia']
Chromosome ref|NC_001136|, position in reference 1154129, indel of length 5392 removed in relation to reference - ['Argentina']
Chromosome ref|NC_001136|, position in reference 513680, indel of length 5964 removed in relation to reference - ['Africa_do_Sul']
Chromosome ref|NC_001137|, position in reference 492645, indel of length 5753 removed in relation to reference - ['Africa_do_Sul']
Chromosome ref|NC_001143|, position in reference 322088, indel of length 16925 added in relation to reference - ['Africa_do_Sul']
Chromosome ref|NC_001141|, position in reference 205352, indel of length 5088 removed in relation to reference - ['Inglaterra']
Chromosome ref|NC_001143|, position in reference 322088, indel of length 16924 added in relation to reference - ['Inglaterra']
Chromosome ref|NC_001145|, position in reference 184164, indel of length 5919 removed in relation to reference - ['Australia']
Chromosome ref|NC_001147|, position in reference 594792, indel of length 5934 removed in relation to reference - ['Australia', 'Espanha']
Chromosome ref|NC_001136|, position in reference 1206697, indel of length 5923 removed in relation to reference - ['Senegal']
Chromosome ref|NC_001148|, position in reference 56450, indel of length 6266 removed in relation to reference - ['Belgica', 'Republica_Tcheca', 'Romenia']
Chromosome ref|NC_001147|, position in reference 969301, indel of length 7056 removed in relation to reference - ['Republica_Tcheca']
Chromosome ref|NC_001139|, position in reference 205412, indel of length 6180 added in relation to reference - ['Espanha']
```
