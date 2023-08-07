CIRCshapemapper v2 
============================================
CIRCshapemapper2  is a pipeline to analyze the circSHAPE-MaP data. Using this pipeline, you could product a SHAPE reactivity file of circRNA, which was used for modeling RNA secondary structure.

A schematic flow shows the pipeline
-------------------------------------

![image](https://github.com/YangLab/CIRCshapemapper/blob/master/CIRCshapemapper%20v2/CIRCshapemapper%20v2%20pipeline.jpg)


Requirement
------------------------------------
* ShapeMapper-2.1.3   
The requirement of ShapeMapper-2.1.3 installation please read the README.md in ShapeMapper2 package. 
* python3

Install
------------------------------------
Install ShapeMapper-2.1.3 first.
please add PATH of shapemapper_mutation_counter and shapemapper_mutation_parser from ’/shapemapper-2.1.3/bin/‘ folder in your .bashrc or .profile

All add script of CIRCshapemapper v2 in your PATH. The file 'normalize_profiles.py' ,'make_reactivity_profiles.py' , 'tab_to_shape.py' and 'render_figures.py' are from  shapemapper-2.1.3 .



Usage
----------------------------------
#### CIRCshapemapper_V2.sh
Description:
	CIRCshapemapper version 2 pipeline. Merge two primer product sequence. Suit for single circRNA. multipal circRNAs should run manully.

Usage:

    /data/rnomics6/nanfang/script/CIRCshapemapper_V2.sh -f circRNA.fa -t circRNA -p  circPrimer.txt -P primer_pair -N NAI_R1.fastq -n NAI_R2.fastq -D DC_R1.fastq -d DC_R2.fastq -U DMSO_R1.fastq -u DMSO_R2.fastq -M 1000 -o outDir
Option:

 -f     -- circRNA fa file
 -t     -- target name in shape
 -p     -- circRNA primer for trim;txt formate
 -P		-- circRNA primer for generate new reference;txt formate
 -N     -- NAI samples R1
 -n     -- NAI samples R2
 -D     -- Denatured samples R1
 -d     -- Denatured samples R2
 -U     -- Untreated sample R1
 -u     -- Untreated sample R2
 -M     -- min-deep defualt 1000
 -o     -- output document
 -v     -- verbose mode, output command but out execute
 -h     -- help information

Version: 2.0




#### SM2_Trim_primer_v2_py3.py

		usage: Trim primers from alignment file [-h] [-i INPUT] [-p PRIMER] [-o OUT]
		optional arguments:
		-h, --help            show this help message and exit
		-i INPUT, --input INPUT sam file
		-p PRIMER, --primer PRIMER primer file 
		-o OUT, --out OUT    output file
![image](https://github.com/YangLab/circSHAPEmapper/blob/master/CIRCshapemapper%20v2/003.jpeg)


File formate
----------------------------------

-p  filename   -- circRNA primer for trim;txt formate
each primer in one line, include forward primers and reverse primers
>TGTGATGGGGACATTGTTAT
>AACTTGCACCTGCCACAGTC
>ACAGTGTGAGTGTGTCCTGCACAATA
>ACAGTGGACCATGGGAGAATGCGGAC

ACAGTG is barcode for seperate data. If your primer don't contain barcode, just recode your primer is ok. 



-P	primer_pair	-- circRNA primer for genarate new reference;txt formate 
include two files

primer_pair_F.txt
>TGTGATGGGGACATTGTTAT
>AACTTGCACCTGCCACAGTC

primer_pair_R.txt
>TGAGTGTGTCCTGCACAATA
>GACCATGGGAGAATGCGGAC
 
 
 

## Citation

## License
Copyright (C) 2019 YangLab. Licensed GPLv3 for open source use or contact YangLab (yanglab@@picb.ac.cn) for commercial use.
