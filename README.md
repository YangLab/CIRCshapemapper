circSHAPEmapper
============================================
circSHAPEmapper is a pipeline to analyze the circSHAPE-MaP data. Using this pipeline, you could product  a SHAPE reactivity file of circRNA, which was used for modeling RNA secondary structure.

A schematic flow shows the pipeline
-------------------------------------



Requirement
------------------------------------
* ShapeMapper2   http://www.chem.unc.edu/rna/software.html

The requirement of ShapeMapper2 installation please read the README.md in ShapeMapper2 package. 

Usage
----------------------------------
#### circSHAPEsplit
		usage: Split junction reads [-h] [-s SAM] [-o OUT]	
		optional arguments:
		-h, --help         show this help message and exit
		-s SAM, --sam SAM  input the mapped sam file
		-o OUT, --out OUT  output the single-end fastq file


#### circSHAPEtrim

		usage: Trim primers from alignment file [-h] [-i INPUT] [-p PRIMER] [-o OUT]
		optional arguments:
		-h, --help            show this help message and exit
		-i INPUT, --input INPUT sam file
		-p PRIMER, --primer PRIMER primer file 
		-o OUT, --out OUT    output file
