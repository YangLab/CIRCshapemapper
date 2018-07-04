#!/usr/bin/env python
#########################################################################
# File Name: circSHAPEsplit.py
# Author: FangNan
# mail: nanfang@picb.cn
#########################################################################

from string import maketrans
import argparse
import re

parser = argparse.ArgumentParser("Split junction reads")
parser.add_argument("-s", "--sam", type=str, help="input the mapped sam file")
parser.add_argument("-o", "--out", type=str, help="output file name, out_split.fastq. Output a single-end fastq file")
args = parser.parse_args()


def revSeq(seq):
	nuc1="atcgATCG"
	nuc2="tagcTAGC"
	trans=maketrans(nuc1,nuc2)
	compSeq=seq.translate(trans)
	revComp=compSeq[::-1]
	return revComp
											

									
def main():
	fileOut=open(args.out+"_split.fastq",'w')
	fileIn=open(args.sam)
	for line in fileIn:
		if line[0] == "@":
			continue
		
		read=line.rstrip().split()
		if read[3] == '*':
			continue
		
		flag=bin(int(read[1]))[2:].zfill(12)
		if flag[-5]=='1':
			mate='_2'
		else:
			mate='_1'
		
		seq=[read[5],read[9],read[10],mate,read[0]]    
		if 'S' in seq[0] :                             
			item=re.compile(r"[0-9]+[M|S|D|I]")
			splitCigar=item.findall(seq[0])
			if 'S' in splitCigar[0] and 'S' not in splitCigar[-1]:
				leftLength=int(splitCigar[0].strip("S"))
				leftSeq=seq[1][:leftLength]
				leftQua=seq[2][:leftLength]
				rightSeq=seq[1][leftLength:]
				rightQua=seq[2][leftLength:]
				fileOut.write("@"+seq[-1]+seq[-2]+"1\n"+leftSeq+"\n+\n"+leftQua+"\n")
				fileOut.write("@"+seq[-1]+seq[-2]+"3\n"+rightSeq+"\n+\n"+rightQua+"\n")
			if 'S' in splitCigar[-1] and 'S' not in splitCigar[0]:
				rightLength=int(splitCigar[-1].strip("S"))
				rightSeq=seq[1][-rightLength:]
				rightQua=seq[2][-rightLength:]
				leftSeq=seq[1][:-rightLength]
				leftQua=seq[2][:-rightLength]
				fileOut.write("@"+seq[-1]+seq[-2]+"3\n"+rightSeq+"\n+\n"+rightQua+"\n")
				fileOut.write("@"+seq[-1]+seq[-2]+"1\n"+leftSeq+"\n+\n"+leftQua+"\n")
			if 'S' in splitCigar[0] and 'S' in splitCigar[-1]:
				leftLength=int(splitCigar[0].strip("S"))
				rightLength=int(splitCigar[-1].strip("S"))
				leftSeq=seq[1][:leftLength]
				leftQua=seq[2][:leftLength]
				rightSeq=seq[1][-rightLength:]
				rightQua=seq[2][-rightLength:]
				middleSeq=seq[1][leftLength:-rightLength]
				middleQua=seq[2][leftLength:-rightLength]                
				fileOut.write("@"+seq[-1]+seq[-2]+"3\n"+rightSeq+"\n+\n"+rightQua+"\n")
				fileOut.write("@"+seq[-1]+seq[-2]+"1\n"+leftSeq+"\n+\n"+leftQua+"\n")
				fileOut.write("@"+seq[-1]+seq[-2]+"2\n"+middleSeq+"\n+\n"+middleQua+"\n")
		
        else:
			fileOut.write("@"+read[0]+"\n"+read[9]+"\n+\n"+read[10]+"\n")
	fileOut.close()
	fileIn.close()
																																										
main()
																																                 
