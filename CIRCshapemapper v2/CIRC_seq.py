#!/usr/bin/env python
import argparse
import re



parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",type=str, help="fa file")
parser.add_argument("-p","--primer",type=str, help="primer name, primerName_F.txt and primerName_R.txt")
parser.add_argument("-o","--output",type=str, help="output.fa")
args = parser.parse_args()


def revSeq(seq):
        nuc1="atcguATCGU"
        nuc2="tagcaTAGCA"
        trans=str.maketrans(nuc1,nuc2)
        compSeq=seq.translate(trans)
        revComp=compSeq[::-1]
        return revComp

def parsePrimer(primerFile):
        primerIn=open(primerFile)
        primerList=[]
        for line in primerIn:
                primer=line.strip().upper()
                revPrimer=revSeq(primer)
                primerList.append(revPrimer)
        return primerList


def main():
        primerF=[]
        f=open(args.primer+"_F.txt")
        for i in f:
            primerF.append(i.strip())
        f.close()
        
        primerR=parsePrimer(args.primer+"_R.txt")

        fastaIn=open(args.input)
        for line in fastaIn:
            name=line.strip()
            tmp=next(fastaIn).strip().upper()
            seq=tmp+tmp
        L=int(len(seq)/2)
        #print seq
        Flist=[]
        for i in primerF:
            a=seq.find(i)
            if a != -1 :
                Flist.append(a)
        
        Rlist=[]
        for i in primerR:
            a=seq[L:].find(i)
            if a != -1:
                Rlist.append(a+L+len(i))
        #print primerF,Flist
        #print primerR,Rlist
        Out=open(args.output+"_shape.fa","w")
        Out.write(name+"\n")
        Out.write(seq[min(Flist):max(Rlist)])
        Out.close()
        #print(name)
        #print(seq[min(Flist):max(Rlist)])
        #name,subseq=subPrimer(args.input,rePattern)
        #name,subseq=subPrimer("POLR2A.fa",rePattern)
        #print name,subseq
        
main()

