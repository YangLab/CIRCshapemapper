#!/usr/bin/env python
import argparse
import re


parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",type=str, help="sam file")
parser.add_argument("-p","--primer",type=str, help="primer")
parser.add_argument("-o","--out",type=str, help="output file")
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
                primerList.append("^"+primer)
                primerList.append("^"+revPrimer)
                primerPattern=re.compile(r"|".join(primerList))
        primerIn.close()       
        return primerPattern
    
    

def subPrimer(Line,rePattern):        
        tmp=Line.strip().split()
        new=tmp[:] #a new list 
        line2=tmp[9]
        subLine2=re.sub(rePattern,"",line2)
        L=len(line2)
        subL=len(subLine2)
        cutL=L-subL
        quality=tmp[10][cutL:]
        if len(subLine2)==0:
            new[2]="*"
        # define regular expression patterns for parsing CIGAR string
        item = re.compile(r"[0-9]+[M|I|D|S|N|H|P|=|X]")
        number = re.compile(r"([0-9]+)")
        #print "CIGAR: %s\nsplit CIGAR: %s"%(cigarString,str(cigarList))
        splitCigarString = item.findall(tmp[5].strip())
        cigarList = [number.split(x)[1:] for x in splitCigarString]
        newCigarList=cigarList[:]
        
        pos=tmp[3] # first matched base position 
        new[9]=subLine2 #seq
        new[10]=quality #seq quality
        MD=tmp[17]
        for i in tmp:   
            if "MD:Z" in i:
                MDsite=tmp.index(i)            
        MDlist=tmp[MDsite].split(":")   # to locat delection
        MDstr=re.findall(r'[0-9]+|[A|T|C|G|N|^]+',MDlist[-1])
        newMDstr=MDstr[:]
        newMDstr[0]=str(int(MDstr[0])-cutL)
        newMDstr="".join(newMDstr)
        MDlist[-1]=newMDstr
        newMD=":".join(MDlist)
        new[MDsite]=newMD
 
        if cutL >0:
            if cigarList[0][1]=="S" : # have no influence  on mutation parse
                new=tmp[:]
                newMD=tmp[MDsite]
                new[6]=tmp[6]
            elif cigarList[0][1]=="M" : # when match if effect mutation parse
                if int(cigarList[0][0]) >= cutL and int(MDstr[0])>= cutL:
                    if int(cigarList[0][0]) == cutL :
                        if len( cigarList)==2 and cigarList[1][1] == "S":
                            new[2]="*"
                        else:
                            new[ 3]=str(int(tmp[3])+cutL)
                            newCigarList=cigarList[1:]
                            new[7]=tmp[7]
                            new[MDsite]=newMD
                    
                    elif int(cigarList[0][0]) > cutL :
                        new[3]=str(int(tmp[3])+cutL)
                        newCigarList[0][0]=str(int(cigarList[0][0])-cutL)
                        new[7]=tmp[7]
                        new[MDsite]=newMD           
                else:
                    new=tmp[:]
                    newMD=tmp[MDsite]                   
                    new[6]=tmp[6]
    
                
            newCigar=""
            for i in newCigarList:
                newCigar=newCigar+"".join(i)    
            new[5]=newCigar
            new[MDsite]=newMD
        else:
            new=tmp[:]
        #print "\t".join(new)
        fastqOut.write("\t".join(new)+"\n")
        
def main():
        fastqIn=open(args.input)      
        #rePattern=parsePrimer("primer_circ_Seq.txt")
        rePattern=parsePrimer(args.primer)

        for line in fastqIn:
            if "@" in line:
                fastqOut.write(line)
            else:
                subPrimer(line,rePattern)
                
        fastqIn.close()
        fastqOut.close()
       
fastqOut=open(args.out,'w')
        
main()
