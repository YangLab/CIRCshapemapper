#!/usr/bin/env bash
helpdoc(){
    cat <<EOF
Description:
	CIRCshapemapper version 2 pipeline. Merge two primer product sequence. Suit for single circRNA. multipal circRNAs should run respectively. 

Usage:

    $0 -f circRNA.fa -t circRNA -p  circPrimer -P primer_pair -N NAI_R1.fastq -n NAI_R2.fastq -D DC_R1.fastq -d DC_R2.fastq -U DMSO_R1.fastq -u DMSO_R2.fastq -M 1000 -o outDir  
Option:

-f     -- circRNA fa file
-t     -- target name in shape
-p     -- circRNA primer for trim;txt formate
-P     -- circRNA primer to generate new reference; txt formate
-N     -- NAI samples R1
-n     -- NAI samples R2
-D     -- Denatured samples R1
-d     -- Denatured samples R2
-U     -- Untreated sample R1
-u     -- Untreated sample R2
-M     -- min-deep default 1000
-o     -- output document
-v     -- verbose mode, output command but out execute
-h     -- help information

Version: 2.0

Author: NanFang 
EOF
}

BOLD_s="\x1b[1;39;49m"
RED_s="\x1b[0;31;107m"
END="\x1b[0m"

declare verbose=0

if [ $# = 0 ]
then
    helpdoc
    exit 1
fi

while getopts ":f:t:p:P:N:n:D:d:U:u:M:o:v" opt; do
  case ${opt} in
    f )
       Fasta=${OPTARG};
      ;;
  	t )
      target=${OPTARG};
      ;;
	p )
      primer=${OPTARG};
	   ;;
	P )
      Pair=${OPTARG};
	  ;;
    N )
        NAI1=${OPTARG};
      ;;
    n )
        NAI2=${OPTARG};
      ;;
    D )
        DC1=${OPTARG};
      ;;
    d )
        DC2=${OPTARG};
      ;;
    U )
        DMSO1=${OPTARG};
      ;;
    u )
        DMSO2=${OPTARG};
      ;;
    M )
		deepth=1000;
        deepth=${OPTARG};
      ;;
    o )
        OutDir=${OPTARG};
      ;;
    v )
        verbose=1
      ;;
    \? )
        helpdoc;
      ;;
  esac
done

shift $((OPTIND -1))

if [[ "Fastqa" == "" ]]; then
    echo -e "${RED_s}Error: Specify -f${END}";
    helpdoc;
    exit;
elif [[ "target" == "" ]]; then
    echo -e "${RED_s}Error: Specify -t${END}";
    helpdoc;
    exit;
elif [[ "primer" == "" ]]; then
    echo -e "${RED_s}Error: Specify -p${END}";
    helpdoc;
    exit;
elif [[ "Pair" == "" ]]; then
    echo -e "${RED_s}Error: Specify -p${END}";
    helpdoc;
    exit;
elif [[ "${NAI1}" == "" ]]; then
    echo -e "${RED_s}Error: Specify -N${END}";
    helpdoc;
    exit;
elif [[ "${NAI2}" == "" ]]; then
    echo -e "${RED_s}Error: Specify -n${END}";
    helpdoc;
    exit;
elif [[ "${DC1}" == "" ]]; then
    echo -e "${RED_s}Error: Specify -D${END}";
    helpdoc;
    exit;
elif [[ "${DC2}" == "" ]]; then
    echo -e "${RED_s}Error: Specify -d${END}";
    helpdoc;
    exit;
elif [[ "${DMSO1}" == "" ]]; then
    echo -e "${RED_s}Error: Specify -U${END}";
    helpdoc;
    exit;
elif [[ "${DMSO2}" == "" ]]; then
    echo -e "${RED_s}Error: Specify -u${END}";
    helpdoc;
    exit;
elif [[ "${deepth}" == "" ]]; then
    echo -e "${RED_s}Error: Specify -M${END}";
    helpdoc;
    exit;
elif [[ "${OutDir}" == "" ]]; then
    echo -e "${RED_s}Error: Specify -o${END}";
    helpdoc;
    exit;
fi



CMD="mkdir ${OutDir}"
echo -e "  >> ${CMD}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi

## correct fa 
CMD="CIRC_seq.py -i ${Fasta} -p ${Pair} -o ${OutDir}/${Fasta%.*}"
echo -e "  >> ${CMD}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi



## First round


CMD="shapemapper --name ${target} --target  ${OutDir}/${Fasta%.*}_shape.fa --out TempRound --verbose --serial  --nproc 16 --min-depth ${deepth} --modified --R1 ${NAI1} --R2 ${NAI2} --untreated --R1 ${DMSO1} --R2 ${DMSO2} --denatured --R1 ${DC1} --R2 ${DC2}"
echo -e "  >> ${CMD}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi



## Trim primer
CMD="SM2_Trim_primer_v2_py3.py -i shapemapper_temp/${target}/Modified/BowtieAligner/${target}_Modified_BowtieAligner_aligned.sam -p ${primer} -o shapemapper_temp/${target}/Modified/BowtieAligner/${target}_Modified.sam"
echo -e "  >> ${CMD}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi
	
CMD="SM2_Trim_primer_v2_py3.py -i shapemapper_temp/${target}/Denatured/BowtieAligner/${target}_Denatured_BowtieAligner_aligned.sam -p ${primer} -o shapemapper_temp/${target}/Denatured/BowtieAligner/${target}_Denatured.sam"
echo -e "  >> ${CMD}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi

CMD="SM2_Trim_primer_v2_py3.py -i shapemapper_temp/${target}/Untreated/BowtieAligner/${target}_Untreated_BowtieAligner_aligned.sam -p ${primer} -o shapemapper_temp/${target}/Untreated/BowtieAligner/${target}_Untreated.sam"
echo -e "  >> ${CMD}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi







## CIRC shape

CMD="mkdir shapemapper_temp/${target}_trimed";
echo -e "  >> ${CMD}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi

Len=`awk '{if(NR==2)print length($0)}' ${OutDir}/${Fasta%.*}_shape.fa` 
echo -e "  >> ${Len}";

CMD=" shapemapper_mutation_parser -i shapemapper_temp/${target}/Modified/BowtieAligner/${target}_Modified.sam -o shapemapper_temp/${target}_trimed/${target}_MutationParser_Modified_parsed_mutations.mut -w -m 10 "
echo -e "  >> ${CMD}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi

CMD=" shapemapper_mutation_counter -i shapemapper_temp/${target}_trimed/${target}_MutationParser_Modified_parsed_mutations.mut  -w -c shapemapper_temp/${target}_trimed/${target}_MutationCounter_Modified_mutations.txt --exclude_3prime 1 --max_internal_match 5 --min_qual 30 --length ${Len} "
echo -e "  >> ${CMD}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi


CMD=" shapemapper_mutation_parser -i shapemapper_temp/${target}/Untreated/BowtieAligner/${target}_Untreated.sam -o shapemapper_temp/${target}_trimed/${target}_MutationParser_Untreated_parsed_mutations.mut -w -m 10 "
echo -e "  >> ${CMD}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi

CMD=" shapemapper_mutation_counter -i shapemapper_temp/${target}_trimed/${target}_MutationParser_Untreated_parsed_mutations.mut  -w -c shapemapper_temp/${target}_trimed/${target}_MutationCounter_Untreated_mutations.txt --exclude_3prime 1 --max_internal_match 5 --min_qual 30 --length ${Len} "
echo -e "  >> ${CMD}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi


CMD=" shapemapper_mutation_parser -i shapemapper_temp/${target}/Denatured/BowtieAligner/${target}_Denatured.sam -o shapemapper_temp/${target}_trimed/${target}_MutationParser_Denatured_parsed_mutations.mut -w -m 10 "
echo -e "  >> ${CMD}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi

CMD=" shapemapper_mutation_counter -i shapemapper_temp/${target}_trimed/${target}_MutationParser_Denatured_parsed_mutations.mut  -w -c shapemapper_temp/${target}_trimed/${target}_MutationCounter_Denatured_mutations.txt --exclude_3prime 1 --max_internal_match 5 --min_qual 30 --length ${Len} "
echo -e "  >> ${CMD}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi



CMD="python3 make_reactivity_profiles.py --fa ${OutDir}/${Fasta%.*}_shape.fa --rna ${Fasta%.*} --counts shapemapper_temp/${target}_trimed/${target}_MutationCounter_Modified_mutations.txt  shapemapper_temp/${target}_trimed/${target}_MutationCounter_Untreated_mutations.txt shapemapper_temp/${target}_trimed/${target}_MutationCounter_Denatured_mutations.txt --out shapemapper_temp/${target}_trimed/${target}_ProfileHandler_CalcProfile_profile_temp.txt --mindepth ${deepth} --maxbg 0.05 " 
echo -e "  >> ${CMD}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi


## change profile to right sequence 
CMD="Circ_NewStratgy.py -i shapemapper_temp/${target}_trimed/${target}_ProfileHandler_CalcProfile_profile_temp.txt --RNA ${Fasta}  -o shapemapper_temp/${target}_trimed/${target}_ProfileHandler_CalcProfile_profile --mindepth ${deepth} --maxbg 0.05 "
echo -e "  >> ${CMD}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi

	


	
CMD=" python3 normalize_profiles.py --warn-on-error --tonorm shapemapper_temp/${target}_trimed/${target}_ProfileHandler_CalcProfile_profile_pre-profile.txt --normout ${OutDir}/${target}_profile.txt"
echo -e "  >> ${CMD}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi
	
CMD="python3 tab_to_shape.py --infile ${OutDir}/${target}_profile.txt --shape ${OutDir}/${target}.shape  --map  ${OutDir}/${target}.map"
echo -e "  >> ${CMD}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi

CMD="python3 render_figures.py --infile ${OutDir}/${target}_profile.txt --mindepth ${deepth} --maxbg 0.05 --plot ${OutDir}/${target}_profiles.pdf --hist ${OutDir}/${target}_histograms.pdf --title RNA:${Fasta%.*}  > ${OutDir}/${target}.log "
echo -e "  >> ${CMD}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi


