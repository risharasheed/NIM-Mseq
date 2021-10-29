#!/bin/bash
echo Current Date and Time is: `date +"%Y-%m-%d %T"`
args=("$@")
outputDir=${args[0]}
projDir=${args[1]}

porchop="${outputDir}/porechop_out/BC*.fastq"
bwa_human="${outputDir}/BWA_Human/*_combined_human_un_count.file"
bwa_silva="${outputDir}/BWA_Silva/*_combined_silva_un_count.file"
kraken="${outputDir}/Kraken_out/*_report"
blastn_f="${outputDir}/blastn_out/*_blastout.txt"
blastn="${outputDir}/blastn_out/*_blastdraft.txt"

count_report="${outputDir}/count_report"
len=04
rm -rf $count_report
mkdir $count_report

if [ "$(ls -A $porchop)" ]; then
 for eachfile in $porchop;
   do
     fname="$(basename ${eachfile})"
	 barcode=${fname:0:$len}
	 	 cat ${eachfile} |	 awk '{s++}END{print s/4}' ${eachrec} |awk -v z=${barcode} '{sub(/$/,"\t"z)} {print}' >> ${count_report}/porchop_cnt.txt
	done
fi	

if [ "$(ls -A $bwa_human)" ]; then
 for eachfile in $bwa_human;
   do
     fname="$(basename ${eachfile})"
	 barcode=${fname:0:$len}
	 awk '{print $1/4"\t"}' ${eachfile} |awk '{rec = rec $0} END{print rec}' | awk -v z=${barcode} '{sub(/$/,z)} {print}' >> ${count_report}/bwa_human_cnt.txt
	echo ${barcode} >>  ${count_report}/reffile
	done
fi

if [ "$(ls -A $bwa_silva)" ]; then
 for eachfile in $bwa_silva;
   do
     fname="$(basename ${eachfile})"
	 barcode=${fname:0:$len}
	 awk '{print $1/4"\t"}' ${eachfile} |awk '{rec = rec $0} END{print rec}' | awk -v z=${barcode} '{sub(/$/,z)} {print}' >> ${count_report}/bwa_silva_cnt.txt
	done
fi


if [ "$(ls -A $kraken)" ]; then
 for eachfile in $kraken;
   do
     fname="$(basename ${eachfile})"
	 barcode=${fname:0:$len}
	 head -2 ${eachfile} | awk 'NR==1 { if ($6 =="unclassified" ) print $2"\t"; else if ($6 =="root" ) print "0\t\n"$2"\t"; else 	print "0\t0\t";} NR==2 {if  ($6=="root") print$2"\t";}'   |awk '{rec = rec $0} END{print rec}' | awk -v z=${barcode} '{sub(/$/,z)} {print}' >> ${count_report}/kraken_cnt.txt
	done
fi

if [ "$(ls -A $blastn)" ]; then
 for eachfile in $blastn;
   do
     fname="$(basename ${eachfile})"
	 barcode=${fname:0:$len}
	 awk '{b++}END{print b}' ${eachfile} |awk -v z=${barcode} '{sub(/$/,"\t"z)} {print}' >> ${count_report}/blast_full_cnt.txt
	done
fi	

if [ "$(ls -A $blastn_f)" ]; then
 for eachfile in $blastn_f;
   do
     fname="$(basename ${eachfile})"
	 barcode=${fname:0:$len}
	 awk '{bo++}END{print bo}' ${eachfile} |awk -v z=${barcode} '{sub(/$/,"\t"z)} {print}' >> ${count_report}/blast_cnt.txt
	done
fi

python ${projDir}/countReport.py -init ${count_report}/porchop_cnt.txt -kra ${count_report}/kraken_cnt.txt -silva ${count_report}/bwa_silva_cnt.txt  -human ${count_report}/bwa_human_cnt.txt -blast ${count_report}/blast_cnt.txt -blastf ${count_report}/blast_full_cnt.txt  -inref ${count_report}/reffile -outDir ${count_report}
python ${projDir}/reportBarChart.py -inFile ${count_report}/countReport -outDir ${count_report} -prefix 'countBarChart'
