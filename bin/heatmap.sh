#!/bin/bash
echo Current Date and Time is: `date +"%Y-%m-%d %T"`
args=("$@")
outputDir=${args[0]}
projDir=${args[1]}
#outputDir="/home/risha_rasheed88/output/nextflow/run1"
#projDir="/home/risha_rasheed88/nextflow/NIM-Mseq/bin"
#
map_name="heatMap"
heatmap_report="${outputDir}/heatmap"
mapfl="$outputDir/pathogen_report/BC*_r.heatmapM"
len=04
#
rm -rf $heatmap_report
mkdir $heatmap_report
#
echo " ------------------------------------------------------------------------------------------"
echo "HeatMap pre file process is in progrss..... "
#
if [ "$(ls -A $mapfl)" ]; then
 for eachfile in $mapfl;
   do
     fname="$(basename ${eachfile})"
	 barcode=${fname:0:$len}
	 cat ${eachfile} |	 sort -k1,1 -k3r | sort -t$'\t' -u -k1,1 > ${heatmap_report}/${barcode}_uniq_spname.txt
	done
fi
#
cat 	${heatmap_report}/*_uniq_spname.txt>  ${heatmap_report}/heatmapM.txt
cat ${heatmap_report}/heatmapM.txt | sort  -k1,1 > ${heatmap_report}/heatmapM_sorted.txt
#
python ${projDir}/heatmapPlot.py -inFile ${heatmap_report}/heatmapM_sorted.txt -outDir ${heatmap_report} -prefix ${map_name}
rm -rf  ${heatmap_report}/*_uniq_spname.txt
