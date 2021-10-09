#!/usr/bin/env nextflow

println "NIM-Mseq(Nanopore Informatics for Metagenomics of Microbial sequences) - Taxonomic classification, similarity search and reference-based assembly of microbial sequences"
println " "
println "Risha rasheed                             "
println "*********************************************************************"


inputfile = params.inputDir
demultiplex_ch1 = Channel
    .fromPath(inputfile)

if(params.Demultiplex == true){
process Demultiplex {
	
	publishDir "$outputDir/porechop_out", mode: "copy", pattern: "*"
	 
    output:
    file ("BC*.fastq")  into demultiplex_ch
   
     script:
     """
	
		porechop -i ${inputfile} -b ./ -t ${threads}
		
     """
 } 

 demultiplex_Out = demultiplex_ch
            .flatten()
            .map {file -> tuple(file.simpleName, file)}
}

if (params.Demultiplex == false){
	demultiplex_Out = demultiplex_ch1
            .flatten()
            .map {file -> tuple(file.simpleName, file)}
}

demultiplex_Out.into { inFastq_NanoPlot; inFastq_nanoFilt}

if(params.quality_process == true){
 process Quality_Nanoplot {
	 label 'NIMMseq'
	 
	 publishDir "$outputDir/quality_process", mode: "copy", pattern: "*"
     tag {datasetID}
	 
	 input:
	 set datasetID, file(fastqFl) from inFastq_NanoPlot
	 
	 output:
	 set datasetID, file("${datasetID}_*") into quality_process
	 
	 script:
     """
	 NanoPlot --fastq ${fastqFl} -p ${datasetID}_  -t ${threads}
	 """
	 
 }
 
}

if(params.filter_rocess == true){
process Filter_NanoFilt {
     label 'NIMMseq'
	 
	 
     tag {datasetID}
	 
	 input:
	 set datasetID, file(fastqFl) from inFastq_nanoFilt
	 
	 output:
	 set datasetID, file("${datasetID}_*") into nanoFilt_trimmed
	 
	 script:
     """
	 cat ${fastqFl} | NanoFilt -q 7 -l 75  > ${datasetID}_trimmed.fastq
	 """
}
}

if (params.filter_rocess == false){
	nanoFilt_trimmed = inFastq_nanoFilt
}


if(params.bwa_process == true){
process BWA_Human {
	 label 'NIMMseq'

     publishDir "$outputDir/BWA_Human", mode: "copy", pattern: "${datasetID}_combined_human_*.f*"
     tag {datasetID}
	 maxForks 1
	 input:
	 set datasetID, file(fastqFl) from nanoFilt_trimmed
	 
	 output:
	 set datasetID, file("${datasetID}_combined_human_unmapped.fastq") into Bwa_out
	 set datasetID, file("${datasetID}_combined_human_un_count.file")  into Bwa_out_sub
	 
	 script:
     """
		bwa mem -t ${threads} ${human_g1k} ${fastqFl} > ${datasetID}_combined_human.sam
		samtools view -b ${datasetID}_combined_human.sam --threads ${threads} > ${datasetID}_combined_human.bam
		samtools view -b -f 4 ${datasetID}_combined_human.bam --threads ${threads} > ${datasetID}_combined_human_unmapped.bam
		samtools view -b -F 4 ${datasetID}_combined_human.bam --threads ${threads} > ${datasetID}_combined_human_mapped.bam
		samtools fastq ${datasetID}_combined_human_unmapped.bam --threads ${threads} > ${datasetID}_combined_human_unmapped.fastq
		samtools fastq ${datasetID}_combined_human_mapped.bam --threads ${threads} > ${datasetID}_combined_human_mapped.fastq
		wc -l ${datasetID}_combined_human_unmapped.fastq >> ${datasetID}_combined_human_un_count.file
		wc -l ${datasetID}_combined_human_mapped.fastq >> ${datasetID}_combined_human_un_count.file
		wc -l ${fastqFl} >> ${datasetID}_combined_human_un_count.file
		
	 """
}

}

if(params.bwa_process == false){
	Bwa_out = nanoFilt_trimmed
}

if(params.silva_process == true){
process BWA_Silva_db {
	 label 'NIMMseq'

     publishDir "$outputDir/BWA_Silva", mode: "copy", pattern: "${datasetID}_combined_silva_*.f*"
     tag {datasetID}
	 maxForks 1
	 input:
	 set datasetID, file(inFastq) from Bwa_out
	 
	 output:
	 set datasetID, file("${datasetID}_combined_silva_unmapped.fastq") into Silva_out
	 set datasetID, file("${datasetID}_combined_silva_un_count.file"), file("${datasetID}_combined_silva_mapped.fastq")  into Silva_out_sub
	 
	 script:
     """
		bwa mem -t ${threads} ${silva_db} ${inFastq} > ${datasetID}_combined_silva.sam
		samtools view -bS ${datasetID}_combined_silva.sam --threads ${threads} > ${datasetID}_combined_silva.bam
		samtools view -b -f 4 ${datasetID}_combined_silva.bam --threads ${threads} > ${datasetID}_combined_silva_unmapped.bam
		samtools view -b -F 4 ${datasetID}_combined_silva.bam --threads ${threads} > ${datasetID}_combined_silva_mapped.bam
		samtools fastq ${datasetID}_combined_silva_unmapped.bam --threads ${threads} > ${datasetID}_combined_silva_unmapped.fastq
		samtools fastq ${datasetID}_combined_silva_mapped.bam --threads ${threads} > ${datasetID}_combined_silva_mapped.fastq
		wc -l ${datasetID}_combined_silva_unmapped.fastq >> ${datasetID}_combined_silva_un_count.file
		wc -l ${datasetID}_combined_silva_mapped.fastq >> ${datasetID}_combined_silva_un_count.file
		wc -l ${inFastq} >> ${datasetID}_combined_silva_un_count.file
		
		
	 """
}

}

if(params.silva_process == false){
	Silva_out = Bwa_out
}

 
Silva_out.into { kraken_silva; kraken_parasite}

if(params.kraken2_process == true ){
process Classification_Kraken2_DB1 {
	 label 'NIMMseq'
	 
     publishDir "$outputDir/Kraken_out", mode: "copy", pattern: "${datasetID}_*"
	 
	 tag {datasetID}
	 
	 maxForks 1
	 
	 input:
	 set datasetID, file(inFastq) from  kraken_silva
	 
	 output:
	 set datasetID,  file("${datasetID}_krona.html"), file("${datasetID}_report"), file("${datasetID}_kraken"), file("${datasetID}_classified"), file("${datasetID}_unclassified") into Kraken_fasta_sub
	 set datasetID, file("${inFastq}") into Kraken_blast
	 	 
	 script:
     """
		
		kraken2 --db ${KRAKEN2_DB} --threads ${threads} --fastq-input ${inFastq} > ${datasetID}_kraken --report ${datasetID}_report  --classified-out ${datasetID}_classified --unclassified-out ${datasetID}_unclassified  
		cat ${datasetID}_kraken | cut -f 2,3 >  ${datasetID}_krona.input 
		ktImportTaxonomy  ${datasetID}_krona.input  -o ${datasetID}_krona.html
	    
	 """
}

}

if(params.kraken2_process == false ){
	Kraken_blast = kraken_silva
}



if(params.blastn_process == true){
process Similarity_Search_Blastn  {

	 label 'NIMMseq'

     publishDir "$outputDir/blastn_out", mode: "copy", pattern: "${datasetID}_*.txt"
	 
     tag {datasetID}
	 
	 maxForks 1
	 
	 input:
	 set datasetID, file(inFastq) from  Kraken_blast
	 
	 output:
	 set datasetID, file("${inFastq}"), file("${datasetID}_blastout.txt") into blastn_out
	 set datasetID, file("${datasetID}_blastdraft.txt") into blastb_out_sub
	 
	 script:
     """
	    
	    seqtk seq -a ${inFastq} > ${datasetID}_combined_silva_unmapped.fasta
		blastn -query ${datasetID}_combined_silva_unmapped.fasta -db ${BLAST_db_nrnt} -evalue 1e-10 -max_target_seqs 1 -num_threads ${threads} -perc_identity 80 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs staxids sscinames" > ${datasetID}_blastdraft.txt
        cat ${datasetID}_blastdraft.txt | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge | tr '|' '\t' | cut -f1,5- | awk '!/(uncultured|Homo sapiens)/' > ${datasetID}_blastout.txt
	   	   
	 """
}

}

if(params.blast_ref_search == true){
process Blast_Reference_Search  {
	 label 'NIMMseq'

     publishDir "$outputDir/blastref_out", mode: "copy", pattern: "${datasetID}_n*"
	 
     tag {datasetID}
	 
	 maxForks 1
	 
	 input:
	 set datasetID, file(inFastq), file("${datasetID}_blastout.txt") from blastn_out
	 
	 output:
	 set datasetID, file("${inFastq}"),  file("${datasetID}_blastout.txt"), file("${datasetID}_ntref.fasta") into blastref_out
	 set datasetID, file("${datasetID}_nrefAccNumber") into blastref_out_sub
	 
	 script:
     """
	   cat ${datasetID}_blastout.txt | cut -f 2 | sort | uniq > ${datasetID}_nrefAccNumber
	   seqtk  subseq $BLAST_Fasta_nt ${datasetID}_nrefAccNumber > ${datasetID}_ntref.fasta
			   
	 """

}

}

if (params.mapping_process == true){
process Alignment_Minimap2  {
	 label 'NIMMseq'

     publishDir "$outputDir/Alignment_Out", mode: "copy", pattern: "${datasetID}.*"
	 
     tag {datasetID}
	 
	 maxForks 1
	 
	 input:
	 set datasetID, file(inFastq), file("${datasetID}_blastout.txt"),  file("${datasetID}_ntref.fasta") from blastref_out
	
	 
	 output:
	 set datasetID, file("${inFastq}"), file("${datasetID}_ntref.fasta"),  file("${datasetID}_blastout.txt"),  file("${datasetID}.mapStats"), file("${datasetID}.newsort.bam") into mapping_out
	 set datasetID, file("${datasetID}.mapping"), file("${datasetID}.mapA.Report"), file("${datasetID}.mapB.Report") into mapping_out_sub
	 script:
     """
   	    minimap2 -ax map-ont --frag=yes --secondary=no ${datasetID}_ntref.fasta ${inFastq} --split-prefix ${datasetID} -t ${threads} > ${datasetID}.sam
		samtools view  -b ${datasetID}.sam --threads ${threads} > ${datasetID}.bam 
		samtools sort ${datasetID}.bam -o ${datasetID}.newsort.bam --threads ${threads}
        samtools index ${datasetID}.newsort.bam -@ ${threads}
        python ${projDir}/bamReport.py -inFile ${datasetID}.newsort.bam -splID ${datasetID}
		
		samtools flagstat ${datasetID}.newsort.bam --threads ${threads} > ${datasetID}.mapA.Report
		samtools stats ${datasetID}.newsort.bam --threads ${threads} > ${datasetID}.mapB.Report
		   
	 """

}

}
if(params.report_pathogen == true){
process Report {
	 label 'NIMMseq'

     publishDir "$outputDir/pathogen_report", mode: "copy", pattern: "${datasetID}_r.*"
	 
     tag {datasetID}
	 
	 input:
	 set datasetID, file(inFastq), file("${datasetID}_ntref.fasta"), file("${datasetID}_blastout.txt"), file("${datasetID}.mapStats"), file("${datasetID}.newsort.bam") from mapping_out
	
	 
	 output:
	 set datasetID, file("${inFastq}"), file("${datasetID}_ntref.fasta"),   file("${datasetID}.newsort.bam") into pathogen_report
	 set datasetID, file("${datasetID}_r.pathogenReport"), file("${datasetID}_r.paf"),  file("${datasetID}_r.detailed.Report.txt"), file("${datasetID}_r.heatmapM"), file("${datasetID}_r.ident_removed.txt") into pathogen_report_sub
	 script:
     """
   	    minimap2 -cx map-ont --frag=yes  --secondary=no ${datasetID}_ntref.fasta ${inFastq} --split-prefix ${datasetID} -t ${threads} > ${datasetID}_r.paf
		cat  ${datasetID}_r.paf | cut -f 1,2,6,7,8,9,10,11,12 >  ${datasetID}_paf_filtered.txt
		python ${projDir}/detailedReport.py -inRef ${datasetID}_paf_filtered.txt -splID ${datasetID}
		cat  ${datasetID}_blastout.txt | cut -f 2,15,16 | sort -u -k1,1 > ${datasetID}.refAcc
		Rscript ${projDir}/reportIdentity.R ${datasetID}_r.ident_reseq.txt  ${datasetID}
		python ${projDir}/pathogenReport.py -inIdent ${datasetID}.identSummary -inMap ${datasetID}.mapStats -inRef ${datasetID}.refAcc -splID ${datasetID} 
		
					   
	 """
}

}


if(params.consensus_process == true){
process Consensus {
	 label 'NIMMseq'

     publishDir "$outputDir/consensus_out", mode: "copy", pattern: "${datasetID}_con*"
	 
     tag {datasetID}
	 
	 input:
	 set datasetID, file(inFastq), file("${datasetID}_ntref.fasta"), file("${datasetID}.newsort.bam") from pathogen_report
	 
	 output:
	 set datasetID, file("${inFastq}"), file("${datasetID}_cons_mpiled.vcf"), file("${datasetID}_consensus.fasta") into consensus_out
	 set datasetID, file("${inFastq}") into Kraken_db2
	 script:
     """
   	   samtools mpileup -uf ${datasetID}_ntref.fasta ${datasetID}.newsort.bam | bcftools call -c -o ${datasetID}_cons_mpiled.vcf
	   cat ${datasetID}_cons_mpiled.vcf  | vcfutils.pl vcf2fq > ${datasetID}_consensus.fastq
	   seqtk seq -a ${datasetID}_consensus.fastq > ${datasetID}_consensus.fasta
			   
	 """
}

}
if(params.kraken2_DB2 == true ){
process Classification_Kraken2_DB2 {
	 label 'NIMMseq'
	 
     publishDir "$outputDir/Kraken_parasite", mode: "copy", pattern: "*"
	 
	 tag {datasetID}
	 
	 maxForks 1
	 
	 input:
	 set datasetID, file(inFastq) from Kraken_db2
	 
	 output:
	 set datasetID,  file("${datasetID}_krona.html"), file("${datasetID}.report"), file("${datasetID}.kraken"), file("${datasetID}.classified"), file("${datasetID}.unclassified") into Kraken_fasta_parasite
	
	  
	 script:
     """
		
		kraken2 --db ${KRAKEN2_DB2} --threads ${threads} --fastq-input ${inFastq} > ${datasetID}.kraken --report ${datasetID}.report  --classified-out ${datasetID}.classified --unclassified-out ${datasetID}.unclassified  
		cat ${datasetID}.kraken | cut -f 2,3 >  ${datasetID}.krona.input 
		ktImportTaxonomy  ${datasetID}.krona.input  -o ${datasetID}_krona.html
	    
	 """
}

}

println "*********************************************************************"

/* 
 * Display information about the completed run
 */
 
workflow.onComplete {
    log.info "                                          "
    log.info "                                          "
    log.info "****************************************************************"
    log.info "Nextflow Version: $workflow.nextflow.version"
    log.info "Command line: $workflow.commandLine"
    log.info "Container:    $workflow.container"
    log.info "Duration: $workflow.duration"
	log.info "****************************************************************"
	
   }


