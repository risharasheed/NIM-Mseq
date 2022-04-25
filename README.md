## NIM-MSeq (Nanopore Informatics for Metagenomics of Microbial Sequences)
The NIM-Mseq is a NextFlow based, open-source bioinformatics pipeline that enables detection of microbes from basecalled Oxford Nanopore data. This pipeline is academically based on opensource tools and bioinformatic algorithms and is optimized for scalable multi thread processor in any cloud platform or standalone processor desktops. This pipeline enables parallel processing of multiple barcoded datasets. All the software used for this pipeline are packed in conda environment that can be executed on wide range of computing systems. 

The pipeline and its dependent softwares can be easily installed from GitHub and conda. The user can straightforwardly run the pipeline by updating the configuration file and input dataset. This pipeline enables the user to monitor the progress of each process. The pipeline output displays the detailed logs like usage of CPU, memory and process informations.  Once the pipeline completes the analysis, report will be generated for each dataset(barcode). Minimum system requirement is 16 CPU and 64+ GB ram to run the pipeline. 
 
## Build environment 

This pipeline uses the nextflow manager and conda environment, please install nextflow and conda.
* Install nextflow from [Nextflow Documentation](https://www.nextflow.io/docs/latest/getstarted.html)
* Install conda from [conda installation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
## Get NIM-Mseq code from GitHub
* Create a working directory where you want to save the code.
* Clone the code from GitHub to working directory
 ```
git clone https://github.com/risharasheed/NIM-Mseq.git
```

## Build Conda environment 
* Use the [link](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) to build the conda environment and install the all packages based on the dependent list from the software/tool table below, or follow the below steps to build the conda environment from file.  
- Use the terminal or an Anaconda Prompt for the following steps:
     - Copy the environmnet **environment.yml**  file to your working directory
     - Create the environment from the **environment.yml** file
     -  ```   conda env create -f environment.yml   ```
     -  Activate the new environment: 
     -  ``` conda activate NIM-Mseq ```
     -  Verify that the new environment was installed correctly:  
     -  ``` conda env list  ``` 
     -   it will show the newly created **NIM-Mseq** environment 
## Build databases for NIM-Mseq
* Human - GRCh37 human genome/proteins 
     * Download the GRCh37 genome data from NCBI
     * unzip the fasta file
     * Build the BWA index using BWA command 
     * ``` bwa index -a bwtsw  hg19_reference/hg19_exome.fa  ```   
     * hg19_reference/hg19_exome.fa  is your downloaded GRCh37.
     * host enviroment databse, in this process we used human databse.
* Silva DB
     * Download the silva db from [here](https://www.arb-silva.de/download/arb-files/)
     * unzip the fasta file
     * Build the BWA index using BWA command
     *  ``` bwa index -a bwtsw  silvadb/silvadb.fa   ``` 
     * silvadb/silvadb.fa  is your downloaded silva db.      
* kraken2 (archaea, bacteria, viral, plasmid, human, UniVec_Core, protozoa & fungi)
     * Build the kraken database based on your required pathogen from NCBI using kraken2-build command, refer [here](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown) to build DB. It required more RAM and CPU to build the database.
     * Other option is to use the existing database which is already built on the cloud. Please refer [this link](https://benlangmead.github.io/aws-indexes/k2) to get the kraken database which already built for different pathogen.
     * in this pipeline we have used two kraken step from two datbases whihc is copied from above link, you can make a single database and set false in config file to skip  one kraken process.
* BLAST (similarity search blastn) DB
     * Download the required blast data from [NCBI FTP site.](https://ftp.ncbi.nlm.nih.gov/blast/db/)
     * Download the reference fsata file from [NCBI](https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/)
## nextflow configuration
The nextflow config file will be used to setup the variable used in nextflow to run the pipeline.
open the netflow.config file and update the below variable based on your directory/data.
```
inputDir = "/../run1/fastq_combined/"     - Directory contains input sequencing reads. Specify the directory to read all reads or specify the full path for specific reads   
outputDir = "/../output/nextflow/run1"    - Specify the directory in which the output need to be stored
threads=30                                - Number of parallel process or number of processors
human_g1k="../human/human_g1k_v37.fasta”  - Specify the directory in which contains human genome index build using BWA.
silva_db="../silva_db/SILVA.fasta"        - Specify the directory in which contains silva fasta file 
KRAKEN2_DB="../kraken2_db"                - Specify the directory in which contains kraken2 database
BLAST_db_nrnt="../blast_db/nt/nt"         - Specify the directory in which contains the blastn database
BLAST_Fasta_nt="../blast_fasta/nt"        - Specify the directory in which contains blast reference fasta file
projDir = "../nextflow/bin"               - This is the directory in which the bin folder is present after git clone.
BLASTDB="$BLASTDB:/../blast_db/nt/"       - This is same as blast database except the last nt to get details for taxomic information from blast
```
## Run the pipeline
Once the next flow configure is completed, use the below command to run the pipeline. to bypass any process from the pipeline, use the config file to set *false* for that process. In this pipeline we have used Porechop as the demultiplexer/barcode.  Likewise based on the requirement alternate demultiplexer can be integrated to the pipeline and used. In our study we have used both Porechop and Guppy.   To integrate Guppy please use [nanopore community]( https://nanoporetech.com/) and [download](  https://community.nanoporetech.com/downloads) the software.  
```
nextflow run NIM-Mseq.nf -with-report
```
once the process is completed you will see the output like below.

![Image of run output](https://github.com/risharasheed/NIM-Mseq/blob/main/images/NIM-Mseq-final-stscreen.png)


## Software/tools used in NIM-Mseq
|Software	|Version or higher	|Link
|---------|  ------|----
|Porechop	|5.0.11	|[porechop](https://github.com/rrwick/Porechop)
|NanoPlot	|1.32.1	|[nanoplot](https://github.com/wdecoster/NanoPlot)
|nanofilt	|2.8	|[nanofilt](https://github.com/wdecoster/nanofilt)
|bwa	|0.7.17	|[BWA](https://github.com/lh3/bwa)
|samtools	|1.7	|[samtools](http://www.htslib.org/)
|kraken2	|2.0.7	|[kraken](https://github.com/DerrickWood/kraken2)
|seqtk	|1.3	|[seqtk](https://github.com/lh3/seqtk)
|Blast	|2.10.1	|[blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
|minimap2	|2.21	|[minimap](https://github.com/lh3/minimap2)
|Krona	|2.8	|[krona](https://github.com/marbl/Krona/wiki)
|bcftools|	1.9	|[bcftools](http://samtools.github.io/bcftools/bcftools.html)
|nextflow	|30	|[nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
|Python	|3.7	|[python](https://www.python.org/)
|R	|3.5	|[R](https://www.r-project.org/)


#### Common errors  
 * There is a possibility to fail quality process due to python seaborn package issue. Hence please verify the seaborn package is 0.10.1. if it uses latest one then replace the latest one with 0.10.1 in specific conda environment python library.
 * The krona report steps may give error like #Loading taxonomy...  #Taxonomy not found in /home/user/. miniconda3/envs/NIM-Mseq/opt/krona/taxonomy.  In this case, please  run the command updaeTaxonomy.sh   from your krona directory. 
```
bash /home/user/anoconda3/envs/NIM-Mseq/opt/krona/updateTaxonomy.sh  
```
* #If above command is failed with make not found error, then do below command and then rerun.
```|
sudo apt-get install build-essential
```
* Then run the updateTaxonomy.sh again.
* Any step is failing without proper error then try to increase the RAM and rerun. 
* The pipeline is setup to run one file at a time in a single process due to space constrain, if you have more RAM then change the value of maxForks in nextflow scripts  


#### References
1.	De Coster, W. et al. (2018) ‘NanoPack: Visualizing and processing long-read sequencing data’, Bioinformatics, 34(15), pp. 2666–2669. doi: 10.1093/bioinformatics/bty149.
2.	Li, H. (2013) ‘Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM’, 00(00), pp. 1–3. Available at: http://arxiv.org/abs/1303.3997.
3.	Li, H. (2018) ‘Minimap2: Pairwise alignment for nucleotide sequences’, Bioinformatics, 34(18), pp. 3094–3100. doi: 10.1093/bioinformatics/bty191.
4.	Quast, C. et al. (2013) ‘The SILVA ribosomal RNA gene database project: Improved data processing and web-based tools’, Nucleic Acids Research, 41(D1), pp. 590–596. doi: 10.1093/nar/gks1219.
5.	Sayers, E. W. et al. (2021). Database resources of the National Center for Biotechnology Information. Nucleic acids research, 49(D1), D10–D17. https://doi.org/10.1093/nar/gkaa892
6.	Wood, D. E., Lu, J. and Langmead, B. (2019) ‘Improved metagenomic analysis with Kraken 2’, bioRxiv. Genome Biology, pp. 1–13. doi: 10.1101/762302.

