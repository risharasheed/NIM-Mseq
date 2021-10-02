#! /usr/bin/env python3
import argparse
import os
import pysam
import subprocess

def checkFile(inFile):
    try:
        bamfile = pysam.AlignmentFile(inFile, "rb")
        bamfile.close()
        return 0
    except:
        return 1
def emptyFile(inFile, splID):
    fofn= "refAcc\treadsMapped\tsupplementaryReads\tbaseTotal\tbaseMapped\n"
    fofn_map= "readID\tref\tquality\treadLen\talnLen\tALNstart\tALNend\tsuppMapping\n"
    with open("{0}.mapping".format(splID), 'w') as outFA:
        outFA.write(fofn_map)
    with open("{0}.mapStats".format(splID), 'w') as outF:
        outF.write(fofn)
        
def reportMapReads(inFile, splID):
	bamfile = pysam.AlignmentFile(inFile, "rb")
#	newBamfile = pysam.AlignmentFile("{0}_new.bam".format(splID), "wb", template=bamfile)
	fofn= "refAcc\treadsMapped\tsupplementaryReads\tbaseTotal\tbaseMapped\n"
	fofn_map= "readID\tref\tquality\treadLen\talnLen\tALNstart\tALNend\tsuppMapping\n"
	
	refs = bamfile.references
	for ref in refs:
		readNum = 0
		baseMapped = 0
		baseTotal = 0
		supplementaryNum = 0
		
		for read in bamfile.fetch(reference=ref):
			read_baseMapped = read.query_alignment_length
			read_baseTotal = read.query_length
#			newBamfile.write(read)
			baseMapped = baseMapped + read.query_alignment_length
			baseTotal = baseTotal + read.query_length
			if read.is_supplementary:
				supplementaryNum = supplementaryNum + 1
				fofn_map += '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(read.query_name, ref, read.mapping_quality, read.query_length, read.query_alignment_length, read.reference_start, read.reference_end, "1")

			else:
				readNum = readNum + 1
				fofn_map += '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(read.query_name, ref, read.mapping_quality, read.query_length, read.query_alignment_length, read.reference_start, read.reference_end, "0")
	
		fofn += '{0}\t{1}\t{2}\t{3}\t{4}\n'.format(ref, readNum, supplementaryNum, baseTotal, baseMapped)
		
	bamfile.close()
#	newBamfile.close()
	
	with open("{0}.mapping".format(splID), 'w') as outFA:
		outFA.write(fofn_map)
		
	with open("{0}.mapStats".format(splID), 'w') as outF:
		outF.write(fofn)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='write aligned reads')
    
    parser.add_argument("-inFile", dest="inFile", required= True, help = "Sorted BAM file")
    parser.add_argument("-splID", dest="splID", required= True, help = "Sample ID")

    params = parser.parse_args()
    bamch = checkFile(inFile = params.inFile)
    if bamch == 0:
        reportMapReads(inFile = params.inFile, splID = params.splID)
            
    else:
        emptyFile(inFile = params.inFile, splID = params.splID)
