#! /usr/bin/env python3
import sys
from Bio import SeqIO
import gzip
import os
import subprocess
import argparse
import csv


def dictIdnt(inIdent):
	dictIdent={}

	identFile = open(inIdent, 'r')
	identFile.readline()
	for line in identFile:
		line = line.rstrip()
		Ctemp = line.split("\t")
		accID = Ctemp[0]
		
		identText = "\t".join((Ctemp[1],Ctemp[2])) 
		dictIdent.setdefault(accID, identText)

	identFile.close()
	return dictIdent

def dictM(inMap):
    dictMap={}

    mapFile = open(inMap, 'r')
    mapFile.readline()
    for line in mapFile:
        line = line.rstrip()
        Mtemp = line.split("\t")

        accID = Mtemp[0]
        mapText = "\t".join((Mtemp[1], Mtemp[4]))
        dictMap.setdefault(accID, mapText)

    mapFile.close()
    return dictMap


def combineReport(inRef, dictIdent, dictMap, splID):
    fofn_ref=''
    fofn=''
    headline = "\t".join(("Barcode",	"Accession ID", "Scientific Name", "Per.Identity", "% Ref-Seq", "Mapped reads", "Total Mapped Bases"))
    fofn+='{0}\n'.format(headline)
    refFile = open(inRef, 'r')
    lines = refFile.readlines()
    for Rline in lines:
        Rline = Rline.rstrip()
        temp = Rline.split("\t")
        RaccID = temp[0]
        Rname = temp[2]
        try:
            Rident = dictIdent[RaccID]
            Rmap = dictMap[RaccID]
            fofn+='{0}\t{1}\t{2}\t{3}\t{4}\n'.format(splID, RaccID, Rname, Rident, Rmap)
            Rref = Rident.rstrip()
            temp = Rref.split("\t")
            Rrefseq = temp[1]
            fofn_ref+='{0}\t{1}\t{2}\n'.format(Rname,splID,Rrefseq)
        except KeyError:
            try: 
                Rident = dictIdent[RaccID]
                Rmap = dictMap[RaccID]
                fofn+='{0}\t{1}\t{2}\t{3}\t{4}\n'.format(splID, RaccID, Rname, Rident, Rmap)
                Rref = Rident.rstrip()
                temp = Rref.split("\t")
                Rrefseq = temp[1]
                fofn_ref+='{0}\t{1}\t{2}\n'.format(Rname,splID,Rrefseq)
            except KeyError:
                pass
			
    return fofn, fofn_ref

def writeFofn(splID, fofn, fofn_ref):
    with open("{0}_r.pathogenReport".format(splID),'wt') as outhandle:
        outhandle.write(fofn)
    with open("{0}_r.heatmapM".format(splID),'wt') as outhandle:
        outhandle.write(fofn_ref)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Produce report for identified viral species')
    parser.add_argument("-inIdent", dest="inIdent", required= True, help = "Ident stats")
    parser.add_argument("-inMap", dest="inMap", required= True, help = "Mapping stats")
    parser.add_argument("-inRef", dest="inRef", required= True, help = "Reference sequences")
    parser.add_argument("-splID", dest="splID", required= True, help = "sample ID")
    

    params= parser.parse_args()

    dictIdent= dictIdnt(params.inIdent)
    dictMap= dictM(params.inMap)
  

    fofn, fofn_ref= combineReport(params.inRef, dictIdent, dictMap, params.splID)
    writeFofn(params.splID, fofn, fofn_ref)
