import sys
import os
import pandas as pd
import argparse

def initialCnt(initfile):
	initcnt={}

	initfil = open(initfile, 'r')
	lines = initfil.readlines()
	for line in lines:
		line = line.rstrip()
		Ctemp = line.split("\t")
		accID = Ctemp[1]
		
		inittext = Ctemp[0]
		initcnt.setdefault(accID, inittext)

	initfil.close()
	return initcnt

def blastfull(blastfl):
    blast={}

    blastfile = open(blastfl, 'r')
    lines = blastfile.readlines()
    for line in lines:
        line = line.rstrip()
        Mtemp = line.split("\t")

        accID = Mtemp[1]
        blasttext = Mtemp[0]
        blast.setdefault(accID, blasttext)

    blastfile.close()
    return blast

def blastout(blastoutfl):
    blastf={}

    blastof = open(blastoutfl, 'r')
    lines = blastof.readlines()
    for line in lines:
        line = line.rstrip()
        Btemp = line.split("\t")

        accID = Btemp[1]
        blastotext = Btemp[0]
        blastf.setdefault(accID, blastotext)

    blastof.close()
    return blastf
 
def humanCnt(humanf):
    human={}

    humanfl = open(humanf, 'r')
    lines = humanfl.readlines()
    for line in lines:
        line = line.rstrip()
        Htemp = line.split("\t")

        accID = Htemp[3]
        humantxt = "\t".join((Htemp[2],Htemp[0],Htemp[1]))
        human.setdefault(accID, humantxt)

    humanfl.close()
    return human 
def silvaCnt(silvaf):
    silva={}

    silvafl = open(silvaf, 'r')
    lines = silvafl.readlines()
    for line in lines:
        line = line.rstrip()
        Stemp = line.split("\t")

        accID = Stemp[3]
        silvatxt = "\t".join((Stemp[0],Stemp[1]))
        silva.setdefault(accID, silvatxt)

    silvafl.close()
    return silva
    
def krakenCnt(kraken):
    krakend={}

    krakenfl = open(kraken, 'r')
    lines = krakenfl.readlines()
    for line in lines:
        line = line.rstrip()
        Ktemp = line.split("\t")
         
        accID = Ktemp[2]
        krakentxt = "\t".join((Ktemp[1],Ktemp[0]))
        krakend.setdefault(accID, krakentxt)

    krakenfl.close()
    
    return krakend
def combineReport(inRef, krakend, silva, human,blastf,blast,initcnt):
    fofn_ref=''
    fofn=''
	
    headline = "\t".join(("SampleID", "fastq count", "Nanofilt count", "Human unmapped","Human mapped", "Silva mapped", "Silva unmapped","Blast result count", "Blast count after filter","Kraken Classified","Kraken Unclassified"))
    fofn+='{0}\n'.format(headline)
    refFile = open(inRef, 'r')
    lines = refFile.readlines()
    
    for Rline in lines:
        Rline = Rline.rstrip()
        temp = Rline.split("\t")
        RaccID = temp[0]
        
        try:
            incnt = initcnt[RaccID]
            human_cnt = human[RaccID]
            silva_cnt = silva[RaccID]
            blastfcnt = blastf[RaccID]
            blastocnt = blast[RaccID]
            kraken_cnt = krakend[RaccID]
            fofn+='{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(RaccID, incnt, human_cnt, silva_cnt,blastfcnt,blastocnt,kraken_cnt)
            

        except KeyError:
            try: 
                incnt = initcnt[RaccID]
                human_cnt = human[RaccID]
                silva_cnt = silva[RaccID]
                blastfcnt = blastf[RaccID]
                blastocnt = blast[RaccID]
                kraken_cnt = krakend[RaccID]
                fofn+='{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(RaccID, incnt, human_cnt, silva_cnt,blastfcnt,blastocnt,kraken_cnt)
                 
            except KeyError:
                print("key error ",RaccID)
                print(initcnt)
                pass
			
    return fofn

def writeFofn(fofn):
    with open("countReport",'wt') as outhandle:
        outhandle.write(fofn)
        
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Produce report for identified viral species')
    parser.add_argument("-init", dest="init", required= True, help = "Ident stats")
    parser.add_argument("-kra", dest="kra", required= True, help = "kraken report")
    parser.add_argument("-silva", dest="silva", required= True, help = "silva count")
    parser.add_argument("-human", dest="human", required= True, help = "human count")
    parser.add_argument("-blast", dest="blast", required= True, help = "blast count")
    parser.add_argument("-blastf", dest="blastf", required= True, help = "blast full")
    parser.add_argument("-inref", dest="inref", required= True, help = "ref file")
    parser.add_argument("-outDir", dest="outDir", required= True, help = "output folder")
    params= parser.parse_args()
    path=params.outDir
    os.chdir(path)
    krakend= krakenCnt(params.kra)
    silva= silvaCnt(params.silva)
    human=humanCnt(params.human)
    blastf=blastfull(params.blastf)
    blast=blastout(params.blast)
    initcnt=initialCnt(params.init)
    fofn= combineReport(params.inref, krakend, silva, human,blastf,blast,initcnt)
    writeFofn(fofn)

    
