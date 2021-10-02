import pandas as pd
import argparse


def combineReport(inRef, splID):
    fofn=''
    fofn_idr=''
    headline = "\t".join(("Query ID", "Accession ID", "Query Length", "Accession Length", "Matching Bases", "Number Of Bases", "Mapping Quality",  "% Ref-Seq", "% Ref-Seq / Read", "Per.Identity"))
    fofn+='{0}\n'.format(headline)
    fofn_id=''
    headline_idn = "\t".join(("ref","refseq","identity"))
    fofn_id+='{0}\n'.format(headline_idn)
    fofn_idr+='{0}\n'.format(headline_idn)
    refFile = open(inRef, 'r')
    lines = refFile.readlines()
    for Rline in lines:
        Rline = Rline.rstrip()
        temp = Rline.split("\t")
        Qname = temp[0]
        Qlength = temp[1]
        RefName = temp[2]
        RefLength = temp[3]
        RefStart = temp[4]
        RefEnd = temp[5]
        Match = temp[6]
        MissGap = temp[7]
        Quality = temp[8]
		
        PrefSeq = round(100 * ((int(RefEnd) - int(RefStart)) / int(RefLength)),2)
        RefSeqP = round(100 * ((int(RefEnd) - int(RefStart)) / int(Qlength)),2)
        Pidentity = round(100 * (int(Match) / int(MissGap)),2)
        
        fofn+='{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(Qname, RefName, Qlength, RefLength, Match,MissGap,Quality,PrefSeq,RefSeqP,Pidentity)
        if float(PrefSeq) > 1.0:
            fofn_id+='{0}\t{1}\t{2}\n'.format(RefName,PrefSeq,Pidentity)
        else:
            fofn_idr+='{0}\t{1}\t{2}\n'.format(RefName,PrefSeq,Pidentity)
            
    return fofn, fofn_id,fofn_idr

def writeFofn(splID, fofn, fofn_id,fofn_idr):
    with open("{0}_r.detailed.Report.txt".format(splID),'wt') as outhandle:
        outhandle.write(fofn)
    with open("{0}_r.ident_reseq.txt".format(splID),'wt') as outhandleSummary:
        outhandleSummary.write(fofn_id)
    with open("{0}_r.ident_removed.txt".format(splID),'wt') as outhandleSummary:
        outhandleSummary.write(fofn_idr)
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Produce report for identified viral species')
    parser.add_argument("-inRef", dest="inRef", required= True, help = "INput file stats")
    parser.add_argument("-splID", dest="splID", required= True, help = "sample ID")
    
    params= parser.parse_args()

    fofn, fofn_id, fofn_idr= combineReport(params.inRef, params.splID)
    writeFofn(params.splID, fofn, fofn_id, fofn_idr)
    