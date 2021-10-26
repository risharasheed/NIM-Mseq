import sys
import os
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

def drawReportChart(inFile,prefix):
    with open(inFile) as reads:                                                                                          
    	dataset = pd.read_csv(reads, delimiter='\t')
    fastq_cnt    = pd.DataFrame({'total':dataset['fastq count'],'label':'Total Read','barcode':dataset['SampleID']})
    nan_filt     = pd.DataFrame({'total':dataset['Nanofilt count'],'label':'Nanofilt count','barcode':dataset['SampleID']})
    human_un     = pd.DataFrame({'total':dataset['Human_unmap'],'label':'Human unmapped'})
    human_map    = pd.DataFrame({'total':dataset['Human mapped'],'label':'Human mapped','barcode':dataset['SampleID']})
    silva_map    = pd.DataFrame({'total':dataset['Silva mapped'],'label':'Silva mapped','barcode':dataset['SampleID']})
    silva_unmap  = pd.DataFrame({'total':dataset['Silva unmapped'],'label':'Silva unmapped','barcode':dataset['SampleID']} )
    kraken_class = pd.DataFrame({'total':dataset['Kraken Classified'],'label':'Kraken Classified','barcode':dataset['SampleID']})
    krakne_uncla = pd.DataFrame({'total':dataset['Kraken Unclassified'],'label':'Kraken Unclassified','barcode':dataset['SampleID']})

    df = pd.concat([fastq_cnt,nan_filt,human_un,human_map,silva_map,silva_unmap,kraken_class,krakne_uncla])
    results_path = prefix
    results_path+=".png"
    g = sns.catplot(x="label",y="total", col="barcode", col_wrap=4,
                data=df,
                kind="bar", height=3, aspect=1,sharey=False)
    (g.set_axis_labels("", "Read Count")
      .set_xticklabels(rotation=90)
      .set_titles("{col_name}")
      .despine(left=True))
    plt.savefig(results_path, dpi=700)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Produce Heatmap for viral species')
    parser.add_argument("-inFile", dest="inFile", required= True, help = "in put file name")
    parser.add_argument("-outDir", dest="outDir", required= True, help = "output folder")
    parser.add_argument("-prefix", dest="prefix", required= True, help = "prefix ")
    
    params= parser.parse_args()
    path=params.outDir
    os.chdir(path)
    drawReportChart(params.inFile,params.prefix)
