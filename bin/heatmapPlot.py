import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import argparse

def drawHeatmap(inFile,prefix):
    with open(inFile) as reads:                                                                                          
    	df = pd.read_csv(reads, delimiter='\t',names=["Scientific Name","Barcode","% Refseq"])
    df = df.pivot("Scientific Name","Barcode","% Refseq")
    results_path = prefix
    results_path+=".png"
#    sb.set(font_scale=2)
    fig, ax = plt.subplots(figsize=(32, 20))
    sb.heatmap(df,cmap="YlGnBu",linewidths=.01,linecolor='gray')
#    sb.heatmap(df,cmap="YlGnBu")
    plt.yticks(rotation=20)
#    plt.title("", fontsize=30)
    plt.savefig(results_path, dpi=700)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Produce Heatmap for viral species')
    parser.add_argument("-inFile", dest="inFile", required= True, help = "in put file name")
    parser.add_argument("-outDir", dest="outDir", required= True, help = "output folder")
    parser.add_argument("-prefix", dest="prefix", required= True, help = "prefix ")
    
    params= parser.parse_args()
    path=params.outDir
    os.chdir(path)
    drawHeatmap(params.inFile,params.prefix)