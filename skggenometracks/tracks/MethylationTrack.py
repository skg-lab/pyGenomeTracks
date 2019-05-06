# -*- coding: utf-8 -*-
from . GenomeTrack import GenomeTrack

class MethylationTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.mr']  #mr=methylation rate
    TRACK_TYPE = 'methylation_rate_graph'
    OPTIONS_TXT = """
height = 3
title =
text =
# x position of text in the plot (in bp)
x position =
"""

    import pandas as pd
    import numpy as np

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import pylab
    #%matplotlib inline

    input_df=pd.read_table(file,header=None)
    #test_dfにcolumnsを追加
    input_df.columns = ["chr", "start", "end","methyl_sum","de_methyl_sum","per_methyl","low","ratio","high"]

    #startとendの中央をとる
    input_df["end"]=input_df["end"]-1000

    chr_input_df=input_df[input_df['chr'] == chrom]
    range_chr_df=chr_input_df[(chr_input_df['end'] >= region_start) & (chr_input_df['end'] <= region_end)]


    pylab.figure(figsize=(10, 4), dpi=200)
    plt.plot(range_chr_df["ratio"])
    plt.fill_between(range_chr_df.index,range_chr_df["low"], range_chr_df["high"], color="blue", alpha=0.2)

    # X軸の数字をオフセットを使わずに表現する
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)

    #グラフに凡例をつける
    plt.legend(["ratio"])
    plt.ylim([0, 1.0])
