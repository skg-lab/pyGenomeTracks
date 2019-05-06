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
    """.format(TRACK_TYPE)

    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)

        import pandas as pd
        import numpy as np

        import matplotlib as mpl
        import matplotlib.pyplot as plt
        import pylab

        #from skggenometracks.plotTracks import get_region(region_string)

        #self.bw = pyBigWig.open(self.properties['file'])
        input_df=pd.read_table(self.properties['file'],header=None)
        #test_dfにcolumnsを追加
        input_df.columns = ["chr", "start", "end","methyl_sum","de_methyl_sum","per_methyl","low","ratio","high"]

    def plot(self, ax, chrom_region, start_region, end_region):
        formated_region = "{}:{}-{}".format(chrom_region, start_region, end_region)

        import pandas as pd
        import numpy as np

        import matplotlib as mpl
        import matplotlib.pyplot as plt
        import pylab

        input_df=pd.read_table(self.properties['file'],header=None)
        input_df.columns = ["chr", "start", "end","methyl_sum","de_methyl_sum","per_methyl","low","ratio","high"]
        #startとendの中央をとる
        input_df["end"]=input_df["end"]-1000

        chr_input_df=input_df[input_df['chr'] == chrom_region]
        range_chr_df=chr_input_df[(chr_input_df['end'] >= start_region) & (chr_input_df['end'] <= end_region)]


        pylab.figure(figsize=(10, 4), dpi=200)
        ax.plot(range_chr_df["ratio"])
        ax.fill_between(range_chr_df.index,range_chr_df["low"], range_chr_df["high"], color="blue", alpha=0.2)

        # X軸の数字をオフセットを使わずに表現する
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)

        #グラフに凡例をつける
        ax.legend(["ratio"])
        plt.ylim([0, 1.0])

        return ax
        #plt.show()
