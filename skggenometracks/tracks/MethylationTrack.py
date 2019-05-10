# -*- coding: utf-8 -*-

from . BedGraphTrack import BedGraphTrack
from . GenomeTrack import GenomeTrack
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class MethylationTrack(GenomeTrack):
    # this track class extends a BedGraphTrack that is already part of
    # pyGenomeTracks. The advantage of extending this class is that
    # we can re-use the code for reading a bedgraph file
    SUPPORTED_ENDINGS = ['mr']
    TRACK_TYPE = 'methylation_rate_graph'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + """
        # bismark file which describe 1bp resolution methylation rate.
        # E.g.
        # chr1    10563   10563   90.9090909090909        10      1
        # chr1    10571   10571   83.3333333333333        10      2
        file_type = {}
        """.format(TRACK_TYPE)

    def __init__(self, properties_dict):
        self.properties = properties_dict

        if 'alpha' not in self.properties:
            self.properties['alpha'] = 1
        if 'color' not in self.properties:
            self.properties['color'] = '#FF0080'
        if 'size' not in self.properties:
            self.properties['size'] = 10


    def plot(self, ax, chrom_region, start_region, end_region):
        """
        Args:
           ax: matplotlib axis to plot
           chrom_region: chromosome name
           start_region: start coordinate of genomic position
           end_region: end coordinate
        """

        # the BedGraphTrack already has methods to read files
        # in which the first three columns are chrom, start,end
        # here we used the get_scores method inherited from the
        # BedGraphTrack class

        # print(chrom_region, start_region, end_region)

        df = pd.read_csv(self.properties['file'], header=None, sep='\t')
        df.columns =  ["chr", "start", "end","methyl_sum","de_methyl_sum","per_methyl","low","ratio","high"]
        #startとendの中央をとる
        df["center"] = df["end"]-1000

        df_q = df[df['chr'] == 'chr'+str(chrom_region)]
        df_s= df_q.query('@start_region <= center & center <= @end_region ')
        #fig = plt.figure()
        #ax = fig.add_subplot(1,1,1)
        ax.plot(df_s['ratio'])
        #ax.plot(df_s['ratio'],alpha=self.properties['alpha'],
        #        color=self.properties['color'])
        ax.set_ylim(-5,105)

    def plot_y_axis(self, ax, plot_axis):
        """
        Plot the scale of the y axis with respect to the plot_axis
        Args:
            ax: axis to use to plot the scale
            plot_axis: the reference axis to get the max and min.
        Returns:
        """
        if 'show data range' in self.properties and self.properties['show data range'] == 'no':
            return

        def value_to_str(value):
            # given a numeric value, returns a
            # string that removes unneeded decimal places
            if value % 1 == 0:
                str_value = str(int(value))
            else:
                str_value = "{:.1f}".format(value)
            return str_value

        ymin, ymax = 0, 1

        ymax_str = value_to_str(ymax)
        ymin_str = value_to_str(ymin)
        # plot something that looks like this:
        # ymax ┐
        #      │
        #      │
        # ymin ┘

        # the coordinate system used is the ax.transAxes (lower left corner (0,0), upper right corner (1,1)
        # this way is easier to adjust the positions such that the lines are plotted complete
        # and not only half of the width of the line.
        x_pos = [0, 0.5, 0.5, 0]
        y_pos = [0.01, 0.01, 0.99, 0.99]
        ax.plot(x_pos, y_pos, color='black', linewidth=1, transform=ax.transAxes)
        ax.text(-0.2, -0.01, ymin_str, verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes)
        ax.text(-0.2, 1, ymax_str, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes)
        ax.patch.set_visible(False)
