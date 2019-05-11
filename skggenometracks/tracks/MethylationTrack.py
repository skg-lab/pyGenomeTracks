# -*- coding: utf-8 -*-

from . BedGraphTrack import BedGraphTrack
from . GenomeTrack import GenomeTrack
import pandas as pd
import matplotlib.pyplot as plt

class MethylationTrack(GenomeTrack):
    # this track class extends a BedGraphTrack that is already part of
    # pyGenomeTracks. The advantage of extending this class is that
    # we can re-use the code for reading a bedgraph file
    SUPPORTED_ENDINGS = ['.mr']
    TRACK_TYPE = 'methylation_rate_graph'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + """
        # input methylation rate file which describe per 2000 window resolution methylation rate.
        # E.g.
        # chr start end methyl_sum de_methyl_sum per_methyl low ratio high
        # chrX    7265601 7267600 90      0       0.9891304347826086      0.9676158352635656      0.9924119338606987      0.9994364962522191
        # chrX    7266001 7268000 82      0       0.9880952380952381      0.9645504318541149      0.9916836033081886      0.9993821994187705
        half_window_step = 1000
        file_type = {}
        """.format(TRACK_TYPE)

    def __init__(self, properties_dict):
        self.properties = properties_dict

        if 'alpha' not in self.properties:
            self.properties['alpha'] = 1
        if 'base_color' not in self.properties:
            self.properties['base_color'] = '#FF0080'
        if 'fill_between_color' not in self.properties:
            self.properties['fill_between_color'] = 'blue'
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

        df = pd.read_csv(self.properties['file'], header=None, sep='\t')
        df.columns =  ["chr", "start", "end","methyl_sum","de_methyl_sum","per_methyl","low","ratio","high"]
        #startとendの中央をとる
        #half_window_step makes the 'center' of 'start' and 'end'
        df["center"] = df["end"]-int(self.properties['half_window_step'])

        #E.g. chrom_region=chrX
        df_q = df[df['chr'] == str(chrom_region)]
        df_s= df_q.query('@start_region <= center <= @end_region')

        ax.plot(df_s['center'], df_s['ratio'],
                c=self.properties['base_color'], alpha=self.properties['alpha'])
        ax.fill_between(df_s['center'],df_s["low"],
                df_s["high"], color=self.properties['fill_between_color'], alpha=0.2)

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
