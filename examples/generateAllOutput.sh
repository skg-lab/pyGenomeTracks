# This bash script can be used to regenerate some examples output when a global change happens.
# It needs to be launched from the root folder.
sgt --tracks ./examples/bigwig_track.ini --region X:2,500,000-3,000,000 -o ./examples/bigwig.png
sgt --tracks ./examples/bigwig_with_genes.ini --region X:2,800,000-3,100,000 -o ./examples/bigwig_with_genes.png
sgt --tracks ./examples/bigwig_with_genes_and_vlines.ini --region X:2,800,000-3,100,000 -o ./examples/bigwig_with_genes_and_vlines.png
sgt --tracks ./examples/hic_track.ini -o hic_track.png --region chrX:2500000-3500000 -o ./examples/hic_track.png
sgt --tracks ./examples/epilogos_track.ini --region X:3100000-3150000 -o ./examples/epilogos_track.png
sgt --tracks ./examples/epilogos_track2.ini --region X:3100000-3150000 -o ./examples/epilogos_track2.png
sgt --tracks ./examples/bedgraph_matrix_lines.ini --region X:2000000-3500000 -o ./examples/bedgraph_matrix_lines.png

# The following examples require a modification of skggenometracks (adding a new track class)
# sgt --tracks ./examples/new_track.ini --region X:3000000-3200000 -o ./examples/new_track.png
# sgt --tracks ./examples/bedgraph_matrix.ini --region X:2000000-3500000 -o ./examples/bedgraph_matrix.png
