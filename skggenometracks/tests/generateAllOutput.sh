# This bash script can be used to regenerate all test output when a global change happens.
# It needs to be launched from the root folder.
./bin/sgt --tracks ./skggenometracks/tests/test_data/bed_and_gtf_tracks.ini --region X:3000000-3300000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o ./skggenometracks/tests/test_data/master_bed_and_gtf.png
./bin/sgt --tracks ./skggenometracks/tests/test_data/bedgraph.ini --region X:2850000-3150000 --trackLabelFraction 0.2 --dpi 130 -o ./skggenometracks/tests/test_data/master_bedgraph.png
./bin/sgt --tracks ./skggenometracks/tests/test_data/bigwig.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./skggenometracks/tests/test_data/master_bigwig.png
./bin/sgt --tracks ./skggenometracks/tests/test_data/alpha.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./skggenometracks/tests/test_data/master_alpha.png
./bin/sgt --tracks ./skggenometracks/tests/test_data/epilogos.ini --region X:3100000-3150000 --trackLabelFraction 0.2 --dpi 130 -o ./skggenometracks/tests/test_data/master_epilogos.png
./bin/sgt --tracks ./skggenometracks/tests/test_data/browser_tracks_hic.ini --region X:2500000-3500000 --trackLabelFraction 0.23 --width 38 --dpi 130 -o ./skggenometracks/tests/test_data/master_plot_hic.png
./bin/sgt --tracks ./skggenometracks/tests/test_data/narrow_peak.ini --region X:2760000-2802000 --trackLabelFraction 0.2 --dpi 130 -o ./skggenometracks/tests/test_data/master_narrowPeak.png
./bin/sgt --tracks ./skggenometracks/tests/test_data/browser_tracks.ini --region X:3000000-3500000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./skggenometracks/tests/test_data/master_plot.png

