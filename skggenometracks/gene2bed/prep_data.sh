#! bin/bash
set -xeu

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gtf.gz

python prep_data.py gencode.v30.annotation.gtf.gz data/human.bed
python prep_data.py gencode.vM21.annotation.gtf.gz data/mouse.bed

cd data
gzip human.bed
gzip mouse.bed

cd ..
