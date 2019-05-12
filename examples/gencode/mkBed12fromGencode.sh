wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.basic.annotation.gtf.gz
gtfToGenePred gencode.v30.basic.annotation.gtf.gz gencode.v30.basic.annotation.gp -geneNameAsName2 -ignoreGroupsWithoutExons -genePredExt
awk '{ print $12"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"}' gencode.v30.basic.annotation.gp > gencode.v30.basic.annotation.gene_name.gp
genePredToBed gencode.v30.basic.annotation.gene_name.gp gencode.v30.basic.annotation.bed12
bedtools sort -i gencode.v30.basic.annotation.bed12 > gencode.v30.basic.annotation.sort.bed12
pigz gencode.v30.basic.annotation.sort.bed12
rm gencode.v30.basic.annotation.gp
rm gencode.v30.basic.annotation.gene_name.gp

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.basic.annotation.gtf.gz
gtfToGenePred gencode.vM21.basic.annotation.gtf.gz gencode.vM21.basic.annotation.gp -geneNameAsName2 -ignoreGroupsWithoutExons -genePredExt
awk '{ print $12"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"}' gencode.vM21.basic.annotation.gp > gencode.vM21.basic.annotation.gene_name.gp
genePredToBed gencode.vM21.basic.annotation.gene_name.gp gencode.vM21.basic.annotation.bed12
bedtools sort -i gencode.vM21.basic.annotation.bed12 > gencode.vM21.basic.annotation.sort.bed12
pigz gencode.vM21.basic.annotation.sort.bed12
rm gencode.vM21.basic.annotation.gp
rm gencode.vM21.basic.annotation.gene_name.gp
