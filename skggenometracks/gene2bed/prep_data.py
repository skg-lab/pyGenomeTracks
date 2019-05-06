import pandas as pd
import sys

f_gtf = sys.argv[1]
f_out = sys.argv[2]

df_gtf = pd.read_csv(f_gtf, sep='\t', comment='#', header=None)


def get_gene_name(s):
    l = [x.split() for x in s.split(';')[:-1]]
    return [x[1].replace('"', '') for x in l if x[0] == 'gene_name'][0]


df_gtf['gene_name'] = df_gtf.loc[:,8].map(get_gene_name)

l = []
for g in df_gtf.gene_name.unique():
#     print(g)
    df_q = df_gtf[df_gtf.gene_name==g]
    chrom = list(df_q[0])[0]
    start = df_q[3].min()
    end = df_q[4].max()
    strand = list(df_q[6])[0]
    l.append([chrom, start, end, g, 0, strand])
    
df_bed = pd.DataFrame(l)
df_bed.to_csv(f_out, header=None, index=None, sep='\t')
