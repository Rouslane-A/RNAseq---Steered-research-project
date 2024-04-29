import pandas as pd

featurecounts_df = pd.read_csv('output_data.csv', index_col=0)

sum_reads = featurecounts_df.sum()


features_100k = featurecounts_df.loc[:, sum_reads >= 100000]

num_genes = (features_100k >0).sum()

min_genes = 2346
max_genes = 9884

filtered_df = features_100k.loc[:,(num_genes >= min_genes) & (num_genes <= max_genes)]

filtered_df.to_csv('filtered_counts.csv')

