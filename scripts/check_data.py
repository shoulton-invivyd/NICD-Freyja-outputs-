import pandas as pd

agg_results = '../agg_demixed.tsv'
meta = '../sample_metadata.csv'

agg_df = pd.read_csv(agg_results, skipinitialspace=True, sep='\t', index_col=0)#.drop_duplicates()
times_df = pd.read_csv(meta, skipinitialspace=True).drop_duplicates().reset_index(drop=True)

#check for empty entries

print(agg_df[agg_df.duplicated(keep='first')].shape[0])

if agg_df[agg_df.duplicated(keep='first')].shape[0] >0:
    agg_df = agg_df[~agg_df.duplicated(keep='first')]
    agg_df.to_csv(agg_results,sep='\t')


# times_df[times_df.isnull().any(axis=1)]

#check for  duplicated seq ids
print(times_df[times_df['Sequence_ID'].duplicated(keep='first')].groupby('Sequence_ID').agg({"SiteName": lambda x: x.nunique()}).shape[0])

#check for duplicated lab numbers
print(times_df[times_df['LabNumber'].duplicated(keep=False)].dropna())


print(times_df[times_df['Sequence_ID'].duplicated(keep='first')].shape[0])
print(times_df[times_df['LabNumber'].duplicated(keep=False)].dropna())
# agg_df.to_csv(agg_results,sep='\t')
times_df.to_csv(meta)

