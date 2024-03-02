import pandas as pd
import pickle
import json
import requests
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.dates as mdates
from scipy import signal
from datetime import date,timedelta
import yaml
import copy
import numpy as np 

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

with open( "plot_config.yml", "r" ) as f :
    plot_config = yaml.load( f, Loader=yaml.FullLoader )

#--- borrowed from SEARCH wastewater surveillance dashboard, coordinated by Scripps Research.--#
def convert_rbg_to_tuple( rgb ):
    rgb = rgb.lstrip( "#" )
    return tuple( int( rgb[i :i + 2], 16 ) for i in (0, 2, 4) )
def convert_tuple_to_rgb( r, g, b ):
    return '#%02x%02x%02x' % (int(r), int(g), int(b))
def lighten_field( value, alpha, gamma=2.2 ):
    return pow( pow(255, gamma) * (1 - alpha) + pow( value, gamma ) * alpha, 1 / gamma)
def lighten_color( r, g, b, alpha, gamma=2.2 ):
    return lighten_field(r, alpha, gamma ), lighten_field( g, alpha, gamma ), lighten_field( b, alpha, gamma )

children_dict = dict()
delta = 0.15
for key in reversed( list( plot_config.keys() ) ):
    for value in ["name", "members"]:
        assert value in plot_config[key], f"YAML entry {key} is not complete. Does not contain '{value}' entry."
    if "color" not in plot_config[key]:
        assert "parent" in plot_config[key], f"YAML entry {key} is incomplete. Must specify either a 'color' or 'parent' entry."
        if plot_config[key]["parent"] in children_dict:
            children_dict[plot_config[key]["parent"]] += 1
        else:
            children_dict[plot_config[key]["parent"]] = 1
        child_idx = children_dict[plot_config[key]["parent"]]
        parent_color = plot_config[plot_config[key]["parent"]]["color"]
        parent_color = convert_rbg_to_tuple( parent_color )
        plot_config[key]["color"] = convert_tuple_to_rgb( *lighten_color( *parent_color, alpha=1.0-(delta*child_idx) ) )

#----#

with open('lineages.yml', 'r') as f:
        try:
            lineages_yml = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            raise ValueError('Error in lineages.yml file: ' + str(exc))

lineage_info = {}
for lineage in lineages_yml:
    lineage_info[lineage['name']] = {'children': lineage['children']}

#Load json file

# agg_df = pd.read_csv('agg-masters.tsv', skipinitialspace=True, sep='\t',index_col=0) 
agg_df = pd.read_csv('../agg_demixed.tsv', skipinitialspace=True, sep='\t',index_col=0) 
agg_df = agg_df[agg_df['coverage']>50.0] 


from freyja.utils import prepLineageDict, prepSummaryDict
agg_df = prepSummaryDict(agg_df)
agg_df = prepLineageDict(agg_df,thresh=0.0000000001,config=plot_config,lineage_info=lineage_info)

#group other recombinants, move everything else to "Other"
groupNames = set(list(plot_config.keys()))
processed_linDictMod = []
for j, sampLabel in enumerate(agg_df.index):
    linDictMod = copy.deepcopy(agg_df.loc[sampLabel, 'linDict'])
    linDictModCopy = copy.deepcopy(agg_df.loc[sampLabel, 'linDict'])
    for rInd in linDictModCopy.keys():
        if rInd.startswith('X') and (rInd not in groupNames):
            if "Recombinants" in linDictMod.keys():
                linDictMod["Recombinants"] += linDictMod.pop(rInd)
            else:
                linDictMod["Recombinants"] = linDictMod.pop(rInd)
        elif (rInd not in groupNames):
            if "Other" in linDictMod.keys():
                linDictMod["Other"] += linDictMod.pop(rInd)
            else:
                linDictMod["Other"] = linDictMod.pop(rInd)
    processed_linDictMod.append(linDictMod)
agg_df.loc[:, 'linDict'] = processed_linDictMod
# move everything else to "Other"
# agg_df.index = [adi.replace('_','-') for adi in agg_df.index]

agg_df.index = [adi.replace('_','-').replace('ENV-','').replace('ENC-','').replace('.tsv','') for adi in agg_df.index]
agg_df.index  = ['-'.join(adi.split('-')[0:3]) if (adi[0:3]=='NIC' or adi[0:3]=='COV') else '-'.join(adi.split('-')[0:2]) for adi in agg_df.index ]

agg_df = agg_df[['linDict']]

times_df = pd.read_csv('../sample_metadata.csv', skipinitialspace=True)
times_df['LabNumber'] = times_df['LabNumber'].apply(lambda x:x.replace('ENV-',''))
times_df = times_df.set_index('LabNumber')
times_df['SampleCollectionDate'] = pd.to_datetime(times_df['SampleCollectionDate'])
#if exact duplicate, drop
times_df = times_df.loc[~times_df.duplicated(keep='first')]


merged_df = times_df.merge(agg_df,left_index=True,right_index=True)
merged_df['Lineages']= merged_df['linDict'].apply(lambda x: ' '.join(x.keys()))
merged_df['Abundances']= merged_df['linDict'].apply(lambda x: list(x.values()))
merged_df = merged_df.drop(columns=['linDict'])
merged_df.to_csv('merged_data.tsv',sep='\t')

# merged_df.to_csv()
# asdf
# meta_df = pd.read_excel('Metadata.xlsx')
# meta_df['SampleCollectionDate'] = pd.to_datetime(meta_df['SampleCollectionDate'])
# # meta_df = meta_df[meta_df['SampleCollectionDate']>'2022-01-01']
# meta_df['LabNumber'] = [adi.replace('ENV-','') for adi in meta_df['LabNumber']]
# meta_df = meta_df.set_index('LabNumber')
# meta_df = meta_df[~meta_df.index.duplicated(keep='first')] #drop duplicates rows in the metadata file
# times_df = meta_df


h0 = [agi for agi in agg_df.index if agi in times_df.index]

hMissing = [agi for agi in agg_df.index if agi not in times_df.index]
if len(hMissing)>0:
    print('Not all samples are in metadata.')
    print("Missing ", hMissing)

duplicates = agg_df[agg_df.index.duplicated(keep=False)]
if len(duplicates)>0:
    print('Samples are duplicated.')
    print(duplicates)

queryType = "linDict"
# build a dataframe 
df_abundances = pd.DataFrame()
for i, sampLabel in enumerate(agg_df.index):
    dat = agg_df.loc[sampLabel, queryType]
    if isinstance(dat, list):
        df_abundances = pd.concat([
            df_abundances,
            pd.Series(
                agg_df.loc[sampLabel, queryType][0],
                name=times_df.loc[sampLabel,
                                  'SampleCollectionDate'])
        ], axis=1)
    else:
        df_abundances = pd.concat([
            df_abundances,
            pd.Series(
                agg_df.loc[sampLabel, queryType],
                name=times_df.loc[sampLabel,
                                  'SampleCollectionDate'])
        ], axis=1)

# fill nans, group data by the appropriate interval
df_abundances = df_abundances.fillna(0).T
# df_abundances = df_abundances.loc[df_abundances.index>="2022-04-01"]

df2 = pd.melt(df_abundances.reset_index(), id_vars='index')
df2.columns = ['collection_date', 'lineage', 'abundance']
df2.to_csv('aggregated_dated.csv')
### in case of round off error, we'll add the leftover bits to Other. 

df_abundances['Other'] += 1.- df_abundances.sum(axis=1)
### now prepare the plot. 
ordering = [v['name'] for v in plot_config.values()] + ['Recombinants','Other']
colors0 = [v['color'] for v in plot_config.values()] + ["#6A6C6E","#DDDDDD"]

presentInds = [i for i in range(len(ordering)) if ordering[i] in df_abundances]
ordering = [ordering[ind] for ind in presentInds]
colors0 = [colors0[ind] for ind in presentInds]
colorDict = {l:c for l,c in zip(ordering,colors0)}
with open('color_map.json', 'w') as f0:
    json.dump(colorDict, f0)
colors0 = colors0[::-1]


df = df_abundances[ordering[::-1]]
### group by month. 
df_ = df.groupby(pd.Grouper(freq='MS')).mean()
df_.to_csv('NICD_monthly.csv')
fig,ax = plt.subplots(figsize=(10.5,5))

for i in range(0, df_.shape[1]):
    ax.bar(df_.index, df_.iloc[:, i],
           width=20, bottom=df_.iloc[:, 0:i].sum(axis=1),
           label=df_.columns[i], color=colors0[i])
    # ax.stackplot(df_.index,df_.T,labels=df_.columns,colors=colors0)
ax.set_xlim(df_.index.min()-timedelta(days=15),df_.index.max()+timedelta(days=15))
ax.set_ylim(0,1)
locator = mdates.MonthLocator(bymonthday=1)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1],loc='center left', bbox_to_anchor=(1, 0.5))
fig.tight_layout()

plt.savefig('../figures/NICD_stackplot_Monthly.pdf')



### group by week. 
df_ = df.groupby(pd.Grouper(freq='W-SAT')).mean()

fig,ax = plt.subplots(figsize=(10.5,5))

for i in range(0, df_.shape[1]):
    ax.bar(df_.index, df_.iloc[:, i],
           width=7, bottom=df_.iloc[:, 0:i].sum(axis=1),
           label=df_.columns[i], color=colors0[i])
# ax.stackplot(df_.index,df_.T,labels=df_.columns,colors=colors0)
ax.set_xlim(df_.index.min()-timedelta(days=15),df_.index.max()+timedelta(days=15))
ax.set_ylim(0,1)
locator = mdates.MonthLocator(bymonthday=1)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1],loc='center left', bbox_to_anchor=(1, 0.5))
fig.tight_layout()
plt.savefig('../figures/NICD_stackplot_weekly.pdf')


### average across days, rolling 
df_ = df.groupby(pd.Grouper(freq='D')).mean()
windowSize=28
df_ = df_.rolling(windowSize, center=True,min_periods=0).mean()

fig,ax = plt.subplots(figsize=(10.5,5))
ax.stackplot(df_.index,df_.T,labels=df_.columns,colors=colors0)
locator = mdates.MonthLocator(bymonthday=1)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1],loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_xlim(df_.index.min(),df_.index.max())
ax.set_ylim(0,1)
ax.set_ylabel('Lineage prevalence')
fig.tight_layout()

plt.savefig('../figures/NICD_stackplot_daily.pdf')
plt.close('all')



