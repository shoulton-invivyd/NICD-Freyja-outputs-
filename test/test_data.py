import pytest
import pandas as pd
import pandas.api.types as ptypes
import yaml

agg_results = 'agg_demixed.tsv'
meta = 'sample_metadata.csv'
@pytest.fixture
def agg_df():
    agg_ = pd.read_csv(agg_results, skipinitialspace=True, sep='\t', index_col=0)
    return agg_

@pytest.fixture
def times_df():
    times_ = pd.read_csv(meta, skipinitialspace=True, index_col=0)
    return times_


def test_agg_dups(agg_df):
    assert agg_df[agg_df.index.duplicated(keep=False)].shape[0] == 0,'Aggregate has duplicate entries'


def test_dups_meta(times_df):
    # check for  duplicated seq ids
    print(times_df)
    assert times_df[times_df['Sequence_ID'].duplicated(keep='first')].shape[0] ==0, 'Duplicate seq IDs'
    # check for duplicated lab numbers
    assert times_df[times_df['LabNumber'].duplicated(keep=False)].dropna().shape[0] == 0, 'Duplicate lab numbers'
