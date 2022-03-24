"""
#This script was made by @roryflynn in colaboration with myself providing knowledge of output and the tools used.


The rules should be:
If "this_study" only and more than 1 in cluster then = novel genera
If "this_study" plus "input_db" then = novel genera clustering with input_db
If *nothing in header" plus "this_study" plus can contain  "any other flag" then = known taxonomy

And the output should then be a subset of the entire csv, where you pull out only clusters that have "this_study" flag, and the label from the rule in a separate column+ the list of the different viruses that are in that cluster in case more manual parsing is needed
"""
import re
import pandas as pd
import numpy as np
import click


# assert re.findall('[A-z,0-9]*_', "dkjkd_4334_djkfdj3_bla~bla") = "dkjkd_4334_djkfdj3_"
def check_for_other(df):
    for i in df['Genome'].unique():
        j = i.replace('this_study_', '').replace('input_db', '')
        if len(re.findall('[A-z,0-9]*_', j)) > 1:
            return True
    return False

def check_rule(df:pd.DataFrame):
    this_study = ['this_study' in i for i in df['Genome'].unique()]
    if not np.any(this_study):
        return pd.DataFrame()
    if np.all(this_study):
        df['Auto_type'] = 'novel genera'
        return df
    input_db = ['input_db' in i for i in df['Genome'].unique()]
    if np.all(np.logical_or(input_db, this_study)):
        df['Auto_type'] = 'novel genera clustering with input_db'
        return df
    #if check_for_other(df)
    df['Auto_type'] = 'known taxonomy'
    return df
    raise ValueError("Reaching this error should not be possible based on rules.")


@click.command('distill')
@click.option('-i', '--input', help="A vContact2 file, to parse",
              required=True)
@click.option('-o', '--output', help="output name, optional",
              required=False)
def parse(input, output:str = None):
    data = pd.read_csv(input, index_col=0)
    data = data[data['Size'] > 1]
    data = data.groupby('VC').apply(check_rule)
    if output is None:
        output = input.replace('.csv', "_auto_typed.csv")
    data.to_csv(output)

if __name__ == '__main__':
    parse()
"""
import os
os.system("python ./vcontact2_parser.py -i genome_by_genome_overview_rorytest.csv")
exit


"""
