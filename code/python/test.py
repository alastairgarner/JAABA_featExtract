

### Import Dependencies ###

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob, os
import yaml

### Set Parameters



### FUNCTIONS ### 
def get_params(filename="params.yaml"):
    files = glob.glob('*.yaml')
    if filename in files:
        file_open = files[files.index(filename)]
    else:
        file_open = files[0]

    print("Opening {} parameter file".format(file_open))

    with open(file_open) as f:
        params = yaml.load(f, Loader=yaml.FullLoader)
    return params


def get_filelist(data_type="salam",file_type='.txt'):
    search_dir = os.path.join(params["directories"][data_type],'**',('*'+file_type))
    f = []
    for file in glob.glob(search_dir,recursive=True):
        f.append(file)
    return f


### SCRIPTING ###

if __name__ == '__main__':
    params = get_params()

    f = get_filelist(data_type="chore_output",file_type="compiledChore.txt")

    df = pd.read_csv(f[0], sep=" ", header=0)
    df.set_index(["timestamp","id"],inplace=True)
    ave = df.groupby(['time']).mean()
    sem = df.groupby(['time']).sem()
    tracked = df.groupby(['time']).count()["speed"]
    df["time"].max()

    bins = np.arange(0,np.ceil(df["time"].max()),0.5)
    labels = np.arange(0,bins.size-1)
    df["bins"] = pd.cut(df["time"],bins, labels=labels)
    ave = df.groupby(['bins']).mean()
    sem = df.groupby(['bins']).sem()

    sns.lineplot(data=ave["speed"])
    ax = sns.lineplot(x=df["bins"],y=df["curve"])
    plt.show(ax)



    
    

    print("Finished")