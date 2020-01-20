### main_plots.py


### Dependencies

import numpy as np
import pandas as pd
import scipy.io as spio

import glob, os, re
import yaml

import seaborn as sns
import matplotlib.pyplot as plt


### FUNCTIONS

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

def get_mat_filelist(directory="./data_compiled"):
    search_dir = os.path.join(directory,'**',('*'+".mat"))
    f = []
    for file in glob.glob(search_dir,recursive=True):
        f.append(file)
    return f

def get_file_details(filepath):
    expressions = []
    expressions.append(r"(?P<date>\d{8})[_](?P<time>\d{6})")
    expressions.append(r"[@\/](?P<driver>(?<=\d{8}_\d{6}[@\/])[\w.]*)"
                r"[@](?P<effector>\w*)")
    expressions.append(r"[@](?P<tracker>[t]\d{2})")
    expressions.append(r"[@#]*(?P<protocol1>[\w\d_]*)"
                r"[#](?P<protocol2>[\w\d_]*)"
                r"[#](?P<protocol3>[\w\d_]*)"
                r"[#](?P<protocol4>[\w\d_]*)")

    file_deets = dict()
    for expr in expressions:
        regex = re.compile(expr)
        result = regex.search(filepath)
        if result is None:
            file_deets.update(dict(regex.groupindex))
            continue
        file_deets.update(result.groupdict())
    file_deets["display"] = "{driver}@{effector}@{protocol1}@{protocol2}".format(**file_deets)
    return file_deets

def group_genotypes(filelist):
    deets = list(map(get_file_details,filelist))
    exps = [deet["display"] for deet in deets]
    un_exps, indices = np.unique(exps, return_inverse=True)
    groups = [filelist[ii==indices] for ii in range(len(un_exps))]
    return groups

def mat_check_keys(dict):
    """
    sdfsdf
    """
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = mat_to_dict(dict[key])
    return dict 

def mat_to_dict(matobj):
    """
    A recursive function which constructs from matobjects nested dictionaries
    """
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = mat_to_dict(elem)
        # elif (hasattr(elem,"__iter__")) and ([isinstance(el, spio.matlab.mio5_params.mat_struct) for el in elem]):
        #     for ii in range(len(elem)):
        #         if isinstance(elem[ii], spio.matlab.mio5_params.mat_struct):
        #             # print(elem[ii])
        #             elem[ii] = mat_to_dict(elem[ii])
        #     dict[strg] = elem
        else:
            dict[strg] = elem
    return dict


### SCRIPTING

if __name__ == "__main__":
    print("Yeah Boi")

    searchterm = "2019082"

    params = get_params(filename="params_AG.yaml")
    f = np.array(get_mat_filelist("./data_processed"))
    groups = group_genotypes(f)

    filename = groups[0][0]
    data = spio.loadmat(filename, squeeze_me=True, struct_as_record=False,variable_names="temp")
    data = data["temp"].data_chore
    data = [mat_to_dict(dct) for dct in data]
    data[0]
    ndict = dict()
    for fname in data[0].keys():
        temp = np.array()
        ndict[fname] = [dat[fname] for dat in data]
    
    ndict["et"]
    pd.DataFrame(ndict)

    
    
    
    
    
    
    




matfile.keys()

fig,ax = plt.subplots()
ax.plot(x,y)
plt.show()

