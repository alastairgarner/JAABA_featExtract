# JAABA_featExtract

This is the code I used to generate analyses for my MSc thesis. It's a mess, but I'll try and outline a few of the key scripts.

**code/params_AG.yaml**

Configuration file where you can set the location of directories that contain different data types (MWT, Choreography, Salam, JAABA, JB). Other scripts rely on this file to know where to find the data.

**code/MATLAB/main_process.m**

Takes raw data from each of the pipelines used in the Ohyama lab (MWT, choreography, Salam, JAABA, JB) and transforms them to a single, common format.

**code/MATLAB/main_setup.m**

> Must be run before *individual_plots.m*

Searches the config directories for valid files. A selection dialogue box allows you to choose which experiment to load, prompting the relevant files to be copied to a temporary local directory. This was implemented because I kept most of my data on an external hard drive, and read/write operations from external drives are slow.

**code/MATLAB/indiviaul_plots.m**

> Must be run after *main_setup.m*

Generates plots & stats for an experiment. Given the experiment selected in the previous step, you now define which genotypes to analyse. This is done either through a GUI file selector, or by using regular expressions (regex, defined at the top of the script.

**./blacklist.txt**

A list of timestamps that identifies experiments which ran incorrectly (stimulus didn't occur, timing was wrong, software crashed, etc). These timestamps are suppressed when performing analyses, to prevent anomalous results. This file can be added to.

---

## Third Party packages

I'd like to express my thanks to the following third party packages, which were essential for this work and are included in this repo.

- [cbrewer](https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab?s_tid=srchtitle)
- [Perceptually Uniform Colormaps](https://www.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps?s_tid=srchtitle)
- [plotboxpos](https://github.com/kakearney/plotboxpos-pkg)
- [plotSpread](https://www.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points-beeswarm-plot)
- [Violinplot](https://github.com/bastibe/Violinplot-Matlab)
- [yamlmatlab](https://github.com/ewiger/yamlmatlab)