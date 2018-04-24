# The MIT License
# 
# Copyright (c) 2017 Piero Dalle Pezze
# 
# Permission is hereby granted, free of charge, 
# to any person obtaining a copy of this software and 
# associated documentation files (the "Software"), to 
# deal in the Software without restriction, including 
# without limitation the rights to use, copy, modify, 
# merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom 
# the Software is furnished to do so, 
# subject to the following conditions:
# 
# The above copyright notice and this permission notice 
# shall be included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR 
# ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


"""
Author: Piero Dalle Pezze
Affiliation: The Babraham Institute (Cambridge, UK)
Aim: A Snakemake workflow to process time course data of generic autophagy using ATG13 as readout.
Date: 30 Nov 2017
Run: snakemake -s Snakefile --configfile config.yaml
"""


import sys
import os


# ---- DATA ---- #

csv_ext  = ".csv"
plot_ext = ".png"
data_dir = "data"

# import from yaml config file
datasets = config['datasets']
thres_lv = config['thres_lv']  
thres_hv = config['thres_hv']
readout = config['readout']
models_dir = config['models_dir']
observable = config['observable']
treatment = config['treatment']

# remove file extension
datasets = [os.path.splitext(ds)[0] for ds in datasets]


# ---- FUNCTIONS ---- #

# Use these functions to get the corresponding threshold values for to use with the data set (wildcard)
def get_thres_lv(wildcards):
    wildcards = os.path.splitext(str(wildcards))[0]
    if wildcards in datasets:
        index = datasets.index(wildcards)
        return thres_lv[index]
    return 0

def get_thres_hv(wildcards):
    wildcards = os.path.splitext(str(wildcards))[0]
    if wildcards in datasets:
        index = datasets.index(wildcards)
        return thres_hv[index]
    return 0

def get_treatment(wildcards):
    wildcards = os.path.splitext(str(wildcards))[0]
    if wildcards in datasets:
        index = datasets.index(wildcards)
        return treatment[index]
    return 'FALSE'
    
# ---- RULES ---- #

rule all:
    input:
        ds1 = expand(os.path.join(data_dir, "{dataset}" + "_sync" + csv_ext), dataset=datasets),
        ds2 = expand(os.path.join("scripts", "2_upper_lower_bound_analysis", "{dataset}" + "_sync_peaks_stats" + csv_ext), dataset=datasets),
        ds3 = expand(os.path.join(models_dir, "{dataset}" + "_copasi" + csv_ext), dataset=datasets)
        
rule synchronise_by_maxval:
    input:
        file = os.path.join(data_dir, "{dataset}" + csv_ext)
    output:
        file = os.path.join(data_dir, "{dataset}" + "_sync" + csv_ext)
    params:
        readout = readout
    shell:
        "Rscript scripts/1_synchronise_by_maxval/synchronise_by_maxval.R {input.file} {params.readout}"

rule upper_lower_bound_analysis:
    input:
        file = os.path.join(data_dir, "{dataset}" + "_sync" + csv_ext)
    output:
        file = os.path.join("scripts", "2_upper_lower_bound_analysis", "{dataset}" + "_sync_peaks_stats" + csv_ext)
    params:
        thres_lv = get_thres_lv,    
        thres_hv = get_thres_hv 
    shell:
        "Rscript scripts/2_upper_lower_bound_analysis/upper_lower_bound_analysis.R {input.file} {params.thres_lv} {params.thres_hv}"        

rule time_courses_data_for_copasi:
    input:
        file = os.path.join(data_dir, "{dataset}" + "_sync" + csv_ext)
    output:
        file = os.path.join(models_dir, "{dataset}" + "_copasi" + csv_ext)
    params:
        models_dir = models_dir,
        outputfile = "{dataset}" + "_copasi",
        observable = observable,
        treatment = get_treatment
    shell:
        "Rscript scripts/3_time_courses_data_for_copasi/time_courses_data_for_copasi.R {input.file} {params.models_dir} {params.outputfile} {params.observable} {params.treatment}"
