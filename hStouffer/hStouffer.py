#=============================================================================#
# Author    : DaeHee Kim
# Date      : 2025-05-23
# Usage     : Perform meta-analysis by integrating DESeq2 result files
# Example   : hStouffer.py -d dataset/low_fat/ -o meta_low_fat
# Description   : hStouffer is a meta-analysis method that combines p-values from multiple studies to identify differentially expressed genes (DEGs) across datasets.
#=============================================================================#

import argparse
import glob
import time
from collections import defaultdict
from scipy.stats import combine_pvalues
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm
import random

parser = argparse.ArgumentParser()
parser.add_argument('-d', action='store', help='Directory containing DESeq2 result files')
parser.add_argument('-p', action='store', help='p-value threshold (e.g. 3 = 1e-3)', default=False)
parser.add_argument('-o', action='store', help='Output file name')
parser.add_argument('-m', action='store', help='Group Size', default= False)
parser.add_argument('-r', action='store', help='Number of repetitions', default= 1000)
parser.add_argument('-n', action='store', help='Apply scaling (True/False)', default= False)
parser.add_argument('-t', action='store', help="Method for p-value combination {'fisher', 'pearson', 'tippett', 'stouffer', 'mudholkar_george', 'median', 'percentile_70'}", default= "stouffer")
parser.add_argument('-c', action='store', help='Number of CPU cores', default= 32)
parser.add_argument('-w', action='store', help='Apply p-value capping (True/False)', default= True)
parser.add_argument('-l', action='store', help='Max or min cutoff threshold', default= False)

args = parser.parse_args()
directory = args.d
output = args.o
dup = args.n
is_winsored = args.w
max_min = args.l

DESeq2_files = glob.glob(f"{directory}/*deg.txt")
if dup:
    dup = int(dup)
    random_dup = [random.choice(DESeq2_files) for _ in range(dup-len(DESeq2_files))]
    DESeq2_files += random_dup

repeat = int(args.r)
combination_method = args.t
winsor = 1e-3
core = int(args.c)
combi = len(DESeq2_files)

x = np.array([10,20,30,40,50,60,70,80,90,100,
              110,120,130,140,150,160,170,180,190,200])
average_y = np.array([1.301029996, 2.650514998, 4.650514998, 7, 10.5, 13, 15.5, 18.5, 22, 25,
              27.5, 30.5, 34.5, 36.5, 40.5, 43.5, 46.5, 50, 53.5, 57])
min_y = np.array([1.301029996, 4, 8, 12, 18, 22, 27, 33, 39, 44, 49, 54, 61, 65, 72, 77, 83, 89, 95, 101])
max_y = np.array([1.301029996, 1.301029996, 1.301029996, 2, 3, 4, 4, 4, 5, 6, 6, 7, 8, 8, 9, 10, 10, 11, 12, 13
])

no_bagged_y = np.array([1.650514998, 7.5, 11, 14.5, 18, 22, 25.5, 29.5, 33, 36.5,
              40.5, 44, 48, 51.5, 55.5, 58.5, 62.5, 66, 70, 73.5])

def linear_interpolation(x, y, predict_x):
    if predict_x < x[0] or predict_x > x[-1]:
        raise ValueError("Out of range.")

    for i in range(len(x)-1):
        if x[i] <= predict_x <= x[i+1]:
            x1, x2 = x[i], x[i+1]
            y1, y2 = y[i], y[i+1]
            
            predicted_y = y1 + (y2 - y1) * (predict_x - x1) / (x2 - x1)
            return predicted_y
        
    raise ValueError("Not proper size.")

if args.p:
    p = str(args.p)
    pval_threshold = float("1e-"+p)
else:
    pval_threshold = linear_interpolation(x, average_y, combi)

if max_min:
    if max_min == 'max':
        pval_threshold = linear_interpolation(x, max_y, combi)
    elif max_min == 'min':
        pval_threshold = linear_interpolation(x, min_y, combi)

pval_threshold = 10**(-float(pval_threshold))  
print(len(DESeq2_files),combi, pval_threshold)


start_time = time.time()

def benjamini_hochberg(p_values, alpha=0.05):
    p_values = np.asarray(p_values)
    n = len(p_values)
    sorted_indices = np.argsort(p_values)
    sorted_p_values = p_values[sorted_indices]
    ranks = np.arange(1, n + 1)
    thresholds = (ranks / n) * alpha
    sorted_adjusted_p_values = np.minimum.accumulate(sorted_p_values[::-1] * n / ranks[::-1])[::-1]
    adjusted_p_values = np.empty(n)
    adjusted_p_values[sorted_indices] = np.minimum(sorted_adjusted_p_values, 1.0)
    return adjusted_p_values

gses = []
basemean = defaultdict(list)
log2Fold = defaultdict(list)
pval = defaultdict(list)

for i, dsq in enumerate(sorted(DESeq2_files)):
    print(f"{i+1} / {len(DESeq2_files)} DESeq2 file loading           ", end="\r")

    df = pd.read_csv(dsq, sep="\t")
    for _, row in df.iterrows():
        gene = row.iloc[0]
        if is_winsored == True:
            if row.iloc[5] < winsor:
                pval[gene].append(winsor)
            else:
                pval[gene].append(row.iloc[5])
        else:
            pval[gene].append(row.iloc[5])

print("\nDESeq2 data loaded!")

data_nums = len(DESeq2_files)

combinations = set()
while len(combinations) < repeat:
    combination = tuple(sorted(random.choices(range(combi), k=data_nums)))
    combinations.add(combination)
combinations = [list(c) for c in combinations]
    
def median_p(pvals):
    pvals = np.array(pvals)
    pvals = pvals[~np.isnan(pvals)]
    if len(pvals) == 0:
        return np.nan
    return np.median(pvals)
def percentile_70(pvals):
    pvals = np.array(pvals)
    pvals = pvals[~np.isnan(pvals)]
    if len(pvals) == 0:
        return np.nan
    return np.percentile(pvals, 30)

def process_combination(combi, pval):
    stouffer_arr, genes_arr = [], []
    epsilon = np.finfo(float).eps 
    for gene in sorted(pval.keys()):
        org_pval = np.array(pval[gene], dtype=float)
        sub_pval = org_pval[np.array(combi)]
        pvals = sub_pval[~np.isnan(sub_pval)]
        pvals[pvals == 0] = epsilon
        if len(pvals) > 0:
            genes_arr.append(gene)
            
            if combination_method in ['fisher', 'pearson', 'tippett', 'stouffer', 'mudholkar_george']:
                _, combined_p = combine_pvalues(pvals, method= combination_method)
            elif combination_method == 'median':
                combined_p = median_p(pvals)
            elif combination_method == 'percentile_70':
                combined_p = percentile_70(pvals)
            stouffer_arr.append(combined_p)
            

    adjusted_stouffer = benjamini_hochberg(stouffer_arr)
    return genes_arr, adjusted_stouffer

results = []
parallel = Parallel(n_jobs=core, backend="loky")
results = parallel(delayed(process_combination)(combi, pval) for combi in tqdm(combinations, desc="Combination Processing"))

fisher_dict = defaultdict(list)
stouffer_dict = defaultdict(list)

for genes_arr, adjusted_stouffer in results:
    for n, gene in enumerate(genes_arr):
        stouffer_dict[gene].append(adjusted_stouffer[n])
        
final_result = defaultdict(list)
for gene in stouffer_dict.keys():
    
    total = len(stouffer_dict[gene])
    stouffer_deg_num = len([p for p in stouffer_dict[gene] if p < pval_threshold])
    final_result[gene] = stouffer_deg_num/total


number_of_degs = 0
for sto in final_result.values():
    if (sto >= 0.95):
        number_of_degs += 1
print(number_of_degs)

with open(output, "w") as f:
    f.write("Gene\tStouffer\n")
    f.write(f"Total DEGs\t{number_of_degs}\n")
    for gene in sorted(final_result.keys()):
        f.write(f"{gene}\t{final_result[gene]}\n")


end_time = time.time()

print(f"Start time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}")
print(f"End time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}")
print(f"Total duration: {end_time - start_time:.2f} seconds")