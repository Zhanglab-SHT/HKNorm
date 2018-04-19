#!/usr/bin/python
"""
The main purpose of the script is to  select the hk_genes and caculate size_factor based on those selected hk_genes

=============================
Usage: 
-h:	help

-i: input the raw_expression_matrix                                [No default]

-n:	the number group what you want to split the the hk_gene.csv    [default value 4]

-p:	the threshold value that select hk_genes                       [default value 0.5]

-o: output the file name of hk_gene_corr.csv                  		   [No default]

-N:	output the file name of normalization data.csv            		   [No default]

-s: the file name of size_factor                                   [default file name :size_factor.csv]

-m: output the file name of hk_genes_mean.csv                      [No default]

-k: input the option and the 'h' represents human_hk_genes ,the 'm' represents 'mouse_hk_genes'   [No default]


===================

input description:

input files:

1.raw_expression_matrix

2.alternative reference hk_gene_names file

======================

output files and figures

1.size_factor.csv

2.corr.csv

3.optimal_hk_gene_corr.figure

4.hk_gene_mean.figure

5.alternative normalizd_data.csv

6.alternative hk_mean.csv


"""

''' 
command line

'''



import sys
import argparse
import pandas as pd
from pandas import DataFrame
from pandas import Series
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
parser = argparse.ArgumentParser('selecting hk_genes to estimate size_factors')
parser.add_argument('I',type = str, help = 'input the raw_data')
parser.add_argument('S',type = str, help = 'input the file name of  the size_factor data',const = 1 , nargs = "?", default = 'size_factor.csv')
parser.add_argument('-b',type = str, help = 'choose the built-in reference hk_genes.csv')
parser.add_argument('-lo',type = str, help = 'choose yourself local hk_genes.csv')
parser.add_argument('-n',type = int, help = 'input the number group what you want to split the the hk_gene.csv',const = 1 , nargs = "?", default = 4)
parser.add_argument('-p',type = float, help = 'the threshold value that select hk_genes',const = 1 , nargs = "?", default = 0.5)
parser.add_argument('-o',type = str, help = 'input the filename of the optimal_hk_gene_corr data',const = 1 , nargs = "?", default = 'corr.csv')
parser.add_argument('-N',type = str, help = 'input the file name of  the normalization data')
parser.add_argument('-m',type = str, help = 'input the file name of  the hk_mean.csv')
args = parser.parse_args()

### define some global variables
CORR = DataFrame()
SF=[]

### define some functions

def corr(x):
	corr_b=hk_genes_matrix
	corr_a=DataFrame(x)
	select_gene_list=[]
	max_pearson_list=[]
	max_pearson=1
	while abs(max_pearson)>= args.p:
		c_frame=pd.Series()
		for i in corr_b.columns:
			f1=lambda y:y.corr(corr_b[i])
			d_frame=corr_a.apply(f1)
			c_frame=pd.concat([c_frame,d_frame],axis=1)
		len_col = len(c_frame.columns)
		c_frame = c_frame.iloc[:,1:len_col]
		c_frame = c_frame.fillna(0)
		c_frame.columns=corr_b.columns
		e_series=c_frame.min()
		select_gene=e_series.idxmax()
		max_pearson=e_series.max()
		if abs(max_pearson) >= args.p:
			corr_a[select_gene]=corr_b[select_gene]
			corr_b=corr_b.drop([select_gene],axis=1)
			select_gene_list.append(select_gene)
			max_pearson_list.append('%.4f' % max_pearson)
		else:
			break
	c={'pc':max_pearson_list,'gene_name':select_gene_list}
	c=DataFrame(c)
	print('.')
	global CORR
	CORR=pd.concat([c,CORR],axis=1)


def plot_corr(x):
	plot_corr_frame = pd.read_csv(x,header = 0)
	red_dot_frame = plot_corr_frame.iloc[0:20]
	black_dot_frame = plot_corr_frame.iloc[20:]
	x1=red_dot_frame.index
	y1 = red_dot_frame['pc']
	x2=black_dot_frame.index
	y2 = black_dot_frame['pc']
	plt.figure()
	ax = plt.subplot()
	ax.scatter(x1,y1,c = 'red')
	ax.scatter(x2,y2,c = 'black')
	plt.savefig("corr.png")
	plt.ylabel('corr')
	plt.xlabel('hk_gene')


def plot_mean (x,y):
	plt.figure()
	plt.scatter(x,y,c= 'black')
	plt.ylabel('counts mean')
	plt.xlabel('hk_gene')
	plt.savefig("mean.png")


def sf(x):
	select_gene_rawdata_gene_mean = x
	select_gene_rawdata_gene_mean_array = select_gene_rawdata_gene_mean.values
	select_gene_rawdata_gene_mean_array = select_gene_rawdata_gene_mean_array.T
	model = LinearRegression(fit_intercept=False)
	c = len(select_gene_rawdata_gene_mean_array) -1
	x = select_gene_rawdata_gene_mean_array[c]                           
	x = x.reshape(-1,1)
	b = 0
	for i in select_gene_rawdata_gene_mean_array:
		y = i.reshape(-1,1)                                                                                                 
		model.fit(x,y)                                              
		a = model.coef_[0,0]
		global SF
		SF.append(a)




### code begins
if args.I is None:
	print(
		'Please ensure that raw data file has been entered'
		 )
	sys.exit(0)
if args.b is None and args.lo is None:
	print(
		'Please input raw data and options'
		 )
	sys.exit(0)
if args.b is not None and args.lo is not None:
	print(
		'Please input raw data and options'
		 )
	sys.exit(0)
'''if args.lo and args.b is None:
	print('Please ensure that files and options have been entered')
	sys.exit(0)'''
'''if args.I and args.b is not None:
	print('Please ensure that files and options have been entered')
	sys.exit(0)'''
if args.b is not None:
	hk_list = pd.read_csv('hk_genes.csv',header = 0 ,index_col = False)
	if args.b == 'he':
		hk_list = hk_list['ENS_ID'].tolist()
	if args.b == 'hg':
		hk_list = hk_list['ENS_GENE'].tolist()
	if args.b == 'me':
		hk_list = hk_list['ENSMUS_ID'].tolist()
	if args.b == 'mg':
		hk_list = hk_list['ENSMUS_GENE'].tolist()
	if args.b == 'ze':
		hk_list = hk_list['ENSDAR_ID'].tolist()
	if args.b == 'zg':
		hk_list = hk_list['ENSDAR_GENE'].tolist()
	if args.b == 'fe':
		hk_list = hk_list['FB_ID'].tolist()
	if args.b == 'fg':
		hk_list = hk_list['FB_GENE'].tolist()
if args.lo is not None:
	hk_list = pd.read_csv(args.lo, index_col = 0)
	hk_list = hk_list.index.tolist()

# partI: return the hk_genes_matrix
raw_data = pd.read_csv(args.I, header = 0, index_col = 0)
raw_data_list = raw_data.index.tolist()
raw_data_list = set(raw_data_list)
hk_list = set(hk_list)
hk_list = raw_data_list & hk_list
hk_genes_matrix = raw_data.loc[hk_list]

# partII: caculate hk_genes_matrix ~~ hk_select corr
hk_genes_matrix = hk_genes_matrix.T
mean_sort = hk_genes_matrix.median()
mean_sort = mean_sort.sort_values()
mean_sort = mean_sort.index.tolist()
len_sample = len(mean_sort)
step_size = len_sample // args.n
mean_sort_list = []
i = 0
while i < len_sample:
	mean_sort_list.append(mean_sort[i])
	i += step_size
hk_select = hk_genes_matrix[mean_sort_list]
hk_select.apply(corr)

# partIII: return the best one hk_gene and plot corr figure
CORR = CORR.fillna(0)
len_col = len(CORR.columns.tolist())
list_gene_name = []
list_pcc = []
for i in xrange(0,len_col,2):
	list_gene_name.append(i)
dict_gene_name = CORR.iloc[0,list_gene_name]
for i in xrange(1,len_col,2):
	list_pcc.append(i)
dict_pcc = CORR.iloc[19,list_pcc] 
dict_int = dict(zip(dict_pcc,dict_gene_name))
dict_int_sort_keys = sorted(dict_int.keys())
if float(dict_int_sort_keys[-1]) > 0.5 :
	optimal_hk_gene = dict_int[dict_int_sort_keys[-1]]
else:
	print('dont find appropriate housekeeping gene')
	sys.exit(0)
CORR.columns = CORR.iloc[0]
hk_frame_columns = list(CORR.columns)
optimal_hk_gene_index = hk_frame_columns.index(optimal_hk_gene)
optimal_hk_gene_pcc = CORR.iloc[:,[optimal_hk_gene_index,optimal_hk_gene_index+1]]
optimal_hk_gene_pcc.columns = ['optimal_hk','pc']
optimal_hk_gene_pcc.to_csv(args.o,index = False, header = True)
plot_corr(args.o) # apply the plot_corr function
	

# partIV: plot mean figure 
###: caculate the hk_mean
select_hk_gene_list = optimal_hk_gene_pcc['optimal_hk'].tolist()
optimal_hk_gene_pcc_frame = hk_genes_matrix[select_hk_gene_list]
optimal_hk_gene_pcc_frame = optimal_hk_gene_pcc_frame.T
sum_series = optimal_hk_gene_pcc_frame.apply(lambda x: x.sum())
sum_frame = DataFrame(sum_series)
sum_frame = sum_frame.T
sum_frame.index = ['col_sum']
optimal_hk_gene_pcc_frame = optimal_hk_gene_pcc_frame.append(sum_frame)
sum_series = optimal_hk_gene_pcc_frame.loc['col_sum']
sum_median = sum_series.median()
mean_fig_frame= (optimal_hk_gene_pcc_frame / sum_series) * sum_median
row_mean_series = mean_fig_frame.mean(axis = 1)
optimal_hk_gene_pcc_frame['row_mean'] = row_mean_series
optimal_hk_gene_pcc_frame.drop(['col_sum'],axis = 0)
###:plot mean figure
row_mean_series_sorted = row_mean_series.sort_values()
row_mean_series_sorted = row_mean_series_sorted.drop('col_sum') ### y1
row_mean_series_sorted_reset_index = row_mean_series_sorted.reset_index(drop=True) ### x1
x1 = row_mean_series_sorted_reset_index.index
y1 = row_mean_series_sorted.tolist()
a_list=row_mean_series_sorted.index.tolist()
b_list = row_mean_series_sorted.tolist()
c_dict={'hk_gene_name':a_list,'counts_mean':b_list}
c_frame = DataFrame(c_dict)
if args.m is not None:
	columns = ['hk_gene_name','counts_mean']
	#counts_format = c_frame.apply(lambda x: '%.f' % x)
	c_frame.to_csv(args.m,header = True ,index=False,columns = columns)
plot_mean(x1,y1) #apply the plot_mean function

# partV: caculate size_factor
sf(optimal_hk_gene_pcc_frame) # apply sf function and return the SF
sf_series=Series(SF)
a = len(sf_series) -1
sf_series = sf_series.drop(a)
sf_frame = DataFrame(sf_series)
sf_frame = sf_frame.T
sf_frame.columns=raw_data.columns
if args.N is not None:
	join_frame = pd.concat([raw_data,sf_frame],ignore_index=True)
	b  = len(join_frame.index) -1
	sf = join_frame.iloc[b]
	raw_frame_normalize = raw_data / sf
	raw_frame_normalize = raw_frame_normalize.applymap(lambda x: '%.f' % x)
	raw_frame_normalize.to_csv(args.N,header = True,index=True)
sf_frame.index = ['size_factor']
sf_frame = sf_frame.T
sf_frame = sf_frame.applymap(lambda x: '%.2f' % x)
sf_frame.to_csv(args.S,header=True,index=True)




