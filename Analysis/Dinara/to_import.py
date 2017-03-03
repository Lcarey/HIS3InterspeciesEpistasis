#http://scipy.github.io/old-wiki/pages/Cookbook/Matplotlib/Show_colormaps
#http://matplotlib.org/1.4.0/examples/images_contours_and_fields/interpolation_methods.html
#http://chrisalbon.com/python/set_the_color_of_a_matplotlib.html
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg 
mpl.rcParams['pdf.fonttype'] = 42

import numpy as np
from StringIO import StringIO
import pandas as pd
import os
from IPython.display import display
from collections import defaultdict
import csv
import scipy 

from IPython.display import display
from collections import defaultdict
from collections import Counter
from itertools import combinations
from scipy import stats
from scipy import optimize
import os
import pandas as pd

import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import numpy as np
from random import randint
from matplotlib.colors import LogNorm

import matplotlib.cm as cm

import time
pd.set_option('display.max_columns', None)
from scipy.optimize import curve_fit
from scipy import interpolate

def hamdist(str1, str2, ignoregaps=False):
    diffs = 0
    if len(str1) != len(str2):
          return max(len(str1),len(str2))
    for ch1, ch2 in zip(str1, str2):
        if (ignoregaps):
            if ((ch1!='-') & (ch2!='-') & (ch1 != ch2)):
                diffs+=1
        else:
              if ch1 != ch2:
                    diffs += 1
    return diffs

def linear_fit_min_res(x,y):
    x=np.array(x)
    y=np.array(y)
    xy=np.mean(x*y)
    x2=np.mean(x*x)
    y2=np.mean(y*y)
    a=xy
    b=x2-y2
    c=-xy
    k1=(-b+np.sqrt(b**2-4*a*c))/2/a
    k2=(-b-np.sqrt(b**2-4*a*c))/2/a
    r1=sum((y-k1*x)**2)
    r2=sum((y-k2*x)**2)
    if (r2>r1):
        return [k1,r1]
    else:
        return [k2,r2]
    return np.sqrt(np.nanmean(((np.array(ydata)-ypred)/sigmas)**2))
    
# function used for fit
def func_single_exponent(x,y0i,si):
    return y0i*np.exp(si*x)
def func_calc_residuals(ydata,xdata,sigmas,popt):
    ypred=map(lambda x: func_single_exponent(x,popt[0],popt[1]),xdata)
    return np.sqrt(np.nanmean(((np.array(ydata)-ypred)/sigmas)**2))    
def func_fit_single_exponents(data,data_sigma,s12,k,fit_all=True):
    tmp_residuals=[]
    tmp_y0=[]
    tmp_s=[]
    tmp_y0_std=[]
    tmp_s_std=[]
    for i in range(len(data)):
        ydata=data[i].copy()
        ydata[2]=ydata[2]*np.exp(s12*k) 
        sigmas=data_sigma[i].copy()
        sigmas[2]=sigmas[2]*np.exp(s12*k)
        xdata=[0,1,k]

        try:
            popt,pcov=scipy.optimize.curve_fit(func_single_exponent,xdata,list(ydata),sigma=sigmas,p0=[ydata[0],0],absolute_sigma=False)
        except RuntimeError:
            popt=[np.nan,np.nan]
            pcov=[[np.nan,np.nan],[np.nan,np.nan]]
        tmp_residuals.append(func_calc_residuals(ydata,xdata,sigmas,popt)) 
        tmp_y0.append(popt[0])
        tmp_s.append(popt[1])
        tmp_y0_std.append(np.sqrt(np.diag(pcov))[0])
        tmp_s_std.append(np.sqrt(np.diag(pcov))[1])
        
    if (fit_all==True):
        return tmp_y0, tmp_s, tmp_y0_std, tmp_s_std, tmp_residuals
    else:
        bad_residuals_limit=np.sort(tmp_residuals)[int(0.95*len(tmp_residuals))]
        residuals_mean=np.nanmean([item for item in tmp_residuals if item<bad_residuals_limit])
        return residuals_mean 
    #return np.nanmean([item for item in tmp_residuals if item<bad_residuals_limit])
    
def int_array_from_string(string):
    for ch in ["[","]","\r\n","'",","]:
        string=string.replace(ch,'')
    array = [int(s) for s in string.split()] 
    return array

def func_extract_sequencies_from_alignmnet(file_name):
    d_seqs={}
    infile=open(file_name).readlines()
    n=0
    while (infile[n]!="List of extant and reconstructed sequences\n" ):
        n+=1
    for i in range(41):
        tmp=infile[n+4+i].split()
        tmp_seq=''
        if (i<21):
            for j in tmp[1:]:
                tmp_seq=tmp_seq+j
            d_seqs[tmp[0]]=tmp_seq
        else:
            for j in tmp[2:]:
                tmp_seq=tmp_seq+j
            d_seqs[tmp[0]+'_'+tmp[1][1:]]=tmp_seq
    return d_seqs

def make_mut_list(seq,wtaa,length):
    tmp=[]
    if (len(seq)==length):
        for i in range(length):
            if (seq[i]!=wtaa[i]):
                tmp.append(str(i)+seq[i])
    return(tmp)


def func_extract_info_about_natural_and_library_variants(wt_data):
    natural_allele =pd.DataFrame(columns=map(lambda x: 'S'+str(x),range(1,13)),index=['nat','lib','nat_lib'])
    for s in range(1,13):
        natural_allele.loc['nat','S'+str(s)]=[]
        natural_allele.loc['lib','S'+str(s)]=[]
        wt_mask=wt_data.loc[s]
        pos1=int(wt_mask.iloc[0].loc['position'][1:])
        len_right=int(wt_mask.iloc[-1].loc['position'][1:])-pos1+1
        for i in range(0,len(wt_mask)):
            pos=int(wt_mask.iloc[i].loc['position'][1:])

            aa_list=list(wt_mask.iloc[i]['aa'])[0::2]
            for j in aa_list:
                natural_allele.loc['nat','S'+str(s)].append(str(pos-pos1)+j)

            aa_list=list(wt_mask.iloc[i]['lib'])[0::2]
            for j in aa_list:
                natural_allele.loc['lib','S'+str(s)].append(str(pos-pos1)+j)

            natural_allele.loc['nat_lib','S'+str(s)]=list(set(natural_allele.loc['lib','S'+str(s)])&set(natural_allele.loc['nat','S'+str(s)]))
    return natural_allele

def is_natural(mut_list,allele):
    for i in mut_list:
        if (i in allele):
            pass
        else:
            return 0
    return 1