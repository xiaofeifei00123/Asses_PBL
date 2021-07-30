#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读micpas的探空资料
-----------------------------------------
Time             :2021/07/26 09:48:05
Author           :Forxd
Version          :1.0
'''

import pandas as pd
import numpy as np
import xarray as xr

# %%
flnm = '/mnt/zfm_18T/fengxiang/DATA/UPAR/Upar_2016/16050108.000'

col_names = ['pressure', 'height', 't', 'td', 'wind_direct', 'wind_speed']
# %%

with open(flnm, 'r') as f:
    whole_text = f.readlines()
    




# df = pd.read_table(
#     flnm,
#     sep='\\s+',
#     skiprows=1,
#     usecols=[0,1,2,3,4,5],
#     names=col_names,
# )

# # %%
# print(df)


# %%
