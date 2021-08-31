#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
不同区域月总降水量平均值的
柱状图
目的是为了验证数据的正确性
-----------------------------------------
Time             :2021/08/29 09:31:42
Author           :Forxd
Version          :1.0
'''
# %%
import xarray as xr
from read_data import TransferData, GetData
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler
from global_variable import station_dic
from data_process_landause import get_landmask
# import datetime
from draw_time_station_land import Data_process
import seaborn as sns


# %%
##### 处理数据
area = {"lat1": 24.875, "lat2": 40.125, "lon1": 69.875, "lon2": 105.125}
month = 'Jul'
time_flag = 'all'
gd = GetData(month)
rain = gd.get_rain_hourly()
dp = Data_process(rain)
d3 = dp.get_rain_land_area()
d3_sum = d3.sum(dim='time')
ds = d3_sum
# %%
# %%
##### 画图
# da = d3_sum.sel(landtype='grass').to_array()
fig = plt.figure(figsize=(14, 9), dpi=300)  # 创建页面
ax = fig.add_axes([0.1,0.1, 0.8,0.8])
land_type = ds['landtype'].values
labels = ds.data_vars
# x = np.arange(len(labels))
x = np.arange(len(land_type))
width = 0.14
# da = ds.sel(landtype='grass').to_array()
# ds.sel(landtype='grass').to_array().values
# %%
# labels.values
# for i in labels:
    # print(i)
# len(land_type)
# x

# %%
# da
# ds['QNSE']

rects1 = ax.bar(x-width*2, ds['obs'], width, label='OBS')
rects2 = ax.bar(x-width*1, ds['ACM2'], width, label='ACM2')
rects3 = ax.bar(x, ds['YSU'], width, label='YSU')
rects4 = ax.bar(x+width*1, ds['QNSE'], width, label='QNSE')
rects5 = ax.bar(x+width*2, ds['QNSE_EDMF'], width, label='QNSE_EDMF')
rects6 = ax.bar(x+width*3, ds['TEMF'], width, label='TEMF')



# rects1 = ax.bar(x-width, ds.sel(landtype='grass').to_array().values, width, label='grass')
# rects2 = ax.bar(x, ds.sel(landtype='bush').to_array().values, width, label='bush')
# rects3 = ax.bar(x+width, ds.sel(landtype='bare').to_array().values, width, label='bare')
ax.bar_label(rects1, padding=0, fmt='%d')
ax.bar_label(rects2, padding=0, fmt='%d')
ax.bar_label(rects3, padding=0, fmt='%d')
ax.bar_label(rects4, padding=0, fmt='%d')
ax.bar_label(rects5, padding=0, fmt='%d')
ax.bar_label(rects6, padding=0, fmt='%d')

# labels = ['OBS' if i == 'obs' else i for i in labels]
ax.set_xticks(x)
ax.set_xticklabels(land_type, fontsize=15)
ax.set_yticks(np.arange(0,351,50))
ax.set_yticklabels(np.arange(0,351,50), fontsize=15)
ax.set_ylim(0,360)
ax.legend(fontsize=16, edgecolor='white', loc='upper right')
plt.show()
fig.savefig('/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/picture/rain_bar.png')