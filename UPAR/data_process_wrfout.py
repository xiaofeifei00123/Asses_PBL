#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
将wrfout的格点数据插值到站点上
计算相关的诊断量
wrfout数据是个矩阵数据，没有按照气压来分布
所以用相同的气压坐标会有偏差
这种算法可能会带来巨大的偏差
TEMF方案的气压值和YSU的能一样吗
水平插值是对的
垂直插值可能存在问题
prc那块重新给就行

# 不同试验 
# 不同变量
# 不同站点
各个站点所在的经纬度不一样，各格点对应的pressure高度不一样，
所以对每个站点分别处理，比较慢
-----------------------------------------
Time             :2021/07/30 08:52:47
Author           :Forxd
Version          :1.0
'''
# %%
import xarray as xr
import numpy as np
import pandas as pd
import os
# from netCDF4 import Dataset
# from wrf import getvar, vinterp, interplevel
from data_process_main import GetData
from global_variable import station_dic
# from functools import reduce
from multiprocessing import Pool

# %%

def get_data_one_station(station):
    """多个模式的某一变量数据，统一读取
        所有模式，所有变量聚合
    """
    # model_dic = {}

    model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
    # model_list = ['ACM2', ]
    month = 'Jul'
    path_wrfout = '/mnt/zfm_18T/fengxiang/Asses_PBL/data/wrfout_data/'

    da_list = []
    for model in model_list:
        # print(model)
        ds = xr.Dataset()
        for var in ['temp', 'td', 'height_agl', 'U', 'V', 'pressure']:
            file_name = str(var) + "_" + str(month) + "_" + str(model) + "_latlon"
            flnm_var = os.path.join(path_wrfout, file_name)
            ds_var = xr.open_dataset(flnm_var)
            ds[var] = ds_var[var]
        ds_concat = ds.sel(lat=station['lat'], lon=station['lon'], method='nearest')
        ds_concat = ds_concat.rename({'height_agl':'height'})  # 保持变量名一致性
        # print(ds['height_agl'].max())
        # print(ds_concat)
        print("读完%s"%model)
        da_list.append(ds_concat.to_array())
    print("开始聚合不同模式的数据")
    da_return = xr.concat(da_list, pd.Index(model_list, name='model'))
    print("聚合完成")
    return da_return 


def get_one(station):
    pass
    pressure_level = np.arange(800,100,-1)
    da = get_data_one_station(station)
    # return da
    ### 给站点数据设置统一的气压维度
    dds = da.to_dataset(dim='variable')
    pressure_mean = dds['pressure'].mean(dim=['time', 'model']).squeeze()
    dds.drop_vars('pressure')
    ## 设置副轴
    dds = dds.assign_coords({'pressure':('bottom_top',pressure_mean.values)})
    dds_swap = dds.swap_dims({'bottom_top':'pressure'})
    
    ds_input = dds_swap.interp(pressure=pressure_level, kwargs={'fill_value':'extrapolate'})
    
    gd = GetData()
    ds_diagnostic = gd.caculate_diagnostic(ds_input)
    ds_return = xr.merge([ds_input, ds_diagnostic])
    da_return = ds_return.to_array()
    return da_return


## 测试  ##
# station = station_dic['DuLan']
# aa = get_one(station)
# aa
# aa.rename({'height_alg':'height'})
## 测试结束  ##



# %%
########## 单进程开始 ##########
# ds = xr.Dataset()
# for key in station_dic:
#     print("###读[%s]站点的数据"%key)
#     # gw = 
#     ds[key] = get_one(station_dic[key])

# dds = ds.where(ds.sel(variable='height')>0, drop=True)
########## 单进程结束  ##########

# %%

############# 多进程开始  ############
pool = Pool(8)
result = []

for key in station_dic:
    print("###读[%s]站点的数据"%key)
    tr = pool.apply_async(get_one, args=(station_dic[key],))
    result.append(tr)

pool.close()
pool.join()

ds_return = xr.Dataset()
for i,j in zip(result,station_dic):
        ds_return[j] = i.get()

ds_return = ds_return.where(ds_return.sel(variable='height')>0, drop=True)
# ds_return = ds_return.rename_vars({'height_agl':'height'})
print(ds_return)
ds_return.to_netcdf('/mnt/zfm_18T/fengxiang/Asses_PBL/data/wrfout_Jul_station2.nc')
##########  多进程结束  ##########
