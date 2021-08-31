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
from netCDF4 import Dataset
from wrf import getvar, vinterp, interplevel

from data_process_main import GetData
from global_variable import station_dic

from functools import reduce

from multiprocessing import Pool

# %%

def get_data_one_station(station):
    """多个模式的某一变量数据，统一读取
        所有模式，所有变量聚合
    """
    # model_dic = {}

    model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
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
    return da
    # da_swap = da.swap_dims({'bottom_top':'pressure'})
    # da_interp = da_swap.interp(pressure=pressure_level, kwargs={'fill_value':'extrapolate'})
    # ds_input = da_interp.to_dataset(dim='variable')
    # gd = GetData()
    # ds_diagnostic = gd.caculate_diagnostic(ds_input)
    # ds_return = xr.merge([ds_input, ds_diagnostic])
    # da_return = ds_return.to_array()
    # return da_return


## 测试
station = station_dic['DuLan']
aa = get_one(station)

# %%
ddds = aa.to_dataset(dim='variable')
ddds['height_agl'].max()

# %%
# ddds['height_agl'].max()


# %%



# %%
    # print(dd)
## 单进程
ds = xr.Dataset()
for key in station_dic:
    print("###读[%s]站点的数据"%key)
    # gw = 
    ds[key] = get_one(station_dic[key])


# %%

### 多进程
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
    # print(ds_nc)

print(ds_return)
ds_return.to_netcdf('/mnt/zfm_18T/fengxiang/Asses_PBL/data/wrfout_Jul_station2.nc')

# %%
# ds_return




# %%
# ddd['td'].max()                
# ds
# ddd

# %%
                
                

                # model_dic[model] = self.get_data_single_once(var, model)
                # ds[model] = self.get_data_single_once(var, model)

        #     elif var == 'wind_s':
        #         U = self.get_data_single_once('U', model)
        #         V = self.get_data_single_once('V', model)
        #         # model_dic[model] = t - td
        #         ds[model] = xr.ufuncs.sqrt(U**2 + V**2)
        # return ds
        

# ds_station = xr.Dataset()
# for key in station_dic:
# station = station_dic['GaiZe']
#     station = station_dic[key]
# ds_var = get_station_one(station)
#     # print(ds_var)
#     ds_var_diag = gd.caculate_diagnostic(ds_var)
#     ds_var_return = xr.merge([ds_var, ds_var_diag])
#     da_var = ds_var_return.to_array()
#     ds_station[station['name']] = da_var

# # %%
# flnm_save = '/mnt/zfm_18T/fengxiang/Asses_PBL/data/'+'wrfout_'+str(month)+"_station.nc"
# # ds_station.to_netcdf(flnm_save)
        
        
        
        
        
        
        
# %%
        


class SaveData():
    
    def get_station_one(self, station):
        """获得一个站点的所有数据"""
        gw = GetWrfout(station, month)
        var_list = ['temp', 'td', 'wind_s']

        ds_var = xr.Dataset()  # 不同变量的聚合
        for var in var_list:
            ds_model = gw.get_data_var(var)  # 不同模式的
            da_model = ds_model.to_array()
            da_model = da_model.rename({'variable':'model'}) # 完成不同模式的聚合
            # print(da_model)
            ds_var[var] = da_model

        # time_index = reduce(np.intersect1d,
        # time_index_temp = ds_var['temp'].sel(model='TEMF').dropna(dim='time', how='all').time.values 
        # time_index_td = ds_var['td'].sel(model='TEMF').dropna(dim='time', how='all').time.values 
        # time_index_wind_s = ds_var['wind_s'].sel(model='TEMF').dropna(dim='time', how='all').time.values 

        # time_index1 = np.intersect1d(time_index_temp, time_index_td)
        # time_index = np.intersect1d(time_index1, time_index_wind_s)
        # ds_return = ds_var.sel(time=time_index)
        ds_return = ds_var
        return ds_return

    def save_station_nc(self, month):
        """将不同站点的数据保存为一个文件
        """
            
        ds_station = xr.Dataset()
        gd = GetData()
        for key in station_dic:
            station = station_dic[key]
            ds_var = self.get_station_one(station)
            print("读[%s]站的数据" %key)
            ds_var_diag = gd.caculate_diagnostic(ds_var)
            ds_var_return = xr.merge([ds_var, ds_var_diag])
            da_var = ds_var_return.to_array()
            ds_station[station['name']] = da_var

        flnm_save = '/mnt/zfm_18T/fengxiang/Asses_PBL/data/'+'wrfout_'+str(month)+"_station.nc"
        ds_station.to_netcdf(flnm_save)

# if __name__ == '__main__':
    
#     # for month in ['May', 'Jul']:
#     for month in ['Jul']:
#         sd = SaveData()
#         sd.save_station_nc(month)
    
    


    # %%
        
    # ds_station = xr.Dataset()
    # gd = GetData()
    # for key in station_dic:
    #     station = station_dic[key]
    #     ds_var = get_station_one(station)
    #     # print(ds_var)
    #     ds_var_diag = gd.caculate_diagnostic(ds_var)
    #     ds_var_return = xr.merge([ds_var, ds_var_diag])
    #     da_var = ds_var_return.to_array()
    #     ds_station[station['name']] = da_var

    # # %%
    # flnm_save = '/mnt/zfm_18T/fengxiang/Asses_PBL/data/'+'wrfout_'+str(month)+"_station.nc"
    # # ds_station.to_netcdf(flnm_save)
    


    # --------------------------------------
    #####          测试
    # --------------------------------------
    # %%
    # month = 'Jul'
    # station = station_dic['GaiZe']
    # # %%
    # gw = GetWrfout(station, month)
    # # da = gw.get_data_single_once('td', 'TEMF')
    # ds = get_station_one(station)
    # # %%
    # # da = ds_model['TEMF']


# %%
