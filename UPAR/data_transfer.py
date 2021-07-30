#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
都处理成站点数据了，选取就可以了

# 需要什么样的数据
1. 给定的
站点
月
变量
2. 变化的
模式
观测
fnl
micaps


# 对时间维度的处理
分别有哪些时间维度的数据






-----------------------------------------
Time             :Mon Jun 28 20:02:27 CST 2021
Author           :Forxd
Version          :0.3
'''

import xarray as xr
import os
import numpy as np
import pandas as pd

import datetime

# from netCDF4 import Dataset
# from wrf import getvar, vinterp, interplevel

# from metpy.units import units
# from metpy.calc import specific_humidity_from_dewpoint
# from metpy.calc import mixing_ratio_from_specific_humidity
# from metpy.calc import virtual_potential_temperature
# from metpy.calc import potential_temperature
# from metpy.calc import relative_humidity_from_dewpoint
from global_variable import station_dic



class TransferData():
    """将获得的数据传递出去
    """
    def __init__(self, station, month, var) -> None:
        pass
        self.station = station
        self.month = month
        self.var = var
        self.path_fnl = '/mnt/zfm_18T/fengxiang/Asses_PBL/data/fnl_station.nc'
        self.path_micaps = '/mnt/zfm_18T/fengxiang/DATA/UPAR/upar_2016_all_station.nc'
        if self.month== 'Jul':
            self.path_gps = '/mnt/zfm_18T/fengxiang/DATA/UPAR/GPS_Upar_2016_Jul.nc'
            self.path_wrfout = '/mnt/zfm_18T/fengxiang/Asses_PBL/data/wrfout_Jul_station.nc'
        elif self.month == 'May':
            self.path_gps = '/mnt/zfm_18T/fengxiang/DATA/UPAR/GPS_Upar_2016_May.nc'
            self.path_wrfout = '/mnt/zfm_18T/fengxiang/Asses_PBL/data/wrfout_May_station.nc'
        else:
            print("输入的月份有误")

    def get_data_one(self, var):
        """获取单个原始变量的值, 比如temp各试验+观测的值
        """
        ds_wrfout = xr.open_dataset(self.path_wrfout)
        ds_micaps = xr.open_dataset(self.path_micaps)
        ds_gps = xr.open_dataset(self.path_gps)
        ds_fnl = xr.open_dataset(self.path_fnl)

        ds_list = [ds_wrfout, ds_micaps, ds_gps, ds_fnl]
        ## 取时间交集, 初步筛选时间
        time_index = ds_gps.time.values
        for ds in ds_list:
            time_index1 = ds.time.values
            time_index = np.intersect1d(time_index, time_index1)

        station_name = self.station['name']
        ds_wrfout = ds_wrfout[station_name].sel(time=time_index, variable=var)
        ds_micaps=ds_micaps[station_name].sel(time=time_index, variable=var)
        ds_fnl = ds_fnl[station_name].sel(time=time_index, variable=var)
        ds_gps = ds_gps[station_name].sel(time=time_index, variable=var)


        dds = xr.concat([ds_micaps, ds_fnl, ds_gps], pd.Index(['micaps', 'fnl', 'gps'], name='model'))
        ds_return = xr.concat([ds_wrfout, dds], dim='model')

        print(ds_return)

        return ds_return






        
        
        
        


    def transfer_data(self, var):
        """传递数据、存储数据
        对数据进行统一插值处理
        """

        if var in ['temp', 'td',  'wind_s', 'q', 'theta_v', 'rh']:
            var_ds = self.get_data_one(var)
        else:
            var_ds = None

        return var_ds
        
        



if __name__ == '__main__':

    # %%
    # station = station_dic['TingRi']
    station = station_dic['ShenZha']
    tr = TransferData(station, 'May', 'theta_v')
    # %%
    # print(aa)
    ttt = tr.transfer_data('wind_s')  # 多试验的某变量(temp, rh..)数据
    # print(model_dic)  
    # # TODO 别传model_dic了，全部换成ds
    # # %%
    # aa = model_dic['obs']
    # # aa = model_dic['fnl']
    # for i in range(len(aa.time)):
    #     print(aa.isel(time=i).values)
    # # print(aa)
# %%

# %%
