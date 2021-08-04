#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
都处理成站点数据了，选取就可以了

# dimensions
站点: 
    GaiZe, ShiQuanhe, TingRi, ChangDu, LaSa, NaQu, LinZhi, ShenZha, 'TuoTuohe'
月: 
    May, Jul
变量: 
    Theta_v, q, td, temp, wind_s, rh
模式: 
    观测, fnl, micaps

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
import time

import sys
sys.path.append('/mnt/zfm_18T/fengxiang/Asses_PBL/Rain')  # 调用Rain目录下文件的读程序
from read_data import TransferData as rain_tr
from read_data import GetData as rain_gd

# from netCDF4 import Dataset
# from wrf import getvar, vinterp, interplevel

# from metpy.units import units
# from metpy.calc import specific_humidity_from_dewpoint
# from metpy.calc import mixing_ratio_from_specific_humidity
# from metpy.calc import virtual_potential_temperature
# from metpy.calc import potential_temperature
# from metpy.calc import relative_humidity_from_dewpoint
from global_variable import station_dic

class TransferMain():
    """这里父类就充当一个写公共属性和变量的类
    """
    pass
    def __init__(self, station, month, time_hour) -> None:
        pass
        self.station = station
        self.month = month
        self.time_hour = time_hour
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
        self.area = {"lat1": 24.875, "lat2": 40.125, "lon1": 69.875, "lon2": 105.125}

class TimeProcess(TransferMain):
    pass

    def get_time_index_rain(self,):
        """获得有降水的时刻
        随站点，时次和月份而变化
        """

        gd = rain_gd(self.month)
        rain = gd.get_rain_hourly()  # 所有时次，所有格点的降水数据
        
        time_index = rain['obs'].time.sel(time=datetime.time(int('12')))  
        rain = rain.sel(time=time_index)
        time_index_12 = rain['obs'].time.sel(time=datetime.time(int('12')))  
        time_index_00 = rain['obs'].time.sel(time=datetime.time(int('00')))  
        time_index = np.union1d(time_index_00.values, time_index_12.values)
        rain = rain.sel(time=time_index)

        # 根据格点降水数据得到站点降水
        tr = rain_tr(ds=rain, area=self.area, time_flag='all')
        # station = station_dic['TingRi']  # 不同站点有降水的时次是不一样的
        dd = tr.rain_station(self.station)
        if dd['obs'].values.max() < 0.1:
            return None
        # os.system('pause')
        ## 满足条件的值保留，不满足条件的赋值为np.nan
        cc = dd.where(dd['obs']>0.1 , drop=True)   # 站点数据>0.1的时次

        time_index_return = cc.time.values
        return time_index_return


    def get_time_index_upar(self, condition, var):
        pass
        """包括了降水的时间和实际数据有的时间
        时间筛选的条件:
            降水
            upar(micaps)
            upar(TEMF)
        """
        ## 该站点所有的00时和12时中,降水>0.1的时次
        time_index_rain = self.get_time_index_rain()  

        ## 得到探空的数据,这里写任意一个变量都可以(主要看这个时次有无数据), 以q代替
        tr = TransferData(self.station, self.month, self.time_hour, var)
        # tr = TransferData(self.station, self.month, self.time_hour, 'q')
        ds = tr.get_data_one()  # 该站点的各模式探空数据, 多时次
        time_index_upar = ds.time.values  # 已经是筛选出来的Upar数据的共有时间

        ## 如果没有探空数据, 直接返回None
        if time_index_upar is None:
            return None  # 不管什么条件，都没有数据返回

        ## 有降水时间的条件
        if not time_index_rain is None:
            pass
            time_index_precip = np.intersect1d(time_index_rain,
                                               time_index_upar)
            if condition == 'rain':
                pass
                time_index_return = time_index_precip
                ## 返回有降水的条件下的探空数据时次
            elif condition == 'dry':
                pass
                time_index_return = np.setdiff1d(time_index_upar,
                                              time_index_precip)
            else:
                ## 所有的时次
                time_index_return = time_index_upar
        else:
            ## 无降水时次
            pass
            if condition == 'rain':
                pass
                time_index_return = None
                ## 返回有降水的条件下的探空数据时次
            elif condition == 'dry':
                pass
                time_index_return = time_index_upar
            else:
                ## 所有的时次
                time_index_return = time_index_upar

        return time_index_return
        

        ## 将降水和探空资料的时间取交集
        # time_index_rain = time_index_rain.values
        # time_index_return = np.intersect1d(time_index_upar, time_index_rain)
        # return time_index_return  # 包括00时、12时、甚至是06时在内的所有数据


class TransferData(TransferMain):
    """将获得的数据传递出去
    原有的传递思路没变, 
    提供站点
    """
    def __init__(self, station, month, time_hour, var) -> None:
        pass
        self.var = var
        super().__init__(station, month, time_hour)

    def get_data_one(self):
        """获取单个原始变量的值, 比如temp各试验+观测的值
            这里没有进行nan值的丢弃
        """
        ## 获取站点数据
        da_wrfout = xr.open_dataset(self.path_wrfout)[self.station['name']]
        da_micaps = xr.open_dataset(self.path_micaps)[self.station['name']]
        da_fnl = xr.open_dataset(self.path_fnl)[self.station['name']]
        if self.station['name'] in ['GaiZe', 'ShenZha', 'ShiQuanhe']:
            da_gps = xr.open_dataset(self.path_gps)[self.station['name']]
            da_list = [da_wrfout, da_micaps, da_gps, da_fnl]
        else:
            da_list = [da_wrfout, da_micaps, da_fnl]

        ## 取时间交集, 初步筛选时间, 这里是将矩阵时间不一致的给灭了
        # time_index = ds_gps.time.values
        time_index = da_micaps.time.values
        for da in da_list:
            time_index1 = da.time.values
            time_index = np.intersect1d(time_index, time_index1)
        ## 三个数据都有的时间交集

        station_name = self.station['name']
        da_wrfout = da_wrfout.sel(time=time_index, variable=self.var)
        da_micaps=da_micaps.sel(time=time_index, variable=self.var)
        da_fnl = da_fnl.sel(time=time_index, variable=self.var)

        if self.station['name'] in ['GaiZe', 'ShenZha', 'ShiQuanhe']:
            da_gps = da_gps.sel(time=time_index, variable=self.var)
            dds = xr.concat([da_micaps, da_fnl, da_gps], 
                            pd.Index(['micaps', 'fnl', 'gps'], name='model'))
        else:
            dds = xr.concat([da_micaps, da_fnl], 
                            pd.Index(['micaps', 'fnl'], name='model'))
        ds_return = xr.concat([da_wrfout, dds], dim='model')
        # print(ds_return)
        ds_return = ds_return.dropna(dim='time', how='any')
        ds_return = ds_return.to_dataset(dim='model')
        return ds_return   # 返回DataSet

    def transfer_data(self, condition):
        """传递数据
        """
        tp = TimeProcess(self.station, self.month, self.time_hour)
        # time_index_rain = tp.get_time_index_upar(condition)  # 最后的既有降水, 观测和模式共有数据的时间
        time_index_rain = tp.get_time_index_upar(condition, self.var)  # 最后的既有降水, 观测和模式共有数据的时间
        # if time_index_rain is None:  # time_index_rain是一个numpy数组
        #     return None

        var_ds = self.get_data_one()  # 各模式数据的站点ds
        # ## 所有时次的探空数据
        # ds_all = var_ds.sel(time=datetime.time(int(self.time_hour))) 
        ## 有降水的时次的探空数据
        if not time_index_rain is None:  # time_index_rain是一个numpy数组

            # time_index_var = var_ds.time.values            
            # time_index = np.intersect1d(time_index_rain, time_index_var)
            # ds_return = var_ds.sel(time=time_index)

            ds_return = var_ds.sel(time=time_index_rain)
            ds_return = ds_return.sel(
                time=datetime.time(int(self.time_hour)))  
            tt_time = ds_return.time.values
            if len(tt_time) == 0:
                return None
            else:
                print("最后筛选的时间是  %s"%str(ds_return.time.values))
                return ds_return
        else:
            print("没有符合条件的时间")
            return None

        



if __name__ == '__main__':
    pass
    ## 申扎站7月弃用
    # %%
    month = 'Jul'
    station = station_dic['LaSa']
    time_hour = '12'
    
    # %%
    tr = TransferData(station, month, time_hour, 'wind_s')
    bb = tr.transfer_data('all')
    print(bb)
    
