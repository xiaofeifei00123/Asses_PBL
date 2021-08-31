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
import xarray as xr
import numpy as np
import pandas as pd
import os
from netCDF4 import Dataset
from wrf import getvar, vinterp, interplevel

from data_process_main import GetData
from global_variable import station_dic

from functools import reduce

class GetWrfout(GetData):

    def __init__(self, station, month):
        self.station = station
        # self.pressure_level = np.arange(610, 100, -5)
        self.pressure_level = np.arange(570, 280, -5)
        self.path = '/mnt/zfm_18T/fengxiang/Asses_PBL/'
        self.path_wrfout = '/mnt/zfm_18T/fengxiang/Asses_PBL/data/wrfout_data/'
        self.model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        self.month = month

        if self.month == 'Jul':
            self.month_num = '07'
            self.time_first = '2016-07-01 13:00'
            # self.flnm_height_obs = '/mnt/zfm_18T/fengxiang/DATA/GPS_Upar_2016/SCEX_TIPEX3_UPAR_GPS_MUL_55228-201607/upar_G_55228_2016070206.txt'
            # self.rh_file = '/mnt/zfm_18T/fengxiang/DATA/FNL/FNL_2016/fnl_rh_201607'
            # self.fnl_file = '/mnt/zfm_18T/fengxiang/DATA/FNL/fnl_201607.nc'
        elif self.month == 'May':
            self.month_num = '05'
            self.time_first = '2016-05-01 13:00'
            # self.flnm_height_obs = '/mnt/zfm_18T/fengxiang/DATA/GPS_Upar_2016/SCEX_TIPEX3_UPAR_GPS_MUL_55228-201605/upar_G_55228_2016051112.txt'
            # self.rh_file = '/mnt/zfm_18T/fengxiang/DATA/FNL/FNL_2016/fnl_rh_201605'
            # self.fnl_file = '/mnt/zfm_18T/fengxiang/DATA/FNL/fnl_201605.nc'
        else:
            print("%s这个月份不在数据集内"%self.month)
            
    """
    获取wrfout数据
    """
    def regrid_wrfout(self, da_temp):
        """对数据进行水平插值和垂直插值,
        先给数据添加prc属性，再插值到特定层上
        """
        def get_pressure_lev():
            # 得到各层的pressure值
            # 不同时刻各层气压值，差别可以忽略不计,
            # 后面还要对气压层进行插值, 这里不对它做过高精度要求
            # 不同的方案和时间，在同一站点，各层气压值相差小于1度
            # 故不作分开考虑
            # path = '/mnt/zfm_18T/Asses_PBL/wrfout_data/'
            # path = os.path.join(self.path, '/data/wrfoutdata')
            path = self.path_wrfout
            # flnm_pressure = os.path.join(path, 'pressure_Jul_YSU_latlon')
            pressure_name = 'pressure_'+self.month+'_YSU_latlon'
            flnm_pressure = os.path.join(path, pressure_name)
            # flnm_pressure = os.path.join(path, 'pressure_Jul_YSU_latlon')
            ds_pressure = xr.open_dataset(flnm_pressure)
            pr = ds_pressure.pressure
        
            # prb = pr.sel(time='2016-07-01 13:00')
            prb = pr.sel(time=self.time_first)
            lat = self.station['lat']
            lon = self.station['lon']
            # prc = prb.sel(lat=32.13, lon=92.5, method='nearest')
            prc = prb.sel(lat=lat, lon=lon, method='nearest')
            # prc = prc[0] - prc  # 离地气压高度
            return prc


        def regrid():
            """对wrfout数据进行插值的
            需要水平插值和垂直插值两项
            Args:
                da_temp (DataArray): 需要插值的变量

            Returns:
                DataArray: 插值后的变量
            """
            # 将bottom_top坐标换成气压坐标和高度坐标
            time_coord = da_temp.time.values
            lat_coord = da_temp.lat.values
            lon_coord = da_temp.lon.values
            prc = get_pressure_lev()
            pressure_coord = prc.values
            da_temp_reset = da_temp.values
            da = xr.DataArray(
                da_temp_reset,
                coords=[time_coord, pressure_coord, lat_coord, lon_coord],
                dims=['time', 'pressure', 'lat', 'lon'])
            # 水平插值
            da = da.sel(lat=self.station['lat'],
                        lon=self.station['lon'],
                        method='nearest')
            # 垂直插值
            da_return = da.interp(pressure=self.pressure_level)
            return da_return
        return regrid()

    def get_data_single_once(self, var, model):
        """读一个模式一个变量的数据
            模块尽可能的小
        """
        file_name = str(var) + "_" + str(
            self.month) + "_" + str(model) + "_latlon"
        flnm_var = os.path.join(self.path_wrfout, file_name)
        ds_var = xr.open_dataset(flnm_var)
        da_var = ds_var[var]
        da_var = self.regrid_wrfout(da_var)
        # da_return = da_var.dropna(dim='time',how='all')
        return da_var

    def get_data_var(self, var):
        """多个模式的某一变量数据，统一读取
        """
        # model_dic = {}
        ds = xr.Dataset()
        for model in self.model_list:
            if var in ['temp', 'td', 'height_agl']:
                # model_dic[model] = self.get_data_single_once(var, model)
                ds[model] = self.get_data_single_once(var, model)

            # elif var == 't_td':
            #     t = self.get_data_single_once('temp', model)
            #     td = self.get_data_single_once('td', model)
            #     # model_dic[model] = t - td
            #     ds[model] = t - td

            elif var == 'wind_s':
                U = self.get_data_single_once('U', model)
                V = self.get_data_single_once('V', model)
                # model_dic[model] = t - td
                ds[model] = xr.ufuncs.sqrt(U**2 + V**2)
        return ds
        

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

if __name__ == '__main__':
    
    # for month in ['May', 'Jul']:
    for month in ['Jul']:
        sd = SaveData()
        sd.save_station_nc(month)
    
    


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
