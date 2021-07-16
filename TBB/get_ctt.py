#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
画月平均降水
-----------------------------------------
Time             :2021/05/29 12:31:34
Author           :Forxd
Version          :1.0
'''

import xarray as xr
import cmaps
import numpy as np
import xesmf as xe
import os
import wrf
import netCDF4 as nc


class Summer_mean:
    
    def __init__(self,):
        # self.path = '/mnt/zfm/TPV_statistics/ERA5_V_850hpa/'
        # self.year = year
        pass

    def month_mean(self,flnm):
        # filein = self.path+'ERA5_V_850_'+str(self.year)+'_0'+str(month)+'.nc'
        ds = xr.open_dataset(flnm)
        # r = ds.RAINNC
        r = ds.ctt

        bb = ds.mean(dim='time')
        # print(bb.values)
        
        # return bb.values
        return bb.ctt
    
def get_total_rain():

    model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','MYJ','YSU']
    path = '/home/fengxiang/Project/Asses_pbl_July/Data/Process/'
    file_list = []
    for i in model_list:
        # flnm = path + "RAINNC_Jul_"+i+'_latlon'
        flnm = path + "ctt_Jul_"+i+'_latlon'
        file_list.append(flnm)

    rain = {}   # 各模式降水和
    su = Summer_mean()
    # for fl in file_list[0:1]:
    for i in range(len(file_list)):
        # ds = xr.open_dataset(fl)
        fl = file_list[i]
        # print(fl)
        rain_total = su.month_mean(fl)
        rain[model_list[i]] = rain_total
    # print(rain['QNSE'])
    return rain  # 返回的是一个各试验tbb和的列表



if __name__ == '__main__':

    ctt = get_total_rain()
    print(type(ctt))
    # ds = xr.Dataset(ctt)
    # ds.tonetcdf('/home/fengxiang/Project/Asses_pbl_July/Data/ctt.nc')
    # print(ds)

    # # model_list = ['MYJ', 'QNSE', 'QNSE_EDMF', 'TEMF', 'YSU']
    # model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','MYJ','YSU']
    # path = '/home/fengxiang/Project/Asses_pbl_July/Data/Process/'
    # file_list = []
    # for i in model_list:
    #     flnm = path + "RAINNC_Jul_"+i+'_latlon'
    #     file_list.append(flnm)

    # rain = {}   # 各模式降水和
    # su = Summer_mean()
    # # for fl in file_list[0:1]:
    # for i in range(len(file_list)):
    #     # ds = xr.open_dataset(fl)
    #     fl = file_list[i]
    #     # print(fl)
    #     rain_total = su.month_mean(fl)
    #     rain[model_list[i]] = rain_total
    # print(rain['QNSE'])
        
