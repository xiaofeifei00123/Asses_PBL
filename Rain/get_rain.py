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
        r = ds.RAINNC
        # print(ds.time.values)
        r_list = []
        for i in range(15):
            j = i*48
            k = (i+1)*48
            # print(j,k)
            rr = r.loc[j:k-1,:,:]  # 每一个试验区段降水

            rrr = rr[-1] - rr[0]
            r_list.append(rrr)
        rain_sum= sum(r_list)
        return rain_sum
    
def get_total_rain():

    model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','MYJ','YSU']
    path = '/home/fengxiang/Project/Asses_pbl_July/Data/Process/'
    file_list = []
    for i in model_list:
        flnm = path + "RAINNC_Jul_"+i+'_latlon'
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
    return rain



if __name__ == '__main__':

    rain = get_total_rain()
    print(rain)

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
        
