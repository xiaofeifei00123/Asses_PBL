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

# %%
from typing import Pattern
import pandas as pd
import numpy as np
import xarray as xr
import re
import sys
import os
from io import StringIO
from global_variable import station_dic

from data_process_main import GetData

class GetMicaps(GetData):
    
    def __init__(self, station_number):
        pass
        self.station_number = station_number
        self.path_micaps = '/mnt/zfm_18T/fengxiang/DATA/UPAR/Upar_2016/'
        self.pressure_level = np.arange(570, 280, -5)

    def read_data_once(self, flnm):
        """根据初始的txt文件
        筛选出需要的站点数据
        """
        with open(flnm, 'r', encoding='gbk') as f:
            whole_text = f.read()
        pattern = '[0-9]+\s+\S+\.\d+\s+\S+\.\d+\s+\d+\s+\d+'
        data = re.split(pattern, whole_text)
        dd = data.pop(0)  # 去掉首行的信息
        seg_data = re.findall(pattern, whole_text)
        length = len(seg_data)
        station = seg_data[0][0:4]

        for i in range(length):
            # if 
            station = seg_data[i]
            station_number = station.split()[0]
            # print(station)
            if str(station_number) == str(self.station_number):
                pass
                df_data = data[i]
                data_file = StringIO(df_data)  # 将这一个站点的数据，模拟到文件中去
                return data_file

    def transpose_data_once(self, data_file):
        """将站点数据文件转换为df格式, 
        再转为DataArray
        data_file, 字符串生成的虚假文件
        """
        # data_file = self.read_data_once()
        col_names = ['pressure', 'height', 'temp', 'td', 'wind_d', 'wind_s']
        df = pd.read_table(
            data_file,  # 由字符串虚拟的文件
            sep='\\s+',
            skiprows=1,
            usecols=[0,1,2,3,4,5],
            names=col_names,
        )
        df = df.where(df < 9999, np.nan)  # 将缺省值赋值为NaN
        df = df.set_index(['pressure'])  # 将pressure这一列设为index
        df.columns.name = 'variable'
        da = xr.DataArray(df)
        return da

    def data_micaps(self, ):
        
        aa = os.listdir(self.path_micaps)  # 文件名列表
        aa.sort()  # 排一下顺序，这个是对列表本身进行操作
        da_time = []  # 每个时次变量的列表
        ttt = []  # 时间序列列表

        for flnm in aa:
            fl_time = '20'+flnm[0:-4]
            # print(fl_time)
            tt = pd.to_datetime(fl_time, format='%Y%m%d%H')
            tt = tt - pd.Timedelta(hours=8)
            ## 这时间是不规则的
            file_name = os.path.join(self.path_micaps, flnm)
            # print(file_name)
            data_file = self.read_data_once(file_name)
            if not data_file:
                # return None
                continue
            ttt.append(tt)
            da = self.transpose_data_once(data_file)
            dda = da.interp(pressure=self.pressure_level)
            da_time.append(dda)  # 很多时次都是到595hPa才有值, 气压和高度的对应关系会随着时间发展而变化, 气压坐标和高度坐标不能通用
        da_return = xr.concat(da_time, pd.Index(ttt, name='time'))
        return da_return

if __name__ == '__main__':
    pass
    # %%
    ds_nc = xr.Dataset()
    for key in station_dic:

        print("读取 %s 站的数据"%key)
        station_number = station_dic[key]['number']
        # print(station_number)
        gd = GetMicaps(station_number=station_number)    
        da = gd.data_micaps()

        da = da.transpose(*(...,'pressure'))
        ds = da.to_dataset(dim='variable')
        ds_diagnostic = gd.caculate_diagnostic(ds)
        ## 将原来的变量和计算的诊断变量合并为一个DataArray
        ds_return = xr.merge([ds, ds_diagnostic])
        dda = ds_return.to_array()
        ## 不同站点的数据组合为一个Dataset
        ds_nc[key] = dda

    ds_nc.to_netcdf('/mnt/zfm_18T/fengxiang/DATA/UPAR/upar_2016_all_station.nc')

    #-------------------------------
    ## 测试
    #-------------------------------