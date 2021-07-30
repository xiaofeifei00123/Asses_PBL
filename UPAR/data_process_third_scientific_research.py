#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
处理第三次科考的GPS探空数据
这些数据是按站点给的
-----------------------------------------
Time             :2021/07/29 11:30:09
Author           :Forxd
Version          :1.0
'''
from metpy.calc.thermo import dewpoint
import pandas as pd
import numpy as np
import xarray as xr
import os

from metpy.units import units
from metpy.calc import specific_humidity_from_dewpoint
from metpy.calc import mixing_ratio_from_specific_humidity
from metpy.calc import virtual_potential_temperature
from metpy.calc import potential_temperature
from metpy.calc import relative_humidity_from_dewpoint
from xarray.core import variable

from global_variable import station_dic
from data_process_main import GetData


class GetThirdScientific(GetData):
    """获取观测数据, 几种变量统一读取
    """
    def __init__(self, station, month):
        self.station = station
        # self.pressure_level = np.arange(610, 100, -5)
        self.pressure_level = np.arange(570, 280, -5)
        self.path = '/mnt/zfm_18T/fengxiang/DATA/UPAR'
        self.path_wrfout = '/mnt/zfm_18T/fengxiang/Asses_PBL/data/wrfout_data/'
        self.model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        self.month = month

        # self.month = 'May'
        if self.month == 'Jul':
            self.month_num = '07'
            self.time_first = '2016-07-01 13:00'
        elif self.month == 'May':
            # self.month = 'May'
            self.month_num = '05'
            self.time_first = '2016-05-01 13:00'
        else:
            print("%s这个月份不在数据集内"%self.month)

    def read_single(self, flnm):
        """
        读取单个文件的观测资料
        因为变量的缺省不同，所以这里是分别对每个变量进行插值的
        这里出现的这些问题，都是由于自己对pandas库不熟悉所导致的，
        读取单个时次的观测数据
        """
        col_names = [
            'pressure', 'height', 'temp', 'td', 't_td', 'wind_d', 'wind_s'
        ]
        ## 按列数据, 要哪几列的数据, 再命名
        df = pd.read_table(
            flnm,
            sep=' ',
            #  skiprows=0,
            usecols=[26, 27, 29, 30, 31, 32, 33],
            names=col_names,
        )
        df = df.where(df < 9999, np.nan)  # 将缺省值赋值为NaN
        df = df.sort_values('pressure', ascending=False)  # 按照某一列排序
        
        # -------------------------------------------------------
        # 将含有缺省值的行删掉
        df1 = df.dropna(axis=0, subset=['pressure'])
        df1 = df.dropna(axis=0, subset=['temp'])
        # print(df1)
        da_list = []
        var_list = ['temp', 'td', 't_td', 'wind_s']
        for var in var_list:
            df2 = df1  # 每一次重新传入这个df1,尽可能多的保留数据
            df2 = df2.dropna(axis=0, subset=[var])
            # 坐标处理
            df2 = df2.drop_duplicates('pressure', 'first')  # 返回副本
            # print(df2)
            # 改变pressure坐标为相对的还是绝对的
            df2 = df2.set_index(['pressure'])  # 将pressure这一列设为index
            df3 = df2
            # print(df3)

            # -------------------------------------------------------
            # 对某些时次td没有值的，补充NaN值
            if var in ['td', 't_td']:
                if df3['t_td'].size == 0:
                    da = xr.DataArray([np.nan, np.nan],
                                      coords=[[0, 1]],
                                      dims='pressure')
                    df3['t_td'] = da.to_series()
                elif df3['td'].size == 0:
                    da = xr.DataArray([np.nan, np.nan],
                                      coords=[[0, 1]],
                                      dims='pressure')
                    df3['td'] = da.to_series()
            # -------------------------------------------------------
            ## 垂直坐标处理, 允许有缺省值的出现
            da = xr.DataArray.from_series(df3[var])
            dda = da.interp(pressure=self.pressure_level)
            da_list.append(dda)
        # da_return = xr.concat(da_list, dim=var_list)
        # da_return = xr.concat(da_list, dim=var_list)
        da_return = xr.concat(da_list, pd.Index(var_list, name='variable'))
        return da_return

    def read_obs(self, ):
        number = self.station['number']
        path1 = os.path.join(self.path, "GPS_Upar_2016/SCEX_TIPEX3_UPAR_GPS_MUL_")
        path = path1 + str(number) + "-2016"+self.month_num
        aa = os.listdir(path)  # 文件名列表
        aa.sort()  # 排一下顺序，这个是对列表本身进行操作
        ds_time = []  # 每个时次变量的列表
        ttt = []  # 时间序列列表
        for flnm in aa:
            fl_time = flnm[-14:-4]
            tt = pd.to_datetime(fl_time, format='%Y%m%d%H')
            ttt.append(tt)
            # 这时间是不规则的
            flnm = os.path.join(path, flnm)
            aa = self.read_single(flnm)
            ds_time.append(aa)  # 很多时次都是到595hPa才有值, 气压和高度的对应关系会随着时间发展而变化, 气压坐标和高度坐标不能通用
        da = xr.concat(ds_time, dim='time')
        da.coords['time'] = ttt
        # print(da)
        return da


def third_scientific_main():
    pass
    ## 只有三个站
    key_list = ['GaiZe', 'ShenZha', 'ShiQuanhe']
    # month = 'May'
    month = 'Jul'
    ds = xr.Dataset()
    for key in key_list:
        station = station_dic[key]
        # print(station)
        gb = GetThirdScientific(station=station, month=month)
        da = gb.read_obs()
        ## 将pressure这一维度，放到最后一个，方便后面的矩阵计算
        da = da.transpose(*(...,'pressure'))
        da_diagnostic = gb.caculate_diagnostic(da)
        ## 将原来的变量和计算的诊断变量合并为一个DataArray
        dda = xr.concat([da, da_diagnostic], dim='variable')
        ## 不同站点的数据组合为一个Dataset
        ds[key] = dda

    ## 保存文件
    flnm ='/mnt/zfm_18T/fengxiang/DATA/UPAR/GPS_Upar_2016_'+str(month)+".nc" 
    ds.to_netcdf(flnm)


if __name__ == '__main__':

    # %%
    key_list = ['GaiZe', 'ShenZha', 'ShiQuanhe']
    # month = 'May'
    month = 'Jul'
    ds = xr.Dataset()
    # %%
    for key in key_list:
        station = station_dic[key]
        # print(station)
        gb = GetThirdScientific(station=station, month=month)
        da = gb.read_obs()
        ## 将pressure这一维度，放到最后一个，方便后面的矩阵计算
        da = da.transpose(*(...,'pressure'))
        ds = da.to_dataset(dim='variable')
        ds_diagnostic = gb.caculate_diagnostic(ds)
        ds_return = xr.merge([ds, ds_diagnostic])
        dda = ds_return.to_array()
        ds[key] = dda
    flnm  = '/mnt/zfm_18T/fengxiang/DATA/UPAR/GPS_Upar_2016_'+str(month)+".nc"
    # ds.to_netcdf(flnm)

    # --------------------------------------
    ### 测试        
    # --------------------------------------
    # # %%
    # key_list = ['GaiZe', 'ShenZha', 'ShiQuanhe']
    # # month = 'May'
    # month = 'Jul'
    # ds = xr.Dataset()

    # # for key in key_list:
    # station = station_dic['GaiZe']
    # # print(station)
    # gb = GetThirdScientific(station=station, month=month)
    # da = gb.read_obs()
    # da = da.transpose(*(...,'pressure'))
    # ds = da.to_dataset(dim='variable')
    # ds_diagnostic = gb.caculate_diagnostic(ds)
    # ds_return = xr.merge([ds, ds_diagnostic])
    # dda = ds_return.to_array()

    