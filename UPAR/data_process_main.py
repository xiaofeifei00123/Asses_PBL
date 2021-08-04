#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
资料分类：
    fnl资料
        - 计算好了rh, q, u, v, t, td 都计算好了
    micaps资料
        - 原始的t, td, u, v, 整合为一个文件
    加密探空资料
        - 暂时处理不到，原始数据未做处理
    wrfout资料
        - 各要素单独存放为一个文件

目标数据格式：
    某站点(筛选出的经纬度)
    多时次
    多高度
    多要素(t, td, rh, u, v, wind_s, theta_v, theta)
    ds


获取探空观测图所需要的数据
需要哪些资料(站点的，随高度和时间变化的):
    1. T -- 温度
    2. Td -- 露点温度
    3. T-Td
    4. wind_s -- 风速

    5. Theta_v -- 虚位温
    6. Theta -- 位温
    7. q -- 比湿

    以及它们随高度变化的梯度

    获取不同月份的数据，5月和7月

注意：
    垂直坐标的选取(h, p)
    缺省数据的处理
    所有时次的数据都计算好，不进行时间平均等操作

数据保存格式:
    某一观测站，某一变量，
    所有模式和观测数据统一保存成一个xr.dataset的格式

要求, 'May'：
    处理好缺省值
    插值好
-----------------------------------------
Time             :Mon Jun 28 20:02:27 CST 2021
Author           :Forxd
Version          :0.3
'''

import xarray as xr
import os
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from wrf import getvar, vinterp, interplevel

from metpy.units import units
from metpy.calc import specific_humidity_from_dewpoint
from metpy.calc import mixing_ratio_from_specific_humidity
from metpy.calc import virtual_potential_temperature
from metpy.calc import potential_temperature
from metpy.calc import relative_humidity_from_dewpoint
from metpy.calc import dewpoint_from_relative_humidity
from xarray.core import variable

from global_variable import station_dic


class GetData():
    """获取数据的公共变量
    """

    def grads_data(self, model_dic):
        """ 子程序
        根据输入的垂直文件，求它做梯度后的文件
        Args:
            model_dic : {'model':da[pressure, time]}
        """
        pass
        dic_return = {}
        for key in model_dic:
            ser = []
            da = model_dic[key]
            pressure = da.pressure.values
            height = da.height_coord.values
            del da['height_coord']
            # 第0层是600hPa
            for i in range(len(pressure)):
                if i == 0:
                    # 边界上使用前差或后差
                    sr = (da[:, i + 1] - da[:, i]) / \
                        (pressure[i + 1] - pressure[i])
                    ser.append(sr)
                # 中间的使用中央差分
                elif i > 0 and i < len(pressure) - 1:
                    sr = (da[:, i + 1] - da[:, i - 1]) / \
                        (pressure[i] - pressure[i - 1]) / 2
                    ser.append(sr)
                elif i == len(pressure) - 1:
                    pass
                    sr = (da[:, i] - da[:, i - 1]) / \
                        (pressure[i] - pressure[i - 1])
                    ser.append(sr)
            da = xr.concat(ser, 'pressure')
            da.coords['pressure'] = pressure
            da = da.transpose()
            dic_return[key] = da
        return dic_return


    def caculate_diagnostic(self, ds):
        """计算比湿，位温等诊断变量
        根据td或者rh计算q,theta_v
        返回比湿q, 虚位温theta_v, 相对湿度rh

        Args:
            ds (Dataset): 包含有temp ,td的多维数据
            这里传入Dataset合理一点
        """
        pass        
        ## 获得温度和露点温度
        dims_origin = ds['temp'].dims  # 这是一个tuple, 初始维度顺序
        ds = ds.transpose(*(...,'pressure'))

        var_list = ds.data_vars
        t = ds['temp']

        ## 转换单位
        pressure = units.Quantity(t.pressure.values, "hPa")

        ## 针对给的rh 或是td做不同的计算
        ## 需要确定t的单位
        if 'td' in var_list:
            """探空资料的温度单位大多是degC"""
            td = ds['td']
            dew_point = units.Quantity(td.values, "degC")
            temperature = units.Quantity(t.values, "degC")
        elif 'rh' in var_list:
            """FNL资料的单位是K"""
            # rh = da.sel(variable='rh')
            rh = ds['rh']
            rh = units.Quantity(rh.values, "%")
            temperature = units.Quantity(t.values, "K")
            dew_point = dewpoint_from_relative_humidity(temperature, rh)
        else:
            print("输入的DataArray中必须要有rh或者td中的一个")
        
        ## 记录维度坐标
        # time_coord = t.time.values
        # pressure_coord = t.pressure.values

        ## 计算诊断变量
        q = specific_humidity_from_dewpoint(pressure, dew_point)
        w = mixing_ratio_from_specific_humidity(q)
        theta_v = virtual_potential_temperature(pressure, temperature, w)

        if 'td' in var_list:
            rh = relative_humidity_from_dewpoint(temperature, dew_point)
            var_name_list = ['q', 'rh', 'theta_v']
            var_data_list = [q, rh, theta_v]
        elif 'rh' in var_list:
            pass
            var_name_list = ['q', 'td', 'theta_v']
            var_data_list = [q, dew_point, theta_v]

        ## 融合各物理量为一个DataArray

        ds_return = xr.Dataset()

        for var_name, var_data in zip(var_name_list, var_data_list):
            pass
            ## 为了去除单位
            dda = xr.DataArray(
                var_data, 
                # coords=[time_coord,pressure_coord],
                coords = t.coords,
                dims=t.dims)
                # dims=['time', 'pressure'])

            ds_return[var_name] = xr.DataArray(
                dda.values, 
                # coords=[time_coord,pressure_coord],
                coords=t.coords,
                dims=t.dims)
        # da_return  = ds_return.to_array()

        ## 转换维度顺序        
        # da_return = da_return.transpose(*dims_origin)
        ds_return = ds_return.transpose(*dims_origin)
        return ds_return

class SaveData():
    pass

if __name__ == '__main__':
    pass



    # %%
    # station = station_dic['TingRi']
#     station = station_dic['ShiQuanhe']
#     tr = TransferData(station, 'May')
#     # %%
#     # aa = tr.get_data_one('temp')
#     # print(aa)
#     model_dic = tr.transfer_data('theta_v')  # 多试验的某变量(temp, rh..)数据
#     print(model_dic)  
#     # TODO 别传model_dic了，全部换成ds
#     # %%
#     # aa = model_dic['obs']
#     aa = model_dic['fnl']
#     for i in range(len(aa.time)):
#         print(aa.isel(time=i).values)
#     # print(aa)
# # %%
