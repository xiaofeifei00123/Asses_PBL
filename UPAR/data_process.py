#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
-----------------------
放弃height变量的添加
-----------------------
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


class GetData():
    """获取数据的公共变量
    """
    def __init__(self, station, month):
        self.station = station
        # self.pressure_level = np.arange(610, 100, -5)
        self.pressure_level = np.arange(570, 280, -5)
        self.path = '/mnt/zfm_18T/Asses_PBL/'
        self.model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        
        self.month = month

        # self.month = 'May'
        if self.month == 'Jul':
            self.month_num = '07'
            self.time_first = '2016-07-01 13:00'
            self.flnm_height_obs = '/mnt/zfm_18T/Asses_PBL/GPS_Upar_2016/SCEX_TIPEX3_UPAR_GPS_MUL_55228-201607/upar_G_55228_2016070206.txt'
            self.rh_file = '/mnt/zfm_18T/Asses_PBL/FNL/fnl_rh_201607'
        elif self.month == 'May':
            # self.month = 'May'
            self.month_num = '05'
            self.time_first = '2016-05-01 13:00'
            self.flnm_height_obs = '/mnt/zfm_18T/Asses_PBL/GPS_Upar_2016/SCEX_TIPEX3_UPAR_GPS_MUL_55228-201605/upar_G_55228_2016051112.txt'
            self.rh_file = '/mnt/zfm_18T/Asses_PBL/FNL/fnl_rh_201605'
        else:
            print("%s这个月份不在数据集内"%self.month)
            
    

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

    def caculate_diagnostic(self, dic_t, dic_td, var):
        # -------------获取数据-------------------
        """这个只能计算一个时次的
        必须要求每一个点上都没有缺省值
        计算，获取诊断变量
        """

        q = {}  # specific humidity,比湿，水汽质量/气团总质量(g/g,g/kg)
        w = {}  # mixing ratio
        theta_v = {}  # virtual potential temperature
        theta = {}  # potential temperature
        dic_return = {}
        rh = {}

        for model in dic_t:
            ## 删除完全缺测列(时间), 主要是TEMF方案温度缺测
            t = dic_t[model].dropna(dim='time', how='all')
            time_coord = t.time.values
            td = dic_td[model].sel(time=time_coord)  # 保证时间一致

            ## 删除有缺测的行, 露点温度缺测较多
            # td = td.dropna(dim='pressure')
            td = td.dropna(dim='pressure', how='any')
            # td = td.interpolate_na(dim='time', method='linear',fill_value='extrapolate')
            # td = td.interpolate_na(dim='time')
            prc = td.pressure.values
            # height = td.height_coord.values

            # for i in range(len(prc) - 1):
            #     if prc[i] > prc[i + 1]:
            #         break
            # prc = prc[0:i + 1]
            # height = height[0:i + 1]
            t = t.sel(pressure=prc).dropna(dim='pressure')
            td = td.sel(pressure=prc)

            # 添加单位
            pressure = units.Quantity(td.pressure.values, "hPa")
            temperature = units.Quantity(t.values, "degC")
            dew_point = units.Quantity(td.values, "degC")

            ## 计算维度
            time_coord = t.time.values
            pressure_coord = t.pressure.values

            if var == 'theta_v':
                q[model] = specific_humidity_from_dewpoint(
                    pressure, dew_point)  # 比湿
                w[model] = mixing_ratio_from_specific_humidity(q[model])
                theta_v[model] = virtual_potential_temperature(
                    pressure, temperature, w[model])
                da = xr.DataArray(theta_v[model], coords=[
                    time_coord, pressure_coord], dims=['time', 'pressure'])
                # da['height_coord'] = ('pressure', t.height_coord)
            elif var == 'rh':  # 相对湿度
                rh[model] = relative_humidity_from_dewpoint(temperature, dew_point)
                da = xr.DataArray(rh[model], coords=[
                    time_coord, pressure_coord], dims=['time', 'pressure'])
                # dic_return[model] = da
            elif var == 'theta':
                theta[model] = potential_temperature(pressure, temperature)
                da = xr.DataArray(theta[model], coords=[
                    time_coord, pressure_coord], dims=['time', 'pressure'])
                # da['height_coord'] = ('pressure', t.height_coord)
                # dic_return[model] = da
            elif var == 'q':
                q[model] = specific_humidity_from_dewpoint(
                    pressure, dew_point)  # 比湿
                # q[model] = q[model]*10**3  # 改变单位为g/kg
                da = xr.DataArray(q[model]*1000, coords=[
                    time_coord, pressure_coord], dims=['time', 'pressure'])
            else:
                print("没有这个诊断变量, 请检查")
            # da['height_coord'] = ('pressure', t.height_coord)
            dic_return[model] = da

        return dic_return

class GetObs(GetData):
    """获取观测数据, 几种变量统一读取
    """
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
        da_return = xr.concat(da_list, dim=var_list)
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
        ds = xr.concat(ds_time, dim='time')
        ds.coords['time'] = ttt
        return ds


class GetWrfout(GetData):
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
            path = '/mnt/zfm_18T/Asses_PBL/wrfout_data/'
            # 不同的方案和时间，在同一站点，各层气压值相差小于1度
            # 故不作分开考虑
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

        def get_height_lev():
            # 得到各层的pressure值
            # 不同时刻各层气压值，差别可以忽略不计,
            # 后面还要对气压层进行插值, 这里不对它做过高精度要求
            path = '/mnt/zfm_18T/Asses_PBL/wrfout_data/'
            # 不同的方案和时间，在同一站点，各层气压值相差小于1度
            # 故不作分开考虑
            height_name = 'height_agl_'+self.month+'_YSU_latlon'
            # flnm_pressure = os.path.join(path, 'height_agl_Jul_YSU_latlon')
            flnm_pressure = os.path.join(path, height_name)
            ds = xr.open_dataset(flnm_pressure)
            pr = ds.height_agl
            # prb = pr.sel(time='2016-07-01 13:00')
            prb = pr.sel(time=self.time_first)
            lat = self.station['lat']
            lon = self.station['lon']
            # prc = prb.sel(lat=32.13, lon=92.5, method='nearest')
            height = prb.sel(lat=lat, lon=lon, method='nearest')
            prc = get_pressure_lev()
            z = xr.DataArray(height.values, coords=[prc], dims=['pressure'])
            z = z.interp(pressure=self.pressure_level)
            return z.values

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
            # da_return = da
            # z = get_height_lev()
            # da_return['height_coord'] = ('pressure', z)
            return da_return

        return regrid()

    def get_data_single_once(self, var, model):
        """读一个模式一个变量的数据
            模块尽可能的小
        """
        file_name = str(var) + "_" + str(
            self.month) + "_" + str(model) + "_latlon"
        flnm_var = os.path.join(self.path,'wrfout_data', file_name)
        ds_var = xr.open_dataset(flnm_var)
        da_var = ds_var[var]
        da_var = self.regrid_wrfout(da_var)
        return da_var

    def get_data_var(self, var):
        model_dic = {}
        for model in self.model_list:
            if var in ['temp', 'td', 'height_agl']:
                model_dic[model] = self.get_data_single_once(var, model)

            elif var == 't_td':
                t = self.get_data_single_once('temp', model)
                td = self.get_data_single_once('td', model)
                model_dic[model] = t - td

            elif var == 'wind_s':
                U = self.get_data_single_once('U', model)
                V = self.get_data_single_once('V', model)
                # model_dic[model] = t - td
                model_dic[model] = xr.ufuncs.sqrt(U**2 + V**2)
        return model_dic


class TransferData(GetObs, GetWrfout):
    """将获得的数据传递出去
    """

    def get_data_one(self, var):
        """获取单个原始变量的值, 比如temp各试验+观测的值
        """
        ds_obs = self.read_obs()  # 观测数据
        if var in ['temp', 'td', 't_td', 'wind_s', 'height_agl']:
            var_dic = self.get_data_var(var)
            var_dic['obs'] = ds_obs.sel(concat_dim=var)
            # var_dic['fnl'] = self.get_data_fnl(var)
        else:
            print("该变量需要计算, 不能直接获取")
        return var_dic

    def transfer_data(self, var):
        """传递数据、存储数据
        对数据进行统一插值处理
        """
        # ds_obs = self.read_obs()  # 观测数据
        if var in ['temp', 'td', 't_td', 'wind_s', 'height_agl']:
            var_dic = self.get_data_one(var)
        elif var[-5:] == 'grads':  # 所有的变量都可以求梯度
            var = var[0:-6]
            var_dic = self.get_data_one(var)
            var_dic = self.grads_data(var_dic)
        elif var in ['theta_v', 'q', 'theta', 'rh']:
            var_dic_t = self.get_data_one('temp')
            var_dic_td = self.get_data_one('td')
            var_dic = self.caculate_diagnostic(var_dic_t,
                                            var_dic_td, var)
            ## 读所有格点的fnl数据, 再插值
            if var == 'rh':
                # ds = xr.open_dataset('/mnt/zfm_18T/Asses_PBL/FNL/fnl_rh_201607')
                ds = xr.open_dataset(self.rh_file)
                da = ds.rh
                ## 水平插值
                da = da.sel(lat=self.station['lat'], 
                                lon=self.station['lon'],
                                method='nearest')
                ## 垂直插值
                var_dic['fnl'] = da.interp(pressure=self.pressure_level)
        else:
            print("输入错误")
            var_dic = None
        return var_dic


if __name__ == '__main__':

    station_dic = {
        'GaiZe': {
            'lat': 32.3,
            'lon': 84.0,
            'name': 'GaiZe',
            'number': '55248',
            'height': 4400,
        },
        'ShenZha': {
            'lat': 30.9,
            'lon': 88.7,
            'name': 'ShenZha',
            'number': '55472',
            'height': 4672
        },
        'ShiQuanhe': {
            'lat': 32.4,
            'lon': 80.1,
            'name': 'ShiQuanhe',
            'number': '55228',
            'height': 4280
        },
        'LaSa': {
            'lat': 29.66,
            'lon': 91.14,
            'name': 'LaSa',
            'number': '55591',
            'height': 3648.8999,
        },
        'TingRi': {
            'lat': 28.63,
            'lon': 87.08,
            'name': 'TingRi',
            'number': '55664',
            'height': 4302,
        },
        'NaQU': {
            'lat': 31.48,
            'lon': 92.06,
            'name': 'NaQu',
            'number': '55299',
            'height': 4508,
        },
        'LinZhi': {
            'lat': 29.65,
            'lon': 94.36,
            'name': 'LinZhi',
            'number': '56312',
            'height': 2991.8,
        },
        'ChangDu': {
            'lat': 31.15,
            'lon': 97.17,
            'name': 'ChangDU',
            'number': '56137',
            'height': 3315,
        },
    }

    station = station_dic['ShenZha']
    tr = TransferData(station, 'May')
    model_dic = tr.transfer_data('rh')  # 多试验的某变量(temp, rh..)数据
    print(model_dic)