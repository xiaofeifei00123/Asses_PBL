#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取fnl数据，
将相对湿度转为比湿保存
5月和7月所有湿度数据保存为一个文件
这个fnl资料里面有
gh, t, r, u, v

同时将其插值到站点上
-----------------------------------------
Author           :Forxd
Version          :1.0
Time：2021/07/12/ 15:28
'''
import os

import xarray as xr
from metpy.units import units
from metpy.calc import (dewpoint_from_relative_humidity,
                        specific_humidity_from_dewpoint,
                        mixing_ratio_from_specific_humidity,
                        virtual_potential_temperature,
                        )

from data_process_main import GetData
from global_variable import station_dic
import numpy as np

class GetFnl():
    
    def __init__(self,) -> None:
        pass
        # self.path = '/mnt/zfm_18T/fengxiang/DATA/FNL/FNL_2016/*20160701*00.grib2'  

        self.path = '/mnt/zfm_18T/fengxiang/DATA/FNL/FNL_2016/*20160[5,7]*00.grib2'  

        self.fnl_file = "/mnt/zfm_18T/fengxiang/DATA/FNL/fnl_2016.nc"
        self.pressure_level = np.arange(570, 280, -5)

    def concat_fnl(self,):
        """获取相对湿度等fnl变量
        将它们聚合成一个ds文件
        """
        ## 这个支持正则表达式
        # path = '/mnt/zfm_18T/fengxiang/DATA/FNL/FNL_2016/*20160701*00.grib2'  
        # path = '/mnt/zfm_18T/fengxiang/DATA/FNL/FNL_2016/*201605*00.grib2'  
        # path = '/mnt/zfm_18T/fengxiang/DATA/FNL/FNL_2016/*20160[5,7]*00.grib2'  
        path = self.path
        fl_list = os.popen('ls {}'.format(path))  # 打开一个管道
        fl_list = fl_list.read().split()
        # print(fl_list)

        rh_list = []
        for fl in fl_list:
            ds = xr.open_dataset(fl, engine='cfgrib',
                                backend_kwargs={'filter_by_keys':
                                    {'typeOfLevel': 'isobaricInhPa'}})
            # ds = ds.rename({'r':'rh', 't':'temp'})  # 统一变量名称
            # rh = ds[var]
            rh = ds['r']
            t = ds['t']
            u = ds['u']
            v = ds['v']
            dds = xr.Dataset()
            dds['rh'] = rh
            dds['temp'] = t
            dds['u'] = u
            dds['v'] = v
            
            # rh = ds
            rh_list.append(dds)
        ds = xr.concat(rh_list, dim='time')

        ds = ds.rename({'isobaricInhPa': 'pressure',
                        'latitude': 'lat', 'longitude': 'lon'})
        ds.attrs['description'] = 'the combine of all time rh, full grid'
        # print(ds)
        return ds

    def get_station(self):
        """插值到指定气压高度
           指定站点位置
        """
        pass
        ds_station = xr.Dataset()
        ds_origin = xr.open_dataset(self.fnl_file)  # fnl格点文件
        for key in station_dic:
            station = station_dic[key]
            # da = ds.rh
            ## 水平插值
            ds = ds_origin.sel(lat=station['lat'], 
                            lon=station['lon'],
                            method='nearest')
            ## 垂直插值
            dds = ds.interp(pressure=self.pressure_level)
            dda = dds.to_array()
            ds_station[station['name']] = dda
        return ds_station






if __name__ == '__main__':
    pass
    
    # 对fnl数据进行聚合和计算
    gf = GetFnl()
    # ds = gf.concat_fnl()
    # gd = GetData()
    # ds_diag = gd.caculate_diagnostic(ds)
    # ds_return = xr.merge([ds, ds_diag])
    # ds_return['wind_s'] = xr.ufuncs.sqrt(ds['u']**2+ds['v']**2)
    # ds_return.to_netcdf("/mnt/zfm_18T/fengxiang/DATA/FNL/fnl_2016.nc")
    
    # 对fnl数据进行插值
    dd = gf.get_station()
    fnl_station_path = '/mnt/zfm_18T/fengxiang/Asses_PBL/data/fnl_station.nc'
    dd.to_netcdf(fnl_station_path)




    #-----------------------------------
    #### 测试
    #-----------------------------------
    # fnl_file = "/mnt/zfm_18T/fengxiang/DATA/FNL/fnl_2016.nc"
    # pressure_level = np.arange(570, 280, -5)

    # # %%
    # # ds_station = xr.Dataset()
    # # for key in station_dic:
    # station = station_dic['GaiZe']
    # ds = xr.open_dataset(fnl_file)  # fnl格点文件

    # # %%
    # # da = ds.rh
    # ## 水平插值
    # ds = ds.sel(lat=station['lat'], 
    #                 lon=station['lon'],
    #                 method='nearest')
    # ## 垂直插值
    # dds = ds.interp(pressure=pressure_level)

    # # %%
    # dda = dds.to_array()
    

# %%
