#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取fnl数据，
将相对湿度转为比湿保存
5月和7月所有湿度数据保存为一个文件
这个fnl资料里面有
gh, t, r, u, v
-----------------------------------------
Author           :Forxd
Version          :1.0
Time：2021/07/12/ 15:28
'''
import os
from metpy.calc.thermo import dry_static_energy 
import xarray as xr
from metpy.units import units
from metpy.calc import dewpoint_from_relative_humidity
from metpy.calc import specific_humidity_from_dewpoint

class GetFnl():
    
    def __init__(self) -> None:
        pass

    def concat_fnl(self,):
        """获取相对湿度等fnl变量
        将它们聚合成一个ds文件
        """
        ## 这个支持正则表达式
        path = '/mnt/zfm_18T/fengxiang/DATA/FNL/FNL_2016/*201607*00.grib2'  
        # path = '/mnt/zfm_18T/fengxiang/DATA/FNL/FNL_2016/*201605*00.grib2'  
        # path = '/mnt/zfm_18T/fengxiang/DATA/FNL/FNL_2016/*20160[5,7]*00.grib2'  
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

    def caculate_fnl(self, ds):
        """根据已有的ds文件内的数据，计算出q等其他变量

        Args:
            ds ([type]): 合并时间后的fnl文件
        """
        t = ds['temp'].squeeze(drop=True)
        rh = ds['rh'].squeeze(drop=True)
        pressure = units.Quantity(t.pressure.values, "hPa")
        temperature = units.Quantity(t.values, "degC")
        td = dewpoint_from_relative_humidity(t, rh)

        time_coord = t.time.values
        pressure_coord = t.pressure.values
        lat_coord = t.lat.values
        lon_coord = t.lon.values
        ## 转换维度顺序
        dew_point = td.transpose(*('lon', 'lat', 'time', 'pressure'))
        q = specific_humidity_from_dewpoint(pressure, dew_point)
        q = q.transpose(*('time', 'pressure', 'lat', 'lon'))
        ds_return = xr.Dataset()

        ## 因为metpy和DataArray的结构不一致，需要保存的话，必须转换数据格式

        ds_return['q'] = xr.DataArray(
            q, 
            coords=[time_coord,pressure_coord,lat_coord, lon_coord],
            dims=['time', 'pressure', 'lat', 'lon'])

        ds_return['td'] = xr.DataArray(
            td, 
            coords=[time_coord,pressure_coord,lat_coord, lon_coord],
            dims=['time', 'pressure', 'lat', 'lon'])
            
        ds_return['u'] = xr.DataArray(
            ds['u'], 
            coords=[time_coord,pressure_coord,lat_coord, lon_coord],
            dims=['time', 'pressure', 'lat', 'lon'])
            
        ds_return['v'] = xr.DataArray(
            ds['v'], 
            coords=[time_coord,pressure_coord,lat_coord, lon_coord],
            dims=['time', 'pressure', 'lat', 'lon'])

        ds_return['rh'] = xr.DataArray(
            ds['rh'], 
            coords=[time_coord,pressure_coord,lat_coord, lon_coord],
            dims=['time', 'pressure', 'lat', 'lon'])
            
        return ds_return


if __name__ == '__main__':
    gf = GetFnl()
    # %%
    ds = gf.concat_fnl()

    # %%
    ds_return = gf.caculate_fnl(ds)

    
    # %%
    print(ds_return)
    # ds_return.to_netcdf("/mnt/zfm_18T/fengxiang/DATA/FNL/fnl_201605.nc")
    ds_return.to_netcdf("/mnt/zfm_18T/fengxiang/DATA/FNL/fnl_201607.nc")
    # %%
    # def get_fnl(ds):
    #     t = ds['temp'].squeeze(drop=True)
    #     rh = ds['rh'].squeeze(drop=True)
    #     pressure = units.Quantity(t.pressure.values, "hPa")
    #     temperature = units.Quantity(t.values, "degC")
    #     td = dewpoint_from_relative_humidity(t, rh)

    #     ## 转换维度顺序
    #     dew_point = td.transpose(*('lon', 'lat', 'time', 'pressure'))
    #     q = specific_humidity_from_dewpoint(pressure, dew_point)
    #     q = q.transpose(*('time', 'pressure', 'lat', 'lon'))
    #     ds_return = xr.Dataset()
    #     ds_return['q'] = q
    #     ds_return['td'] = td
    #     ds_return['temp'] = t
    #     ds_return['u'] = ds['u']
    #     ds_return['v'] = ds['v']
    #     return ds_return
    # aa = get_fnl(ds)

# %%
