#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
用来读取需要的数据的模块
改用站点插值和cmorph融合的降水数据
-----------------------------------------
Time             :2021/06/04 14:47:01
Author           :Forxd
Version          :1.0
'''
import xarray as xr
import pandas as pd
# import salem  # 过滤高原外的数据, 滤出地形外的数据
# import geopandas
import xesmf as xe  # 插值
import numpy as np
import os
## 读grd文件的库
from xgrads import CtlDescriptor
from xgrads import open_CtlDataset

from global_variable import station_dic


class GetData():

    def __init__(self,):
        # self.flag = flag
        # self.area = area
        self.var = 'RAINNC'
        self.path_wrf = '/mnt/zfm_18T/fengxiang/Asses_PBL/data/wrfout_data'
        self.path_obs = '/mnt/zfm_18T/fengxiang/DATA/PRECIPTATION/CMORPH_STATION_RAIN'

    def get_rain_wrf(self, ):
        """循环出不同试验的结果, 返回空间分布的累计降水
            所有时次降水累积和的空间分布 
        """
        model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','ACM2','YSU']
        file_list = []
        for i in model_list:
            flnm = str(self.var)+"_May_"+i+'_latlon'
            flnm = os.path.join(self.path_wrf, flnm)
            file_list.append(flnm)
        rain = xr.Dataset()
        for i in range(len(file_list)):
            fl = file_list[i]
            ## 数据已经处理过了，ds就是逐小时的降水
            ds = xr.open_dataset(fl)  # 一定要防止变量串台啊
            rain[model_list[i]] = ds.RAINNC  # 构建ds
        return rain

    def get_rain_obs(self,):
        """格点降水和CMORPH降水融合数据
        """
        flnm = '/mnt/zfm_18T/fengxiang/DATA/PRECIPTATION/CMORPH_STATION_RAIN/05/CHN_PRCP_HOUR_MERG_DISPLAY_0.1deg.lnx.ctl'
        ds = open_CtlDataset(flnm)
        da = ds.crain.squeeze(drop=True)
        da = da.where(da.values>=0, np.nan)

        ds_in = da.to_dataset()
        ds_in = ds_in.drop_dims('time')

        ## 插值器构建，前面必须是一个二维的dataset格式的数组
        # ds_out = xe.util.grid_2d(69.75, 105.125, 0.25, 24.75, 45.125, 0.25)
        # ds_out = xe.util.grid_2d(69.75, 105.125, 0.25, 24.75, 45.125, 0.25)
        ds_out = xr.Dataset(
            {'lat':(['lat'],np.arange(24.875, 45.125+0.25, 0.25)),
             'lon':(['lon'], np.arange(69.875, 105.125+0.25, 0.25))}
            )
        regridder = xe.Regridder(ds_in, ds_out, 'bilinear')  

        ### regrid, 插值
        da = regridder(da)  
        
        # da = regridder(da)  
        # xtime = da.time.values
        # lat = np.arange(24.875, 45.125+0.25,0.25)
        # lon = np.arange(69.875, 105.125+0.25,0.25)

        # da_return = xr.DataArray(da.values, coords=[xtime, lat, lon], dims=['time', 'lat', 'lon'])
        # print(da)
        # print(da_return)
        # return da_return
        return da

    def get_rain_hourly(self):
        
        rain = self.get_rain_wrf()  # 这是个ds
        aa = self.get_rain_obs()
        # rain['obs'] = self.get_rain_obs()  # 这个返回da
        rain['obs'] = aa
        return rain

class TransferData():

    def __init__(self, ds, area, time_flag):
        pass
        self.ds = ds
        self.area = area
        self.flag = time_flag  # 时间标签

    def get_time_list(self):
        time_list = self.ds.time.values
        # time_list = r.time.values
        day_list = []
        night_list = []
        for i in time_list:
            hour = str(i)[11:13]
            if int(hour) <= 13:  # 00时--13时时白天, 14--23时是夜间
                day_list.append(i)
            else:
                night_list.append(i)
        print(night_list)
        time_index = {'all':time_list, 'day':day_list, 'night':night_list}
        return time_index

    def rain_hourly(self,):
        """取固定空间范围内，逐小时的数据
        Returns:
            [type]: [description]
        """
        area = self.area
        lat_index = np.arange(area['lat1'], area['lat2']+0.1, 0.25)
        lon_index = np.arange(area['lon1'], area['lon2']+0.1, 0.25)
        # ds_return = rain.sel(lat=lat_index, lon=lon_index, method='nearest')
        ds_return = self.ds.sel(lat=lat_index, lon=lon_index, method='nearest')
        return ds_return

    def rain_time_total(self, ):
        """求所有时次的降水和
            白天、夜间
        """
        ds = self.rain_hourly()  # 筛选过区域的数据
        time_list = ds.time.values
        day_list = []
        night_list = []
        for i in time_list:
            hour = str(i)[11:13]
            if int(hour) <= 13:  # 00时--13时时白天, 14--23时是夜间
                day_list.append(i)
            else:
                night_list.append(i)
        
        ## 根据白天还是夜间的需求，求总降水量
        if self.flag == 'all':
            rrr = ds.sum(dim='time')
        elif self.flag == 'day':
            rain_day = ds.sel(time=day_list)
            rrr = rain_day.sum(dim='time')
        else:
            rain_night = ds.sel(time=night_list)
            rrr = rain_night.sum(dim='time')
        return rrr

    def rain_station(self, station):
        """返回完整区域数据，需要进行数据筛选，自己根据dataset进行数据mask
        """
        # ds = xr.open_dataset(flnm)
        r = self.rain_hourly()
        # shp_file = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_tp/Tibet.shp'
        # shp = geopandas.read_file(shp_file)
        r = r.sel(lat=station['lat'], lon=station['lon'], method='nearest')
        # r = r.mean(dim='lat')
        # r = r.mean(dim='lon')
        # r = r.salem.roi(shape=shp)  # mask高原外的值
        return r

            

if __name__ == '__main__':

    # %%
    pass
    # %%
    ### 获取数据
    time_flag = 'all' ## 白天还是夜间
    # time_flag = 'night' ## 白天还是夜间
    area = {"lat1":24.875, "lat2":45.125, "lon1":70.875, "lon2":105.125}
    gd = GetData()
    rain = gd.get_rain_hourly()
    # aa = gd.get_rain_obs()
    # print(aa.values)
    # print(rain)

    # %%
    ### 筛选数据
    tr = TransferData(ds=rain, area=area, time_flag=time_flag)
    # a = tr.rain_hourly()
    # b = tr.rain_time_total()
    station = station_dic['GaiZe']
    c = tr.rain_station(station)
    print(c)
    # c = tr.rain_space_average()
    

# %%
