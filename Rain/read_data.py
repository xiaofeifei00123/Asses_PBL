#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
用来读取需要的数据的模块
-----------------------------------------
Time             :2021/06/04 14:47:01
Author           :Forxd
Version          :1.0
'''
import xarray as xr
import pandas as pd

import salem  # 过滤高原外的数据
import geopandas
import numpy as np

class Get_data(object):
    # class Summer_mean(object):
        
    def __init__(self, flag,area):
        self.flag = flag
        self.area = area
        pass
    def month_sum_whole(self,ds,rain_type):
        """求月降水之和,
            白天降水之和, 
            夜间降水之和
         """
        r = ds[self.var]
        shp_file = '/home/fengxiang/Data/shp_tp/Tibet.shp'
        shp = geopandas.read_file(shp_file)
        r = r.salem.roi(shape=shp)  # mask高原外的值
        area = self.area
        if rain_type == 'obs':
            pass
            r = r.loc[:,area['lat2']:area['lat1'],area['lon1']:area['lon2']]  # 一天8次的降水
        else:
            r = r.loc[:,area['lat1']:area['lat2'],area['lon1']:area['lon2']]  # 一天8次的降水
        ## 判断某时次是白天还是夜间, 获得白天和夜间时次列表
        time_list = r.time.values
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
            rrr = r.sum(dim='time')
        elif self.flag == 'day':
            rain_day = r.sel(time=day_list)
            rrr = rain_day.sum(dim='time')
        else:
            rain_night = r.sel(time=night_list)
            rrr = rain_night.sum(dim='time')
        return rrr

    def get_rain_spacea_average(self,ds, rain_type):
        pass
        # ds = xr.open_dataset(flnm)
        r = ds[self.var]
        shp_file = '/home/fengxiang/Data/shp_tp/Tibet.shp'
        shp = geopandas.read_file(shp_file)
        r = r.salem.roi(shape=shp)  # mask高原外的值
        # area = {"lat1":24.875, "lat2":45.125, "lon1":69.875, "lon2":105.125, "rain_type":'obs'}
        # r = r.loc[:,45.125:24.875,69.875:105.125]  # 一天8次的降水
        area = self.area
        if rain_type == 'obs':
            pass
            r = r.loc[:,area['lat2']:area['lat1'],area['lon1']:area['lon2']]  # 一天8次的降水
        else:
            r = r.loc[:,area['lat1']:area['lat2'],area['lon1']:area['lon2']]  # 一天8次的降水

        ## 判断某时次是白天还是夜间, 获得白天和夜间时次列表
        time_list = r.time.values
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
            # rrr = r.sum(dim='time')
            rr = r.mean(dim='lat')
            rrr = rr.mean(dim='lon')
        elif self.flag == 'day':
            rain_day = r.sel(time=day_list)
            rr = rain_day.mean(dim='lat')
            rrr = rr.mean(dim='lon')
            # rrr = rain_day.sum(dim='time')
        else:
            rain_night = r.sel(time=night_list)
            # rrr = rain_night.sum(dim='time')
            rr = rain_night.mean(dim='lat')
            rrr = rr.mean(dim='lon')
        return rrr

class Get_wrf(Get_data):
    def __init__(self, flag,area, var,):
        super().__init__(flag, area)
        self.var = var

    def get_rain_hourly(self, ):
        """循环出不同试验的结果, 返回空间分布的累计降水
            所有时次降水累积和的空间分布 
        """

        model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','MYJ','YSU']
        #model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','YSU']
        path = '/home/fengxiang/Project/Asses_PBL/Data/standard_data/'
        file_list = []
        for i in model_list:
            # flnm = path + "RAINNC_Jul_"+i+'_latlon'
            flnm = path +str(self.var)+"_Jul_"+i+'_latlon'
            file_list.append(flnm)

        rain = {}   # 各模式降水和, 各模式的降水字典
        for i in range(len(file_list)):
            fl = file_list[i]
            ## 数据已经处理过了，ds就是逐小时的降水
            ds = xr.open_dataset(fl)  # 一定要防止变量串台啊
            # rain_total = self.month_sum_whole(ds, 'model')
            rain[model_list[i]] = ds
        return rain
        
    def get_rain_wrf(self, ):
        """循环出不同试验的结果, 返回空间分布的累计降水
            所有时次降水累积和的空间分布 
        """

        model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','MYJ','YSU']
        #model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','YSU']
        path = '/home/fengxiang/Project/Asses_PBL/Data/standard_data/'
        file_list = []
        for i in model_list:
            # flnm = path + "RAINNC_Jul_"+i+'_latlon'
            flnm = path +str(self.var)+"_Jul_"+i+'_latlon'
            file_list.append(flnm)

        rain = {}   # 各模式降水和, 各模式的降水字典
        for i in range(len(file_list)):
            fl = file_list[i]
            ## 数据已经处理过了，ds就是逐小时的降水
            ds = xr.open_dataset(fl)  # 一定要防止变量串台啊
            rain_total = self.month_sum_whole(ds, 'model')
            rain[model_list[i]] = rain_total
        return rain

    def get_rain_wrf_times(self,):
        """返回区域平均降水的时间变化序列, 
            区域平均 
        """
        pass
        shp_file = '/home/fengxiang/Data/shp_tp/Tibet.shp'
        shp = geopandas.read_file(shp_file)
        model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','MYJ','YSU']
        #model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','YSU']
        path = '/home/fengxiang/Project/Asses_PBL/Data/standard_data/'
        file_list = []
        for i in model_list:
            # flnm = path + "RAINNC_Jul_"+i+'_latlon'
            flnm = path +str(self.var)+"_Jul_"+i+'_latlon'
            file_list.append(flnm)

        rain = {}   # 各模式降水和, 各模式的降水字典
        for i in range(len(file_list)):
            fl = file_list[i]
            ds = xr.open_dataset(fl)
            rain_times = self.get_rain_spacea_average(ds, 'model')
            rain[model_list[i]] = rain_times
        return rain
    

class Get_obs(Get_data):
    def __init__(self, flag,area, var):
        super().__init__(flag,area)
        self.var = var

    def get_file_list(self,):
        """"获得文件路径+文件名列表"""
        path = '/mnt/Disk4T_5/fengxiang_file/Data/CMORPH/'
        # flnm = 'cmorph.3hr-025deg.20160501.nc'
        ## 获得时间列表
        time_start = '20160701'
        time_end= '20160731'
        time_list = []
        time_array = pd.date_range(start=time_start, end=time_end, freq='D')
        for i in time_array:
            time_list.append(i.strftime('%Y%m%d'))

        ## 获得文件列表
        file_list = []
        for i in time_list:
            flnm = path+'cmorph.3hr-025deg.'+i+'.nc'
            file_list.append(flnm)
        return file_list

    def get_rain_hourly(self,):
        """日降水数据文件,一天8次降水，分别对应00-03时，03-06时...,
        这样几个时间段内的降水"""
        ## 第0时次的降水，应该是00-01时，01-02时，02-03时的观测降水
        ## 以观测时次来看，应该是01,02,03时次的降水
        file_list = self.get_file_list()  # 日数据文件列表
        # rain_oneday_list = []
        # rrt = 0
        rain_list = []
        # time_start = '20160701 0000'
        # time_end= '20160731 2300'
        time_start = '20160701 0100'
        time_end= '20160801 0000'
        time_array = pd.date_range(start=time_start, end=time_end, freq='1H')
        for i in range(len(file_list)):
            ds = xr.open_dataset(file_list[i])  # 打开一天的降水数据
            da = ds.cmorph_precip
            da = da.loc[:,45.125:24.875,69.875:105.125]  # 一天8次的降水
            for j in range(8):
                r_hour = da[j]  # 小时降水量的值
                r_hour = r_hour.reset_coords()
                r_hour = r_hour.drop_vars('time')
                rain_list.append(r_hour)
                rain_list.append(r_hour)
                rain_list.append(r_hour)
        ds = xr.concat(rain_list, dim='time')
        ds.coords['time'] = time_array
        # ds_return = ds.isel(time=np.arange(13,733,1))
        ds_return = ds.isel(time=np.arange(12,732,1))
        return ds_return

    def get_rain_obs(self,):
        """所有时次的降水累积和
        """
        ds = self.get_rain_hourly()
        rain_total = self.month_sum_whole(ds, 'obs')
        return rain_total
    
    def get_rain_obs_times(self, ):
        """返回区域平均降水的时间变化序列"""
        pass
        ds = self.get_rain_hourly()
        rain_times = self.get_rain_spacea_average(ds, 'obs')
        return rain_times
    def get_rain_obs_daily(self,):
        pass
        ds = self.get_rain_hourly()




class Return_data():

    def __init__(self, flag='all', area={"lat1":24.875, "lat2":45.125, "lon1":69.875, "lon2":105.125}) -> None:
        pass
        self.flag = flag   # 白天还是夜间, 还是all
        self.area = area


    def get_rain_total(self,):
        """所有时次的区域降水和
            月总降水量 
        """
        pass
        model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','MYJ','YSU']
        rain = {}
        ## 得到所有试验的降水
        rain_modle_RAINNC = Get_wrf(self.flag, self.area,'RAINNC').get_rain_wrf()
        for i in model_list:
            rain_modle = rain_modle_RAINNC[i]  
            rain[i] = rain_modle/30

        gb = Get_obs(self.flag, self.area, 'cmorph_precip')
        rain_obs = gb.get_rain_obs()/30
        rain['obs'] = rain_obs
        return rain


    def get_rain_specific_hour(self, sp_time):
        """得到特定时次的降水

        Args:
            sp_time [str]: '2016/07/02 12' 精确到小时
        """

        model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','MYJ','YSU']
        rain = {}
        ## 得到所有试验的降水
        rain_modle_RAINNC = Get_wrf(self.flag, self.area,'RAINNC').get_rain_hourly()
        for i in model_list:
            rain_modle = rain_modle_RAINNC[i]['RAINNC']
            select_time = pd.to_datetime(sp_time)
            rain_tm = rain_modle.sel(time=select_time)
            rain[i] = rain_tm.RAINNC

        gb = Get_obs(self.flag, self.area, 'cmorph_precip')
        rain_obs = gb.get_rain_hourly()
        rain['obs'] = rain_obs.sel(time=select_time).cmorph_precip
        ### 返回特定时次的，所有模式和观测的降水
        return rain

    def get_rain_specific_day(self, sp_time):
        """特定日期的降水和
            这一天的降水是从00时--24时之间所有的降水
            即01时的降水+....+到下一日00时的降水
            没有区域限制，因为就是要获得所有区域

        Args:
            sp_time [str]: '2016/07/02', 尽量不要使用1日和31日
        """
        select_day = pd.to_datetime(sp_time)
        aa = select_day.strftime('%d')
        if aa in ['01', '31']:
            print("输入的日期没有完整的数据，请输入其他的日期")
        else:
            pass
            start = select_day.strftime('%Y/%m/%d')+' 01'
            time_series = pd.date_range(start, periods=24, freq='1H')
        
        model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','MYJ','YSU']
        rain = {}
        ## 得到所有试验的降水
        rain_modle_RAINNC = Get_wrf(self.flag, self.area,'RAINNC').get_rain_hourly()
        for i in model_list:
            rain_modle = rain_modle_RAINNC[i]['RAINNC']
            select_time = pd.to_datetime(sp_time)
            rain_hour = rain_modle.sel(time=time_series)
            rain[i] = rain_hour.sum(dim='time')

        gb = Get_obs(self.flag, self.area, 'cmorph_precip')
        rain_obs = gb.get_rain_hourly().cmorph_precip
        rain_obs = rain_obs.sel(time=time_series)
        rain_obs = rain_obs.sum(dim='time')
        rain['obs'] = rain_obs
        # ### 返回特定时次的，所有模式和观测的降水
        # return rain
        # print(rain['obs'])
        return rain
        


    def get_rain_times(self,):
        """不同时次(逐小时降水)的区域平均值
            所有时次降水的时间序列 
        """
        pass
        model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','MYJ','YSU']
        rain = {}
        for i in model_list:
            rain_modle_RAINNC = Get_wrf(self.flag, self.area,'RAINNC').get_rain_wrf_times()
            rain_modle = rain_modle_RAINNC[i]
            rain[i] = rain_modle

        gb = Get_obs(self.flag, self.area,'cmorph_precip')
        rain_obs = gb.get_rain_obs_times()
        rain['obs'] = rain_obs
        return rain

    def get_rain_daily(self,):
        """得到日降水量的区域平均值, 的时间序列

        Args:
            flag (str, optional): [是白天还是夜间]. Defaults to 'all'.
            area (dict, optional): [求某个区域]. Defaults to {"lat1":24.875, "lat2":45.125, "lon1":69.875, "lon2":105.125}.

        Returns:
            dic: 分别为每日的区域平均降水
        """
        ds = self.get_rain_times(flag='all', area=self.area)
        time_index = np.arange(0,720+1,24)  # 不能有小数

        model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','MYJ','YSU', 'obs']

        dic_model = {}
        for j in model_list:  # 循环出每种降水
            rain_list = []
            for i in range(len(time_index)-1):  # 循环出日降水文件
                time_array = np.arange(time_index[i], time_index[i+1])  # 每天的时间列表
                ds_return = ds[j].isel(time=time_array)
                rain_day = ds_return.sum(dim='time')  # 日降水, 区域平均的
                rain_list.append(rain_day)

            ## 聚合，添加时间属性
            dic_model[j] = xr.concat(rain_list, dim='time')
            time_array = pd.date_range(start='20160701', end='20160730', freq='1D')
            dic_model[j].coords['time'] = time_array
        return dic_model

    def get_rain_24hour_sequence(self, ):
        """得到日降水量的区域平均值, 的时间序列, 重构成0-23时
        逐小时的一个时间序列, 观测的降水是3小时平均值
        这里的flag=all, 表示不区分白天还是黑夜，因为这个就是算每天的时间变化的, 
        """

        ds = self.get_rain_times(flag=self.flag, area=self.area)
        model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','MYJ','YSU', 'obs']

        dic_model = {}
        time_start = '20160701 1300'
        time_start_end = '20160702 1200'

        time_array_start= pd.date_range(start=time_start, end=time_start_end, freq='1H')
        time_array_end= '20160731 1200'
        for j in model_list:  # 循环出每种试验降水
            rain_list = []
            for i in time_array_start:
                # timearray = 
                time_array = pd.date_range(start=i, end=time_array_end, freq='1D')
                rain_hour = ds[j].sel(time=time_array).mean(dim='time')
                rain_list.append(rain_hour)

            ## 聚合，添加时间属性
            dic_model[j] = xr.concat(rain_list, dim='time')
            # time_array = time_array_start
            dic_model[j].coords['time'] = time_array_start.hour

            # dic_model[j] = xr.CFTimeIndex.sort_values(return_indexer=True, )
            dic_model[j] = dic_model[j].sortby('time')

        return dic_model

    def get_rain_3hour_sequence(self,):

        dic = self.get_rain_24hour_sequence()
        model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','MYJ','YSU', 'obs']
        dic_return = {}
        for model in model_list:
            rain_3h_list = []
            rain_24h = dic[model]
            for i in range(1, 23,3):
                j = i+1
                k = i+2
                if k == 24:
                    k = 0

                rain_3h = rain_24h[i]+rain_24h[j]+rain_24h[k]
                rain_3h_list.append(rain_3h)
            time_array = pd.date_range(start='20160701 03', end='20160702 00', freq='3H')
            rain_3h_array = xr.concat(rain_3h_list, dim='time')
            rain_3h_array.coords['time'] = time_array.hour
            dic_return[model] = rain_3h_array
            dic_return[model] = dic_return[model].sortby('time')
            # dic_re

        return dic_return




if __name__ == '__main__':
    flag = 'all' ## 白天还是夜间
    # Get_obs = Get_obs(flag, 'cmorph_precip')
    # gd = Get_obs.get_rain_obs()
    area = {"lat1":24.875, "lat2":46.125, "lon1":70.875, "lon2":105.125}
    # gd = Get_obs.get_rain_obs_times()
    # aa = get_rain_total(flag, area)
    # aa = get_rain_times(flag, area)
    # aa = get_rain_daily(area)
    # aa = get_rain_3hour_sequence(area)
    # print(aa)
    rd = Return_data(flag, area)
    aa = rd.get_rain_specific_day('2016/07/02 22')