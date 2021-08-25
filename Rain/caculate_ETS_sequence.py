#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
计算并绘制逐小时的
画SAL评分的图
二分类评分的图
计算时间序列评分
-----------------------------------------
Time             :2021/06/04 14:32:20
Author           :Forxd
Version          :1.0
'''

# %%
from numpy.lib.shape_base import column_stack
from pandas.core.frame import DataFrame
import xarray as xr
import meteva.method as mem
import meteva.base as meb
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from read_data import TransferData
from read_data import GetData

from global_variable import station_dic

# %%

class Caculate():

    def __init__(self, rain,) -> None:
        # self.flag = flag
        # self.area = area
        self.rain = rain
        # self.threshold = threshold

    def get_two_scale(self, threshold):
        # flag = 'all'
        """计算ETS评分这些值
            这个阈值是要变化的
        """
        rain = self.rain  # 空间分布的总降水
        model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']

        ### 计算二分类评分
        hfmc_dic = {}
        Accuracy = {}
        ETS = {}
        TS = {}
        POFD = {}
        HK = {}
        ODR = {}
        BIAS = {}
        for model in model_list:
            # hfmc = mem.hfmc(rain['obs'].values/30, rain[model].values/30, grade_list=[0.25])

            ## mask高原外的数据
            rain_obs = rain['obs'].values
            rain_model = rain[model].values
            # print(rain_model.max())
            print("计算 %s" %model)
            hfmc = mem.hfmc(rain_obs, rain_model, grade_list=[threshold])   # 这里算的都是平均值的评分, 我想要的是评分的平均值, 也就是要算每个时次的评分
            # ETS[model] = mem.ets_hfmc(hfmc) 
            TS[model] = mem.ts_hfmc(hfmc) 
            Accuracy[model] = mem.pc_hfmc(hfmc)
            # POFD[model] = mem.pofd_hfmc(hfmc)
            # POFD[model] = mem.far_hfmc(hfmc)
            # HK[model] = mem.hk_yesorno_hfmc(hfmc)
            # ODR[model] = mem.odds_ratio_hfmc(hfmc)
            BIAS[model] = mem.bias_hfmc(hfmc)

        # grade_list = [Accuracy, ETS, TS, POFD, HK, ODR, BIAS]
        # grade_list_name = ['Accuracy', 'ETS', 'TS', 'POFD', 'HK', 'ODR', 'BIAS']
        grade_list = [Accuracy, TS,  BIAS]
        grade_list_name = ['Accuracy', 'TS', 'BIAS']

        b = []
        for i in grade_list:
            a = pd.Series(i)
            b.append(a)
            # print(a)
        df = pd.concat(b, axis=1)
        df.columns = grade_list_name
        df = df.T
        # print(df)
        # df.to_csv('/home/fengxiang/Project/Asses_PBL/Data/tt.csv')
        return df



    def get_space_scale(self, rain_threshold):
        """计算空间评分
        """
        rain = self.rain
        """空间降水"""
        model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']

        sal_list = []
        for model in model_list:
            pass
            print("计算 %s" %model)
            # tttime = pd.date_range()
            rain_model = rain[model]
            rain_obs  = rain['obs']
            grid_obs = meb.xarray_to_griddata(rain_obs)
            grid_model = meb.xarray_to_griddata(rain_model)

            look_ff = mem.mode.feature_finder(grid_obs, grid_model, smooth=10, threshold=rain_threshold/15, minsize=5)
            look_match = mem.mode.centmatch(look_ff)
            look_merge = mem.mode.merge_force(look_match)
            sal = mem.sal(look_ff)
            ssal = pd.Series(sal)
            sal_list.append(ssal)
        
        df = pd.concat(sal_list, axis=1)
        df.columns = model_list
        return df

    def get_time_scale(self, rain):
        """计算时间序列的评分
        rain是时间序列的降水, 就是一个一维的数组
        """
        # rd = rd(self.flag, self.area)
        # rain = rd.get_rain_total()
        # rain = rd.get_rain_times()
        # print(rain)
        model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']

        MAE = {} # 平均绝对误差
        RMSE = {} # 均方根误差
        SD = {}  # 预报标准差/观测标准差
        CORR = {}
        for model in model_list:
            pass
            rain_model = rain[model].values
            rain_obs  = rain['obs'].values
            # MAE[model] = mem.mae(rain_obs, rain_model)
            sd = mem.ob_fo_std(rain_obs, rain_model)
            SD[model] = sd[1]/sd[0]
            RMSE[model] = mem.rmse(rain_obs, rain_model)
            CORR[model] = mem.corr(rain_obs, rain_model)
        # grade_list = [MAE, RMSE, SD, CORR]
        # grade_list_name = ['MAE', 'RMSE', 'SD', 'CORR']
        # grade_list = [RMSE, SD, CORR]
        grade_list = [RMSE, SD]
        # grade_list_name = ['RMSE', 'SD[Model]/SD[OBS]']
        grade_list_name = ['RMSE', 'SD[Fcst]/SD[Obs]']
        # grade_list_name = ['RMSE', 'SD', 'CORR']
        # grade_list_name = ['Accuracy', 'TS', 'POFD', 'HK', 'ODR', 'BIAS']

        b = []
        for i in grade_list:
            a = pd.Series(i)
            b.append(a)
            # print(a)
        df = pd.concat(b, axis=1)
        df.columns = grade_list_name
        df = df.T
        print(df)
        # df.to_csv('/home/fengxiang/Project/Asses_PBL/Data/time_score.csv')
        return df
        
# %%
class Pretreat():
    """
    继承读取数据的类, 直接用它的方法
    获取原始数据
    预处理成需要的数据
    """
    def __init__(self, month, time_flag) -> None:
        pass
        self.month = month
        self.time_flag = time_flag
        self.area = {"lat1": 24.875, "lat2": 40.125, "lon1": 69.875, "lon2": 105.125}

    def get_rain_sum(self, gd):
        """计算月总降水
        """
        rain = gd.get_rain_hourly()
        tr = TransferData(ds=rain, area=self.area, time_flag=self.time_flag)
        rain = tr.rain_time_total()
        return rain

    def caculate_threshold(self, rain_sum):
        """"计算降水阈值"""
        rain = rain_sum['obs'].values  # 区域内观测降水值
        rain_dim1 = rain.flatten()  # 二维降水数据变为1维
        # tt = rain_dim1[~np.isin(rain_dim1, 0.)]  # 删除0降水的点
        # rain_dim1 = np.sort(tt)  # 排序
        rain_dim1 = np.sort(rain_dim1)  # 这里有没有计算0降水量影响在2mm内
        max_num = len(rain_dim1) - 1  # 取最大降水序号
        # 最大降水序号*0.95为降水阈值
        rain_threshold = rain_dim1[int(max_num * 0.95)]
        # rain_dim1
        # rain_threshold/15
        return rain_threshold


# %%
def get_grade_time():
    """获得日降水多个站点的评分
    """
    pass
    month = 'Jul'
    gd = GetData(month)
    rain = gd.get_rain_hourly()  # 直接获取小时降水数据
    rain = rain.resample(time='D').sum()   # 日降水时间序列
    # pr = Pretreat(month, 'all') # 只能是all, 这个不区分白天和夜间
    # rain_day = pr.get_rain_day(gd)
    rain_station_list = []
    for key in station_dic:  
        station = station_dic[key]
        # rain1 = tr.rain_station(station_dic_dic[key])  
        r = rain.sel(lat=station['lat'], lon=station['lon'], method='nearest')
        rain_station_list.append(r)  # 站点降水的列表
    rain_station_array = xr.concat(rain_station_list, pd.Index(station_dic.keys(), name='station'))
    # rain_24h_mean = rain_station_array.mean(dim='station')  # 多站的平均值
    return rain_station_array

def concat_grade_time():
    rain_station_array = get_grade_time()

    ## 计算每个站点的时间评分
    time_score_list = []
    for key in station_dic:
        print("计算 %s 站的时间评分" %key)
        rain_station = rain_station_array.sel(station=key)
        ca = Caculate(rain_station)
        rain_station_score  = ca.get_time_scale(rain_station)
        time_score_list.append(rain_station_score)
    df = pd.concat(time_score_list)

    # 增加计算平均降水的评分
    rain_24h_mean = rain_station_array.mean(dim='station')  # 多站的平均值
    ca = Caculate(rain_24h_mean)
    rain_station_mean_score  = ca.get_time_scale(rain_24h_mean)
    df_return = pd.concat([df, rain_station_mean_score])
    df_return = df_return.round(2)

    ## 给各评分增加站点名称
    aa = station_dic.keys()
    bb = []
    for i in aa:
        bb.append(i)
        bb.append(None)
        # bb.append(i)
    bb.append('Mean')
    bb.append(None)

    # bb.append('mean')
    df_return['Station'] = bb
    # df = df_return

    # 将station 变为第一列    
    df_return = df_return.reset_index()
    new_columns = df_return.columns.to_list()
    new_columns.insert(0,'Station')
    new_columns.pop(-1)
    df_return = df_return.reindex(columns=new_columns)
    df_return = df_return.rename(columns={'index':'Grade'})
    
    ## 保存
    df_return.to_excel('/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/time_score.xlsx')
    return df_return

# %%
# dd    
# new_columns
# dc = dd.set_index('station', inplace=False)
# dc = dd.insert(0, 'QNSE')
# new_columns.insert(0,'station')
# new_columns.pop(-1)
# new_columns
# dd
# dc = ddc.reindex(columns=new_columns)
# ddd = dc.rename(columns={'index':'Grade'})
# ddd.to_excel('/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/test.xlsx')
# ddd
# dc.reset_index(level='station')
# ddc


# %%


def single_space_grade(time_flag):
    """获取单个时间段的空间评分

    Args:
        time_flag ([type]): [description]

    Returns:
        [type]: [description]
    """
    pass
    # 获取数据
    month = 'Jul'
    # time_flag = 'all'
    gd = GetData(month)
    pr = Pretreat(month, time_flag)
    rain_sum = pr.get_rain_sum(gd)
    rain_threshold = pr.caculate_threshold(rain_sum)

    ## 计算评分
    ca = Caculate(rain_sum)
    threshold1 = 30
    threshold2 = rain_threshold
    ts1 = ca.get_two_scale(threshold1)
    ts2 = ca.get_two_scale(threshold2)
    sal = ca.get_space_scale(rain_threshold)
    df = pd.concat([ts1, ts2, sal])
    # print(df)
    return df  # 一个时间段的评分


def get_grade_space():
    """获取全部时间段的空间评分
    """
    time_flag_list = ['all', 'day', 'night']

    ## 循环计算评分
    df_list = []
    for time_flag in time_flag_list:
        pass
        df = single_space_grade(time_flag)
        df_list.append(df.round(2))
    ## 存储数据
    writer = pd.ExcelWriter('/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/space.xlsx')
    num = len(time_flag_list)    
    for i in range(num):  # 正好是3个时间段的数据
        df_list[i].to_excel(writer, sheet_name = time_flag_list[i])
    writer.save()
    writer.close()

if __name__ == '__main__':
    concat_grade_time()
    # get_grade_space()
        
