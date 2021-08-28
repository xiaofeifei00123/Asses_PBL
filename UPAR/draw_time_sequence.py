#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取单个站点的所有时次，高度，方案的数据
绘制它的时间高度廓线
填色图
-----------------------------------------
Time             :2021/06/17 18:59:18
Author           :Forxd
Version          :2.0
'''

# %%
# from Rain.draw_distribution import draw_dual
import xarray as xr
import os
import numpy as np
import pandas as pd
import wrf
import netCDF4 as nc
from netCDF4 import Dataset
from wrf import getvar, vinterp, interplevel
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import cmaps
from global_cmap import get_cmap_temp, get_cmap_q
import datetime
# from data_process import TransferData
from data_transfer import TransferData
from global_variable import station_dic

# %%
def get_var(station, var):
    """获取一个站点,某个变量的所有时次的数据

    Args:
        station ([type]): 站点的字典
        var ([type]): 变量名

    Returns:
        [DataArray]: 该站点所有时次的数据
    """
    month = 'Jul'
    time_hour = '12'
    tr = TransferData(station, month, time_hour, 'q')
    bb = tr.get_data_all()
    var_return = bb.sel(variable=var)
    return var_return

# %%
#### 测试
# station = station_dic['GaiZe']
# var = 'theta_v'
# aa = get_var(station, var)
# # str(aa.variable.values)
# # time_index = aa.time.sel(time=datetime.time(int(12)))  
# # cc = aa.sel(time=time_index)
# # dd = cc['QNSE']
# # str(dd.coords['variable'].values)
# # aa.dims
# ds = aa
# dims_origin = ds['QNSE'].dims  # 这是一个tuple, 初始维度顺序
# # dims_origin
# ds2 = ds.transpose(*(...,'pressure'))
# # ds2['ACM2']
# # ds['ACM2']
# ds.pressure.values
#### 测试结束

# %%
def grads_data(ds):
    """ 子程序
    根据输入的垂直文件，求它做梯度后的文件
    Args:
        model_dic : {'model':da[pressure, time]}
    """
    # pass
    # dic_return = {}
    # for key in model_dic:
    #     ser = []
    # da = model_dic[key]
    da = ds
    pressure = da.pressure.values
    # height = da.height_coord.values
    # del da['height_coord']
    # 第0层是600hPa
    ser = []
    for i in range(len(pressure)):
    
        if i == 0:
            sr = (da.isel(pressure=i+1) - da.isel(pressure=i)) / (pressure[i + 1] - pressure[i])
            ser.append(sr)
        # 中间的使用中央差分
        elif i > 0 and i < len(pressure) - 1:
            # sr = (da[:, i + 1] - da[:, i - 1]) / \
            #     (pressure[i] - pressure[i - 1]) / 2
            sr = (da.isel(pressure=i+1) - da.isel(pressure=i-1)) / (pressure[i + 1] - pressure[i-1])
            ser.append(sr)
        elif i == len(pressure) - 1:
            pass
            # sr = (da[:, i] - da[:, i - 1]) / \
            #     (pressure[i] - pressure[i - 1])
            sr = (da.isel(pressure=i) - da.isel(pressure=i-1)) / (pressure[i] - pressure[i-1])
            ser.append(sr)
    da = xr.concat(ser, 'pressure')
    da.coords['pressure'] = pressure
    variable_origin = str(da.variable.values)
    variable_new = variable_origin+"_"+'grads'
    # da = da.assign_coords(variable=variable_new)*10**3
    da = da.assign_coords(variable=variable_new)*10**3*(-1)
    da = da.transpose(*(...,'pressure'))  # 将pressure方位最后一个坐标
    return da
    # return dic_return
# cc = grads_data(ds)
# %%
# cc
# ds['ACM2']
# dd = cc.assign_coords(variable='theta_v_grads')
# dd.variable
## 最小是-410， 最大是68
# cc.max()

# %%
# %%
# %%
class Draw():
    def __init__(self, ):
        # self.da_var = da_var
        # self.name = name
        # self.station = station
        pass

    def draw_main(self, station_dic, var):

        for key in station_dic:
            pass
            station = station_dic[key]

            # 获得数据
            tr = TransferData(station)
            # model_dic = gd.get_data_main(var, station)
            model_dic = tr.transfer_data(var)  # 循环获得变量
            # print("yes")

            # 画图
            bb = self.combine_fig(var, model_dic, station['name'])

    def draw_contourf_single(self, ax, val):
        """[summary]

        Args:
            val (DataArray): 需要换图的变量
        """

        y = val.pressure.values
        # x = time_index
        x = val.time.values
        # time_index = time_index['data']
        # val = val.sel(time=time_index)
        var = str(val.coords['variable'].values)
        val = val.values.swapaxes(0, 1)  # 转置矩阵
        color_map = get_cmap_temp()


        if var == 'temp':
            level = np.arange(-20, 30, 2.5)
        elif var == 't_td':
            level = np.arange(0, 30, 1)
        elif var == 'wind_s':
            level = np.arange(0, 25, 2.5)
        elif var == 'temp_grads':
            level = np.arange(-0.3, 0.3, 0.02)
        elif var == 't_td_grads':
            # level = np.arange(0, 0.25, 0.01)
            level = np.arange(-0.25, 0.26, 0.01)
        elif var == 'wind_grads':
            level = np.arange(-0.25, 0.26, 0.05)
        elif var == 'q':
            level = np.arange(0.1, 10, 0.5)*10**(-3)
            color_map = get_cmap_q()
        elif var == 'theta_v':
            level = np.arange(320, 360, 2.5)
        elif var == 'theta_v_grads':
            # level = np.arange(-400, 100, 10)
            level = np.arange(-100, 104, 10)

        CS = ax.contourf(x,
                         y,
                         val,
                         levels=level,
                         cmap=color_map,
                         extend='max')
                        #  extend='neither')
        # cb = fig.colorbar(CS, orientation='horizontal', shrink=0.8, pad=0.14, fraction=0.14) # 这里的cs是画填色图返回的对象
        # 设置标签大小
        # ax.set_xticks(x[::2])  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        # xlabel = x[::2].dt.strftime('%m%d').values
        ax.set_xticks(x[::4])  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        # xlabel = x.dt.strftime('%d').values
        aa = pd.to_datetime(x[::4])
        xlabel = aa.strftime('%d%H')
        # xlabel = aa.strftime('%d').values
        ax.set_xticklabels(xlabel, rotation=45)
        ax.invert_yaxis()

        # plt.gca().invert_yaxis()
        ax.set_ylim(570, 400)
        # ax.set_title(title, fontsize=14)
        # ## 设置标签名称
        ax.set_xlabel("Time(UTC, day)", fontsize=14)
        ax.set_ylabel("Pressure (hPa)", fontsize=14)
        # fig.savefig(fig_name,bbox_inches = 'tight')
        return CS

    def draw_dual(self, ds, picture_dic):
        """
        画一个多图, 6张， 一个站点的
        Args:
            根据ds来绘图
            model_dic ([type]): 各试验数据的字典
            负责绘图，就只负责绘图就好了，可以不用负责保存图片
        """

        fig = plt.figure(figsize=(10, 15), dpi=200)  # 创建页面
        grid = plt.GridSpec(4,
                            2,
                            figure=fig,
                            left=0.10,
                            right=0.98,
                            bottom=0.15,
                            top=0.93,
                            wspace=0.3,
                            hspace=0.35)

        num = 7
        axes = [None] * num
        # axes = [None] * 14  # 设置一个维度为8的空列表
        for i in range(num):
            axes[i] = fig.add_subplot(grid[i])

        model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF', 'micaps', 'fnl']
        ax6 = fig.add_axes([0.2, 0.06, 0.7, 0.02])  # 重新生成一个新的坐标图

        var = str(ds.variable.values)
        time_select = str(ds.time[0].dt.hour.values)

        for i in range(len(model_list)):
            CS = None
                    # time_index = model_dic[model_list[i]].time.sel(time=datetime.time(12))
            CS = self.draw_contourf_single(axes[i],
                                            ds[model_list[i]], 
                                            )

            # title = str(station['name']) + "_" + model_list[i] + \
            #                 "_" + str(var)+"_"+str(time_select)
            # title = picture_dic['title'] + "_" + model_list[i]
            title = model_list[i]
            if model_list[i] == 'micaps':
                title = 'OBS'
            elif model_list[i] == 'fnl':
                title = 'FNL'
            axes[i].set_title(title, fontsize=14, loc='left')
        # plt.show()
            cb = fig.colorbar(CS,
                              cax=ax6,
                              orientation='horizontal',
                              shrink=0.8,
                              pad=0.14,
                              fraction=0.14)  # 这里的cs是画填色图返回的对象
        # plt.show()
        return fig

        # return ds

# station = station_dic['LaSa']
# station = station_dic['ShiQuanhe']
# station = station_dic['DuLan']
# # station = station_dic['GaiZe']
# # var = 'wind_s'
# var = 'theta_v'
# # var = 'q'
# Dr = Draw()
# aa = Dr.draw_dual(var,station, 12)

# Dr = Draw()
# picture_dic = {}
# time_select = '12'
# dds = cc
# time_index = ds.time.sel(time=datetime.time(int(time_select)))  
# dds = dds.sel(time=time_index)
# picture_dic['title'] = str(station['name']) + "_" + str(var)+"_"+str(time_select)
# fig = Dr.draw_dual(dds,picture_dic)







# %%
    
def main_station():
    """绘制单个站，单个要素的时间高度廓线
    一个站一个站的绘制
    """

    path = '/mnt/zfm_18T/fengxiang/Asses_PBL/UPAR/picture_time_sequence/'
    var_list = ['wind_s', 'theta_v', 'q']
    # time_select = '12'
    time_select_list = ['00', '12']

    picture_dic = {}
    
    for time_select in time_select_list:
        for var in var_list:
            for key in station_dic:
                
                ## 获取数据
                station = station_dic[key]
                ds = get_var(station, var)
                time_index = ds.time.sel(time=datetime.time(int(time_select)))  
                ds = ds.sel(time=time_index)
                
                ## 画图
                Dr = Draw()
                picture_dic['title'] = str(station['name']) + "_" + str(var)+"_"+str(time_select)
                fig = Dr.draw_dual(ds,picture_dic)
                fig_name = os.path.join(path, str(var) + "_" + str(station['name'])+"_"+str(time_select))
                fig.savefig(fig_name)


def get_station_land(var, time_select):
    """获得多个站点的平均数据，时间高度廓线

    Args:
        var (str): 不同的变量, q, theta_v
        time_select (str): 不同的时次

    Returns:
        [dict]: 三种下垫面类型的数据
    """
    
    low = ['ShiQuanhe', 'MangYa', 'GeErmu']
    medium = ['GaiZe']
    high  = ['TingRi', 'ShenZha', 'LaSa', 'NaQu', 'YuShu', 
              'DaRi', 'BaTang', 'LinZhi', 'ChangDu']
    Mean = ['ShiQuanhe', 'MangYa', 'GeErmu', 'TingRi', 'ShenZha', 'LaSa', 'NaQu', 'YuShu', 
              'DaRi', 'BaTang', 'LinZhi', 'ChangDu']

    landuse_list = [low, medium, high, Mean]
    rain_list = []

    for station_list in landuse_list:
        rain_station_list = []
        for key in station_list:  
            station = station_dic[key]
            ds = get_var(station, var)  # 单个站点，单个变量的数据
            time_index = ds.time.sel(time=datetime.time(int(time_select)))  
            ds = ds.sel(time=time_index)  # 12时或者00时的数据
            rain_station_list.append(ds)
            
        rain_station_array = xr.concat(rain_station_list, pd.Index(station_list, name='station'))
        rain_station_mean = rain_station_array.mean(dim='station')  # 多站的平均值
        rain_list.append(rain_station_mean)

    rain_land_dic = {}        
    rain_land_dic['low'] = rain_list[0]
    rain_land_dic['medium'] = rain_list[1]
    rain_land_dic['high'] = rain_list[2]
    rain_land_dic['mean'] = rain_list[3]
    return rain_land_dic


def main_land():
    """获得多个站点的平均时间高度廓线, 画图
    """
    pass
    path = '/mnt/zfm_18T/fengxiang/Asses_PBL/UPAR/picture_time_sequence/'
    # var_list = ['wind_s', 'theta_v', 'q']
    var_list = ['theta_v']
    # time_select = '00'
    # time_select = '12'
    time_select_list = ['00', '12']

    
    for time_select in time_select_list:
        picture_dic = {}
        for var in var_list:
            pass
            dic = get_station_land(var, time_select)
            for key in dic:
                ds = dic[key]
                ## 求梯度
                ds = grads_data(ds)
                print(ds)
                ##
                Dr = Draw()
                picture_dic['title'] = str(key) + "_" + str(var)+"_"+str(time_select)
                fig = Dr.draw_dual(ds,picture_dic)
                fig_name = os.path.join(path, str(var) + "_" + str(key)+"_"+str(time_select))

                ## 设置图片标题
                tt_tile = str(key)                
                if str(key) == 'low':
                    tt_tile = 'BareLand'
                elif str(key) == 'medium':
                    tt_tile = 'Bush'
                elif str(key) == 'high':
                    tt_tile = 'GrassLand'
                
                var = str(ds.variable.values)
                fig_title = tt_tile+"_" + str(var)+"_"+str(time_select)
                fig.suptitle(fig_title, fontsize=26)
                fig.savefig(fig_name)




# %%

if __name__ == '__main__':
    # main_station()
    main_land()

    pass

    # 最终画图
    # var_list = ['temp', 't_td', 'wind', 'temp_grads',
    #             't_td_grads', 'wind_grads', 'theta_v', 'q']
    # var = var_list[-1]

    # Dr = Draw()
    # Dr.draw_main(station_dic, var)
