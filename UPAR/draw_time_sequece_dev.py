#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取单个站点的所有时次，高度，方案的数据
绘制它的时间高度廓线
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
    month = 'Jul'
    time_hour = '12'
    tr = TransferData(station, month, time_hour, 'q')
    bb = tr.get_data_all()
    var_return = bb.sel(variable=var)
    return var_return

# station = station_dic['GaiZe']
# var = 'wind_s'
# aa = get_var(station, var)
















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

    def draw_contourf_single(self, var, ax, val, title):
        """[summary]

        Args:
            val (DataArray): 需要换图的变量
        """

        y = val.pressure.values
        # x = time_index
        x = val.time.values
        # time_index = time_index['data']
        # val = val.sel(time=time_index)
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
            level = np.arange(330, 360, 2.5)

        CS = ax.contourf(x,
                         y,
                         val,
                         levels=level,
                         cmap=color_map,
                         extend='both')
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
        ax.set_ylim(570, 300)
        ax.set_title(title, fontsize=14)
        # ## 设置标签名称
        ax.set_xlabel("Time(UTC, day)", fontsize=14)
        ax.set_ylabel("Pressure (hPa)", fontsize=14)
        # fig.savefig(fig_name,bbox_inches = 'tight')
        return CS

    def draw_dual(self, var, station, time_select):
        """画一个多图, 6张， 一个站点的
        Args:
            model_dic ([type]): 各试验数据的字典
        """
        # fig = plt.figure(figsize=(8, 5), dpi=400)  # 创建页面
        # ax = fig.add_axes([0.1, 0.13, 0.85, 0.8])  # 重新生成一个新的坐标图
        # print(area_dic['all'])


        fig = plt.figure(figsize=(10, 15), dpi=200)  # 创建页面
        grid = plt.GridSpec(4,
                            2,
                            figure=fig,
                            left=0.10,
                            right=0.98,
                            bottom=0.15,
                            top=0.95,
                            wspace=0.3,
                            hspace=0.35)

        num = 7
        axes = [None] * num
        # axes = [None] * 14  # 设置一个维度为8的空列表
        for i in range(num):
            axes[i] = fig.add_subplot(grid[i])

        # station = station_dic['LaSa']
        # var = 'wind_s'
        ds = get_var(station, var)
        time_index = ds.time.sel(time=datetime.time(int(time_select)))  
        # time_index_00 = ds.time.sel(time=datetime.time(int('00')))  
        # time_index = np.union1d(time_index_00.values, time_index_12.values)
        # time_index = time_index_12
        ds = ds.sel(time=time_index)
        

        model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF', 'micaps', 'fnl']

        ax6 = fig.add_axes([0.2, 0.06, 0.7, 0.02])  # 重新生成一个新的坐标图
        # keys = ['00', '06', '12']
        # for model in model_list:
        # for key in keys:
        for i in range(len(model_list)):
            CS = None
                    # time_index = model_dic[model_list[i]].time.sel(time=datetime.time(12))
            title = str(station['name']) + "_" + model_list[i] + \
                            "_" + str(var)+"_"+str(time_select)
            CS = self.draw_contourf_single(var, axes[i],
                                            ds[model_list[i]], title,
                                            )

        # plt.show()
            cb = fig.colorbar(CS,
                              cax=ax6,
                              orientation='horizontal',
                              shrink=0.8,
                              pad=0.14,
                              fraction=0.14)  # 这里的cs是画填色图返回的对象

        path = '/mnt/zfm_18T/fengxiang/Asses_PBL/UPAR/picture_time_sequence/'
        fig_name = os.path.join(
            path,
            str(var) + "_" + str(station['name'])+"_"+str(time_select))
        # fig.savefig('/mnt/zfm_18T/fengxiang/Asses_PBL/tt.png')
        fig.savefig(fig_name)
        plt.show()
        return ds

# station = station_dic['LaSa']
# station = station_dic['ShiQuanhe']
# station = station_dic['DuLan']
# # station = station_dic['GaiZe']
# # var = 'wind_s'
# var = 'theta_v'
# # var = 'q'
# Dr = Draw()
# aa = Dr.draw_dual(var,station, 12)


def main():
    # station = station_dic['DuLan']
    # station = station_dic['GaiZe']
    # var = 'wind_s'
    var_list = ['wind_s', 'theta_v', 'q']
    time_select = '00'
    
    for var in var_list:

        for key in station_dic:
            station = station_dic[key]
            # var = 'theta_v'
            # var = 'q'
            # var = 'wind_s'
            Dr = Draw()
            aa = Dr.draw_dual(var,station, time_select)
    




# %%

if __name__ == '__main__':
    main()

    pass

    # 最终画图
    # var_list = ['temp', 't_td', 'wind', 'temp_grads',
    #             't_td_grads', 'wind_grads', 'theta_v', 'q']
    # var = var_list[-1]

    # Dr = Draw()
    # Dr.draw_main(station_dic, var)
