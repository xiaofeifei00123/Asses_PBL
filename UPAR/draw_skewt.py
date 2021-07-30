#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
关注相对湿度的廓线
-----------------------------------------
Time             :2021/06/17 18:59:18
Author           :Forxd
Version          :2.0
'''

import os
import numpy as np
from netCDF4 import Dataset
from wrf import getvar, vinterp, interplevel
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import cmaps
from global_cmap import get_cmap_temp
import datetime

from data_process import TransferData
from global_variable import station_dic


class Draw_skewt():

    def __init__(self, station, month):
        pass
        self.station = station
        self.month = month
        if self.month == 'Jul':
            if station['name'] == 'ShenZha':
                self.model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'obs']
                # self.color_list = ['orange', 'red', 'cyan', 'blue', 'black']
                self.color_list = [
                    'red', 'red', 'blue', 'blue', 'green', 'black'
                ]
            else:
                self.model_list = [
                    'ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF', 'obs'
                ]
                self.color_list = [
                    'red', 'red', 'blue', 'blue', 'green', 'black'
                ]
        else:
            self.model_list = [
                'ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF', 'obs'
            ]
            self.color_list = [
                'red', 'red', 'blue', 'blue', 'green', 'black'
            ]

        self.line_style_list = ['solid', 'dashed','solid', 'dashed','solid', 'solid', 'solid']

        self.path = '/mnt/zfm_18T/fengxiang/Asses_PBL/UPAR/picture_upar_q'


    def draw_upar_single(self, model_dic, title_dic):
        """画探空的廓线
        """
        pass
        # print()
        fig = plt.figure(figsize=(4, 7), dpi=400)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.85])  # 左下起始位置，宽度和高度
        # ax2 = fig.add_axes([0.42, 0.15, 0.22, 0.8])
        # ax3 = fig.add_axes([0.75, 0.15, 0.22, 0.8])

        model_list = self.model_list
        # model_q = {}
        # model_wind = {}
        model_var = {}

        ## 筛选时间
        time_index1 = model_dic['obs'].time.sel(
            time=datetime.time(int(title_dic['time'])))  

        time_index2 = model_dic['TEMF'].time.sel(
            time=datetime.time(int(title_dic['time'])))  

        ## 求两个时间的交集
        time_index = np.intersect1d(time_index1.values,
                                    time_index2.values)
        for i in range(len(model_list)):

            model_var[model_list[i]] = \
                model_dic[model_list[i]].sel(
                    time=time_index).mean(dim='time', skipna=True)


            model_var['fnl'] = \
                model_dic['fnl'].sel(
                    time=time_index).mean(dim='time', skipna=True)


        ## 画不同的变量
        # model_t = model_var
        # model_var = model_q

        for [model, i, j] in zip(model_list, self.color_list, self.line_style_list):

            # y = model_t[model].pressure
            # x = model_t[model].values
            # ax1.plot(x, y, color=i, label=model, linestyle=j)

            y2 = model_var[model].pressure
            x2 = model_var[model].values
            ax.plot(x2, y2, color=i, label=model, linestyle=j)

            # y3 = model_wind[model].pressure
            # x3 = model_wind[model].values
            # ax3.plot(x3, y3, color=i, label=model, linestyle=j)
            # ax3.set_xlim(0, 25)

        fts = 14
        # ax.set_ylabel('Pressure/(hPa)', fontsize=fts)
        # ax.set_xlabel('Theta_v/(K)', fontsize=fts)
        # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        # ax1.set_xticks(np.arange(320, 355, 10))
        # ax1.set_xlim(320, 360)

        x2_fnl = model_var['fnl'].values
        y2_fnl = model_var['fnl'].pressure
        ax.plot(x2_fnl, y2_fnl, color='cyan', label='fnl', linestyle='solid')
        ax.set_xticks(np.arange(0, 10, 1))
        # ax.legend(loc='lower center',
        #            bbox_to_anchor=(0.35, -0.18, 0.2, 0.2),
        #            ncol=4)
        ax.set_xlabel('q  (g/kg)', fontsize=fts)

        # ax3.set_xlabel('Wind_speed (m/s)', fontsize=fts)
        # ax3.set_xticks(np.arange(0, 26, 5))

        # ax1.invert_yaxis()
        ax.invert_yaxis()
        # ax3.invert_yaxis()
        # ax1.set_ylim(570, 300)
        ax.set_ylim(570, 300)
        # ax3.set_ylim(570, 300)

        fig.suptitle(
            title_dic['station_name'] + "_" + self.month+"_"+title_dic['time'],
            fontsize=fts * 1.3)

        fig_name = title_dic['station_name'] + "_" +self.month+"_"+ title_dic['time']+'_single'
        fgnm = os.path.join(self.path, fig_name)
        fig.savefig(fgnm)


    def draw_upar(self, model_dic, model_dic_q,
                  model_dic_wind, title_dic):
        """画探空的廓线
        """
        pass
        # print()
        fig = plt.figure(figsize=(9, 7), dpi=400)
        ax1 = fig.add_axes([0.1, 0.15, 0.22, 0.8])  # 左下起始位置，宽度和高度
        ax2 = fig.add_axes([0.42, 0.15, 0.22, 0.8])
        ax3 = fig.add_axes([0.75, 0.15, 0.22, 0.8])

        model_list = self.model_list
        model_q = {}
        model_wind = {}
        model_var = {}

        ## 筛选时间
        time_index1 = model_dic['obs'].time.sel(
            time=datetime.time(int(title_dic['time'])))  

        time_index2 = model_dic['TEMF'].time.sel(
            time=datetime.time(int(title_dic['time'])))  

        ## 求两个时间的交集
        time_index = np.intersect1d(time_index1.values,
                                    time_index2.values)
        for i in range(len(model_list)):

            model_var[model_list[i]] = \
                model_dic[model_list[i]].sel(
                    time=time_index).mean(dim='time', skipna=True)

            model_q[model_list[i]] = \
                model_dic_q[model_list[i]].sel(
                    time=time_index).mean(dim='time', skipna=True)

            model_q['fnl'] = \
                model_dic_q['fnl'].sel(
                    time=time_index).mean(dim='time', skipna=True)

            model_wind[model_list[i]] = \
                model_dic_wind[model_list[i]].sel(
                    time=time_index).mean(dim='time', skipna=True)

        ## 画不同的变量
        model_t = model_var
        model_var = model_q

        for [model, i, j] in zip(model_list, self.color_list, self.line_style_list):

            y = model_t[model].pressure
            x = model_t[model].values
            ax1.plot(x, y, color=i, label=model, linestyle=j)

            y2 = model_var[model].pressure
            x2 = model_var[model].values
            ax2.plot(x2, y2, color=i, label=model, linestyle=j)

            y3 = model_wind[model].pressure
            x3 = model_wind[model].values
            ax3.plot(x3, y3, color=i, label=model, linestyle=j)
            ax3.set_xlim(0, 25)

        fts = 14
        ax1.set_ylabel('Pressure/(hPa)', fontsize=fts)
        ax1.set_xlabel('Theta_v/(K)', fontsize=fts)
        # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax1.set_xticks(np.arange(320, 355, 10))
        ax1.set_xlim(320, 360)

        x2_fnl = model_var['fnl'].values
        y2_fnl = model_var['fnl'].pressure
        ax2.plot(x2_fnl, y2_fnl, color='cyan', label='fnl', linestyle='solid')
        ax2.set_xticks(np.arange(0, 10, 1))
        ax2.legend(loc='lower center',
                   bbox_to_anchor=(0.25, -0.18, 0.2, 0.2),
                   ncol=4)
        ax2.set_xlabel('q  (g/kg)', fontsize=fts)

        ax3.set_xlabel('Wind_speed (m/s)', fontsize=fts)
        ax3.set_xticks(np.arange(0, 26, 5))

        ax1.invert_yaxis()
        ax2.invert_yaxis()
        ax3.invert_yaxis()
        ax1.set_ylim(570, 300)
        ax2.set_ylim(570, 300)
        ax3.set_ylim(570, 300)

        fig.suptitle(
            title_dic['station_name'] + "_" + self.month+"_"+title_dic['time'],
            fontsize=fts * 1.3)

        fig_name = title_dic['station_name'] + "_" +self.month+"_"+ title_dic['time']
        fgnm = os.path.join(self.path, fig_name)
        fig.savefig(fgnm)

    def draw_main(self, ):

        # time_index = ['00', '06', '12']
        time_index = ['00', '12']
        # time_index = ['12']

        # 循环出一个时次一个站点
        for time_select in time_index:
            ## 获得数据
            tr = TransferData(self.station, self.month)
            # model_dic = tr.transfer_data('theta_v')
            model_dic_q = tr.transfer_data('q')
            # model_dic_wind = tr.transfer_data('wind_s')
            ## 画图
            title = {'time': time_select, 'station_name': self.station['name']}
            self.draw_upar_single(
                # model_dic,
                model_dic_q,
                # model_dic_wind,
                title,
            )
            # self.draw_upar(
            #     model_dic,
            #     model_dic_q,
            #     model_dic_wind,
            #     title,
            # )

if __name__ == '__main__':

    ## 测试

    # %%
    pass
    # %%
    station = station_dic['GaiZe']
    month = 'Jul'
    time_select = '00'
    # %%
    tr = TransferData(station, month)
    model_dic_q = tr.transfer_data('q')
    # %%
    time_select = '00'
    title = {'time': time_select, 'station_name': station['name']}
    dr = Draw_skewt(station, month)
    dr.draw_upar_single(
        model_dic_q,
        title,
        )


'''     ## 画图
    for key in station_dic:
        station = station_dic[key]
        # print(station)
        Dr = Draw_skewt(station, 'Jul')
        # Dr = Draw_skewt(station, 'May')
        Dr.draw_main() '''

# %%
