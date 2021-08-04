#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
画时间序列图, 各区域降水量24小时的变化
-----------------------------------------
Time             :2021/06/04 14:32:20
Author           :Forxd
Version          :1.0
'''
    

from re import T
import xarray as xr
from read_data import TransferData, GetData
import meteva.method as mem
import meteva.base as meb
import numpy as np
import pandas as pd

import salem  # 过滤高原外的数据
import geopandas

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.pyplot import savefig

from global_variable import station_dic
import datetime



# %%





class Draw():
    
    def __init__(self, month):
        self.fontsize = 10
        self.month = month

    def draw_time_sequence(self,ax, dic):

        ccolor = ['black','red', 'cyan', 'green', 'blue', 'orange', ]
        QNSE = dic['QNSE']
        x_label = QNSE.coords['time']
        # x_label = str(x_label.dt.strftime('%d%H').values).split()
        x_label = x_label.dt.strftime('%d%H')


        # x_label = x_label.dt.strftime("%H")  # 转换时间维字符串格式
        y = QNSE.values
        module_list = ['obs', 'ACM2','YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        x = np.arange(len(x_label))
        # print(x)

        # fig = plt.figure(figsize=(8, 5), dpi=400)  # 创建页面
        # ax = fig.add_axes([0.1, 0.13, 0.85, 0.8])  # 重新生成一个新的坐标图
        j = 0
        for i in module_list:
            # y = dr.loc[i,:].values
            y = dic[i].values
            y = np.around(y,2)
            ax.plot(x_label, y, label=i, color=ccolor[j])
            j +=1 

        ax.set_xticks(x_label[::2])  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        # ax.xaxis.set_tick_params(labelsize=15)
        ax.xaxis.set_tick_params(labelsize=self.fontsize*1.8, rotation=45)
        ax.tick_params(which='major',length=8,width=1.0) # 控制标签大小 
        ax.tick_params(which='minor',length=4,width=0.5)  #,colors='b')
        ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
        # ax.set_yticks(np.arange(0, 5.01, 0.1))  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax.yaxis.set_tick_params(labelsize=self.fontsize*1.8)
        # ax.set_ylim(0,5.1)
        # ax.set_ylim(0.0,0.6)
        ax.set_xlabel("Time(UTC)", fontsize=self.fontsize*2.0)
        ax.set_ylabel("Precipitation (mm)", fontsize=self.fontsize*2.0)
        # ax.legend()
        # ax.set_title("201607", fontsize=18)
        # fig.savefig('/home/fengxiang/Project/Asses_pbl_July/Draw/Rain/time_sequecnce.png')
    
    def combine_fig(self, tr):

        fig = plt.figure(figsize=(16, 20), dpi=200)  # 创建页面
        grid = plt.GridSpec(4,
                            1,
                            figure=fig,
                            left=0.05,
                            right=0.98,
                            bottom=0.1,
                            top=0.97,
                            wspace=0.2,
                            hspace=0.25)

        axes = [None] * 6  # 设置一个维度为8的空列表
        axes[0] = fig.add_subplot(grid[0])
        axes[1] = fig.add_subplot(grid[1])
        axes[2] = fig.add_subplot(grid[2])
        axes[3] = fig.add_subplot(grid[3])

        time_flag=all
        dic = {}

        station_list = ['GaiZe', ]

        station_dic1 = {
                'ShiQuanhe':station_dic['ShiQuanhe'],
                'GaiZe':station_dic['GaiZe'], 
                'ShenZha':station_dic['ShenZha'],
                'TingRi':station_dic['TingRi'],
                    }        
        
        station_dic2 = {
                'TuoTuohe':station_dic['TuoTuohe'],
                'NaQu':station_dic['NaQu'], 
                'LaSa':station_dic['LaSa'],
                'LinZhi':station_dic['LinZhi'],
                    }        
        station_dic3 = {
                'ChangDu':station_dic['ChangDu'],
                'MangYa':station_dic['MangYa'], 
                'GeErmu':station_dic['GeErmu'],
                'DuLan':station_dic['DuLan'],
                    }        

        i = 0
        station_dic_dic = station_dic3  # 这里要改
        for key in station_dic_dic:  
            rain1 = tr.rain_station(station_dic_dic[key])  
            dic[key] = rain1
            axes[i].set_title(key, fontsize=self.fontsize*2.0, loc='left', y=0.88, x=0.05)
            i += 1

        for i,j in zip(range(4),dic):
            self.draw_time_sequence(axes[i], dic[j])
            axes[i].set_ylim(0.0, 10.0)
            axes[i].set_yticks(np.arange(0, 10.1, 2.0))
            
        axes[3].legend(ncol=3 ,bbox_to_anchor=(0.5,-0.55) ,loc='lower center',fontsize=self.fontsize*2.0, edgecolor='white')
        # fig.suptitle("May", fontsize=self.fontsize*2.5)
        fig.suptitle(self.month, fontsize=self.fontsize*2.5)
        flnm = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/rain_staion3_'+self.month+'.png'   # 这里要改
        fig.savefig(flnm)



if __name__ == '__main__':
    area = {"lat1": 24.875, "lat2": 40.125, "lon1": 69.875, "lon2": 105.125}
    month = 'Jul'
    # month = 'May'
    gd = GetData(month)
    rain = gd.get_rain_hourly()

    time_flag = 'all'
    flag = 'all'
    da_obs = rain['obs']
    ## 选取最大范围和全部时间，保证能取出所需要的数据
    time_index_12 = da_obs.time.sel(time=datetime.time(int('12')))  
    time_index_00 = da_obs.time.sel(time=datetime.time(int('00')))  
    time_index = np.union1d(time_index_00.values, time_index_12.values)
    rain2 = rain.sel(time=time_index)


    # %%
    tr = TransferData(ds=rain2, area=area, time_flag=time_flag)
    Dr = Draw(month)
    Dr.combine_fig(tr)

# %%
