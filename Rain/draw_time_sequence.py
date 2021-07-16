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
    

import xarray as xr
from read_data  import (get_rain_total, get_rain_times,  
        get_rain_24hour_sequence, get_rain_3hour_sequence)
import meteva.method as mem
import meteva.base as meb
import numpy as np
import pandas as pd

import salem  # 过滤高原外的数据
import geopandas

# import matplotlib.pyplot as plot
import matplotlib.pyplot as plt
from matplotlib.pyplot import savefig
# get_data = Get_data('night')


class Draw():

    def draw_time_sequence(self,ax, dic):

        ccolor = ['black','red', 'cyan', 'green', 'blue', 'orange', ]
        QNSE = dic['QNSE']
        x_label = QNSE.coords['time']

        # x_label = x_label.dt.strftime("%H")  # 转换时间维字符串格式
        y = QNSE.values
        module_list = ['obs', 'MYJ','YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        x = np.arange(len(x_label))

        # fig = plt.figure(figsize=(8, 5), dpi=400)  # 创建页面
        # ax = fig.add_axes([0.1, 0.13, 0.85, 0.8])  # 重新生成一个新的坐标图
        j = 0
        for i in module_list:
            # y = dr.loc[i,:].values
            y = dic[i].values
            y = np.around(y,2)
            ax.plot(x_label, y, label=i, color=ccolor[j])
            j +=1 

        ax.set_xticks(x_label)  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax.xaxis.set_tick_params(labelsize=15)
        ax.set_yticks(np.arange(0, 5.01, 0.5))  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax.yaxis.set_tick_params(labelsize=15)
        ax.set_ylim(0,5.1)
        ax.set_xlabel("Time(UTC)", fontsize=18)
        ax.set_ylabel("Precipitation (mm)", fontsize=18)
        # ax.legend()
        # ax.set_title("201607", fontsize=18)
        # fig.savefig('/home/fengxiang/Project/Asses_pbl_July/Draw/Rain/time_sequecnce.png')
    
    def combine_fig(self, area_dic):
        # fig = plt.figure(figsize=(8, 5), dpi=400)  # 创建页面
        # ax = fig.add_axes([0.1, 0.13, 0.85, 0.8])  # 重新生成一个新的坐标图
        # print(area_dic['all'])
        fig = plt.figure(figsize=(10, 8), dpi=200)  # 创建页面
        grid = plt.GridSpec(2,
                            2,
                            figure=fig,
                            left=0.10,
                            right=0.98,
                            bottom=0.18,
                            top=0.95,
                            wspace=0.2,
                            hspace=0.21)

        axes = [None] * 6  # 设置一个维度为8的空列表
        axes[0] = fig.add_subplot(grid[0, 0:1])
        axes[1] = fig.add_subplot(grid[0, 1:2])
        axes[2] = fig.add_subplot(grid[1, 0:1])
        axes[3] = fig.add_subplot(grid[1, 1:2])


        axes[0].set_title("all", fontsize=18, loc='left', y=0.88, x=0.05)
        axes[1].set_title("north", fontsize=18, loc='left', y=0.88, x=0.05)
        axes[2].set_title("south west", fontsize=18, loc='left', y=0.88, x=0.05)
        axes[3].set_title("south east", fontsize=18, loc='left', y=0.88, x=0.05)
        # axes[4] = fig.add_subplot(grid[2, 0:1], projection=proj)
        # axes[5] = fig.add_subplot(grid[2, 1:2], projection=proj)
        # # area = area_list[0]
        # dic0 = area_dic['all']
        # dic0 = get_rain_24hour_sequence(area_dic['all'])
        # dic1 = get_rain_24hour_sequence(area_dic['north'])
        # dic2 = get_rain_24hour_sequence(area_dic['south left'])
        # dic3 = get_rain_24hour_sequence(area_dic['south right'])

        dic0 = get_rain_3hour_sequence(area_dic['all'])
        dic1 = get_rain_3hour_sequence(area_dic['north'])
        dic2 = get_rain_3hour_sequence(area_dic['south left'])
        dic3 = get_rain_3hour_sequence(area_dic['south right'])
        # ax = fig.add_axes([0.1, 0.13, 0.85, 0.8])  # 重新生成一个新的坐标图
        # print(dic0)

        # self.draw_time_sequence(ax, dic0)
        self.draw_time_sequence(axes[0], dic0)
        self.draw_time_sequence(axes[1], dic1)
        self.draw_time_sequence(axes[2], dic2)
        self.draw_time_sequence(axes[3], dic3)
        axes[2].legend(ncol=3 ,bbox_to_anchor=(1.0,-0.55) ,loc='lower center',fontsize=16,edgecolor='white')
        fig.savefig('/home/fengxiang/Project/Asses_PBL/Draw/Rain/time_sequecnce.png')



if __name__ == '__main__':
    area0 = {"lat1":24.875, "lat2":46.125, "lon1":70.875, "lon2":105.125}
    area1 = {"lat1":33.5, "lat2":40, "lon1":80, "lon2":90}  # north
    area2 = {"lat1":28, "lat2":33, "lon1":83, "lon2":94}  # south left 
    area3 = {"lat1":26, "lat2":33, "lon1":95, "lon2":103}  # south right
    area4 = {"lat1":27, "lat2":28, "lon1":86, "lon2":92}  # bottom
    
    ## 选择时间范围和计算区域
    flag = 'all'

    area_dic = {}
    area_dic['all'] = area0
    area_dic['north'] = area1
    area_dic['south left'] = area2
    area_dic['south right'] = area3
    # print(area_dic)
    Dr = Draw()
    # Dr.draw_time_sequence()
    Dr.combine_fig(area_dic)
    # pass




    


