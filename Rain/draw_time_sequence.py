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

# import matplotlib.pyplot as plot
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.pyplot import savefig
# get_data = Get_data('night')

from global_variable import station_dic


class Draw():
    
    def __init__(self):
        self.fontsize = 10

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

        ax.set_xticks(x_label[::24])  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        # ax.xaxis.set_tick_params(labelsize=15)
        ax.xaxis.set_tick_params(labelsize=self.fontsize*1.8)
        # ax.xaxis.set_major_formatter(x_label[::24])
        ax.xaxis.set_minor_locator(plt.MultipleLocator(4))
        ax.tick_params(which='major',length=8,width=1.0) # 控制标签大小 
        ax.tick_params(which='minor',length=4,width=0.5)  #,colors='b')
        # ax.xaxis.set_minor_locator(x_label.values[::12])
        # ax.set_yticks(np.arange(0, 5.01, 0.1))  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax.yaxis.set_tick_params(labelsize=self.fontsize*1.8)
        # ax.set_ylim(0,5.1)
        # ax.set_ylim(0.0,0.6)
        ax.set_xlabel("Time(UTC)", fontsize=self.fontsize*2.0)
        ax.set_ylabel("Precipitation (mm)", fontsize=self.fontsize*2.0)
        # ax.legend()
        # ax.set_title("201607", fontsize=18)
        # fig.savefig('/home/fengxiang/Project/Asses_pbl_July/Draw/Rain/time_sequecnce.png')
    
    def combine_fig(self, area_dic, rain):
        fig = plt.figure(figsize=(20, 16), dpi=200)  # 创建页面
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

        axes[0].set_title("all", fontsize=self.fontsize*2.0, loc='left', y=0.88, x=0.05)
        axes[1].set_title("north", fontsize=self.fontsize*2.0, loc='left', y=0.88, x=0.05)
        axes[2].set_title("south west", fontsize=self.fontsize*2.0, loc='left', y=0.88, x=0.05)
        axes[3].set_title("south east", fontsize=self.fontsize*2.0, loc='left', y=0.88, x=0.05)
        time_flag=all
        dic = {}
        for key in area_dic:
            # print(key)
            tr = TransferData(ds=rain, area=area_dic[key], time_flag=time_flag)
            rain1 = tr.rain_space_average()
            dic[key] = rain1

        for i,j in zip(range(4),dic):
            self.draw_time_sequence(axes[i], dic[j])
        axes[0].set_ylim(0.0, 0.6)
        axes[1].set_ylim(0.0, 0.6)
        axes[2].set_ylim(0.0, 0.6)
        axes[3].set_ylim(0.0, 2.0)
        axes[0].set_yticks(np.arange(0, 0.7, 0.1)) 
        axes[1].set_yticks(np.arange(0, 0.7, 0.1))
        axes[2].set_yticks(np.arange(0, 0.7, 0.1))
        axes[3].set_yticks(np.arange(0, 2.1, 0.2))
            
        axes[3].legend(ncol=3 ,bbox_to_anchor=(0.5,-0.55) ,loc='lower center',fontsize=self.fontsize*2.0, edgecolor='white')
        fig.suptitle("May", fontsize=self.fontsize*2.5)
        fig.savefig('/mnt/zfm_18T/fengxiang/Asses_PBL/Rain_May/picture/rain_time.png')



if __name__ == '__main__':
    # area0 = {"lat1":24.875, "lat2":46.125, "lon1":70.875, "lon2":105.125}
    area0 = {"lat1": 24.875, "lat2": 40.125, "lon1": 69.875, "lon2": 105.125}
    area1 = {"lat1":33.5, "lat2":40, "lon1":80, "lon2":90}  # north
    area2 = {"lat1":28, "lat2":33, "lon1":83, "lon2":94}  # south left 
    area3 = {"lat1":26, "lat2":33, "lon1":95, "lon2":103}  # south right
    area4 = {"lat1":27, "lat2":28, "lon1":86, "lon2":92}  # bottom
    
    ## 选择时间范围和计算区域
    flag = 'all'
    # %%

    area_dic = {}
    area_dic['all'] = area0
    area_dic['north'] = area1
    area_dic['south left'] = area2
    area_dic['south right'] = area3
    # print(area_dic)
    gd = GetData()
    rain = gd.get_rain_hourly()

    time_flag = 'all'
    # tr = TransferData()    
    # tr = TransferData(ds=rain, area=area_dic['north'], time_flag=time_flag)
    # aa = tr.rain_space_average()
    # print(aa)
    # %%
    

    Dr = Draw()
    # Dr.draw_time_sequence()
    Dr.combine_fig(area_dic, rain)
    # pass




    



# %%
