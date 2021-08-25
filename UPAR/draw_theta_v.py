#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
关注虚位温的廓线
-----------------------------------------
Time             :2021/06/17 18:59:18
Author           :Forxd
Version          :2.0
'''

# %%
import os
from threading import Condition
import numpy as np
from netCDF4 import Dataset
from wrf import getvar, vinterp, interplevel
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from cycler import cycler
import cmaps
from global_cmap import get_cmap_temp
import datetime

# from data_process import TransferData
from data_transfer import TransferData
from global_variable import station_dic



# %%
# # %%
# station = station_dic['ShiQuanhe']
# # station = station_dic['GeErmu']
# month = 'Jul'
# # condition = 'rain'
# condition = 'all'
# time_select = '12'
# tr1 = TransferData(station, month, time_select, 'theta_v')
# ds_q = tr1.transfer_data(condition)
# print(ds_q)
# # %%





# %%
def get_all_data():        
    """获取所有站点的数据

    Returns:
        [type]: [description]
    """
    i = 0
    # station_dic_dic = station_dic  # 这里要改
    dic = {}
    month = 'Jul'
    condition = 'all'
    time_select = '12'
    # tr = TransferData(station, month, time_select, 'q')
    fontsize = 10
    for key in station_dic:  
        station = station_dic[key]
        # rain1 = tr.rain_station(station_dic_dic[key])  
        # r = rain.sel(lat=station['lat'], lon=station['lon'], method='nearest')
        tr = TransferData(station, month, time_select, 'theta_v')
        ds_q = tr.transfer_data(condition)

        dic[key] = ds_q
        # axes[i].set_title(key, fontsize=fontsize*2.0, loc='left', y=0.88, x=0.05)
        i += 1
    return dic
        
dic = get_all_data()
dic_theta_v = dic


# %%
# dic_theta_v = ds

# %%
class DrawSkewt():

    def __init__(self) -> None:
        pass
        self.path = '/mnt/zfm_18T/fengxiang/Asses_PBL/UPAR/picture_upar_q/'
    
    def draw_upar_single(self, ax, ds):
        """画探空的廓线
        """
        ccolor = ['black', 'cyan','red', 'red', 'blue', 'blue', 'green']
        model_list = ['micaps', 'fnl', 'ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        # ccolor = ['black','red', 'cyan', 'green', 'blue', 'orange', ]
        lline_style = ['-', '-', '--', '-', '--', '-.', '-']
        mmarker = ['o', '^', '^', '*', '*', '+', 'o']
        custom_cycler = (
            cycler(color=ccolor) +
            cycler(linestyle=lline_style) +           
            # cycler(label=module_list)
            cycler(marker=mmarker)
                        )
        # ds = ds*10**3
        # fig = plt.figure(figsize=(4, 7), dpi=400)
        # ax = fig.add_axes([0.1, 0.1, 0.8, 0.85])  # 左下起始位置，宽度和高度

        model_var = ds.mean(dim='time')
        ## 画不同的变量
        ax.set_prop_cycle(custom_cycler)
        for model in model_list:
            y2 = model_var[model].pressure
            x2 = model_var[model].values
            if model == 'micaps':
                model = 'OBS'
            ax.plot(x2[::2], y2[::2], label=model, lw=2, ms=4)
        fts = 18
        ax.invert_yaxis()
        ax.set_ylim(570, 300)
        # ax.set_xlim(0, 10)

        ax.xaxis.set_tick_params(labelsize=fts*1.3)
        ax.yaxis.set_tick_params(labelsize=fts*1.3)
        ax.tick_params(which='major',length=8,width=1.0) # 控制标签大小 
        ax.tick_params(which='minor',length=4,width=0.5)  #,colors='b')
        ax.yaxis.set_minor_locator(plt.MultipleLocator(10))
        ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
        ax.set_xlabel("theta (K)", fontsize=fts*2.0)
        ax.set_ylabel("Pressure (hPa)", fontsize=fts*2.0)
        # ax.legend()
        # plt.show()
        
        
        

    def combine_fig(self, dic):
        fig = plt.figure(figsize=(18, 25), dpi=200)  # 创建页面
        grid = plt.GridSpec(5,
                            3,
                            figure=fig,
                            left=0.05,
                            right=0.98,
                            bottom=0.12,
                            top=0.97,
                            wspace=0.2,
                            # hspace=0.25)
                            hspace=0.3)

        # num = 15
        num = len(station_dic)
        axes = [None] * num # 设置一个num维度为8的空列表
        # axes = [None] * 14  
        for i in range(num):
            axes[i] = fig.add_subplot(grid[i])

        fontsize  = 16
        for i,j in zip(range(num),dic):
            if j == 'micaps':
                j = 'OBS'
            axes[i].set_title(j, fontsize=fontsize*2.0, loc='left', y=1.01, x=0.05)
            self.draw_upar_single(axes[i], dic[j])
            print(dic[j])

            
        # i = 0            
        # for key in station_dic:  
        #     station = station_dic[key]
        #     if j == 'micaps':
        #         j = 'OBS'
        #     axes[i].set_title(j, fontsize=fontsize*2.0, loc='left', y=1.01, x=0.05)
        #     self.draw_upar_single(axes[i], dic[key])
        #     # print(dic[j])
        #     i += 1

            

        # # plt.show()
            
        # axes[7].legend(ncol=4 ,bbox_to_anchor=(0.5,-0.55) ,loc='lower center',fontsize=fontsize*2.0, edgecolor='white')
        # # flnm = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/rain_staion1_'+self.month+'.png'   # 这里要改
        # flnm = self.path+'theta_v'+'.png'   # 这里要改
        # fig.savefig(flnm)
        

dr = DrawSkewt()
dr.combine_fig(dic_theta_v)



# %%
