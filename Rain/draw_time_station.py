#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
日降水变化曲线
画时间序列图, 各区域降水量24小时的变化
逐日降水的时间变化曲线
这里的时间是从前一天的
-----------------------------------------
Time             :2021/06/04 14:32:20
Author          :Forxd
Version          :1.0
'''

# %%
# from re import T
import xarray as xr
from read_data import TransferData, GetData
# import meteva.method as mem
# import meteva.base as meb
import numpy as np
import pandas as pd

# import salem  # 过滤高原外的数据
# import geopandas

import matplotlib.pyplot as plt
# import matplotlib.ticker as mticker
# from matplotlib.pyplot import savefig
from cycler import cycler

from global_variable import station_dic
import datetime



# %%
area = {"lat1": 24.875, "lat2": 40.125, "lon1": 69.875, "lon2": 105.125}
month = 'Jul'
# month = 'May'
gd = GetData(month)
rain = gd.get_rain_hourly()

time_flag = 'all'
rain = rain.resample(time='D').sum()   # 日降水时间序列


def get_station_mean(rain):
    """获得多个站点的平均值
    """
    pass
    rain_station_list = []
    for key in station_dic:  
        station = station_dic[key]
        # rain1 = tr.rain_station(station_dic_dic[key])  
        r = rain.sel(lat=station['lat'], lon=station['lon'], method='nearest')
        rain_station_list.append(r)
    rain_station_array = xr.concat(rain_station_list, pd.Index(station_dic.keys(), name='station'))
    # rain_station_array
    rain_24h_mean = rain_station_array.mean(dim='station')  # 多站的平均值
    return rain_24h_mean

rain_mean = get_station_mean(rain)

# %%
# rain_mean


def get_station_land(rain):
    """获得多个站点的平均值
    """
    pass
    # low = ['ShiQuanhe', 'MangYa', 'GeErmu', 'DuLan']
    low = ['ShiQuanhe', 'MangYa', 'GeErmu']
    # medium = ['GaiZe', 'TuoTuohe']
    medium = ['GaiZe']
    high  = ['TingRi', 'ShenZha', 'LaSa', 'NaQu', 'YuShu', 
              'DaRi', 'BaTang', 'LinZhi', 'ChangDu']

    landuse_list = [low, medium, high]
    rain_list = []

    for station_list in landuse_list:
        rain_station_list = []
        for key in station_list:  
            station = station_dic[key]
            r = rain.sel(lat=station['lat'], lon=station['lon'], method='nearest')
            rain_station_list.append(r)
        rain_station_array = xr.concat(rain_station_list, pd.Index(station_list, name='station'))
        rain_station_mean = rain_station_array.mean(dim='station')  # 多站的平均值
        rain_list.append(rain_station_mean)

    rain_land_dic = {}        
    rain_land_dic['low'] = rain_list[0]
    rain_land_dic['medium'] = rain_list[1]
    rain_land_dic['high'] = rain_list[2]

    return rain_land_dic
rain_land_dic = get_station_land(rain)
# %%
rain_land_dic['high']










## 这样就能做到数据处理和画图按顺序进行了吗
# %%

class Draw():
    
    def __init__(self, month):
        self.fontsize = 10
        self.month = month
        self.path = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/picture/'   # 这里要改

    def draw_time_sequence(self,ax, dic):

        QNSE = dic['QNSE']
        x_label = QNSE.coords['time']
        # x_label = str(x_label.dt.strftime('%d%H').values).split()
        x_label = x_label.dt.strftime('%d')


        # x_label = x_label.dt.strftime("%H")  # 转换时间维字符串格式
        # y = QNSE.values
        module_list = ['obs', 'ACM2','YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        # x = np.arange(len(x_label))
        # print(x)


        ccolor = ['black','red', 'red', 'blue', 'blue', 'green', ]
        # ccolor = ['black','red', 'cyan', 'green', 'blue', 'orange', ]
        lline_style = ['-', '-', '--', '-', '--', '-']
        # mmarker = ['o', '^', '^', '*', '*', '+']
        mmarker = ['o', 'o', 'o', 'o', 'o', 'o']
        custom_cycler = (
            cycler(color=ccolor) +
            cycler(linestyle=lline_style))            
            # # cycler(label=module_list)
            # + cycler(marker=mmarker)
            #             )
        
        j = 0
        ax.set_prop_cycle(custom_cycler)
        for i in module_list:
            # y = dr.loc[i,:].values
            y = dic[i].values
            y = np.around(y,2)
            # ax.plot(x_label, y, label=i, color=ccolor[j])
            # ax.plot(x_label, y, label=i, lw=4)
            if i == 'obs':
                i = 'OBS'
            ax.plot(x_label, y, label=i, lw=2, markersize=5)
            j +=1 

        # ax.set_xticks(x_label[::24])  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax.set_xticks(x_label[::2])  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        # ax.xaxis.set_tick_params(labelsize=15)
        ax.xaxis.set_tick_params(labelsize=self.fontsize*1.8, rotation=45)
        ax.tick_params(which='major',length=8,width=1.0) # 控制标签大小 
        ax.tick_params(which='minor',length=4,width=0.5)  #,colors='b')
        ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
        # ax.set_yticks(np.arange(0, 5.01, 0.1))  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax.yaxis.set_tick_params(labelsize=self.fontsize*1.8)
        ax.set_ylim(0,50)
        # ax.set_ylim(0.0,0.6)
        ax.set_xlabel("Time(UTC)", fontsize=self.fontsize*2.0)
        ax.set_ylabel("Precipitation (mm)", fontsize=self.fontsize*2.0)
        # ax.legend()
        # ax.set_title("201607", fontsize=18)
        # fig.savefig('/home/fengxiang/Project/Asses_pbl_July/Draw/Rain/time_sequecnce.png')
    
    def combine_fig(self, rain):
        """[summary]

        Args:
            rain ([type]): 所有模式的日降水
        """

        # fig = plt.figure(figsize=(18, 15), dpi=200)  # 创建页面
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

        # axes = [None] * 9  # 设置一个维度为8的空列表
        num = 15
        axes = [None] * num
        # axes = [None] * 14  # 设置一个维度为8的空列表
        for i in range(num):
            axes[i] = fig.add_subplot(grid[i])
        
        i = 0
        # station_dic_dic = station_dic  # 这里要改
        dic = {}

        for key in station_dic:  
            station = station_dic[key]
            # rain1 = tr.rain_station(station_dic_dic[key])  
            r = rain.sel(lat=station['lat'], lon=station['lon'], method='nearest')
            dic[key] = r
            axes[i].set_title(key, fontsize=self.fontsize*2.0, loc='left', y=0.88, x=0.05)
            i += 1

        for i,j in zip(range(num),dic):
            self.draw_time_sequence(axes[i], dic[j])
            # axes[i].set_ylim(0.0, 40.0)
            # axes[i].set_yticks(np.arange(0, 40.1, 5.0))
            axes[i].set_ylim(0.0, 50.0)
            axes[i].set_yticks(np.arange(0, 50.1, 5.0))
            
        axes[13].legend(ncol=3 ,bbox_to_anchor=(0.5,-0.45) ,loc='lower center',fontsize=self.fontsize*2.0, edgecolor='white')
        # flnm = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/rain_staion1_'+self.month+'.png'   # 这里要改
        flnm = self.path+'rain_staion_'+self.month+'.png'   # 这里要改
        fig.savefig(flnm)

    def draw_single(self, rain, title):
        """[summary]

        Args:
            rain ([DataArray]): 一个模式的降水
        """
        fig = plt.figure(figsize=(12, 9), dpi=200)  # 创建页面
        ax = fig.add_axes([0.12, 0.25, 0.8, 0.7])

        self.draw_time_sequence(ax, rain)
        ax.set_yticks(np.arange(0, 51, 5))
            
        ax.legend(ncol=3 ,bbox_to_anchor=(0.5,-0.31) ,loc='lower center',fontsize=self.fontsize*2.0, edgecolor='white')
        ax.set_title(title, fontsize=22)
        # # flnm = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/rain_staion1_'+self.month+'.png'   # 这里要改
        flnm = self.path+'rain_staion_mean_'+self.month+"_"+title+'.png'   # 这里要改
        # fig.suptitle('The mean time sequence', fontsize=self.fontsize*2.0)
        fig.savefig(flnm)


Dr = Draw(month)
# Dr.combine_fig(rain)


title_dic = {'low':'BareLand', 'medium':'Bush', 'high':'GrassLand'}
Dr.draw_single(rain_mean, 'Mean')
for i in title_dic:
    Dr.draw_single(rain_land_dic[i], title_dic[i])
# plt.show()




# time_index_12 = da_obs.time.sel(time=datetime.time(int('12')))  
# time_index_00 = da_obs.time.sel(time=datetime.time(int('00')))  
# time_index = np.union1d(time_index_00.values, time_index_12.values)