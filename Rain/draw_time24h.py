#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
逐日降水的时间变化曲线
1. 这里的时间是从当天的08时到第二天的08时(BJT), 也就是世界时(00-23)

2. 24小时降水变化曲线(平均), 各站点分开
3. 多站点求平均，上面两张图应该都要有
站点降水的时间序列
-----------------------------------------
Time             :2021/06/04 14:32:20
Author          :Forxd
Version          :1.0
'''

# %%
import xarray as xr
from read_data import TransferData, GetData
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from cycler import cycler

from global_variable import station_dic
import datetime



# %%
area = {"lat1": 24.875, "lat2": 40.125, "lon1": 69.875, "lon2": 105.125}
month = 'Jul'
# month = 'May'
gd = GetData(month)
rain = gd.get_rain_hourly()

def get_rain24h(rain):
    rain_list = []
    ttime_index = pd.date_range(start='20160701 00', end='20160701 23', freq='H')
    for i in range(24):
        time_index = rain.time.sel(time=datetime.time(int(i)))  
        rain_hour = rain.sel(time=time_index).mean(dim='time')
        rain_list.append(rain_hour)
    # rain_24h = xr.concat(rain_list, pd.Index(np.arange(24), name='time'))
    rain_24h = xr.concat(rain_list, pd.Index(ttime_index, name='time'))
    return rain_24h

rain = get_rain24h(rain)  # 各测站单独的降水

# %%

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

rain_24h_mean = get_station_mean(rain)
    

# cc
# %%
# rain_hour
rain_24h_mean

## 这样就能做到数据处理和画图按顺序进行了吗
# %%

class Draw():
    
    def __init__(self, month):
        self.fontsize = 10
        self.month = month
        self.path = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/picture/'   # 这里要改
        self.module_list = ['obs', 'ACM2','YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']

    def draw_time_sequence(self,ax, dic):

        QNSE = dic['QNSE']
        x_label = QNSE.coords['time']
        # x_label = str(x_label.dt.strftime('%d%H').values).split()
        x_label = x_label.dt.strftime('%H')


        # x_label = x_label.dt.strftime("%H")  # 转换时间维字符串格式
        # y = QNSE.values
        # module_list = ['obs', 'ACM2','YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        # x = np.arange(len(x_label))
        # print(x)


        ccolor = ['black','red', 'red', 'blue', 'blue', 'green', ]
        # ccolor = ['black','red', 'cyan', 'green', 'blue', 'orange', ]
        lline_style = ['-', '-', '--', '-', '--', '-.']
        mmarker = ['o', '^', '^', '*', '*', '+']
        custom_cycler = (
            cycler(color=ccolor) +
            cycler(linestyle=lline_style) +           
            # cycler(label=module_list)
            cycler(marker=mmarker)
                        )
        
        j = 0
        ax.set_prop_cycle(custom_cycler)
        for i in self.module_list:
            # y = dr.loc[i,:].values
            y = dic[i].values
            y = np.around(y,2)
            # ax.plot(x_label, y, label=i, color=ccolor[j])
            # ax.plot(x_label, y, label=i, lw=4)
            if i == 'obs':
                i = 'OBS'
            ax.plot(x_label, y, label=i, lw=2)
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
        # ax.set_ylim(0,5.1)
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

        fig = plt.figure(figsize=(18, 15), dpi=200)  # 创建页面
        grid = plt.GridSpec(3,
                            3,
                            figure=fig,
                            left=0.05,
                            right=0.98,
                            bottom=0.12,
                            top=0.97,
                            wspace=0.2,
                            # hspace=0.25)
                            hspace=0.3)

        axes = [None] * 9  # 设置一个维度为8的空列表
        axes[0] = fig.add_subplot(grid[0])
        axes[1] = fig.add_subplot(grid[1])
        axes[2] = fig.add_subplot(grid[2])
        axes[3] = fig.add_subplot(grid[3])
        axes[4] = fig.add_subplot(grid[4])
        axes[5] = fig.add_subplot(grid[5])
        axes[6] = fig.add_subplot(grid[6])
        axes[7] = fig.add_subplot(grid[7])
        axes[8] = fig.add_subplot(grid[8])


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

        for i,j in zip(range(9),dic):
            self.draw_time_sequence(axes[i], dic[j])
            axes[i].set_ylim(0.0, 1)
            # axes[i].set_yticks(np.arange(0, 40.1, 5.0))
            axes[i].set_yticks(np.arange(0, 1.1, 0.1))
            
        axes[7].legend(ncol=3 ,bbox_to_anchor=(0.5,-0.55) ,loc='lower center',fontsize=self.fontsize*2.0, edgecolor='white')
        # flnm = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/rain_staion1_'+self.month+'.png'   # 这里要改
        flnm = self.path+'rain_staion_24h'+self.month+'.png'   # 这里要改
        fig.savefig(flnm)

    def draw_single(self, rain):
        """[summary]

        Args:
            rain ([DataArray]): 一个模式的降水
        """
        fig = plt.figure(figsize=(12, 9), dpi=200)  # 创建页面
        ax = fig.add_axes([0.12, 0.25, 0.8, 0.7])

        self.draw_time_sequence(ax, rain)
        ax.set_yticks(np.arange(0, 1.1, 0.1))
            
        ax.legend(ncol=3 ,bbox_to_anchor=(0.5,-0.31) ,loc='lower center',fontsize=self.fontsize*2.0, edgecolor='white')
        # # flnm = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/rain_staion1_'+self.month+'.png'   # 这里要改
        flnm = self.path+'rain_staion_24h_station_mean'+self.month+'.png'   # 这里要改
        # fig.suptitle('The mean time sequence', fontsize=self.fontsize*2.0)
        fig.savefig(flnm)



Dr = Draw(month)
# Dr.combine_fig(rain)
Dr.draw_single(rain_24h_mean)
plt.show()




# time_index_12 = da_obs.time.sel(time=datetime.time(int('12')))  
# time_index_00 = da_obs.time.sel(time=datetime.time(int('00')))  
# time_index = np.union1d(time_index_00.values, time_index_12.values)