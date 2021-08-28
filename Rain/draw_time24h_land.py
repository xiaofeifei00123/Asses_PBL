#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
逐日降水的时间变化曲线
1. 这里的时间是从当天的08时到第二天的08时(BJT), 也就是世界时(00-23)
2. 24小时降水变化曲线(平均), 各站点分开
3. 多站点求平均，上面两张图应该都要有
站点降水的时间序列
4. 不同下垫面类型区域的平均降水
-----------------------------------------
Time             :2021/06/04 14:32:20
Author          :Forxd
Version          :1.0
'''

# %%
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler

from read_data import TransferData, GetData
from data_process_landause import get_landmask
from global_variable import station_dic
import datetime



# %%
area = {"lat1": 24.875, "lat2": 40.125, "lon1": 69.875, "lon2": 105.125}
month = 'Jul'
# month = 'May'
gd = GetData(month)
rain = gd.get_rain_hourly()  # 原始的小时降水数据



def get_rain24h(rain):
    """[summary]

    Args:
        rain ([type]): [description]

    Returns:
        rain_24h: (lat:82, lon:142, time:24)
    """
    rain_list = []
    ttime_index = pd.date_range(start='20160701 00', end='20160701 23', freq='H')
    for i in range(24):
        time_index = rain.time.sel(time=datetime.time(int(i)))  
        rain_hour = rain.sel(time=time_index).mean(dim='time')
        rain_list.append(rain_hour)
    # rain_24h = xr.concat(rain_list, pd.Index(np.arange(24), name='time'))
    rain_24h = xr.concat(rain_list, pd.Index(ttime_index, name='time'))
    return rain_24h  

rain_24h = get_rain24h(rain)  # 所有格点的吧

# %%
# rain_24h['obs'][0].plot()
# rain

# %%

# %%
# rain
def get_rain_land_area():
    land_mask = get_landmask()
    rain_land_area_list = []

    land_list = ['grass', 'bare', 'bush']

    for i in land_list:
        # rain_land[i] = rain_24h*land_mask[i]
        rain_land_grid = rain_24h*land_mask[i]
        rain_land_area_list.append(rain_land_grid.mean(dim=['lat','lon']))
        # rain_land_area[i] = rain_land_grid.mean(dim=['lat','lon'])
    # rain_land_area['mean'] = rain_24h.mean(dim=['lat', 'lon'])
    rain_land_area = xr.concat(rain_land_area_list,pd.Index(land_list, name='landuse') )
        # rain_station_array = xr.concat(rain_station_list, pd.Index(station_list, name='station'))
    return rain_land_area

# %%
aa = get_rain_land_area()

# %%
# aa
import dask.array as da
da_ac = aa['ACM2']
# y = da.from_array(aa, chunks=100)
# y
# da_ac
# da_ac
y = da.from_array(da_ac, chunks=1)


# %%
aa.sel(landuse='grass')
# aa['ACM2']
# cc = xr.Dataset()
# cc['aa'] = 'aa'
# cc



# %%


# %%
# rain_24h

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
        lline_style = ['-', '-', '--', '-', '--', '-']
        # mmarker = ['o', '^', '^', '*', '*', '+']
        mmarker = ['o', 'o', 'o', 'o', 'o', 'o']
        custom_cycler = (
            cycler(color=ccolor) +
            cycler(linestyle=lline_style))           
            # cycler(label=module_list)
            # cycler(marker=mmarker)
                        # )
        
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

        # num = 15
        num = len(station_dic)
        axes = [None] * num # 设置一个num维度为8的空列表
        # axes = [None] * 14  
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
            axes[i].set_ylim(0.0, 1.5)
            axes[i].set_yticks(np.arange(0, 1.6, 0.2))
            if key == 'LinZhi':
                axes[i].set_ylim(0.0, 3)
                axes[i].set_yticks(np.arange(0, 4.1, 0.5))
            i += 1

        for i,j in zip(range(num),dic):
            self.draw_time_sequence(axes[i], dic[j])
            # axes[i].set_yticks(np.arange(0, 40.1, 5.0))
            # axes[i].set_yticks(np.arange(0, 1.1, 0.1))
            
        axes[13].legend(ncol=3 ,bbox_to_anchor=(0.5,-0.55) ,loc='lower center',fontsize=self.fontsize*2.0, edgecolor='white')
        # flnm = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/rain_staion1_'+self.month+'.png'   # 这里要改
        flnm = self.path+'rain_staion_24h'+self.month+'.png'   # 这里要改
        fig.savefig(flnm)

    def draw_single(self, rain, title):
        """[summary]

        Args:
            rain ([DataArray]): 一个模式的降水
        """
        fig = plt.figure(figsize=(12, 9), dpi=200)  # 创建页面
        ax = fig.add_axes([0.12, 0.25, 0.8, 0.7])

        self.draw_time_sequence(ax, rain)
        ax.set_yticks(np.arange(0, 1.1, 0.1))
            
        ax.legend(ncol=3 ,bbox_to_anchor=(0.5,-0.31) ,loc='lower center',fontsize=self.fontsize*2.0, edgecolor='white')
        ax.set_title(title, fontsize=22)
        # # flnm = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/rain_staion1_'+self.month+'.png'   # 这里要改
        flnm = self.path+'rain_staion_24h_station_mean'+self.month+"_"+title+'.png'   # 这里要改
        # fig.suptitle('The mean time sequence', fontsize=self.fontsize*2.0)
        fig.savefig(flnm)



Dr = Draw(month)
# Dr.combine_fig(rain)
# Dr.draw_single(rain_24h_mean)

# rain_land_dic = get_station_land(rain)  # 站点的数据

# title_dic = {'low':'BareLand', 'medium':'Bush', 'high':'GrassLand'}
# Dr.draw_single(rain_24h_mean, 'Mean')
# for i in title_dic:
#     Dr.draw_single(rain_land_dic[i], title_dic[i])

# plt.show()


rain_land_area = get_rain_land_area()
land_list = ['grass', 'bare', 'bush']

for i in land_list:
    pass
    Dr.draw_single(rain_land_area.sel(landuse=i), i+'land')
Dr.draw_single(rain_land_area.mean(dim='landuse'), 'mean'+'_land')
    







# time_index_12 = da_obs.time.sel(time=datetime.time(int('12')))  
# time_index_00 = da_obs.time.sel(time=datetime.time(int('00')))  
# time_index = np.union1d(time_index_00.values, time_index_12.values)