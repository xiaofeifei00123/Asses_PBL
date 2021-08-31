#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
日降水变化曲线, 每天的降水量连成的曲线
画时间序列图, 各区域降水量24小时的变化
逐日降水的时间变化曲线
这里的时间是从前一天的


选择性的在于:
    降水时间序列
    1. 不同站点   
    2. 不同下垫面站点的均值
    3. 不同下垫面区域的均值

流程:
    获取数据
    处理数据
    画图

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
from data_process_landause import get_landmask
import datetime

# %%
#######################################################
## 获取数据+预处理(数据筛选)
area = {"lat1": 24.875, "lat2": 40.125, "lon1": 69.875, "lon2": 105.125}
month = 'Jul'
time_flag = 'all'
gd = GetData(month)
rain = gd.get_rain_hourly()
#######################################################
def get_rain_day(rain):
    """获取数据
    """
    rain_day = rain.resample(time='D').sum()   # 日降水时间序列
    return rain_day
def get_rain24h(rain):
    """每小时的(每天这个小时的平均值)降水
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

rain_day = get_rain_day(rain)
rain_24h = get_rain24h(rain)  # 所有格点的吧

# %%
#######################################################
## 处理数据
class Data_process():
    """处理数据
    返回值统一控制成Dataset
    """
    def __init__(self, rain_day) -> None:
        pass
        self.rain_day = rain_day  # 逐日的降水


    def get_station(self,):
        """ 1. 不同站点
        """
        pass
        rain_station_list = []
        for key in station_dic:  
            station = station_dic[key]
            # rain1 = tr.rain_station(station_dic_dic[key])  
            r = self.rain_day.sel(lat=station['lat'], lon=station['lon'], method='nearest')
            rain_station_list.append(r)
        rain_station = xr.concat(rain_station_list, pd.Index(station_dic.keys(), name='station'))
        # rain_station_array
        # rain_24h_mean = rain_station_array.mean(dim='station')  # 多站的平均值
        # return rain_24h_mean
        return rain_station

    def get_station_mean(self,):
        """2. 不同站点均值
        按照下垫面类型划分
        """
        pass
        land_dic = {}

        ## 过渡区
        land_dic['transition'] = ['GaiZe', 'TuoTuohe', 'ShenZha', 'DuLan', 'TingRi']
        land_dic['dry'] = ['ShiQuanhe', 'MangYa', 'GeErmu']
        land_dic['wet'] = ['ShenZha','LaSa', 'NaQu', 'YuShu', 'DaRi', 'BaTang', 'ChangDu', 'LinZhi']
        land_dic['all'] = land_dic['transition']+land_dic['dry']+land_dic['wet']

        ## 下垫面类型区域
        land_dic['bare'] = ['ShiQuanhe', 'MangYa', 'GeErmu']
        land_dic['bush'] = ['GaiZe']
        land_dic['grass'] = ['ShenZha','LaSa', 'NaQu', 'YuShu', 'DaRi', 'BaTang', 'ChangDu']

        rain_list = []

        for key in land_dic:
            station_list = land_dic[key]
            rain_station_list = []
            for i in station_list:  # 循环出站点
                station = station_dic[i]
                r = self.rain_day.sel(lat=station['lat'], lon=station['lon'], method='nearest')
                rain_station_list.append(r)
            rain_station_array = xr.concat(rain_station_list, pd.Index(station_list, name='station'))
            rain_station_mean = rain_station_array.mean(dim='station')  # 多站的平均值
            rain_list.append(rain_station_mean)
        rain_land = xr.concat(rain_list, pd.Index(land_dic.keys(), name='landtype'))
        return rain_land

    def get_rain_land_area(self,):
        """ 3. 不同下垫面区域均值
        """
        land_mask = get_landmask()
        rain_land_area_list = []

        land_list = ['grass', 'bare', 'bush']  # 主要就是这三个区域，区域也没有分界线的说法
        for i in land_list:
            rain_land_grid = self.rain_day*land_mask[i]
            rain_land_area_list.append(rain_land_grid.mean(dim=['lat','lon']))
            # rain_land_area[i] = rain_land_grid.mean(dim=['lat','lon'])
        # rain_land_area['mean'] = rain_24h.mean(dim=['lat', 'lon'])
        rain_land_area = xr.concat(rain_land_area_list,pd.Index(land_list, name='landtype') )
            # rain_station_array = xr.concat(rain_station_list, pd.Index(station_list, name='station'))
        rain_land_area_mean = rain_land_area.mean(dim='landtype')
        rain_land_area_mean.coords['landtype'] = 'all'
        rain_land_area_return = xr.concat([rain_land_area, rain_land_area_mean], dim='landtype')
        return rain_land_area_return

###################  结束  ##########################





# %%
class Draw():
    
    def __init__(self, month):
        self.fontsize = 10
        self.month = month
        self.path = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/picture/'   # 这里要改

    def draw_time_sequence(self,ax, ds, title_dic):

        QNSE = ds['QNSE']
        x_label_origin = QNSE.coords['time']
        # x_label = str(x_label.dt.strftime('%d%H').values).split()
        if title_dic['time_type'] == '31':
            x_label = x_label_origin.dt.strftime('%d')
            ax.set_ylim(0,50)
            ax.set_yticks(np.arange(0, 51, 5))  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        elif title_dic['time_type'] == '24':
            x_label = x_label_origin.dt.strftime('%H')
            ax.set_ylim(0,1)
            ax.set_yticks(np.arange(0, 1.1, 0.1))  # 这个是选择哪几个坐标画上来的了,都有只是显不显示


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
            y = ds[i].values
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
        ax.set_xlabel("Time(UTC)", fontsize=self.fontsize*2.0)
        ax.set_ylabel("Precipitation (mm)", fontsize=self.fontsize*2.0)
    
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

    def draw_single(self, rain, title_dic):
        """[summary]

        Args:
            rain ([DataArray]): 一个模式的降水
        """
        fig = plt.figure(figsize=(12, 9), dpi=200)  # 创建页面
        ax = fig.add_axes([0.12, 0.25, 0.8, 0.7])
        # x_label = rain['QNSE'].coords['time']

        # if title[-2:] == '31':
        #     ax.set_yticks(np.arange(0, 51, 5))
        #     x_label = x_label.dt.strftime('%d')
        #     ax.set_ylim(0,50)
        # elif title[-2:] == '24':
        #     ax.set_yticks(np.arange(0, 1, 0.1))
        #     x_label = x_label.dt.strftime('%H')
        #     ax.set_ylim(0,1)
        # else:
        #     print("画图时时间出了问题")
            
        self.draw_time_sequence(ax, rain, title_dic)
        # ax.set_xticks(x_label[::2])  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax.legend(ncol=3 ,bbox_to_anchor=(0.5,-0.31) ,loc='lower center',fontsize=self.fontsize*2.0, edgecolor='white')
        ax.set_title(title_dic['landtype'], fontsize=22)
        # # flnm = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/rain_staion1_'+self.month+'.png'   # 这里要改
        # flnm = self.path+self.month+"_"+title+'.png'   # 这里要改
        flnm = self.path + title_dic['area_type']+ "_"+title_dic['time_type']+ "_"+title_dic['landtype']+".png"
        # fig.suptitle('The mean time sequence', fontsize=self.fontsize*2.0)
        fig.savefig(flnm)


if __name__ == '__main__':
    pass
    # rain_draw = rain_day
    for rain_draw in [rain_day, rain_24h]:
        time_type = str(len(rain_draw.time))

        dap = Data_process(rain_draw)
        # d1 = dap.get_station()
        d2 = dap.get_station_mean()
        d3 = dap.get_rain_land_area()
        
        month = 'Jul'
        Dr = Draw(month)

        # rain_land_area = get_rain_land_area()
        # land_list = ['grass', 'bare', 'bush']
        land_list = ['transition', 'dry', 'wet', 'all']
        for i in land_list:
            pass
            title_dic = {
                'landtype':i,
                'time_type':str(time_type),
                'area_type':'station'
            }
            Dr.draw_single(d2.sel(landtype=i), title_dic)
            
        land_list = ['bare', 'bush', 'grass', 'all']
        for i in land_list:
            pass
            title_dic = {
                'landtype':i,
                'time_type':str(time_type),
                'area_type':'space'
            }
            Dr.draw_single(d3.sel(landtype=i), title_dic)
            Dr.draw_single(d2.sel(landtype=i), title_dic)
        
        

# time_index_12 = da_obs.time.sel(time=datetime.time(int('12')))  
# time_index_00 = da_obs.time.sel(time=datetime.time(int('00')))  
# time_index = np.union1d(time_index_00.values, time_index_12.values)