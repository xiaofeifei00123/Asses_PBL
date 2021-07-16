#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取模式和观测的TBB数据，并画图
-----------------------------------------
Time             :2021/06/02 12:04:38
Author           :Forxd
Version          :1.0
'''
import xarray as xr
import cmaps
import numpy as np
import xesmf as xe
import os
import wrf
import netCDF4 as nc
import pandas as pd

import salem
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.io.shapereader import Reader, natural_earth
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
import geopandas
import cmaps

from multiprocessing import Pool
from multiprocessing import Manager
import multiprocessing

import struct

class Get_data(object):
    pass
    def __init__(self, flag) -> None:
        super().__init__()
        self.flag = flag
        # self.var = 'tbb'

class Get_wrf(Get_data):

    def month_mean(self,flnm):
        """求月降水之和 """
        ds = xr.open_dataset(flnm)
        # r = ds[self.var]
        r = ds.ctt

        ## 判断某时次是白天还是夜间, 获得白天和夜间时次列表
        time_list = r.time.values
        day_list = []
        night_list = []
        for i in time_list:
            hour = str(i)[11:13]
            if int(hour) <= 13:  # 00时--13时时白天, 14--23时是夜间
                day_list.append(i)
            else:
                night_list.append(i)
        
        ## 根据白天还是夜间的需求，求总降水量
        if self.flag == 'all':
            rrr = r.mean(dim='time')
        elif self.flag == 'day':
            rain_day = r.sel(time=day_list)
            rrr = rain_day.mean(dim='time')
        else:
            rain_night = r.sel(time=night_list)
            rrr = rain_night.mean(dim='time')
        return rrr
    
    def get_total_rain(self,):

        model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','MYJ','YSU']
        # path = '/home/fengxiang/Project/Asses_pbl_July/Data/Process/'
        path = '/home/fengxiang/Project/Asses_pbl_July/Data/standard_data/'
        file_list = []
        for i in model_list:
            # flnm = path + "RAINNC_Jul_"+i+'_latlon'
            flnm = path + "ctt_Jul_"+i+'_latlon'
            file_list.append(flnm)

        rain = {}   # 各模式降水和
        # su = Summer_mean()
        # for fl in file_list[0:1]:
        for i in range(len(file_list)):
            # ds = xr.open_dataset(fl)
            fl = file_list[i]
            # print(fl)
            rain_total = self.month_mean(fl)
            rain[model_list[i]] = rain_total
        # print(rain['QNSE'])
        return rain  # 返回的是一个各试验tbb和的列表

class Get_obs(Get_data):

    def get_tbb_obs(self):
        pass
        file_path = '/home/fengxiang/Project/Asses_pbl_July/Data/standard_data/'
        # flnm = 'tbb_Jul_latlon'
        flnm = 'tbb_Jul_test'
        flnm = os.path.join(file_path,flnm)
        r = xr.open_dataarray(flnm)

        time_list = r.time.values
        day_list = []
        night_list = []
        for i in time_list:
            hour = str(i)[11:13]
            print(hour)
            if int(hour) <= 13:  # 00时--13时时白天, 14--23时是夜间
                day_list.append(i)
            else:
                night_list.append(i)
        
        ## 根据白天还是夜间的需求，求总降水量
        if self.flag == 'all':
            rrr = r.mean(dim='time')
        elif self.flag == 'day':
            rain_day = r.sel(time=day_list)
            rrr = rain_day.mean(dim='time')
        else:
            rain_night = r.sel(time=night_list)
            rrr = rain_night.mean(dim='time')
        return rrr


class Draw(Get_data):
    def draw_station(self,ax):
        """在地图上标注站点

        Args:
            ax ([type]): [description]
        """
        station = {
            # 'TR': {
            #     'lat': 28.6,
            #     'lon': 87.0
            # },
            'NQ': {
                'lat': 31.4,
                'lon': 92.0
            },
            'LS': {
                'lat': 29.6,
                'lon': 91.1
            },
            'TTH': {
                'lat': 34.2,
                'lon': 92.4
            },
            'GZ': {
                'lat': 32.3,
                'lon': 84.0
            },
            'SZ': {
                'lat': 30.9,
                'lon': 88.7
            },
            'SQH': {
                'lat': 32.4,
                'lon': 80.1
            },
            # 'JinChuan': {
            #     'lat': 31.29,
            #     'lon': 102.04
            # },
            # 'JinLong': {
            #     'lat': 29.00,
            #     'lon': 101.50
            # },
        }
        values = station.values()
        station_name = list(station.keys())
        print(type(station_name[0]))
        # print(station_name[0])
        x = []
        y = []
        for i in values:
            y.append(float(i['lat']))
            x.append(float(i['lon']))

        ## 标记出站点
        ax.scatter(x,
                y,
                color='black',
                transform=ccrs.PlateCarree(),
                alpha=1.,
                linewidth=0.2,
                s=10)
        ## 给站点加注释
        for i in range(len(x)):
            print(x[i])
            ax.text(x[i] - 1,
                    y[i] + 0.5,
                    station_name[i],
                    transform=ccrs.PlateCarree(),
                    alpha=1.,
                    fontdict={
                        'size': 9,
                    })


    def create_map(self,ax):
        """创建地图对象
        ax 需要添加底图的画图对象

        Returns:
            ax: 添加完底图信息的坐标子图对象
        """
        proj = ccrs.PlateCarree()
        # --设置地图属性
        # 画省界
        provinces = cfeat.ShapelyFeature(Reader(
            '/mnt/Disk4T_5/fengxiang_file/Data/Map/cn_shp/Province_9/Province_9.shp'
        ).geometries(),
                                        proj,
                                        edgecolor='k',
                                        facecolor='none')

        # 画青藏高原
        Tibet = cfeat.ShapelyFeature(
            Reader('/home/fengxiang/Data/shp_tp/Tibet.shp').geometries(),
            proj,
            edgecolor='k',
            facecolor='none')

        # ax.add_feature(provinces, linewidth=0.6, zorder=2)
        ax.add_feature(Tibet, linewidth=0.6, zorder=2)  # 添加青藏高原区域

        # --设置图像刻度
        # ax.set_xticks(np.arange(80, 102 + 2, 4))
        # ax.set_yticks(np.arange(26, 38 + 2, 2))
        ax.set_xticks(np.arange(70, 105 + 2, 4))
        ax.set_yticks(np.arange(25, 45 + 2, 2))
        ax.xaxis.set_major_formatter(LongitudeFormatter())
        ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
        ax.yaxis.set_major_formatter(LatitudeFormatter())
        ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
        ax.tick_params(axis='both', labelsize=5, direction='out')
        # -- 设置图像范围
        # ax.set_extent([78, 98, 26, 38], crs=ccrs.PlateCarree())
        return ax


    def draw_contourf(self, data, ax, cmap, title):
        """根据具体的数据，画子图

        Args:
            data ([type]): [description]
            ax ([type]): [description]
            cmap ([type]): [description]
            title ([type]): [description]

        Returns:
            [type]: [description]
        """
        # levels = np.arange(0, 66, 5)  # 设置colorbar分层
        # levels = [200, 205, 210, 215, 220, 225, 235, 245, 250]  # 需要画出的等值线
        # levels = [200, 205, 210, 215, 220, 225, 235, 245, 250, 260]  # 需要画出的等值线
        # levels = [220, 225, 235, 245, 250, 260, 270, 280, 290, 300]  # 需要画出的等值线
        # levels = [240, 250, 260, 270, 280, 290, 300,310]  # 需要画出的等值线
        levels = [220,230, 240, 250, 260, 270, 280, 290, 300]  # 需要画出的等值线
        ax = self.create_map(ax)
        print(title[0])

        x = data.lon
        y = data.lat
        crx = ax.contourf(x,
                        y,
                        data,
                        cmap=cmap,
                        extend='both',
                        levels=levels,
                        transform=ccrs.PlateCarree())
        ## 画时间标签
        # ax.set_title(title[0], loc='left', y=0.82, fontsize=12)
        ax.set_title(title[0], loc='left', fontsize=12)
        # ax.set_extent([78, 102, 26, 38], crs=ccrs.PlateCarree())
        ax.set_extent([70, 105, 25, 45], crs=ccrs.PlateCarree())
        ax.xaxis.set_tick_params(labelsize=10)
        ax.yaxis.set_tick_params(labelsize=10)
        ## 画(a), (b)..
        # ax.text(99.5, 37, title[1])
        # ax.text(100, 38, title[1])
        ax.set_title(title[1], loc='right', fontsize=12)
        self.draw_station(ax)
        return crx


    def draw_obs(self, ctt):
        """画观测的TBB云顶亮温，主程序
        Args:
            year ([type]): 
        """
        # flnm_2005 = "./Data/TBB_2005.nc"
        # flnm_2014 = "./Data/TBB_2014.nc"

        # # ftime_2005 = get_time_list('2005')
        # # ftime_2014 = get_time_list('2014')
        # ds_2005 = xr.open_dataset(flnm_2005)
        # ds_2014 = xr.open_dataset(flnm_2014)

        ## 获取地形文件
        shp_file = '/home/fengxiang/Data/shp_tp/Tibet.shp'
        shp = geopandas.read_file(shp_file)
        # print(ctt['QNSE'])

        # time_list_2005 = ['0728_1200', '0728_2000', '0729_0000']
        # time_list_2014 = ['0819_0600', '0819_1200', '0819_1900']

        # ## 生成每张图上的标签
        # tbb_2005 = {}
        # tbb_2014 = {}
        # #     dic_title = {}
        # # title_list = [None]*6
        # for i in time_list_2005:
        #     str_title = i[-4:-2] + "00UTC" + " " + i[2:4] + "Jul 2005"
        #     tbb_val = ds_2005[i]
        #     tbb_2005[i] = {'title': str_title, 'val': tbb_val}
        # for i in time_list_2014:
        #     str_title = i[-4:-2] + "00UTC" + " " + i[2:4] + "Aug 2014"
        #     tbb_val = ds_2014[i]
        #     tbb_2014[i] = {'title': str_title, 'val': tbb_val}

        ## --->画图
        proj = ccrs.PlateCarree()  # 创建坐标系
        fig = plt.figure(figsize=(8, 7), dpi=400)  # 创建页面

        ## 设置子图位置和大小
        grid = plt.GridSpec(3,
                            2,
                            figure=fig,
                            left=0.07,
                            right=0.96,
                            bottom=0.12,
                            top=0.96,
                            wspace=0.3,
                            hspace=0.3)
        axes = [None] * 6  # 设置一个维度为5的空列表
        # # print(axes)
        axes[0] = fig.add_subplot(grid[0, 0:1], projection=proj)
        axes[1] = fig.add_subplot(grid[0, 1:2], projection=proj)
        axes[2] = fig.add_subplot(grid[1, 0:1], projection=proj)
        axes[3] = fig.add_subplot(grid[1, 1:2], projection=proj)
        axes[4] = fig.add_subplot(grid[2, 0:1], projection=proj)
        axes[5] = fig.add_subplot(grid[2, 1:2], projection=proj)
        # # 设置colorbar位置和名称等
        ccc = cmaps.precip3_16lev_r
        colors = mpl.cm.get_cmap(ccc)
        col = colors(np.linspace(0, 1, 18))
        cccc = mpl.colors.ListedColormap([
            col[0],
            col[1],
            col[2],
            (231 / 250, 177 / 250, 22 / 250),
            col[4],
            col[6],
            '#85f485',
            '#16c516',
            'white',
        ])

        cmap = cccc
        label_list_2005 = ['(a)', '(b)', '(c)']

        # ctt = get_total_rain()

        shp_file = '/home/fengxiang/Data/shp_tp/Tibet.shp'
        shp = geopandas.read_file(shp_file)
        ## 画子图
        self.draw_contourf(ctt['QNSE'].salem.roi(shape=shp), axes[2], cmap, ['QNSE', '(a)'])
        self.draw_contourf(ctt['QNSE_EDMF'].salem.roi(shape=shp), axes[3], cmap, ['QNSE_EDMF', '(b)'])
        self.draw_contourf(ctt['TEMF'].salem.roi(shape=shp), axes[4], cmap, ['TEMF', '(c)'])
        self.draw_contourf(ctt['YSU'].salem.roi(shape=shp), axes[1], cmap, ['YSU', '(d)'])
        cf = self.draw_contourf(ctt['MYJ'].salem.roi(shape=shp), axes[0], cmap, ['MYJ', '(e)'])
        cf = self.draw_contourf(ctt['obs'].salem.roi(shape=shp), axes[5], cmap, ['obs', '(f)'])
        # cf = draw_contourf(tbb_2014[time_list_2014[2]]['val'], axes[5], cmap, [tbb_2014[time_list_2014[2]]['title'], '(f)'])


        # ## 画色标
        ax6 = fig.add_axes([0.18, 0.04, 0.7, 0.03])  # 重新生成一个新的坐标图
        # bounds = [205, 210, 215, 220, 225, 235, 245]

        bounds = [220,230, 240,250, 260,270,280,290]  # 需要画出的等值线
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')
        cb = fig.colorbar(
            cf,
            cax=ax6,
            orientation='horizontal',
            ticks=bounds,
            fraction=0.05,  # 色标大小
            pad=0.1,  # 透明度
            # extend='both'
        )

        fig_name = "TBB"
        fig_name = 'TBB_'+self.flag
        fig.savefig(fig_name)
if __name__ == '__main__':
    pass
    # flag = 'all'
    flag = 'day'
    # flag = 'night'
    
    # gd = Get_data(flag)
    gw = Get_wrf(flag)
    ctt = gw.get_total_rain()
    get_obs = Get_obs(flag)
    ctt_obs = get_obs.get_tbb_obs()
    ctt['obs'] = ctt_obs
    print(ctt_obs)
    dr = Draw(flag)
    dr.draw_obs(ctt)



    ## Test

    # file_list = get_obs.get_file_list()
    # flnm = file_list[0]
    # # t = pd.to_datetime('20160701_1015', format='%Y%m%d_%H%M')
    # cc = get_obs.get_tbb_obs(flnm)
    # # aa = get_obs.regrid(flnm)
    # ctt['obs'] = aa


    # # ff = go.get_dataset(file_list)
    # # # print(flnm)
    # # cc = go.get_tbb_obs(flnm)
    # # print(cc)
    # # # print(ctt)


    # # # print(ctt)


    # # # # print(cc)
    # # # # print(cc.lat.values)