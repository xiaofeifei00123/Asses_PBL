#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
画降水空间分布图
白天、夜间、所有时次的总降水分布
-----------------------------------------
Time             :2021/06/15 17:21:14
Author           :Forxd
Version          :1.0
'''
# %%

import sys,os
BASE_PATH = os.path.dirname(os.path.abspath(__file__))
sys.path.append(BASE_PATH)
# from pandas.core import frame
# from shapely.ops import transform
# import xarray as xr
# import pandas as pd
import numpy as np
# import xesmf as xe
import os
# import netCDF4 as nc

# from datetime import datetime
# import datetime as Datetime
# 
import salem
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.io.shapereader import Reader, natural_earth
import matplotlib as mpl
from matplotlib.path import Path
# import matplotlib.patches as patches
import matplotlib.pyplot as plt
# import matplotlib.ticker as mticker
# import matplotlib.gridspec as gridspec
import geopandas
# import cmaps
# import meteva.base as mb
from read_data import TransferData
from read_data import GetData
from get_cmap import get_cmap_rain2
from global_variable import station_dic


# %%
######################################
## 测试单个图的
######################################
# area = {"lat1": 24.875, "lat2": 40.125, "lon1": 69.875, "lon2": 105.125}
# month = 'Jul'
# time_flag = 'all'
# gd = GetData(month)
# rain = gd.get_rain_hourly()
# tr = TransferData(ds=rain, area=area, time_flag=time_flag)
# rain = tr.rain_time_total()


# %%

class Draw(object):

    def __init__(self, fig_name, flag, levels) -> None:
        self.fig_name = fig_name
        self.flag = flag
        self.levels = levels
        ## 所有路径推荐使用绝对路径
        self.path_province = '/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/Province_9/Province_9.shp'
        self.path_tibet = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_tp/Tibet.shp'
        self.picture_path = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/picture'


        # self.levels = [10, 25,  50, 75, 100, 200,
        #                250, 300, 400, 600, 800, 1000,]
        # self.levels = [1, 5, 10, 25, 50, 100, 200,
        #                250, 300, 400, 600, 800, 1000,]
        # self.levels = [0.1, 0.5, 1, 2.5, 5.0, 10.0, 20.0,
        #                25.0, 30.0, 40.0, 60.0, 80.0, 100.0,]
    

    def draw_station(self, ax):
        # station = {
        #     'NQ': {'lat': 31.4, 'lon': 92.0},
        #     'LS': {'lat': 29.6, 'lon': 91.1},
        #     'TTH': {'lat': 34.2, 'lon': 92.4},
        #     'GZ': {'lat': 32.3, 'lon': 84.0},
        #     'SZ': {'lat': 30.9, 'lon': 88.7},
        #     'SQH': {'lat': 32.4, 'lon': 80.1},
        # }
        station = station_dic
        values = station.values()
        # station_name = list(station.keys())
        station_name = []
        x = []
        y = []
        for i in values:
            y.append(float(i['lat']))
            x.append(float(i['lon']))
            station_name.append(i['abbreviation'])

        # 标记出站点
        ax.scatter(x,
                   y,
                   color='black',
                   transform=ccrs.PlateCarree(),
                   alpha=1.,
                   linewidth=0.2,
                   s=10)
        # 给站点加注释
        for i in range(len(x)):
            # print(x[i])
            ax.text(x[i] - 1,
                    y[i] + 0.5,
                    station_name[i],
                    transform=ccrs.PlateCarree(),
                    alpha=1.,
                    fontdict={
                        'size': 9,
            })

    def create_map(self, ax):
        """创建地图对象
        ax 需要添加底图的画图对象

        Returns:
            ax: 添加完底图信息的坐标子图对象
        """
        proj = ccrs.PlateCarree()
        # --设置地图属性
        # 画省界
        provinces = cfeat.ShapelyFeature(
            Reader(self.path_province).geometries(),
            proj,
            edgecolor='k',
            facecolor='none')

        # 画青藏高原
        Tibet = cfeat.ShapelyFeature(
            Reader(self.path_tibet).geometries(),
            proj,
            edgecolor='k',
            facecolor='none')

        # ax.add_feature(provinces, linewidth=0.6, zorder=2)
        ax.add_feature(Tibet, linewidth=0.6, zorder=2)  # 添加青藏高原区域

        # --设置图像刻度
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

        norm = mpl.colors.BoundaryNorm(self.levels, cmap.N, extend='both')
        ax = self.create_map(ax)  # 创建坐标图像

        # if title[0] == 'obs':
        #     pass
        #     self.draw_patch(ax)
        # #     print("yes")

        x = data.lon
        y = data.lat
        crx = ax.contourf(x,
                          y,
                          data,
                          cmap=cmap,
                          norm=norm,
                          extend='both',
                          # extend='max',
                          #   levels=levels,
                          levels=self.levels,
                          transform=ccrs.PlateCarree())
        # ax.set_title(title[2], loc='left',  fontsize=12)
        ax.set_title(title[0], loc='left',  fontsize=12)
        ax.set_extent([70, 105, 25, 41], crs=ccrs.PlateCarree())
        ax.xaxis.set_tick_params(labelsize=10)
        ax.yaxis.set_tick_params(labelsize=10)
        # ax.text(99.5, 38.5, title[1])
        ax.set_title(title[1], loc='right',  fontsize=12)
        # ax.set_tilte(title[1], loc='right')
        # ax.text(78,26.2,title[0], fontsize=12)
        self.draw_station(ax)
        return crx

    def draw_patch(self, ax):
        """画矩形框
        """

        area = [None] * 3  # 设置一个维度为8的空列表
        area[0] = {"lat1": 33.5, "lat2": 40, "lon1": 80, "lon2": 90}  # north
        area[1] = {"lat1": 28, "lat2": 33,
                   "lon1": 83, "lon2": 94}  # south left
        area[2] = {"lat1": 26, "lat2": 33,
                   "lon1": 95, "lon2": 103}  # south right
        for i in range(3):
            lon = np.empty(4)
            lat = np.empty(4)
            lon[0], lat[0] = area[i]['lon1'], area[i]['lat1']
            lon[1], lat[1] = area[i]['lon2'], area[i]['lat1']
            lon[2], lat[2] = area[i]['lon2'], area[i]['lat2']
            lon[3], lat[3] = area[i]['lon1'], area[i]['lat2']
            x, y = lon, lat
            xy = list(zip(x, y))
            poly = plt.Polygon(xy, edgecolor="red", fc="none", lw=.9, alpha=1)
            ax.add_patch(poly)

    def draw_main(self, rain):

        shp = geopandas.read_file(self.path_tibet)
        model_list = ['QNSE', 'QNSE_EDMF', 'TEMF', 'ACM2', 'YSU']

        # --->画图
        proj = ccrs.PlateCarree()  # 创建坐标系
        fig = plt.figure(figsize=(9, 9), dpi=400)  # 创建页面
        grid = plt.GridSpec(3,
                            2,
                            figure=fig,
                            left=0.07,
                            right=0.96,
                            bottom=0.12,
                            top=0.96,
                            wspace=0.2,
                            hspace=0.1)

        axes = [None] * 6  # 设置一个维度为8的空列表
        axes[0] = fig.add_subplot(grid[0, 0:1], projection=proj)
        axes[1] = fig.add_subplot(grid[0, 1:2], projection=proj)
        axes[2] = fig.add_subplot(grid[1, 0:1], projection=proj)
        axes[3] = fig.add_subplot(grid[1, 1:2], projection=proj)
        axes[4] = fig.add_subplot(grid[2, 0:1], projection=proj)
        axes[5] = fig.add_subplot(grid[2, 1:2], projection=proj)

        cmap = get_cmap_rain2()

        time_2005 = '2800_2906 Jul 2005'
        time_2014 = '1900_2000 Aug 2014'

        self.draw_contourf(rain['ACM2'].salem.roi(
            shape=shp), axes[0], cmap, ['ACM2', '(a)', time_2005])
        self.draw_contourf(rain['ACM2'].salem.roi(
            shape=shp), axes[1], cmap, ['YSU', '(b)', time_2005])
        self.draw_contourf(rain['QNSE'].salem.roi(
            shape=shp), axes[2], cmap, ['QNSE', '(c)', time_2005])
        self.draw_contourf(rain['QNSE_EDMF'].salem.roi(
            shape=shp), axes[3], cmap, ['QNSE_EDMF', '(d)', time_2005])
        self.draw_contourf(rain['TEMF'].salem.roi(
            shape=shp), axes[4], cmap, ['TEMF', '(e)', time_2005])
        cf = self.draw_contourf(rain['obs'].salem.roi(
            shape=shp), axes[5], cmap, ['OBS', '(f)', time_2005])
        ax6 = fig.add_axes([0.18, 0.06, 0.7, 0.02])  # 重新生成一个新的坐标图

        cb = fig.colorbar(
            cf,
            cax=ax6,
            orientation='horizontal',
            ticks=self.levels,
            # fraction = 0.1,  # 色标大小
            pad=0.1,  #
        )
        path = os.path.join(self.picture_path, self.fig_name)
        fig.savefig(path)



######################################
## 测试一个图的
# levels = [10, 25,  50, 75, 100, 200, 250, 300, 400, 600, 800, 1000,]
# levels = [10, 25,  50, 75, 100, 125, 150, 175, 200, 250, 300, 500,]
# levels = [10, 50, 100, 150, 200, 225, 250, 275, 300,325, 350, 500]
# levels = [10, 40, 70, 100, 130, 160, 190, 220, 250, 300, 500,]

# # rain = rain-rain['obs']
# fig_name = 'rain_'+month + "_"+time_flag
# print("画  %s" %fig_name)
# dr = Draw(fig_name, time_flag, levels=levels)
# dr.draw_main(rain)
######################################




# %%

if __name__ == '__main__':
    def draw_hourly():
        """每个小时的降水画一张图
        """
        time_flag = 'all'
        area = {"lat1": 24.875, "lat2": 40.125, "lon1": 69.875, "lon2": 105.125}
        month = 'Jul'
        gd = GetData(month)
        rain = gd.get_rain_hourly()
        tr = TransferData(ds=rain, area=area, time_flag=time_flag)
        rain = tr.rain_hourly()
        fl_time = rain.time

        for i in fl_time:
            rain_hour = rain.sel(time=i)
            flnm_time = str(i.dt.strftime('%Y%m_%d%H').values)
            fig_name = 'rain_'+flnm_time
            dr = Draw(fig_name, time_flag)
            dr.draw_main(rain_hour)

    def draw_dual():
        # 一次画多个需要的图
        time_flag = 'night'
        area = {"lat1": 24.875, "lat2": 40.125, "lon1": 69.875, "lon2": 105.125}
        # levels = [10, 50, 100, 150, 200, 225, 250, 275, 300,325, 350, 500]
        levels = [10, 40, 70, 100, 130, 160, 190, 220, 250, 300, 500,]
        month_list = ['May', 'Jul']
        time_list = ['all', 'day', 'night']
        for month in month_list:
            for time_flag in time_list:

                gd = GetData(month)
                rain = gd.get_rain_hourly()

                # %%
                tr = TransferData(ds=rain, area=area, time_flag=time_flag)
                rain = tr.rain_time_total()
                # rain = rain-rain['obs']
                fig_name = 'rain_'+month + "_"+time_flag
                print("画  %s" %fig_name)
                dr = Draw(fig_name, time_flag, levels)
                dr.draw_main(rain)

    draw_dual()




## --------------------------------------------
''' ## 按照白天、夜间和月份画降水
    # time_flag = 'all'
    # time_flag = 'day'
    time_flag = 'night'
    area = {"lat1": 24.875, "lat2": 40.125, "lon1": 69.875, "lon2": 105.125}
    # month = 'Jul'
    month = 'May'

    # %%

    # %%
    gd = GetData(month)
    rain = gd.get_rain_hourly()

    # %%
    tr = TransferData(ds=rain, area=area, time_flag=time_flag)
    rain = tr.rain_time_total()
    # rain = rain-rain['obs']
    fig_name = 'rain_'+month + "_"+time_flag
    dr = Draw(fig_name, time_flag)
    dr.draw_main(rain) '''
## --------------------------------------------



## --------------------------------------------
    ## 测试逐小时降水画
    # rain = tr.rain_hourly()
    # rain = rain-rain['obs']
    # print(rain['obs'])
    # %%
    # fl_time = rain.time

    # for i in fl_time:
    #     rain_hour = rain.sel(time=i)
    #     flnm_time = str(i.dt.strftime('%Y%m_%d%H').values)
    #     fig_name = 'rain_'+flnm_time

    #     dr = Draw(fig_name, time_flag)
    #     dr.draw_main(rain_hour)
## --------------------------------------------
