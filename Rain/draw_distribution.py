#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
画月平均降水分布图
-----------------------------------------
Time             :2021/06/15 17:21:14
Author           :Forxd
Version          :1.0
'''

from pandas.core import frame
from shapely.ops import transform
import xarray as xr
import pandas as pd
import cmaps
import numpy as np
import xesmf as xe
import os
import wrf
import netCDF4 as nc

from datetime import datetime
import datetime as Datetime

import salem
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.io.shapereader import Reader, natural_earth
import matplotlib as mpl
from matplotlib.path  import Path
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
import geopandas
import cmaps
import meteva.base as mb

from read_data import Return_data as rd
# from read_data import  get_rain_total, get_rain_specific_day
from get_cmap import get_cmap_rain

class Draw(object):
    def __init__(self,fig_name, flag) -> None:
        self.fig_name = fig_name
        self.flag = flag

    def draw_station(self,ax):
        station = {
            'NQ': { 'lat': 31.4, 'lon': 92.0 },
            'LS': { 'lat': 29.6, 'lon': 91.1 },
            'TTH': { 'lat': 34.2, 'lon': 92.4 },
            'GZ': { 'lat': 32.3, 'lon': 84.0 },
            'SZ': { 'lat': 30.9, 'lon': 88.7 },
            'SQH': { 'lat': 32.4, 'lon': 80.1 },
        }
        values = station.values()
        station_name = list(station.keys())
        # print(type(station_name[0]))
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
            # print(x[i])
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


    def draw_contourf(self,data, ax, cmap, title):
        # levels = np.arange(0, 66, 5)  # 设置colorbar分层
        # levels = [200, 205, 210, 215, 220, 225, 235, 245, 250]  # 需要画出的等值线
        # levels = np.arange(0, 66, 5)  # 设置colorbar分层
        # levels = np.arange(0, 110, 5)  # 设置colorbar分层
        # levels = np.arange(60, 500, 5)  # 设置colorbar分层
        # levels = np.arange(0, 481, 30)  # 设置colorbar分层
        # levels = np.arange(0, 481, 40)  # 设置colorbar分层
        levels = np.arange(0, 31, 2.5)  # 设置colorbar分层

        ax = self.create_map(ax)  # 创建坐标图像
        # print(title[0])
        # if title[0] == 'obs' and self.flag == 'all':
        if title[0] == 'obs':
            pass
            self.draw_patch(ax)
            print("yes")

        x = data.lon
        y = data.lat
        crx = ax.contourf(x,
                        y,
                        data,
                        cmap=cmap,
                        # extend='both',
                        extend='max',
                        levels=levels,
                        # levels=12,
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

        area = [None] * 3  # 设置一个维度为8的空列表
        area[0] = {"lat1":33.5, "lat2":40, "lon1":80, "lon2":90}  # north
        area[1] = {"lat1":28, "lat2":33, "lon1":83, "lon2":94}  # south left 
        area[2] = {"lat1":26, "lat2":33, "lon1":95, "lon2":103}  # south right
        for i in range(3):
            lon = np.empty(4)
            lat = np.empty(4)
            lon[0], lat[0] = area[i]['lon1'],area[i]['lat1']   # lower left (ll)
            lon[1], lat[1] = area[i]['lon2'],area[i]['lat1']   # lower right (ll)
            lon[2], lat[2] = area[i]['lon2'],area[i]['lat2']   # upper right(ll)
            lon[3], lat[3] = area[i]['lon1'],area[i]['lat2']   # upper left (ll)
            x, y = lon, lat
            xy = list(zip(x, y))
            # print(xy)
            poly = plt.Polygon(xy, edgecolor="red", fc="none", lw=.9, alpha=1)
            ax.add_patch(poly)
            # plt.gca().add_patch(poly)

    def draw_obs(self,rain):

        ## 获取地形文件
        shp_file = '/home/fengxiang/Data/shp_tp/Tibet.shp'
        shp = geopandas.read_file(shp_file)
        # module_list = ['YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        model_list = ['QNSE', 'QNSE_EDMF', 'TEMF','MYJ','YSU']

        ## --->画图
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

        cmap = get_cmap_rain()

        time_2005 = '2800_2906 Jul 2005'
        time_2014 = '1900_2000 Aug 2014'
        # rain = self.get_total_rain()


        # print(rain['obs'])
        # --------------------------------------------------
        # #### test
        # ds = xr.open_dataset('/mnt/Disk4T_5/fengxiang_file/Data/ERA5/temp/gpcp_v02r03_monthly_d201607.nc')
        # da = ds.precip
        # da = da.rename({'longitude':'lon', 'latitude':'lat'})
        # da = da.squeeze()
        # print(da)
        # print(rain['obs'])
        # rain['obs'] = da*30
        # --------------------------------------------------

        # rain_QNSE = ds_QNSE.totrain.salem.roi(shape=shp)
        self.draw_contourf(rain['MYJ'].salem.roi(shape=shp), axes[0], cmap, ['MYJ', '(a)',time_2005])
        self.draw_contourf(rain['MYJ'].salem.roi(shape=shp), axes[1], cmap, ['YSU', '(b)',time_2005])
        self.draw_contourf(rain['QNSE'].salem.roi(shape=shp), axes[2], cmap, ['QNSE', '(c)',time_2005])
        self.draw_contourf(rain['QNSE_EDMF'].salem.roi(shape=shp), axes[3], cmap, ['QNSE_EDMF', '(d)',time_2005])
        self.draw_contourf(rain['TEMF'].salem.roi(shape=shp), axes[4], cmap, ['TEMF', '(e)',time_2005])
        cf = self.draw_contourf(rain['obs'].salem.roi(shape=shp), axes[5], cmap, ['obs', '(f)',time_2005])
        # cf = self.draw_contourf(rain['obs'], axes[5], cmap, ['obs', '(f)',time_2005])
        ax6 = fig.add_axes([0.18, 0.06, 0.7, 0.02])  # 重新生成一个新的坐标图

        ### 在obs上画边框
        # mpath.Path.unit_re
        # axes[5].add_patch()
        # self.draw_patch(axes[5])

        cb = fig.colorbar(
            cf,
            cax=ax6,
            orientation='horizontal',
            # ticks=bounds,
            # fraction = 0.1,  # 色标大小
            pad=0.1,  # 
            # extend='both'
        )
        # fig.savefig('comapre')
        path = '/home/fengxiang/Project/Asses_PBL/Draw/Rain/'
        fig.savefig(path+self.fig_name)
        

if __name__ == '__main__':
    pass



    # #### test
    # ds = xr.open_dataset('/mnt/Disk4T_5/fengxiang_file/Data/ERA5/temp/gpcp_v02r03_monthly_d201607.nc')
    # da = ds.precip
    # da = da.rename({'longitude':'lon', 'latitude':'lat'})
    # da = da.squeeze()
    # print(da)




    flag = 'night'
    # flag = 'day'
    # flag = 'all'
    area = {"lat1":24.875, "lat2":40.125, "lon1":69.875, "lon2":105.125}
    # rain = get_rain_total(flag, area)  # 这个和area, flag也应该没有关系
    rd = rd(flag, area)
    # rain = rd.get_rain_specific_day('2016/07/14 00')
    rain = rd.get_rain_total()
    
    print(rain['obs'])


    # #### 画图
    fig_name = 'rain_'+flag
    dr = Draw(fig_name, flag)
    dr.draw_obs(rain)