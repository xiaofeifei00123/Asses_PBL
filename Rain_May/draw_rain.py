#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
将200507月28日20时的的试验结果放在左边
将201408月19日12时的的试验结果放在右边
画多个试验降水在一张图上
难点在于子图的摆放位置等等设置
和上一个的改变，是可以自定义的多了，逻辑简单了
-----------------------------------------
Time             :2021/03/24 11:38:10
Author           :Forxd
Version          :2.0
'''

import xarray as xr
import numpy as np
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


from get_rain import get_total_rain


def draw_station(ax):
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

def create_map(ax):
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
    ax.set_xticks(np.arange(80, 102 + 2, 4))
    ax.set_yticks(np.arange(26, 38 + 2, 2))
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.tick_params(axis='both', labelsize=5, direction='out')
    # -- 设置图像范围
    # ax.set_extent([78, 98, 26, 38], crs=ccrs.PlateCarree())
    return ax


def draw_contourf(data, ax, cmap, title):
    # levels = np.arange(0, 66, 5)  # 设置colorbar分层
    # levels = [200, 205, 210, 215, 220, 225, 235, 245, 250]  # 需要画出的等值线
    # levels = np.arange(0, 66, 5)  # 设置colorbar分层
    levels = np.arange(60, 500, 5)  # 设置colorbar分层
    ax = create_map(ax)
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
    # ax.set_title(title[2], loc='left',  fontsize=12)
    ax.set_title(title[0], loc='left',  fontsize=12)
    ax.set_extent([78, 102, 26, 38], crs=ccrs.PlateCarree())
    ax.xaxis.set_tick_params(labelsize=10)
    ax.yaxis.set_tick_params(labelsize=10)
    # ax.text(99.5, 38.5, title[1])
    ax.set_title(title[1], loc='right',  fontsize=12)
    # ax.set_tilte(title[1], loc='right')
    # ax.text(78,26.2,title[0], fontsize=12)
    draw_station(ax)
    return crx

def draw_obs():

    # flnm_2005 = "./Data/total_rain2005.nc"
    # flnm_2014 = "./Data/total_rain2014.nc"

    # # ftime_2005 = get_time_list('2005')
    # # ftime_2014 = get_time_list('2014')
    # ds_2005 = xr.open_dataset(flnm_2005)
    # ds_2014 = xr.open_dataset(flnm_2014)

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
    axes = [None] * 5  # 设置一个维度为8的空列表
    # # print(axes)
    axes[0] = fig.add_subplot(grid[0, 0:1], projection=proj)
    axes[1] = fig.add_subplot(grid[0, 1:2], projection=proj)
    axes[2] = fig.add_subplot(grid[1, 0:1], projection=proj)
    axes[3] = fig.add_subplot(grid[1, 1:2], projection=proj)
    axes[4] = fig.add_subplot(grid[2, 0:1], projection=proj)
    # axes[5] = fig.add_subplot(grid[2, 1:2], projection=proj)
    # axes[6] = fig.add_subplot(grid[3, 0:1], projection=proj)
    # axes[7] = fig.add_subplot(grid[3, 1:2], projection=proj)
    # axes[8] = fig.add_subplot(grid[4, 0:1], projection=proj)
    # axes[9] = fig.add_subplot(grid[4, 1:2], projection=proj)

    cmap = cmaps.precip3_16lev

    time_2005 = '2800_2906 Jul 2005'
    time_2014 = '1900_2000 Aug 2014'
    rain = get_total_rain()

    # rain_QNSE = ds_QNSE.totrain.salem.roi(shape=shp)
    draw_contourf(rain['MYJ'].salem.roi(shape=shp), axes[0], cmap, ['MYJ', '(a)',time_2005])
    draw_contourf(rain['MYJ'].salem.roi(shape=shp), axes[1], cmap, ['YSU', '(b)',time_2005])
    draw_contourf(rain['QNSE'].salem.roi(shape=shp), axes[2], cmap, ['QNSE', '(c)',time_2005])
    draw_contourf(rain['QNSE_EDMF'].salem.roi(shape=shp), axes[3], cmap, ['QNSE_EDMF', '(d)',time_2005])
    cf = draw_contourf(rain['TEMF'].salem.roi(shape=shp), axes[4], cmap, ['TEMF', '(e)',time_2005])
    # draw_contourf(ds_2014.OBS.salem.roi(shape=shp), axes[1], cmap, ['OBS', '(f)',time_2014])
    # draw_contourf(ds_2014.YSU.salem.roi(shape=shp), axes[3], cmap, ['YSU', '(g)',time_2014])
    # draw_contourf(ds_2014.QNSE.salem.roi(shape=shp), axes[5], cmap, ['QNSE', '(h)',time_2014])
    # draw_contourf(ds_2014.QNSE_EDMF.salem.roi(shape=shp), axes[7], cmap, ['QNSE_EDMF', '(i)',time_2014])
    # cf = draw_contourf(ds_2014.TEMF.salem.roi(shape=shp), axes[9], cmap, ['TEMF', '(j)',time_2014])
    ax6 = fig.add_axes([0.18, 0.06, 0.7, 0.02])  # 重新生成一个新的坐标图

    # ## 画色标
    bounds = np.arange(5, 65, 5)
    bounds = np.arange(60, 500, 50)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

    cb = fig.colorbar(
        cf,
        cax=ax6,
        orientation='horizontal',
        ticks=bounds,
        # fraction = 0.1,  # 色标大小
        pad=0.1,  # 
        # extend='both'
    )
    fig.savefig('comapre')

if __name__ == '__main__':
    draw_obs()
    # draw_obs('2014')
