#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
read in namelist.wps , draw wrf domain and plot some station
-----------------------------------------
Time             :2021/03/28 17:28:59
Author           :Forxd
Version          :1.0
'''

# %%
import xarray as xr
import numpy as np
import salem

import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.io.shapereader import Reader, natural_earth
import cartopy.feature as cf
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import geopandas
import cmaps

from matplotlib.path import Path
import matplotlib.patches as patches

# from global_variable import station_dic
from global_station import station_dic


def draw_screen_poly(lats, lons):
    '''
    lats: 纬度列表
    lons: 经度列表
    purpose:  画区域直线
    '''
    x, y = lons, lats
    xy = list(zip(x, y))
    print(xy)
    poly = plt.Polygon(xy, edgecolor="black", fc="none", lw=.8, alpha=1)
    plt.gca().add_patch(poly)


def create_map(info):
    """创建一个包含青藏高原区域的Lambert投影的底图

    Returns:
        ax: 坐标图对象
    """
    ## --创建画图空间

    ref_lat = info['ref_lat']
    ref_lon = info['ref_lon']
    true_lat1 = info['true_lat1']
    true_lat2 = info['true_lat2']
    false_easting = (info['e_we'][0] - 1) / 2 * info['dx']
    false_northing = (info['e_sn'][0] - 1) / 2 * info['dy']

    # print(true_lat1)

    proj_lambert = ccrs.LambertConformal(
        central_longitude=ref_lon,
        central_latitude=ref_lat,
        standard_parallels=(true_lat1, true_lat2),
        cutoff=-30,
        false_easting=false_easting,
        false_northing=false_northing,
    )
    # proj = ccrs.PlateCarree(central_longitude=ref_lon)  # 创建坐标系
    proj = ccrs.PlateCarree()  # 创建坐标系
    ## 创建坐标系
    fig = plt.figure(figsize=(4, 3.4), dpi=300)  # 创建页面
    ax = fig.add_axes([0.12, 0.02, 0.85, 0.95], projection=proj_lambert)

    ## 读取青藏高原地形文件
    Tibet = cfeat.ShapelyFeature(
        Reader('/home/fengxiang/Data/shp_tp/Tibet.shp').geometries(),
        proj,
        edgecolor='black',
        lw = 1.,
        linewidth = 1.,
        facecolor='none',
        alpha=1.)

    ## 将青藏高原地形文件加到地图中区
    ax.add_feature(Tibet, linewidth=0.5, zorder=2)

    ## --设置网格属性, 不画默认的标签
    gl = ax.gridlines(draw_labels=True,
                      dms=True,
                      linestyle=":",
                      linewidth=0.3,
                      x_inline=False,
                      y_inline=False,
                      color='k')
    # # gl=ax.gridlines(draw_labels=True,linestyle=":",linewidth=0.3 , auto_inline=True,x_inline=False, y_inline=False,color='k')

    ## 关闭上面和右边的经纬度显示
    gl.top_labels = False  #关闭上部经纬标签
    # gl.bottom_labels = False
    # # gl.left_labels = False
    gl.right_labels = False
    ## 这个东西还挺重要的，对齐坐标用的
    gl.rotate_labels = None

    gl.xformatter = LONGITUDE_FORMATTER  #使横坐标转化为经纬度格式
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(75, 105 + 5, 5))
    gl.ylocator = mticker.FixedLocator(np.arange(20, 50, 5))
    gl.xlabel_style = {'size': 10}  #修改经纬度字体大小
    gl.ylabel_style = {'size': 10}
    ax.spines['geo'].set_linewidth(0.6)  #调节边框粗细
    # ax.set_extent([60, 120, 10, 60], crs=proj)
    # ax.set_extent([0, 2237500*2, 0, 1987500*2], crs=proj_lambert)
    ax.set_extent([0, false_easting * 2, 0, false_northing * 2],
                  crs=proj_lambert)
    # print(false_northing)
    # print(false_easting)
    return ax


def get_information(flnm):
    """根据namelist.wps文件，获取地图的基本信息

    Args:
        flnm ([type]): [description]

    Returns:
        [type]: [description]
    """
    ## getting namelist.wps domain information
    name_dict = {}

    with open(flnm) as fr:
        for line in fr:
            if "=" in line:  # 这里没有考虑注释的那些行吧, 不过wps一般也没人注释就是了
                line = line.replace("=", "").replace(",", "")
                name_dict.update({line.split()[0]:
                                  line.split()[1:]})  # 这个字典直接可以更新

    dx = float(name_dict["dx"][0])  # 转换为公里
    dy = float(name_dict["dy"][0])
    max_dom = int(name_dict["max_dom"][0])
    # print(max_dom)
    parent_grid_ratio = list(map(int, name_dict["parent_grid_ratio"]))
    i_parent_start = list(map(int, name_dict["i_parent_start"]))
    j_parent_start = list(map(int, name_dict["j_parent_start"]))
    e_sn = list(map(int, name_dict["e_sn"]))
    e_we = list(map(int, name_dict["e_we"]))
    ref_lat = float(name_dict["ref_lat"][0])  # 模式区域中心位置
    ref_lon = float(name_dict["ref_lon"][0])
    truelat1 = float(name_dict["truelat1"][0])  # 和投影相关的经纬度
    truelat2 = float(name_dict["truelat2"][0])

    cenlon = np.arange(max_dom)
    cenlat = np.arange(max_dom)
    cenlon_model = dx * (e_we[0] - 1) / 2.0  # 中心点偏离边界的距离
    cenlat_model = dy * (e_sn[0] - 1) / 2.0

    dict_return = {
        "dx": dx,
        "dy": dy,
        "max_dom": max_dom,
        "parent_grid_ratio": parent_grid_ratio,
        "j_parent_start": j_parent_start,
        "i_parent_start": i_parent_start,
        "e_sn": e_sn,
        "e_we": e_we,
        'ref_lat': ref_lat,
        'ref_lon': ref_lon,
        'true_lat1': truelat1,
        'true_lat2': truelat2,
        'parent_grid_ratio': parent_grid_ratio,
    }
    return dict_return


def draw_d02(info):
    """绘制domain2

    Args:
        info ([type]): [description]
    """
    max_dom = info['max_dom']
    dx = info['dx']
    dy = info['dy']
    i_parent_start = info['i_parent_start']
    j_parent_start = info['j_parent_start']
    parent_grid_ratio = info['parent_grid_ratio']
    e_we = info['e_we']
    e_sn = info['e_sn']

    if max_dom >= 2:
        ### domain 2
        # 4 corners 找到四个顶点和距离相关的坐标
        ll_lon = dx * (i_parent_start[1] - 1)
        ll_lat = dy * (j_parent_start[1] - 1)
        ur_lon = ll_lon + dx / parent_grid_ratio[1] * (e_we[1] - 1)
        ur_lat = ll_lat + dy / parent_grid_ratio[1] * (e_sn[1] - 1)

        lon = np.empty(4)
        lat = np.empty(4)

        lon[0], lat[0] = ll_lon, ll_lat  # lower left (ll)
        lon[1], lat[1] = ur_lon, ll_lat  # lower right (lr)
        lon[2], lat[2] = ur_lon, ur_lat  # upper right (ur)
        lon[3], lat[3] = ll_lon, ur_lat  # upper left (ul)

        draw_screen_poly(lat, lon)  # 画多边型

        ## 标注d02
        # plt.text(lon[0] * 1+100000, lat[0] * 1. - 225000, "d02", fontdict={'size':14})
        plt.text(lon[0] * 1+100000, lat[0] * 1.+50000, "d02", fontdict={'size':14})
        print(lon[3])


def draw_station():
    """画站点
    """
    # station = {
    #     # 'TR': {
    #     #     'lat': 28.6,
    #     #     'lon': 87.0
    #     # },
    #     'NQ': {
    #         'lat': 31.4,
    #         'lon': 92.0
    #     },
    #     'LS': {
    #         'lat': 29.6,
    #         'lon': 91.1
    #     },
    #     'TTH': {
    #         'lat': 34.2,
    #         'lon': 92.4
    #     },
    #     'GZ': {
    #         'lat': 32.3,
    #         'lon': 84.0
    #     },
    #     'SZ': {
    #         'lat': 30.9,
    #         'lon': 88.7
    #     },
    #     'SQH': {
    #         'lat': 32.4,
    #         'lon': 80.1
    #     },
    #     # 'JinChuan': {
    #     #     'lat': 31.29,
    #     #     'lon': 102.04
    #     # },
    #     # 'JinLong': {
    #     #     'lat': 29.00,
    #     #     'lon': 101.50
    #     # },
    #     'TR': {
    #         'lat': 28.63,
    #         'lon': 87.08,
    #         'name': 'TingRi',
    #         'number': '55664',
    #         'height': 4302
    #     },
    #     'LZ': {
    #         'lat': 29.65,
    #         'lon': 94.36,
    #         'name': 'LinZhi',
    #         'number': '56312',
    #         'height':2991.8 
    #     },
    #     'CD': {
    #         'lat': 31.15,
    #         'lon': 97.17,
    #         'name': 'ChangDu',
    #         'number': '56137',
    #         'height': 3315
    #     },
    # }
    station = station_dic
    values = station.values()
    # station_name = list(station.keys())
    # station_name = list(station.values()['abbreviation'])
    station_name = []
    # print(type(station_name[0]))
    # print(station_name[0])
    x = []
    y = []
    for i in values:
        y.append(float(i['lat']))
        x.append(float(i['lon']))
        station_name.append(i['abbreviation'])

    ## 标记出站点
    ax.scatter(x,
               y,
               color='black',
               transform=ccrs.PlateCarree(),
               linewidth=0.2,
               s=12)
    ## 给站点加注释
    for i in range(len(x)):
        print(x[i])
        if station_name[i] == 'LS':
            # x[i] = x[i]
            y[i] = y[i] - 2.0 
        if station_name[i] == 'SQH':
            x[i] = x[i] - 0.5
        if station_name[i] == 'GEM':
            x[i] = x[i] - 2.0
        if station_name[i] == 'TTH':
            x[i] = x[i] - 1.5
        if station_name[i] == 'CD':
            y[i] = y[i] - 2.0
            x[i] = x[i] - 0.5
        if station_name[i] == 'BT':
            x[i] = x[i] + 0.5
        plt.text(x[i] - 1,
                 y[i] + 0.5,
                 station_name[i],
                 transform=ccrs.PlateCarree(),
                 fontdict={
                     'size': 9,
                 })

    plt.text(62,
             45,
             'd01',
             transform=ccrs.PlateCarree(),
             fontdict={
                 'size': 14,
             })


if __name__ == '__main__':
    file_folder = "./"
    file_name = "namelist.wps"
    flnm = file_folder + file_name

    info = get_information(flnm)  # 获取namelist.wps文件信息
    ax = create_map(info)  # 在domain1区域内，添加地理信息，创建底图
    draw_d02(info)  # 绘制domain2区域
    draw_station()  # 将站点位置绘制到图上
    # plt.title('d01', loc='left')
    # plt.savefig("domain.tiff")
    plt.savefig("domain.png")