#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
import xarray as xr
import numpy as np
import salem


import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter,LatitudeFormatter
from cartopy.io.shapereader import Reader, natural_earth
import cartopy.feature as cf
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import geopandas
import cmaps

from matplotlib.path import Path
import matplotlib.patches as patches

def draw_screen_poly( lats, lons):
    '''
    lats: 纬度列表
    lons: 经度列表
    purpose:  画区域直线
    '''
    x, y =  lons, lats 
    xy = list(zip(x,y))
    print(xy)
    poly = plt.Polygon( xy, edgecolor="blue",fc="none", lw=2, alpha=1)
    plt.gca().add_patch(poly)

def create_map():
    """创建一个包含青藏高原区域的Lambert投影的底图

    Returns:
        ax: 坐标图对象
    """
    ## --创建画图空间
    proj = ccrs.PlateCarree()  # 创建坐标系
    proj_lambert = ccrs.LambertConformal(
                    central_longitude=90.0,
                    central_latitude=33.0,
                    standard_parallels=(30,60,),
                    cutoff=-30,
                    false_easting=2237500,
                    false_northing=1987500,
                )

    ## 创建坐标系
    fig = plt.figure(figsize=(4, 4), dpi=500)  # 创建页面
    ax = fig.add_axes([0.1,0.1,0.8,0.8], projection=proj_lambert)  

    ## 读取青藏高原地形文件
    Tibet = cfeat.ShapelyFeature(
        Reader('/home/fengxiang/Data/shp_tp/Tibet.shp').geometries(),
        proj, edgecolor='k',
        facecolor='none', alpha=0.9
        )

    ## 将青藏高原地形文件加到地图中区
    ax.add_feature(Tibet, linewidth=0.5, zorder=2)

    ## --设置网格属性, 不画默认的标签
    gl=ax.gridlines(draw_labels=True,linestyle=":",linewidth=0.3 ,x_inline=False, y_inline=False,color='k')

    ## 关闭上面和右边的经纬度显示
    gl.top_labels=False #关闭上部经纬标签                                  
    # gl.bottom_labels = False
    # gl.left_labels = False

    gl.right_labels=False
    gl.xformatter = LONGITUDE_FORMATTER  #使横坐标转化为经纬度格式            
    gl.yformatter = LATITUDE_FORMATTER                                        
    gl.xlocator=mticker.FixedLocator(np.arange(60,120,10))                              
    gl.ylocator=mticker.FixedLocator(np.arange(10,60,10)) 
    gl.xlabel_style={'size':4}#修改经纬度字体大小                             
    gl.ylabel_style={'size':4}
    ax.spines['geo'].set_linewidth(0.6)#调节边框粗细
    # ax.set_extent([60, 120, 10, 60], crs=proj)
    ax.set_extent([0, 1700000*3, 0, 1700000*3], crs=proj_lambert)
    return ax
    
# ax = create_map()


def get_information(flnm, ax):
    """根据namelist.wps文件，获取地图的基本信息

    Args:
        flnm ([type]): [description]

    Returns:
        [type]: [description]
    """
    ## getting namelist.wps domain information
    name_dict={}

    with open(flnm) as fr:
        for line in fr:
            if "=" in line:   # 这里没有考虑注释的那些行吧, 不过wps一般也没人注释就是了
                line=line.replace("=","").replace(",","")
                name_dict.update({line.split()[0]: line.split()[1:]})  # 这个字典直接可以更新

    # dx = float(name_dict["dx"][0])/1000  # 转换为公里
    # dy = float(name_dict["dy"][0])/1000
    dx = float(name_dict["dx"][0])  # 转换为公里
    dy = float(name_dict["dy"][0])
    max_dom = int(name_dict["max_dom"][0])
    # print(max_dom)
    parent_grid_ratio = list(map(int, name_dict["parent_grid_ratio"]))
    i_parent_start = list(map(int, name_dict["i_parent_start"]))
    j_parent_start = list(map(int, name_dict["j_parent_start"]))
    e_sn = list(map(int, name_dict["e_sn"]))
    e_we = list(map(int, name_dict["e_we"]))
    ref_lat=  float(name_dict["ref_lat"][0])  # 模式区域中心位置
    ref_lon=  float(name_dict["ref_lon"][0])
    truelat1 = float(name_dict["truelat1"][0])  # 和投影相关的经纬度
    truelat2 = float(name_dict["truelat2"][0])

    # print(info['dx'])
    ## plot center position 画中心点
    ## 默认一个纬度是110公里，经度90km
    dlat = 110
    dlon = 90

    print(ref_lat)
    print(ref_lon)
    

    cenlon= np.arange(max_dom) 
    cenlat=np.arange(max_dom)
    # print(dx)
    # print(e_we)
    cenlon_model=dx*(e_we[0]-1)/2.0  # 中心点偏离边界的距离
    cenlat_model=dy*(e_sn[0]-1)/2.0
    print(cenlon_model)
    print(cenlat_model)
    # cenlon[0], cenlat[0]=m(cenlon_model, cenlat_model, inverse=True)
    # print(cenlon_model)
    
    if max_dom >= 2:
        ### domain 2
        # 4 corners
        ll_lon = dx*(i_parent_start[1]-1)
        ll_lat = dy*(j_parent_start[1]-1)
        ur_lon = ll_lon + dx/parent_grid_ratio[1] * (e_we[1]-1)
        ur_lat = ll_lat + dy/parent_grid_ratio[1] * (e_sn[1]-1)

        ## lower left (ll)
        lon = [1,2,3,4]
        lat = [1,2,3,4]
        lon[0],lat[0] = ll_lon, ll_lat

        ## lower right (lr)
        lon[1],lat[1] = ur_lon, ll_lat

        ## upper right (ur)
        lon[2],lat[2] = ur_lon, ur_lat

        ## upper left (ul)
        lon[3],lat[3] = ll_lon, ur_lat
        print(lat)
        print(lon)
        draw_screen_poly(lat, lon)   # 画多边型

        ## 标注d02
        plt.text(lon[3]*1, lat[3]*1., "d02")

        ### 区域2画多边形中点
        # cenlon_model = ll_lon + (ur_lon-ll_lon)/2.0
        # cenlat_model = ll_lat + (ur_lat-ll_lat)/2.0
        # cenlon[1], cenlat[1]=m(cenlon_model, cenlat_model, inverse=True)
        # # plt.plot(cenlon_model, cenlat_model,marker="o")  # 这个画的是区域2的中点
        # # plt.text(cenlon_model*0.8, cenlat_model*1.01, "({cenlat}, {cenlon})".format(cenlat=round(cenlat[1],2), cenlon=round(cenlon[1],2)))
        # plt.plot(cenlon_model,cenlat_model, marker="o", color="gray")
        # plt.text(cenlon_model*0.8, cenlat_model*1.01, "({cenlat}, {cenlon})".format(
        #                 cenlat=round(cenlat[0],2), cenlon=round(cenlon[0],2)))
        # plt.savefig('111')

    dict_return = {
                    "dx":dx,
                    "dy":dy,
                    "max_dom":max_dom,
                    "parent_grid_ratio":parent_grid_ratio,
                    "j_parent_start":j_parent_start,
                    "i_parent_start":i_parent_start,
                    "e_sn":e_sn,
                    "e_we":e_we,
                    'ref_lat':ref_lat,
                    'ref_lat':ref_lon,
                    'truelat1':truelat1,
                    'truelat2':truelat2,
                }
    return dict_return
    


if __name__ == '__main__':
    proj = ccrs.PlateCarree()  # 创建坐标系
    ax = create_map()
    # ax.scatter(50,90)
    # ax.scatter(x=10,y=0,s=1,c=1, transform=proj)
    # ax.scatter(x=90,y=33,s=1,c=1, transform=proj)

    # proj_lambert = ccrs.LambertConformal(
    #                 central_longitude=90,
    #                 central_latitude=30,
    #             )
    # proj_lambert = ccrs.LambertConformal(
    #                 central_longitude=90,
    #                 central_latitude=30,
    #             )
    # def draw_screen_poly( lats, lons):
    #     '''
    #     lats: 纬度列表
    #     lons: 经度列表
    #     purpose:  画区域直线
    #     '''
    #     x, y =  lons, lats 
    #     xy = list(zip(x,y))
    #     # print(xy)
    #     poly = plt.Polygon( xy, edgecolor="blue",fc="none", lw=2, alpha=1)
    #     plt.gca().add_patch(poly)

    # # lats =[23, 37, 37, 23]
    # lats =[14, 51, 51, 14]
    # lons = [69, 57, 122, 110]
    # # lons = [80, 80, 90, 90]
    # x, y =  lons, lats 
    # xy = list(zip(x,y))
    # # print(xy)
    # # poly = ax.Polygon(xy)
    # poly = plt.Polygon( xy, edgecolor="blue",fc="none", lw=1, alpha=1,transform=proj)
    # # poly = plt.Polygon( xy, edgecolor="blue",fc="none", lw=1, alpha=1,transform=proj_lambert)
    # # poly = plt.Polygon( xy, edgecolor="blue",fc="none", lw=1, alpha=1)
    # plt.gca().add_patch(poly)
    




    file_folder="./"
    file_name="namelist.wps"
    flnm=file_folder+file_name


    info = get_information(flnm, ax)

    plt.plot(850000, 925000,marker="o")  # 这个画的是区域2的中点
    plt.plot(0,0,marker='o')
    plt.savefig("ddddddd")
    # # print(info['dx'])
    # ## plot center position 画中心点
    # cenlon= np.arange(info['max_dom']) 
    # print(cenlon)
    # cenlat=np.arange(info['max_dom'])
    # cenlon_model=dx*(e_we[0]-1)/2.0
    # cenlat_model=dy*(e_sn[0]-1)/2.0
    # # cenlon[0], cenlat[0]=m(cenlon_model, cenlat_model, inverse=True)
    # ## 画区域1的中点和标注
    # # plt.plot(cenlon_model,cenlat_model, marker="o", color="gray")
    # # plt.text(cenlon_model*0.8, cenlat_model*1.01, "({cenlat}, {cenlon})".format(cenlat=round(cenlat[0],2), cenlon=round(cenlon[0],2)))