#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
'''
    File name: draw_wrf_domain.py
    Author: Liang Chen
    E-mail: chenliang@tea.ac.cn
    Date created: 2016-12-22
    Date last modified: 2021-3-3
    ##############################################################
    Purpos:
    this function reads in namelist.wps and plot the wrf domain
'''

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from matplotlib.colors import LinearSegmentedColormap
import shapefile
from matplotlib.collections import LineCollection
import matplotlib.colors
import sys
import numpy as np


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


sShapeFiles="/home/fengxiang/Data/shp_tp/"
shape_line=['Tibet.shp',]

## setting namelist.wps domain information
file_folder="./"
file_name="namelist.wps_l"
sfile=file_folder+file_name
name_dict={}

with open(sfile) as fr:
    for line in fr:
        if "=" in line:   # 这里没有考虑注释的那些行吧, 不过wps一般也没人注释就是了
           line=line.replace("=","").replace(",","")
           name_dict.update({line.split()[0]: line.split()[1:]})  # 这个字典直接可以更新

dx = float(name_dict["dx"][0])
dy = float(name_dict["dy"][0])
max_dom = int(name_dict["max_dom"][0])
print(max_dom)
parent_grid_ratio = list(map(int, name_dict["parent_grid_ratio"]))
i_parent_start = list(map(int, name_dict["i_parent_start"]))
j_parent_start = list(map(int, name_dict["j_parent_start"]))
e_sn = list(map(int, name_dict["e_sn"]))
e_we = list(map(int, name_dict["e_we"]))
ref_lat=  float(name_dict["ref_lat"][0])  # 模式区域中心位置
ref_lon=  float(name_dict["ref_lon"][0])
truelat1 = float(name_dict["truelat1"][0])  # 和投影相关的经纬度
truelat2 = float(name_dict["truelat2"][0])



# # ## draw map

fig = plt.figure(figsize=(7,6))   # 设置画板大小
#Custom adjust of the subplots
plt.subplots_adjust(left=0.05,right=0.97,top=0.9,bottom=0.1)  # 调整画布大小
ax = plt.subplot(111)
m = Basemap(resolution="l", projection="lcc", rsphere=(6370000.0, 6370000.0), lat_1=truelat1, lat_2=truelat2, lat_0=ref_lat, lon_0=ref_lon, width=dx*(e_we[0]-1), height=dy*(e_sn[0]-1))

# m.drawcoastlines()
#m.drawcountries(linewidth=2)
#m.drawcountries()
#m.fillcontinents()
#m.fillcontinents(color=(0.8,1,0.8))
#m.drawmapboundary()
#m.fillcontinents(lake_color="aqua")
#m.drawmapboundary(fill_color="aqua")



### 根据地形文件，画底图
ii=0  # 控制变量
for sr in shape_line:
    # print(sr)
    r = shapefile.Reader(sShapeFiles+sr)  # 读地形文件
    shapes = r.shapes()
    records = r.records()

    for record, shape in zip(records,shapes):
        lons,lats = zip(*shape.points)
        data = np.array(m(lons, lats)).T
        if len(shape.parts) == 1:
            segs = [data,]
        else:
            segs = []
            for i in range(1,len(shape.parts)):
                index = shape.parts[i-1]
                index2 = shape.parts[i]
                segs.append(data[index:index2])
            segs.append(data[index2:])


        lines = LineCollection(segs,antialiaseds=(1,))
#       lines.set_facecolors(cm.jet(np.random.rand(1)))

        if ii==0:
            lines.set_edgecolors('black')
            lines.set_linewidth(2)
        else:
            lines.set_edgecolors('k')
            lines.set_linewidth(1)
        ax.add_collection(lines)
    ii=ii+1

## 画标签
m.drawparallels(np.arange(-90, 90, 10), labels = [1,0,0,0], fontsize=16,dashes=[1,1])
# m.drawmeridians(np.arange(-180, 180, 10), labels = [0,0,0,1], fontsize=16,dashes=[1,1])
print(ref_lat, ref_lon)

## plot center position 画中心点
cenlon= np.arange(max_dom); cenlat=np.arange(max_dom)
cenlon_model=dx*(e_we[0]-1)/2.0
cenlat_model=dy*(e_sn[0]-1)/2.0
cenlon[0], cenlat[0]=m(cenlon_model, cenlat_model, inverse=True)
## 画区域1的中点和标注
plt.plot(cenlon_model,cenlat_model, marker="o", color="gray")
plt.text(cenlon_model*0.8, cenlat_model*1.01, "({cenlat}, {cenlon})".format(cenlat=round(cenlat[0],2), cenlon=round(cenlon[0],2)))


#### draw nested domain rectangle
#### 区域2
#### 画多边形
lon=np.arange(4); lat=np.arange(4)

if max_dom >= 2:
    ### domain 2
    # 4 corners
    ll_lon = dx*(i_parent_start[1]-1)
    ll_lat = dy*(j_parent_start[1]-1)
    ur_lon = ll_lon + dx/parent_grid_ratio[1] * (e_we[1]-1)
    ur_lat = ll_lat + dy/parent_grid_ratio[1] * (e_sn[1]-1)

    ## lower left (ll)
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
    cenlon_model = ll_lon + (ur_lon-ll_lon)/2.0
    cenlat_model = ll_lat + (ur_lat-ll_lat)/2.0
    cenlon[1], cenlat[1]=m(cenlon_model, cenlat_model, inverse=True)
    # plt.plot(cenlon_model, cenlat_model,marker="o")  # 这个画的是区域2的中点
    # plt.text(cenlon_model*0.8, cenlat_model*1.01, "({cenlat}, {cenlon})".format(cenlat=round(cenlat[1],2), cenlon=round(cenlon[1],2)))

if max_dom >= 3:
    ### domain 3
    ## 4 corners
    ll_lon += dx/parent_grid_ratio[1]*(i_parent_start[2]-1)
    ll_lat += dy/parent_grid_ratio[1]*(j_parent_start[2]-1)
    ur_lon = ll_lon +dx/parent_grid_ratio[1]/parent_grid_ratio[2]*(e_we[2]-1)
    ur_lat =ll_lat+ dy/parent_grid_ratio[1]/parent_grid_ratio[2]*(e_sn[2]-1)

    ## ll
    lon[0],lat[0] = ll_lon, ll_lat
    ## lr
    lon[1],lat[1] = ur_lon, ll_lat
    ## ur
    lon[2],lat[2] = ur_lon, ur_lat
    ## ul
    lon[3],lat[3] = ll_lon, ur_lat

    draw_screen_poly(lat, lon)
    plt.text(lon[0]-lon[0]/10,lat[0]-lat[0]/10,"({i}, {j})".format(i=i_parent_start[2], j=j_parent_start[2]))
    #plt.plot(lon,lat,linestyle="",marker="o",ms=10)
    cenlon_model = ll_lon + (ur_lon-ll_lon)/2.0
    cenlat_model = ll_lat + (ur_lat-ll_lat)/2.0

#    plt.plot(cenlon,cenlat,marker="o",ms=15)
    #print m(cenlon, cenlat)cenlon, cenlat, ll_lon, ll_lat, ur_lon, ur_lat
    #print m(cenlon, cenlat,inverse=True)

    cenlon[2], cenlat[2]=m(cenlon_model, cenlat_model, inverse=True)

if max_dom >= 4:

    ### domain 3
    ## 4 corners
    ll_lon += dx/parent_grid_ratio[1]/parent_grid_ratio[2]*(i_parent_start[3]-1)
    ll_lat += dy/parent_grid_ratio[1]/parent_grid_ratio[2]*(j_parent_start[3]-1)
    ur_lon = ll_lon +dx/parent_grid_ratio[1]/parent_grid_ratio[2]/parent_grid_ratio[3]*(e_we[3]-1)
    ur_lat =ll_lat+ dy/parent_grid_ratio[1]/parent_grid_ratio[2]/parent_grid_ratio[3]*(e_sn[3]-1)
    
    ## ll
    lon[0],lat[0] = ll_lon, ll_lat
    ## lr
    lon[1],lat[1] = ur_lon, ll_lat
    ## ur
    lon[2],lat[2] = ur_lon, ur_lat
    ## ul
    lon[3],lat[3] = ll_lon, ur_lat
    draw_screen_poly(lat, lon)
    #plt.plot(lon,lat,linestyle="",marker="o",ms=10)

    cenlon_model = ll_lon + (ur_lon-ll_lon)/2.0
    cenlat_model = ll_lat + (ur_lat-ll_lat)/2.0
#    plt.plot(cenlon,cenlat,marker="o",ms=15)
    #print m(cenlon, cenlat)cenlon, cenlat, ll_lon, ll_lat, ur_lon, ur_lat
    #print m(cenlon, cenlat,inverse=True)
    cenlon[3], cenlat[3]=m(cenlon_model, cenlat_model, inverse=True)

## 标注站点

plt.plot(cenlon_model, cenlat_model,marker="o")  # 这个画的是区域2的中点
print(cenlon_model/25000, cenlat_model/25000)
# plt.text(cenlon_model*0.8, cenlat_model*1.01, "({cenlat}, {cenlon})".format(cenlat=round(cenlat[1],2), cenlon=round(cenlon[1],2)))

cenlon_model=dx*(e_we[0]-1)/2.0
print(dx)
print(dy)
Tingri={'lat':28.6,'lon':87.0,'name':'Tingri'}
plt.plot(Tingri['lon']*25000, Tingri['lat']*25000,marker="o") 
# plt.text(Tingri['lon']*0.8, Tingri['lat']*1.01, "({cenlat}, {cenlon})".format(cenlat=round(cenlat[1],2), cenlon=round(cenlon[1],2)))


plt.savefig("tttt.png")
