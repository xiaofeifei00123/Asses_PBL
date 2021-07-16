#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
对tbb亮温是0.的缺省处理
读卫星的星下点经纬度的处理
np.arange() 对步长为0.1的数据处理产生问题，
使用np.around() 限制小数位
-----------------------------------------
Time             :2021/04/12 08:53:15
Author           :Forxd
Version          :1.0
'''


import numpy as np
import xarray as xr
from pandas import Series, DataFrame
import pandas as pd

import cmaps  # 设置色标的
import xesmf as xe  # 插值的
import os  # 操作文件的

def read_tbb_bin(flmn):
    """读取卫星的TBB资料，筛选出范围内的数据，返回成DataSet的形式
    Args:
        flmn,观测的TBB的小时数据
    Returns:
        [type]: 返回的是带有坐标属性的DataSet格式的原始卫星TBB数据
    """
    # fx = "./NOM_ITG_2288_2288(0E0N)_LE.dat"  # 该卫星云顶亮温的描述文件
    fx = '/home/fengxiang/Project/Asses_pbl_July/Data/Process/NOM_ITG_2288_2288(0E0N)_LE/NOM_ITG_2288_2288(0E0N)_LE.dat'

    ## ----------------------------------------------
    ## -- 获取卫星资料的描述文件
    ## -- 就是对应的格点数据的经纬度
    ## ----------------------------------------------
    with open(fx, 'r') as f:
        data = np.fromfile(f, dtype=np.float32, count=2288*2288)  # 读经度lon原数据
        data1 = data.reshape(2288,2288)  # 根据卫星数据TBB的格点分布,重新分成这个维度
        data = np.fromfile(f, dtype=np.float32, count=2288*2288)  # 读纬度lat
        data2 = data.reshape(2288,2288) 

    # 根据星下点经度，对原始经纬度数据进行校正, 
    lonCenter                    = 104.5 # footpoint of FY-2C\E\F
    latCenter                    = 0
    # 这里的lat2d周边有的点是没有数据的，也不是规则的
    lon2d                        = data1+lonCenter  # 二维的经度坐标
    lat2d                        = data2+latCenter

    # print(type(data))
    # print(data1)
    ## ----------------------------------------------
    ## 读TBB数据
    # TBB = xr.open_dataarray(target_dir+target_file)
    TBB = xr.open_dataarray(flmn)
    tbb = TBB.values  # 将云顶亮温转为numpy的DataArray格式

    print(lat2d[500][500])
    print(lat2d[500][501])
    print(lat2d[500][502])
    ## 设置区域范围
    # lat_s = 26
    # lat_n = 38
    # lon_w = 80
    # lon_e = 98

    # lat_s = 20
    # lat_n = 45 
    # lon_w = 75
    # lon_e = 120
    lat_s = 10
    lat_n = 45 
    lon_w = 80
    lon_e = 120
    ## 筛选范围内的云顶亮温, 直接使用numpy切片
    ## 这一步是比较没想到的
    index = np.where((lon2d<lon_e)&(lon2d>lon_w)&(lat2d>lat_s)&(lat2d<lat_n))  # 这个是numpy的where语法，返回的是各维度坐标，组成的元组
    # print(index)
    # tbb[loc] = None
    tbb_reset = tbb[index[0][0]:index[0][-1], index[1][-1]:index[1][0]]  # 筛选的范围内的云顶亮温TBB
    lat_reset = lat2d[index[0][0]:index[0][-1], index[1][-1]:index[1][0]]  # 筛选的范围内的云顶亮温TBB
    lon_reset = lon2d[index[0][0]:index[0][-1], index[1][-1]:index[1][0]]  # 筛选的范围内的云顶亮温TBB
    # print(lon_reset)
    # print(lat_reset)
    # 将数组变成连续的，加快读写速度
    tbb_reset = np.ascontiguousarray(tbb_reset)
    lat_reset = np.ascontiguousarray(lat_reset)
    lon_reset = np.ascontiguousarray(lon_reset)

    # loc = np.where(tbb_reset<100)
    # tbb_reset[loc] = None
    # print(tbb_reset[loc])
    # 为TBB赋坐标属性
    tbb_return = xr.Dataset(
                    {
                        'tbb':(['south_north', 'west_east'], tbb_reset)
                    },
                    coords={
                        'lat':(['south_north','west_east'], lat_reset),
                        'lon':(['south_north','west_east'], lon_reset),
                        },
                        attrs={'variable':'tbb' },)

    # print(tbb_return)
    # # print(tbb_return['tbb'])
    return tbb_return

def regrid(dataset):
    """利用xESMF库，将非标准格点的数据，插值到标准格点上去
    Args:
        dataset ([type]): Dataset格式的数据, 这里是由TBB构成的DataSet
    读的是80-102度的数据
    """

    ## 创建ds_out, 利用函数创建,这个东西相当于掩膜一样
    # ds_out = xe.util.grid_2d(79.875,98,0.25,25.875,38,0.25) # 80-98E,26-38N
    ds_out = xe.util.grid_2d(79.875,102,0.25,25.875,38,0.25) # 80-98E,26-38N
    # ds_out = xe.util.grid_2d(79.95,98,0.1,25.95,38,0.1) # 80-98E,26-38N
    # ds_out = xe.util.grid_2d(79.95,102,0.1,25.95,38,0.1) # 80-98E,26-38N
    regridder = xe.Regridder(dataset, ds_out, 'bilinear')  # 好像是创建了一个掩膜一样
    dp = dataset['tbb']  # 获取变量
    dp_out = regridder(dp)  # 返回插值后的变量

    dp_val = dp_out.values
    loc = np.where(dp_val<100)
    dp_val[loc] = None


    ## 将二维的附属变量(经纬度), 变成一维的, 为了画图方便
    # np.around(data,n) 保留n位小数
    # lat = np.around(np.arange(26, 38+0.01, 0.1), 1)
    # lon = np.around(np.arange(80, 102+0.01, 0.1), 1)
    lat = np.around(np.arange(26, 38+0.01, 0.25), 2)
    lon = np.around(np.arange(80, 102+0.01, 0.25), 2)
    # print(lon[-1])
    # dp_return_dataarray = xr.DataArray(dp_out.values, coords=[lat, lon], dims=('lat', 'lon'))
    dp_return_dataarray = xr.DataArray(dp_val, coords=[lat, lon], dims=('lat', 'lon'))
    # # print(dp_return)
    # dp_return_dataset = xr.Dataset({"tbb": dp_return_dataarray})
    # # return dp_return_dataset
    return dp_return_dataarray


def get_tbb_obs(fnm):
    """得到插值过后的tbb

    Args:
        fnm ([type]): tbb原始文件

    Returns:
        dp_dataset: 插值后的tbb
    """

                  
    tbb_origin = read_tbb_bin(fnm)  # 从原始数据中读出来的
    print(tbb_origin)
    
    tbb_regrid = regrid(tbb_origin)  # 插值到标准经纬度上

    dp_dataset = xr.Dataset({"tbb": tbb_regrid})
    return dp_dataset

if __name__ == '__main__':
    
    ## 2014年8月19日09-15时TBB平均值
    target_dir="/mnt/zfm/Fengx/Assess_pbl_data/TBB_data/2014/"
    # target_file = "FY2C_TBB_IR1_NOM_20050729_0300.hdf"  # 卫星云顶亮温
    str_file = "FY2E_TBB_IR1_NOM_2014"

    file_time = [ '0819_1000',]
    # fnm = target_dir+str_file+file_time[0]+".hdf"
    flnm = '/home/fengxiang/Project/Asses_pbl_July/Data/TBB_FY2E_201607/'
    tbb_dataset = get_tbb_obs(flnm)
    # tbb_array = tbb_dataset.tbb
    # tbb = tbb_array.values
    # loc = np.where(tbb<100)

    # aa = read_tbb_bin(fnm)
    # tbb_dataset.to_netcdf("./Data/TBB_Obs_2812_15.nc")
    # bb = aa.tbb
    # bb = tbb.tbb.values
    # # print(bb)
    # loc = np.where(bb<100)
    # bb[loc] = None
    # # print(bb)
    # print(bb[0,-1])
    # # print(bb)
    # cc = bb.sel(lat=37, lon=102)
    # print(cc)
    # for i in cc.lon:
    #     print(i.values)
    # print(cc.lon)
    # dd = cc.sel(lon=101.9)
    # print(dd)
    # 取特定经纬度的值
    # print(bb[dict(lat=32, lon=102)])
    # (aa.tbb.sel(lat=30,lon=100))
    # print(aa.lat)
    # print(aa.lat)
