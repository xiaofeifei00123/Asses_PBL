#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取lambert投影的metem文件中的下垫面类型数据
转换投影至lat-lon投影
获得一个下垫面地表类型数据的掩膜，即
是该地形的赋值为1， 不是该地形的为0
-----------------------------------------
Time             :2021/08/27 09:56:15
Author           :Forxd
Version          :1.0
'''


# %%
import xarray as xr
import xesmf as xe
import numpy as np
import salem
import geopandas

# %%
def get_landmask():
    """获得根据下垫面类型区分的掩膜文件
    是该地形则为1， 否则是np.nan值
    Returns:
        [type]: [description]
    """
    ## 读取下垫面类型文件
    flnm = '/mnt/zfm_18T/fengxiang/Asses_PBL/data/geo_em.d01.nc'
    ds = xr.open_dataset(flnm)
    land = ds['LU_INDEX'].squeeze().values
    lon = ds['XLONG_M'].squeeze().values
    lat = ds['XLAT_M'].squeeze().values

    ## 生成二维的nc文件, 下垫面类型的原始数据
    ds_in = xr.Dataset(
        {
            'LU_INDEX':(['x','y'],land)
        },
        coords={
            'lon':(['x','y'], lon),
            'lat':(['x','y'], lat),
        },
    )


    ## 对下垫面的文件进行转换投影和插值
    ds_out = xr.Dataset(
        {'lat':(['lat'],np.arange(24.875, 45.125+0.25, 0.25)),
            'lon':(['lon'], np.arange(69.875, 105.125+0.25, 0.25))}
        )
    regridder = xe.Regridder(ds_in, ds_out, 'bilinear')  

    dds = regridder(ds_in)  
    ds_out = dds.round(0)
    ## 去除高原外的数据
    shp_file = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_tp/Tibet.shp'
    shp = geopandas.read_file(shp_file)
    ds_tibet = ds_out.salem.roi(shape=shp)  # 青藏高原区域上的下垫面类型数据
    da = ds_tibet['LU_INDEX']

    mask_dic = {}
    mask_dic['bare'] = xr.where(da==16, 1,np.nan)
    mask_dic['bush'] = xr.where(da==7, 1,np.nan)
    mask_dic['grass'] = xr.where(da==10, 1,np.nan)
    return mask_dic

# %%
# tt = mask_dic










