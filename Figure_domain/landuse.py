#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:

-----------------------------------------
Time             :2021/08/24 20:39:29
Author           :Forxd
Version          :1.0
'''


# %%
import xarray as xr
import  matplotlib.pyplot as plt
# %%
flnm = '/mnt/zfm_18T/fengxiang/DATA/LANDUSE/geo_em.d01.nc'
ds = xr.open_dataset(flnm)
land = ds['LU_INDEX'].squeeze()
# print(land)
# %%
# land.dims
# ds['XLAT_M'].shape
# land['XLAT_M']
# land
# ds['XLAT_M'].max()
ds['XLONG_M'].max()