import xarray as xr

flnm = '/mnt/zfm/Fengx/Assess_pbl_data/200507_new/QNSE/wrfout_d02_2005-07-27_22:00:00'

ds = xr.open_dataset(flnm)
print(ds.RAINC.values.max())