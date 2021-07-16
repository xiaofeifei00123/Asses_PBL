# %%
import xarray as xr
import os
path = '/mnt/zfm_18T/Asses_PBL/FNL/201607/*_00_00.grib2'
fl_list = os.popen('ls {}'.format(path))
fl_list = fl_list.read().split()

# %%
rh_list = []
for fl in fl_list:
    ds = xr.open_dataset(fl, engine='cfgrib',
                         backend_kwargs={'filter_by_keys':
                            {'typeOfLevel': 'isobaricInhPa'}})
    rh = ds.r
    rh_list.append(rh)
da_rh = xr.concat(rh_list, dim='time')
da = da_rh.sel(latitude=33, longitude=95,
               method='nearest')
da = da.rename({'isobaricInhPa': 'pressure',
                'latitude': 'lat', 'longitude': 'lon'})
print(da)
