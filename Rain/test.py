# %%
import xarray as xr
from xgrads import CtlDescriptor
from xgrads import open_CtlDataset
import numpy as np

# %%
flnm = '/mnt/zfm_18T/fengxiang/DATA/PRECIPTATION/CMORPH_STATION_RAIN/07/CHN_PRCP_HOUR_MERG_DISPLAY_0.1deg.lnx.ctl'
ds = open_CtlDataset(flnm)
da = ds.crain.squeeze(drop=True)
da = da.where(da.values>=0, np.nan)

ds_in = da.to_dataset()
ds_in = ds_in.drop_dims('time')
# %%
ds_in