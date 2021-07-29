#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:

-----------------------------------------
Time             :2021/07/15 10:00:54
Author           :Forxd
Version          :1.0
'''

from xgrads import CtlDescriptor
from xgrads import open_CtlDataset

flnm = '/mnt/zfm_18T/Asses_PBL/CMORPH_RAIN/07/CHN_PRCP_HOUR_MERG_DISPLAY_0.1deg.lnx.ctl'
dest = open_CtlDataset(flnm)
print(dest)