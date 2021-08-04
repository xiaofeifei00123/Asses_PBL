#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
绘制探空曲线的主程序
-----------------------------------------
Author           :Forxd
Version          :1.0
Time：2021/07/28/ 12:04
'''


from threading import Condition
from global_variable import station_dic
from draw_skewt import Draw_skewt
from data_transfer import TransferData

import time

def draw_single(month, condition):
    """画单独的比湿廓线图
    """

    # month = 'May'
    # month = 'Jul'
    # condition = 'rain'
    # condition = 'dry'
    for key in station_dic:
        station = station_dic[key]
        Dr = Draw_skewt(station, month)
        time_index = ['00', '12']
        for time_select in time_index:
            print("画[%s]站， [%s]时的图"%(key, time_select))
            tr = TransferData(station, month, time_select, 'q')
            ds_q = tr.transfer_data(condition)
            # print(ds_q)
            # time.sleep(100)
            if not ds_q is None:  #  判断不是None
            # if not ds_q is None :  #  判断不是None
                # ds_q = ds_q.to_dataset(dim='model')
                title = {'time': time_select, 
                        'station_name': station['name'],
                        'condition':condition}
                dr = Draw_skewt(station, month)
                dr.draw_upar_single(
                    ds_q,
                    title,
                    )

def draw_dual():
    
    month_list = ['May', 'Jul']                 
    Condition_list = ['dry', 'rain', 'all']
    for i in month_list:
        for j in Condition_list:
            draw_single(i, j)


if __name__ == '__main__':
    # draw_single()
    draw_dual()
    


