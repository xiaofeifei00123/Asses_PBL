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


from global_variable import station_dic
from skewt_q import Draw_skewt
from data_process import TransferData

def draw_single():
    """画单独的比湿廓线图
    """

    month = 'May'
    # month = 'Jul'
    for key in station_dic:
        station = station_dic[key]
        Dr = Draw_skewt(station, month)
        time_index = ['00', '12']
        for time_select in time_index:
            tr = TransferData(station, month)
            model_dic_q = tr.transfer_data('q')
            title = {'time': time_select, 'station_name': station['name']}
            dr = Draw_skewt(station, month)
            dr.draw_upar_single(
                model_dic_q,
                title,
                )

if __name__ == '__main__':
    draw_single()
    


