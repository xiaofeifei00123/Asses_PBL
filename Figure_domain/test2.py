import numpy as np
import math


def get_information(flnm):
    """根据namelist.wps文件，获取地图的基本信息

    Args:
        flnm ([type]): [description]

    Returns:
        [type]: [description]
    """
    ## getting namelist.wps domain information
    name_dict={}

    with open(flnm) as fr:
        for line in fr:
            if "=" in line:   # 这里没有考虑注释的那些行吧, 不过wps一般也没人注释就是了
                line=line.replace("=","").replace(",","")
                name_dict.update({line.split()[0]: line.split()[1:]})  # 这个字典直接可以更新

    dx = float(name_dict["dx"][0])  # 转换为公里
    dy = float(name_dict["dy"][0])
    max_dom = int(name_dict["max_dom"][0])
    # print(max_dom)
    parent_grid_ratio = list(map(int, name_dict["parent_grid_ratio"]))
    i_parent_start = list(map(int, name_dict["i_parent_start"]))
    j_parent_start = list(map(int, name_dict["j_parent_start"]))
    e_sn = list(map(int, name_dict["e_sn"]))
    e_we = list(map(int, name_dict["e_we"]))
    ref_lat=  float(name_dict["ref_lat"][0])  # 模式区域中心位置
    ref_lon=  float(name_dict["ref_lon"][0])
    truelat1 = float(name_dict["truelat1"][0])  # 和投影相关的经纬度
    truelat2 = float(name_dict["truelat2"][0])
    # print(truelat2)
    # print(e_sn)

    # print(dy)

    # print(e_sn[0])
    dlat = (e_sn[0]-1)/2*dy/1000/110  # 每一度是110km
    latn_angle = ref_lat+dlat
    lats_angle = ref_lat-dlat
    latn = (ref_lat+dlat)*math.pi/180
    lats = (ref_lat-dlat)*math.pi/180

    dlonn = (e_we[0]-1)/2*dx/1000/(math.pi*6371/180*math.cos(latn))
    dlons = (e_we[0]-1)/2*dx/1000/(math.pi*6371/180*math.cos(lats))
    
    pos1 = [lats_angle, ref_lon-dlons]
    pos2 = [latn_angle, ref_lon-dlonn]
    pos3 = [latn_angle, ref_lon+dlonn]
    pos4 = [lats_angle, ref_lon+dlons]
    # print(dlonn)
    # print(dlons)
    print(pos1)
    print(pos2)
    print(pos3)
    print(pos4)






    # print(ref_lat)
    # print(latn)
    # print(lats)


if __name__ == '__main__':

    file_folder="./"
    file_name="namelist.wps"
    flnm=file_folder+file_name
    info = get_information(flnm)
    