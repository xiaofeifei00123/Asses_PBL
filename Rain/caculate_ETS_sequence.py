#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
计算并绘制逐小时的
画SAL评分的图
二分类评分的图
计算时间序列评分
-----------------------------------------
Time             :2021/06/04 14:32:20
Author           :Forxd
Version          :1.0
'''
import xarray as xr
import meteva.method as mem
import meteva.base as meb
import numpy as np
import pandas as pd

import salem  # 过滤高原外的数据
import geopandas

# import matplotlib.pyplot as plot
import matplotlib.pyplot as plt
from matplotlib.pyplot import savefig
# get_data = Get_data('night')

from read_data import Return_data  as rd
from read_data import GetData, TransferData
# from read_data  import get_rain_total, get_rain_times, get_rain_daily


class Caculate():

    def __init__(self, rain) -> None:
        # self.flag = flag
        # self.area = area
        self.rain = rain

    def get_two_scale(self,):
        # flag = 'all'
        """计算ETS评分这些值
        """
        # rd = rd(self.flag, self.area)
        # rain = rd.get_rain_total()
        rain = self.rain
        model_list = ['MYJ', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']

        ### 计算二分类评分
        hfmc_dic = {}
        Accuracy = {}
        ETS = {}
        TS = {}
        POFD = {}
        HK = {}
        ODR = {}
        BIAS = {}
        for model in model_list:
            # hfmc = mem.hfmc(rain['obs'].values/30, rain[model].values/30, grade_list=[0.25])

            ## mask高原外的数据
            rain_obs = rain['obs'].values
            rain_model = rain[model].values
            print(rain_model.max())

            # hfmc = mem.hfmc(rain['obs'].values/30, rain[model].values/30, grade_list=[0.1])
            hfmc = mem.hfmc(rain_obs/30, rain_model/30, grade_list=[0.1])   # 这里算的都是平均值的评分, 我想要的是评分的平均值, 也就是要算每个时次的评分
            # ets = mem.ets_hfmc(hfmc)
            # ETS[model] = ets
            ETS[model] = mem.ets_hfmc(hfmc) 
            TS[model] = mem.ts_hfmc(hfmc) 
            Accuracy[model] = mem.pc_hfmc(hfmc)
            POFD[model] = mem.pofd_hfmc(hfmc)
            # POFD[model] = mem.far_hfmc(hfmc)
            HK[model] = mem.hk_yesorno_hfmc(hfmc)
            ODR[model] = mem.odds_ratio_hfmc(hfmc)
            BIAS[model] = mem.bias_hfmc(hfmc)

        grade_list = [Accuracy, ETS, TS, POFD, HK, ODR, BIAS]
        grade_list_name = ['Accuracy', 'ETS', 'TS', 'POFD', 'HK', 'ODR', 'BIAS']

        b = []
        for i in grade_list:
            a = pd.Series(i)
            b.append(a)
            # print(a)
        df = pd.concat(b, axis=1)
        df.columns = grade_list_name
        df = df.T
        print(df)
        df.to_csv('/home/fengxiang/Project/Asses_PBL/Data/tt.csv')
        return df


    def get_time_scale(self,):
        """计算时间序列的评分
        """
        # rd = rd(self.flag, self.area)
        # rain = rd.get_rain_total()
        rain = rd.get_rain_times()
        # print(rain)
        model_list = ['MYJ', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']

        MAE = {} # 平均绝对误差
        RMSE = {} # 均方根误差
        SD = {}  # 预报标准差/观测标准差
        CORR = {}
        for model in model_list:
            pass
            rain_model = rain[model].values
            rain_obs  = rain['obs'].values
            MAE[model] = mem.mae(rain_obs, rain_model)
            sd = mem.ob_fo_std(rain_obs, rain_model)
            SD[model] = sd[1]/sd[0]
            RMSE[model] = mem.rmse(rain_obs, rain_model)
            CORR[model] = mem.corr(rain_obs, rain_model)
        grade_list = [MAE, RMSE, SD, CORR]
        grade_list_name = ['MAE', 'RMSE', 'SD', 'CORR']
        # grade_list_name = ['Accuracy', 'TS', 'POFD', 'HK', 'ODR', 'BIAS']

        b = []
        for i in grade_list:
            a = pd.Series(i)
            b.append(a)
            # print(a)
        df = pd.concat(b, axis=1)
        df.columns = grade_list_name
        df = df.T
        print(df)
        df.to_csv('/home/fengxiang/Project/Asses_PBL/Data/time_score.csv')
        return df

    def get_space_scale(self):
        """计算空间评分
        """
        # rd = rd(self.flag, self.area)
        rain = rd.get_rain_total()
        model_list = ['MYJ', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']

        sal_list = []
        for model in model_list:
            pass
            # tttime = pd.date_range()
            rain_model = rain[model]
            rain_obs  = rain['obs']
            grid_obs = meb.xarray_to_griddata(rain_obs)
            grid_model = meb.xarray_to_griddata(rain_model)

            look_ff = mem.mode.feature_finder(grid_obs, grid_model, smooth=10, threshold=20, minsize=5)
            look_match = mem.mode.centmatch(look_ff)
            look_merge = mem.mode.merge_force(look_match)
            sal = mem.sal(look_ff)
            ssal = pd.Series(sal)
            # print(ssal)
            # print(ssal)
            sal_list.append(ssal)
            # print(sal)
        
        df = pd.concat(sal_list, axis=1)
        df.columns = model_list
        print(df)
        df.to_csv('/home/fengxiang/Project/Asses_PBL/Data/sal.csv')
        return df
        # print(tt)

class Draw():

    # def __init__(self, flag_list, area_dic) -> None:
    #     self.flag_list = flag_list
    #     self.area_dic = area_dic

    # dic = get_rain_daily()

    def draw_SAL(self,df, ax, title):

        labels = ["Structure", "Amplitude", "Location"]
        x = np.arange(len(labels)) * 2
        width = 0.3

        colors = ['red', 'green', 'blue', 'orange', 'cyan']
        # 画柱状图

        rects1 = ax.bar(x - width * 2, df['YSU'][0:3], width, label='YSU', color=colors[0])
        rects2 = ax.bar(x - width, df['QNSE'][0:3], width, label='QNSE', color=colors[1])
        rects3 = ax.bar(x, df['QNSE_EDMF'][0:3], width, label='QNSE_EDMF', color=colors[2])
        rects4 = ax.bar(x + width, df['TEMF'][0:3], width, label='TEMF', color=colors[3])
        rects5 = ax.bar(x + width*2, df['MYJ'][0:3], width, label='MYJ', color=colors[4])

        # ax.set_ylabel('SAL')
        ax.set_xticks(x)
        ax.xaxis.set_tick_params(labelsize=16)
        ax.set_yticks(np.arange(-1.5, 1.51, 0.5))  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax.set_ylim(-1.5, 1.5)  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax.yaxis.set_tick_params(labelsize=16)
        ax.set_xticklabels(labels)
        ax.set_title(title, fontsize=20)
        ax.axhline(y=0, color='black') # 画0线

    def SAL_combine(self, flag_list, area_dic):
        df_list = []
        title_list = []
        for key in area_dic:
            for flag in flag_list:
                pass
                area = area_dic[key]
                # print(area)
                ca = Caculate(flag, area)
                df = ca.get_space_scale()
                df_list.append(df)
                # title = str(flag)+'_'+str(key)
                title = str(key)+'_'+str(flag)
                title_list.append(title)

        fig = plt.figure(figsize=(12, 15), dpi=400)  # 创建页面
        grid = plt.GridSpec(4,
                            3,
                            figure=fig,
                            left=0.07,
                            right=0.96,
                            bottom=0.1,
                            top=0.96,
                            wspace=0.5,
                            hspace=0.3)

        axes = [None] * 12  # 设置一个维度为8的空列表
        axes[0] = fig.add_subplot(grid[0, 0:1])
        axes[1] = fig.add_subplot(grid[0, 1:2])
        axes[2] = fig.add_subplot(grid[0, 2:3])
        axes[3] = fig.add_subplot(grid[1, 0:1])
        axes[4] = fig.add_subplot(grid[1, 1:2])
        axes[5] = fig.add_subplot(grid[1, 2:3])
        axes[6] = fig.add_subplot(grid[2, 0:1])
        axes[7] = fig.add_subplot(grid[2, 1:2])
        axes[8] = fig.add_subplot(grid[2, 2:3])
        axes[9] = fig.add_subplot(grid[3, 0:1])
        axes[10] = fig.add_subplot(grid[3, 1:2])
        axes[11] = fig.add_subplot(grid[3, 2:3])


        for i in range(12):
            self.draw_SAL(df_list[i], axes[i], title_list[i])
            print(title_list[i])

        axes[10].legend(ncol=5 ,bbox_to_anchor=(0.5,-0.5) ,loc='lower center',fontsize=16,edgecolor='white')
        flnm = '/home/fengxiang/Project/Asses_PBL/Draw/Rain/SAL_201607.png'
        fig.savefig(flnm)


    def draw_performance(self, flag, area, ax, title):
        # pass
        # flag = flag_list[0]
        # area = area_dic['all']
        rain = get_rain_total(flag, area)
        rain_obs = rain['obs'].values
        ## 聚合模式降水数据成dataset
        model_list = ['MYJ', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        rain_list = []
        for i in model_list:
            precip = rain[i]
            rain_list.append(precip)
        rain_model = xr.concat(rain_list, dim='model')
        rain_model = rain_model.values
        # print(da)

        mem.performance(
            rain_obs, 
            rain_model,
            member_list=model_list,
            grade_list=[0.0, 0.1, 5, 10, 20, 50, 100],
            # grade_list=[0.1, 0.5, 1, 2, 5],
            # x_y='far_mr',
            # title='综合评分表现图',
            sup_fontsize =15,
            title=title,
            save_path='/home/fengxiang/Project/Asses_PBL/Draw/Rain/performance.png',
            ax = ax,
        )

    def add_legend(self, ax1):
        pass
        model_list = ['MYJ', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        grade_list=[0.0, 0.1, 5, 10, 20, 50, 100]

        legend_num = len(model_list)

        label = []
        label.extend(model_list)
        sup_fontsize =16
        # for line in range(legend_num):
        #     for i in range(len(grade_list)):
        #         ax1.plot(new_sr[line, i], new_pod[line, i], marker[i],label = i*line, color=colors[line], markersize=6)
        #         a_list.append(i*line)
        lines,label1 = ax1.get_legend_handles_labels()
        legend2 = ax1.legend(lines[0:len(lines):len(grade_list)],label,loc="lower left",
                            bbox_to_anchor=(0.9, -0.5),ncol=3, fontsize=sup_fontsize * 0.9)
        legend1=ax1.legend(lines[:len(grade_list)],['grade:'+str(i)for i in grade_list],loc="lower right",
                        bbox_to_anchor=(0.9, -0.5), ncol=4, fontsize=sup_fontsize * 0.9,
                        # ncol=5 ,bbox_to_anchor=(0.5,-0.5) ,loc='lower center',fontsize=16,edgecolor='white'
                        )
        ax1.add_artist(legend1)
        ax1.add_artist(legend2)

    def combine_performance(self, flag_list, area_dic):
        pass
        flag_l = []
        area_l = []
        title_l = []
        for key in area_dic:
            for flag in flag_list:
                pass
                area = area_dic[key]
                area_l.append(area)
                flag_l.append(flag)
                title = str(key)+'_'+str(flag)
                title_l.append(title)

        print(flag_l)
        fig = plt.figure(figsize=(12, 16), dpi=400)  # 创建页面
        grid = plt.GridSpec(4,
                            3,
                            figure=fig,
                            left=0.07,
                            right=0.96,
                            bottom=0.1,
                            top=0.96,
                            # wspace=0.2,
                            # hspace=0.3,
                            wspace=0.3,
                            hspace=0.3,
                            )

        axes = [None] * 12  # 设置一个维度为8的空列表
        axes[0] = fig.add_subplot(grid[0, 0:1])
        axes[1] = fig.add_subplot(grid[0, 1:2])
        axes[2] = fig.add_subplot(grid[0, 2:3])
        axes[3] = fig.add_subplot(grid[1, 0:1])
        axes[4] = fig.add_subplot(grid[1, 1:2])
        axes[5] = fig.add_subplot(grid[1, 2:3])
        axes[6] = fig.add_subplot(grid[2, 0:1])
        axes[7] = fig.add_subplot(grid[2, 1:2])
        axes[8] = fig.add_subplot(grid[2, 2:3])
        axes[9] = fig.add_subplot(grid[3, 0:1])
        axes[10] = fig.add_subplot(grid[3, 1:2])
        axes[11] = fig.add_subplot(grid[3, 2:3])


        for i in range(12):
            self.draw_performance(flag_l[i], area_l[i], axes[i], title_l[i])
            # print(title_list[i])


        font = 15
        axes[0].set_ylabel("Hit ratio", fontsize= font)
        axes[3].set_ylabel("Hit ratio", fontsize= font)
        axes[6].set_ylabel("Hit ratio", fontsize= font)
        axes[9].set_ylabel("Hit ratio", fontsize= font)
        axes[9].set_xlabel("Success ratio", fontsize= font)
        axes[10].set_xlabel("Success ratio", fontsize= font)
        axes[11].set_xlabel("Success ratio", fontsize= font)
        # axes[10].legend(ncol=5 ,bbox_to_anchor=(0.5,-0.5) ,loc='lower center',fontsize=16,edgecolor='white')
        self.add_legend(axes[10])
        flnm = '/home/fengxiang/Project/Asses_PBL/Draw/Rain/performance_201607.png'
        fig.savefig(flnm)

    


if __name__ == '__main__':
    area0 = {"lat1":24.875, "lat2":46.125, "lon1":70.875, "lon2":105.125}
    area1 = {"lat1":33.5, "lat2":40, "lon1":70, "lon2":105}  # north
    area2 = {"lat1":28, "lat2":33, "lon1":83, "lon2":94}  # south left 
    area3 = {"lat1":26, "lat2":33, "lon1":95, "lon2":103}  # south right
    # area4 = {"lat1":27, "lat2":28, "lon1":86, "lon2":92}  # bottom
    
    ## 选择时间范围和计算区域
    flag = 'all'
    # flag = 'day'
    # flag = 'night'
    # area = area4

    
    
    area_dic = {}
    area_dic['all'] = area0
    area_dic['north'] = area1
    area_dic['south left'] = area2
    area_dic['south right'] = area3
    month = 'May'    
    gd = GetData(month)
    rain = gd.get_rain_hourly()
    print(rain)
    area = area_dic['all']

    tr = TransferData(rain, area, flag)
    rain1 = tr.rain_hourly()
    print(rain1)
    

    #### 计算评分
    # ca = Caculate(rain)
    # ca.get_two_scale()    
    # ca.get_time_scale()    
    # ca.get_space_scale()    

    #### 画图
    # Dr = Draw()
    # Dr.draw_time_sequence()
    # Dr.SAL_combine(flag_list, area_dic)
    # Dr.combine_performance(flag_list, area_dic)





    


