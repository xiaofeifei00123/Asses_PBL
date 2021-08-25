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
