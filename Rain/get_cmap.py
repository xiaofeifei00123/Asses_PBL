#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
标准化设置色标
-----------------------------------------
Time             :2021/07/19 15:36:15
Author           :Forxd
Version          :1.0
'''

import os
import sys
import matplotlib as mpl
import cmaps
import numpy as np
import meteva.base as mb
import matplotlib.pyplot as plt


def get_cmap_rain():
    ccc = cmaps.precip3_16lev_r
    ccc, clev = mb.def_cmap_clevs(mb.cmaps.rain_1h)
    colors = mpl.cm.get_cmap(ccc)
    col = colors(np.linspace(0, 1, 18))
    cccc = mpl.colors.ListedColormap([
        'white',
        (165/250, 243 / 250, 141 / 250),  # RGB withn 0-1 range
        (153/250, 210 / 250, 202 / 250),  # RGB withn 0-1 range
        (155/250, 188 / 250, 232 / 250),  # RGB withn 0-1 range
        (107/250, 157 / 250, 225 / 250),  # RGB withn 0-1 range
        (59/250, 126 / 250, 219 / 250),  # RGB withn 0-1 range
        (43/250, 92 / 250, 194 / 250),  # RGB withn 0-1 range
        (28/250, 59 / 250, 169 / 250),  # RGB withn 0-1 range
        (17/250, 44 / 250, 144 / 250),  # RGB withn 0-1 range
        (7/250, 30 / 250, 120 / 250),  # RGB withn 0-1 range
        (70/250, 25 / 250, 129 / 250),  # RGB withn 0-1 range
        (134/250, 21 / 250, 139 / 250),  # RGB withn 0-1 range
        (200/250, 17 / 250, 169 / 250),  # RGB withn 0-1 range
        (129/250, 0 / 250, 64 / 250),  # RGB withn 0-1 range
    ])
    cmap = cccc
    return cmap


def get_temp():
    cccc = mpl.colors.ListedColormap([
        (2 / 250, 12 / 250,  100/250),
        (7 / 250, 30 / 250,  120/250),
        (17 / 250, 49 / 250,  139/250),
        (27 / 250, 68 / 250,  159/250),
        (38 / 250, 87 / 250,  179/250),
        (48 / 250, 106 / 250,  199/250),
        (59 / 250, 126 / 250,  219/250),
        (78 / 250, 138 / 250,  221/250),
        (97 / 250, 150 / 250,  224/250),
        (116 / 250, 163 / 250,  226/250),
        (135 / 250, 175 / 250,  229/250),
        (155 / 250, 188 / 250,  232/250),
        (154 / 250, 196 / 250,  220/250),
        (153 / 250, 205 / 250,  208/250),
        (152 / 250, 214 / 250,  196/250),
        (151 / 250, 232 / 250,  173/250),
        (215 / 250, 222 / 250,  126/250),
        (234 / 250, 219 / 250,  112/250),
        (244 / 250, 217 / 250,  99 / 250),
        (250 / 250, 204 / 250,  79 / 250),
        (247 / 250, 180 / 250,  45 / 250),
        (242 / 250, 155 / 250,  0 / 250),
        (241 / 250, 147 / 250,  3 / 250),
        (240 / 250, 132 / 250,  10 / 250),
        (239 / 250, 117 / 250,  17 / 250),
        (238 / 250, 102 / 250,  24 / 250),
        (238 / 250, 88 / 250,  31 / 250),
        (231 / 250, 75 / 250,  26 / 250),
        (224 / 250, 63 / 250,  22 / 250),
        (217 / 250, 51 / 250,  18 / 250),
        (208 / 250, 36 / 250,  14 / 250),
        (194 / 250, 0 / 250,  3 / 250),
        (181 / 250, 1 / 250,  9 / 250),
        (169 / 250, 2 / 250,  16 / 250),
        (138 / 250, 5 / 250,  25 / 250),
        (111 / 250, 0 / 250,  21 / 250),
        (80 / 250, 0 / 250,  15 / 250),
    ])

    temp_level = [
        -30, -28, -26, -24, -22, -20, -18, -16, -14, -12, -10,
        -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
        22, 24, 26, 28, 30, 32, 34, 35, 37, 40,
    ]

    cmap = cccc
    return cmap, temp_level


def get_cmap_temp():
    # ccc = mb.cmaps.rain_1h
    ccc, clev = mb.def_cmap_clevs(mb.cmaps.temp_2m)
    # cmap = ccc
    cmap = cmaps.temp_19lev
    return cmap


def draw():
    fig, ax = plt.subplots(figsize=(6, 1))
    fig.subplots_adjust(bottom=0.5)
    cmap, clevel = get_temp()
    # print(clevel)
    # cmap = get_cmap_rain()
    # cmap = get_cmap_temp()
    # print(type(cmap))
    # norm = mpl.colors.Normalize(vmin=5, vmax=10)
    # norm = mpl.colors.Normalize(clevel)
    bounds = clevel  # 分隔点的值
    norm = mpl.colors.BoundaryNorm(
        boundaries=bounds, ncolors=len(bounds)+1, extend='both')

    # print(norm)
    ticks = [-10, 0, 2,10]
    fig.colorbar(
        mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
        # cax=ax,   # 指定colorbar是画在哪个ax内
        ax = ax,  # 指定ax是哪几个子图的
        orientation='horizontal',
        label='Test colorbar',
        ticks=ticks,
        spacing='uniform',
        # fraction=100,
        # shrink=1.0,
        aspect=30,  # 长短边比例
        pad=10  # 离图片的距离
        )
    path = sys.path[0]
    fig_name = os.path.join(path, 'cmap_test')
    fig.savefig(fig_name)


if __name__ == '__main__':

    # aa = get_cmap_rain()
    # print(aa)
    draw()
