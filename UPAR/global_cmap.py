
# %%
import matplotlib as mpl
import cmaps
import numpy as np
import meteva.base as mb
import matplotlib.pyplot as plt


# %%
def get_cmap_rain():
    ccc = cmaps.precip3_16lev_r
    ccc,clev = mb.def_cmap_clevs(mb.cmaps.rain_1h)
    colors = mpl.cm.get_cmap(ccc)
    col = colors(np.linspace(0, 1, 18))
    # print(col)
    cccc = mpl.colors.ListedColormap([
        # col[0],
        # col[1],
        # col[2],
        # col[3], 
        # col[4], 
        # # col[5], 
        # col[6], 
        # col[7], 
        # col[8], 
        # col[9], 
        # col[10], 
        # col[11], 
        # col[12], 
        # col[13], 
        # col[14], 
        # (231 / 250, 177 / 250, 22 / 250),
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
        # col[4],
        # col[6],
        # '#85f485',
        # '#16c516',
        # 'white',
    ])
    cmap = cccc
    return cmap

def get_cmap_q():
    ccc = cmaps.precip3_16lev_r
    ccc,clev = mb.def_cmap_clevs(mb.cmaps.rain_1h)
    colors = mpl.cm.get_cmap(ccc)
    col = colors(np.linspace(0, 1, 18))
    # print(col)
    cccc = mpl.colors.ListedColormap([
        (151/250, 232 / 250, 173 / 250),  # RGB withn 0-1 range
        (153/250, 210 / 250, 202 / 250),  # RGB withn 0-1 range
        (155/250, 188 / 250, 232 / 250),  # RGB withn 0-1 range
        (107/250, 157 / 250, 225 / 250),  # RGB withn 0-1 range
        (59/250, 126 / 250, 219 / 250),  # RGB withn 0-1 range
        (43/250, 92 / 250, 194 / 250),  # RGB withn 0-1 range
        (28/250, 59 / 250, 169 / 250),  # RGB withn 0-1 range
        (17/250, 44 / 250, 144 / 250),  # RGB withn 0-1 range
        (7/250, 30 / 250, 120 / 250),  # RGB withn 0-1 range
        (0/250, 15 / 250, 80 / 250),  # RGB withn 0-1 range
    ])
    cmap = cccc
    return cmap
def get_cmap_temp():
    # ccc = mb.cmaps.rain_1h
    # ccc,clev = mb.def_cmap_clevs(mb.cmaps.temp_2m, vmin=-30, vmax=35)
    # print(clev)

    # cccc, clev = mb.def_cmap_clevs(cmap=ccc, clevs=np.arange(300))

    ccc = cmaps.GMT_panoply
    # cmap = cmaps.temp_19lev
    # ccc = cmaps.circular_1
    # ccc = cmaps.precip3_16lev
    # print(ccc)
    colors = mpl.cm.get_cmap(ccc)
    col = colors(np.linspace(0, 1, 15))
    # # print(col)
    # col = colors(np.arange(0,1,11))
    # # # print(col)
    cccc = mpl.colors.ListedColormap([
        col[0],
        col[1],
        col[2],
        col[3],
        col[4],
        col[5],
        # col[6],
        # col[7],
        # 'white',
        col[8],
        col[9],
        col[10],
        col[11],
        col[12],
        col[13],
        col[14],
    ])

    # ## 在原来的基础上，多搞几个细分的颜色
    cccc, clev = mb.def_cmap_clevs(cmap=cccc, clevs=np.arange(40))

    # # # print(len(cmap.values))
    cmap = cccc
    return cmap

def get_cmap_landuse():
    pass
    # cccc = mpl.colors.ListedColormap([
    #             ('#FF0000'),   
    #             ('#FF4500'),   
    #             ('#FF8C00'),   
    # ])
    # cmap = cccc
    cmap = cmaps.vegetation_modis
    return cmap

def draw():
    fig, ax = plt.subplots(figsize=(6, 1))
    fig.subplots_adjust(bottom=0.5)
    # cmap = get_cmap_rain()
    # cmap = get_cmap_rain()
    # cmap = get_cmap_q()
    cmap = get_cmap_landuse()
    # norm = mpl.colors.Normalize(vmin=5, vmax=10)
    fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap),
             cax=ax, orientation='horizontal', label='Some Units')
    fig.savefig('test')

        
    # aa = get_cmap_rain()
    # print(aa)
draw()
# %%
