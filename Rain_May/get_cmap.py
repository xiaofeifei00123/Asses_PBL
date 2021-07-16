
import matplotlib as mpl
import cmaps
import numpy as np
import meteva.base as mb
import matplotlib.pyplot as plt

def get_cmap_rain():
    ccc = cmaps.precip3_16lev_r
    ccc,clev = mb.def_cmap_clevs(mb.cmaps.rain_1h)
    colors = mpl.cm.get_cmap(ccc)
    col = colors(np.linspace(0, 1, 18))
    # print(col)
    # cccc = mpl.colors.ListedColormap([
    #     col[0],
    #     col[1],
    #     col[2],
    #     col[3], 
    #     col[4], 
    #     # col[5], 
    #     col[6], 
    #     col[7], 
    #     col[8], 
    #     col[9], 
    #     col[10], 
    #     col[11], 
    #     # col[12], 
    #     col[13], 
    #     # col[14], 
    #     # (231 / 250, 177 / 250, 22 / 250),
    #     # col[4],
    #     # col[6],
    #     # '#85f485',
    #     # '#16c516',
    #     # 'white',
    # ])
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

def get_cmap_temp():
    # ccc = mb.cmaps.rain_1h
    ccc,clev = mb.def_cmap_clevs(mb.cmaps.temp_2m)
    # cmap = ccc
    cmap = cmaps.temp_19lev
    return cmap


def draw():
    fig, ax = plt.subplots(figsize=(6, 1))
    fig.subplots_adjust(bottom=0.5)
    cmap = get_cmap_rain()
    # cmap = get_cmap_temp()
    # print(type(cmap))
    # norm = mpl.colors.Normalize(vmin=5, vmax=10)
    fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap),
             cax=ax, orientation='horizontal', label='Some Units')
    fig.savefig('test')

if __name__ == '__main__':
        
    # aa = get_cmap_rain()
    # print(aa)
    draw()