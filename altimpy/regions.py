"""
Module with definition of ice-shelf "boxes" and other regions.

WARNING: if a different grid with more cells is used, the boundaries
should be reviewed.

"""

import numpy as np

def where_isnan(name, lon, lat):

    if lon.ndim == 1:
        lon, lat = np.meshgrid(lon, lat)
    reg = {}

    # regional ice shelves
    #---------------------

    reg['qml'] = reg['Queen Maud'] = \
            (lon > 62) & (lon < 330) | (lat < -76.5) | (lat > -66)

    reg['ws'] = reg['West-Shackleton'] = \
            (lon < 80) | (lon > 105) | (lat < -68) | (lat > -64.5)

    reg['wvl'] = reg['Wilkes-Victoria'] = \
            (lon < 80) | (lon > 174) | (lat < -76) | (lat > -64)

    reg['ris'] = reg['Ross'] = \
            (lon < 156) | (lon > 215) | (lat < -85.8) | (lat > -77) | \
            (lon > 190) & (lat > -77.7) | \
            (lon > 208) & (lat > -78.5)

    reg['as'] = reg['Amundsen'] = \
            (lon < 202) | (lon > 262.5) | (lat < -78) | (lat > -73.2)

    reg['bs'] = reg['Bellingshausen'] = \
            (lon < 256) | (lon > 295) | (lat < -74) | (lat > -69.5) | \
            (lon < 262) & (lat < -73.2)

    reg['lis'] = reg['Larsen'] = \
            (lon < 294) | (lon > 300.5) | (lat < -73.2) | (lat > -65.3)

    reg['ap'] = reg['Antarctic Peninsula'] = \
            (lon < 280) | (lon > 305) | (lat < -74.5) | (lat > -64) | \
            (lon > 295) & (lat < -73)

    reg['fris'] = reg['Filchner-Ronne'] = \
            (lon < 275) | (lon > 338) | (lat < -84) | (lat > -74.7) | \
            (lon > 330) & (lat > -78)

    reg['eais'] = reg['East Antarctica'] = \
            (lon > 179) & (lon < 313.5) & (lat > -80.6) | \
            (lon > 179) & (lon < 307) & (lat > -80.7) | \
            (lon > 179) & (lon < 306) & (lat > -81.5) | \
            (lon > 179) & (lon < 304.5) & (lat > -82)

    reg['wais'] = reg['West Antarctica'] = \
            (lon < 179) | (lon > 313) | (lat < -90) | (lat > -60) | \
            (lon > 306) & (lat < -80.6) | (lon > 304) & (lat < -81.5)

    reg['ais'] = reg['All Antarctica'] = \
            (lat < -90) | (lat > -60)

    # individual ice shelves
    #-----------------------

    reg['brunt'] = reg['Brunt'] = \
            (lon < 332) | (lon > 342.5) | (lat < -76.2) | (lat > -74) | \
            (lat > -74.5) & (lon > 339.5)

    reg['riiser'] = reg['Riiser'] = \
            (lon < 339.5) | (lon > 349) | (lat > -71.6) | (lat < -74.5)

    reg['fimbul'] = reg['Fimbul'] = \
            (lon > 7.5) & (lon < 357.4) | (lat > -69) | (lat < -71.8)

    reg['lazarev'] = reg['Lazarev'] = \
            (lon < 12.3) | (lon > 16) | (lat < -70.7) | (lat > -69) | \
            (lon < 13) & (lat < -70)

    reg['baudouin'] = (lon < 24) | (lon > 33.3) | (lat > -69) | (lat < -71.2)

    reg['harald'] = (lon < 33.3) | (lon > 37) | (lat > -68.5) | (lat < -70)

    reg['amery'] = reg['Amery'] = (lon < 66) | (lon > 74.5) | (lat > -68) | \
            (lat < -73.5)

    reg['west'] = reg['West'] = \
            (lon < 81) | (lon > 89.5) | (lat < -68) | (lat > -66)

    reg['shackleton'] = reg['Shackleton'] = \
            (lon < 94.8) | (lon > 103) | (lat < -67) | (lat > -64.8)

    reg['totten'] = reg['Totten'] = (lon < 113.5) | (lon > 118) | \
            (lat < -67.6) | (lat > -66.3)

    reg['moscow'] = reg['Moscow'] = (lon < 118.5) | (lon > 123) | \
            (lat < -67.5) | (lat > -66.3)

    reg['holmes'] = reg['Holmes'] = (lon < 126.3) | (lon > 128) | \
            (lat < -67.2) | (lat > -66.3)

    reg['dibble'] = reg['Dibble'] = (lon < 133.5) | (lon > 135.5) | \
            (lat < -66.6) | (lat > -65.8)

    reg['mertz'] = reg['Mertz'] = (lon < 143) | (lon > 147) | \
            (lat < -68.3) | (lat > -66.5)

    reg['cook'] = reg['Cook'] = (lon < 151) | (lon > 154) | \
            (lat < -69) | (lat > -68)

    reg['rennick'] = reg['Rennick'] = (lon < 160.7) | (lon > 162.7) | \
            (lat < -71.3) | (lat > -70)

    reg['drygalski'] = reg['Drygalski'] = (lon < 160) | (lon > 166) | \
            (lat < -75.8) | (lat > -75.2)

    reg['mariner'] = reg['Mariner'] = (lon < 166) | (lon > 170) | \
            (lat < -73.7) | (lat > -72.8)

    reg['rosse'] = reg['Ross EAIS'] = (lon < 156) | (lon > 179) | \
            (lat < -85.8) | (lat > -77.1)

    reg['rossw'] = reg['Ross WAIS'] = (lon < 179) | (lon > 215) | \
            (lat < -85.8) | (lat > -77.7) | (lon > 208) & (lat > -78.5)

    reg['sulzberger'] = reg['Sulzberger'] = (lon < 207.5) | (lon > 215) | \
            (lat < -78) | (lat > -76.3) | (lon > 211.2) & (lat > -76.4)

    reg['nickerson'] = reg['Nickerson'] = (lon < 210.5) | (lon > 218) | \
            (lat < -76.4) | (lat > -75.2)

    reg['getz'] = reg['Getz'] = (lon < 225) | (lon > 245.5) | \
            (lat < -75.3) | (lat > -73.2)

    reg['dotson'] = reg['Dotson'] = (lon < 245.2) | (lon > 248.7) | \
            (lat < -75.3) | (lat > -73.7) | (lon > 247.5) & (lat < -75)

    reg['crosson'] = reg['Crosson'] = (lon < 249.3) | (lon > 251) | \
            (lat < -75.5) | (lat > -74.5) | (lon < 249) & (lat > -75)

    reg['thwaites'] = reg['Thwaites'] = (lon < 251.3) | (lon > 257) | \
            (lat < -75.5) | (lat > -74.5)

    reg['pig'] = reg['Pine Island'] = (lon < 257) | (lon > 262.5) | \
            (lat < -75.8) | (lat > -74)

    reg['cosgrove'] = reg['Cosgrove'] = (lon < 257) | (lon > 262) | \
            (lat < -74) | (lat > -73.2)

    reg['abbot'] = reg['Abbot'] = (lon < 256.2) | (lon > 269) | \
            (lat < -73.5) | (lat > -72) | (lon < 261) & (lat < -73.3)

    reg['venable'] = reg['Venable'] = (lon < 271.2) | (lon > 275) | \
            (lat < -73.5) | (lat > -72.7)

    reg['stange'] = reg['Stange'] = (lon < 281) | (lon > 285) | \
            (lat < -74) | (lat > -72.5)

    reg['bach'] = reg['Bach'] = (lon < 286) | (lon > 290) | \
            (lat < -72.5) | (lat > -71.6)

    reg['wilkins'] = reg['Wilkins'] = (lon < 285) | (lon > 291) | \
            (lat < -71.5) | (lat > -69.5) | (lon > 290) & (lat > -70.5)

    reg['georgevi'] = reg['George VI'] = (lon < 285.5) | (lon > 293.8) | \
            (lat < -73.8) | (lat > -70) | (lon < 290.7) & (lat > -72.5)

    reg['larsenb'] = reg['Larsen B'] = (lon < 297.5) | (lon > 300) | \
            (lat < -66.1) | (lat > -65.3)

    reg['larsenc'] = reg['Larsen C'] = (lon < 294) | (lon > 300.5) | \
            (lat < -69.5) | (lat > -66.1)

    reg['larsend'] = reg['Larsen D'] = (lon < 297) | (lon > 300.5) | \
            (lat < -73) | (lat > -69.5)

    reg['ronne'] = reg['Ronne'] = (lon < 275) | (lon > 313) | (lat < -84) | \
            (lat > -74.7) | (lon > 306) & (lat < -80.6) | \
            (lon > 304) & (lat < -81.5)

    reg['filchner'] = reg['Filchner'] = (lon < 300) | (lon > 338) | \
            (lat < -84) | (lat > -78) | (lon < 313.5) & (lat > -80.6) | \
            (lon < 307) & (lat > -80.7) | \
            (lon < 306) & (lat > -81.5) | \
            (lon < 304.5) & (lat > -82)

    return np.where(reg[name])


# DEPRECATED Use conditions instead
# general regions
ais = allantarctica = antarctica = (0, 360, -90, -60)
wais = westantarctica = (179, 311.5, -90, -60)
eais1 = eastantarctica1 = (0, 179, -90, -60)         # lon is 0/360!
eais2 = eastantarctica2 = (311.5, 360, -90, -60)
eais1_no_rf = (0, 179, -77, -60)                     # w/o Ross and Filchner 
eais2_no_rf = (330, 360, -77, -60)
queenmaud = qml = (0, 38, -71.5, -68.5)              # FIXME, it's not complete!
fris = (272, 332, -90, -74.5)
antpen = ap = (276, 310, -74.5, -64)
larsen = easternap = lis = (294, 302, -74, -64)
belling = bellingshausen =  (260.8, 294.4, -74, -69)  # revise boundaries!
amundsen = (220, 263, -76, -73.1)
ross = ris = (155, 213.3, -85, -77.3)                 # revise boundaries!
aws = (62, 106, -74, -63)
wilkes = (112, 167, -72, -65)
westshac = (80, 105, -68, -64.5)
tottmosc = (113, 124, -67.5, -66)

# DEPRECATED Use conditions instead
# individual ice shelves (clock-wise starting at Brunt)
brunt = stancomb = (331, 340, -76.5, -73.2)
riiser = (339, 350, -74.5, -71.5)
fimbulw = (349, 360, -72, -68)     # revise boundaries!
fimbule = (0, 13, -72, -68)        # revise boundaries!
lazarev = (13.2, 34.5, -72, -68)   # revise boundaries!
amery = (62, 81, -74, -67)
west = (80, 90, -68, -65)
shackleton = (94.5, 101.5, -67, -64)
conger = (101.5, 105.5, -67, -65)
totten = (112, 117.5, -68, -66)
moscow = (117.5, 123, -68, -66)
holmes = (126.3, 128, -67.2, -66.3)
rosse = (155, 179, -85, -77)    # revise boundaries!
rossw = (178.3, 213.3, -85, -77.5)
sulzberger = (207, 216, -78, -76)
getz = (224, 245, -76, -73)
dotson = (245, 249, -76, -73.5)    # revise boundaries!
crosson = (248.5, 251, -76, -73.5) # revise boundaries!
thwaites = (251, 256.5, -76, -74)
pig = (256.5, 263, -76, -74)
abbot = (255, 272, -73.3, -71)
stange = (275, 285, -75, -72)
bach = (286, 290, -72.8, -71.5) 
wilkins = (285, 290, -71.8, -68)
georges = (284.8, 290, -74, -72.5)
georgen = (289.2, 295, -74, -69)
larsena = (297, 302, -65.5, -64)
larsenb = (295, 301.5, -66.5, -65)
larsenc = (294, 301.5, -69.2, -66.2)
larsend = (295.5, 302, -73, -69.2)
ronne = (275, 312, -82, -74.2)
filchner = (312, 332, -82, -77)


