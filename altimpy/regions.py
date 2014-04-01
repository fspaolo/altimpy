"""
Module with definition of ice-shelf "boxes" and other regions.

"""
# general regions
ais = (0, 360, -90, -60),
wais = (179, 311.5, -90, -60)
eais1 = (0, 179, -85, -60)         # lon is 0/360!
eais2 = (311.5, 360, -85, -60)
eais1_no_rf = (0, 179, -77, -60)   # w/o Ross and Filchner 
eais2_no_rf = (330, 360, -85, -60)
antpen = (276, 310, -74.5, -64)
amundsen = (220, 276, -76, -71)

# individual ice shelves (clock-wise starting at lon=0)
fimbulw = (349, 360, -72, -68)     # revise boundaries!
fimbule = (0, 13, -72, -68)        # revise boundaries!
lazarev = (13.2, 34.5, -72, -68)   # revise boundaries!
amery = (62, 81, -74, -67)
west = (80, 90, -68, -65)
shackleton = (94.5, 101.5, -67, -64)
conger = (101.5, 105.5, -67, -65)
totten = (112, 117.5, -68, -66)
moscow = (117.5, 123, -68, -66)
ross = (150, 225, -82, -77)
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
larsenb = (294.5, 302, -66.5, -65)
larsenc = (292, 302, -69.2, -66.2)
larsend = (295.5, 302, -73, -69.2)
ronne = (275, 312, -82, -74.2)
filchner = (312, 332, -82, -77)
brunt = (331, 340, -76.5, -73.2)
riiser = (339, 350, -74.5, -71.5)
