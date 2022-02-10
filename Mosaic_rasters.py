#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 12:43:08 2022

@author: jianmin.wang
"""


#### using mosaic_overlapp
###################           HLS RASTER MOSAIC                ##############################
#'/gpfs/data/xyz/jianmin/Rodman/spatial_data/Spring/Spring.tif', 
#import os
#import re
from osgeo import gdal
#from osgeo import osr
import numpy as np
import sys
import math
sys.path.append("/home/jianmin.wang/codes/self_packages/python3")
# from file import *
# from JWraster import *
import JWraster
dir(JWraster)


FILL=32767
pathin = '/gpfs/data/xyz/jianmin/HLSdata/phenology_nofusion/'
pathout = '/gpfs/data/xyz/jianmin/HLSdata/mosaics/'
phens = ['gri', 'giMD', 'gre', 'sei', 'seMD', 'see', 
         'griVI', 'giMDVI', 'greVI', 'seiVI', 'seMDVI', 'seeVI', 
         'minVI', 'maxVI', 'grprate', 'setrate', 'growlength', 'growarea']
nbs = [3]*len(phens)
nb = 3
tiles = ['T13SEA', 'T13SEB', 'T13SDA', 'T13SDB']
years = ['2018']
#years = ['2016', '2017', '2018']

for phen in phens: 
    for year in years:
        inputs = ['%sz1%s_HLS_%s_%s' % (pathin, phen, tile, year) for tile in tiles]
        output = '%sz1%s_HLS_mosaic_%s' % (pathout, phen, year) 
        JWrasetr.mosaic_overlapp(inputs, output, nb, "ENVI")
        outF = open(output+'.hdr', 'a')
        outF.write("data ignore value = %s" % FILL)
        outF.close()
        os.remove(output+".aux.xml")
        print("phenology metric %s in year %s was done!" % (phen, year))
####################           HLS RASTER MOSAIC                ##############################   



#%% using mosaic_gdal
#### MOSAIC NAIP IMAGES from a text file#####
import sys
sys.path.append("/home/jianmin.wang/codes/self_packages/python3")
import JWraster
dir(JWraster)
pathin = '/gpfs/data/xyz/jianmin/NAIPdata/Ponil_Complex/'
pathout = '/gpfs/data/xyz/jianmin/NAIPdata/Ponil_Complex/'
file = open(pathin+'file_list.txt', 'r')
files = file.readlines()
file.close()
files = [pathin+file.rstrip("\n") for file in files]
mosaic_gdal(files, pathout+'Ponil_Complex_2018.tif', 255, 'GTiff')

