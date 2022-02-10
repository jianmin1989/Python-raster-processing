#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 12:46:20 2022

@author: jianmin.wang
"""

###Convert hdf file to geotiff file using gdal.translate
import re, os
from osgeo import gdal
file =  '/scratch/jianmin.wang/data/LST_MODIS/MOD11A1.A2017365.h11v04.006.2018002063958.hdf'
fileout = '/gpfs/scratch/jianmin.wang/ARD_MODIS_LST/test_h11v04'

info = gdal.Info(file)
## SELECT A BAND TO DO THE TRANSLATE
## Here by defaul use the first one

band1 = re.search('SUBDATASET_1_NAME=(.*)\n', info).group(1) #group number is the order of () in the pattern
tmp = gdal.Translate(fileout, band1, format='ENVI')
tmp = None

##ONLY IF THE OUTPUT IF ENVI AND THE INPUT HDF IS SINUSOIDAL PROJECTION.
if os.path.isfile(fileout+'.hdr'):
    outF = open(fileout+'.hdr', 'r')
    lines = outF.readlines()
    outF.close()
       
    new_file = open(fileout+'.hdr', 'w+')
    for line in lines:
        if re.match("coordinate system string = {", line, re.IGNORECASE):   
            new_file.write('projection info = {16, 6371007.181, 0.000000, 0.0, 0.0, D_Unknown, Sinusoidal, units=Meters}\n')
            new_file.write('coordinate system string = {PROJCS["Sinusoidal",GEOGCS["GCS_Unknown",DATUM["D_Unknown",SPHEROID["S_Unknown",6371007.181,0.0]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],UNIT["Meter",1.0]]}\n')
        else:
            new_file.write(line)
    #new_file.write("data ignore value = %s" % -1.0)
    new_file.close()
###Convert hdf file to geotiff file to get the projection information  




## convert hdf4 to GTiff with a specified projection. (convert and reproject)
import re
from osgeo import gdal


## Parameters
file =  '/scratch/jianmin.wang/data/LST_MODIS/MOD11A1.A2017365.h11v04.006.2018002063958.hdf'
fileout = '/gpfs/scratch/jianmin.wang/temp/test_h11v04.tif'
band = 1
reso = 0.01
dstNodata = 0
resampleAlg = gdal.GRA_NearestNeighbour
dstSRS = 'EPSG:4326'


#### process
info = gdal.Info(file)
band = re.search('SUBDATASET_%s_NAME=(.*)\n' % band, info).group(1) #group number is the order of () in the pattern
op = gdal.WarpOptions(format='GTiff', errorThreshold=0,
                            xRes= reso, yRes =reso, 
                            targetAlignedPixels=True, dstSRS = dstSRS, 
                            resampleAlg=resampleAlg, dstNodata=dstNodata)
ds = gdal.Warp(fileout, band, options=op)
ds = None



########################## Convert HDF5 to GeoTiff
import sys
sys.path.append("/home/jianmin.wang/codes/self_packages/python3")
# from file import *
# from JWraster import *
import JWraster
dir(JWraster)

outdir = "/gpfs/scratch/jianmin.wang/temp/"
file = "/scratch/yongchang.ye/product/viirs/VNP22Q2/2019/001/VNP22Q2.A2019001.h11v04.001.2020225123025.h5"
# file = "/gpfs/scratch/jianmin.wang/temp/VNP21A1D.A2017121.h33v07.001.2019170075505.h5"


# proj = JWraster.get_sinu_projection()
h5 = JWraster.HDF5(file)
gt = h5.get_geotransformation()
layernames = h5.get_layernames()
layers = range(len(layernames))
for layer in layers:
    obj = h5.open_layer(layer)
    # get layer attributes
    attrs = h5.LayerAttrs(obj)

    # get filled value
    fill_value = attrs.get_fill_value()
    print("layer %d fill %d " % (layer, fill_value ))
    print("layer %d - %s" % (layer, layernames[layer]))
    h5.write_to_geotiff(outdir, layer, proj, gt)