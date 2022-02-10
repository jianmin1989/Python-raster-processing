#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 12:39:49 2022

@author: jianmin.wang
"""

### REPROJECT MULTIPLE FILES TO A SINGLE FILE TO MATCH ANOTHER FILE ######
from osgeo import gdal 
FILL = 32767
path = "/gpfs/scratch/jianmin.wang/ARD_MODIS_LST/composite_MODIS/"
tiles = ["h10v04", "h11v04"]
year = 2017
dname = 'LST' #LSTQA LST
infiles = ["%sMODIS_%s.%s.%s.year1.BIP" % (path, dname, year, tile) for tile in tiles]

fileout = "%sMODIS_%s.%s.ARD016006.tif" % (path, dname, year)


nrow, ncol, minx, miny, maxx, maxy, resox, resoy, srs = get_extent("/gpfs/scratch/jianmin.wang/ARD_MODIS_LST/fused_fulltile_nokshift/z1LSTfuse_016006_2016_tif.tif")
ds =  gdal.Warp(fileout, infiles, format='GTiff', outputBounds=[minx, miny, maxx, maxy], errorThreshold=0,
                            xRes= 960, yRes =960, dstSRS = srs, #, transfomerOptions=True, targetAlignedPixels=True, 
                            outputType = gdal.GDT_Int16, workingType = gdal.GDT_Int16,
                            resampleAlg=gdal.GRA_NearestNeighbour, srcNodata=FILL, dstNodata=FILL)  #gdal.GRA_Bilinear
ds = None
### REPROJECT MULTIPLE FILES TO A SINGLE FILE TO MATCH ANOTHER FILE ######


