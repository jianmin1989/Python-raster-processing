#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 12:54:01 2022

@author: jianmin.wang
"""

###Quick reprojection NLCD
infile = '/hunter/data1/wangj/NLCD/nlcd_2011_landcover_2011_edition_2014_10_10'
infile = '/hunter/data1/wangj/NLCD/NLCD_2016_Land_Cover_L48_20190424.img'
outfile = '/hunter/data1/wangj/LSP_30m/nlcd/nlcd_18TXM_2016.tif'  #_UTM18N wgs84 
infile = '/hunter/data1/wangj/NLCD/NLCD_2013_Land_Cover_L48_20190424.img' 
outfile = '/hunter/data1/wangj/LSP_30m/nlcd/nlcd_18TXM_2013'  #_UTM18N wgs84 

proj_code = 'EPSG:32618'
resox =30 
resoy = 30
ns = 3660
nl =3660
minx = 600000
maxy = 4700040
miny = maxy-nl*resoy
maxx = minx + ns*resox
outputBounds = [minx, miny, maxx, maxy] #[minx, miny, maxx, maxy]
tmp =  gdal.Warp(outfile, infile, format='ENVI',  errorThreshold=0,
                        xRes= resox, yRes =resoy, outputBounds=outputBounds,
#                            coordinateOperation = srs,
                        dstSRS = proj_code, #srs, #, transfomerOptions=True, targetAlignedPixels=True, 
                        resampleAlg=gdal.GRA_NearestNeighbour, 
                        outputType = gdal.GDT_Byte, workingType = gdal.GDT_Byte,
                        srcNodata=0, dstNodata=0)
tmp = None

