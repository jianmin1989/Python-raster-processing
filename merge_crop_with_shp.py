#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 12:54:37 2022

@author: jianmin.wang
"""

###merge rasters   
from osgeo import gdal, osr, ogr
area = 'large'
pathin = '/scratch/jianmin.wang/HLS_VIIRS_Kyle/%s/' % area
os.chdir(pathin)


# driver = ogr.GetDriverByName('ESRI Shapefile')
# dataset = driver.Open(shapefile)
# spatialRef = dataset.GetLayer().GetSpatialRef()
# spatialRef.ExportToWkt()
# dataset =None

pathout = '/scratch/jianmin.wang/HLS_VIIRS_Kyle/'
phens_out = ['Greenup', 'Senesce', 'Peak', 'MaxVI', 'AmpVI', 'VISD']
phens_out = ['VISD']
# phens_out = ['Greenup', 'Senesce', 'Peak', 'MaxVI', 'VISD']
# phen = 'VISD'
for phen in phens_out:
    print("generating phenology %s" % phen)
    if phen == 'VISD':
        fillvalue = -9999
    else:
        fillvalue = 32767
    files = ["%s_VIIRS_HLS_B4_%s_%s.tif" % (phen, tile, year) for tile in tiles]        
    # print(files)       
    outfile =   "mosaic_%s_VIIRS_HLS_B4_%s_%s.tif" % (phen, area, year)           
    # mosaic_gdal_merge(files, outfile, fillvalue, 'GTiff')   

    # cropfile= os.path.join(pathout, "%s_VIIRS_HLS_B4_%s_%s.tif" % (phen, area, year)) 
    ds_out = gdal.Warp(outfile, files, format='GTiff',  errorThreshold=0,
                            xRes= 30, yRes =30, #outputBounds=outputBounds,
    #                            coordinateOperation = srs,
                            dstSRS = 'EPSG:32613', targetAlignedPixels=True, #srs, #, transfomerOptions=True, targetAlignedPixels=True, 
                            resampleAlg=gdal.GRA_NearestNeighbour, 
                            #outputType = gdal.GDT_Byte, workingType = gdal.GDT_Byte,
                            # cutlineDSName = shapefile, 
                            # cutlineLayer = shape, 
                            # cropToCutline=True,
                            srcNodata=fillvalue, dstNodata=fillvalue)        
    ds_out = None        




# shape_ds = ogr.Open(shapefile)
# shape = shape_ds.GetLayer()
for phen in phens_out:
    print("generating phenology %s" % phen)
    if phen == 'VISD':
        fillvalue = -9999
    else:
        fillvalue = 32767
        
    outfile =   "mosaic_%s_VIIRS_HLS_B4_%s_%s.tif" % (phen, area, year)   
    cropfile= os.path.join(pathout, "%s_VIIRS_HLS_B4_%s_%s.tif" % (phen, area, year)) 
    ds_out = gdal.Warp(cropfile, outfile, format='GTiff',  errorThreshold=0,
                            #xRes= resox, yRes =resoy, outputBounds=outputBounds,
    #                            coordinateOperation = srs,
                            #dstSRS = proj_code, #srs, #, transfomerOptions=True, targetAlignedPixels=True, 
                            resampleAlg=gdal.GRA_NearestNeighbour, 
                            #outputType = gdal.GDT_Byte, workingType = gdal.GDT_Byte,
                            cutlineDSName = shapefile, 
                            # cutlineLayer = shape, 
                            cropToCutline=True,
                            srcNodata=fillvalue, dstNodata=fillvalue)        
    ds_out = None
