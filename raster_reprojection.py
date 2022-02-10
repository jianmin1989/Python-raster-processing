#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 12:36:24 2022

@author: jianmin.wang
"""

#%% ####################           RASTER REPROJECTION                ##############################
#'/gpfs/data/xyz/jianmin/Rodman/spatial_data/Spring/Spring.tif', 
#import os
#import re
from osgeo import gdal
import glob
import os
#Fires: Borrego  H12  High_Park  Ponil_Complex  Saw  Spring  Waldo_Canyon  West_Fork  York
Fires=[ 'H12', 'Ponil_Complex',  'Saw',  'Spring' , 'West_Fork',  'York']
bounds = {'H12':[483900, 4087680, 488700, 4094640],
          'Ponil_Complex':[484620, 4046880, 509040, 4074570],
          'Saw':[517770, 4080750, 523530, 4085520],
          'Spring':[494220, 4090800, 510240, 4101570],
          'West_Fork':[510540, 4076340, 517920, 4084170],
          'York':[500790, 4086600, 505980, 4091130]
        }
pathin = '/gpfs/data/xyz/jianmin/Rodman/spatial_data/'
pathout = '/gpfs/data/xyz/jianmin/Rodman/'
pathmtbsin = '/gpfs/data/xyz/jianmin/Fires/MTBS/'
pathmtbsout = '/gpfs/data/xyz/jianmin/Fires/'

#[(value[3]-value[1])/30 for value in bounds.values()]


fire = 'Saw'

infiles = [glob.glob(pathin+fire+'/*_PSME_Model.tif')[0], 
           glob.glob(pathin+fire+'/*_PIPO_Model.tif')[0], 
           glob.glob(pathmtbsin + fire+ '/*_dnbr6.tif')[0]]

outfiles = [pathout + fire+'_PSME_Model_UTMWGS84', 
            pathout + fire+'_PIPO_Model_UTMWGS84', 
            pathmtbsout + fire+'_mtbs_UTMWGS84.tif']
outputBounds=bounds[fire]
#print("for %s" % fire)
#print('inputs are :' )
#print(infiles)
#print('Output are :' )
#print(outfiles)

ops = []
ops.append(gdal.WarpOptions(format='ENVI', outputBounds=outputBounds,  errorThreshold=0,
                            xRes= 30, yRes =30, targetAlignedPixels=True, dstSRS = 'EPSG:32613', #, transfomerOptions=True
                            resampleAlg=gdal.GRA_NearestNeighbour, dstNodata=-1.0))
ops.append(gdal.WarpOptions(format='ENVI', outputBounds=outputBounds, errorThreshold=0,
                            xRes= 30, yRes =30, targetAlignedPixels=True, dstSRS = 'EPSG:32613', 
                            resampleAlg=gdal.GRA_NearestNeighbour, dstNodata=-1.0))
#ops.append(gdal.WarpOptions(format='ENVI', outputBounds=outputBounds,  #errorThreshold=0,
#                            xRes= 1, yRes =1, targetAlignedPixels=True, dstSRS = 'EPSG:32613', 
#                            resampleAlg=gdal.GRA_NearestNeighbour, dstNodata=255))
ops.append(gdal.WarpOptions(format='GTiff', outputBounds=outputBounds,  errorThreshold=0,
                            xRes= 30, yRes =30, targetAlignedPixels=True, dstSRS = 'EPSG:32613', 
                            resampleAlg=gdal.GRA_NearestNeighbour, dstNodata=255))



for infile, outfile, op in zip(infiles, outfiles, ops):
#    infile=infiles[1]
#    outfile = outfiles[1]+'_2'
#    op = gdal.WarpOptions(format='ENVI', outputBounds=outputBounds, errorThreshold=0,
#                            xRes= 480, yRes =480, dstSRS = 'EPSG:32613', 
#                            resampleAlg=gdal.GRA_NearestNeighbour, dstNodata=-1.0)
    print("reprojecting file %s to %s" % (infile, outfile)) 
    tmp = gdal.Warp(outfile, infile, options=op)
    del tmp
    if os.path.isfile(outfile+'.hdr'):
        outF=open(outfile+'.hdr', 'a')
        outF.write("data ignore value = %s" % -1.0)
        outF.close()
    try :
        os.remove(outfile+".aux.xml")
    except OSError:
        pass
#    gdal.Warp(Filename, vsicurl_url, cutlineDSName = Crop_file, cropToCutline = True)  ###Change here the projection using dstSRS
    
####################           RASTER REPROJECTION                ##############################   


