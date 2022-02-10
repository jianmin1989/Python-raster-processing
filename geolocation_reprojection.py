#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 12:53:33 2022

@author: jianmin.wang
"""

###Covert a point from a coordinate system to another (reproject)
from osgeo import ogr, osr
import sys
sys.path.append("/home/jianmin.wang/codes/self_packages/python3")

InSR = osr.SpatialReference()
InSR.ImportFromEPSG(4326)       # WGS84/Geographic
OutSR = osr.SpatialReference()
# OutSR.ImportFromWkt(JWraster.get_sinu_projection())
file = "/scratch/jianmin.wang/data/MCD12Q1_bip/test_h08v04"
OutSR.ImportFromWkt(gdal.Open(file).GetProjection())

Point = ogr.Geometry(ogr.wkbPoint)
Point.AddPoint( 38.43091667,-120.9658861) # use your coordinates here (lat, longitude)
Point.AssignSpatialReference(InSR)    # tell the point what coordinates it's in
Point.TransformTo(OutSR)              # project it to the out spatial reference
print('{0},{1}'.format(Point.GetX(),Point.GetY())) # output projected X and Y coordinates
# 39.232,-121.2972
