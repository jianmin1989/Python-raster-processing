#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 23:31:29 2020

@author: wangj
"""

### This code is study area targeted, and thus no tiles. The tiles is actually the study area code. Different areas even with different projection
Pathin = '/hunter/data1/wangj/planet_data/HLS18TXM/Reflectance/'
Tile = 'HLS18TXM'
Pathsnow = '/gpfs/data/xyz/jianmin/MODIS1km/MOD09GA/HLS18TXM_planet/'
#snowTiles = ['h12v04']
Pathout = '/hunter/data1/wangj/planet_data/HLS18TXM/3dcomposite/'
## ['T13SDB', 'T13SEB', 'T13SDA', 'T13SEA'] "13SDV"
Years = list(range(2017, 2018))   ### Only 2017 
Outstyle = 0  ### IF 0, the output are separated files of each is a composite (day);;;  if 1, the output is the one BIP file for a tile year
Nyear = 2 #### 1 for the current year and 2 for the previous and proceeding half year and the current year. 
Nday = 3  ## How many days composite
FILL = 32767


Reso = 3
#DstSRS = 'EPSG:32610'  #UTM 10N
#OutputBounds = [677355, 4253628, 679155, 4255428] #[minx, miny, maxx, maxy]
DstSRS = 'EPSG:32618'  #UTM 18N
OutputBounds = None #[minx, miny, maxx, maxy]  #needs to be multiply of reso
#Can be created by http://geojson.io/#map=10/33.7392/-107.9874
Regionfile = '/hunter/data/jianmin/codes/planet/Planet_api_J1/HLS18TXM/interested_area/POLYGON.shp'
Cut2region = True ## False   True #This is just an initial value. the real value would be determined in the code
#Regionfile and OutputBounds can be used together. 
#The spatial polygons will be used to check whether the imagery is overlapped with the regionfile
#OutputBounds is used to set the spatial range of the output
# Regionfile could be irregular shape while OutputBounds must be a rectangle. 
# At least one of OutputBounds and Regionfile is needed. 
# One can be calculated based on another.  


#### This line below is to get a hdr file just for test
#ds = gdal.Warp('/hunter/data1/wangj/planet_data/HLS18TXM/3dcomposite/test', 
##                       '/hunter/data1/wangj/planet_data/Cal_Liu/Reflectance/20190629_182839_0f4e_3B_AnalyticMS_SR.tif',  
#               '/hunter/data1/wangj/planet_data/HLS18TXM/Reflectance/20160703_145144_0e16_3B_AnalyticMS_SR.tif',
#               format = 'ENVI', errorThreshold=0,
#               cutlineLayer = Region, cropToCutline = Cut2region, #Can specify a regionfile here
#               outputBounds = OutputBounds, #check transformerOptions for tie point loc, not sure.
#               xRes= Reso, yRes = Reso, dstSRS = DstSRS, targetAlignedPixels =True, 
#               resampleAlg=gdal.GRA_NearestNeighbour, srcNodata=0, dstNodata=32767)
#ds = None




def main():
    import sys
    import numpy as np 
    import os
    import time
    import re
    import geopandas as gpd
    from shapely.geometry import Polygon

    if not os.path.exists(Pathout):
        os.mkdir(Pathout)
    
    
    #### This block is to define the output grids and spatial reference
    ## There three 
    global Nrow, Ncol, Region, Cut2region, Years, OutputBounds
    if (OutputBounds is None) and (Regionfile is None) :
        sys.exit('Please at least setup one of OutputBounds or Regionfile')
    elif (OutputBounds is not None) and (Regionfile is None):
        OutputBounds = [round(x/Reso)*Reso for x in OutputBounds]  #Make sure it is multiply of Reso
        longs = [OutputBounds[0], OutputBounds[0], OutputBounds[2], OutputBounds[2], OutputBounds[0]]
        lats = [OutputBounds[3], OutputBounds[1], OutputBounds[1], OutputBounds[3], OutputBounds[3]]
        Region = gpd.GeoDataFrame(index=[0], crs={'init':DstSRS}, geometry=[Polygon(zip(longs, lats))])
        Cut2region = True #in fact it doesn't matter to be T or F
    elif (OutputBounds is None) and (Regionfile is not None):
        Region = gpd.read_file(Regionfile)
#        print(Region.crs)
        Region = Region.to_crs({'init':DstSRS} )
        OutputBounds = Region.bounds.values.tolist()[0]
        OutputBounds = [round(x/Reso)*Reso for x in OutputBounds]
    else:
        Cut2region = True
    
    Nrow=round((OutputBounds[3]-OutputBounds[1])/Reso)
    Ncol=round((OutputBounds[2]-OutputBounds[0])/Reso)
    print(OutputBounds)
    print("Nrow = %s  Ncol = %s" % (Nrow, Ncol))    
    #### This block is to define the output grids and spatial reference    
        
    
    #### Collect the necessary input files.       
    time1 = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    Years = sorted(Years) ### or years.sort()
    tile = Tile
    allfiles = os.listdir(Pathin)
    allfiles = [file for file in allfiles if re.match(".*_SR\.tif", file)]
    if Pathsnow is not None:
        allfiles_snow = os.listdir(Pathsnow)
#        pattern = '(' + '|'.join(snowTiles) + ')'
#        allfiles = [file for file in allfiles if re.match("^MOD09GA\.A\d{7}\.%s\.006.hdf" % pattern, file)] #MOD09GA.A2014190.h09v05.006.hdf
        allfiles_snow = [file for file in allfiles_snow if re.match("^MOD09GA\.A\d{7}\.%s.*tif" % tile, file)] #MOD09GA.A2014190.h09v05.006.hdf
        allfiles_snow = sorted(allfiles_snow)
    for year in Years:
        print('generating tile %s in year %s' % (tile, year))
        
        dds, yys = dayandyear(year, Years, Nyear, Nday)
            
        for dd, yy, in zip(dds, yys):
            ### To generate the file names in each year 
            dates = [str(yy) + str(item).zfill(3) for item in dd+np.arange(Nday)] 
            mmdds = [date_transfer(yy, item, -1) for item in dd+np.arange(Nday)]
#                    outfile = pathout + 'VIIRS_VI.'+dates[0]+'.' + tile + '.BIP.gz'
            outfile = '%sPlanet_VI.%s.%s.BIP.gz' % (Pathout, dates[0], tile)
            
            if (os.path.isfile(outfile)) and (os.path.getsize(outfile) > 50): 
                print("      %s exists and file size is > 50 kb" % outfile)
                continue
            else:
                print("generating %s" % outfile)
                
            #files_ref = [pathref + 'VNP43IA4.A' + item + '.'+ tile + '.001.h5' for item in dates]
            #[pathref + 'VNP43IA4.A' + item1 + '.h27v07.001.h5' + item2 for item1, item2 in zip(dates, a)]
            #files_qa = [pathqa + 'VNP43IA2.A' + item + '.'+ tile + '.001.h5' for item in dates]
            #files_ref = ["%sHLS.%s.%s.%s.v1.4.tar" % (pathref, sensor, tile, item) for sensor in Sensors for item in dates ]
            files = [file for file in allfiles for mmdd in mmdds if re.match(mmdd+"_.*", file)]
            if Pathsnow is not None:
                pattern = '(' + '|'.join(dates) + ')'
                files_snow = [Pathsnow + file for file in allfiles_snow if re.match('MOD09GA.A%s' % pattern, file)]
            else :
                files_snow = None
                
            start = time.process_time()
            Success = manage_composite(Pathin, files, files_snow, outfile)
            end = time.process_time()
            print(end-start)
            
            if Success != 1 :
                print('Compositing has problems with day %s in year %s' % (dd, yy))
                sys.exit(1)
                    
    if Outstyle == 1:
        #reognize the output files
        for year in Years:
            print(" for year %s" % year)
            #outfile = Pathout + 'VIIRS_VI.'+str(year)+'.' + tile + '.year' + str(nyear) + '.BIP.gz'
            outfile = '%sAHI_VI.%s.%s.year%s.BIP' % (Pathout, year, tile, Nyear)
            if os.path.isfile(outfile) : 
                continue
            dds, yys = dayandyear2(year, Nyear, Nday)                  
            
            # files=[pathout + 'VIIRS_VI.'+ str(yy) + str(dd).zfill(3) +'.' + tile + '.BIP.gz' for yy, dd in zip(yys, dds)]
            files=['%sAHI_VI.%s%s.%s.BIP.gz' % (Pathout, yy, str(dd).zfill(3), tile) for yy, dd in zip(yys, dds)]
            Success = stackgzfiles(files, outfile)
                    
    time2 = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())      
    print('started at %s' % time1)
    print('ended at %s' % time2)
    


# dd <0 returns yyyymmdd , or mode=1 returns yyyydoy
def date_transfer(year, mm, dd):
    mondint=[0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
    if year%4 == 0 and (year%100 !=0 or year%400 ==0)  :
        mondint=[0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
        
    if dd>0 :
        doy=mondint[mm-1]+dd
        return(str(year) + str(doy).zfill(3) )
    elif dd<=0:
        doy=mm
        mm=0
        while mm<=11 and mondint[mm] <= doy:
            mm=mm+1
        dd=doy-mondint[mm-1]
        return(str(year) + str(mm).zfill(2) + str(dd).zfill(2) )
        



#def date_transfer_manage(dates):
#    import sys
#    
#    mode = len(dates[0])
#    if mode==7:
#        years=[int(date[0:4]) for date in dates]
#        doys=[int(date[4:7]) for date in dates]
#        
#        
#    elif mode==8:
#        years=[date[0:4] for date in dates]
#    else:
#        sys.exit('Please input dates as "2000001" or "20000101"!')
    
    
    

    
def dayandyear(year, years, nyear, nday):
    import numpy as np
    
    DDs=np.arange(1, 367, nday)
    YYs=np.array([year]*len(DDs))
    dds=DDs
    yys=YYs
    if nyear == 2:
        Nfile = len(DDs)
        Nfilehalf = Nfile//2
        if year-1 not in years:
            dds = np.concatenate((DDs[Nfilehalf:Nfile], DDs), axis=0)
            yys = np.concatenate((np.array([year-1]*(Nfile-Nfilehalf)), YYs), axis=0)
            
        if year+1 not in years:
            dds = np.concatenate((dds, DDs[0:Nfilehalf]), axis=0)
            yys = np.concatenate((yys, np.array([year+1]*Nfilehalf)), axis=0)   
        
    return(dds, yys)


def dayandyear2(year, nyear, nday):
    import numpy as np
    
    DDs=np.arange(1, 366, nday)
    YYs=np.array([year]*len(DDs))
    dds=DDs
    yys=YYs
    if nyear == 2:
        Nfile = len(DDs)
        Nfilehalf = Nfile//2
        dds = np.concatenate((DDs[Nfilehalf:Nfile], DDs), axis=0)
        yys = np.concatenate((np.array([year-1]*(Nfile-Nfilehalf)), YYs), axis=0)
            
        dds = np.concatenate((dds, DDs[0:Nfilehalf]), axis=0)
        yys = np.concatenate((yys, np.array([year+1]*Nfilehalf)), axis=0)   
        
    return(dds, yys)   


def manage_composite(pathin, files, files_snow, outfile):
    import gzip 
#    import struct 
    import numpy as np
    import time
    
    with gzip.open(outfile, 'wb') as f:
    #RES = np.zeros(Nrow, Ncol, 4)     
        if len(files)==0:
            #bytes_write = struct.pack('hhhh', int(FILL), int(9), int(FILL), int(FILL))
            #[f.write(bytes_write) for i in range(Nrow) for j in range(Ncol)]
            res=np.full((Nrow, Ncol, 3), FILL, dtype=np.int16)
            res[:,:,1]=9
        
        else:
            start1 = time.process_time()   
            EVI2, NDVI, QA = index_calculate(pathin, files, files_snow)
            end1 = time.process_time()
            print("                           index calculation finished, spent %f s" % (end1-start1))
            if EVI2 is None: 
                res=np.full((Nrow, Ncol, 3), FILL, dtype=np.int16)
                res[:,:,1]=9
            else: 
                start2 = time.process_time()
                evi2, qa, ndvi = COMPOSITE(EVI2, NDVI, QA)
                end2 = time.process_time()
                print("                           composition finished, spent %f s" % (end2-start2))
                res=np.stack((evi2, qa, ndvi), axis=2)
                #[f.write(struct.pack('hhhh', tmpevi2, tmpqa, tmpndvi)) for tmpevi2, tmpqa, tmpndvi in zip(np.nditer(evi2), np.nditer(qa), np.nditer(ndvi))]

        f.write(res.tostring())


    return(1)
    

def check_for_overlap(file, bound):
    from osgeo import gdal
    ds = gdal.Open(file)
    rows = ds.RasterYSize
    cols = ds.RasterXSize
    gt = ds.GetGeoTransform()
    minx = gt[0]
    maxy = gt[3]
    maxx = gt[0] + gt[1] * cols
    miny = gt[3] + gt[5] * rows
    rect = [minx, miny, maxx, maxy]
    ds = None
    
    overlap = 1
    if(rect[3]<bound[1] or rect[2]<bound[0]):
        overlap = 0
    #the blue color                          or   green 
    elif(rect[1]>bound[3] or rect[0]>bound[2]):
        overlap = 0
    
    return(overlap)


    
        
def index_calculate(pathin, files, files_snow) : 
    import numpy as np
    import os
    from osgeo import gdal
    import re
#    from PIL import Image
    
    REF_files = files #[os.path.join(pathin, file) for file in files]
    QA_files = [file.replace('SR.tif', 'DN_udm.tif') for file in REF_files]
    ids = [iid for iid, file_qa in enumerate(QA_files) if os.path.isfile(os.path.join(pathin, file_qa)) ]
#    ids = [iid for iid, (file_qa, file_ref) in enumerate(zip(QA_files, REF_files)) if os.path.isfile(file_qa) and os.path.isfile(file_ref)]
#    QA_files=[QA_files[i] for i in ids]
    REF_files=[REF_files[i] for i in ids]
    
    REF_files=[file for file in REF_files if check_for_overlap(os.path.join(pathin, file), OutputBounds)==1]
    REF_files=sorted(REF_files)
    QA_files = [file.replace('SR.tif', 'DN_udm.tif') for file in REF_files]
    
    
    [print('%s selected' % file) for file in REF_files]
     
    n=len(REF_files)
    if n<=0:
        return(None, None, None)
        
    EVI2 = np.full((Nrow, Ncol, n), FILL, dtype=np.int16)
    NDVI = np.full((Nrow, Ncol, n), FILL, dtype=np.int16)
    QA = np.full((Nrow, Ncol, n), 4, dtype=np.int16)
    for k in range(n):
        file_ref = os.path.join(pathin, REF_files[k])
        file_qa = os.path.join(pathin, QA_files[k])
#        file_ref = '/hunter/data1/wangj/planet_data/Cal_Liu/20190629_182839_0f4e_3B_AnalyticMS_SR.tif'
#        file_qa = file_ref.replace('SR.tif', 'DN_udm.tif') 
        ds = gdal.Warp('', file_ref, format = 'VRT', errorThreshold=0,
               cutlineLayer = Region, cropToCutline = Cut2region, #Can specify a regionfile here
               outputBounds = OutputBounds, #check transformerOptions for tie point loc, not sure.
               xRes= Reso, yRes = Reso, dstSRS = DstSRS, targetAlignedPixels =True, 
               resampleAlg=gdal.GRA_NearestNeighbour, srcNodata=0, dstNodata=32767)
        
        ref_I1 = np.array(ds.GetRasterBand(3).ReadAsArray())
        ref_I2 = np.array(ds.GetRasterBand(4).ReadAsArray())
        ds = None
        
        ds = gdal.Warp('', file_qa, format = 'VRT', errorThreshold=0,
               cutlineLayer = Region, cropToCutline = Cut2region, #Can specify a regionfile here
               outputBounds = OutputBounds, #check transformerOptions for tie point loc, not sure.
               xRes= Reso, yRes = Reso, dstSRS = DstSRS, targetAlignedPixels =True, 
               resampleAlg=gdal.GRA_NearestNeighbour, srcNodata=255, dstNodata=255)
        qar = np.array(ds.GetRasterBand(1).ReadAsArray())
        ds = None
        
        
        if files_snow is not None:
            if k==0: 
                olddate = 'nodate'
            date = date_transfer(int(REF_files[k][0:4]), int(REF_files[k][4:6]), int(REF_files[k][6:8]))
            if date!=olddate:
                olddate = date
                file_snow = [file for file in files_snow if re.match('%sMOD09GA\.A%s' % (Pathsnow, date), file)]
                if len(file_snow) >0 :
#                    vrtds = [gdal.Translate('', re.search('SUBDATASET_2_NAME=(.*)\n', gdal.Info(file)).group(1), format='VRT') for file in file_snow]   
#                    ds = gdal.Warp('', vrtds, format = 'VRT', errorThreshold=0,
#                               cutlineLayer = Region, cropToCutline = Cut2region, #Can specify a regionfile here
#                               outputBounds = OutputBounds, #check transformerOptions for tie point loc, not sure.
#                               xRes= Reso, yRes = Reso, dstSRS = DstSRS, targetAlignedPixels =True, 
#                               resampleAlg=gdal.GRA_NearestNeighbour, srcNodata=65535, dstNodata=65535)
                    ds = gdal.Open(file_snow[0])
#                    ds = gdal.Warp('', file_snow[0], format = 'VRT', errorThreshold=0,
#                               cutlineLayer = Region, cropToCutline = Cut2region, #Can specify a regionfile here
#                               outputBounds = OutputBounds, #check transformerOptions for tie point loc, not sure.
#                               xRes= Reso, yRes = Reso, dstSRS = DstSRS, targetAlignedPixels =True, 
#                               resampleAlg=gdal.GRA_NearestNeighbour, srcNodata=65535, dstNodata=65535)
                    snow = np.array(ds.GetRasterBand(1).ReadAsArray())
                    ds = None
                else:
                    snow = None
            
            
         
#            state_1km_1 (720, 8)
#    16-bit unsigned integer,    1200 x 1200
#    Number of attributes = 6
#        long_name = 1km Reflectance Data State QA - first layer
#        units = bit field
#        valid_range = 0,57335
#        _FillValue = 65535
#        Nadir Data Resolution = 1km
#        QA index = 
#	Bits are listed from the MSB (bit 15) to the LSB (bit 0):
#	Bit    Description
#	15     internal snow algorithm flag; 
#	       1 -- yes, 0.00%
#	       0 -- no, 100.00%
#	14     Salt Pan;
#	       1 -- yes, 0.00%
#	       0 -- no, 100.00%
#	13     Pixel is adjacent to cloud;
#	       1 -- yes, 8.19%
#	       0 -- no, 91.81%
#	12     MOD35 snow/ice flag;
#	       1 -- yes, 0.04%
#	       0 -- no, 99.96%
#	11     internal fire algorithm flag;
#	       1 -- fire, 0.00%
#	       0 -- no fire, 100.00%
#	10     internal cloud algorithm flag;
#	       1 -- cloud, 63.34%
#	       0 -- no cloud, 36.66%
#	8-9    cirrus detected;
#	       00 -- none, 81.05%
#	       01 -- small, 2.92%
#	       10 -- average, 3.98%
#	       11 -- high, 12.06%
#	6-7    aerosol quantity;
#	       00 -- climatology, 65.58%
#	       01 -- low, 10.68%
#	       10 -- average, 14.51%
#	       11 -- high, 9.23%
#	3-5    land/water flag;
#	       000 -- shallow ocean, 9.15%
#	       001 -- land, 42.88%
#	       010 -- ocean coastlines and lake shorelines, 28.46%
#	       011 -- shallow inland water, 12.64%
#	       100 -- ephemeral water, 0.00%
#	       101 -- deep inland water, 3.28%
#	       110 -- continental/moderate ocean, 2.34%
#	       111 -- deep ocean, 1.24%
#	2      cloud shadow;
#	       1 -- yes, 4.20%
#	       0 -- no, 95.80%
#	0-1    cloud state;
#	       00 -- clear, 41.10%
#	       01 -- cloudy, 53.27%
#	       10 -- mixed, 5.64%
#	       11 -- not set, assumed clear, 0.00%
               
            
        evi2 = 25000.0 * (ref_I2 - ref_I1)/ (ref_I2 + 2.4*ref_I1 + 10000)
        ndvi = 10000.0*(ref_I2-ref_I1)/(ref_I2+ref_I1)
            
        
        #mask = (ref_I1 <=0) | (ref_I1 >= 10000) | (ref_I2 <=0) | (ref_I2 >= 10000) | (ref_I3 <=0) | (ref_I3 >= 10000)
        #ref_I1_valid = np.ma.MaskedArray(ref_I1, mask)
        
        tmpind = (ref_I1 <=0) | (ref_I1 >= 10000) | (ref_I2 <=0) | (ref_I2 >= 10000)  
        evi2[tmpind] = FILL
        ndvi[tmpind] = FILL
            
            
            
        qa = np.full((Nrow, Ncol), 0, dtype=np.int16)
        qa[(qar // 2)%2 ==1] = 4
        if snow is not None:
            qa[(snow // 32768)%2==1]= 1
        qa[(qar // 32)%2 ==1] = 9
        qa[(qar // 16)%2 ==1] = 9
        qa[(evi2 == FILL) | (ndvi==FILL) | (qar%2 == 1)] = 9 ##FILL
        
        EVI2[:,:,k] = evi2.astype(np.int16)
        NDVI[:,:,k] = ndvi.astype(np.int16)
        QA[:,:,k] = qa
            
    #            VI_vec = np.vectorize(VI, otypes=[np.int16, np.int16, np.int16, np.int16])
    #            evi2, ndvi, ndwi, qa = VI_vec(ref_I1, ref_I2, ref_I3, qa_I1, qa_I2, water, snow)
    
    return(EVI2, NDVI, QA)


#def VI(r1, r2, r3, q1, q2, wr, sn):
#    if (r1 <=0) or (r1 >= 10000) or (r2 <=0) or (r2 >= 10000):
#        eevi2=FILL
#        nndvi=FILL
#    else:
#        eevi2 = 25000.0 * (r2 - r1)/ (r2 + 2.4*r1 + 10000)
#        nndvi = 10000.0*(r2-r1)/(r2+r1)
#        
#    if (r3 <=0) or (r3 >= 10000) or (r2 <=0) or (r2 >= 10000):
#        nndwi=FILL
#    else:
#        nndwi = 10000.0*(r2-r3)/(r2+r3)
#    
#    qqa=0
#    if (q1 !=0) or (q2 != 0):
#        qqa=4
#    if (wr < 1) or (wr > 2):
#        qqa = 100
#    if (sn ==1):
#        qqa=1
#    if (eevi2 == FILL) or (nndvi == FILL) or (q1 >100) or (q2 >100):
#        qqa=9
#    
#    return(eevi2, nndvi, nndwi, qqa)



#a = np.zeros([3, 3, 4])
#def func(x):
#    x[x==0]=1
#func(a[0, 0, :])

def COMPOSITE(EVI2, NDVI, QA):
    import numpy as np
    
    n=np.shape(EVI2)[2]
    outEVI2 = np.full([Nrow, Ncol], FILL, dtype=np.int16)
    outNDVI = np.full([Nrow, Ncol], FILL, dtype=np.int16)
    outQA = np.full([Nrow, Ncol], 9, dtype=np.int16)
        
    index = (QA == 0) & ((EVI2 < -3000) | (EVI2 > 10000))
    QA[index] = 9
    EVI2[index]=FILL
    
    index = np.any(QA==100, axis = 2)
    outQA[index] = 100

    index = (np.logical_not(index)) & (np.any(QA<9, axis =2))
    outQA[index] = np.nanmin(QA[index], axis = 1)
    
    outQAn = np.repeat(outQA[:, :, np.newaxis], n, axis=2)
    indexqa = (QA != outQAn)
    indexqa[np.logical_not(index)] = True
    
    maskedEVI2 = np.ma.array(EVI2, mask = indexqa)
    maskedNDVI = np.ma.array(NDVI, mask = indexqa)
    tmpEVI2 = np.amax(maskedEVI2, axis=2)
    tmpNDVI = np.amax(maskedNDVI, axis=2)
    
#    EVI2[indexqa]=np.nan
#    NDVI[indexqa]=np.nan
#    NDWI[indexqa]=np.nan    
#    tmpEVI2 = np.nanmax(EVI2, axis=2)
#    tmpNDVI = np.nanmax(NDVI, axis=2)
#    tmpNDWI = np.nanmax(NDWI, axis=2)
    
    outEVI2[index] = tmpEVI2[index]
    outNDVI[index] = tmpNDVI[index]
    
    
    outQA[(outQA!=100) & (np.any(QA==1, axis = 2))] = 1
    
    return(outEVI2, outQA, outNDVI)
    
    
        
def stackgzfiles(files, outfile):
    import gzip 
    
#    with gzip.open(outfile, 'wb') as f:
    with open(outfile, 'wb') as f:
        try: 
            fins = [gzip.open(file, 'rb') for file in files]
        finally:
            k=0
            comp=True
            while comp:
                for fin in fins:
                    bytes_write = fin.read(8)
                    if not bytes_write:
                        comp=False
                        break
                    f.write(bytes_write)
                    
                k=k+1
                if(k % Ncol == 0):
                    irow=k/Ncol
                    if(irow % 100 ==0):
                        print("irow = %d" % irow)
                        
                
            for fin in fins:
                fin.close()
                
                
main()    

