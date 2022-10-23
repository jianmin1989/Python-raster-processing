    #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 23:31:29 2020

@author: wangj
"""


Pathin = '/gpfs/data/xyz/jianmin/LandsatARD/reflectance/' 
Pathout = '/gpfs/data/xyz/jianmin/LandsatARD/3dcomposite/'
Tiles = ['011011'] ## ['T13SDB', 'T13SEB', 'T13SDA', 'T13SEA'] "13SDV"
Years = list(range(2000, 2019))   ### Only 2017 
Outstyle = 0  ### IF 0, the output are separated files of each is a composite (day);;;  if 1, the output is the one BIP file for a tile year
Nyear = 2 #### 1 for the current year and 2 for the previous and proceeding half year and the current year. 
Nday = 3  ## How many days composite
Nrow = 5000
Ncol = 5000
FILL = 32767
Sensors=['LT04', 'LT05', 'LE07', 'LC08']






def main(pathin, pathout, tiles, years, outstyle, nyear, nday):
    import sys
    import numpy as np 
    import os
    import time
    import re

    if not os.path.exists(pathout):
        os.mkdir(pathout)
               
    time1 = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    years = sorted(years) ### or years.sort()
    allfiles = os.listdir(pathin)
    allfiles = [file for file in allfiles if re.match(".*_SR\.tar", file)]
    for tile in tiles:
        allfiles_t = [file for file in allfiles if re.match(".*CU_"+tile+"_.*_SR\.tar", file)]
        for year in years:
            print('generating tile %s in year %s' % (tile, year))
            
            dds, yys = dayandyear(year, years, nyear, nday)
                
            for dd, yy, in zip(dds, yys):
                ### To generate the file names in each year 
                dates = [str(yy) + str(item).zfill(3) for item in dd+np.arange(nday)] 
                mmdds = [date_transfer(yy, item, -1) for item in dd+np.arange(nday)]
#                    outfile = pathout + 'VIIRS_VI.'+dates[0]+'.' + tile + '.BIP.gz'
                outfile = '%sLARD_VI.%s.%s.BIP.gz' % (pathout, dates[0], tile)
                
                if (os.path.isfile(outfile)) and (os.path.getsize(outfile) > 50): 
                    print("      %s exists and file size is > 50 kb" % outfile)
                    continue
                else:
                    print("generating %s" % outfile)
                    
                #files_ref = [pathref + 'VNP43IA4.A' + item + '.'+ tile + '.001.h5' for item in dates]
                #[pathref + 'VNP43IA4.A' + item1 + '.h27v07.001.h5' + item2 for item1, item2 in zip(dates, a)]
                #files_qa = [pathqa + 'VNP43IA2.A' + item + '.'+ tile + '.001.h5' for item in dates]
                #files_ref = ["%sHLS.%s.%s.%s.v1.4.tar" % (pathref, sensor, tile, item) for sensor in Sensors for item in dates ]
                files=[file for file in allfiles_t for mmdd in mmdds if re.match(".*CU_"+tile+"_"+mmdd+"_.*", file)]
                
                
                
                start = time.process_time()
                Success = manage_composite(pathin, files, outfile)
                end = time.process_time()
                print(end-start)
                
                if Success != 1 :
                    print('Compositing has problems with day %s in year %s' % (dd, yy))
                    sys.exit(1)
                        
        if outstyle == 1:
            #reognize the output files
            for year in years:
                print(" for tile %s year %s" % (tile, year))
                #outfile = pathout + 'VIIRS_VI.'+str(year)+'.' + tile + '.year' + str(nyear) + '.BIP.gz'
                outfile = '%sAHI_VI.%s.%s.year%s.BIP' % (pathout, year, tile, nyear)
                if os.path.isfile(outfile) : 
                    continue
                dds, yys = dayandyear2(year, nyear, nday)                  
                
                # files=[pathout + 'VIIRS_VI.'+ str(yy) + str(dd).zfill(3) +'.' + tile + '.BIP.gz' for yy, dd in zip(yys, dds)]
                files=['%sAHI_VI.%s%s.%s.BIP.gz' % (pathout, yy, str(dd).zfill(3), tile) for yy, dd in zip(yys, dds)]
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
    
    DDs=np.arange(1, 366, nday)
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


def manage_composite(pathin, files, outfile):
    import gzip 
#    import struct 
    import numpy as np
    import time
    
    with gzip.open(outfile, 'wb') as f:
    #RES = np.zeros(Nrow, Ncol, 4)     
        if len(files)==0:
            #bytes_write = struct.pack('hhhh', int(FILL), int(9), int(FILL), int(FILL))
            #[f.write(bytes_write) for i in range(Nrow) for j in range(Ncol)]
            res=np.full((Nrow, Ncol, 4), FILL, dtype=np.int16)
            res[:,:,1]=9
        
        else:
            start1 = time.process_time()   
            EVI2, NDVI, NDWI, QA = index_calculate(pathin, files)
            end1 = time.process_time()
            print("                           index calculation finished, spent %f s" % (end1-start1))
            
            start2 = time.process_time()
            evi2, qa, ndvi, ndwi = COMPOSITE(EVI2, NDVI, NDWI, QA)
            end2 = time.process_time()
            print("                           composition finished, spent %f s" % (end2-start2))
            res=np.stack((evi2, qa, ndvi, ndwi), axis=2)
            #[f.write(struct.pack('hhhh', tmpevi2, tmpqa, tmpndvi, tmpndwi)) for tmpevi2, tmpqa, tmpndvi, tmpndwi in zip(np.nditer(evi2), np.nditer(qa), np.nditer(ndvi), np.nditer(ndwi))]

        f.write(res.tostring())


    return(1)
    
        
        
def index_calculate(pathin, files) : 
    import numpy as np
    import os
    import shutil
    import re
    import tarfile
    from osgeo import gdal
#    from PIL import Image
    
    
    tmppath = os.path.join(pathin, 'tmp')
    if not os.path.exists(tmppath):
        os.mkdir(tmppath)
    #print('%s has been created' % tmppath)
     
    n=len(files)
    EVI2 = np.full((Nrow, Ncol, n), FILL, dtype=np.int16)
    NDVI = np.full((Nrow, Ncol, n), FILL, dtype=np.int16)
    NDWI = np.full((Nrow, Ncol, n), FILL, dtype=np.int16)
    QA = np.full((Nrow, Ncol, n), 4, dtype=np.int16)
    
    for k in range(n):
        print('openning %s' % os.path.join(pathin, files[k]))
        with tarfile.open(os.path.join(pathin, files[k])) as tar:
            def is_within_directory(directory, target):
                
                abs_directory = os.path.abspath(directory)
                abs_target = os.path.abspath(target)
            
                prefix = os.path.commonprefix([abs_directory, abs_target])
                
                return prefix == abs_directory
            
            def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
            
                for member in tar.getmembers():
                    member_path = os.path.join(path, member.name)
                    if not is_within_directory(path, member_path):
                        raise Exception("Attempted Path Traversal in Tar File")
            
                tar.extractall(path, members, numeric_owner=numeric_owner) 
                
            
            safe_extract(tar, tmppath)
        
        if re.search("LC08_", files[k]):
            bred='SRB4.tif'
            bnir='SRB5.tif'
            bswir='SRB4.tif'
            bqa='PIXELQA.tif'
        else:
            bred='SRB3.tif'
            bnir='SRB4.tif'
            bswir='SRB5.tif'
            bqa='PIXELQA.tif'
        
        
        file_red=os.path.join(tmppath, files[k][0:-6]+bred)
        file_nir=os.path.join(tmppath, files[k][0:-6]+bnir)
        file_swir=os.path.join(tmppath, files[k][0:-6]+bswir)
        file_qa=os.path.join(tmppath, files[k][0:-6]+bqa)
    
        ds=gdal.Open(file_red)
        ref_I1 = np.array(ds.GetRasterBand(1).ReadAsArray())
        ds = None
        ds=gdal.Open(file_nir)
        ref_I2 = np.array(ds.GetRasterBand(1).ReadAsArray())
        ds = None
        ds=gdal.Open(file_swir)
        ref_I3 = np.array(ds.GetRasterBand(1).ReadAsArray())
        ds = None
        ds=gdal.Open(file_qa)
        qar = np.array(ds.GetRasterBand(1).ReadAsArray())
        ds = None
        
#        ds=Image.open(file_red)
#        ref_I1=np.array(ds)
#        ds.close()
#        ds=Image.open(file_nir)
#        ref_I2=np.array(ds)
#        ds.close()
#        ds=Image.open(file_swir)
#        ref_I3=np.array(ds)
#        ds.close()
#        ds=Image.open(file_qa)
#        qar=np.array(ds)
#        ds.close()
            
        
        ### Calculate index
        evi2 = 25000.0 * (ref_I2 - ref_I1)/ (ref_I2 + 2.4*ref_I1 + 10000)
        ndvi = 10000.0*(ref_I2-ref_I1)/(ref_I2+ref_I1)
        ndwi = 10000.0*(ref_I2-ref_I3)/(ref_I2+ref_I3)
        #mask = (ref_I1 <=0) | (ref_I1 >= 10000) | (ref_I2 <=0) | (ref_I2 >= 10000) | (ref_I3 <=0) | (ref_I3 >= 10000)
        #ref_I1_valid = np.ma.MaskedArray(ref_I1, mask)
        
        tmpind = (ref_I1 <=0) | (ref_I1 >= 10000) | (ref_I2 <=0) | (ref_I2 >= 10000)  
        evi2[tmpind] = FILL
        ndvi[tmpind] = FILL
        tmpind = (ref_I3 <=0) | (ref_I3 >= 10000) | (ref_I2 <=0) | (ref_I2 >= 10000)
        ndwi[tmpind] = FILL
            
            
            ### Record QA
#7-6    Cloud Confidence:
#00 - None
#01 - Low
#10 - Medium
#11 - High
#5      Cloud
#4      snow/ice
#3      cloud shadow
#2      Water
#1      Clear
#0      FILL     
            
            
        qa = np.full((Nrow, Ncol), 0, dtype=np.int16)
        qa[(qar // 32)%2 ==1] = 4
        qa[(qar // 4)%2 ==1] = 100
        qa[(qar // 16)%2 ==1] = 1
        qa[(evi2 == FILL) | (ndvi==FILL) | (qar%2 == 1)] = 9 ##FILL
        
        EVI2[:,:,k] = evi2.astype(np.int16)
        NDVI[:,:,k] = ndvi.astype(np.int16)
        NDWI[:,:,k] = ndwi.astype(np.int16)
        QA[:,:,k] = qa
            
    #            VI_vec = np.vectorize(VI, otypes=[np.int16, np.int16, np.int16, np.int16])
    #            evi2, ndvi, ndwi, qa = VI_vec(ref_I1, ref_I2, ref_I3, qa_I1, qa_I2, water, snow)
    
    shutil.rmtree(tmppath)
    return(EVI2, NDVI, NDWI, QA)


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

def COMPOSITE(EVI2, NDVI, NDWI, QA):
    import numpy as np
    
    n=np.shape(EVI2)[2]
    outEVI2 = np.full([Nrow, Ncol], FILL, dtype=np.int16)
    outNDVI = np.full([Nrow, Ncol], FILL, dtype=np.int16)
    outNDWI = np.full([Nrow, Ncol], FILL, dtype=np.int16)
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
    maskedNDWI = np.ma.array(NDWI, mask = indexqa)
    tmpEVI2 = np.amax(maskedEVI2, axis=2)
    tmpNDVI = np.amax(maskedNDVI, axis=2)
    tmpNDWI = np.amax(maskedNDWI, axis=2)
    
#    EVI2[indexqa]=np.nan
#    NDVI[indexqa]=np.nan
#    NDWI[indexqa]=np.nan    
#    tmpEVI2 = np.nanmax(EVI2, axis=2)
#    tmpNDVI = np.nanmax(NDVI, axis=2)
#    tmpNDWI = np.nanmax(NDWI, axis=2)
    
    outEVI2[index] = tmpEVI2[index]
    outNDVI[index] = tmpNDVI[index]
    outNDWI[index] = tmpNDWI[index]
    
    
    outQA[(outQA!=100) & (np.any(QA==1, axis = 2))] = 1
    
    return(outEVI2, outQA, outNDVI, outNDWI)
    
    
        
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
                
                
main(Pathin, Pathout, Tiles, Years, Outstyle, Nyear, Nday)        


