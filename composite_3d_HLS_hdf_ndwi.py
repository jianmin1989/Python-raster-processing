#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 15:29:04 2020

@author: wangj
"""

Pathref = '/gpfs/data/xyz/jianmin/HLSdata/reflectance/' 
Pathqa = '/gpfs/data/xyz/jianmin/HLSdata/reflectance/' 
Pathout = '/gpfs/data/xyz/jianmin/HLSdata/3dcomposite/'
Tiles =  ['T13SDB', 'T13SEB'] ## ['T13SDB', 'T13SEB', 'T13SDA', 'T13SEA'] "13SDV"
Years = list(range(2018, 2019))   ### Only 2017 
Outstyle = 0  ### IF 0, the output are separated files of each is a composite (day);;;  if 1, the output is the one BIP file for a tile year
Nyear = 2 #### 1 for the current year and 2 for the previous and proceeding half year and the current year. 
Nday = 3  ## How many days composite
Nrow = 3660
Ncol = 3660
FILL = 32767
Sensors=['L30', 'S30']






def main(pathref, pathqa, pathout, tiles, years, outstyle, nyear, nday):
    import sys
    import numpy as np 
    import os
    import time

    if not os.path.exists(pathout):
        os.mkdir(pathout)
    time1 = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    years = sorted(years) ### or years.sort()
    for tile in tiles:
        for year in years:
            dds, yys = dayandyear(year, years, nyear, nday)
                
            for dd, yy, in zip(dds, yys):
                ### To generate the file names in each year 
                dates = [str(yy) + str(item).zfill(3) for item in dd+np.arange(nday)] 
#                    outfile = pathout + 'VIIRS_VI.'+dates[0]+'.' + tile + '.BIP.gz'
                outfile = '%sHLS_VI.%s.%s.BIP.gz' % (pathout, dates[0], tile)
                
                if (os.path.isfile(outfile)) and (os.path.getsize(outfile) > 50): 
                    print("      %s exists and file size is > 50 kb" % outfile)
                    continue
                else:
                    print("generating %s" % outfile)
                    
                #files_ref = [pathref + 'VNP43IA4.A' + item + '.'+ tile + '.001.h5' for item in dates]
                #[pathref + 'VNP43IA4.A' + item1 + '.h27v07.001.h5' + item2 for item1, item2 in zip(dates, a)]
                #files_qa = [pathqa + 'VNP43IA2.A' + item + '.'+ tile + '.001.h5' for item in dates]
                files_ref = ["%sHLS.%s.%s.%s.v1.4.hdf" % (pathref, sensor, tile, item) for sensor in Sensors for item in dates ]
                files_qa = files_ref
                
                start = time.process_time()
                Success = manage_composite(files_ref, files_qa, outfile)
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


def manage_composite(files_ref, files_qa, outfile):
    import gzip 
#    import struct 
    import numpy as np
    import time
    
    with gzip.open(outfile, 'wb') as f:
    #RES = np.zeros(Nrow, Ncol, 4)
        start1 = time.process_time()
        EVI2, NDVI, NDWI, QA = index_calculate(files_ref, files_qa)
        end1 = time.process_time()
        print("                           index calculation finished, spent %f s" % (end1-start1))
        if EVI2 is None:
            #bytes_write = struct.pack('hhhh', int(FILL), int(9), int(FILL), int(FILL))
            #[f.write(bytes_write) for i in range(Nrow) for j in range(Ncol)]
            res=np.full((Nrow, Ncol, 4), FILL, dtype=np.int16)
            res[:,:,1]=9
        
        else:
            start2 = time.process_time()
            evi2, qa, ndvi, ndwi = COMPOSITE(EVI2, NDVI, NDWI, QA)
            end2 = time.process_time()
            print("                           composition finished, spent %f s" % (end2-start2))
            res=np.stack((evi2, qa, ndvi, ndwi), axis=2)
            #[f.write(struct.pack('hhhh', tmpevi2, tmpqa, tmpndvi, tmpndwi)) for tmpevi2, tmpqa, tmpndvi, tmpndwi in zip(np.nditer(evi2), np.nditer(qa), np.nditer(ndvi), np.nditer(ndwi))]

        f.write(res.tostring())


    return(1)
    
        
        
def index_calculate(REF_files, QA_files) : 
    import numpy as np
    import os
    from pyhdf.SD import SD, SDC
    import re
    
    [print("%s is missed" % file) for file in REF_files if not os.path.isfile(file)]
#    QA_files=[file for file, file_ref in zip(QA_files, REF_files) if os.path.isfile(file_ref)]
    REF_files=[file for file in REF_files if os.path.isfile(file)]
    
    n=len(REF_files)
    if n > 0:  
        EVI2 = np.full((Nrow, Ncol, n), FILL, dtype=np.int16)
        NDVI = np.full((Nrow, Ncol, n), FILL, dtype=np.int16)
        NDWI = np.full((Nrow, Ncol, n), FILL, dtype=np.int16)
        QA = np.full((Nrow, Ncol, n), 4, dtype=np.int16)
        k=0
        for k in range(n):
            file_ref = REF_files[k]
            if re.search("HLS.L30.", file_ref):
                bred='band04'
                bnir='band05'
                bswir='band06'
                bqa='QA'
            else:
                bred='B04'
                bnir='B8A'
                bswir='B11'
                bqa='QA'
            ### Read Files
            hand_ref = SD(file_ref, SDC.READ)
            ref_I1 = hand_ref.select(bred).get()
            ref_I2 = hand_ref.select(bnir).get()
            ref_I3 = hand_ref.select(bswir).get()
            qar = hand_ref.select(bqa).get() 
            hand_ref.end()
            
            
            
            #(ref_I2[0,0] - ref_I3[0,0])/ (ref_I2[0,0] + 2.4*ref_I3[0,0] + 10000.0)
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
#7-6    aerosol:
#00 - climatology
#01 - low
#10 - average
#11 - high
#5      water
#4      snow/ice
#3      cloud shadow
#2      adjacent to cloud
#1      cloud
#0      cirrus     
            
            
            qa = np.full((Nrow, Ncol), 0, dtype=np.int16)
            qa[(qar % 64 > 0) | (qar//64 == 3)] = 4
            qa[(qar // 32)%2 ==1] = 100
            qa[(qar // 16)%2 ==1] = 1
            qa[(evi2 == FILL) | (ndvi==FILL)] = 9 ##FILL
            
            EVI2[:,:,k] = evi2.astype(np.int16)
            NDVI[:,:,k] = ndvi.astype(np.int16)
            NDWI[:,:,k] = ndwi.astype(np.int16)
            QA[:,:,k] = qa
            
    #            VI_vec = np.vectorize(VI, otypes=[np.int16, np.int16, np.int16, np.int16])
    #            evi2, ndvi, ndwi, qa = VI_vec(ref_I1, ref_I2, ref_I3, qa_I1, qa_I2, water, snow)
    
        return(EVI2, NDVI, NDWI, QA)
    else:
        return(None, None, None, None)


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
                
                
main(Pathref, Pathqa, Pathout, Tiles, Years, Outstyle, Nyear, Nday)        


