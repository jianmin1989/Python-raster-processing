#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 12:58:20 2022

@author: jianmin.wang
"""

import sys, warnings, os, re
import subprocess
import math
import numpy as np
from osgeo import gdal, ogr, osr

#### write a header file for MODIS or VIIRS sinusoid projection
def writeMODISheader(files, ns=None, nl=None, nb=1, reso=None, FILL=32767, dtype=2, xstart=1, ystart=1, overwrite=True):
    """
    Date and Author: Sep. 21th, 2021 by Jianmin Wang
    Aim:
        write a header for MODIS/VIIRS binary tile.
    :param files: the files need header
    :param ns: number of samples
    :param nl: number of lines
    :param nb: number of bands
    :param reso: resolution, 250, 500, or 1000
    :param FILL: Fill value in header file
    :param dtype: data type in ENVI header file
    :param xstart: x start in ENVI header file
    :param ystart: y start in ENVI header file
    :param overwrite: if True, overwrite the header file
    """
    
    
    if type(files)==str: files = [files]
    
    
    ### PREPARE THE REFERENCE
    # reffile = "/home/jianmin.wang/codes/h10v04_1km_dt12_band1.hdr" 
    # infile = open(reffile, 'r')
    # lines = infile.readlines()
    # infile.close()   
    lines = ['ENVI\n',
             'description = {h10v04_1km_dt12_band1}\n',
             'samples = 1200\n',
             'lines = 1200\n',
             'bands = 1\n',
             'header offset = 0\n',
             'file type = ENVI Standard\n',
             'data type = 12\n',
             'interleave = bsq\n',
             'byte order = 0\n',
             'x start = 1\n',
             'y start = 1\n',
             'map info = {Sinusoidal, 1, 1, -8895604.158132, 5559752.598833, 926.625433138333, 926.625433139167}\n',
             'projection info = {16, 6371007.181, 0.000000, 0.0, 0.0, D_Unknown, Sinusoidal, units=Meters}\n',
             'coordinate system string = {PROJCS["Sinusoidal",GEOGCS["GCS_Unknown",DATUM["D_Unknown",SPHEROID["S_Unknown",6371007.181,0.0]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],UNIT["Meter",1.0]]}\n',
             'data ignore value = 32767\n'
             ]
    ## can also try readhdr function in raster_read_plot.py
    RESO = 926.625433138333
    NS = 1200
    NL = 1200
    UPLEFT = np.array([-8895604.158132, 5559752.598833])

    # PREPARE THE INPUTS
    if ns==None or nl==None:
        if reso==None: 
            sys.exit("Must at least assign a value to ns/nl or reso")
        elif reso!=250 and reso!=500 and reso!=1000:
            sys.exit("reso should be one of 250, 500, 1000")
        else:
            nscale = 1000/reso
            reso = RESO/nscale
            ns = NS * nscale
            nl = NL * nscale
    else:
        if ns!=nl or (ns!=1200 and ns!=2400 and ns!=4800):
            sys.exit("ns and nl must be equal, should be one of 1200, 2400, 4800")
        nscale = ns/NS
        if reso!=None and reso!=1000/nscale:
            warnings.warn("reso is not matching with ns/nl and changed to %d" % 1000/nscale)
        reso = RESO/nscale
        
    #### PROCESS
    for file in files:
        fileout = file+".hdr"
        if not overwrite and os.path.isfile(fileout):  continue
        print("writing file hdr %s" % fileout)
        
        outfile = open(fileout, "w")
        tmp = re.match(".*_(h(\d{2})v(\d{2}))", file)
        # tmp = re.match("(MCD12Q1)_LC_500m_byte_(h(\d{2})v(\d{2}))", file)
        if tmp:
            tile = tmp.group(1)  
            htile = int(tmp.group(2))
            vtile = int(tmp.group(3))
            upleft = UPLEFT+np.array([(htile-10)*NS*RESO, (4-vtile)*NL*RESO])
        else :
            sys.exit("Please check file %s" % file)
          

        newlines = lines
        for line in newlines:
            if re.match("description = ", line, re.IGNORECASE):   
                outfile.write('description = {%s}\n' % file)
            elif re.match('samples =', line, re.IGNORECASE):
                outfile.write("samples = %d\n" % ns )
            elif re.match('lines =', line, re.IGNORECASE):
                outfile.write("lines = %d\n" % nl )
            elif re.match('bands =', line, re.IGNORECASE):
                outfile.write("bands = %d\n" % nb )
            elif re.match('data type', line, re.IGNORECASE):
                outfile.write("data type = %s\n" % dtype)
            elif re.match('x start', line, re.IGNORECASE):
                outfile.write("x start = %d\n" % xstart)
            elif re.match('y start', line, re.IGNORECASE):
                outfile.write("y start = %d\n" % ystart)
            elif re.match('interleave', line, re.IGNORECASE):
                outfile.write('interleave = bip\n')
            elif re.match('map info =', line, re.IGNORECASE):
                outfile.write('map info = {Sinusoidal, 1, 1, %f, %f, %f, %f}\n' % (upleft[0], upleft[1], reso, reso))
            elif re.match('data ignore value = ', line, re.IGNORECASE):
                outfile.write('data ignore value = %s\n' % FILL)
            else:
                outfile.write(line)
        outfile.close()



def readbinary(file, nl, ns, nb, datatype, asdatatype = None, FILL=None):
    """
    Date and Author: Nov. 11th, 2021 by Jianmin Wang
    Aim:
        read binary file to numpy array
    :param file: the file to read
    :param nl: number of lines in the input file
    :param ns: number of samples in the input file
    :param nb: number of bands
    :param datatype: float, np.int16, np.uint8, np.float32 
    :param asdatatype: convert the default datatype to this type
    :param FILL: 
    """
    
    
    res = np.fromfile(file, dtype = datatype).reshape((nl, ns, nb))
    if asdatatype != None:
        res = res.astype(asdatatype)
        if FILL != None:
            res[res==FILL] = np.nan
    return(res)






def readgdal(file, bands, asdatatype = None, FILL=None):
    """
    Date and Author: Nov. 11th, 2021 by Jianmin Wang
    Aim:
        read raster file to numpy array
    :param bands: bands to read, start from 1. if set as -1, read all the bands
    :param datatype: float, np.int16, np.uint8, np.float32 
    :param asdatatype: convert the default datatype to this type
    :param FILL: 
    """
    
    ds = gdal.Open(file)
    nb = ds.RasterCount
    nl = ds.RasterYSize
    ns = ds.RasterXSize
    # dtype = gdal.GetDataTypeName(ds.GetRasterBand(1).DataType)
    
    if isinstance(bands, list):
        for band in bands:
            if band<1 or band >nb:
                print("Please set up bands with in the range [1, %d]" % nb)
                print("Otherwise use -1 to read all bands")
                sys.exit(1)
    elif isinstance(bands, int):
        if bands == -1:
            bands = range(1, nb+1) 
        elif bands >=1 and bands <=nb:
            bands = [bands]
        else:
            print("Please set up bands with in the range [1, %d]" % nb)
            print("Otherwise use -1 to read all bands")
            sys.exit(1)
    else:
        print("Please set up bands as integers with in the range [1, %d]" % nb)
        print("Otherwise use -1 to read all bands")
        sys.exit(1)
    
    nbands = len(bands)
    for i, band in enumerate(bands):
        tmp = ds.GetRasterBand(band).ReadAsArray()
        if i==0: res = np.zeros((nl, ns, nbands),  dtype=tmp.dtype)
        if nbands > 1: res[:,:,i] = tmp
        else: res = tmp
        
    ds = None        
    
    
    if asdatatype != None:
        res = res.astype(asdatatype)
        if FILL != None:
            res[res==FILL] = np.nan
    return(res)






def readhdr(file):
    """
    Date and Author: Nov. 11th, 2021 by Jianmin Wang
    Aim:
        read binary ENVI header file
        returns (rows, cols, minx, miny, maxx, maxy, xreso, yreso, src, nb, FILL)
    :param file: the file to read, can be a file with/o hdr HDR 
    """
    
    filehdr = None
    #https://docs.python.org/3/library/re.html#search-vs-match
    if re.search("(\.hdr$|\.HDR$)", file, re.IGNORECASE):
        filehdr = file
        # print("a")
    else:
        if os.path.isfile(file+'.hdr'): 
            filehdr = file+'.hdr'
        elif os.path.isfile(file+'.HDR'): 
            filehdr = file+'.HDR'
            
    if filehdr is None:
        sys.exit("Please check the header file")
        
    infile = open(filehdr, 'r')
    lines = infile.readlines()
    infile.close()       
    # lines = lines[0:-2]
    
    
    cols = None
    rows = None
    minx = None
    miny=None
    maxx = None
    maxy = None
    xreso = None
    yreso = None
    src =None
    nb = None
    FILL = None
    for line in lines:
        if re.match("^samples", line, re.IGNORECASE):   
            cols = int(line.split("=")[1].strip())
        elif re.match("^lines", line, re.IGNORECASE):   
            rows = int(line.split("=")[1].strip())
        elif re.match("^bands", line, re.IGNORECASE):   
            nb = int(line.split("=")[1].strip())
        elif re.match("^data ignore value", line, re.IGNORECASE):   
            FILL = int(line.split("=")[1].strip())
        elif re.match("^map info", line, re.IGNORECASE):  
            tmp = re.search("\{(.*)\}", line).group(1).split(",")
            xreso = float(tmp[5].strip())
            yreso = -float(tmp[6].strip())
            minx = float(tmp[3].strip()) - (float(tmp[1].strip())-1)*xreso
            maxy = float(tmp[4].strip()) - (float(tmp[2].strip())-1)*yreso ### yreso <0
            maxx = minx + cols*xreso 
            miny = maxy + rows*yreso ###yreso is <0
            transform = [minx, xreso, 0, maxy, 0, yreso]  
        elif re.match("^coordinate system string", line, re.IGNORECASE): 
            # print("Here")
            src = re.search("\{(.*)\}", line).group(1)
            
         
    return(rows, cols, minx, miny, maxx, maxy, xreso, yreso, src, nb, FILL)




#### Save each band in a ENVI binary file to individual bands.
def savetif():
    import numpy as np
    from osgeo import gdal
    path = "/scratch/jianmin.wang/ARD_MODIS_LST/fused_fulltile_nokshift2/"
    pathout = "/scratch/jianmin.wang/ARD_MODIS_LST/fused_fulltile_nokshift2/tifs2/"
    tile = "016006"
    year = 2016
    
    
    #row, col, minx, miny, maxx,  maxy, resx, resy, src = get_extent("/scratch/jianmin.wang/ARD_MODIS_LST/fused_fulltile_nokshift2/z1LSTfuse_016006_2016001.tif")
    row, col, minx, miny, maxx,  maxy, resx, resy, src, Nb, FILL = JWraster.readhdr("/scratch/jianmin.wang/ARD_MODIS_LST/fused_fulltile_nokshift2/z1LSTfuse_016006_2016")
    year = 2017
    ftype = "pred"
    for ftype in ["npairs"]: #["pred", "qa", "fuse", "npairs"]:
        f = open("%sz1LST%s_%s_%s" % (path, ftype, tile, year))
        FILL = 32767
        Nb = 366
        if ftype == "qa":
            data = np.fromfile(f, dtype = np.ubyte).reshape((row, col, Nb))
            FILL = 255
        else: 
            if ftype == "npairs":
                Nb = 1
            else:
                Nb = 366
            data = np.fromfile(f, dtype = np.int16).reshape((row, col, Nb))     
        f.close()
        
        
        for ib in range(Nb):
            # ib = 0
            fileout = "%sz1LST%s_%s_%s%03d.tif" % (pathout, ftype, tile, year, ib+1)
            driver = gdal.GetDriverByName("GTiff")
            if ftype == "qa":
                ds_o = driver.Create(fileout, row, col, 1, gdal.GDT_Byte)
            else:
                ds_o = driver.Create(fileout, row, col, 1, gdal.GDT_Int16)
            ds_o.SetGeoTransform([minx, resx, 0, maxy, 0, resy])
            ds_o.SetProjection(src)
        
            ds_o.GetRasterBand(1).WriteArray(data[:,:,ib])
            ds_o.GetRasterBand(1).SetNoDataValue(FILL)
            ds_o = None



#### Mosaic subsets to a large binary file
def mosaic_binary(files, fileout, rowblock, colblock, nb, nl, ns, datatype):
    """
    Date and Author: Sep. 21th, 2021 by Jianmin Wang
    Aim:
        mosaic a bunch of binary files (wall-to-wall but no overlapps) 
        into a larger binary file.
    :param files: the files (without header file) to be merge
    :param fileout: the output file
    :param rowblock: number of rows in each block (ether a number or a list)
    :param colblock: number of cols in each block (ether a number or a list)
    :param nb: number of bands
    :param nl: number of lines in the output file
    :param ns: number of samples in the output file
    :param datatype: number of bytes (datatype)
    """
    
    
    if type(rowblock)==int: 
        if (nl%rowblock)!=0:  sys.exit("Please check rowblock, nl!!!")
        nll = nl//rowblock
        rowblock = [rowblock]*nll
    else:
        if sum(rowblock) != nl: sys.exit("Please check rowblock and nl")
        nll = len(rowblock)
    if type(colblock)==int:
        if (ns%colblock)!=0 :  sys.exit("Please check colblock and ns!!!")
        nss = ns//colblock
        colblock = [colblock]*nss
    else:
        if sum(colblock) != ns: sys.exit("Please check colblock and ns!!!")
        nss = len(colblock)
    
    ###Chekc if file size if correct
    index = 0 
    for i, row in enumerate(rowblock):
        for j, col in enumerate(colblock):
            file = files[index]
            fsize = os.path.getsize(file)
            if fsize!=row*col*nb*datatype:
                sys.exit("Please check the size (%d) of file %s which has a nl=%d, ns=%d, nb=%d, dtype=%d" % (fsize, file, row, col, nb, datatype))
            index= index+1
    
    
    ### Process
    fout = open(fileout, 'wb')
    k = 0
    for i in range(nll):
        tmpfiles = files[(i*nss):(i*nss+nss)]
        
        try:
            fins = [open(file, 'rb') for file in tmpfiles]
        finally:
            comp = True
            while comp:
                for fin, col in zip(fins, colblock):
                    bytes_write = fin.read(col*nb*datatype)
                    if not bytes_write:
                        comp=False
                        break
                    
                    fout.write(bytes_write)
                    
                # print("i=%d k = %d" % (i, k))
                if comp: 
                    k = k + 1 #why 2405 but not 2399
    
            for fin in fins:
                fin.close()   
    
    fout.close()        







def mosaic_gdal(inputs, output, fill, oformat="GTiff"):
    """
    Date and Author: Sep. 21th, 2021 by Jianmin Wang
    Aim:
        mosaic a bunch of gdal files into a large gdal file.
        the projection of first file would be used
        if the projection or resolutions are different among inputs.
        can deal with Overlapps
        
    :param inputs: the files to be mosaic
    :param output: the output file
    :param fill: FILL value
    :param oformat: output format (NOTE please use GTiff, otherwise FILL is not correctly assigned)
    """
    
    try:
        os.remove(output)
    except OSError:
        pass
    
    # merge_command=' '.join(['source /hunter/data1/wangj/anaconda3/etc/profile.d/conda.sh &&', ##!!!IMPORTANT to activate conda activate!!!
    #                     'conda activate python37 &&', 
    #                     'gdal_merge.py -n', str(fill), '-a_nodata', str(fill), '-init', str(fill), 
    #                     '-o', output, '-of', ft, ' '.join(inputs), 
    #                     '&& conda deactivate'])
    merge_command=' '.join(['module load python/3.7 &&', ##!!!IMPORTANT to activate conda activate!!!
                        'module load gdal &&', 
                        'gdal_merge.py -n', str(fill), '-a_nodata', str(fill), '-init', str(fill), 
                        '-o', output, '-of', oformat, ' '.join(inputs)])
    subprocess.call(merge_command,shell=True)
#    os.system(merge_command)















def mosaic_overlapp(inputs, output, nb, oformat):
    """
    Date and Author: Sep. 21th, 2021 by Jianmin Wang
    Aim:
        mosaic a bunch of gdal files into a large gdal file.
        requring all the files are in same projection and resolution.
        it is faster than mosaic_gdal. 
        it can deal with Overlapps (e.g., HLS data), compared to mosaic_binary
        
    :param inputs: the files with header file to be mosaic
    :param output: the output file
    :param nb: number of band
    :param oformat: output format (NOTE please use GTiff, otherwise FILL is not correctly assigned)
    """
    
    nin = len(inputs)
    ROW = np.zeros(nin)
    COL = np.zeros(nin)
    RESX = np.zeros(nin)
    RESY = np.zeros(nin)
    MINX = np.zeros(nin)
    MINY = np.zeros(nin)
    MAXX = np.zeros(nin)
    MAXY = np.zeros(nin)
    for i, file in enumerate(inputs):
        res = get_extent(file)
        
        ROW[i] = res['rows']
        COL[i] = res['cols']
        MINX[i] = res['minx']
        MINY[i] = res['miny']
        MAXX[i] = res['maxx']
        MAXY[i] = res['maxy']
        RESX[i] = res['xreso'] 
        RESY[i] = -res['yreso'] 
        src = res['proj']
        if i>0 and (RESX[i]!=RESX[i-1] or RESY[i]!=RESY[i-1]):
            sys.exit('The resolution is not the same for %s and %s' % (file, inputs[i-1]))
        
        
    resy = RESY[0]
    resx = RESX[0]        
    minx = min(MINX)
    miny = min(MINY)
    maxx = max(MAXX)
    maxy = max(MAXY)
    ns=math.floor((maxx-minx)/resx+0.5)
    nl=math.floor((miny-maxy)/resy+0.5)
    
    MINC=((MINX-minx)/resx+0.5).astype(np.int16)
    MAXC=((MAXX-minx)/resx+0.5).astype(np.int16)
    MINR=((MAXY-maxy)/resy+0.5).astype(np.int16)
    MAXR=((MINY-maxy)/resy+0.5).astype(np.int16)
    
    
    
    ### Write data into numpy array
    result = np.zeros([nl, ns, nb], dtype=np.int16)+FILL
    for i, file in enumerate(inputs):
        ds=gdal.Open(file)
        for ib in range(nb):
            data = np.array(ds.GetRasterBand(ib+1).ReadAsArray()) #ReadAsArray(0, 0, MAXC[i]-MINC[i], MAXR[i]-MINR[i])
            result[MINR[i]:MAXR[i], MINC[i]:MAXC[i], ib]=data
        ds = None
        
    
    #Save to ENVI file
    fileformat = oformat
    driver = gdal.GetDriverByName(fileformat)
    metadata = driver.GetMetadata()
    if metadata.get(gdal.DCAP_CREATE) != "YES":
        print("Driver {} does NOT support Create() method.".format(fileformat))
#    if gdal.DCAP_CREATE in metadata and metadata[gdal.DCAP_CREATE] == "YES":
#        print("Driver {} supports CreateCopy() method.".format(fileformat))
     
    ds = driver.Create(output, xsize=ns, ysize=nl, bands=nb, eType=gdal.GDT_Int16) #INTERLEAVE
    ds.SetGeoTransform([minx, resx, 0, maxy, 0, resy])
    ds.SetProjection(src)
    for ib in range(nb):
        ds.GetRasterBand(ib+1).WriteArray(result[:,:,ib])
        ds.GetRasterBand(ib+1).SetNoDataValue(FILL)
    ds=None












###
# file = "/scratch/jianmin.wang/realtime/Anci0_shift1_over336_BKy1_noehance/z1gri_VNP43_RT_h11v04_2017LSPref_p15_row0000t0099_col0000t2399"
# file = "/scratch/jianmin.wang/realtime/Anci0_shift3_over366_BK15y1_enhanceno_1styear_bestr/merged/z1gri_VNP43_RT_h11v04_2017gri_0"
def get_extent(file):
    """
    Date and Author: Sep. 21th, 2021 by Jianmin Wang
    Aim:
        mosaic a bunch of gdal files into a large gdal file.
        the projection of first file would be used
        if the projection or resolutions are different among inputs.
        
    :param file: the input file
    :output is a dictionary
    """
    
    
    ds = gdal.Open(file)

    rows = ds.RasterYSize
    cols = ds.RasterXSize
    gt = ds.GetGeoTransform()
    src = ds.GetProjection()
    minx = gt[0]
    maxy = gt[3]
    maxx = gt[0] + gt[1] * cols
    miny = gt[3] + gt[5] * rows
    ds = None
    
    result = {"rows": rows,
              "cols": cols, 
              "minx": minx, "miny":miny, "maxx":maxx, "maxy":maxy, 
              "xreso": gt[1], 
              "yreso": -gt[5], 
              "gt": gt, 
              "proj": src}
    return(result)
    # return (rows, cols, minx, miny, maxx, maxy, gt[1], gt[5], src)

####class HDF5 




# #### Convert h5 to geotiffs
# def h5_to_tif(filein):



def get_sinu_projection():
    """
    Date and Author: by Jianmin Wang
    Aim:
        get file projection
    :return: a string of projection information
    Note: currently, no projection information was attached with layers in .h5 file
    """
    # set projection, VIIRS is same with MODIS

    
    prj = 'PROJCS["Sinusoidal",\
                GEOGCS["GCS_Undefined",\
                    DATUM["Undefined",\
                        SPHEROID["User_Defined_Spheroid",6371007.181,0.0]],\
                    PRIMEM["Greenwich",0.0],\
                    UNIT["Degree",0.0174532925199433]],\
                PROJECTION["Sinusoidal"],\
                PARAMETER["False_Easting",0.0],\
                PARAMETER["False_Northing",0.0],\
                PARAMETER["Central_Meridian",0.0],\
                UNIT["Meter",1.0]]'
    return(prj)
### End function ###


    

def get_wgs84_projection():
    """
    Date and Author: by Jianmin Wang
    Aim:
        get file projection
    :return: a string of projection information
    """

    InSR = osr.SpatialReference()
    InSR.ImportFromEPSG(4326)       # WGS84/Geographic
    OutSR = osr.SpatialReference()