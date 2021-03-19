#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 10:38:51 2019

@author: nanfeng
"""

from ENVI import read_envi_header, write_envi_header
import numpy as np
import os, glob

def create_sca(in_dir):
    merge_dir = os.path.join(in_dir, "merge")
    vnir_dir = os.path.join(in_dir, "vnir")
    
    old_sca_file = glob.glob(os.path.join(merge_dir, "*_SCA"))
    if old_sca_file == []:
        print("Cannot find any SCA file!")
        return
    else:
        old_sca_file = old_sca_file[0]
        
    bg_mask_file = glob.glob(os.path.join(merge_dir, "*_BackgroundMask"))
    if bg_mask_file == []:
        print("Cannot find any BackgroundMask file!")
        return
    else:
        bg_mask_file = bg_mask_file[0]
        
    imugps_file = glob.glob(os.path.join(vnir_dir, "*_ProcessedIMUGPS.txt"))
    if imugps_file == []:
        print("Cannot find any Processed IMUGPS file!")
        return
    else:
        imugps_file = imugps_file[0]
        
    dem_image_file = glob.glob(os.path.join(merge_dir, "*_DEM"))
    if dem_image_file == []:
        print("Cannot find any DEM file!")
        return
    else:
        dem_image_file = dem_image_file[0]
        
    #Read old SCA image data.
    old_sca_header = read_envi_header(old_sca_file+'.hdr')
    old_sca_image = np.memmap(old_sca_file,
                              mode='r',
                              dtype='float32',
                              shape=(old_sca_header['bands'],
                                     old_sca_header['lines'],
                                     old_sca_header['samples']))
    view_zenith = np.copy(old_sca_image[0,:,:])
    view_azimuth = np.copy(old_sca_image[1,:,:])
    old_sca_image.flush()
    del old_sca_image
    
    #Read background mask image data.
    bg_mask_header = read_envi_header(bg_mask_file+'.hdr')
    bg_mask_image = np.memmap(bg_mask_file,
                              mode='r',
                              dtype='bool',
                              shape=(bg_mask_header['lines'],
                                     bg_mask_header['samples']))

    #Get sun azimuth.
    sun_azimuth = float(old_sca_header['sun azimuth'])

    #Recalculate sensor angles.
    relative_azimuth = sun_azimuth-view_azimuth
    index = relative_azimuth>0.0
    view_zenith[index] = -view_zenith[index]
    view_zenith = view_zenith*100
    view_zenith[bg_mask_image] = 9100
    view_zenith = np.int16(view_zenith)

    relative_azimuth[relative_azimuth<0] += 360.0
    relative_azimuth[relative_azimuth>180] = 360.0-relative_azimuth[relative_azimuth>180]
    relative_azimuth = relative_azimuth*10
    relative_azimuth[bg_mask_image] = -1
    relative_azimuth = np.int16(relative_azimuth)

    #Get DEM.
    dem_header = read_envi_header(dem_image_file+'.hdr')
    dem_image = np.memmap(dem_image_file,
                          mode='r',
                          dtype='float32',
                          shape=(dem_header['lines'],
                                 dem_header['samples']))
    avg_dem = dem_image[dem_image>0.0].mean()
    del dem_image

    #Get GPS and IMU
    imugps = np.loadtxt(imugps_file) # ID, X, Y, Z, R, P, H, R_Offset, P_Offset, H_Offset, Grid_Convergence
    abg_altitude = np.ones(relative_azimuth.shape, dtype='int16')*int(imugps[:,3].mean()-avg_dem)
    abg_altitude[bg_mask_image] = -1

    #Write angle data.
    new_sca_file = old_sca_file+"_New"
    fid = open(new_sca_file, 'wb')
    fid.write(view_zenith.tostring())
    fid.write(relative_azimuth.tostring())
    fid.write(abg_altitude.tostring())
    fid.close()
    del view_zenith, view_azimuth, relative_azimuth, abg_altitude, bg_mask_image
    
    old_sca_header['bands'] = 3
    old_sca_header['data type'] = 2
    old_sca_header['GPS long-lat-alt'] = [imugps[:,12].mean(), imugps[:,13].mean(), imugps[:,3].mean()]
    old_sca_header['heading[deg]'] = imugps[:,6].mean()
    old_sca_header['DEM height[m]'] = avg_dem
    old_sca_header['band names'].append('ALT')
    write_envi_header(new_sca_file+'.hdr', old_sca_header)

flight_dirs = ["Z:/townsenduser-rw/hyspex_processed/Hancock_ARS/20180706/HANCOCK_20180706_1200/HANCOCK-1200_20180706_13",
               "Z:/townsenduser-rw/hyspex_processed/Hancock_ARS/20180718/HANCOCK_20180718/HANCOCK_20180718_09"]
for flight_dir in flight_dirs:
    print("Flight: %s." %flight_dir)
    create_sca(flight_dir)
"""
#in_dir = "Z:/townsenduser-rw/HyspexPro/Output/Hancock_Potato/20190625-2"
#date_dirs = glob.glob(os.path.join(in_dir, "*"))
date_dirs = ["Z:/townsenduser-rw/HyspexPro/Output/Hancock_Potato/20190625-2",
             "Z:/townsenduser-rw/HyspexPro/Output/Hancock_Potato/20190708-2"]
for date_dir in date_dirs:
    flight_dirs = glob.glob(os.path.join(date_dir, "*"))
    for flight_dir in flight_dirs:
        print("Flight: %s." %flight_dir)
        create_sca(flight_dir)
"""
