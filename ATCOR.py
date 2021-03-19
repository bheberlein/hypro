import glob, os
from collections import OrderedDict
from datetime import datetime
from ENVI import read_envi_header, write_envi_header
import numpy as np

def make_atcor_inn(params, atcor_fn):
    try:
        fid = open(atcor_fn, 'w')
    except:
        raise IOError('Cannot open \n%s!' %(atcor_fn))
    for field, value in params.items():
        if field in ['sensor name',
                     'cali_file',
                     'sca_file',
                     'dem_file',
                     'slope_file',
                     'aspect_file',
                     'skyview_file',
                     'shadow_file',
                     'refl_atm_file',
                     'therm_atm_file',
                     ]:
            fid.write(value+'\n')
        elif field == 'cali_fn':
            fid.write(value+'\n')
        else:
            fid.write('%s %s\n' %(value, field))
    fid.close()
    
def search_file(the_dir, keyword):
    file = glob.glob(os.path.join(the_dir, keyword))
    if len(file)==0:
        print('Cannot find any files in %s with the keyword %s.' %(the_dir, keyword))
        return None
    elif len(file)>1:
        print('Mulitple files found in %s with the keyword %s.' %(the_dir, keyword))
        return None
    else:
        return file[0]
    
def process_sca_file(raw_sca_image_file, processed_sca_image_file):
    # Read raw scan angle data.
    raw_sca_header = read_envi_header(os.path.splitext(raw_sca_image_file)[0]+'.hdr')
    raw_sca_image = np.memmap(raw_sca_image_file,
                              mode='r',
                              dtype='float32',
                              shape=(raw_sca_header['bands'],
                                     raw_sca_header['lines'],
                                     raw_sca_header['samples']))

    # Initialize processed scan angle data.
    processed_sca_image = np.memmap(processed_sca_image_file,
                                    mode='w+',
                                    dtype='int16',
                                    shape=(raw_sca_header['bands'],
                                           raw_sca_header['lines'],
                                           raw_sca_header['samples']))
    # zenith*100
    angle = raw_sca_image[0,:,:]*100
    angle[angle<0] = 9100
    processed_sca_image[0,:,:] = angle.astype('int16')
    del angle
    
    # zenith*10
    angle = raw_sca_image[1,:,:]*10
    angle[angle<0] = -1
    processed_sca_image[1,:,:] = angle.astype('int16')
    del angle
    
    # Clear data
    processed_sca_image.flush()
    raw_sca_image.flush()
    del processed_sca_image, raw_sca_image
    
    # Write header
    raw_sca_header['description'] = 'SCA'
    raw_sca_header['data type'] = 2
    raw_sca_header['band names'] = ['Sensor Zenith [100*deg]', 'Sensor Azimuth [10*deg]']
    raw_sca_header['data ignore value'] = None
    write_envi_header(os.path.splitext(processed_sca_image_file)[0]+'.hdr', raw_sca_header)

def get_avg_elev(dem_image_file):
    """ Get the average elevation of a DEM image.
    Notes:
        Pixels with a negative elevation are excluded from averaging.
    Arguments:
        dem_image_file: str
            DEM image filename.
    Returns:
        avg_elev: float
            Average elevation.
    """
    import gdal
    
    ds = gdal.Open(dem_image_file, gdal.GA_ReadOnly)
    dem_image = ds.GetRasterBand(1).ReadAsArray()
    avg_elev = dem_image[dem_image>0.0].mean() # Negative values are ignored.
    ds = None
    del dem_image
    return avg_elev

def main(in_dir, atcor_dir, idl_file):
    cmd_file = os.path.join(in_dir, 'atcor_cmd.txt')
    fid = open(cmd_file, 'w')
    
    flight_dirs = glob.glob(os.path.join(in_dir, '*'))
    for flight_dir in flight_dirs:
        
        # Skip if the flight diretory is not a diretory.
        if not os.path.isdir(flight_dir):
            continue
        print(flight_dir)
        
        # Get the merge diretory.
        merge_dir = os.path.join(flight_dir, 'merge')
        if not os.path.isdir(merge_dir):
            continue
        
        # Get radiance and scan angles filenames.
        rdn_image_file = search_file(merge_dir, '*Rdn')
        raw_sca_image_file = search_file(merge_dir, '*SCA')
        dem_image_file = search_file(merge_dir, '*DEM')
        if (rdn_image_file is None) or (raw_sca_image_file is None) or (dem_image_file is None):
            continue
        if (not os.path.exists(rdn_image_file+'.hdr')) or (not os.path.exists(raw_sca_image_file+'.hdr')) or (not os.path.exists(dem_image_file+'.hdr')):
            continue
        
        # Read radiance header file data.
        rdn_header = read_envi_header(os.path.splitext(rdn_image_file)[0]+'.hdr')

        # Read sca header file data.
        raw_sca_header = read_envi_header(os.path.splitext(raw_sca_image_file)[0]+'.hdr')

        # Get ATCOR parameters.
        params = OrderedDict()

        # 1: date
        when = datetime.strptime(rdn_header['acquisition time'], '%Y-%m-%dT%H:%M:%S.%f')
        params['Date (dd-mm-year)'] = when.date().strftime('%d-%m-%Y')#dd-mm-yyyy
        del when

        # 2: reflectance scale factor
        params['Scale factor reflectance'] = '1.0'

        # 3: pixel size
        params['Pixel size [m]'] = rdn_header['map info'][5]
        
        # 4: sensor name
        params['sensor name'] = 'HyspexPro, HyspexPro,'

        # 5: gain setting
        params['gainset'] = '1.0'

        # 6: calibration file
        cali_fn = glob.glob(os.path.join(atcor_dir, 'sensor', 'HyspexPro', '*.cal'))
        params['cali_fn'] = cali_fn[0]
        
        # 7: scan angle file
        processed_sca_image_file = os.path.join(merge_dir, os.path.basename(raw_sca_image_file)+'_INT16')
        process_sca_file(raw_sca_image_file, processed_sca_image_file)
        params['sca_fn'] = processed_sca_image_file

        # 8: emissivity and dem unit
        params['iemiss, dem_unit  (0=[m], 1=[dm], 2=[cm])'] = '0.98  0'

        # 9: dem file
        params['dem_fn'] = dem_image_file
        
        # 10: slope
        params['slope_fn'] = ''

        # 11: aspect
        params['aspect_fn'] = ''

        # 12: skyview
        params['skyview_fn'] = ''

        # 13: shadow
        params['shadow_fn'] = ''

        # 14: atm filename
        imugps_file = glob.glob(os.path.join(flight_dir, 'vnir', '*ProcessedIMUGPS.txt'))[0]
        altitude = np.loadtxt(imugps_file)[:,3].mean()
        altitudes = np.array([1000, 2000, 3000, 4000, 5000])
        atm_altitude = altitudes[(np.abs(altitudes-altitude)).argmin()]
        basename = ('h%05d' %atm_altitude) +'_wv20_rura.atm'
        params['refl_atm_fn'] = os.path.join(atcor_dir, 'atm_lib', 'HyspexPro', basename)
        params['refl_atm_fn'] = basename
        del basename, altitudes, atm_altitude

        # 15: atmospheric LUT filename in the thermal region
        params['therm_atm_fn'] = ''

        # 16: adjacency range[km]
        avg_elev = get_avg_elev(dem_image_file)
        params['Adjacency range [km]'] = '%.3f' %((altitude-avg_elev)*0.1/1000.0) #meter to km
        
        # 17: visibility
        params['Visibility [km]'] = '30'

        # 18: mean ground elevation [km]
        params['Mean ground elevation [km]'] = '%.3f' %(avg_elev/1000.0)
        del avg_elev
        
        # 19: solar zenith and azimuth angle
        params['Solar zenith, azimuth [degree]'] = '%s  %s' %(raw_sca_header['sun zenith'], raw_sca_header['sun azimuth'])

        # 20: Flight altitude [km], heading[deg]
        altitude = altitude/1000.0 #meter to km
        heading = np.loadtxt(imugps_file)[:,6].mean()
        params['Flight altitude [km], heading[deg]'] = '%.1f  %.1f' %(altitude, heading)

        # npref, iwaterwv, ihaze, iwat_shd, ksolflux, ishadow, icl_shadow
        """
        Notes:
            npref=1: variable visibility
            iwaterwv=1: water vapor correction using bands in the 940 nm region
            iwaterwv=2: water vapor correction using bands in the 1130 nm region
            ihaze=0: no haze correction
            iwat_shd=0: water pixels are exclueded from de-shadowing (land)
            ksolflux=0: file with value added channels not calculated
            ishadow=0; no DEM cast shadow map is used
            icl_shadow=0: no cloud/building shadow correction
        """
        params['npref, iwaterwv, ihaze, iwat_shd, ksolflux, ishadow, icl_shadow'] = '1 2 0 0 0 0 0'

        # itriang, ratio_red_NIR , ratio_blue_red
        """
        Notes:
            itriang: average vis. index of reference areas is employed for non-reference pixels
            ratio_red_NIR: red-to-NIR ratio
            ratio_blu_red: blue-to-red ratio
        """
        params['itriang, ratio_red_NIR , ratio_blue_red'] = '0 0.1 0.5'

        # BRDF parameters
        params['ibrdf, beta_thr, thr_g'] = '0 65.0 0.250'

        # LAI parameters
        params['lai_model, a0_vi, a1_vi, a2_vi'] = '1 0.820 0.780 0.600'

        # FPAR model: C, A, B
        params['FPAR model: C, A, B'] = '0.900 0.950 0.380'

        # air temperature and emissivity
        params['t_air, emiss_air   (flat terrain)'] = '20.0 0.830'

        # t_air, z0_ref, tgradient, p_wv, zh_pwv  (rugged terrain)
        params['t_air, z0_ref, tgradient, p_wv, zh_pwv  (rugged terrain)'] = '20.0 0.50 0.65 15.0 6.3'

        # ihot_mask, ihot_dynr
        params['ihot_mask, ihot_dynr'] = '2 2'

        # iclshad_mask, thr_shad, phi_unscl_max, phi_scl_min, istretch_type
        params['iclshad_mask, thr_shad, phi_unscl_max, phi_scl_min, istretch_type'] = '2 -999.00 -999.00 0.08 1'

        # bands for  820 nm water vapor retrieval
        params['bands for 820 nm water vapor retrieval'] = '0 0 0 0 0 0'

        # bands for 1130 nm water vapor retrieval
        params['bands for 1130 nm water vapor retrieval'] = '212 212 218 223 240 240'
        
        # bands for thermal water vapor retrieval (night data)
        params['bands for thermal water vapor retrieval (night data)'] = '0 0 0 0'

        # surface esmissitive (adjusted NEM, channel with Tmax)
        params['surface esmissitive (adjusted NEM, channel with Tmax)'] = '0.990 0.980 0.975 0.965'

        # iwv_model (1=wv without regression, 2=regression with several bands)
        params['iwv_model (1=wv without regression, 2=regression with several bands)'] = '2'

        # irrad0  (1=solar flux at ground, 2=+surface reflected radiance)
        params['irrad0  (1=solar flux at ground, 2=+surface reflected radiance)'] = '0'

        atcor_file = os.path.join(merge_dir, os.path.basename(rdn_image_file)+'.inn')
        make_atcor_inn(params, atcor_file)
        fid.write('%s -rt="%s/bin/atcor4.sav" -args "%s" F "%s_atm"\n' %(idl_file, atcor_dir, rdn_image_file, rdn_image_file))
    fid.close()
    
atcor_dir = r'S:\nanfeng\rese\atcor_4'
idl_file = r'S:/nanfeng/rese/IDL84/bin/bin.x86_64/idlrt.exe'
in_dir = r'Z:\townsenduser-rw\Potato\Hyspex_Ref'
main(in_dir, atcor_dir, idl_file)
