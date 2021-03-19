#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" 
Preparing ATM look up table and clipping DEM for each flightline. 
Contains Part 1 to Part 3 of HyspexPro workflow. 

"""

import argparse, glob, json, logging, os, re

from Geography        import get_map_crs, get_sun_angles
from IMUGPS           import prepare_imugps_Hyspex
from SensorModel      import determine_if_rotated, make_sensor_model
from DEM              import prepare_dem
from Boresight        import boresight_calibration
from GeoReferencing   import calculate_igm, calculate_sca, build_glt
from Figure           import plot_image_area, plot_angle_geometry, make_quicklook, plot_avg_rdn, plot_wvc_model, plot_smile_effect
from AtmLUT           import build_atm_lut, resample_atm_lut
from Radiometry       import make_radio_cali_file_Hyspex, dn2rdn_Hyspex, resample_rdn
from Classification   import pre_classification
from SmileEffect      import average_rdn, detect_smile_effect
from WVC              import build_wvc_model, estimate_wvc
from GeoRectification import orthorectify_rdn, orthorectify_dem, orthorectify_sca
from ImageMerging     import merge_dem_sca, merge_rdn
from VIS              import estimate_vis
from AtmCorr          import atm_corr_image

def get_flight_indices(config):
    """ Get Hyspex flight indices.
    Arguments:
        config: dict
            Configurations.
    Returns:
        flight_indicies: list of strs
            Flight indices.
    """

    sensor_name = config['Sensors'][list(config['Sensors'].keys())[0]]['id']
    flight_indices = []
    dn_image_files = glob.glob(os.path.join(config['Data']['input_dir'], '*%s*.hyspex' %sensor_name))
    for dn_image_file in dn_image_files:
        basename = os.path.basename(dn_image_file)
        span = re.search('%s' %sensor_name, basename).span()
        flight_indices.append(basename[:span[0]-1])
        del basename, span
    del dn_image_files, sensor_name

    return flight_indices

def create_flight_log(output_dir, log_file_basename):
    """ Create a Hyspex flight processing log.
    Arguments:
        output_dir: str
            Output directory.
        log_file_basename: str
            Log file basename.
    Returns:
        flight_log: logging object
            Flight log.
    """

    log_file = os.path.join(output_dir, '%s.log' %log_file_basename)
    logging.basicConfig(filename=log_file,
                        level=logging.DEBUG,
                        format="%(asctime)s %(funcName)25s: %(message)s",
                        datefmt='%Y-%m-%dT%H:%M:%S',
                        filemode='w')
    flight_log = logging.getLogger()

    return flight_log

def initialize_flight_dict(config, flight_index):
    """ Initialize a Hyspex flight dictionary.
    Arguments:
        config: dict
            User-defined configurations.
        flight_index: str
            Flight index.
    Returns:
        flight_dict: dict
            Flight dictionary.
    """

    # flight dictionary
    flight_dict = dict()

    # flight output directory
    flight_dict['output_dir'] = os.path.join(config['Data']['output_dir'], flight_index)
    if not os.path.exists(flight_dict['output_dir']):
        os.mkdir(flight_dict['output_dir'])

    # flight atmospheric lookup table directory
    flight_dict['atm_dir'] = os.path.join(flight_dict['output_dir'], 'atm')
    if not os.path.exists(flight_dict['atm_dir']):
        os.mkdir(flight_dict['atm_dir'])

    # flight merged image directory.
    flight_dict['merge_dir'] = os.path.join(flight_dict['output_dir'], 'merge')
    if not os.path.exists(flight_dict['merge_dir']):
        os.mkdir(flight_dict['merge_dir'])

    # atmospheric correction parameters
    flight_dict['atm_database_dir'] = config['Atmospheric_Correction']['atm_database_dir']
    flight_dict['atm_mode'] = config['Atmospheric_Correction']['atm_mode']
    flight_dict['vis_retrieval'] = config['Atmospheric_Correction']['vis_retrieval']
    flight_dict['wvc_retrieval'] = config['Atmospheric_Correction']['wvc_retrieval']

    # raw dem
    flight_dict['dem'] = config['DEM']

    # boresight offsets
    flight_dict['boresight_options'] = config['Geometric_Correction']['boresight']['options']

    # sensor dictionary
    flight_dict['sensors'] = dict()
    for sensor_index in config['Sensors'].keys():
        # sensor parameters
        sensor_dict = config['Sensors'][sensor_index].copy()

        # digital number (DN) image
        dn_image_file = search_file(config['Data']['input_dir'], '%s_%s*.hyspex' %(flight_index, sensor_dict['id']))
        sensor_dict['dn_image_file'] = dn_image_file

        # radiometric calibration parameters
        sensor_dict['setting_file'] = config['Radiometric_Calibration']['setting_file'][sensor_index]

        # geometric correction parameters
        sensor_dict['pixel_size'] = config['Geometric_Correction']['pixel_size'][sensor_index]
        if config['Geometric_Correction']['boresight']['gcp_file'] is None or config['Geometric_Correction']['boresight']['gcp_file'][sensor_index] is None:
            sensor_dict['gcp_file'] = None
        else:
            if dn_image_file in config['Geometric_Correction']['boresight']['gcp_file'][sensor_index].keys():
                sensor_dict['gcp_file'] = config['Geometric_Correction']['boresight']['gcp_file'][sensor_index][dn_image_file]
            else:
                sensor_dict['gcp_file'] = None

        # boresight offsets
        sensor_dict['boresight_offsets'] = config['Geometric_Correction']['boresight']['offsets'][sensor_index]

        # sensor model
        sensor_dict['sensor_model_file'] = config['Geometric_Correction']['sensor_model_file'][sensor_index]

        # imugps
        raw_imugps_file = search_file(config['Data']['input_dir'], '%s_%s*.txt' %(flight_index, sensor_dict['id']))
        sensor_dict['raw_imugps_file'] = raw_imugps_file

        # output directory
        sensor_dict['output_dir'] = os.path.join(flight_dict['output_dir'], sensor_index)
        if not os.path.exists(sensor_dict['output_dir']):
            os.mkdir(sensor_dict['output_dir'])

        flight_dict['sensors'][sensor_index] = sensor_dict
        del sensor_dict, dn_image_file, raw_imugps_file

    return flight_dict

def search_file(in_dir, keyword):
    """ Search a specific file with the keyword.
    Arguments:
        in_dir: str
            Hyspex data input directory.
        keyword: str
            Searching keyword.
    Returns:
        file: str
            Found filename.
    """

    file = glob.glob(os.path.join(in_dir, keyword))

    if len(file) == 0:
        raise IOError('Cannot find any file in %s with the keyword: %s.' %(in_dir, keyword))
    elif len(file) > 1:
        raise IOError('Multiple files are found in %s with the keyword: %s.' %(in_dir, keyword))
    else:
        return file[0]

def get_center_lon_lat(raw_imugps_file):
    """ Get Hyspex image center longitude and latitude.
    Arguments:
        raw_imugps_file: str
            Hyspex raw imugps filename.
    Returns:
        [lon, lat]: list of floats
            Image center longitude and latitude.
    """

    import numpy as np

    imugps = np.loadtxt(raw_imugps_file)
    lon, lat = imugps[:,1].mean(), imugps[:,2].mean()
    del imugps

    return [lon, lat]

def get_acquisition_time(dn_header_file, raw_imugps_file):
    """ Get Hyspex image acquistion time.
    Notes:
        (1) This code is adapted from Brendan Heberlein (bheberlein@wisc.edu).
    Arguments:
        header_file: str
            Hyspex DN image header filename.
        imugps_file: str
            Hyspex raw imugps filename.
    Returns:
        when: datetime object
            Image acquisition time.
    """

    from datetime import datetime, timedelta
    from ENVI     import read_envi_header

    import numpy  as np

    header = read_envi_header(dn_header_file)
    week_start = datetime.strptime(f"{header['acquisition date']} 00:00:00", "%Y-%m-%d %H:%M:%S")
    week_seconds = np.loadtxt(raw_imugps_file)[:,7].mean()
    epoch = datetime(1980, 1, 6, 0, 0)
    gps_week = (week_start-epoch).days//7
    time_elapsed = timedelta(days=gps_week*7, seconds=week_seconds)
    when = epoch+time_elapsed

    return when

def get_sun_earth_distance(when):
    """ Get sun-earth distance of a day.
    Arguments:
        when: datetime object
            Date and time.
    Returns:
        d: float
            Sun-Earth distance.
    """

    import numpy as np
    cur_dir = os.path.dirname(os.path.realpath(__file__))
    d = np.loadtxt(os.path.join(cur_dir, 'data/sun-earth-distance.dat'))
    doy = when.timetuple().tm_yday
    d = d[doy-1]
    
    return d

def HyspexPro(config_file):
    print(config_file)
    # Load configurations.
    config = json.load(open(config_file, 'r'))
    
    # Make an output directory.
    if not os.path.exists(config['Data']['output_dir']):
        os.mkdir(config['Data']['output_dir'])

    # Get flight indices.
    flight_indices = get_flight_indices(config)

    # Create a processing log file.
    flight_log = create_flight_log(config['Data']['output_dir'],
                                   os.path.splitext(os.path.basename(config['Data']['input_dir']))[0])

    # Process each flight.
    for flight_index in flight_indices:
        print(flight_index)
        #----------------------------------------Part 0----------------------------------------#
        flight_log.info('%sFlight: %s%s' %('='*20, flight_index, '='*20))

        # Initialize the flight dictionary.
        flight_dict = initialize_flight_dict(config, flight_index)

        #----------------------------------------Part 1----------------------------------------#
        flight_log.info('%sPART 1: Extract flight information.' %('-'*10))

        # center longitude and latitude
        sensor_dict = flight_dict['sensors'][list(flight_dict['sensors'].keys())[0]]
        flight_dict['center_lon_lat'] = get_center_lon_lat(sensor_dict['raw_imugps_file'])
        flight_log.info('Image center longitude and latitude [deg]: %.6f, %.6f' %(flight_dict['center_lon_lat'][0], flight_dict['center_lon_lat'][1]))

        # image acquisition time
        flight_dict['acquisition_time'] = get_acquisition_time(os.path.splitext(sensor_dict['dn_image_file'])[0]+'.hdr', sensor_dict['raw_imugps_file'])
        flight_log.info('Image acquisition time: %s' %flight_dict['acquisition_time'])

        # sun-earth distance
        flight_dict['sun_earth_distance'] = get_sun_earth_distance(flight_dict['acquisition_time'])

        # map coordinate system
        flight_dict['map_crs'] = get_map_crs(flight_dict['dem'], flight_dict['center_lon_lat'][0], flight_dict['center_lon_lat'][1])
        flight_log.info('Map coordinate system: %s' %flight_dict['map_crs'].GetAttrValue('projcs'))

        # sun zenith and azimuth angles
        flight_dict['sun_angles'] = get_sun_angles(flight_dict['center_lon_lat'][0], flight_dict['center_lon_lat'][1], flight_dict['acquisition_time'])
        flight_log.info('Sun zenith and azimuth angle [deg]: %.2f, %.2f' %(flight_dict['sun_angles'][0], flight_dict['sun_angles'][1]))
        del sensor_dict

        #----------------------------------------Part 2----------------------------------------#
        flight_log.info('%sPart 2: Do geo-referencings.' %('-'*10))
        for sensor_index, sensor_dict in flight_dict['sensors'].items():
            flight_log.info('Sensor: %s' %sensor_index)

            # Initialize.
            basename = os.path.basename(sensor_dict['dn_image_file'][:-len('_raw.hyspex')])

            # Process IMUGPS.
            flight_log.info('Prepare the IMU (Inertial Measurement Unit) and GPS (Global Positioning System) data.')
            sensor_dict['processed_imugps_file'] = os.path.join(sensor_dict['output_dir'], basename+'_ProcessedIMUGPS.txt')
            prepare_imugps_Hyspex(sensor_dict['processed_imugps_file'],
                                  sensor_dict['raw_imugps_file'],
                                  sensor_dict['boresight_offsets'],
                                  flight_dict['map_crs'],
                                  flight_dict['boresight_options'])

            # Generate sensor model.
            if sensor_dict['sensor_model_file'] is None:
                flight_log.info('Generate sensor model.')
                sensor_dict['sensor_model_file'] = os.path.join(sensor_dict['output_dir'], basename+'_SensorModel.txt')
                if_rotated = determine_if_rotated(sensor_dict['processed_imugps_file'])
                if if_rotated:
                    flight_log.info('The sensor is 180 degree rotated.')
                make_sensor_model(sensor_dict['sensor_model_file'],
                                  sensor_dict['fov'],
                                  sensor_dict['ifov'][1],
                                  sensor_dict['samples'],
                                  if_rotated)

            # Process DEM.
            flight_log.info('Prepare the DEM (Digital Elevation Model) data.')
            sensor_dict['dem_image_file'] = os.path.join(sensor_dict['output_dir'], basename+'_DEM')
            prepare_dem(sensor_dict['dem_image_file'],
                        flight_dict['dem'],
                        sensor_dict['processed_imugps_file'],
                        sensor_dict['fov'],
                        flight_dict['map_crs'],
                        sensor_dict['pixel_size'])

            # Do boresighting if the gcp file is available.
            if sensor_dict['gcp_file'] is not None:
                flight_log.info('Do boresighting.')
                sensor_dict['boresight_file'] = os.path.join(sensor_dict['output_dir'], basename+'_Boresight')
                boresight_calibration(sensor_dict['boresight_file'],
                                      sensor_dict['gcp_file'],
                                      sensor_dict['processed_imugps_file'],
                                      sensor_dict['sensor_model_file'],
                                      sensor_dict['dem_image_file'],
                                      flight_dict['boresight_options'])

            # Build IGM.
            flight_log.info('Calculate the IGM (Input Geometry).')
            sensor_dict['igm_image_file'] = os.path.join(sensor_dict['output_dir'], basename+'_IGM')
            calculate_igm(sensor_dict['igm_image_file'],
                          sensor_dict['processed_imugps_file'],
                          sensor_dict['sensor_model_file'],
                          sensor_dict['dem_image_file'],
                          flight_dict['boresight_options'])

            # Create SCA.
            flight_log.info('Calculate the SCA (Scan Angle).')
            sensor_dict['raw_sca_image_file'] = os.path.join(sensor_dict['output_dir'], basename+'_RawSCA')
            calculate_sca(sensor_dict['raw_sca_image_file'],
                          sensor_dict['processed_imugps_file'],
                          sensor_dict['igm_image_file'],
                          flight_dict['sun_angles'])

            # Build GLT.
            flight_log.info('Build the GLT (Geographic Lookup Table).')
            sensor_dict['glt_image_file'] = os.path.join(sensor_dict['output_dir'], basename+'_GLT')
            build_glt(sensor_dict['glt_image_file'],
                      sensor_dict['igm_image_file'],
                      sensor_dict['pixel_size']/2.0,
                      flight_dict['map_crs'])

            # Plot image areas.
            flight_log.info('Plot the image area.')
            sensor_dict['image_area_figure_file'] = os.path.join(sensor_dict['output_dir'], basename+'_ImageArea.png')
            plot_image_area(sensor_dict['image_area_figure_file'],
                            sensor_dict['dem_image_file'],
                            sensor_dict['igm_image_file'],
                            sensor_dict['processed_imugps_file'])

            # Build angle geometries.
            flight_log.info('Plot the angle geometry.')
            sensor_dict['angle_geometry_figure_file'] = os.path.join(sensor_dict['output_dir'], basename+'_AngleGeometry.png')
            plot_angle_geometry(sensor_dict['angle_geometry_figure_file'],
                                sensor_dict['raw_sca_image_file'])

            del basename
        del sensor_index, sensor_dict

        #----------------------------------------Part 3----------------------------------------#
        flight_log.info('%sPart 3: Build an ALT (Atmospheric Lookup Table).' %('-'*10))
        flight_dict['raw_atm_lut_file'] = os.path.join(flight_dict['atm_dir'], '%s_RawALT' %flight_index)
        build_atm_lut(flight_dict)

    logging.shutdown()
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config_file")
    args = parser.parse_args()
    flight_dict = HyspexPro(args.config_file)