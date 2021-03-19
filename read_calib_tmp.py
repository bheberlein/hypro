from Radiometry import make_radio_cali_file_Hyspex
from scipy    import constants

    # Read meatadata from the raw Hyspex image.
    header = dict()
    try:
        fid = open(dn_image_file, 'rb')
    except:
        logger.error('Cannot read calibration data from %s.' %dn_image_file)
    header['word'] = np.fromfile(fid, dtype=np.int8, count=8)
    header['hdrSize'] = np.fromfile(fid, dtype=np.int32, count=1)[0]
    header['serialNumber'] = np.fromfile(fid, dtype=np.int32, count=1)[0]
    header['configFile'] = np.fromfile(fid, dtype=np.int8, count=200)
    header['settingFile'] = np.fromfile(fid, dtype=np.int8, count=120)

    header['scalingFactor'] = np.fromfile(fid, dtype=np.float64, count=1)[0]
    header['electronics'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['comsettingsElectronics'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['comportElectronics'] = np.fromfile(fid, dtype=np.int8, count=44)
    header['fanSpeed'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]

    header['backTemperature'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['Pback'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['Iback'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['Dback'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['comport'] = np.fromfile(fid, dtype=np.int8, count=64)

    header['detectstring'] = np.fromfile(fid, dtype=np.int8, count=200)
    header['sensor'] = np.fromfile(fid, dtype=np.int8, count=176)
    header['temperature_end'] = np.fromfile(fid, dtype=np.float64, count=1)[0]
    header['temperature_start'] = np.fromfile(fid, dtype=np.float64, count=1)[0]
    header['temperature_calibration'] = np.fromfile(fid, dtype=np.float64, count=1)[0]

    header['framegrabber'] = np.fromfile(fid, dtype=np.int8, count=200)
    header['ID'] = np.fromfile(fid, dtype=np.int8, count=200)
    header['supplier'] = np.fromfile(fid, dtype=np.int8, count=200)
    header['leftGain'] = np.fromfile(fid, dtype=np.int8, count=32)
    header['rightGain'] = np.fromfile(fid, dtype=np.int8, count=32)

    header['comment'] = np.fromfile(fid, dtype=np.int8, count=200)
    header['backgroundFile'] = np.fromfile(fid, dtype=np.int8, count=200)
    header['recordHD'] = np.fromfile(fid, dtype=np.int8, count=1)
    header['unknownPOINTER'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['serverIndex'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]

    header['comsettings'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['numberOfBackground'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['spectralSize'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['spatialSize'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['binning'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]

    header['detected'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['integrationTime'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['frameperiod'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['defaultR'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['defaultG'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]

    header['defaultB'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['bitshift'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['temperatureOffset'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['shutter'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['backgroundPresent'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]

    header['power'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['current'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['bias'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['bandwidth'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['vin'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]

    header['vref'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['sensorVin'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['sensorVref'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['coolingTemperature'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['windowStart'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]

    header['windowStop'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['readoutTime'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['p'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['i'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['d'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]

    header['numberOfFrames'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['nobp'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['dw'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['EQ'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['lens'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]

    header['FOVexp'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['scanningMode'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['calibAvailable'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['numberOfAvg'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['SF'] = np.fromfile(fid, dtype=np.float64, count=1)[0]

    header['apertureSize'] = np.fromfile(fid, dtype=np.float64, count=1)[0]
    header['pixelSizeX'] = np.fromfile(fid, dtype=np.float64, count=1)[0]
    header['pixelSizeY'] = np.fromfile(fid, dtype=np.float64, count=1)[0]
    header['temperature'] = np.fromfile(fid, dtype=np.float64, count=1)[0]
    header['maxFramerate'] = np.fromfile(fid, dtype=np.float64, count=1)[0]

    header['spectralCalibPOINTER'] = np.fromfile(fid, dtype=np.int32, count=1)[0]
    header['REPOINTER'] = np.fromfile(fid, dtype=np.int32, count=1)[0]
    header['QEPOINTER'] = np.fromfile(fid, dtype=np.int32, count=1)[0]
    header['backgroundPOINTER'] = np.fromfile(fid, dtype=np.int32, count=1)[0]
    header['badPixelsPOINTER'] = np.fromfile(fid, dtype=np.int32, count=1)[0]

    header['imageFormat'] = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    header['spectralVector'] = np.fromfile(fid, dtype=np.float64, count=header['spectralSize'])
    header['QE'] = np.fromfile(fid, dtype=np.float64, count=header['spectralSize'])
    header['RE'] = np.fromfile(fid, dtype=np.float64, count=header['spectralSize']*header['spatialSize'])
    header['RE'].shape = (header['spectralSize'], header['spatialSize'])
    header['background'] = np.fromfile(fid, dtype=np.float64, count=header['spectralSize']*header['spatialSize'])
    header['background'].shape = (header['spectralSize'], header['spatialSize'])

    header['badPixels'] = np.fromfile(fid, dtype=np.uint32, count=header['nobp'])
    if header['serialNumber']>=3000 and header['serialNumber']<=5000:
        header['backgroundLast'] = np.fromfile(fid, dtype=np.float64, count=header['spectralSize']*header['spatialSize'])
        header['backgroundLast'].shape = (header['spectralSize'], header['spatialSize'])

    fid.close()
