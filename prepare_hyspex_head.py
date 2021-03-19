import numpy as np
from ENVI import empty_envi_header, write_envi_header

def read_hyspex_header(dn_image_file):
    header = dict()
    fid = open(dn_image_file, 'rb')
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

    # Convert int8array to string.
    def from_int8array_to_string(int8_array):
        string = ''
        for int8_value in int8_array:
            if int8_value == 0:
                break
            string += chr(int8_value)
        return string
    for key in header:
        if header[key].dtype.name == 'int8':
            header[key] = from_int8array_to_string(header[key])
    return header

dn_image_file = r"Z:\townsenduser-rw\hyspex_raw\2019\20190711\CHEESEHEAD-2_20190711\Flightline_06\CHEESEHEAD-2_20190711_06_VNIR_1800_SN00840_FOVx2_raw.hyspex"
dn_header_file = r"Z:\townsenduser-rw\hyspex_raw\2019\20190711\CHEESEHEAD-2_20190711\Flightline_06\CHEESEHEAD-2_20190711_06_VNIR_1800_SN00840_FOVx2_raw.hdr"
metadata = read_hyspex_header(dn_image_file)
##print("Frameperiod = %s" %metadata["frameperiod"])
##print("Integration time = %s" %metadata["integrationTime"])
##print("Binning = %s" %metadata["binning"])
##print("Number of frames = %s" %metadata["numberOfFrames"])
##print("dw = %s" %metadata["dw"])
##print("EQ = %s" %metadata["EQ"])
##print("FOVexp = %s" %metadata["FOVexp"])
##print("Lens = %s" %metadata["lens"])
##print("NumberOfAvg = %s" %metadata["numberOfAvg"])
##print("CalibAvailable = %s" %metadata["calibAvailable"])
##print("Number of background = %s" %metadata["numberOfBackground"])
##print("Pixelsize x = %s" %metadata["pixelSizeX"])
##print("Pixelsize y = %s" %metadata["pixelSizeY"])
##print("ID = %s" %metadata["ID"])
##print("Comment = %s" %metadata["comment"])
##print("Serialnumber = %s" %metadata["serialNumber"])
##print("Scanningmode = %s" %metadata["scanningMode"])
##print("samples = %s" %metadata["spatialSize"])
##print("lines = %s" %metadata["numberOfFrames"])
##print("bands = %s" %metadata["spectralSize"])
##print("header offset = %s" %metadata["hdrSize"])
##print("acquisition date = Unknown")
##print("acquisition start time = Unknown")
##print("data type = 12")
##print("data ignore value = 65535")
##print("interleave = bil")
##print("default bands = {%s, %s, %s}" %(metadata["defaultR"], metadata["defaultG"], metadata["defaultB"]))
##print("byte order = 0")
##print("wavelength = {%s}" %metadata["spectralVector"])

header = empty_envi_header()
header["description"] = ["Frameperiod = %s" %metadata["frameperiod"],
                         "Integration time = %s" %metadata["integrationTime"],
                         "Binning = %s" %metadata["binning"],
                         "Number of frames = %s" %metadata["numberOfFrames"],
                         "dw = %s" %metadata["dw"],
                         "EQ = %s" %metadata["EQ"],
                         "FOVexp = %s" %metadata["FOVexp"],
                         "NumberOfAvg = %s" %metadata["numberOfAvg"],
                         "CalibAvailable = %s" %metadata["calibAvailable"],
                         "Number of background = %s" %metadata["numberOfBackground"],
                         "Pixelsize x = %s" %metadata["pixelSizeX"],
                         "Pixelsize y = %s" %metadata["pixelSizeY"],
                         "ID = %s" %metadata["ID"],
                         "Comment = %s" %metadata["comment"],
                         "Serialnumber = %s" %metadata["serialNumber"],
                         "Scanningmode = %s" %metadata["scanningMode"]]
header["description"] = "\n"+"\n".join(header["description"])
header["samples"] = metadata["spatialSize"]
header["lines"] = metadata["numberOfFrames"]
header["samples"] = metadata["spatialSize"]
header["bands"] = metadata["spectralSize"]
header["header offset"] = metadata["hdrSize"]
header["acquisition date"] = "2019-07-11"#Give an date here.
header["acquisition start time"] = "13:59:55" #Give a time here.
header["data type"] = 12
header["data ignore value"] = 65535
header["interleave"] = "bil"
header["default bands"] = [metadata["defaultR"], metadata["defaultG"], metadata["defaultB"]]
header["byte order"] = 0
header["wavelength"] = list(metadata["spectralVector"])
write_envi_header(dn_header_file, header)
