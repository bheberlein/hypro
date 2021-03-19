a = None
try:
    if min(None) < 100:
        print('aaa')
except:
    print('bbb')
##import numpy as np
##from ENVI import read_envi_header
##from scipy import interpolate
##
##header_file = 'Z:/townsenduser-rw/HyspexPro/Output/Cheesehead/CHEESE_20190626_tile_01/tmp/CHEESE_20190626_05_VNIR_1800_SN00840_FOVx2_IGM.hdr'
##image_file = 'Z:/townsenduser-rw/HyspexPro/Output/Cheesehead/CHEESE_20190626_tile_01/tmp/CHEESE_20190626_05_VNIR_1800_SN00840_FOVx2_IGM'
##
##header = read_envi_header(header_file)
##lines = header['lines']
##samples = header['samples']
##image = np.memmap(image_file,
##                  mode='r',
##                  dtype='float64',
##                  shape=(header['bands'],
##                         header['lines'],
##                         header['samples']))
##
##igm_image = np.copy(image)
### Inpterpolate IGM.
##nan_lines = []
##nonnan_lines = []
##for line in np.arange(lines):
##    print(line)
##    # Find Nan values.
##    nan_flag = np.isnan(igm_image[0, line, :])
##
##    # If all columns are Nan values;
##    if np.all(nan_flag):
##        del nan_flag
##        nan_lines.append(line)
##        continue
##
##    # If some columns are Nan values;
##    if np.any(nan_flag):
##        nan_samples = np.arange(samples)[nan_flag]
##        nonnan_samples = np.arange(samples)[~nan_flag]
##        f = interpolate.interp1d(nonnan_samples, igm_image[:, line, nonnan_samples], axis=1, fill_value="extrapolate")
##        igm_image[:, line, nan_samples] = f(nan_samples)
##        del nan_samples, nonnan_samples, f, nan_flag
##
##    nonnan_lines.append(line)
##    
##for nan_line in nan_lines:
##    index = np.argmin(np.abs(nonnan_lines-nan_line))
##    nonnan_line = nonnan_lines[index]
##    igm_image[:,nan_line,:] = igm_image[:,nonnan_line,:]
##    del index, nonnan_line
##
##igm_image_file = 'Z:/townsenduser-rw/HyspexPro/Output/Cheesehead/CHEESE_20190626_tile_01/CHEESE_20190626_05/vnir/CHEESE_20190626_05_VNIR_1800_SN00840_FOVx2_IGM'
##fid = open(igm_image_file, 'wb')
##fid.write(igm_image.tostring())
##fid.close()
##del igm_image
##image.close()
##import os
##the_path = os.path.dirname(os.path.realpath(__file__))
##print(the_path)
