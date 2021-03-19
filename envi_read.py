""" ENVI reader tools.
"""
import numpy as np
import os

ext_list = ['', '.img', '.raw', '.hyspex', '.bsq', 'bil', 'bip']

""" ENVI data type
    Reference: https://www.harrisgeospatial.com/docs/ENVIHeaderFiles.html
    The type of data representation:
        1 = Byte: 8-bit unsigned integer
        2 = Integer: 16-bit signed integer
        3 = Long: 32-bit signed integer
        4 = Floating-point: 32-bit single-precision
        5 = Double-precision: 64-bit double-precision floating-point
        6 = Complex: Real-imaginary pair of single-precision floating-point
        9 = Double-precision complex: Real-imaginary pair of double precision floating-point
        12 = Unsigned integer: 16-bit
        13 = Unsigned long integer: 32-bit
        14 = 64-bit long integer (signed)
        15 = 64-bit unsigned long integer (unsigned)
"""
dtype_dict = {1:np.uint8,
              2:np.int16,
              3:np.int32,
              4:np.float32,
              5:np.float64,
              6:np.complex64,
              9:np.complex128,
              12:np.uint16,
              13:np.uint32,
              14:np.int64,
              15:np.uint64}

value_type = {'acquisition time': 'str',
              'band names': 'str_list',
              'bands': 'int', # required
              'bbl': 'float_list',
              'byte order': 'int', # required
              'class lookup': 'int_list',
              'class names': 'str_list',
              'classes': 'int',
              'cloud cover': 'float',
              'complex function': 'str',
              'coordinate system string': 'str',
              'data gain values': 'float_list',
              'data ignore value': 'float',
              'data offset values': 'float_list',
              'data reflectance gain values': 'float_list',
              'data reflectance offset values': 'float_list',
              'data type': 'int', # required
              'default bands': 'float_list',
              'default stretch': 'str',
              'dem band': 'int',
              'dem file': 'str',
              'description': 'str',
              'file type': 'str', # required
              'fwhm': 'float_list',
              'geo points': 'float_list',
              'header offset': 'int',
              'interleave': 'str', # required
              'lines': 'int', # required
              'map info': 'str_list',
              'pixel size': 'float_list',
              'projection info': 'str',
              'read procedures': 'str_list',
              'reflectance scale factor': 'float',
              'rpc info': 'str',
              'samples': 'int', # required
              'security tag': 'str',
              'sensor type': 'str',
              'solar irradiance': 'float',
              'spectra names': 'str_list',
              'sun azimuth': 'float',
              'sun elevation': 'float',
              'wavelength': 'float_list',
              'wavelength units': 'str',
              'x start': 'float',
              'y start': 'float',
              'z plot average': 'int_list',
              'z plot range': 'float_list',
              'z plot titles': 'str'}

def read_header(fn):
    header = dict()
    fid = open(fn, 'r')
    trans_tab = str.maketrans(dict.fromkeys('\n{}'))
    while 1:
        line = fid.readline()
        if '=' in line:
            key, value = line.split('=', 1)
            if ('{' in value) and ('}' not in value):
                while '}' not in line:
                    line = fid.readline()
                    if line.strip()[0] == ';':
                        continue
                    value += line
            key = key.strip()
            value = value.translate(trans_tab).strip()
            if key in value_type.keys():
                if value_type[key] == 'int':
                    value = int(value)
                elif value_type[key] == 'float':
                    value = float(value)
                elif value_type[key] == 'str_list':
                    value = list(map(str.strip, value.split(',')))
                elif value_type[key] == 'int_list':
                    value = list(map(int, value.split(',')))
                elif value_type[key] == 'float_list':
                    value = list(map(float, value.split(',')))
            header[key] = value
        if line == '':
            break
    fid.close()
    # fill unused fields with None
    for key in value_type:
        if key not in header:
            header[key] = None
    if (header['interleave'] is None) or (header['interleave'].lower() not in ['bip','bil','bsq']):
        raise IOError('Unknown ENVI interleave')
    header['interleave'] = header['interleave'].lower()
    if header['header offset'] is None:
        header['header offset'] = 0
    check_required_keys(header)
    return header

def check_required_keys(header):
    required_keys = ['bands',
                     'byte order',
                     'data type',
                     'file type',
                     'interleave',
                     'lines',
                     'samples']
    for key in required_keys:
        if key not in header:
            raise IOError('%s not in header.' %key)
    if ('header offset' not in header) or (header['header offset'] is None):
        header['header offset'] = 0

def read_line(image, line, interleave):
    if interleave == 'bip':
        data = image[line,:,:]
    elif interleave == 'bil':
        data = image[line,:,:]
    elif interleave == 'bsq':
        data = image[:,line,:]
    return data

def read_sample(image, sample, interleave):
    if interleave == 'bip':
        data = image[:,sample,:]
    elif interleave == 'bil':
        data = image[:,:,sample]
    elif interleave == 'bsq':
        data = image[:,:,sample]
    return data

def read_band(image, band, interleave):
    if interleave == 'bip':
        data = image[:,:,band]
    elif interleave == 'bil':
        data = image[:,band,:]
    elif interleave == 'bsq':
        data = image[band,:,:]
    return data

def read_pixel(image, line, sample, interleave):
    if interleave == 'bip':
        data = image[line,sample,:]
    elif interleave == 'bil':
        data = image[line,:,sample]
    elif interleave == 'bsq':
        data = image[:,line,sample]
    return data

class ENVIReader():
    def __init__(self, image_fn=None):
        self.__image_fn = None
        self.__header_fn = None
        self.__image = None
        self.__header = None
        self.__interleave = None
        self.__shape = None
        self.__lines = None
        self.__samples = None
        self.__bands = None
        self.__offset = None
        self.__wavelengths = None
        self.__fwhms = None
        self.__no_data = None
        self.__map_info = None
        self.__crs = None
        self.__ulXY = None
        self.__lrXY = None
        self.__psize = None
        self.__dtype = None
        if image_fn is not None:
            self.image_fn = image_fn

    @property
    def image_fn(self):
        return self.__image_fn

    @property
    def header_fn(self):
        return self.__header_fn

    @property
    def image(self):
        return self.__image

    @property
    def header(self):
        return self.__header

    @property
    def interleave(self):
        return self.__interleave

    @property
    def shape(self):
        return self.__shape

    @property
    def lines(self):
        return self.__lines

    @property
    def samples(self):
        return self.__samples

    @property
    def bands(self):
        return self.__bands

    @property
    def offset(self):
        return self.__offset

    @property
    def wavelengths(self):
        return self.__wavelengths

    @property
    def fwhms(self):
        return self.__fwhms

    @property
    def no_data(self):
        return self.__no_data

    @property
    def map_info(self):
        return self.__map_info

    @property
    def crs(self):
        return self.__crs

    @property
    def ulXY(self):
        return self.__ulXY

    @property
    def lrXY(self):
        return self.__lrXY

    @property
    def psize(self):
        return self.__psize

    @property
    def dtype(self):
        return self.__dtype

    @image_fn.setter
    def image_fn(self, fn):
        self.__image_fn = fn
        self.header_fn = os.path.splitext(fn)[0]+'.hdr'

    @header_fn.setter
    def header_fn(self, fn):
        # read header
        self.__header_fn = fn
        self.read_header(self.header_fn)
        # update image filename
#        image_fn = os.path.splitext(self.header_fn)[0]
#        for ext in ext_list:
#            if os.path.exists(image_fn+ext):
#                self.__image_fn = image_fn
#                return
#        raise IOError('Cannot find the image file of %s.' %self.header_fn)

    @header.setter
    def header(self, header):
        # update header
        check_required_keys(header)
        header['interleave'] = header['interleave'].lower()
        if header['interleave'] not in ['bip', 'bil', 'bsq']:
            raise IOError('Unknown ENVI interleave: %s' %header['interleave'])
        self.__header = header.copy()
        # update other attributes
        self.__interleave = self.header['interleave']
        if self.header['interleave'] == 'bip':
            self.__shape = (self.header['lines'], self.header['samples'], self.header['bands'])
        elif self.header['interleave'] == 'bil':
            self.__shape = (self.header['lines'], self.header['bands'], self.header['samples'])
        elif self.header['interleave'] == 'bsq':
            self.__shape = (self.header['bands'], self.header['lines'], self.header['samples'])

        self.__lines = self.header['lines']
        self.__samples = self.header['samples']
        self.__bands = self.header['bands']
        self.__dtype = dtype_dict[self.header['data type']]
        self.__wavelengths = self.header['wavelength']
        self.__fwhms = self.header['fwhm']
        self.__no_data = self.header['data ignore value']
        self.__crs = self.header['coordinate system string']

        if self.header['map info'] is not None:
            self.__map_info = self.header['map info']
            x_tie_point = float(self.header['map info'][1])
            y_tie_point = float(self.header['map info'][2])
            pixel_easting = float(self.header['map info'][3])
            pixel_northing = float(self.header['map info'][4])
            x_pixel_size = float(self.header['map info'][5])
            y_pixel_size = float(self.header['map info'][6])
            self.__ulXY = (pixel_easting-(x_tie_point-1)*x_pixel_size, pixel_northing+(y_tie_point-1)*y_pixel_size)
            self.__lrXY = (self.__ulXY[0]+(self.__samples-1)*x_pixel_size, self.__ulXY[1]-(self.__lines-1)*y_pixel_size)
            self.__psize = (x_pixel_size, y_pixel_size)
        self.__offset = self.header['header offset']

    def read_header(self, fn):
        self.header = read_header(fn)

    def load_image(self):
        if self.header is None:
            raise IOError('Read header file first.')
        self.__image = np.memmap(self.image_fn, dtype=self.dtype, mode='r', shape=self.shape, offset=self.offset)

    def read_line(self, line):
        return read_line(self.image, line, self.interleave)

    def read_sample(self, sample):
        return read_sample(self.image, sample, self.interleave)

    def read_band(self, band):
        return read_band(self.image, band, self.interleave)

    def read_pixel(self, line, sample):
        return read_pixel(self.image, line, sample, self.interleave)

    def close(self):
        del self.__image
