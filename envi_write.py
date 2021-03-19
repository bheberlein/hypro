"""ENVI writer tools.
"""
import numpy as np
import os

ext_list = ['', '.img', '.raw', '.hyspex']

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

fields = ['acquisition time',
          'band names',
          'bands', # required
          'bbl',
          'byte order', # required
          'class lookup',
          'class names',
          'classes',
          'cloud cover',
          'complex function',
          'coordinate system string',
          'data gain values',
          'data ignore value',
          'data offset values',
          'data reflectance gain values',
          'data reflectance offset values',
          'data type', # required
          'default bands',
          'default stretch',
          'dem band',
          'dem file',
          'description',
          'file type', # required
          'fwhm',
          'geo points',
          'header offset',
          'interleave', # required
          'lines', # required
          'map info',
          'pixel size',
          'projection info',
          'read procedures',
          'reflectance scale factor',
          'rpc info',
          'samples', # required
          'security tag',
          'sensor type',
          'solar irradiance',
          'spectra names',
          'sun azimuth',
          'sun elevation',
          'wavelength',
          'wavelength units',
          'x start',
          'y start',
          'z plot average',
          'z plot range',
          'z plot titles']

def empty_header():
    header = dict()
    for key in fields:
        header[key] = None
    return header

def write_header(header, header_fn):
    fid = open(header_fn, 'w')
    fid.write('ENVI\n')
    for key, value in header.items():
        if value is None:
            continue
        if (key in fields) and (type(value) is list):
            value = '{%s}' %(', '.join(list(map(str, value))))
        else:
            value = str(value)
        fid.write('%s = %s\n' %(key, value))
    fid.close()
    
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

class ENVIWriter():
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
        ext = os.path.splitext(fn)[1]
        if not ext == '.hdr':
            raise IOError('Invalid header file extension: %s.' %ext)
        self.__header_fn = fn
        self.__image_fn = os.path.splitext(self.header_fn)[0]

    @header.setter
    def header(self, header):
        # update header
        check_required_keys(header)
        self.__header = header.copy()
        if not any(list(header.values())):# if all values are None, return.
            return
        header['interleave'] = header['interleave'].lower()
        if header['interleave'] not in ['bip', 'bil', 'bsq']:
            raise IOError('Unknown ENVI interleave: %s' %header['interleave'])
        if 'header offset' not in header or header['header offset'] is None:
            header['header offset'] = 0
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
    
    def load_image(self):
        if None in [self.image_fn, self.dtype, self.shape, self.offset]:
            raise IOError('Invalid values for self.image_fn, self.dtype, self.shape, self.offset.')
        self.__image = np.memmap(self.image_fn, dtype=self.dtype, mode='w+', shape=self.shape, offset=self.offset)
        
    def write_header(self):
        write_header(self.header, self.header_fn)
        
    def write_line(self, line, data):
        if (not self.shape == self.image.shape) or (not self.dtype == self.image.dtype):
            raise IOError('Load image first.')
        if self.interleave == 'bip':
            self.image[line,:,:] = data
        elif self.interleave == 'bil':
            self.image[line,:,:] = data
        elif self.interleave == 'bsq':
            self.image[:,line,:] = data
    
    def write_sample(self, sample, data):
        if (not self.shape == self.image.shape) or (not self.dtype == self.image.dtype):
            raise IOError('Load image first.')
        if self.interleave == 'bip':
            self.image[:,sample,:] = data
        elif self.interleave == 'bil':
            self.image[:,:,sample] = data
        elif self.interleave == 'bsq':
            self.image[:,:,sample] = data
    
    def write_band(self, band, data):
        if (not self.shape == self.image.shape) or (not self.dtype == self.image.dtype):
            raise IOError('Load image first.')
        if self.interleave == 'bip':
            self.image[:,:,band] = data
        elif self.interleave == 'bil':
            self.image[:,band,:] = data
        elif self.interleave == 'bsq':
            self.image[band,:,:] = data
    
    def close(self):
        del self.__image
    