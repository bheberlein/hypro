def linear_percent_stretch(raw_image):
    """ Do linear percent stretch.
    References:
        (1) https://www.harrisgeospatial.com/docs/BackgroundStretchTypes.html
    Arguments:
        raw_image: 2D array
            Raw image data.
    Returns:
        stretched_image: 2D array
            Percent_stretched image.
    """

    stretched_image = np.zeros(raw_image.shape, dtype='uint8')
    low = raw_image.min()
    high = raw_image.max()
    stretched_image = (raw_image.astype('float32')-low)/(high-low)*255

    return stretched_image

def find_nearest_band(wavelengths, wave):
    wavelengths = np.array(wavelengths)
    band = np.argmin(np.abs(wavelengths-wave))
    return band

def make_quickview(qickview_file, reflectance_image_file):
    import os
    # Read reflectance image.
    reflectance_ds = ENVIImage(reflectance_image_file)
    rgb_bands = [find_nearest_band(reflectance_ds.wavelengths, 680),
                    find_nearest_band(reflectance_ds.wavelengths, 550),
                    find_nearest_band(reflectance_ds.wavelengths, 450)]

    # Write RGB image
    driver = gdal.GetDriverByName('GTiff')
    quickview_ds = driver.Create(qickview_file, reflectance_ds.ncolumns, reflectance_ds.nrows, 3, gdal.GDT_Byte)

    for output_band, rgb_band in enumerate(rgb_bands):
        quickview_image = linear_percent_stretch(reflectance_ds.read_band(rgb_band))
        quickview_ds.GetRasterBand(output_band+1).WriteArray(quickview_image.astype("uint8"))
        del quickview_image
    quickview_ds = None
    reflectance_ds = None
