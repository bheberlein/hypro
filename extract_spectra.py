""" extract spectral statistics of potato
"""
import ogr
from envi_read import ENVIReader
import numpy as np
from PIL import Image, ImageDraw
import pandas as pd
import matplotlib.pyplot as plt

def get_shp_info(shp_fn):
    """
    Returns:
        shp_info: dict
            field: fieldorder
            value: N_Rate, Id, Points
    """
    shp_info = dict()
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.Open(shp_fn, 0)
    layer = ds.GetLayer()
    for feature in layer:
        geom = feature.GetGeometryRef()
        ring = geom.GetGeometryRef(0)
        points = ring.GetPoints()
        id = feature.GetField("id")
        shp_info[field_order] = dict()
        shp_info[field_order]['id'] = id
        layer = None
    ds = None
    return(shp_info)

def find_band_id(wavelengths, wav):
    wavelengths = np.array(wavelengths)
    if wavelengths[0]<1.0:
        wavelengths = wavelengths*1000.0
    band_id = np.argmin(np.abs(wavelengths-wav))
    return band_id

def get_spectral_stats(img_fn, shp_fn):
    ndvi_threshold = 0.60
    print(img_fn)
    # read header info
    img = ENVIReader(img_fn)
    img.load_image()
    # read shp info
    info = get_shp_info(shp_fn)
    # make mask
    nir = img.read_band(find_band_id(img.wavelengths, 840))
    red = img.read_band(find_band_id(img.wavelengths, 680))
    ndvi_mask = ((nir-red)/(nir+red+1e-18)<ndvi_threshold).astype(np.uint8)*255

    for key, value in info.items():
        info[key]['avg_spectra'] = []
        info[key]['std_spectra'] = []
        points = np.array(value['Points'])
        points[:,0] = (points[:,0]-img.ulXY[0])/img.psize[0]
        points[:,1] = (img.ulXY[1]-points[:,1])/img.psize[1]
        if (points[:,0].min()>=img.samples or
            points[:,0].max()<0 or
            points[:,1].min()>=img.lines or
            points[:,1].max()<0):
            continue
        for p in points:
            p[0] = max(0, p[0])
            p[0] = min(img.samples-1, p[0])
            p[1] = max(0, p[1])
            p[1] = min(img.lines-1, p[1])
        points = points.astype(np.int)
        sample0 = int(points[:,0].min())
        sample1 = int(points[:,0].max())
        line0 = int(points[:,1].min())
        line1 = int(points[:,1].max())

        mask1 = Image.new("L", (sample1-sample0, line1-line0), 255)
        points[:,0] = points[:,0]-sample0
        points[:,1] = points[:,1]-line0
        points = [tuple(p) for p in points]
        ImageDraw.Draw(mask1).polygon(points, outline=0, fill=0)
        mask1 = np.array(mask1).reshape(line1-line0, sample1-sample0)
        mask2 = ndvi_mask[line0:line1, sample0:sample1]
        mask = Image.fromarray(np.uint8(mask1+mask2))

        """
        plt.imshow(((nir-red)/(nir+red))[line0:line1, sample0:sample1])
        plt.show()
        plt.imshow(mask2)
        plt.show()
        plt.imshow(mask1)
        plt.show()
        plt.imshow(mask1+mask2)
        plt.show()
        mask.show()
        import sys
        sys.exit(1)
        """

        info[key]['mask'] = mask
        info[key]['row/col'] = (line0, line1, sample0, sample1)

    # loop over bands
    print('band: ', end='')
    for band in range(img.bands):
        if band%50 == 0:
            print('%03d, ' %band, end='')
        if not 'mask' in info[key]:
            continue
        data = img.read_band(band)
        for key, value in info.items():
            if img.header['bbl'] is not None and img.header['bbl'][band]==0:
                info[key]['avg_spectra'].append(np.nan)
                info[key]['std_spectra'].append(np.nan)
                continue

            line0, line1, sample0, sample1 = info[key]['row/col']
            mask = info[key]['mask']
            tmp = np.ma.array(data[line0:line1, sample0:sample1], mask=mask)
            info[key]['avg_spectra'].append(tmp.mean())
            info[key]['std_spectra'].append(tmp.std())
            """
            plt.imshow(tmp)
            plt.show()
            import sys
            sys.exit(1)
            """
    for key in info:
        info[key]['avg_spectra'] = np.array(info[key]['avg_spectra'])
        info[key]['std_spectra'] = np.array(info[key]['std_spectra'])
    print('%03d, Done!' %band)
    # wavelength
    wavelengths = np.array(img.wavelengths)#[np.array(img.header['bbl'])!=0]
    return info, wavelengths

import os, glob
shp_fn = 'polygon/filename'
in_img_dir = 'where/you/put/Hyspex/imagery'
out_img_dir = 'where/you/put/Hyspex/spectral/spreadsheet'
hdr_files = glob.glob(os.path.join(in_img_dir, '*.hdr'))

for hdr_file in hdr_files:
	img_fn = hdr_file[:-len('.hdr')]
	out_fn = os.path.join(out_img_dir, os.path.basename(img_fn)+'_spectra.xlsx')
	info, wavelengths = get_spectral_stats(img_fn, shp_fn)
	
	if wavelengths[0]==None:
		wavelengths = np.arange(len(info[list(info.keys())[0]]['avg_spectra']))
	Id = []
	avg_spectra = []
	std_spectra = []
	fill_value = np.full(len(wavelengths), np.nan)
	for value in info.values():
		Id.append(value['id'])
		if 'mask' in value:
			avg_spectra.append(value['avg_spectra'])
			std_spectra.append(value['std_spectra'])
		else:
			avg_spectra.append(fill_value)
			std_spectra.append(fill_value)
	avg_spectra = np.array(avg_spectra)
	std_spectra = np.array(std_spectra)

	df = pd.DataFrame(columns=['field_order', 'id']+list(wavelengths))
	df['field_order'] = list(info.keys())
	df['id'] = Id
	df[df.columns[3:]] = avg_spectra
	df.to_excel(out_fn)