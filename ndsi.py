import ogr
from envi_read import ENVIReader
from envi_write import ENVIWriter
import numpy as np
from copy import deepcopy

def find_band_id(wavelengths, wav):
    wavelengths = np.array(wavelengths)
    if wavelengths[0]<1.0:
        wavelengths = wavelengths*1000.0
    band_id = np.argmin(np.abs(wavelengths-wav))
    return band_id

def get_spectral_stats(img_fn, shp_fn):
#    ndvi_threshold = 0.60
    print(img_fn)
    # read header info
    img = ENVIReader(img_fn)
    img.load_image()
    # read shp info
    info = get_shp_info(shp_fn)
    # make mask
#    nir = img.read_band(find_band_id(img.wavelengths, 840))
#    red = img.read_band(find_band_id(img.wavelengths, 680))
#    ndvi_mask = ((nir-red)/(nir+red+1e-18)<ndvi_threshold).astype(np.uint8)*255

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
#        mask2 = ndvi_mask[line0:line1, sample0:sample1]
        mask = Image.fromarray(np.uint8(mask1))

#        plt.imshow(mask)
#        plt.show()
#        import sys
#        sys.exit(1)

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

def main(in_img_fn, out_img_fn, wav1, wav2, ndvi_threshold=0.6):

    in_img = ENVIReader(in_img_fn)
    in_img.load_image()
    nir = in_img.read_band(find_band_id(in_img.wavelengths, 840))
    red = in_img.read_band(find_band_id(in_img.wavelengths, 680))
    ndvi_mask = (nir-red)/(nir+red+1e-18)<ndvi_threshold

    band1 = in_img.read_band(find_band_id(in_img.wavelengths, wav1))
    band2 = in_img.read_band(find_band_id(in_img.wavelengths, wav2))
    ndsi = (band1-band2)/(band1+band2+1e-18)
    ndsi[ndvi_mask] = -999.0
    out_img = ENVIWriter(out_img_fn)
    header = deepcopy(in_img.header)
    header['band names'] = np.nan
    header['bands'] = 1
    header['bbl'] = np.nan
    header['data ignore value'] = -999.0
    header['fwhm'] = np.nan
    header['wavelength'] = np.nan
    out_img.header = header
    out_img.load_image()
    out_img.write_header()
    out_img.write_band(0, ndsi)
    out_img.close()
    in_img.close()

import os, glob

in_img_dir = 'where/you/put/Hyspex/imagery'
out_img_dir = 'where/you/put/Hyspex/NDVI/imagery'
hdr_files = glob.glob(os.path.join(in_img_dir, '*.hdr'))
wav1, wav2 = 650, 790
ndvi_threshold = 0.6

for hdr_file in hdr_files:
	in_img_fn = hdr_file[:-len('.hdr')]
	out_img_fn = os.path.join(out_img_dir, os.path.basename(in_img_fn)+'_NDVI')
	main(in_img_fn, out_img_fn, wav1, wav2, ndvi_threshold=0.6)