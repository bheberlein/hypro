import os, glob
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

in_dir = "Z:/townsenduser-rw/hyspex_raw/2019/20190819/DEL_20190819"


dn_image_files = glob.glob(os.path.join(in_dir, "*.hyspex"))
for dn_image_file in dn_image_files:
    basename = os.path.basename(dn_image_file)[:-len(".hyspex")]
    dn_header_file = os.path.splitext(dn_image_file)[0]+".hdr"
    raw_imugps_file = glob.glob(os.path.join(in_dir, "%s.txt" %basename))[0]
    acq_time = get_acquisition_time(dn_header_file, raw_imugps_file)
    print(acq_time)
