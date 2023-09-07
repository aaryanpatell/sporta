import ctypes
import sys
import os
import pandas as pd
import openpyxl

lib = ctypes.CDLL("/home/patela4/Documents/database/imager/imager/imager.so")

class Hdf5Header(ctypes.Structure):
    _fields_ = [
        ("distance", ctypes.c_double),
        ("wavelength", ctypes.c_double),
        ("beam_center", ctypes.c_double * 2),
        ("size", ctypes.c_double * 2),
        ("psize", ctypes.c_double),
        ("threshold", ctypes.c_double),
        ("radius_90_percent", ctypes.c_int),
        ("resolution", ctypes.c_double),
        ("num_spots", ctypes.c_int),
        ("intensity", ctypes.c_double),
        ("ice", ctypes.c_int)
    ]

lib.getdirs.restype = ctypes.POINTER(ctypes.c_char_p)
lib.getdirs.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]

lib.getrecentfile.restype = ctypes.c_float
lib.getrecentfile.argtypes = [ctypes.c_char_p]

lib.getResolution.restype = ctypes.c_float
lib.getDistance.restype = ctypes.c_double
lib.getWavelength.restype = ctypes.c_double
lib.getPsize.restype = ctypes.c_double
lib.getThreshold.restype = ctypes.c_double
lib.getRadius90Percent.restype = ctypes.c_int
lib.getNumSpots.restype = ctypes.c_int
lib.getIntensity.restype = ctypes.c_double
lib.getIce.restype = ctypes.c_int

resolution = []
distance = []
wavelength = []
threshold = []
radius_90_percent = []
num_spots = []
intensity = []
ice = []

def getdirs(currentdir, folder):
    count = ctypes.c_int()
    result = lib.getdirs(currentdir.encode(), folder.encode(), ctypes.byref(count))
    dirs = [result[i].decode() for i in range(count.value)]
    lib.free_memory(result)
    
    return dirs

args = sys.argv
argc = len(args)

c_args = (ctypes.c_char_p * (argc + 1))()
c_args[:-1] = [arg.encode('utf-8') for arg in args]  # Assign to all elements except the last one
c_args[-1] = None  # Null-terminate the array

lib.argparse(argc, c_args)

currentdir = os.getcwd()
folder = "data"
dirs = getdirs(currentdir, folder)
for d in dirs:
    print(d)
    header = Hdf5Header()
    lib.getrecentfile(d.encode())
    resolution.append(lib.getResolution(header))
    distance.append(lib.getDistance(header))
    wavelength.append(lib.getWavelength(header))
    threshold.append(lib.getThreshold(header))
    radius_90_percent.append(lib.getRadius90Percent(header))
    num_spots.append(lib.getNumSpots(header))
    intensity.append(lib.getIntensity(header))
    ice.append(lib.getIce(header))

# Your previous code ...

# Make sure all the lists have the same length
if len(resolution) == len(distance) == len(wavelength) == len(threshold) == len(radius_90_percent) == len(num_spots) == len(intensity) == len(ice):
    data_dict = {
        'Resolution': resolution,
        'Distance': distance,
        'Wavelength': wavelength,
        'Threshold': threshold,
        '# Spots': num_spots,
        'Intensity': intensity,
        'Ice': ice,
    }

    df = pd.DataFrame(data_dict, index=dirs)

    # Optionally, you can transpose the DataFrame if you prefer directories as columns and attributes as rows.
    # df = df.T

    # Now you can export the DataFrame to an Excel file
    df.to_excel('test.xlsx', sheet_name='Data Dump')

    print(df)
else:
    print("Error: Data lists have different lengths.")
