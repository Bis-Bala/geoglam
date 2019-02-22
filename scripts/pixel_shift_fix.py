import numpy as np
import gdal
#import scipy.optimize as opt
#import scipy.ndimage
#from PIL import Image
import netCDF4
import json
import os
import sys
import datetime
import argparse
import time
import glob
import calendar


def build_data(fl, output_filename, hdf_file):

    # hdf_file = "/OSM/CDC/LPDAAC/source/lpdaac-tiles/c6/MCD43A4.006/2017.01.01/MCD43A4.A2017001.h31v11.006.2017014051817.hdf"
    
    o_ds = gdal.Open('HDF4_EOS:EOS_GRID:"{}":MOD_Grid_BRDF:Nadir_Reflectance_Band1'.format(hdf_file))
    # import pdb; pdb.set_trace()
    geot = o_ds.GetGeoTransform()
    # geot = (14455356.755667, 463.3127165275005, 0.0, -2223901.039333, 0.0, -463.3127165279167)
    with netCDF4.Dataset(fl) as src:
        with  netCDF4.Dataset(output_filename, "w") as dst:

            # copy all file data except for the excluded
            var_list = ['phot_veg', 'nphot_veg', 'bare_soil']
            xSize = src.dimensions['x'].size
            ySize = src.dimensions['y'].size
        # copy attributes from json file
            with open('nc_metadata.json') as data_file:
                attrs = json.load(data_file)
                for key in attrs:
                    setattr(dst, key, attrs[key])
            #for name in src.ncattrs():
            #    dst.setncattr(name, src.getncattr(name))
        # update glabal attribute comment 
        #    setattr(dst, "comment", "17 Feb 2019: implemented signed datatype(int8) in anomaly files (previously floating point\n19 Feb 2019: Corrected geolocation for all files")
        # copy dimensions
            for name, dimension in src.dimensions.items():
            
                #dst.createDimension(name, (len(dimension) if not dimension.isunlimited else None))
                dst.createDimension(name, (len(dimension)))
        
            for name, variable in src.variables.items():
                # Do the pixel shift logic here
                if name in var_list:
                    x = dst.createVariable(name, variable.datatype, variable.dimensions, fill_value=variable._FillValue, zlib=True, chunksizes=(1, 240, 240))
                else:
                    x = dst.createVariable(name, variable.datatype, variable.dimensions)
                           
                if name == 'x':
                    x[:] = np.arange(0, xSize) * geot[1] + geot[0] + geot[1]/2
                elif name == 'y':
                    x[:] = np.arange(0, ySize) * geot[5] + geot[3] + geot[5]/2
                elif name == 'sinusoidal':
                    dst.variables[name] = src.variables[name]
                else:
                    dst.variables[name][:] = src.variables[name][:]
                # update variable attributes
            
                for attrname in variable.ncattrs():
                    if attrname != '_FillValue':
                        setattr(x, attrname, getattr(variable, attrname))
                    if attrname == 'GeoTransform':
                        setattr(x, attrname, "{} {} {} {} {} {} ".format(*[geot[i] for i in range(6)]))
                       
    print (" output file is %s and hdf_file %s" % (output_filename, hdf_file))
           

def create_netcdf(input_files, output_dir, tile_pat):
    
    output_dir = output_dir + '/'
    hdf_pat = "/g/data2/u39/public/data/modis/lpdaac-tiles-c6/MCD43A4.006/2017.01.01/MCD43A4.A2017001."
    hdf_file = glob.glob(hdf_pat  + tile_pat +'*' + '*.hdf')
    hdf_file = hdf_file[0]
    
    for i, fl in enumerate(input_files):
        output_filename = fl.split('/')[-1]
        output_filename = output_dir + output_filename
        build_data(fl, output_filename, hdf_file)
    
    return None


if __name__ == "__main__":
    # product name can be of SRADIncomingAtmosphericLongwave, SRADNetLongwave, SRADNetRadiation, SRADOutgoingSurfaceLongwave, SRADShortwaveRatio, SRADSkyviewFactor, SRADTotalShortwaveSlopingSurf
    # call this as python create_netcdf.py input_dir="/data/projects/soil/data" product="SRADNetRadiation" output_dir="/data/projects/soil/output" metadata_dir="/data/projects/soil/metadata"
    parser = argparse.ArgumentParser(description="""Solar Radiation raster parser""")
    #parser.add_argument("-i", action='store', dest='input_dir', type=str, default='/g/data2/fr1/modis-fc/native_frac', help="Full path of input file.")
    #parser.add_argument("-i", action='store', dest='input_dir', type=str, default='/short/fr1/jpg599/geoglam_data/data/v310/tiles/monthly/cover', help="Full path of input file.")
    parser.add_argument("-i", action='store', dest='input_dir', type=str, default='/g/data2/tc43/temp', help="Full path of input file.")
    parser.add_argument("-o", action='store', dest='output_dir', type=str, default='/g/data2/fr1/modis-fc/native_frac/output/monthly/cover', help="Full path to destination of the netcdf files")
    parser.add_argument("-t", action='store', dest='tile_pat', type=str, help="date pattern like h00v09 h31v11 etc ")
    args = parser.parse_args()
    print ("inp dir %s" % args.input_dir)
    input_dir = args.input_dir
    output_dir = args.output_dir
    tile_pat = args.tile_pat

    input_files = glob.glob(input_dir +  '/'  + '*' + tile_pat + '*.nc')
    print (input_files)
    create_netcdf(input_files, output_dir, tile_pat)
