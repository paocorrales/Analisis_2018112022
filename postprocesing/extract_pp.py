#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import wrf
import numpy as np
from netCDF4 import Dataset
import xarray as xr

path_ini = "/glade/scratch/jruiz/EXP/E4/FCST/20181122030000/wrfout_d01_2018-11-22_06:00:00"
path_end = "/glade/scratch/jruiz/EXP/E4/FCST/20181122030000/wrfout_d01_2018-11-23_12:00:00"
ncfile_ini = Dataset(path_ini)
ncfile_end = Dataset(path_end)

p = wrf.getvar(ncfile_ini, "RAINNC")
pp_ini = wrf.getvar(ncfile_ini, "RAINNC") + wrf.getvar(ncfile_ini, "RAINC") + wrf.getvar(ncfile_ini, "RAINSH")
pp_end = wrf.getvar(ncfile_end, "RAINNC") + wrf.getvar(ncfile_end, "RAINC") + wrf.getvar(ncfile_end, "RAINSH") 

pp = pp_end - pp_ini
pp = pp.assign_attrs(p.attrs)

def write_xarray_to_netcdf(xarray_array, output_path, mode='w', format='NETCDF4', group=None, engine=None,
                           encoding=None):
    """writes and xarray in a netcdf format outputfile
    Uses the xarray typical for wrf-python. The projection objects are transformed into strings
    to be able to use them as netcdf attributes
    :param xarray_array: xarray.DataArray
    :param output_path: str
    :param format: 'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT' or 'NETCDF3_CLASSIC'
                    default: 'NETCDF4'
    :param group: str, default None
    :param engine: 'netcdf4', 'scipy' or 'h5netcdf'
    :param encoding: dict, default: None
    """
    xarray_array_out = xarray_array.copy(deep=True)
    # coordinates are extracted from variable
    del xarray_array_out.attrs['coordinates']
    # wrf-python projection object cannot be processed
    xarray_array_out.attrs['projection'] = str(xarray_array_out.attrs['projection'])

    xarray_array_out.to_netcdf(path=output_path, mode=mode, format=format, group=group,
                           engine=engine,
                           encoding=encoding)

write_xarray_to_netcdf(pp, "pp_acum_E4_fcst20181122030000.nc", engine ="netcdf4")


