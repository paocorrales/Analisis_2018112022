#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import wrf
import numpy as np
from netCDF4 import Dataset
import xarray as xr

pathANA = "/home/paola.corrales/datosmunin3/EXP/E2/ANA/20181122000000/analysis.ensmean"
#pathANA = "/home/paola.corrales/comGSIv3.7_EnKFv1.3/examples/test/gsi/wrfanl.2018112100"
ncfile = Dataset(pathANA)



pw = wrf.getvar(ncfile, "pw")

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

write_xarray_to_netcdf(pw, "pw_20181122000000_E2.nc", engine ="netcdf4")
