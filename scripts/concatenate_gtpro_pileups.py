#!/usr/bin/env python3

import xarray as xr
import sys
from lib.util import info

if __name__ == "__main__":
    outpath = sys.argv[1]

    data = []
    for path in sys.argv[2:]:
        data.append(xr.open_dataarray(path))
        info(data[-1].sizes)
    data = xr.concat(data, "library_id", fill_value=0)
    info(data.sizes)

    info(f'Writing to {outpath}')
    data.to_dataset(name="tally").to_netcdf(
        outpath, encoding=dict(tally=dict(zlib=True, complevel=6))
    )
