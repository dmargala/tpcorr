# tpcorr

## Dependencies

### Python packages

Most of these are pip installable. `galsim` might be a bit of extra work.

 * numpy
 * scipy
 * matplotlib
 * h5py
 * astropy
 * galsim
 * specsim
 * bossdata

### BOSS data

`bossdata` handles interaction with BOSS data including automatic downloading of required data. One will also need to checkout the platelist repository from the public SDSS SVN.

## Usage

To calculate corrections for the offset fibers on a given plate-mjd and save them to an hdf5 file use:

```
./calc_plate_correction.py --plate 6641 --mjd 56383
```

For more options use `./calc_plate_correction.py --help`.