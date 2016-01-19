# tpcorr

## Dependencies

### Python packages

Most of these are pip installable. `galsim` is a bit of extra work.

 * numpy
 * scipy
 * matplotlib
 * h5py
 * astropy
 * specsim
 * bossdata
 * galsim

### BOSS data

`bossdata` handles interaction with BOSS data including automatic downloading of required data. One will also need to checkout the `speclog` repository from the public SDSS SVN.

## Usage

To calculate corrections for the offset fibers on a given plate-mjd and save them to an hdf5 file use:

```
calc_plate_correction --plate 6641
```

For full list of options use `calc_plate_correction --help`.