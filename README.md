# tpcorr

## Setup

### Installation:

```
python setup.py install
```

Developer mode:

```
python setup.py develop
```

### Requirements

The following packages are used by tpcorr and are installable via pip:

 * numpy
 * scipy
 * matplotlib
 * h5py
 * astropy
 * specsim
 * bossdata

### Optional Dependency

The following are optional:

 * galsim (more complex PSF models)

### BOSS data

bossdata handles interaction with BOSS data including automatic downloading of required data. 

The sdss speclog product is also required. Since this is not available via SAS, checkout a local copy of the repo:

```
svn co https://svn.sdss.org/public/data/sdss/speclog/trunk speclog
```

And remember to set the `BOSS_SPECLOG` env var so bossdata knows how to find it.

## Usage

To calculate corrections for the offset fibers on a given plate-mjd and save them to an hdf5 file use:

```
calc_plate_correction --plate 6641
```

For full list of options use `calc_plate_correction --help`.