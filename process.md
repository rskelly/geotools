# Requirements and Steps for Processing Hyperspectral Data

## Required Data

1) Raw irradiance from the downwelling radiometer (the Flame.)
2) A band map for the spectrometer (the Nano.)
3) Reflectance spectra from the handheld spectrometer (the ASD.)
4) Radiance image from the nano (a raster.)
5) The **imu_gps.txt** file corresponding to the radiance image.
6) The **frameIndex_*n*.txt** file corresponding to the radiance image.

## Data Formats

Each tabular data source must be formatted a certain way. 

### Irradiance Data

Flame output is tab-delimited by default, with several lines of header data that is not formatted as a table. The table
need not be modified for processing, but the `contrem` program must be told where the data begins by setting the `Header Row`
field to 13 and the `First Data Column` to 2. The `Date` and `Time` are at `0` and `1`, respectively.

The format will be changed after convolution. See below.

### Handheld Spectrometer Data

ASD data are produced as individual binary files, which can be converted to text representations, then collated into a
spreadsheet. Because there are different ways of doing this, only the format required for consumption by `convolve` is
given.

By whatever means, convert the ASD data into a table with at least the following columns:

label | date | timestamp | *b<sub>0</sub>* ... *b<sub>n</sub>*
------|------|-----------|------------------------------------
lw2 | 2018-08-28 11:20:47.119000 | 1535480447119 | *v<sub>0</sub> ... v<sub>n</sub>*

Here, *b<sub>0</sub>* ... *b<sub>n</sub>* are the wavelengths, one column for each. The label column refers to the
target for that row and is required to identify the target at a later step.

In `convolve`, set the `Header Row` to `0`, the `First Data Column` to `3` and the `Date` and `Time` to `1` and `2`, respectively.

The format will be changed after convolution. See below.

### Convolved Output

Each tabular dataset (irradiance; handheld reflelctance) must be convolved to match the airborne spectrometer. Use the [convolve](https://github.com/rskelly/contrem/wiki/convolve) program to do this.

The format of convolved data is standardized. The columns are:

date | timestamp | *b<sub>0</sub>* ... *b<sub>n</sub>*
-----|-----------|------------
2018-08-28 11:20:47.119000 | 1535480447119 | *v<sub>0</sub> ... v<sub>n</sub>*

where *b<sub>0</sub>* ... *b<sub>n</sub>* are the band wavelengths, one column for each.

If there are extra columns in the input data, they can be copied and placed into this output manually. For example, the `label` field in the ASD table should be transferred. (The rows will be in the same order as the input.)

## Processing Steps

1) [Convolve](https://github.com/rskelly/contrem/wiki/convolve) the handheld spectrometer (ASD) data to match the airborne spectrometer.
2) [Convolve](https://github.com/rskelly/contrem/wiki/convolve) the irradiance (Flame) data to match the airborne spectrometer.
3) Compute the apparent [reflectance](https://github.com/rskelly/contrem/wiki/reflectance) using the **imu_gps.txt** file, the **frameIndex_*n*.txt** file, the radiance raster and the convolved irradiance.
4) Run the [refl_regress.py](https://github.com/rskelly/contrem/wiki/refl_regress) to produce the coefficients.

