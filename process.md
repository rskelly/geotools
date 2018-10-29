# Requirements and Steps for Processing Hyperspectral Data

## Required Data

1) Raw irradiance from the Flame.
2) Band map for the Nano.
3) Radiance image from the nano.
4) ASD spectra.

## Data Formats

Each tabular data source must be formatted a certain way. 

### Flame Data

Flame data is tab-delimited by default, with several lines of header data that is not formatted as a table. The table
need not be modified for processing, but the `contrem` program must be told where the data begins by setting the `Header Row`
field to 13 and the `First Data Column` to 2. The `Date` and `Time` are at `0` and `1`, respectively.

The format will be changed after convolution. See below.

### ASD Data

ASD data are produced as individual binary files, which can be converted to text representations, then collated into a
spreadsheet. Because there are different ways of doing this, only the format required for consumption by this software is
given.


The format will be changed after convolution. See below.

### Convolved Output

The format of convolved data is standardized. The columns are:

date | timestamp | b_0 ... b_n
------------------------------

where `b_0 ... b_n` are the band wavelengths.

## Preprocessing:

Each tabular
