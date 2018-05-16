# contrem
Application for applying continuum removal using convex hulls.

# Usage

## Direct

contrem can be run directly from the command line. Running it witout arguments will yield the following:

    Usage: contrem [options]
     -d A GDAL-readable data file containing spectral samples; can contain any number of bands >= 2.
     -r An ENVI ROI text file.
     -b A CSV file containing a mapping from wavelength to (1-based) band index.
     -w An integer giving the (0-based) column index in -b which contains wavelengths.
     -i An integer giving the (0-based) column index in -b which contains the band indices.
     -z If given, indicates the presence of a header in the band map that must be skipped.
     -o An output file template. This is a filename with no extension that will be modified as appropriate. Parent directories will be created.
     -l The minimum wavelength to consider.
     -h The maximum wavelength to consider.
     -s The size of the buffer. Default is 256. Larger buffers are possible, but one must consider that multiple buffers may be in memory at once.
     -t The number of threads to use. Default 2.
     -p By default, sample statistics are used. This flag forces the use of population statistics.
     -v The driver to use for output rasters. Defaults to ENVI, but any GDAL-writable format will do.
     -e File extension for raster files. Defaults to .dat for ENVI files.

## Batch

To run batches, run the run.py in the run folder. The run program takes a single argument, the path to a configuration file (CSV) containing configuration parameters.

The csv columns are, in order:

Field Name | Description
---------- | -----------
data_file | A raster file containing spectra. This is any GDAL-readable raster. If the raster has wavelengths embedded in the band metadata or description, these will be loaded as a band map. If they are not present, a separate band map is required. This and roi_file are mutually exclusive.
roi_file | An ENVI ROI file. A band map is required with this format. This and data_file are mutually exclusive.
band_map | A CSV file containing a mapping between the wavelength and band number. If this is given, it will override the band map contained in the raster file.
bm_header | This value indicates whether the band map file has a header that must be skipped. If 't' is given, the header is skipped, otherwise not. If a band map file is given, this is a required field.
bm_wl_col | An integer indicating the column where the wavelength is stored. The first column is zero. If a band map file is given, this is a required field.
bm_idx_col | An integer indicating the column where the band index is stored. Band indices start at 1, which is the first band. The first column is zero. If a band map file is given, this is a required field.
output_template | This is the template for output filenames. If the path has a directory part, it will be created. The output type and extension will be appended to the template to create a full filename. For example, the template "output/result" will result in the filename "output/result_ch.dat" for the convex hull intersection file and the output folder will be created if it does not exist.
low_wl | The lower bound of the range of wavelengths to process. If given, the first wavelength greater than or equal to this one will be used. If empty, the entire range is used. If this is given, high_wl must also be given.
high_wl | The upper bound of the range of wavelengths to process. If given, the last wavelength less than or equal to this one will be used. If empty, the entire range is used. If this is given, low_wl must also be given.
population | If 'f' then sample statistics (using n-1) are calculated, other wise population stats (using n) are computed.
driver | The output driver -- anything that can be produced by GDAL. "GTiff" and "ENVI" are the two most likely choices.
extension | The file extension. The driver will not choose this automatically. ENVI files generally use "" or ".dat", GTiff files use ".tif" 

# Install
Install the program by doing the usual:

## Linux
1) Checkout from git
2) mkdir build
3) cd build
4) cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
5) make
6) sudo make install

The program requires GDAL and GEOS so these must be installed first.

## OSX
1) ./INSTALL_OSX.sh

The program requires GDAL and GEOS so these must be installed first. The install script attempts this
using homebrew.

## Windows
1) Puke and cry
2) Goto 1
