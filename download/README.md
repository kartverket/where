# Necessary Data Files for Where

In general, Where is based on downloading data from original sources, and Where
comes with a mechanism for automatically downloading data files as they are
needed.


## Files Configuration

All files used by Where are listed in a files configuration. The files
configuration is available in the [`config` directory](../config/). Common files
are defined in `files.conf`, while pipeline specific files are defined in files
named `files_pipeline_....conf`.

Within the files configuration are fields

+ `url` - where can the file be downloaded from
+ `directory` - where will the file be stored on the local file system

as well as `origin` which gives a URL which typically should be possible to look
up in your browser for more information about the data file.


## Data Files that are Non-Downloadable

A few data files can not be downloaded from their original source. The main
reason for this is that they are generated for a given set of stations. Examples
of such files are instead made available here.

These files will also be automatically be downloaded by Where when they are
needed (they will be downloaded from the original Where-repository, not your
local copy).

See below for instructions on how to create your own versions of these files. If
you create your own versions of these files, see
[`files.conf`](../config/files.conf) for where they should be stored.


### Ocean Tide Loading

The file `ocean_loading_TPXO72_cmc_off` contains ocean tide loading
coefficients, and is created at the
[Ocean Tide Loading provider at Chalmers](http://holt.oso.chalmers.se/loading/)
by M.S. Bos and H.-G. Scherneck.

If you need to create your own file, either for other stations or with other
models, follow the instructions at the web site. Where uses its `site_id` field
to match stations with what is given in this file. For VLBI `site_id` is the
same as the CDP-number.


### Atmospheric Tide Loading

The files `s1_s2_def_cm.dat` and `s1_s2_def_ce.dat` contain atmospheric tide
loading coefficients, and are created at the
[S1 and S2 Atmospheric Tide Loading Calculator](http://geophy.uni.lu/ggfc-about/tide-loading-calculator.html)
by T. van Dam and R. Ray.

If you need to create your own file follow the instructions at the web site.


### Atmospheric Tide Loading - Center of Mass Corrections

The file `com.dat` contains center of mass corrections for atmospheric tide
loading. The values are available at the website of the
[S1 and S2 Atmospheric Tide Loading Calculator](http://geophy.uni.lu/ggfc-about/tide-loading-calculator.html)
by T. van Dam and R. Ray. However, we have only found the center of mass
corrections in
[PDF format](http://geophy.uni.lu/applications/atm1/download/com_table.pdf),
which is very inconvenient to use.

The `com.dat` file contains the values from this PDF-file in an easier to use
text file.
