# P3 for ICON
## Important Remarks
- needs ComIn 0.4.0 and the corresponding icon branch https://gitlab.dkrz.de/icon/icon/-/tree/comin-0.4.0
- developed and tested only on Horeka cpuonly nodes with Intel oneAPI 2025.1 compilers and Open MPI 5.0
- not tested on Levante or other clusters and not developed for GPU offloading
- developed for a single ICON domain (ICON-LAM, Torus-Grid or Global) and only offline nesting is possible
- as the branch name suggests this is based on P3 v5.4 but will probably be updated to the lastest P3 v5.5 at some time

## Detailed Installation Steps
- clone the icon branch comin-0.4.0 and this repo outside of it and checkout the branch `icon-p3v5.4`
- unpack the P3 lookup tables:
```
cd P3-microphysics-icon/lookup_tables
gunzip *.gz
```
- load Horeka modules needed for building:
```
module purge
module add compiler/intel/2025.1_llvm mpi/openmpi/5.0 lib/netcdf/4.9-serial lib/hdf5/1.14-serial lib/netcdf-fortran/4.6-intel-2025.1_llvm lib/eccodes_nonetcdf/2.32.0 numlib/mkl/2025.1 libfyaml
```
- copy and run special icon-comin config file, then build ICON (includes incomplete building of ComIn):
```
cp P3-microphysics-icon/interfaces/icon/icon-config/hk.cpu.intel-2025-openmpi-5.0_comin_debug my-icon-comin-0.4.0/config/kit
cd my-icon-comin-0.4.0
chmod 755 config/kit/hk.cpu.intel-2025-openmpi-5.0_comin_debug
./config/kit/hk.cpu.intel-2025-openmpi-5.0_comin_debug
make -j8
```
- build ComIn again manually, this time all of it and then export the absolute build path into `ComIn_DIR`:
```
cd externals/comin
rm -rf build
mkdir build
cd build
cmake ..
make
export ComIn_DIR="$(pwd)"
```
- check in `P3-microphysics-icon/interfaces/icon/src/CMakeLists.txt` the compiler options you want (debug or fast) and accordingly leave one of the two target_compile_options() blocks uncommented
- now you can build the P3 plugin outside of but linked to ComIn:
```
cd to P3-microphysics-icon/interfaces/icon
mkdir build
cd build
cmake ../src
make
```
- if this was successfull, this created the dynamic libraries like `libp3plugin.so` which are linked to the ComIn libraries in ICON
- have a look at the example runscript included here and adjust the paths to you P3-microphysics-icon:
```
path_to_plugin="/path/to/myrepo/P3-microphysics-icon/interfaces/icon/build"
lookup_tables_path="/path/to/myrepo/P3-microphysics-icon/lookup_tables"
```
- the following lines in the example are essential to copy over the compiled p3plugin libraries to the workspace dir:
```
(after cd $EXPDIR)
mkdir p3plugin
cp $path_to_plugin/libp3*.so p3plugin
export LD_LIBRARY_PATH="$(pwd)/p3plugin:$LD_LIBRARY_PATH"
```
- and the following namelist comin_nml with these entries needs to exist in you `NAMELIST_NWP`:
```
&comin_nml
plugin_list(1)%name           = 'p3plugin'
plugin_list(1)%plugin_library = 'libp3plugin.so'
plugin_list(1)%comm           = 'p3_comm'
/
```
- when you run the runscript, the libraries in its current state get copied to the ouput dir, this way while the run is in the waiting queue you can make changes to the plugin source code and recompile it. When the submitted simulation finally runs it still uses the original plugin version that was copied over before you changes, this results in a faster workflow while developing in the p3plugin or P3 itself.

## Usage
### Main namelist settings
Apart from the `&comin_nml` group you need to set `inwp_gscp` to -1 (external microphysics) and turn off the saturation adjustment since the condensation is done within P3. Using 3-moments for the ice phase is recommended for all cases since it does barely need more computing time. The number of ice categories has a larger effect on runtime, 2 is a good intermediate setting. When using only 1 larger ice particles will be numerically diluted/shrunk in regions of ice nucleation because the gamma distribution does only represent one size peak that cannot well describe very small and large ice particles.
```
&nwp_phy_nml
 inwp_gscp                   =                         -1         ! -1: external microphysics via ComIn
 inwp_satad                  =                          0         ! saturation adjustment off
 icalc_reff                  =                          0         ! currently only possible setting (no calculation)
 icpl_rad_reff               =                          0         ! currently only possible setting (no calculation)
/
&p3_nml
 n_icecat                    =                          2         ! number of free ice categories (P3 scheme)
 l3mom_ice                   =                      .TRUE.        ! use triple moment ice categories (P3 scheme)
 lliqfrac                    =                      .TRUE.        ! use predicted bulk liquid fraction (P3 scheme)
 itracer_ini                 =                          1         ! 0: dry initialization, 1: read 1M-scheme tracer, 3: read P3 tracer
 tracer_ini_filename         =       '$INILBCDIR/$ININAME'
 lookup_tables_path          =     '${lookup_tables_path}'
/
```
### P3 Initialization
Three options are available set by `itracer_ini` for the initialization of the P3 tracers:
- 0: initialize without clouds and precipitation, "dry" (only qv present)
- 1: initialize from 1M-scheme mass tracers qc, qr, qi, qs, if available also use qg
- 3: initialize from P3 tracers
The tracer_ini_filename file must be in netcdf4 format and also contain the model halve levels of the ini data (named `hhl` or `HHL`) as a 2D-field. The tracer fields can be 2D like e.g. (height, ncells) or 3D like e.g. (time, height, ncells) but the dimensions can be of another order e.g. (height. ncells, time). The fields are interpolated vertically from the ini data levels to the model run levels. This way a typical output file from a 1M-scheme or P3-scheme model run can be used to initialize P3 and therefore "offline" nesting is possible.
In the initialization with P3 tracers, if the ice category number between ini data and model run doesn't match additional categories are ignored or categories left empty. The same applies for differing settings of `l3mom_ice` or `lliqfrac`.

### P3 lateral boundary conditions for ICON-LAM setups
This is not handled at the moment, qc and qi LBC fields are read by ICON and nudged in the nudging zone but qi is not read by P3 and therefore only qc is being nudged in and out of the domain.

### Further settings in P3
P3 is run always with 2-moment cloud droplets (in theory it could run also with 1-moment cloud droplets), the CCN activation and maximum concentration can be set in `microphy_p3.F90` with `iparam_ccn`. The warm-phase autoconversion can be set with `iparam_rain` and the homogeneous freezing level with `t_hom_freeze` (no INP activation implemented currently in P3, only droplet immersion freezing after Bigg, 1953).

### Output variables
All precipitating ice (sum of all categories) updates into `snow_gsp_rate` and `snow_gsp`. Rain is updated into `rain_gsp_rate` and `rain_gsp` and the sum of both into `prec_gsp_rate`, `prec_gsp` and `tot_prec`. 3D equivalent reflectivity (model level) can be output with `ze_p3`. The P3 tracers and number-mean ice diameters can be output with (here for 2 categories):
```
 ml_varlist='qv','qc','qnc','qr','qnr','qi','rho',
 'qitot_1','qnitot_1'[,'qzitot_1'],'qirim_1','birim_1'[,'qiliq_1'],'dmean_i1',
 'qitot_2','qnitot_2'[,'qzitot_2'],'qirim_2','birim_2'[,'qiliq_2'], 'dmean_i2'
 ```
Here `qi` is the sum of all `qitot_x` whereas `qs` is not being used and will always be zero in the output.
