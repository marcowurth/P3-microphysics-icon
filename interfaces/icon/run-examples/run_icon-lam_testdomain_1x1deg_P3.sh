
# ----------------------------------------------------------------------
# SETTINGS: Input/output directories
# ----------------------------------------------------------------------

ICONDIR="/home/hk-project-aci/nw5893/ICON/icon-comin-0.4.0"

# START DATE AND TIME OF THE SIMULATION
STARTDATE="2018-11-10T14:00:00Z"
ENDDATE="2018-11-10T14:15:00Z"

# output directory
EXPDIR=/hkfs/work/workspace/scratch/nw5893-ws5/test_2024122700/domain_1x1deg/output_P3_1Mgrpscheme_exp4_2C

# absolute path to initial conditions and lateral boundary forcing data
INILBCDIR=/hkfs/work/workspace/scratch/nw5893-ws5/test_2024122700/domain_1x1deg/ini_lbc
ININAME="ini_domain_1x1deg_20181110T140000Z_1Mgrpscheme.nc"
LBCNAMEBEG="lbc_domain_1x1deg_"
lbc_grid_filename="grid_lbc_domain_1x1deg.grid.nc"
map_files_dir="/home/hk-project-aci/nw5893/ICON/map_files"

# absolute path to radiation input directory
RADDIR=${ICONDIR}/externals/ecrad

# absolute path to external parameter directory
EXTPARDIR=/hkfs/work/workspace/scratch/nw5893-ws5/test_2024122700/domain_1x1deg/grid_extpar
extpar_filename="external_parameter_icon_domain_1x1deg_DOM01_tiles.nc"

# absolute path to grid files
GRIDDIR=/hkfs/work/workspace/scratch/nw5893-ws5/test_2024122700/domain_1x1deg/grid_extpar
dynamics_grid_filename="domain_1x1deg_DOM01.nc"

# absolute path to model binary, including the executable
MODEL=${ICONDIR}/bin/icon

path_to_plugin="/home/hk-project-aci/nw5893/ICON/p3-microphysics-5.4-icon/interfaces/icon/build"
lookup_tables_path="/home/hk-project-aci/nw5893/ICON/p3-microphysics-5.4-icon/lookup_tables"

# ----------------------------------------------------------------------
# copy input data: grids, external parameters
# ----------------------------------------------------------------------

mkdir -p ${EXPDIR}
script_name=$(basename -- "$0")
#cp $script_name $EXPDIR
cd $EXPDIR
#ulimit -n hard
#ulimit -s unlimited
#ulimit -c unlimited

# limited area grid
ln -sf ${GRIDDIR}/${dynamics_grid_filename} ${dynamics_grid_filename}

# lateral boundary grid
ln -sf ${INILBCDIR}/${lbc_grid_filename} ${lbc_grid_filename}

# external parameter
ln -sf ${EXTPARDIR}/${extpar_filename} ${extpar_filename}

# files needed for radiation
cp -r ${RADDIR}/data ecrad_data

# Dictionary for the mapping: DWD GRIB2 names <-> ICON internal names
ln -sf ${map_files_dir}/ana_varnames_map_file_ICON-EU.txt map_file.ana

# Dictionary for the mapping: GRIB2/Netcdf input names <-> ICON internal names
ln -sf ${map_files_dir}/dict.latbc_ICON-EU map_file.latbc

# configuration
#cp $HOME/.bashrc .
cp ${ICONDIR}/config.log .

mkdir p3plugin
cp $path_to_plugin/libp3*.so p3plugin
export LD_LIBRARY_PATH="$(pwd)/p3plugin:$LD_LIBRARY_PATH"


# ----------------------------------------------------------------------------
# create ICON master namelist
# ----------------------------------------------------------------------------

cat > icon_master.namelist << EOF

! master_nml: ----------------------------------------------------------------
&master_nml
 lrestart                   =                      .FALSE.        ! .TRUE.=current experiment is resumed
/

! master_model_nml: repeated for each model ----------------------------------
&master_model_nml
 model_type                  =                          1         ! identifies which component to run (atmosphere,ocean,...)
 model_name                  =                      "ATMO"        ! character string for naming this component.
 model_namelist_filename     =              "NAMELIST_NWP"        ! file name containing the model namelists
 model_min_rank              =                          1         ! start MPI rank for this model
 model_max_rank              =                      65536         ! end MPI rank for this model
 model_inc_rank              =                          1         ! stride of MPI ranks
/

! time_nml: specification of date and time------------------------------------
&time_nml
 ini_datetime_string = "$STARTDATE"
 end_datetime_string = "$ENDDATE"
/
EOF


# ----------------------------------------------------------------------
# model namelists
# ----------------------------------------------------------------------

cat > NAMELIST_NWP << EOF
&comin_nml
plugin_list(1)%name           = 'p3plugin'
plugin_list(1)%plugin_library = 'libp3plugin.so'
plugin_list(1)%comm           = 'p3_comm'
/

! parallel_nml: MPI parallelization -------------------------------------------
&parallel_nml
 nproma                      =                         16         ! loop chunk length
 p_test_run                  =                     .FALSE.        ! .TRUE. means verification run for MPI parallelization
 num_io_procs                =                          1         ! number of I/O processors
 num_restart_procs           =                          0         ! number of restart processors
 num_prefetch_proc           =                          1         ! number of processors for LBC prefetching
 iorder_sendrecv             =                          3         ! sequence of MPI send/receive calls
 use_omp_input               =                     .FALSE.        ! allows task parallelism for reading atmospheric input data
/

! run_nml: general switches ---------------------------------------------------
&run_nml
 ltestcase                   =                     .FALSE.        ! real case run
 num_lev                     =                        120         ! number of full levels (atm.) for each domain
 dtime                       =                          5.        ! timestep in seconds
 ldynamics                   =                      .TRUE.        ! compute adiabatic dynamic tendencies
 ltransport                  =                      .TRUE.        ! compute large-scale tracer transport
 iforcing                    =                          3         ! forcing of dynamics and transport by parameterized processes
 msg_level                   =                         10         ! detailed report during integration
 ltimer                      =                      .TRUE.        ! timer for monitoring the runtime of specific routines
 timers_level                =                         10         ! performance timer granularity
 check_uuid_gracefully       =                      .TRUE.        ! give only warnings for non-matching uuids
 output                      =                       "nml"        ! main switch for enabling/disabling components of the model output
/

! diffusion_nml: horizontal (numerical) diffusion ----------------------------
&diffusion_nml
 lsmag_3d                    =                     .FALSE.
 lhdiff_vn                   =                      .TRUE.        ! diffusion on the horizontal wind field
 lhdiff_temp                 =                      .TRUE.        ! diffusion on the temperature field
 lhdiff_w                    =                      .TRUE.        ! diffusion on the vertical wind field
 hdiff_order                 =                          5         ! order of nabla operator for diffusion
 itype_vn_diffu              =                          1         ! reconstruction method used for Smagorinsky diffusion
 itype_t_diffu               =                          2         ! discretization of temperature diffusion
 hdiff_smag_fac              =                          0.015     ! scaling factor for Smagorinsky diffusion
/

! dynamics_nml: dynamical core -----------------------------------------------
&dynamics_nml
 iequations                  =                          3         ! type of equations and prognostic variables
 lcoriolis                   =                      .TRUE.        ! Coriolis force
/

! extpar_nml: external data --------------------------------------------------
&extpar_nml
 itopo                       =                          1         ! topography (0:analytical)
 extpar_filename             =        '${extpar_filename}'        ! filename of external parameter input file
/

! initicon_nml: specify read-in of initial state ------------------------------
&initicon_nml
 init_mode                   =                          7         ! start from DWD fg with subsequent vertical remapping 
 lread_ana                   =                     .FALSE.        ! no analysis data will be read
 dwdfg_filename              =       "$INILBCDIR/$ININAME"        ! initial data filename
 ana_varnames_map_file       =              "map_file.ana"        ! dictionary mapping internal names onto GRIB2 shortNames
 ltile_coldstart             =                      .TRUE.        ! coldstart for surface tiles
 ltile_init                  =                     .FALSE.        ! set it to .TRUE. if FG data originate from run without tiles
/

! grid_nml: horizontal grid --------------------------------------------------
&grid_nml
 dynamics_grid_filename      =  '${dynamics_grid_filename}'       ! array of the grid filenames for the dycore
 lredgrid_phys               =                     .FALSE.        ! .true.=radiation is calculated on a reduced grid
 lfeedback                   =                     .FALSE.        ! specifies if feedback to parent grid is performed
 l_limited_area              =                      .TRUE.        ! .TRUE. performs limited area run
 ifeedback_type              =                          2         ! feedback type (incremental/relaxation-based)
 start_time                  =                          0.        ! Time when a nested domain starts to be active [s]
/

! gridref_nml: grid refinement settings --------------------------------------
&gridref_nml
 denom_diffu_v               =                        150.        ! denominator for lateral boundary diffusion of velocity
/

! interpol_nml: settings for internal interpolation methods ------------------
&interpol_nml
 nudge_zone_width            =                         10         ! width of lateral boundary nudging zone
 nudge_max_coeff             =                          0.075
 support_baryctr_intp        =                     .FALSE.        ! barycentric interpolation support for output
/

! io_nml: general switches for model I/O -------------------------------------
&io_nml
 itype_pres_msl              =                          5         ! method for computation of mean sea level pressure
 itype_rh                    =                          1         ! method for computation of relative humidity
 lmask_boundary              =                     .FALSE.        ! mask out interpolation zone in output
/

! limarea_nml: settings for limited area mode ---------------------------------
&limarea_nml
 itype_latbc                 =                          1         ! 1: time-dependent lateral boundary conditions
 dtime_latbc                 =                        900.        ! time difference between 2 consecutive boundary data
 latbc_boundary_grid         =     '${lbc_grid_filename}'         ! Grid file defining the lateral boundary
 latbc_path                  =             '${INILBCDIR}'         ! Absolute path to boundary data
 latbc_varnames_map_file     =           'map_file.latbc'
 latbc_filename              = '${LBCNAMEBEG}<y><m><d>T<h><min>00Z.nc'     ! boundary data inputfilename
 init_latbc_from_fg          =                     .TRUE.         ! .TRUE.: take lbc for initial time from first guess
/

! lnd_nml: land scheme switches -----------------------------------------------
&lnd_nml
 ntiles                      =                          1         ! number of tiles
 lseaice                     =                      .FALSE.       ! .TRUE. for use of sea-ice model
/

! turbdiff_nml: turbulent diffusion -------------------------------------------
&turbdiff_nml
 lconst_z0                   =                         .TRUE.
 const_z0                    =                          0.05      ! roughness length of crops/farmland
 tkhmin                      =                          0.75      ! scaling factor for minimum vertical diffusion coefficient
 tkmmin                      =                          0.75      ! scaling factor for minimum vertical diffusion coefficient
 pat_len                     =                        300.0       ! effective length scale of thermal surface patterns
 tur_len                     =                        500.0       ! asymptotic maximal turbulent distance
 rlam_heat                   =                         10.0
 rat_sea                     =                          0.8       ! controls laminar resistance for sea surface
 frcsmot                     =                          0.2       ! these 2 switches together apply vertical smoothing of the TKE source terms
 imode_frcsmot               =                            2       ! in the tropics (only), which reduces the moist bias in the tropical lower troposphere
 itype_sher                  =                            2       ! type of shear forcing used in turbulence
 ltkeshs                     =                        .TRUE.      ! include correction term for coarse grids in hor. shear production term
 icldm_turb                  =                            1       ! mode of cloud water representation in turbulence
 icldm_tran                  =                            1       ! mode of cloud water representation in turbulence in transfer layer
 itype_wcld                  =                            1       ! type of water cloud diagnosis within the turbulence scheme:1: employing a scheme based on relative humitidy
 q_crit                      =                          1.6       ! critical value for normalized supersaturation
 alpha1                      =                          0.125
 lfreeslip                   =                        .TRUE.      ! use a free-slip lower boundary condition, i.e. neither momentum nor heat/moisture fluxes
/

! nonhydrostatic_nml: nonhydrostatic model -----------------------------------
&nonhydrostatic_nml
 iadv_rhotheta               =                          2         ! advection method for rho and rhotheta
 ivctype                     =                          2         ! type of vertical coordinate
 itime_scheme                =                          4         ! time integration scheme
 ndyn_substeps               =                          5         ! number of dynamics steps per fast-physics step
 exner_expol                 =                          0.333     ! temporal extrapolation of Exner function
 vwind_offctr                =                          0.2       ! off-centering in vertical wind solver
 damp_height                 =                      20000.0       ! height at which Rayleigh damping of vertical wind starts
 rayleigh_coeff              =                          5.0       ! Rayleigh damping coefficient
 divdamp_order               =                         24         ! order of divergence damping 
 divdamp_type                =                         32         ! type of divergence damping
 divdamp_fac                 =                          0.004     ! scaling factor for divergence damping
 igradp_method               =                          3         ! discretization of horizontal pressure gradient
 l_zdiffu_t                  =                      .TRUE.        ! specifies computation of Smagorinsky temperature diffusion
 thslp_zdiffu                =                          0.02      ! slope threshold (temperature diffusion)
 thhgtd_zdiffu               =                        125.0       ! threshold of height difference (temperature diffusion)
 htop_moist_proc             =                      22500.0       ! max. height for moist physics
 hbot_qvsubstep              =                      22500.0       ! height above which QV is advected with substepping scheme
/

! nwp_phy_nml: switches for the physics schemes ------------------------------
&nwp_phy_nml
 inwp_gscp                   =                         -1         ! cloud microphysics and precipitation
 icpl_aero_gscp              =                          0         ! coupling between autoconversion and Tegen aerosol climatology
 inwp_convection             =                          0         ! convection
 icpl_aero_conv              =                          0         ! coupling between autoconversion and Tegen aerosol climatology
 inwp_radiation              =                          4         ! use ecrad
 inwp_cldcover               =                          1         ! cloud cover scheme for radiation
 inwp_turb                   =                          1         ! vertical diffusion and transfer
 inwp_satad                  =                          0         ! saturation adjustment
 inwp_sso                    =                          0         ! subgrid scale orographic drag
 inwp_gwd                    =                          0         ! non-orographic gravity wave drag
 inwp_surface                =                          1         ! surface scheme
 latm_above_top              =                      .TRUE.        ! take into account atmosphere above model top for radiation computation
 lsgs_cond                   =                      .TRUE.
 ldetrain_conv_prec          =                      .TRUE.
 itype_z0                    =                          2         ! type of roughness length data
 icapdcycl                   =                          3         ! apply CAPE modification to improve diurnalcycle over tropical land
 icpl_o3_tp                  =                          1
 icpl_rad_reff               =                          0
 icalc_reff                  =                          0
 dt_conv                     =                          60.       ! time step for convection parameterization
 dt_rad                      =                          60.       ! time step for radiation
 dt_sso                      =                          60.       ! time step for SSO parameterization
/

&p3_nml
 n_icecat                    =                          2         ! number of free ice categories (P3 scheme)
 l3mom_ice                   =                      .TRUE.        ! use triple moment ice categories (P3 scheme)
 lliqfrac                    =                      .TRUE.        ! use predicted bulk liquid fraction (P3 scheme)
 ihydrometeor_ini            =                          1         ! 0: dry initialization, 1: read 1M-scheme tracer, 3: read P3 tracer
 tracer_ini_filename         =       '$INILBCDIR/$ININAME'
 lookup_tables_path          =     '${lookup_tables_path}'
/

! nwp_tuning_nml: additional tuning parameters ----------------------------------
&nwp_tuning_nml
 itune_slopecorr             =                          1
 tune_gkwake                 =                          0.25
 tune_gkdrag                 =                          0.0
 tune_gfrcrit                =                          0.333
 tune_grcrit                 =                          0.25
 tune_minsnowfrac            =                          0.3
 tune_box_liq                =                          0.04
 tune_box_liq_asy            =                          4.0
 tune_gust_factor            =                          7.0
 tune_sgsclifac              =                          1.0
 tune_zvz0i                  =                          0.85
 tune_zceff_min              =                          0.01
 tune_icesedi_exp            =                          0.3
 tune_rcucov                 =                          0.075
 tune_rhebc_land             =                          0.825
 icpl_turb_clc               =                          2
 max_calibfac_clcl           =                          2.0
/

! radiation_nml: radiation scheme ---------------------------------------------
&radiation_nml
 ecrad_data_path             =                'ecrad_data'        ! folder containing ecRad optical properties files
 ecrad_isolver               =                          0
 irad_o3                     =                         79         ! ozone climatology
 irad_aero                   =                          0         ! aerosols
 islope_rad                  =                          0         ! Slope correction for surface radiation
 ecrad_igas_model            =                          0
 ecrad_llw_cloud_scat        =                     .FALSE.
 ecrad_use_general_cloud_optics =                  .FALSE.
 albedo_type                 =                          3         ! fixed surface albedo
 albedo_fixed                =                        0.2
 vmr_co2                     =                    408.e-06
 vmr_ch4                     =                   1850.e-09
 vmr_n2o                     =                   331.0e-09
 vmr_o2                      =                     0.20946
 vmr_cfc11                   =                    240.e-12
 vmr_cfc12                   =                    532.e-12
/

! sleve_nml: vertical level specification -------------------------------------
&sleve_nml
 min_lay_thckn               =                         20.0       ! layer thickness of lowermost layer
 top_height                  =                      30000.0       ! height of model top
 stretch_fac                 =                          0.65      ! stretching factor to vary distribution of model levels
 decay_scale_1               =                       4000.0       ! decay scale of large-scale topography component
 decay_scale_2               =                       2500.0       ! decay scale of small-scale topography component
 decay_exp                   =                          1.2       ! exponent of decay function
 flat_height                 =                      20000.0       ! height above which the coordinate surfaces are flat
/

! transport_nml: tracer transport ---------------------------------------------
&transport_nml
 ivadv_tracer                =                       3, 3         ! tracer specific method to compute vertical advection
 itype_hlimit                =                       3, 4         ! type of limiter for horizontal transport
 ihadv_tracer                =                      52, 2         ! tracer specific method to compute horizontal advection
/

! output = invariant model fields
&output_nml
 filetype                    =                          4         ! output format: 2=GRIB2, 4=NETCDFv2
 dom                         =                          1         ! write all domains
 output_bounds               =               0., 0., 1800.        ! output: start, end, increment
 steps_per_file              =                          1         ! number of output steps in one output file
 steps_per_file_inclfirst    =                     .FALSE. ! first step not counted wrt steps_per_file count (otherwise first file contains 2 times)
 mode                        =                          1         ! 1: relative t-axis, 2: absolute t-axis
 include_last                =                     .FALSE.        ! flag whether to include the last time step
 output_filename             =         'icon_inv_icongrid'
 filename_format             = '<output_filename>_DOM<physdom>' ! file name base
 output_grid                 =                     .FALSE.        ! flag whether grid information is added to output.
 remap                       =                          0         ! 1: remap to lat-lon grid
 ml_varlist                  = 'z_ifc','z_mc','topography_c','fr_land'!,'soiltyp','gz0'
/

! output = single & model level output
&output_nml
 filetype                    =                          4         ! output format: 2=GRIB2, 4=NETCDFv2
 dom                         =                          1         ! write domain 1 only
 output_bounds               =              0., 1800., 60.        ! start, end, increment
 steps_per_file              =                          1         ! number of steps per file
 steps_per_file_inclfirst    =                     .FALSE. ! first step not counted wrt steps_per_file count (otherwise first file contains 2 times)
 mode                        =                          1         ! 1: forecast mode (relative t-axis), 2: climate mode (absolute t-axis)
 include_last                =                     .FALSE.
 output_filename             =         'NWP_LAM_icongrid'
 filename_format             = '<output_filename>_DOM<physdom>_<datetime2>_slml' ! file name base
 output_grid                 =                     .FALSE.
 remap                       =                          0         ! 1: remap to lat-lon grid
 ml_varlist='qv','qc','qr','qnc','qnr','dmean_c','dmean_r','qi','qs','rho','u','v','w','pres','temp','rh',
 'qitot_1','qnitot_1','qzitot_1','qirim_1','birim_1','qiliq_1','dmean_i1',
 'qitot_2','qnitot_2','qzitot_2','qirim_2','birim_2','qiliq_2','dmean_i2'!,
 !'qitot_3','qnitot_3','qzitot_3','qirim_3','birim_3','qiliq_3','dmean_i3'
/

EOF


# ----------------------------------------------------------------------
# run the model!
# ----------------------------------------------------------------------

cp $MODEL ./icon
touch icon
chmod +x icon

cat > job_ICON << EOF
#!/bin/bash -l
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=76
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --time=00:02:00


module purge
#module add compiler/intel/2023.1.0_llvm mpi/openmpi/5.0 lib/netcdf-fortran/4.6-intel-2023.1.0_llvm lib/hdf5/1.14-serial lib/netcdf/4.9-serial eccodes-2.41.0
#module add compiler/intel/2024.0_llvm mpi/openmpi/4.1 lib/netcdf/4.9_serial lib/hdf5/1.12_serial lib/netcdf-fortran/4.6_serial lib/eccodes/2.32.0 numlib/mkl/2022.0.2 libfyaml
#module add compiler/intel/2025.1_llvm mpi/openmpi/4.1 lib/netcdf/4.9-serial lib/hdf5/1.14-serial lib/netcdf-fortran/4.6-intel-2025.1_llvm lib/eccodes_nonetcdf/2.32.0 numlib/mkl/2025.1 libfyaml
module add compiler/intel/2025.1_llvm mpi/openmpi/5.0 lib/netcdf/4.9-serial lib/hdf5/1.14-serial lib/netcdf-fortran/4.6-intel-2025.1_llvm lib/eccodes_nonetcdf/2.32.0 numlib/mkl/2025.1 libfyaml

export OMPI_MCA_coll="^hcoll" # remove mpi-related warnings
export UCX_LOG_LEVEL=error    # remove UCX-related warnings
export UCX_TLS=dc,self        # remove job-failed status
#export OMP_NUM_THREADS=1
mpirun ${EXPDIR}/icon

EOF

#sbatch -A hk-project-aci -p cpuonly --job-name=icon-lam --mail-type=BEGIN,END,FAIL --mail-user=marco.wurth@kit.edu job_ICON
sbatch -A hk-project-aci -p dev_cpuonly --job-name=icon-lamtest-p3 job_ICON
