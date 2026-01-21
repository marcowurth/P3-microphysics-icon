
# ----------------------------------------------------------------------
# SETTINGS: Input/output directories
# ----------------------------------------------------------------------

ICONDIR="/home/hk-project-aci/nw5893/ICON/icon-release-2025.10"

# domain directory
domain_name=torus_40km_100m
STARTDATE="2024-05-08T00:00:00Z"
ENDDATE="2024-05-08T00:30:00Z"
dtime="1.0"

# output directory
EXPDIR=/hkfs/work/workspace/scratch/nw5893-wsx/wk82/${domain_name}/output_P3_test

# absolute path to radiation input directory
RADDIR=${ICONDIR}/externals/ecrad

# paths to external directory
GRID_EXTPARDIR_ORIG="/home/hk-project-aci/nw5893/grid_extpar/${domain_name}"

# filenames
grid_filename="Torus_Triangles_40km_x_40km_res100m.nc"

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
rm -r *.nc

# torus grid
cp ${GRID_EXTPARDIR_ORIG}/${grid_filename} ${grid_filename}

# files needed for radiation
cp -r ${RADDIR}/data ecrad_data

# Dictionary for the mapping: DWD GRIB2 names <-> ICON internal names
#ln -sf ${ICONDIR}/run/ana_varnames_map_file_ICON-EU.txt map_file.ana

# Dictionary for the mapping: GRIB2/Netcdf input names <-> ICON internal names
#ln -sf ${ICONDIR}/run/dict.latbc_ICON-EU map_file.latbc

# configuration
#cp $HOME/.bashrc .
cp ${ICONDIR}/config.log .

mkdir -p p3plugin
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
 num_io_procs                =                          2         ! number of I/O processors
 num_restart_procs           =                          0         ! number of restart processors
 num_prefetch_proc           =                          0         ! number of processors for LBC prefetching
 iorder_sendrecv             =                          3         ! sequence of MPI send/receive calls
 use_omp_input               =                     .FALSE.        ! allows task parallelism for reading atmospheric input data
/

! run_nml: general switches ---------------------------------------------------
&run_nml
 ltestcase                   =                      .TRUE.        ! real case run
 num_lev                     =                        150         ! number of full levels (atm.) for each domain
 lvert_nest                  =                     .FALSE.        ! no vertical nesting
 dtime                       =                     $dtime         ! timestep in seconds
 ldynamics                   =                      .TRUE.        ! compute adiabatic dynamic tendencies
 ltransport                  =                      .TRUE.        ! compute large-scale tracer transport
 iforcing                    =                          3         ! forcing of dynamics and transport by parameterized processes
 msg_level                   =                         10         ! detailed report during integration
 ltimer                      =                      .TRUE.        ! timer for monitoring the runtime of specific routines
 timers_level                =                         10         ! performance timer granularity
 check_uuid_gracefully       =                      .TRUE.        ! give only warnings for non-matching uuids
 output                      =                       "nml"        ! main switch for enabling/disabling components of the model output
/

&nh_testcase_nml
 nh_test_name  = 'wk82'    ! test case identifier
 qv_max_wk     = 0.014
 u_infty_wk    = 0.0
 bub_amp       = 3.
 bub_hor_width = 10000.
 bub_ver_width = 1500.
 bubctr_lon    = 20000.0
 bubctr_lat    = 20000.0
 bubctr_z      = 1500.
/

! diffusion_nml: horizontal (numerical) diffusion ----------------------------
&diffusion_nml
 lsmag_3d                    =                     .FALSE.
 lhdiff_vn                   =                      .TRUE.        ! diffusion on the horizontal wind field
 lhdiff_temp                 =                      .TRUE.        ! diffusion on the temperature field
 lhdiff_w                    =                      .TRUE.        ! diffusion on the vertical wind field
 lhdiff_q                    =                      .TRUE.        ! diffusion on the vertical wind field
 hdiff_order                 =                          5         ! order of nabla operator for diffusion
 itype_vn_diffu              =                          1         ! reconstruction method used for Smagorinsky diffusion
 itype_t_diffu               =                          2         ! discretization of temperature diffusion
/

! dynamics_nml: dynamical core -----------------------------------------------
&dynamics_nml
 lcoriolis                   =                     .FALSE.        ! Coriolis force
/

! extpar_nml: external data --------------------------------------------------
&extpar_nml
 itopo                       =                          0         ! topography (0:analytical)
 n_iter_smooth_topo          =                          1         ! iterations of topography smoother
 heightdiff_threshold        =                        400.        ! height difference between neighb. grid points
 hgtdiff_max_smooth_topo     =                        750.        ! see Namelist doc
 itype_vegetation_cycle      =                          2         ! 2=operational, 1=default
 read_nc_via_cdi             =                      .TRUE.
/

! grid_nml: horizontal grid --------------------------------------------------
&grid_nml
 dynamics_grid_filename      =           '${grid_filename}'       ! array of the grid filenames for the dycore
 lredgrid_phys               =                      .FALSE.       ! .true.=radiation is calculated on a reduced grid
 l_limited_area              =                      .FALSE.       ! .TRUE. performs limited area run
 is_plane_torus              =                      .TRUE.        ! feedback type (incremental/relaxation-based)
/

! ls_forcing_nml: large scale forcing --------------------------------------
!&ls_forcing_nml
! is_ls_forcing               =                     .FALSE.
!/

! gridref_nml: grid refinement settings --------------------------------------
&gridref_nml
 denom_diffu_v               =                        150.        ! denominator for lateral boundary diffusion of velocity
/

! io_nml: general switches for model I/O -------------------------------------
&io_nml
 itype_pres_msl              =                          5         ! method for computation of mean sea level pressure
 itype_rh                    =                          1         ! method for computation of relative humidity
 lmask_boundary              =                     .FALSE.        ! mask out interpolation zone in output
 celltracks_interval         =                        60.
 dt_celltracks               =                        30.
 wshear_uv_heights           =                      6000.
/

! nonhydrostatic_nml: nonhydrostatic model -----------------------------------
&nonhydrostatic_nml
 iadv_rhotheta               =                          2         ! advection method for rho and rhotheta
 ivctype                     =                          2         ! type of vertical coordinate
 itime_scheme                =                          4         ! time integration scheme
 ndyn_substeps               =                          5         ! number of dynamics steps per fast-physics step
 exner_expol                 =                          0.333     ! temporal extrapolation of Exner function
 vwind_offctr                =                          0.2       ! off-centering in vertical wind solver
 damp_height                 =                      15000.0       ! height at which Rayleigh damping of vertical wind starts
 rayleigh_coeff              =                          5.0       ! Rayleigh damping coefficient
 divdamp_order               =                         24         ! order of divergence damping 
 divdamp_type                =                         32         ! type of divergence damping
 divdamp_fac                 =                          0.004     ! scaling factor for divergence damping
 igradp_method               =                          3         ! discretization of horizontal pressure gradient
 l_zdiffu_t                  =                      .TRUE.        ! specifies computation of Smagorinsky temperature diffusion
 thslp_zdiffu                =                          0.02      ! slope threshold (temperature diffusion)
 thhgtd_zdiffu               =                        125.0       ! threshold of height difference (temperature diffusion)
 htop_moist_proc             =                      20000.0       ! max. height for moist physics
 hbot_qvsubstep              =                      20000.0       ! height above which QV is advected with substepping scheme
/

! nwp_phy_nml: switches for the physics schemes ------------------------------
&nwp_phy_nml
 inwp_gscp                   =                         -1         ! cloud microphysics and precipitation
 inwp_convection             =                          0         ! convection
 inwp_radiation              =                          0         ! use ecrad
 inwp_cldcover               =                          0         ! cloud cover scheme for radiation
 inwp_turb                   =                          5         ! vertical diffusion and transfer
 inwp_satad                  =                          0         ! saturation adjustment
 inwp_sso                    =                          0         ! subgrid scale orographic drag
 inwp_gwd                    =                          0         ! non-orographic gravity wave drag
 inwp_surface                =                          0         ! surface scheme
 lsgs_cond                   =                     .FALSE.
 icpl_rad_reff               =                          0
 icalc_reff                  =                          0
 dt_conv                     =                          60.       ! time step for convection parameterization
 dt_ccov                     =                          60.       ! time step for convection parameterization
 dt_rad                      =                          60.       ! time step for radiation
 dt_sso                      =                          60.       ! time step for SSO parameterization
/

&p3_nml
 n_icecat                    =                          2          ! number of free ice categories (P3 scheme)
 l3mom_ice                   =                      .TRUE.         ! use triple moment ice categories (P3 scheme)
 lliqfrac                    =                      .TRUE.         ! use predicted bulk liquid fraction (P3 scheme)
 itracer_ini                 =                          0          ! 0: dry initialization, 1: read 1M-scheme tracer, 3: read P3 tracer
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

! sleve_nml: vertical level specification -------------------------------------
&sleve_nml
 min_lay_thckn               =                         20.0       ! layer thickness of lowermost layer
 top_height                  =                      20000.0       ! height of model top
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

! les_nml: 3d turbulent diffusion -------------------------------------------
&les_nml
 smag_constant               =                            0.23
 turb_prandtl                =                            0.333333
 max_turb_scale              =                          300.
 isrfc_type                  =                            0       ! 0=no sfc fluxes, 1=TERRA, 2=Fixed flux
 vert_scheme_type            =                            2
 ldiag_les_out               =                       .FALSE.
 les_metric                  =                       .FALSE.
/

! turbdiff_nml: turbulent diffusion -------------------------------------------
&turbdiff_nml
 lconst_z0                   =                        .TRUE.
 const_z0                    =                            0.05       ! roughness length of crops/farmland
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
 ml_varlist                  = 'z_ifc','z_mc','topography_c','fr_land','soiltyp','gz0'
/

! output = single & model level output
&output_nml
 filetype                    =                          4         ! output format: 2=GRIB2, 4=NETCDFv2
 dom                         =                          1         ! write domain 1 only
 output_bounds               =             0., 93600., 10.        ! start, end, increment
 steps_per_file              =                          1         ! number of steps per file
 steps_per_file_inclfirst    =                     .FALSE. ! first step not counted wrt steps_per_file count (otherwise first file contains 2 times)
 mode                        =                          1         ! 1: forecast mode (relative t-axis), 2: climate mode (absolute t-axis)
 include_last                =                     .FALSE.
 output_filename             =         'NWP_LAM_icongrid'
 filename_format             = '<output_filename>_DOM<physdom>_<datetime2>_slml' ! file name base
 output_grid                 =                     .FALSE.
 remap                       =                          0         ! 1: remap to lat-lon grid
 !ml_varlist='qv','qc','qr','qnc','qnr','qi','rho','u','v','w','pres','temp','rh',
 !'clcl','clcm','clch','clct','dmean_c','dmean_r','deff_c','dhmax','ze_p3',
 !'qitot_1','qnitot_1','qzitot_1','qirim_1','birim_1','qiliq_1','dmean_i1','deff_i1','rho_i1','vm_i1',
 !'qitot_2','qnitot_2','qzitot_2','qirim_2','birim_2','qiliq_2','dmean_i2','deff_i2','rho_i2','vm_i2'
 !'qitot_3','qnitot_3','qzitot_3','qirim_3','birim_3','qiliq_3','dmean_i3','deff_i3','rho_i3','vm_i3'
 ml_varlist='qv','qc','qr','qnc','qnr','dmean_c','dmean_r','qi','qs','rho','u','v','w','pres','temp','rh'!,
 'qitot_1','qnitot_1','qzitot_1','qirim_1','birim_1','qiliq_1','dmean_i1',
 'qitot_2','qnitot_2','qzitot_2','qirim_2','birim_2','qiliq_2','dmean_i2'!,
/

! output = single & model level output
!&output_nml
! filetype                    =                          4         ! output format: 2=GRIB2, 4=NETCDFv2
! dom                         =                          1         ! write domain 1 only
! output_bounds               =             0., 93600., 10.        ! start, end, increment
! steps_per_file              =                          1         ! number of steps per file
! steps_per_file_inclfirst    =                     .FALSE. ! first step not counted wrt steps_per_file count (otherwise first file contains 2 times)
! mode                        =                          1         ! 1: forecast mode (relative t-axis), 2: climate mode (absolute t-axis)
! include_last                =                     .FALSE.
! output_filename             =           'NWP_LAM_latlon'
! filename_format             = '<output_filename>_DOM<physdom>_<datetime2>_slml' ! file name base
! output_grid                 =                     .FALSE.
! remap                       =                          1  ! 1: remap to lat-lon grid
! reg_lon_def                 =            -180.0,2.0,180.0
! reg_lat_def                 =             -90.0,1.0,90.0
! ml_varlist='z_mc','qc','qr','qnc','qnr','qi','rho','w','pres','temp','ddt_temp_phys',
! 'dmean_c','dmean_r','deff_c','dhmax','ze_p3',
! 'qitot_1','qnitot_1','qirim_1','birim_1','dmean_i1','deff_i1','rho_i1','vm_i1',
! 'qitot_2','qnitot_2','qirim_2','birim_2','dmean_i2','deff_i2','rho_i2','vm_i2'
! !'qitot_3','qnitot_3','qirim_3','birim_3','dmean_i3','deff_i3','rho_i3','vm_i3'
!/

! output = 1min precipitation output
!&output_nml
! filetype                    =                          4         ! output format: 2=GRIB2, 4=NETCDFv2
! dom                         =                          1         ! write domain 1 only
! output_bounds               =             0., 93600., 60.        ! start, end, increment
! steps_per_file              =                          1         ! number of steps per file
! steps_per_file_inclfirst    =                     .FALSE. ! first step not counted wrt steps_per_file count (otherwise first file contains 2 times)
! mode                        =                          1         ! 1: forecast mode (relative t-axis), 2: climate mode (absolute t-axis)
! include_last                =                     .FALSE.
! output_filename             =         'NWP_LAM_icongrid'
! filename_format             = '<output_filename>_DOM<physdom>_<datetime2>_prec' ! file name base
! output_grid                 =                     .FALSE.
! remap                       =                          0         ! 1: remap to lat-lon grid
! ml_varlist='tot_prec','rain_gsp_rate','snow_gsp_rate','dhmax_ground'
!/
!
!! output = 1min-max convective variables output
!&output_nml
! filetype                    =                          4         ! output format: 2=GRIB2, 4=NETCDFv2
! dom                         =                          1         ! write domain 1 only
! output_bounds               =             0., 93600., 60.        ! start, end, increment
! steps_per_file              =                          1         ! number of steps per file
! steps_per_file_inclfirst    =                     .FALSE. ! first step not counted wrt steps_per_file count (otherwise first file contains 2 times)
! mode                        =                          1         ! 1: forecast mode (relative t-axis), 2: climate mode (absolute t-axis)
! include_last                =                     .FALSE.
! output_filename             =         'NWP_LAM_icongrid'
! filename_format             = '<output_filename>_DOM<physdom>_<datetime2>_conv_max' ! file name base
! output_grid                 =                     .FALSE.
! remap                       =                          0         ! 1: remap to lat-lon grid
! ml_varlist='uh_max','uh_max_med','uh_max_low','vorw_ctmax','w_ctmax'
!/
EOF

# ----------------------------------------------------------------------
# run the model!
# ----------------------------------------------------------------------

cp $MODEL ./icon
touch icon
chmod +x icon

cat > job_ICON << EOF
#!/bin/bash -l
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=76
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --time=00:15:00
######SBATCH --mem=200gb


module purge
module add compiler/intel/2025.1_llvm mpi/openmpi/5.0 lib/hdf5/1.14 lib/netcdf/4.9 lib/netcdf-fortran/4.6 lib/eccodes_nonetcdf/2.32.0 numlib/mkl/2025.1 libfyaml

export OMPI_MCA_coll="^hcoll" # remove mpi-related warnings
export UCX_LOG_LEVEL=error    # remove UCX-related warnings
export UCX_TLS=dc,self        # remove job-failed status
#export OMP_NUM_THREADS=1
mpirun ${EXPDIR}/icon

EOF

#sbatch -A hk-project-aci -p dev_cpuonly --job-name=icon-lam --mail-type=BEGIN,END,FAIL --mail-user=marco.wurth@partner.kit.edu job_ICON
sbatch -A hk-project-aci -p dev_cpuonly --job-name=icon-torus-les-wk82-warmbubble-p3 job_ICON
