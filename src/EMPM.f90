                PROGRAM EMPM

! -------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------
!
! EMPM - Explicit Mixing Parcel Model:
! ------------------------------------
!
! 'The EMPM [...] predicts the evolving in-cloud variability of temperature and water vapor mixing
! ratio due to entrainment and finite-rate turbulent mixing using 1D representation of a rising 
! cloudy parcel.' (Krueger, 2008)
!
!
! Authors:     Steven Krueger, Phil Austin, Sonja Drueke
!
! Date:        SEPTEMBER 2017
!
!
! Updates:
! --------
! Sep 28th, 2017
! -- Include temptime and qvtime array in netCDF output file when netCDF flag is chosen
! -- Include spline interpolation for small droplet radii (contains still a bug)
!
! Sep 21st, 2017
! -- Fix all compiler (ifort) warnings and errors (such as undeclared, undefind or unused variables)
!
! Sep 19th, 2017
! -- Allow netCDF output for rtime and xtime variables (possible extension to other variables?)
!
! Aug 23rd, 2017
! -- Allow unformatted FORTRAN file format for rtime.dat and xtime.dat
! 
! July 26th, 2017
! -- Increase array size for writing out individual droplet data (m_index)
!
! July 25th, 2017
! -- Include threshold value (read from namelist) for vertical velocity at which EMPM is stopped
!
! July 24th, 2017
! -- Modify index.f90 to accommodate to new length of incloud vertical velocity profile
!
! July 18th, 2017
! -- Length of incloud vertical velocity is now independent of the environment profile length
! -- When end of incloud vertical velocity is reached, EMPM is stopped.
! -- Correct the units of the molar weight of water (from g/mol to kg/mol)
!
! August 12th, 2016
! -- implementing calculation of LWP
! -- including calculation for effective radius
!
! August 2nd, 2016
! -- including cloud top height as one limit to end the simulation
! -- fixing bug that led to non-constant total water mixing ratio in the case of a gamma-shaped 
!    initial droplet size distribution
!
! June 9th, 2016
! -- fixing bug related to vertical velocity
!    - remove pressure calculation in intpol.f90 and use the in the droplet growth model calculated
!      pressure instead
! -- change buoyancy calculation from actual temperature to virtual potential temperature)

! June 2nd, 2016
! -- including buoyancy calculation and calculate vertical velocity based on buoyancy
!
! May 27th, 2016
! -- implementing option to read incloud vertical velocity from file
!
! May 18th, 2016
! -- fixing a bug in size_dis.f90 which led to an overestimation of liquid water mixing ratio
!
! May 5th, 2016
! -- included the calculation of molecular, temperature and turbulent diffusivity in dependence of
!    dissipation rate of turbulent kinetic energy which has to be set in the namelist
!
! May 4th, 2016
! -- changed all REAL variables to REAL*8
!
! April 22nd, 2016
! -- bug fix: entrained ccn read from file works properly now
!
! April 21st, 2016
! -- fixing the bug that led to a too fast droplet growth 
!
! March 16th, 2016
! -- including the possibility to hold the droplet number constant
! -- entrained ccn can now have either
!      - all the same radius and mass
!      - or the size distribution of ccn can be read from file
!
! March 4th, 2016
! -- implement initial droplet size distribution calculation using gamma function
!
! January 25th, 2016
! -- the old FORTRAN77 code was rewritten in FORTRAN90
!
! March 13th, 1995
! -- revise to let ccn contained in entrained clear air 
!
! -------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------
!
! -- PORCEDURES FOR RUNNING THIS PROGRAM:
!
! 1. Make appropriate changes in Makefile (choose hardware) 
! 2. compile EMPM code by typing make (in src directory)
! 3. Set up your simulation: 
!    - make changes in namelist according to your demands
!    - initial profiles of temperature and humidity and the initial droplet information have to be
!      given in two files (please refer to the example for information about the format)
!    - (optional: ccn data of the entrained air can be read from a third file)
!    - (optional: incloud vertical velocity can be provided in a forth file)
! 4. execute EMPM by typing ./EMPM (in src directory)
!
! -------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------

        USE const
        USE intspline
        USE netcdf

        IMPLICIT NONE

! -------------------------------------------------------------------------------------------------
! domain and time step variables
! -------------------------------------------------------------------------------------------------
! -- grid points in the domain (ngrid)
! -- grid size (dx)
        INTEGER   :: ngrid
        REAL*8    :: dx
! -- m1 - first point in the flow domain
! -- m2 - first point of each entrained blob
! -- m3 - last point of each entrained blob
! -- m4 - last point in the domain (at the right wall boundary)
! -- m6 - last point at one domain distance right to boundary wall
! -- grid point in each entrained blob (xn)
        INTEGER   :: m1
        INTEGER, DIMENSION(:), ALLOCATABLE :: m2, m3
        INTEGER   :: m4, m6
        INTEGER   :: xn
! -- dimensionless domain length (bl), largest eddy size (bl_inte)
! -- smallest eddy size, essentially the model's Komolgorov length scale (eta)
! -- domain volume (volume)
        REAL*8    :: bl, bl_inte, eta 
        REAL*8    :: volume 

! -- flag whether diffusion time step shall be given in namelist (itd)
! -- time step both for turbulent convection (tc) and diffusion (td)
! -- auxiliary variables to calculate the time step for turbulent convection tc (tau_l and tlambda)
! -- time (t), time step determined by comparing tc and td and used in calculation (dt)
! -- start and end time (given in namelist - start_t and end_t, respectively)
        INTEGER   :: itd
        REAL*8    :: tc, td
        REAL*8    :: tau_l, tlambda
        REAL*8    :: t, dt
        REAL*8    :: start_t, end_t
        REAL*8    :: time_step

! -- random length scale for triple mapping (dl)
! -- random position for triplet mapping (p2)
! -- first point of the triplet mapping eddy (i1), grid points inside the random eddy size (n)
! -- auxiliary variable to calculate the grip points (ni)
! -- last point of the triplet mapping eddy (nn)
! -- after how many diffusion 'events' a triplet mapping 'event' (turbulence) is carried out 
! -- (idtrip), how many triplet mapping 'events' (turbulence) per diffusion 'event' (mtrip)
        REAL*8    :: dl
        REAL*8    :: p2
        INTEGER   :: i1, n, ni, nn
        INTEGER   :: idtrip, mtrip

! -- initial time or 'old' time (time0)
! -- time step for entrainment events (dt_entm), entrainment rate (ent_rate)
! -- accumulative time for all entrainment events (t_entm)
! -- size of the entrained parcel in relation to domain size (fraction - psigma)
! -- volume of the entrained air (volume_entm)
        REAL*8    :: time0
        REAL*8    :: dt_entm, ent_rate, t_entm
	REAL*8    :: psigma
        REAL*8    :: volume_entm

! -- number of entrained blobs (n_blob)
! -- final number of blobs (n_blob_final) - it might be n_blob+1 due to periodic boundary condition
! -- maximal number of blobs possible (n_max)
! -- start and end point of entrained parcel (spoint and fpoint, respectively)
! -- variables to encounter split blobs due to periodic B.C. (part and blob_flag)
        INTEGER   :: n_blob, n_blob_final, n_max
        INTEGER, DIMENSION(:), ALLOCATABLE :: spoint, fpoint
        INTEGER   :: part, blob_flag

! -- length of the history arrays (hist_array)
! -- interval (t_interval) for writing out history data
! -- array size for writing out data (has to agree with end_t / t_interval - set in namelist - mip)
        INTEGER   :: hist_array
        REAL*8    :: t_interval
        INTEGER   :: mip

! -- running index (iteration)
! -- start and end iteration index (iteration_s and iteration_e, respectively)
! -- index for calculating statistics (iteration1), maximal number of iteration (iteration2)
! -- index for writing out history data for t_interval (iteration3)
        INTEGER   :: iteration, iteration_s, iteration_e
        INTEGER   :: iteration1, iteration2, iteration3

! -- current number of iteration periods in one realization (ip),
! -- maximal number of iterations in one realization (ipreal)
        INTEGER   :: ip, ipreal

! -- index for realizations (mreal), maximal number of realizations (mrealization)
        INTEGER   :: mreal, mrealization

! -- combining array for dm and dtemp (dA)
! -- molecular (dm), temperature (dtemp) and turbulent (dtur) diffusivity
! -- dissipation rate of turbulent kinetic energy
	REAL*8, DIMENSION(2) :: dA
        REAL*8    :: dm=0, dtemp=0, dtur=0
        REAL*8    :: epsilon
! -- Reynolds number (Re)
        REAL*8    :: Re

! -- array size for droplet properties (idimen, ne, ne1)
        INTEGER   :: idimen, ne, ne1

! -- flag to set the I.C.s and B.C.s in the drop_map subroutine (flag)
        INTEGER   :: flag

! -- flag for switching on/off isobaric mixing (iim)
! -- flag for switching on/off instant mixing (iinm)
! -- flag for switching on/off random variation of calculated timestep for entrainment events (ire)
! -- flag for switching on/off implemtation of ccn in entrained air parcel (iccn)
! -- flag for switching on/off constant droplet number in domain (idc)
! -- flag for switching on/off constant vertical velocity (icw)
! -- flag for switching on/off incloud vertical velocity profile (ivp)
! -- flag for switching on/off cloud base entrainment (cb_ent) - namelist
! -- flag for switching on/off stopping EMPM when a vertical velocity smaller than threshold is 
! -- reached (icaw)
! -- flag for entrainment at cloud base (first time step) (ent_flag) - code

        INTEGER   :: iim, iinm
        INTEGER   :: ire
        INTEGER   :: iccn, idc
        INTEGER   :: icw, ivp
        INTEGER   :: cb_ent
        INTEGER   :: icaw
        LOGICAL   :: ent_flag

! -- running indices
        INTEGER   :: i, j, k

!
! random number generation
! -------------------------------------------------------------------------------------------------
! -- REAL FUNCTION for random number generation (RAND1)
! -- random seeds (iseed1, iseed2, ...)
! -- calculated random numbers (r1, r2)
        REAL*8    :: RAND1
        INTEGER*4 :: iseed1, iseed2, iseed3
        INTEGER*4 :: iseed4, iseed5, iseed6
        REAL*8    :: r1, r2

!
! initial environmental (parcel) profile
! -------------------------------------------------------------------------------------------------
! -- profile length (mi)
! -- profile variables (pressure, mixing ratio, temperature, height)
! -- interpolated mixing ratio (qv_e), temperature (temp_e), height (Z)
! -- potential temperature (pt_e), environmental virtual potential temperature (pt_v_e)
        INTEGER   :: mi
        REAL*8, DIMENSION(:), ALLOCATABLE :: p_array, qv_array, temp_array, z_array
        REAL*8    :: qv_e, temp_e, Z
        REAL*8    :: pt_e, pt_v_e

! -- pressure level (press) and height (height) at each time step
! -- initial pressure and height level (press0 and height0, respectively)
! -- temperature (temp) at each time step
! -- initial temperature (temp0), vertical velocity (w), initial vertical velocity (w0)
! -- liquid water mixing ratio (ql), water vapor mixing ratio (qv)
! -- initial water vapor mixing ratio (qv0), initial liquid water mixing ratio (ql0)
! -- initial volume of liquid water (ql0_vol), initial total water mixing ratio (qt0)
        REAL*8    :: press, height, press0, height0
        REAL*8    :: temp, temp0, w, w0
        REAL*8    :: ql, qv
	REAL*8    :: qv0, ql0, ql0_vol, qt0

! -- pressure and height after the call of the DGM (press_dgm and height_dgm, respectively)
	REAL*8    :: press_dgm, height_dgm

!
! (initial) cloud properties
! -------------------------------------------------------------------------------------------------
! -- flag for initial droplet properties (idp)
        INTEGER   :: idp
! -- initial mass of the solute in the droplet (m_i), initial droplet number concentration (N_i)
! -- initial droplet radius (r_i)
        REAL*8    :: m_i, N_i, r_i
! -- initial liquid water content in cloud (LWC) in kg/m^3, liquid water path (LWP) in kg/m^2
        REAL*8    :: LWC, LWP
! -- initial air density (rho_a) in kg/m^3, domain averaged air density (rho_avg) in kg/m^3
        REAL*8    :: rho_a, rho_avg
! -- vertical velocity inside the cloud (w_cloud), pressure level of the velocity information 
! -- (p_cloud), height level of the velocity information (h_cloud)
        REAL*8, DIMENSION(:), ALLOCATABLE :: w_cloud, p_cloud, h_cloud
! -- profile length for incloud vertical velocity profile
        INTEGER   :: nline_v
        INTEGER   :: count_lines
! -- interpolated vertical velocity (vel) in the case the vertical velocity is read from file
        REAL*8    :: vel
! -- threshold for vertical velocity at which EMPM is stopped (thsw)
        REAL*8    :: thsw

! -- cloud top height (cloudtopheight - set in namelist) to end simulation
        REAL*8 :: cloudtopheight

! -- flag necessay to calculate LWP of entire cloud (LWP_flag)
! -- count of time steps necessary to calculate LWP (count_LWP_
        LOGICAL :: LWP_flag
        REAL*8  :: count_LWP

! -- total water mixing ratio (qt)
! -- moist static energy (hm), liquid water static energy (sl),
! -- domain average of buoyancy (buoy_ave), domain average of potential temperature (pt_ave)
! -- and virtual potential temperature (pt_v_ave), domain average liquid water content (ql_ave)
	REAL*8, DIMENSION(:), ALLOCATABLE  :: qt
	REAL*8, DIMENSION(:), ALLOCATABLE  :: hm, sl
        REAL*8    :: buoy_ave, pt_ave, pt_v_ave, ql_ave

! -- domain average of water vapor mixing ratio (qv_ave) and of temperature (temp_ave)
! -- sum of water vapor mixing ratio (qv_sum) and of temperature (temp_sum)
	REAL*8    :: qv_ave, qv_sum, temp_ave, temp_sum

! -- partial pressure of water vapor (es), mixing ratio (qs)
! -- temperature change (increase) due to rising air (temp_in)
! -- REAL FUNCTION to calculate the saturated partial pressure of water vapor (EW)
	REAL*8    :: es, qs
	REAL*8    :: temp_in
	REAL*8    :: EW

! -- number of scalars (temperature and water vapor mixing ratio - n_SA)
! -- combining array for scalars - mixing ratio (qv/qv_e) and temperature (temp/temp_e) for
! -- whole domain (A), temporay array for execution of diffusion by Euler approximation (B)
! -- seperate arrays for mixing ratio (SA1) and temperature (SA2)
! -- probability density function for SA1 and SA2 (SA1_PDF and SA2_PDF, respectively)
! -- averaged PDF for SA1 and SA2 (SA1_PDF_ave and SA2_PDF_ave, respectively)
! -- maximal and minimal values for SA1 and SA2 for PDF calculation (qv_pdf_max, temp_pdf_max, qv_pdf_min, 
! -- temp_pdf_min, respectively)
        INTEGER   :: n_SA = 2
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: A, B
        REAL*8, DIMENSION(:), ALLOCATABLE   :: SA1, SA2 
        REAL*8, DIMENSION(:), ALLOCATABLE   :: SA1_scalar, SA2_scalar
      	REAL*8, DIMENSION(:,:), ALLOCATABLE :: SA1_PDF, SA2_PDF, SA1_PDF_ave, SA2_PDF_ave
        REAL*8    :: qv_pdf_max, temp_pdf_max, qv_pdf_min, temp_pdf_min

! -- bin number for PDF (scalar_pdf_nbin), array size for scalars (set in namelist - igrid)
        INTEGER   :: scalar_pdf_nbin
        INTEGER*4 :: igrid

! -- number of scalars (n_scalars)
! -- combining array for scalars (qlwater, qt, sl, buoy_ave, hm, SuS, qv and temp - scalar)
! -- domain average of scalars (qlwater, qt, sl, buoy_ave, hm, SuS, qv and temp - scalar_mean)
! -- average  of means of scalars (qlwater, qt, sl, buoy_ave, hm, SuS, qv and temp - 
! -- scalar_mean_ave) of realizations
! -- domain rms of scalars (qlwater, qt, sl, buoy_ave, hm, SuS, qv and temp - scalar_rms)
! -- average  of rms of scalars (qlwater, qt, sl, buoy_ave, hm, SuS, qv and temp - 
! -- scalar_rms_ave) of realizations
! -- probability density function for supersaturation (SuS_pdf)
! -- average pdf of supersaturation of realizations (SuS_pdf_ave)
! -- 'bin values' for pdf of supersaturation (SuS_scalar), dummy variable necessary to call 
! -- subroutine (dummy_SuS)
! -- min and max for pdf calculation of supersaturation (Sus_pdf_max and Sus_pdf_min, respectively)
        INTEGER   :: n_scalar = 8
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: scalar
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: scalar_mean, scalar_mean_ave
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: scalar_rms, scalar_rms_ave
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: SuS_pdf, SuS_pdf_ave
        REAL*8, DIMENSION(:), ALLOCATABLE   :: SuS_scalar, dummy_SuS
        REAL*8    :: Sus_pdf_max, Sus_pdf_min

!
! droplet data
! -------------------------------------------------------------------------------------------------
! -- droplet distribution (dis)
! -- droplet radius (radius_drop), droplet mass (mass_solute) and droplet concentration (con_drop)
! -- initial positions of droplets (x), total water per droplet (qwater)
! -- liquid water mixing ratio per grid cell (qlwater), supersaturation (SuS)
! -- droplet radius (dsize), droplet mass (csize)
! -- terminal velocity of each droplet (v)
! -- auxiliary variable to calculate terminal velocity of droplets (c), droplet radius (radius)
! -- inverse volume of each grid cell (i_grid_vol)
	REAL*8, DIMENSION(:), ALLOCATABLE :: dis
        REAL*8, DIMENSION(:), ALLOCATABLE :: new_con, radius_drop, mass_solute, con_drop
	REAL*8, DIMENSION(:), ALLOCATABLE :: x, qwater
	REAL*8, DIMENSION(:), ALLOCATABLE :: qlwater, SuS
	REAL*8, DIMENSION(:), ALLOCATABLE :: csize, dsize
	REAL*8, DIMENSION(:), ALLOCATABLE :: v
        REAL*8    :: c
	REAL*8    :: radius
        REAL*8    :: i_grid_vol

! -- number of droplet in the domain volume (n_drop), index of droplets (n_drop_index)
! -- number of activated droplets (n_act_drop), index array of all droplets (index_drop)
! -- number of droplets without entrained ccns (n_drop_used)
! -- total number of droplet categories (num_tot)
        INTEGER   :: n_drop, n_drop_index, n_act_drop
        INTEGER, DIMENSION(:), ALLOCATABLE :: index_drop
        INTEGER   :: n_drop_used
        INTEGER   :: num_tot

! -- grid cell index of droplet (j_cell), index for droplet position (jcell_pop)
        INTEGER   :: j_cell
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: jcell_pop

! -- index for triplet mapping (index)
! -- change of droplet position due to triplet mapping (dx_map)
        INTEGER, DIMENSION(:), ALLOCATABLE   :: index
	REAL*8, DIMENSION(:), ALLOCATABLE    :: dx_map

! -- time series of droplet radius (r_time), time series of droplet position (x_time)
! -- control array to check that radius / position are 0 if there are no droplet (time_null)
! -- array size for writing out individual droplet data (m_index)
! -- time at which time series data of droplet radius and position are written out (t_time)
! -- index of time step to write out droplet radius' and positions (t_ip)
	REAL*8, DIMENSION(:), ALLOCATABLE   :: r_time
        REAL*8, DIMENSION(:), ALLOCATABLE   :: x_time
	REAL*8, DIMENSION(:), ALLOCATABLE   :: time_null
        INTEGER   :: m_index
!	REAL*8    :: t_time
!        INTEGER   :: t_ip

! -- auxiliary variables necessary to calculate lenght of individual droplet data
! -- height difference between entrainment events (delta_ent), cloud depth (cloud_depth)
! -- number of entrainment events (number_ent), expected droplet number entrained (n_drop_ent)
! -- expected total number indices (tot_num_index) 
        REAL*8    :: delta_ent, cloud_depth 
        REAL*8    :: number_ent, n_drop_ent, tot_num_index

! -- array of mean droplet radius (r_mean), array of standard deviation of droplets (r_std)
! -- array of maximal droplet radius (r_max)
! -- array of index of maximal droplet radius (r_max_index)
! -- array of number of activated droplets (actndrop)
! -- maximal droplet radius (r_max_drop), sum of standard deviation of droplets (r_std_sum)
! -- sum of droplet radius (r_sum), index of maximal droplet radius (index_r_max)
	REAL*8, DIMENSION(:), ALLOCATABLE  :: r_mean, r_std
	REAL*8, DIMENSION(:), ALLOCATABLE  :: r_max
        INTEGER, DIMENSION(:), ALLOCATABLE :: r_max_index, actndrop
        REAL*8    :: r_max_drop, r_std_sum, r_sum
        INTEGER   :: index_r_max

! -- bin number for droplet radius pdf (r_pdf_nbin)
! -- maximum and minimum value for droplet radius pdf (r_pdf_max and r_pdf_min, respectively)
! -- droplet radius pdf (r_pdf), 'bin value' for droplet radius pdf
        INTEGER   :: r_pdf_nbin
        REAL*8    :: r_pdf_max, r_pdf_min
        REAL*8, DIMENSION(:), ALLOCATABLE  :: r_pdf, r_scalar

! -- mean droplet radius (mean_r), square of mean droplet radius (mean_r2)
! -- mean droplet radius cubed (mean_r3)
        REAL*8, DIMENSION(:), ALLOCATABLE  :: mean_r, mean_r2, mean_r3
! -- effective radius (eff_r) im mu m
        REAL*8    :: eff_r

! -- array of mean supersaturation droplet based (s_mean)
! -- array of standard deviation of supersaturation droplet based (s_std)
! -- sum of standard deviation of supersaturation (s_std_sum)
! -- sum of supersaturation droplet based (s_sum)
! -- maximal supersaturation (s_max), minimal supersaturation (s_min)
! -- difference between maximal and minimal supersaturation (s_diff)
	REAL*8, DIMENSION(:), ALLOCATABLE  :: s_mean, s_std
        REAL*8    :: s_std_sum, s_sum
        REAL*8    :: s_max, s_min, s_diff

! -- type of aerosol (sets the molecular mass of the solute (solute composition, aero)
        INTEGER   :: aero 

!
! entrained ccn
! -------------------------------------------------------------------------------------------------
! -- number of entrained ccn (n_ccn), number of ccn bins given in file (nccn)
! -- number of droplets within the parts of the domain in which the entrained air is mixed (ncut)
! -- mass of entrained ccn (m_ccn), radius of entrained ccn (r_ccn)
! -- radius of entrained ccn (radius_ccn), mass of entrained ccn (mass_ccn), concentration of 
! -- entrained ccn (con_ccn)
        INTEGER   :: n_ccn, nccn, ncut
        REAL*8    :: m_ccn, r_ccn
	REAL*8, DIMENSION(:), ALLOCATABLE :: radius_ccn, mass_ccn, con_ccn

! -- changed index for droplet properties due to entrainment (ichange)
        INTEGER   :: ichange

!
! pdf of total water mixing ratio field
! -------------------------------------------------------------------------------------------------
! -- pdf of total water mixing ratio (qt_pdf), 'bin values' for pdf of qt (qt_scalar)
! -- min and max for pdf calculation of qt (qt_pdf_max and qt_pdf_min, respectively)
! -- bin number for pdf (qt_pdf_nbin)
	REAL*8, DIMENSION(:), ALLOCATABLE   :: qt_pdf, qt_scalar
	REAL*8    :: qt_pdf_max, qt_pdf_min
        INTEGER   :: qt_pdf_nbin

! -- variables for calculating simulation time
        INTEGER*4 :: time, time_start, time_end

! -- switch to set the format for writing xtime.dat and rtime.dat
        INTEGER*4 :: output_format

! -- netCDF variables
        INTEGER   :: ncid, varid_r, varid_x, varid_qv, varid_temp
        REAL*8, DIMENSION(:), ALLOCATABLE  :: qv_out, temp_out

! -- info for the run
        CHARACTER*70 :: run_info

! file names and paths
! -------------------------------------------------------------------------------------------------
! -- ccn properties in entrained air (c_file), drop properties (d_file), 
! -- input file names / paths for initial environmental sounding (i_file), 
! -- input file names / paths for incloud vertical velocity (v_file)
! -- input file names (given in namelist)
! -- input file paths (given in namelist)
! -- output path (given in namelist)
        CHARACTER(LEN=70) :: c_file, d_file, i_file, v_file
        CHARACTER(LEN=30) :: i_profile_file, v_profile_file, drop_file, ccn_file
        CHARACTER(LEN=50) :: i_profile_path, v_profile_path, drop_path, ccn_path
        CHARACTER(LEN=50) :: output_path

! -- file name for radii and position history data (r_filename and x_filename)
        CHARACTER(LEN=100):: r_filename, x_filename
! -- file name for netCDF file (n_filename)
        CHARACTER(LEN=100):: n_filename

! -- characters to number output files for different realizations
        CHARACTER(LEN=4)  :: format_string
        CHARACTER(LEN=2)  :: filenumber

! -------------------------------------------------------------------------------------------------
! beginning
! -------------------------------------------------------------------------------------------------

! -- read namelist
        NAMELIST /inipar/ bl, bl_inte, eta, epsilon,              &
                          ne, ne1, idimen, mi, n_max,             &
                          rho_a, aero,                            &
                          n_blob, psigma, ent_rate,               &
                          run_info,                               &
                          itd, ire, iim, idc, icw, ivp, icaw,     &
                          iinm, iccn, nccn, r_ccn, m_ccn,         &
                          idp, m_i, N_i, r_i, LWC, thsw,          &
                          qv, temp, press, w, height, cb_ent,     &
                          cloudtopheight,                         &
                          r_pdf_min, r_pdf_max, r_pdf_nbin,       &
                          qt_pdf_min, qt_pdf_max, qt_pdf_nbin,    &
                          scalar_pdf_nbin,                        &
                          qv_pdf_min, qv_pdf_max,                 &
                          temp_pdf_min, temp_pdf_max,             &
                          SuS_pdf_min, SuS_pdf_max,               &
                          start_t, end_t, t_interval, time_step,  &
                          mrealization

        NAMELIST /io_file/ output_format,                         &
                           i_profile_file, i_profile_path,        &
                           v_profile_file, v_profile_path,        &
                           drop_file, drop_path,                  &
                           ccn_file, ccn_path,                    &
                           output_path

        OPEN (31, FILE='namelist', FORM='FORMATTED')

        READ (31, NML=inipar)

        READ (31, NML=io_file)

        CLOSE (31)

! -- initialize some variables
        ent_flag       = .FALSE.
        LWP_flag       = .TRUE.
        LWP            = 0.0
        count_LWP      = 0.0


        iteration1 = t_interval/time_step
        iteration2 = end_t/time_step
        mip        = end_t/t_interval
        igrid      = 50000 !end_t/time_step

        time_start = time()

! -------------------------------------------------------------------------------------------------
! allocate arrays
! -------------------------------------------------------------------------------------------------

        allocate( dis(ne) )
        allocate( new_con(ne1), radius_drop(ne1), mass_solute(ne1), con_drop(ne1) ) 
        allocate( con_ccn(ne1), mass_ccn(ne1), radius_ccn(ne1) )
        allocate( dsize(idimen), csize(idimen), x(idimen), qwater(idimen) )
        allocate( v(idimen) )

        allocate( p_array(mi), z_array(mi), qv_array(mi), temp_array(mi) )

        allocate( spoint(n_max), fpoint(n_max), m2(n_max), m3(n_max) )

        allocate( A(0:igrid,2), B(0:igrid,2) )
        allocate( SA1(0:igrid), SA2(0:igrid) )
        allocate( jcell_pop(0:10,igrid) )

        allocate( qt(igrid), sl(igrid) )

        allocate( qlwater(igrid), SuS(igrid), hm(igrid) )
        allocate( dx_map(igrid) ) 

        allocate( index(0:igrid) )

        allocate( r_pdf(r_pdf_nbin), r_scalar(r_pdf_nbin) )
        allocate( mean_r(r_pdf_nbin), mean_r2(r_pdf_nbin), mean_r3(r_pdf_nbin) )

        allocate( qt_pdf(qt_pdf_nbin), qt_scalar(qt_pdf_nbin) )

        allocate( scalar(n_scalar,igrid), dummy_Sus(igrid) )
        allocate( scalar_mean(n_scalar,0:mip), scalar_mean_ave(n_scalar,0:mip) )
        allocate( scalar_rms(n_scalar,0:mip), scalar_rms_ave(n_scalar,0:mip) )
        allocate( SuS_pdf(mip,scalar_pdf_nbin), SuS_pdf_ave(mip,scalar_pdf_nbin) )
        allocate( SuS_scalar(scalar_pdf_nbin) )
        allocate( SA1_PDF(mip,scalar_pdf_nbin), SA2_PDF(mip,scalar_pdf_nbin) )
        allocate( SA1_PDF_ave(mip,scalar_pdf_nbin), SA2_PDF_ave(mip,scalar_pdf_nbin) )
        allocate( SA1_scalar(scalar_pdf_nbin), SA2_scalar(scalar_pdf_nbin) )

        allocate( r_std(mip), r_mean(mip) )
        allocate( s_std(mip), s_mean(mip) )
        allocate( r_max(mip), r_max_index(mip), actndrop(mip) )

        OPEN (UNIT=7,FILE='EMPM.log', FORM='FORMATTED')

        WRITE(7,*) '---------------------------------------------------------------------------'
        WRITE(7,*) '---------------------------------------------------------------------------'
        WRITE(7,*) '      _____             ____           '
        WRITE(7,*) '     |       |\    /|  |    \  |\    /|'
        WRITE(7,*) '     |       | \  / |  |    |  | \  / |'
        WRITE(7,*) '     |---    |  \/  |  |____/  |  \/  |'
        WRITE(7,*) '     |       |      |  |       |      |'
        WRITE(7,*) '     |_____  |      |  |       |      |'
        WRITE(7,*) ''
        WRITE(7,*) '---------------------------------------------------------------------------'
        WRITE(7,*) '---------------------------------------------------------------------------'
        WRITE(7,*) ''
        WRITE(7,*) run_info
        WRITE(7,*) ''
        WRITE(7,*) 'parameter read from namelist'
        WRITE(7,*) '---------------------------------------------------------------------------'
        WRITE(7,*) 'bl      = ', bl
        WRITE(7,*) 'bl_inte = ', bl_inte
        WRITE(7,*) 'eta     = ', eta

! -------------------------------------------------------------------------------------------------
! read spline points from splinepoints.input, the data is then passed to fcnkb and psplint through 
! intspline.f90
! -------------------------------------------------------------------------------------------------

        OPEN (51, file='splinepoints.input', access='sequential', status='old', form='formatted')

        READ(51,*) (((splinearray(i,j,k),i=1,ncats),j=1,nvecs),k=1,npoints)

        CLOSE (51)

! -------------------------------------------------------------------------------------------------
! read initial environmental profile from file
! -------------------------------------------------------------------------------------------------
        i_file = TRIM(i_profile_path)//TRIM(i_profile_file)
        WRITE(7,*) ''
        WRITE(7,*) 'reading initial environmental profile from file:'
        WRITE(7,*) i_file
        WRITE(7,*) '---------------------------------------------------------------------------'
        WRITE(7,*) '  pressure        height          mixing ratio  temperature'
        WRITE(7,*) '   (hPa)           (m)             (kg/kg)       (K)'

        OPEN (51,FILE=i_file)

        DO i = 1,mi
           READ(51,*) p_array(i), z_array(i), qv_array(i), temp_array(i)
           WRITE(7,'(4(3X,F12.6))') p_array(i), z_array(i), qv_array(i), temp_array(i)
        END DO

        CLOSE(51)

        IF (ivp .EQ. 1) THEN

           v_file = TRIM(v_profile_path)//TRIM(v_profile_file)
           nline_v = count_lines(v_file)
           allocate( p_cloud(nline_v), h_cloud(nline_v), w_cloud(nline_v) ) 

           WRITE(7,*) ''
           WRITE(7,*) 'reading incloud vertical velocity profile from file:'
           WRITE(7,*) v_file
           WRITE(7,*) '---------------------------------------------------------------------------'
           WRITE(7,*) '   pressure      height            velocity'
           WRITE(7,*) '    (hPa)         (m)               (m/s) '
           OPEN (52,FILE=v_file)

           DO i = 1,nline_v
              READ(52,*) p_cloud(i), h_cloud(i), w_cloud(i)
              WRITE(7,'(3(3X,F12.6))')  p_cloud(i), h_cloud(i), w_cloud(i)
           END DO

           CLOSE(52)
           w = w_cloud(1)
        END IF

        WRITE(7,*) ''

        WRITE(7,*) 'cloud base properties (given in namelist):'
        WRITE(7,*) 'qv     = ', qv, '   ! water vapor mixing ratio (kg/kg)'
        WRITE(7,*) 'temp   = ', temp, '   ! temperature (K)'
        WRITE(7,*) 'press  = ', press, '   ! pressure (Pa)'
        WRITE(7,*) 'w      = ', w, '   ! cloud base vertical velocity (m/s)'
        WRITE(7,*) 'rho    = ', rho_a, '   ! cloud base air density (kg/m^3)'
        WRITE(7,*) 'height = ', h_cloud(1), '   ! height (m)'
        WRITE(7,*) '---------------------------------------------------------------------------'
        WRITE(7,*) ''

! -- store initial pressure, temperature, height, velocity and water vapor mixing ratio
        press0  = press
        temp0   = temp
        height0 = height
        w0      = w
        qv0     = qv

! -- calculate the volume of the domain (the cross section size is set to roughly the Kolomogorov 
! -- scale)
        volume  = 0.001*0.001*bl

! -- calculating Reynolds number (Re)
        Re      = (bl_inte/eta)**(4.0/3.0)

! -- set time to zero
        t       = 0.0

! -- calculate the number of grid points
! -- (original fomula for fast mixing (small eta))
        ngrid   = NINT(BL_inte/(eta/6.0))*(BL/BL_inte)
! -- (new fomula for slow mixing (big eta))
!       ngrid = NINT(BL_inte/(eta/60.0))*(BL/BL_inte)

        WRITE(7,*) 'ngrid = ', ngrid

! -- calculate the grid size
        dx      = (bl/FLOAT(ngrid))

! -- calculate the turbulent diffusivity
        dtur    = 0.1*BL_inte**(4.0/3.0)*epsilon**(1.0/3.0)

! -- calculate the molecular diffusivity (for temperature (dtemp) and water vapor (dm))
        dm      = dtur/Re
        dtemp   = dm

! -- calculate time step for diffision (td)
        td      = (0.2*(dx**2.0)/dm)

        IF (itd .EQ. 1) THEN
           td = time_step
        END IF

        WRITE(7,*) 'just determined time step for diffusion:'
        WRITE(7,*) 'td = ', td

! -- calculate time step for turbulent convection (tc)
! -- Krueger (1993), 'Linear Eddy Modeling of Entrainment and Mixing in Stratus Clouds', (JAS)
        tau_L   = BL_inte**2.0/dtur
        tlambda = (54/5.0*Re**1.25)/(BL_inte*tau_L)
        tc      = 1.0/(BL*tlambda)

        WRITE(7,*) ''
        WRITE(7,*) 'determined time steps '
        WRITE(7,*) '---------------------------------------------------------------------------'

! -- compare td and tc to decide how many triplet mapping per diffusion
! -- decide which time step to use
        IF (td .GE. tc) THEN
           mtrip  = INT(td/tc)
           idtrip = 1
           dt     = td
           WRITE(7,*) 'td is greater / equal to tc'
        ELSE
           mtrip  = 1
           idtrip = INT(tc/td)+1
           dt     = tc/FLOAT(idtrip)
           WRITE(7,*) 'td is smaller than tc, so let td = tc/idtrip'
        END IF

        WRITE(7,*) 'td     = ', td
        WRITE(7,*) 'tc     = ', tc 
        WRITE(7,*) 'dt     = ', dt
        WRITE(7,*) 'dx     = ', dx
        WRITE(7,*) 'mtrip  = ', mtrip
        WRITE(7,*) 'idtrip = ', idtrip

! -- check whether the interval (t_interval) for writing out history data is n*dt where 
! -- n is an integer
        IF (mod(t_interval/dt, 1.0) .GT. 0) THEN
           WRITE(7,*) mod(t_interval/dt, 1.0)
           STOP 't_interval devide by dt is not an integer of dt'
        END IF

        IF (mod((start_t/t_interval), 1.0) .GT. 0) THEN
           WRITE(7,*) mod((start_t/t_interval), 1.0)
           STOP 'start_t divided by t_interval does not give an integer'
        END IF

        IF (mod((end_t/t_interval), 1.0) .GT. 0) THEN
           WRITE(7,*) mod((end_t/t_interval), 1.0)
           STOP 'end_t divided by t_interval does not give an integer'
        END IF

! -- calculate iteration index for writing out history data for t_interval (iteration3)
        iteration3  = INT( t_interval/dt ) 

! -- calculate start and end iteration index
        iteration_s = start_t/dt
        iteration_e = end_t/dt

! -- find the array length for history data
        hist_array  = (end_t-start_t)/t_interval

        WRITE(7,*) ''
        WRITE(7,*) 'time information for recording history data '
        WRITE(7,*) '---------------------------------------------------------------------------'
        WRITE(7,*) 'start time    = ', start_t
        WRITE(7,*) 'end time      = ', end_t
        WRITE(7,*) 'time interval = ', t_interval
        WRITE(7,*) 'array length of hist. data = ', hist_array
        WRITE(7,*) 'start iteration index      = ', iteration_s
        WRITE(7,*) 'end iteration index        = ', iteration_e   

! -- check if length of the history array exceed the allowed value
        IF (hist_array .GT. mip) THEN
           STOP 'history data array exceeds max.'
        END IF  

        N_i = N_i*1.0E6

! -- define length of r_time, x_time and other arrays
        IF (iccn .EQ. 0) THEN
           m_index = NINT(volume*N_i)
        ELSE 
           delta_ent     = (n_blob/ent_rate)*(psigma/(1.0-psigma))
           cloud_depth   = h_cloud(nline_v)-h_cloud(1)

           number_ent    = cloud_depth/delta_ent

           n_drop_ent    = n_blob*psigma*bl*volume*N_i/bl
           tot_num_index = number_ent*n_drop_ent

           m_index = NINT((volume*N_i)+1.5*tot_num_index)
        END IF

        WRITE(7,*) 'length of array for droplet data = ', m_index

        ALLOCATE( r_time(m_index), x_time(m_index) )
        ALLOCATE( time_null(m_index), index_drop(m_index) )

! -------------------------------------------------------------------------------------------------
! Loop for the realization starts here, the number of realization is set in the namelist
! -------------------------------------------------------------------------------------------------
        DO Mreal=1,Mrealization 

           IF (mreal .LT. 10) THEN
              format_string = "(i1)"
           ELSE IF (mreal .LT. 100) THEN
              format_string = "(i2)"
           ELSE
              WRITE(7,*) 'The number of realizations is too big. It is limited to 99 at the moment.'
              STOP
           END IF

           WRITE(filenumber, format_string) mreal

           allocate ( qv_out(ngrid), temp_out(ngrid) )

! -- files for the whole domain over time / each droplet over time
           r_filename = TRIM(output_path)//'rtime_'//TRIM(filenumber)//'.dat'  
           x_filename = TRIM(output_path)//'xtime_'//TRIM(filenumber)//'.dat'
           IF (output_format .EQ. 1) THEN
              OPEN(61,FILE=TRIM(r_filename),FORM='formatted')
              OPEN(62,FILE=TRIM(x_filename),FORM='formatted')
           ELSE IF (output_format .EQ. 2) THEN
              OPEN(61,FORM='unformatted',FILE=TRIM(r_filename))
              OPEN(62,FORM='unformatted',FILE=TRIM(x_filename))
           ELSE IF (output_format .EQ. 3) THEN
              n_filename = TRIM(output_path)//'EMPM_output_'//TRIM(filenumber)//'.nc'
              CALL netcdf_prep(TRIM(n_filename),m_index,ngrid,mip,dt,dx,iteration1, &
                               varid_r,varid_x,varid_qv,varid_temp,ncid)
           END IF 
           OPEN(63,FILE=TRIM(output_path)//'qvtime_'//TRIM(filenumber)//'.dat',FORM='formatted')
           OPEN(64,FILE=TRIM(output_path)//'temptime_'//TRIM(filenumber)//'.dat',FORM='formatted')
! -- control files / interesting information
           OPEN(71,FILE=TRIM(output_path)//'qt_pdf_time_'//TRIM(filenumber)//'.dat',FORM='formatted')
           OPEN(72,FILE=TRIM(output_path)//'findex_'//TRIM(filenumber)//'.dat',FORM='formatted')
           OPEN(73,FILE=TRIM(output_path)//'entm_'//TRIM(filenumber)//'.dat',FORM = 'formatted')
! -- averaged droplet data
           OPEN(81,FILE=TRIM(output_path)//'r_mean_stats_'//TRIM(filenumber)//'.dat',FORM='formatted')
           OPEN(82,FILE=TRIM(output_path)//'r_mean_time_'//TRIM(filenumber)//'.dat',FORM='formatted')
           OPEN(83,FILE=TRIM(output_path)//'r2_mean_time_'//TRIM(filenumber)//'.dat',FORM='formatted')
           OPEN(84,FILE=TRIM(output_path)//'r3_mean_time_'//TRIM(filenumber)//'.dat',FORM='formatted')
           OPEN(85,FILE=TRIM(output_path)//'r_pdf_time_'//TRIM(filenumber)//'.dat',FORM='formatted')
           OPEN(86,FILE=TRIM(output_path)//'r_eff_time_'//TRIM(filenumber)//'.dat',FORM='formatted')
           OPEN(87,FILE=TRIM(output_path)//'super_mean_'//TRIM(filenumber)//'.dat',FORM='formatted')
           OPEN(88,FILE=TRIM(output_path)//'super_max_min_'//TRIM(filenumber)//'.dat',FORM='formatted')

! -- domain average / average of all realizations
           OPEN(91,FILE=TRIM(output_path)//'ave_'//TRIM(filenumber)//'.dat',FORM = 'formatted')
           OPEN(92,FILE=TRIM(output_path)//'SA1_pdf_ave_'//TRIM(filenumber)//'.dat',FORM='formatted')
           OPEN(93,FILE=TRIM(output_path)//'SA2_pdf_ave_'//TRIM(filenumber)//'.dat',FORM='formatted')
           OPEN(94,FILE=TRIM(output_path)//'temp_qv_ave_'//TRIM(filenumber)//'.dat',FORM='formatted')
           OPEN(95,FILE=TRIM(output_path)//'sus_pdf_'//TRIM(filenumber)//'.dat',FORM='formatted')
           OPEN(96,FILE=TRIM(output_path)//'LWP_'//TRIM(filenumber)//'.dat',FORM='formatted')

           WRITE(7,*) '---------------------------------------------------------------------------'
           WRITE(7,*) 'mreal = ', mreal

! -- create random number seed based on computer time
! -- be careful: time() is an intrinsic which is not fully portable! The values returned by this 
! -- intrinsic might be, or become, negative. (GCC, the GNU Compiler Collection)
           iseed1 = time()
           iseed2 = time()+10000000
           iseed3 = time()+20000000
           iseed4 = time()+30000000
           iseed5 = time()+40000000
           iseed6 = time()+50000000

           WRITE(7,*) ''
           WRITE(7,*) 'random number seeds:'
           WRITE(7,*) '---------------------------------------------------------------------------'
           WRITE(7,*) 'iseed1 = ', iseed1
           WRITE(7,*) 'iseed2 = ', iseed2
           WRITE(7,*) 'iseed3 = ', iseed3
           WRITE(7,*) 'iseed4 = ', iseed4
           WRITE(7,*) 'iseed5 = ', iseed5
           WRITE(7,*) 'iseed6 = ', iseed6

! -- 're'-initialize several variables
           ent_flag  = .FALSE.
           LWP_flag  = .TRUE.
           LWP       = 0.0
           count_LWP = 0.0

           press     = press0
           temp      = temp0
           height    = height0
           qv        = qv0
           w         = w0

           ip        = 0
           num_tot   = 0

           t         = 0.0

! -- calculate the inverse volume of each grid cell
           i_grid_vol = 1.0/(dx*0.001*0.001)

           IF (idp .EQ. 1) THEN
! -- all droplets have the same initial properties (given in namelist)
              DO j=1,INT(BL*N_i*1.0E-6)
                 num_tot        = num_tot+1
                 dsize(num_tot) = r_i
                 csize(num_tot) = m_i
              END DO
              n_drop = num_tot
              WRITE(7,*) ''
              WRITE(7,*) 'All droplets have the same initial properties (given in namelist)'
              WRITE(7,*) '---------------------------------------------------------------------------'
              WRITE(7,*) 'initial radius     = ', r_i
              WRITE(7,*) 'mass of the solute = ', m_i
              WRITE(7,*) 'n_drop             = ', n_drop
              WRITE(7,*) ''
           ELSE IF (idp .EQ. 2) THEN
! -- calculate droplet size distribution using gamma function and number concentration (N_i) and 
! -- LWC given in namelist
              CALL size_dis(N_i,LWC,rho_a,radius_drop,con_drop,mass_solute,volume,m_i,ne)
              WRITE(7,*) '---------------------------------------------------------------------------'
              WRITE(7,*) ' droplet radius  solute mass     droplet concentration'
              WRITE(7,*) '  (m)             (kg)            (#/domain volume)' 
              DO i=1,ne
                 IF (NINT(con_drop(i)) .NE. 0) THEN
                    DO j=1,NINT(con_drop(i))
                       num_tot        = num_tot+1
                       dsize(num_tot) = radius_drop(i)
                       csize(num_tot) = mass_solute(i)
                    END DO
                    WRITE(7,'(3x,3(e12.5,4x))') radius_drop(i), mass_solute(i), con_drop(i)
                 END IF
              END DO
              n_drop = num_tot
              WRITE(7,*) 'n_drop = ', n_drop
           ELSE IF (idp .EQ. 3) THEN
! -- read droplet data from file (the output of Dr. Austin's DGM code)
             d_file = TRIM(drop_path)//TRIM(drop_file)
             WRITE(7,*) ''
             WRITE(7,*) 'reading droplet data from file:'
             WRITE(7,*) d_file
             WRITE(7,*) '---------------------------------------------------------------------------'
             WRITE(7,*) ' droplet radius'
             WRITE(7,*) '  (m)' 

	     OPEN (52,FILE=d_file)

             DO i=1,ne
                READ(52,*) radius_drop(i)
                WRITE(7,*) radius_drop(i)
 	     END DO

	     READ(52,*)(mass_solute(i),con_drop(i),i=1,ne)

	     WRITE(7,*) '---------------------------------------------------------------------------'
             WRITE(7,*) ' solute mass    droplet concentration'
             WRITE(7,*) '  (kg)           (#)' 
	     WRITE(7,'(2(3x,e12.5))')(mass_solute(i),con_drop(i),i=1,ne)

	     CLOSE(52)

! -- determine the new droplet distribution of the present domain size
             CALL accu(volume,dis,con_drop,n_drop,ne,ne1)
             WRITE(7,*) 'n_drop = ', n_drop

             DO i=1,ne
                IF (INT(dis(i)) .NE. 0.0) THEN
                   DO j=1,INT(dis(i))
                      num_tot        = num_tot + 1
                      dsize(num_tot) = radius_drop(i)
                      csize(num_tot) = mass_solute(i)
                   END DO
                END IF
             END DO
           END IF

! -- calculate the initial total liquid water mixing ratio, for given set of droplets
           ql0_vol = 0.0

           DO k=1,n_drop
              ql0_vol = ql0_vol+(pi43*rho_w*dsize(k)**3*i_grid_vol)
           END DO

           ql0     = ql0_vol/ngrid
           qt0     = qv0+ql0

           WRITE(7,*) '---------------------------------------------------------------------------'
           WRITE(7,*) 'initial water vapor, liquid water and total water mixing ratio:'
           WRITE(7,*) 'qv = ', qv0
           WRITE(7,*) 'ql = ', ql0
           WRITE(7,*) 'qt = ', qt0
           WRITE(7,*) ''
	   WRITE(7,*) '---------------------------------------------------------------------------'

           m1 = 1
           m4 = ngrid
           m6 = 2*ngrid

! -- randomly assigns the droplet position
! -- x(i) are the initial positions of the droplets
           CALL assigd(n_drop,ngrid,jcell_pop,dx,m1,x,bl,igrid,idimen,iseed4)

	   DO i=1,n_drop
	      index_drop(i) = i
	   END DO
	   n_drop_index = n_drop

           t_entm    = 0.0                          ! total time elapsed
           iteration = 0

! -- finding environmental temperature and water vapor mixing ratio at cloud base
           CALL intpol(mi,p_array,qv_array,press,qv_e) !dt,w,
           CALL intpol(mi,p_array,temp_array,press,temp_e) !dt,w,

! -- write initial microphysical properties
	   WRITE(94,100) iteration,press0,qv0,temp0,ql0,qt0,w,temp_e,qv_e
100	   format(1x,i5,1x,f12.5,1x,f12.9,1x,f9.4,1x,f16.13,1x,f16.13,1x,f14.9,1x,f9.4,1x,f12.9)

! -- write initial radius statistics (mean and standard deviation)
           r_sum      = 0.0
           r_std_sum  = 0.0
           n_act_drop = n_drop
           DO i=1,n_drop
              r_sum = r_sum+dsize(i)*1.0e6
           END DO

           r_sum = r_sum/FLOAT(n_act_drop)

           DO i=1,n_drop
              r_std_sum = r_std_sum+(dsize(i)*1.0e6-r_sum)**2.0
           END DO

           r_std_sum    = ((r_std_sum/FLOAT(n_act_drop))**0.5)
           WRITE(81,101) t, r_sum, r_std_sum, n_act_drop
101        format(f12.4,2(1x,f12.7),1(1x,i4))

! -------------------------------------------------------------------------------------------------
! Loops for adiabatic lifting of model domain (and eventual entrainment events)
! -------------------------------------------------------------------------------------------------
           DO                                     ! exit if (iteration .ge. iteration2)
              DO                                  ! exit if (mod(iteration,iteration1) .eq. 0)
                 time0 = dt*iteration

! -- end simulation when top of environmental profile is reached
                 IF (press .LT. p_array(SIZE(p_array))*100.0) THEN
                    WRITE(7,*) 'The top of the environmental profile was exceeded. Leaving the code now.'
                    exit
                 END IF
! -- end simulation when top of incloud vertical velocity profile is reached
                 IF (ivp .EQ. 1) THEN
                    IF (press .LT. p_cloud(SIZE(p_cloud))*100.0) THEN
                       WRITE(7,*) 'The top of the incloud vertical velocity profile was exceeded. '
                       WRITE(7,*) 'Leaving the code now.'
                       exit
                    END IF
                 END IF
! -- set flag for entrainment at cloud base (first time step) to TRUE
                 IF (cb_ent .EQ. 1) THEN
                    ent_flag = .TRUE.
                 END IF
! -- interpolate height, mixing ration and temperature to pressure level
! -- interpolated value are stored on Z, qv_e and temp_e
                 IF (ivp .EQ. 1) THEN
	            CALL intpol(nline_v,p_cloud,w_cloud,press,vel) !time0,w,
                    w = vel
                 END IF

	         CALL intpol(mi,p_array,z_array,press,Z) !dt,w,
                 CALL intpol(mi,p_array,qv_array,press,qv_e) !dt,w,
                 CALL intpol(mi,p_array,temp_array,press,temp_e) !dt,w,

                 IF (t .GE. t_entm) THEN
                    WRITE(7,*) 'iteration = ', iteration

	            CALL fdtime(psigma,n_blob,ent_rate,iseed3,w,ire,dt_entm)

	            time0  = t_entm
	            t_entm = t_entm + dt_entm

                    WRITE(73,102) time0, press, Z, qv_e, temp_e, dt_entm
102                 format(6(2x,e15.6))

! -- set velocity to zero after first entrainment event to study isobaric mixing 
                    IF (iim .EQ. 1) THEN
                       IF (time0 .GE. dt_entm) THEN
                          w = 0.0
                       END IF
                    END IF

                    WRITE(7,*) ''
                    WRITE(7,*) 'vertical velocity (w)             = ', w
                    WRITE(7,*) 'interpolated height (Z)           = ', Z
                    WRITE(7,*) 'interpolated mixing ratio (qv_e)  = ', qv_e
                    WRITE(7,*) 'interpolated temperature (temp_e) = ', temp_e
 
! -- find spoint and fpoint to decide locations of m1, m2, m3.....
                    CALL domain(psigma,ngrid,spoint,fpoint,iseed6,n_max,n_blob,n_blob_final,xn)

                    DO i=1,n_blob_final
                       m2(i) = spoint(i)
                       m3(i) = fpoint(i)
                    END DO

                    WRITE(7,*) 'm2(i) = spoint(i) and m3(i) = fpoint(i)'
                    WRITE(7,*) 'm4 = ', m4
                    WRITE(7,*) 'm6 = ', m6
	            
                    DO i=m1,m6
	               SuS(i) = 0.0
	            END DO

! -------------------------------------------------------------------------------------------------
! Begin entrainment
! -------------------------------------------------------------------------------------------------
                    IF (cb_ent .EQ. 0 .AND. iteration .NE. 0) THEN
                       ent_flag = .TRUE.
                    END IF

                    IF (ent_flag) THEN
                       DO j=1,ngrid
                          jcell_pop(0,j) = 0
                       END DO

                       ncut = 0
                       DO k=1,n_drop
                          ichange   = k
                          i         = k-ncut
	                  j         = INT(x(k)/dx)+m1

                          x(i)      = x(ichange)
                          dsize(i)  = dsize(ichange)
                          csize(i)  = csize(ichange)
	                  qwater(i) = qwater(ichange)

	                  index_drop(i) = index_drop(ichange)

                          DO i = 1, n_blob_final
                             IF ((j .GE. m2(i) ) .AND. (j .LE. m3(i)) ) THEN
                                ncut = ncut+1
                                GOTO 622
                             END IF
                          END DO
                          jcell_pop(0,j-m1+1)                   = jcell_pop(0,j-m1+1)+1
                          jcell_pop(jcell_pop(0,j-m1+1),j-m1+1) = i
 622                      CONTINUE
                       END DO
                       n_drop = n_drop-ncut

                       DO k=1,m_index
                          IF (k .GT. n_drop) THEN
                             index_drop(k) = 0
                          END IF
                       END DO

                       IF (iccn .EQ. 0) THEN
! -- insert blob of clean air (contains no ccn)
                          WRITE(7,*) ''
                          WRITE(7,*) 'The entrained air contains no ccn!'
                          WRITE(7,*) '---------------------------------------------------------------------------'
                       ELSE IF (iccn .EQ. 1) THEN
! -- entrained air contains droplets / ccn
! -- all ccn have the same radius and solute mass (given in namelist)
                          WRITE(7,*) ''
                          WRITE(7,*) 'All ccn have the same radius and solute mass (given in namelist)'
                          WRITE(7,*) '---------------------------------------------------------------------------'
                          WRITE(7,*) ' ccn radius (m) = ', r_ccn
                          WRITE(7,*) ' ccn mass (kg)  = ', m_ccn

                          volume_entm = psigma*volume
                          num_tot     = n_drop
	                  n_drop_used = n_drop

                          IF (idc .EQ. 1) THEN
! -- constant droplet number in domain
                             n_ccn    = ncut
                          ELSE
                             n_ccn    = volume_entm*N_i
                          END IF

	                  WRITE(7,*) ' n_ccn          = ', n_ccn
                          WRITE(7,*) '---------------------------------------------------------------------------'
                          DO j=1,INT(n_ccn)
                             num_tot         = num_tot+1

                             n_drop_index    = n_drop_index+1
                             dsize(num_tot)  = r_ccn
                             csize(num_tot)  = m_ccn
                             qwater(num_tot) = pi43*(dsize(num_tot)**3.0)*rho_w*i_grid_vol

	                     index_drop(num_tot) = n_drop_index
                          END DO

                          WRITE(7,103) n_drop_index, dsize(num_tot), qwater(num_tot)
103                       format(i6,2(3x,e12.5))
                          WRITE(7,*) '---------------------------------------------------------------------------'

                          n_drop = n_drop+n_ccn

                          blob_flag = 0
                          IF (n_blob .LT. n_blob_final) blob_flag = 1
                          part = NINT((FLOAT(m3(n_blob))-FLOAT(m2(n_blob)))/FLOAT(xn)*10)

! -- assign the ccn into the portions between m2 and m3
	                  CALL assigd2(n_drop_used,n_ccn,m2,m3,jcell_pop,dx,m1,x,igrid,idimen,iseed5,n_max, & 
                               &       n_blob,part,blob_flag,ngrid)
                       ELSE IF (iccn .EQ. 2) THEN
! -- entrained air contains droplets / ccn
! -- read information about entrained ccns from file
                          c_file = TRIM(ccn_path)//TRIM(ccn_file)
                          WRITE(7,*) ''
                          WRITE(7,*) 'reading ccn data in entrained air:'
                          WRITE(7,*) c_file
                          WRITE(7,*) '---------------------------------------------------------------------------'
                          WRITE(7,*) ' ccn radius'
                          WRITE(7,*) '  (m)' 

	                  OPEN (53,FILE=c_file)

                          DO i=1,nccn
                             READ(53,*) radius_ccn(i)
                             WRITE(7,*) radius_ccn(i)
 	                  END DO
                          READ(53,*)(mass_ccn(i),con_ccn(i),i=1,nccn)

                          WRITE(7,*) '---------------------------------------------------------------------------'
                          WRITE(7,*) ' ccn mass       ccn concentration'
                          WRITE(7,*) '  (kg)           (#/m^3)' 
	                  WRITE(7,'(2(3x,e12.5))')(mass_ccn(i),con_ccn(i),i=1,nccn)

                          CLOSE(53)

	                  n_drop_used = n_drop
                          volume_entm = psigma*volume

! -- determine the ccn distribution of entrained blob
                          CALL accu(volume_entm,dis,con_ccn,n_ccn,ne,ne1)

                          num_tot = n_drop

	                  WRITE(7,*) 'n_ccn = ', n_blob*n_ccn
                          WRITE(7,*) '---------------------------------------------------------------------------'
                          WRITE(7,*) 'index  droplet/ccn    liquid water  '
                          WRITE(7,*) '        radius         mixing ratio '
                          WRITE(7,*) '        (m)            (kg/kg)      ' 

                          dis = dis*n_blob

                          DO i=1,ne
                             IF (INT(dis(i)) .GT. 0) THEN
                                DO j=1,INT(dis(i))
                                   num_tot         = num_tot+1 

                                   n_drop_index    = n_drop_index+1
                                   dsize(num_tot)  = radius_ccn(i)
                                   csize(num_tot)  = mass_ccn(i)
                  	           qwater(num_tot) = pi43*(dsize(num_tot)**3.0)*rho_w*i_grid_vol

	                           index_drop(num_tot) = n_drop_index
                                END DO
                             END IF
                          END DO

                          WRITE(7,104) n_drop_index, dsize(num_tot), qwater(num_tot)
104                       format(i6,2(3x,e12.5))
                          WRITE(7,*) '---------------------------------------------------------------------------'

	                  n_ccn  = n_ccn*n_blob
                          n_drop = n_drop+n_ccn

                          blob_flag = 0
                          IF (n_blob .LT. n_blob_final) blob_flag = 1
                          part = NINT((FLOAT(m3(n_blob))-FLOAT(m2(n_blob)))/FLOAT(xn)*10)

! -- assign the ccn into the portions between m2 and m3
	                  CALL assigd2(n_drop_used,n_ccn,m2,m3,jcell_pop,dx,m1,x,igrid,idimen,iseed5,n_max, & 
                               &       n_blob,part,blob_flag,ngrid)
                       ELSE
                          WRITE(7,*) ''
                          WRITE(7,*) 'iccn has to be either 0,1 or 2! Check namelist!'
                          WRITE(7,*) 'Run aborted'
                          WRITE(7,*) '---------------------------------------------------------------------------'
                          STOP
                       END IF
                    END IF

                    WRITE(7,*) 'n_drop = ', n_drop
                    WRITE(7,*) '---------------------------------------------------------------------------'

                    IF (iteration .EQ. 0) THEN
                       DO i=m1,m4
                          A(i,1) = qv
                          A(i,2) = temp
                       END DO

                       DO i=m4+1,m6
                          A(i,1) = A(i-ngrid,1)
                          A(i,2) = A(i-ngrid,2)
                       END DO
                    ELSE
                       WRITE(7,*) '          k       m2(k)       m3(k)'
	               WRITE(7,*) '---------------------------------------------------------------------------'
                       DO k=1,n_blob_final
                          WRITE(7,*) k, m2(k), m3(k)
                          DO i=m2(k),m3(k)
                             A(i,1) = qv_e
                             A(i,2) = temp_e
                          END DO
                       END DO

                       WRITE(7,*) '---------------------------------------------------------------------------'

                       DO i=m4+1,m6
                          A(i,1) = A(i-ngrid,1)
                          A(i,2) = A(i-ngrid,2)
                       END DO
                    END IF

! -- set initial value for jcell_pop(0,j) 
                    DO j=m1,m4
                       jcell_pop(0,j-m1+1) = 0
                    END DO

                    DO i=1,n_drop
                       j                                     = INT(x(i)/dx+m1)
                       jcell_pop(0,j-m1+1)                   = jcell_pop(0,j-m1+1)+1
                       jcell_pop(jcell_pop(0,j-m1+1),j-m1+1) = i
                    END DO

! -- write out the initial info for the run
!                    IF (iteration .EQ. 0) THEN
!                       WRITE(63,*) run_info
!                       WRITE(63,*) ngrid, dx 
!                       WRITE(63,*) iteration*dt
!                       WRITE(63,*) (A(i,1), i=1,ngrid)
!
!                       WRITE(64,*) run_info
!                       WRITE(64,*) ngrid, dx 
!                       WRITE(64,*) iteration*dt
!                       WRITE(64,*) (A(i,2), i=1,ngrid)
!                    END IF

! -- flag is set to zero to set the I.C.s and B.C.s in the drop_map subroutine
                    flag = 0
                 END IF

                 IF (iteration .EQ. 0) THEN
                    DO i=1,n_drop
                       r_time(index_drop(i)) = dsize(i)
                       x_time(index_drop(i)) = x(i)
                    END DO

                    IF (output_format .EQ. 1) THEN
                       WRITE(61,*) iteration*dt
                       WRITE(61,*) (r_time(i), i=1,m_index)
 
                       WRITE(62,*) iteration*dt
                       WRITE(62,*) (x_time(i), i=1,m_index)
                      
                       WRITE(63,*) run_info
                       WRITE(63,*) ngrid, dx 
                       WRITE(63,*) iteration*dt
                       WRITE(63,*) (A(i,1), i=1,ngrid)

                       WRITE(64,*) run_info
                       WRITE(64,*) ngrid, dx 
                       WRITE(64,*) iteration*dt
                       WRITE(64,*) (A(i,2), i=1,ngrid)
                    ELSE IF (output_format .EQ. 2) THEN
                       WRITE(61) iteration*dt
                       WRITE(61) (r_time(i), i=1,m_index)

                       WRITE(62) iteration*dt
                       WRITE(62) (x_time(i), i=1,m_index)

                       WRITE(63) run_info
                       WRITE(63) ngrid, dx 
                       WRITE(63) iteration*dt
                       WRITE(63) (A(i,1), i=1,ngrid)

                       WRITE(64) run_info
                       WRITE(64) ngrid, dx 
                       WRITE(64) iteration*dt
                       WRITE(64) (A(i,2), i=1,ngrid)
                    ELSE IF (output_format .EQ. 3) THEN
                       ip = iteration/iteration1
                       DO i=m1,m4
                          qv_out(i)   = A(i,1)
                          temp_out(i) = A(i,2)
                       END DO
                       !CALL netcdf_write(r_time,x_time,qv_out,temp_out,m_index,ngrid,ip, &
                       !                  varid_r,varid_x,varid_qv,varid_temp,ncid)
                    END IF

                    CALL radius_pdf(n_drop,r_pdf_min,r_pdf_max,r_pdf_nbin,dsize, &
                                    mean_r,mean_r2,mean_r3,r_pdf,r_scalar)
                    WRITE(82,*) mean_r
                    WRITE(83,*) mean_r2
                    WRITE(84,*) mean_r3
                    WRITE(85,*) r_scalar
                    WRITE(85,*) r_pdf

                    CALL eff_radius(n_drop,r_pdf_nbin,mean_r2,mean_r3,r_pdf,eff_r)
                    WRITE(86,*) t, Z, eff_r

 	            DO i=m1,m4
                       qt(i) = A(i,1)+qlwater(i)
 	            END DO

 	            CALL qt_prob(m1,m4,qt_pdf_min,qt_pdf_max,qt,qt_pdf,qt_scalar,igrid,qt_pdf_nbin)
	            WRITE(71,*) qt_scalar
	            WRITE(71,*) qt_pdf

! -- calculate supersaturation based on droplets
                    s_sum     = 0.0
                    s_std_sum = 0.0

                    DO i=m1,m4
                       es       = EW(A(i,2))*100.0
                       qs       = eps*es/(press-es)
                       SuS(i)   = A(i,1)/qs-1.0
                    END DO

                    DO i=1,n_drop
                       j_cell = INT(x(i)/dx)+m1
                       s_sum  = s_sum+SuS(j_cell)
                    END DO

                    s_sum   = s_sum/FLOAT(n_drop)

                    DO i=1,n_drop
                       j_cell    = INT(x(i)/dx)+m1
                       s_std_sum = s_std_sum+(SuS(j_cell)-s_sum)**2.0
                    END DO

                    s_std_sum  = ((s_std_sum/FLOAT(n_drop))**0.5)
                    WRITE(87,106) t, s_sum, s_std_sum

                 END IF

	         iteration = iteration + 1

! -- perform triplet mapping; note, sometime several mapping occurs per diffusion time step
! -- dl       random length scale for triple mapping
! -- p2       random position for triplet mapping
! -- i1       first point of the triplet mapping eddy
! -- nn       last point of the triplet mapping eddy
! -- n        grid points inside the random eddy size (note that n must be able devided by 
!             3 completely)
	         IF (mod(iteration,idtrip) .EQ. 0) THEN
                    DO j=1,mtrip
	               r1 = RAND1(iseed1)
                       dl = (r1*(BL_inte**(-5./3.)-eta**(-5./3.))+eta**(-5./3.))**(-3./5.)
                       ni = INT(dl/DX)
                       n  = NINT(FLOAT(ni)/3.)*3

	               r2 = RAND1(iseed2)
                       p2 = r2*ngrid

                       i1 = INT(p2)+m1
                       nn = n+i1-1

                       CALL triplet(n,i1,A,igrid)

! -- calulate the droplet new location after triplet mapping and sedimentation i.e., new x(i)
! -- note that n, i1, X are eddy grid points, first point and droplet positions
! -- jcell_pop is the index for droplet position
                       CALL drop_map(n,i1,m1,ngrid,flag,index,0,igrid)
                       flag = 1

! -- reestablish the B.C.s#1, 
! -- during the triplet mapping, the scalar value outside the domain boundary might be changed. 
! -- because of the periodic B.C.s we need to built the B.C.s after each mapping

                       DO i=1,n_SA
                          DO k=i1,nn 
                             IF (k .GT. m4)THEN
                                A(k-ngrid,i) = A(k,i)
	                     ELSE
	                        A(k+ngrid,i) = A(k,i)
                             END IF
                          END DO
                       END DO
                    END DO

! -- calculate dx_map of j_cell
                    DO j=1,ngrid
                       dx_map(index(j)) = (j-index(j)) * dx
                       IF(dx_map(index(j)) .LT. (-1*BL)) THEN
                          WRITE(7,*) j, dx_map(index(j))
                       END IF
                    END DO

! -- dx_map(j) is changed in the location of cell j by triplet mapping.
! -- move drops in each cell by this amount.
                    DO j=m1,m4
                       IF (jcell_pop(0,j-m1+1) .GT. 0) THEN
                          DO k=1,jcell_pop(0,j-m1+1)                        ! for all drops in cell
                             i    = jcell_pop(k,j-m1+1)
                             x(i) = x(i)+dx_map(j-m1+1)
                             IF ((x(i) .GT. BL) .OR. (x(i) .LT. 0.d0)) THEN
                                WRITE(7,*)
		                WRITE(7,*) 'jcell= ', j
                                WRITE(7,105) i, x(i), dx_map(j-m1+1)
105                             format(i6,2(2x,e27.15))
 		                WRITE(7,*) 'x(i) out of domain'
                             END IF
                          END DO
                       END IF
                    END DO

! -- set dx_map to zero again
                    DO j=1,ngrid
                       dx_map(index(j)) = 0.0
                       index(j)         = j
                       index(j+ngrid)   = j
	            END DO
                 END IF

                 dA(1) = dm
                 dA(2) = dtemp

! -- excute diffusion by Euler approximation
                 A(0,1) = A(m4,1)
                 A(0,2) = A(m4,2)
                 DO i=1,n_SA
                    DO k=m1,m4
                       B(k,i) = A(k,i)+dt*dA(i)*(A(k-1,i)-2.0*A(k,i)+A(k+1,i))/(DX**2.0)
                       A(k,i) = B(k,i)
                    END DO
                 END DO

! -- calculate the terminal velocity v(i) for each droplet (formula only valid for droplets 
! -- <40 mu m)
                 DO i=1,n_drop
                    c    = (2.0*g*rho_w)/(9.0*dyn_vis*rho_a)
                    v(i) = c*(dsize(i)**2.0)
                 END DO

! -- including the terminal velocity for each drop
                 DO i=1,n_drop
                    x(i) = x(i)+v(i)*dt
                    IF (x(i) .GT. BL) THEN
	               x(i) = x(i)-BL
                    END IF
	            IF (x(i) .LT. 0.0) THEN
	               x(i) = x(i)+BL
	            END IF 
                 END DO

! -- let initial jcell_pop(0,j)=0
                 DO j=m1,m4
                    jcell_pop(0,j-m1+1) = 0
                 END DO

! -- calculate the cell position for each drop after sedimetation
                 DO i=1,n_drop
	            j   = INT(x(i)/dx+m1)

                    IF((x(i) .GT. BL) .OR. (x(i) .LT. 0)) THEN
                       WRITE(7,*) 'after sedimetation, x(i) = ',x(i)
                       WRITE(7,*) 'j    = ', j
                       WRITE(7,*) 'i    = ', i
                       WRITE(7,*) 'x(i) = ', x(i)
                       WRITE(7,*) 'dx   = ', dx
                       WRITE(7,*) 'v(i) = ', v(i)
                    END IF
                    jcell_pop(0,j-m1+1)                   = jcell_pop(0,j-m1+1)+1
                    jcell_pop(jcell_pop(0,j-m1+1),j-m1+1) = i
                 END DO

                 t = iteration*dt

	         IF (mod(iteration,1) .EQ. 0) THEN
                    DO i=m1,m4
                       qlwater(i)   = 0.0
                    END DO

	            DO i=1,m_index
	               time_null(i) = 0.0
	            END DO

	            DO i=1,n_drop
! -- find in which grid cell the drop is
                       j_cell = INT(x(i)/dx+m1)
 	               radius = dsize(i)

	               IF (iteration .EQ. 1) THEN
	                  press  = press0
	                  height = height0
                       ELSE
                          press  = press_dgm
                          height = height_dgm
                       END IF

! -- call DGM for each drop
                       CALL dgm(radius,A(j_cell,1),A(j_cell,2),SuS(j_cell),press,w,height,qwater(i),&
                             &  csize(i),i_grid_vol,dt,iteration,aero)

                       dsize(i)        = radius
	               qlwater(j_cell) = qlwater(j_cell)+qwater(i)

	               IF (qlwater(j_cell) .LT. 0.0) THEN
                          WRITE(81,*) j_cell, 'qlwater(j_cell)=', qlwater(j_cell), qwater(i), dsize(i)
	               END IF

! -- write down the time series
                       IF ( (MOD(iteration,iteration3) .EQ. 0) .AND. &
                       &    ( t .GE. start_t                 ) .AND. &
                       &    ( t .LE. end_t                   ) &
                       &  ) THEN
                          r_time(index_drop(i)) = radius*1.e6
                          x_time(index_drop(i)) = x(i)
                          time_null(index_drop(i))   = 1.0

                       END IF
                    END DO
                 END IF

! -- write out time history of droplet radius, supersaturation and droplet locations
                 IF ( (MOD(iteration,iteration3) .EQ. 0) .AND. &
                 &    ( t .GE. start_t                 ) .AND. &
                 &    ( t .LE. end_t                   ) &
                 &  ) THEN
	            DO i = 1, m_index
	               IF (time_null(i) .EQ. 0.0) THEN
	                  r_time(i) = 0.0
	                  x_time(i) = 0.0
	               END IF
	            END DO
                 END IF

                 press_dgm  = press
                 height_dgm = height


! -- DT/Dt due to the temperature change (include the effect of rising velocity)
! -- instead of considering this term in fcnkb.f90
                 DO i=m1,m4
                    cpm     = cp*((1.0+cpv/cp*A(i,1))/(1.0+A(i,1)))
                    temp_in = -(g/cpm)*w*dt
                    A(i,2)  = A(i,2)+temp_in
                 END DO

! -- average qv, temp for different gridcells because of homogeneous mixing
! -- summarize domain qlwater 
                 qv_sum   = 0.0
                 temp_sum = 0.0
	         ql       = 0.0

                 DO i=m1,m4
                    qv_sum   = A(i,1)+qv_sum
                    temp_sum = A(i,2)+temp_sum
	            ql       = ql+qlwater(i)
                 END DO

                 qv_ave   = qv_sum/FLOAT(ngrid)
                 temp_ave = temp_sum/FLOAT(ngrid)

! -- domain average liquid water content
                 ql_ave   = ql/FLOAT(ngrid)

! -- instant mixing case   
                 IF (iinm .EQ. 1) THEN
                    DO i=m1,m4
                       A(i,1) = qv_ave
                       A(i,2) = temp_ave
                    END DO
                 END IF

! -- established periodic B.C.s
                 DO i=1,n_SA
                    DO j=m4+1,m6
                       A(j,i)= A(j-ngrid,i)
                    END DO
                 END DO

                 t = dt*iteration

! -- specified time to record data
!                 IF ( (MOD(iteration,iteration3) .EQ. 0) .AND. &
!                    & ( t .GE. start_t                 ) .AND. &
!                    & ( t .LE. end_t                   ) &
!                    & ) THEN

! -- write out field data
!                    WRITE(63,*) t                            ! time
!                    WRITE(63,*) (A(i,1),i=1,ngrid)
!                    WRITE(64,*) t                            ! time
!                    WRITE(64,*) (A(i,2),i=1,ngrid)
!                 END IF

! -- calculate virtual potential temperature for the environment
                 pt_e   = temp_e*(100000/press)**(Rd/cp)
                 pt_v_e = pt_e*(1.0+0.61*qv_e)

! -- calculate virtual potential tempature for the domain average
                 pt_ave   = temp_ave*(100000/press)**(Rd/cp)
                 pt_v_ave = pt_ave*(1.0+0.61*qv_ave-ql_ave)

! -- calculate buoyancy
                 buoy_ave = g*((pt_v_ave-pt_v_e)/pt_v_e)

! -- calculate LWP of the entire cloud (cloud top is assumed where vertical velocity is negativ)
                 IF (w .GT. 0.0 .AND. LWP_flag) THEN
                    rho_avg   = press/(temp_ave*Rd)*((1.0+qv_ave)/(1.0+qv_ave/eps))
                    LWP       = LWP+(ql_ave*rho_avg*w*dt)
                    count_LWP = count_LWP+1.0
                 END IF
                 
! -- vertical velocity changes according to buoyancy
                 IF (icw .EQ. 0) THEN
                    w = w+buoy_ave*dt
                 END IF

! -- change flag for LWP calculation if vertical velocity is negative
                 IF (w .LT. 0.0) THEN
                    LWP_flag = .FALSE.
                 END IF

! -- using height to find the end of the simulation
                 IF (cloudtopheight .GT. 0.0) THEN
                    IF (Z .GT. cloudtopheight) THEN
                       WRITE(7,*) 'Cloud top height is reached. Leaving the code now.'
                       exit
                    END IF
                 END IF

! -- cloud completely evaporates (ql drops below zero)
                 IF (ql .LT. 0.0) THEN
                    WRITE(7,*) 'Cloud completely evaporated. Leaving the code now.'
                    exit
                 END IF
! -- ending simulation when vertical velocity is smaller than threshold value set in namelist
! -- (and according flag is set)
                 IF (icaw .EQ. 1) THEN
                    IF (w .LT. thsw) THEN
                       WRITE(7,*) 'Vertical velocity is smaller than threshold ', thsw, '.'
                       WRITE(7,*) 'According to flag set in namelist (icaw=1) EMPM stops.'
                       exit
                    END IF
                 END IF

                 IF (mod(iteration,iteration1) .EQ. 0) exit
              END DO

!              IF (mod(iteration,iteration1) .EQ. 0.) THEN
! -- write out history data of droplet radius and locations
              IF (output_format .EQ. 1) THEN
                 WRITE(61,*) t
                 WRITE(61,*) (r_time(i), i=1,m_index)

                 WRITE(62,*) t
                 WRITE(62,*) (x_time(i), i=1,m_index)

                 WRITE(63,*) run_info
                 WRITE(63,*) ngrid, dx 
                 WRITE(63,*) iteration*dt
                 WRITE(63,*) (A(i,1), i=1,ngrid)

                 WRITE(64,*) run_info
                 WRITE(64,*) ngrid, dx 
                 WRITE(64,*) iteration*dt
                 WRITE(64,*) (A(i,2), i=1,ngrid)
              ELSE IF (output_format .EQ. 2) THEN
                 WRITE(61) t
                 WRITE(61) (r_time(i), i=1,m_index)

                 WRITE(62) t
                 WRITE(62) (x_time(i), i=1,m_index)

                 WRITE(63) run_info
                 WRITE(63) ngrid, dx 
                 WRITE(63) iteration*dt
                 WRITE(63) (A(i,1), i=1,ngrid)

                 WRITE(64) run_info
                 WRITE(64) ngrid, dx 
                 WRITE(64) iteration*dt
                 WRITE(64) (A(i,2), i=1,ngrid)
              ELSE IF (output_format .EQ. 3) THEN
                 ip = iteration/iteration1 

                 DO i=m1,m4
                    qv_out(i)   = A(i,1)
                    temp_out(i) = A(i,2)
                 END DO
                 CALL netcdf_write(r_time,x_time,qv_out,temp_out,m_index,ngrid,ip, &
                                   varid_r,varid_x,varid_qv,varid_temp,ncid)
                    
              END IF
! -- calculate supersaturation and its pdf
              DO i=m1,m4
                 es       = EW(A(i,2))*100.0
                 qs       = eps*es/(press-es)
                 SuS(i)   = A(i,1)/qs-1.0
              END DO

              s_max = -0.99
              s_min = +0.1
              CALL max_min(s_max,s_min,SuS,m1,m4,igrid)

              s_diff     = s_max-s_min

              WRITE(88,*) s_max, s_min, s_diff

              WRITE(94,100) iteration,press,qv_ave,temp_ave,ql_ave,qv_ave+ql_ave,w,temp_e,qv_e

              r_max_drop = 0.0
              ip         = iteration/iteration1

              CALL rmax(r_max_drop,dsize,n_drop,index_r_max)

              r_max(ip)       = r_max_drop
              r_max_index(ip) = index_drop(index_r_max)

! -- calculate mean radius, mean r^2 and mean r^3 
              CALL radius_pdf(n_drop,r_pdf_min,r_pdf_max,r_pdf_nbin,dsize,mean_r,mean_r2,mean_r3,r_pdf,r_scalar)
              WRITE(82,*) mean_r
              WRITE(83,*) mean_r2
              WRITE(84,*) mean_r3
              WRITE(85,*) r_pdf

              CALL eff_radius(n_drop,r_pdf_nbin,mean_r2,mean_r3,r_pdf,eff_r)
              WRITE(86,*) t, Z, eff_r

! -- calculate pdf of total water
              DO i=m1,m4
                 qt(i) = A(i,1)+qlwater(i)
              END DO

              CALL qt_prob(m1,m4,qt_pdf_min,qt_pdf_max,qt,qt_pdf,qt_scalar,igrid,qt_pdf_nbin)

              WRITE(71,*) qt_pdf
!                END IF

! -- current number of iteration periods in one realization
              ip = INT(iteration/iteration1)
              WRITE(7,*) 'iperiod = ', ip

              r_sum      = 0.0
              r_std_sum  = 0.0
              n_act_drop = n_drop
              DO i=1,n_drop
                 r_sum = r_sum+dsize(i)*1.0e6
              END DO

              r_mean(ip) = r_sum/FLOAT(n_act_drop)

              DO i=1,n_drop
                 r_std_sum = r_std_sum+(dsize(i)*1.0e6-r_mean(ip))**2.0
              END DO

              r_std(ip)    = ((r_std_sum/FLOAT(n_act_drop))**0.5)
              actndrop(ip) = n_act_drop
              WRITE(81,101) t, r_mean(ip), r_std(ip), actndrop(ip)

! -- calculate supersaturation based on droplets
              s_sum     = 0.0
              s_std_sum = 0.0

              DO i=1,n_drop
                 j_cell = INT(x(i)/dx)+m1
                 s_sum  = s_sum+SuS(j_cell)
              END DO

              s_mean(ip)   = s_sum/FLOAT(n_drop)

              DO i=1,n_drop
                 j_cell    = INT(x(i)/dx)+m1
                 s_std_sum = s_std_sum+(SuS(j_cell)-s_mean(ip))**2.0
              END DO

              s_mean(ip) = s_mean(ip)
              s_std(ip)  = ((s_std_sum/FLOAT(n_drop))**0.5)
              WRITE(87,106) t, s_mean(ip), s_std(ip)
106           format(3(1x,f14.7))

! -------------------------------------------------------------------------------------------------
! begin statistics calculation
! -------------------------------------------------------------------------------------------------

! -- separate scalar A(i,1) to SA1(i); and A(i,2) to SA2(i)
              DO i=m1,m4 
                 SA1(i) = A(i,1)
                 SA2(i) = A(i,2)
              END DO

! -- calculate the probability density function of the scalars SA1 and SA2
! -- The number of bins is set in namelist
              CALL prob(m1,m4,qv_pdf_min,qv_pdf_max,SA1,SA1_PDF,SA1_scalar,scalar_pdf_nbin,ip,mip,igrid)
              CALL prob(m1,m4,temp_pdf_min,temp_pdf_max,SA2,SA2_PDF,SA2_scalar,scalar_pdf_nbin,ip,mip,igrid)

! -- calculate the average of the SA1_PDF and SA2_PDF of differnt realizations
              CALL average(SA1_PDF_ave,SA1_PDF,scalar_pdf_nbin,ip,mreal,mip)
              CALL average(SA2_PDF_ave,SA2_PDF,scalar_pdf_nbin,ip,mreal,mip)

! -- other statistical calculations (moist static energy, liquid water static energy)
              DO i=m1,m4
                 qt(i)   = SA1(i)+qlwater(i)
                 sl(i)   = cp*SA2(i)+g*height-Lv_c*qlwater(i)
                 hm(i)   = cp*SA2(i)+g*height+Lv_c*SA1(i)
              END DO

              DO i=m1,m4
                 scalar(1,i) = qlwater(i)
                 scalar(2,i) = qt(i)
                 scalar(3,i) = sl(i)
                 scalar(4,i) = buoy_ave
                 scalar(5,i) = hm(i)
                 scalar(6,i) = SuS(i)
                 scalar(7,i) = SA1(i)
                 scalar(8,i) = SA2(i)
              END DO

! -- calculate mean and root mean quare of scalars
              CALL mean_rms(m1,m4,scalar,scalar_mean,scalar_rms,ip,mip,igrid)

! -- average mean and rms of different realizations
              CALL average_2(scalar_mean_ave,scalar_mean,ip,mreal,mip)
              CALL average_2(scalar_rms_ave,scalar_rms,ip,mreal,mip)

! -- calculate the probability density function of SuS
              dummy_SuS(:) = scalar(6,:)
              CALL prob(m1,m4,Sus_pdf_min,Sus_pdf_max,dummy_SuS,SuS_pdf,SuS_scalar,scalar_pdf_nbin,ip,mip,igrid)

! -- calculate the average of the SA1_PDF and SA2_PDF of differnt realizations
              CALL average(SuS_pdf_ave,SuS_pdf,scalar_pdf_nbin,ip,mreal,mip)

! -- criteria to end the simulation
! -- top of the environmental profiles is reached
              IF (press .LT. p_array(SIZE(p_array))*100.0) exit

! -- cloud top height set in namelist is reached
              IF (cloudtopheight .GT. 0.0) THEN
                 IF (Z .GT. cloudtopheight) exit
              END IF

! -- cloud completely evaporates (ql drops below zero)
              IF (ql .LT. 0.0) exit

! -- ending simulation when vertical velocity is smaller than threshold value set in namelist
! -- (and according flag is set)
              IF (icaw .EQ. 1) THEN
                 IF (w .LT. thsw) exit
              END IF

! -- ending simulation when top of incloud vertical velocity profile is reached
              IF (ivp .EQ. 1) THEN
                 IF (press .LT. p_cloud(SIZE(p_cloud))*100.0) THEN
                    exit
                 END IF
              END IF

              IF (iteration .GE. iteration2) exit
           END DO

           ipreal = INT(iteration/iteration1)
           WRITE(7,*) '---------------------------------------------------------------------------'
           WRITE(7,*) 'ipreal = ', ipreal
           WRITE(7,*) '---------------------------------------------------------------------------'

           DO j=1,ipreal
! -- scalar_pdf_nbin is the bin number for PDF
! -- Note that the array "j" in PDFave(j,k) is the iteration periods in each realization
              IF (j .EQ. 1) THEN
                 WRITE(92,*) SA1_scalar
                 WRITE(93,*) SA2_scalar
              END IF
              WRITE(92,*) (SA1_PDF_ave(j,k),k=1,scalar_pdf_nbin)
              WRITE(93,*) (SA2_PDF_ave(j,k),k=1,scalar_pdf_nbin)
           END DO
 
           DO ip=1,ipreal
              t = ip*dt*iteration1 
              WRITE(91,107) t, scalar_mean_ave(1,ip), scalar_rms_ave(1,ip), &
                    &       scalar_mean_ave(2,ip), scalar_rms_ave(2,ip),    &
                    &       scalar_mean_ave(3,ip), scalar_rms_ave(3,ip),    &
                    &       scalar_mean_ave(4,ip), scalar_rms_ave(4,ip),    &
                    &       scalar_mean_ave(5,ip), scalar_rms_ave(5,ip),    &
                    &       scalar_mean_ave(6,ip), scalar_rms_ave(6,ip),    &
                    &       scalar_mean_ave(7,ip), scalar_rms_ave(7,ip),    &
                    &       scalar_mean_ave(8,ip), scalar_rms_ave(8,ip),    &
                    &       r_mean(ip), r_std(ip), actndrop(ip), r_max(ip), r_max_index(ip)
           END DO
107        format(19(1x,e14.6),1x,1i10,1x,e14.6,1x,1i10)

           DO i=1,n_drop
             WRITE(72,108) i,index_drop(i),dsize(i)*1.e6
           END DO
108        format(2(1x,i6),1x,e12.4)

           DO j=1,ipreal
! -- scalar_pdf_nbin is the bin number for PDF
! -- Note that the array "j" in PDFave(j,k) is the iteration periods in each realization
              IF (j .EQ. 1) THEN
                 WRITE(95,*) SuS_scalar
              END IF

              WRITE(95,*) (SuS_pdf_ave(j,k),k=1,scalar_pdf_nbin)
           END DO

           WRITE(94,100) iteration,press,qv_ave,temp_ave,ql/FLOAT(ngrid),qv_ave+ql/FLOAT(ngrid),w,temp_e,qv_e

           WRITE(96,*) LWP/(count_LWP*dt)

! -- initialize the iteration in order to perform another realization
           iteration = 0

           IF (output_format .EQ. 1 .OR. output_format .EQ. 2) THEN
              CLOSE(unit=61)
              CLOSE(unit=62)
           ELSE IF (output_format .EQ. 3) THEN
              CALL check( NF90_CLOSE(ncid) ) 
           END IF
           CLOSE(unit=63)
           CLOSE(unit=64)

           CLOSE(unit=71)
           CLOSE(unit=72)
           CLOSE(unit=73)

           CLOSE(unit=81)
           CLOSE(unit=82)
           CLOSE(unit=83)
           CLOSE(unit=84)
           CLOSE(unit=85)
           CLOSE(unit=86)
           CLOSE(unit=87)
           CLOSE(unit=88)

           CLOSE(unit=91)
           CLOSE(unit=92)
           CLOSE(unit=93)
           CLOSE(unit=94)
           CLOSE(unit=95)
           CLOSE(unit=96)

! -- continue with another realization
        END DO

        WRITE(7,*) ''
        WRITE(7,*) 'EMPM was executed correctly' 
        WRITE(7,*) '---------------------------------------------------------------------------'
        WRITE(7,*) '---------------------------------------------------------------------------'

        time_end = time()
        WRITE(7,*) 'elpased time in min. ', (time_end - time_start)/60.

        CLOSE(7)

        END
