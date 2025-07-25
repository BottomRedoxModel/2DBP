! This file is part of Bottom RedOx Model (BROM, v.1.1).
! BROM is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the BROM distribution.
!-----------------------------------------------------------------------
! Original author(s): Evgeniy Yakushev, Shamil Yakubov,
!                     Elizaveta Protsenko, Phil Wallhead, Anfisa Berezina, 
!                     Beatriz Arellano-Nava
!-----------------------------------------------------------------------


    module io_netcdf

      use input_mod
      use types_mod
  
      use netcdf
      use fabm_omp, only: type_fabm_model => type_fabm_omp_model
      use fabm_types, only: attribute_length, rk
  
  
      implicit none
      !all is private
      private
      !public functions
      public init_netcdf, input_netcdf_2, save_netcdf, close_netcdf
      !netCDF file id
      integer               :: nc_id
      integer, allocatable  :: parameter_id(:)
      integer, allocatable  :: parameter_fick_id(:)
      integer, allocatable  :: parameter_sink_id(:)
      integer, allocatable  :: parameter_id_diag(:)
  
      integer               :: i_id, z_id, z2_id, time_id, swradWm2_id, hice_id
      integer               :: pH_id, T_id, S_id, Kz_id, Kz_sol_id, Kz_par_id, w_sol_id, w_par_id, u_x_id, gas_air_sea_id
      integer               :: pCO2_id, Om_Ca_id, Om_Ar_id
  
      logical               :: first
  
  
  
      contains
  !=======================================================================================================================
      subroutine input_netcdf_2(z_w, dz_w, hz_w, t_w, s_w, kz_w, use_swradWm2, &
          hice, swradWm2, aice, use_hice, gargett_a0, gargett_q, use_gargett, &
          year, i_max, steps_in_yr, k_wat_bbl, u_x_w)
  
  ! inputs hydrophysical data  (time,depth,t,s,Kz,u,v,ice,light) from netCDF files
  !------------------------------------------------------------------------------------------------------------------------
  
      use io_ascii, only: get_brom_name, get_brom_par
  
  !Input variables
      integer, intent(in)                         :: use_swradWm2, use_hice, year, i_max, steps_in_yr, use_Gargett
      real(rk), intent(in)                        :: gargett_a0, gargett_q
  !Output variables
      real(rk), allocatable, dimension(:), intent(out)        :: z_w, dz_w         ! Layer midpoint depths and spacing between them
      real(rk), allocatable, dimension(:), intent(out)        :: hz_w              ! Layer thicknesses
      real(rk), allocatable, dimension(:,:,:), intent(out)    :: t_w, s_w, kz_w    ! Temperature, salinity, and vertical diffusivity
      real(rk), allocatable, dimension(:,:,:), intent(out)    :: u_x_w             ! Horizontal advection [m/s]
      real(rk), allocatable, dimension(:), intent(out)        :: hice              ! Ice thickness [m]
      real(rk), allocatable, dimension(:), intent(out)        :: aice              ! Ice square share  [-]
      real(rk), allocatable, dimension(:), intent(out)        :: swradWm2          ! 24-hr average surface downwelling shortwave irradiance in air [W/m2]
  !Input/output variables
      integer, intent(out)                        :: k_wat_bbl
  !Local variables
      type(type_input):: input
      class(variable), allocatable:: var
      integer:: i_time, i_depth  !sizes of input arrays
      real(rk), allocatable, dimension(:)         :: time_temp, z_temp, hice_temp, swradWm2_temp, aice_temp, dens !, z_temp2z_w2, z_w_error, Eair_temp,
      real(rk), allocatable, dimension(:,:)       :: t_temp, s_temp, kz_temp, u_temp, v_temp
      character(len=64) :: ncinfile_name !, ncint_name, ncins_name, ncinkz_name, ncinEair_name, ncinhice_name, &
                         ! ncinz_name, ncinz2_name, ncintime_name, ncinlat_name, ncinlon_name, ncinhmix_rate_name
      integer           :: nc_year0, nc_set_k_wat_bbl !, nc_staggered_grid, nc_bottom_to_top, nc_z_increasing_upward,
      integer           :: year0, yeari, nyrdays, ni, nyrdaysm1, i, j, ip, istart, iend, iday !, inext,ilast, ncid
      integer           :: istep! AB, for hourly input
      real(rk)          :: time0, timei, time00
  ! Input and output file names
      ncinfile_name = get_brom_name("ncinfile_name")
  !Other NetCDF parameters
      nc_set_k_wat_bbl = get_brom_par("nc_set_k_wat_bbl")
      nc_year0 = get_brom_par("nc_year0")
  
  !------------------------------------------------------------------------------------------------------------------------
  !---------input data from input file
      input = type_input(ncinfile_name)
      i_time = input%get_1st_dim_length('oc_time')
      i_depth = input%get_1st_dim_length('depth')
  
      call input%get_var('depth', var)
      select type(var)
      class is(alone_variable)
        write(*,*)
      class is(variable_1d)
        allocate(z_temp(i_depth))
        z_temp = abs(var%value)
      class is(variable_2d)
        write(*,*)
      end select
      deallocate(var)
  
      call input%get_var('oc_time', var)
      select type(var)
      class is(alone_variable)
        write(*,*)
      class is(variable_1d)
        allocate(time_temp(i_time))
        time_temp = var%value
      class is(variable_2d)
        write(*,*)
      end select
      deallocate(var)
  
      if (use_hice.eq.1) then
      call input%get_var('hice', var)
      select type(var)
      class is(alone_variable)
        write(*,*)
      class is(variable_1d)
        allocate(hice_temp(i_time))
        hice_temp = var%value
      class is(variable_2d)
        write(*,*)
      end select
      deallocate(var)
      end if
  
      if (use_hice.eq.1) then
      call input%get_var('aice', var)
      select type(var)
      class is(alone_variable)
        write(*,*)
      class is(variable_1d)
        allocate(aice_temp(i_time))
        aice_temp = var%value
      class is(variable_2d)
        write(*,*)
      end select
      deallocate(var)
      end if
  
      if (use_swradWm2.eq.1) then
      call input%get_var('swradWm2', var)
      select type(var)
      class is(alone_variable)
        write(*,*)
      class is(variable_1d)
        allocate(swradWm2_temp(i_time))
        swradWm2_temp = var%value
      class is(variable_2d)
        write(*,*)
      end select
      deallocate(var)
      end if
  
      call input%get_var('salt', var)
      select type(var)
      class is(alone_variable)
        write(*,*)
      class is(variable_1d)
        write(*,*)
      class is(variable_2d)
        allocate(s_temp(i_time,i_depth))
        s_temp = var%value
      end select
      deallocate(var)
  
      call input%get_var('temp', var)
      select type(var)
      class is(alone_variable)
        write(*,*)
      class is(variable_1d)
        write(*,*)
      class is(variable_2d)
        allocate(t_temp(i_time,i_depth))
        t_temp = var%value
      end select
      deallocate(var)
  
      call input%get_var('Kz', var)
      select type(var)
      class is(alone_variable)
        write(*,*)
      class is(variable_1d)
        write(*,*)
      class is(variable_2d)
        allocate(kz_temp(i_time,i_depth))
        kz_temp = var%value
      end select
      deallocate(var)
  
      call input%get_var('u', var)
      select type(var)
      class is(alone_variable)
        write(*,*)
      class is(variable_1d)
        write(*,*)
      class is(variable_2d)
        allocate(u_temp(i_time,i_depth))
        u_temp = var%value
      end select
      deallocate(var)
  
      call input%get_var('v', var)
      select type(var)
      class is(alone_variable)
        write(*,*)
      class is(variable_1d)
        write(*,*)
      class is(variable_2d)
        allocate(v_temp(i_time,i_depth))
        v_temp = var%value
      end select
      deallocate(var)
  
      call input%delete_list()
  !-------------------------------------------------------------------------------------------
  
      !Set the no. grid points in the water column k_wat_bbl using the netCDF input
      if (nc_set_k_wat_bbl.eq.1) then
          k_wat_bbl = i_depth
          write(*,*) "k_wat_bbl set from netCDF input to ", k_wat_bbl
      else
          write(*,*) "k_wat_bbl set from brom.yaml to ", k_wat_bbl
      end if
      !Allocate permanent variables
      allocate(t_w(i_max,k_wat_bbl,steps_in_yr))
      allocate(s_w(i_max,k_wat_bbl,steps_in_yr))
      allocate(kz_w(i_max,k_wat_bbl+1,steps_in_yr))
      allocate(z_w(k_wat_bbl))
      allocate(dens(k_wat_bbl))
      allocate(dz_w(k_wat_bbl))
      allocate(hz_w(k_wat_bbl))
      allocate(hice(steps_in_yr))
      allocate(aice(steps_in_yr))
      allocate(swradWm2(steps_in_yr))
      allocate(u_x_w(i_max,k_wat_bbl,steps_in_yr))
  !-------------------------------------------------------------------------------------------
      !!Establish the initial index istart using the netcdf time variable and the chosen year
      !First calculate the earliest year "year0" and corresponding time "time00" [s] at start of this year
      time0 = time_temp(1) !Assume that netcdf time is time in seconds since nc_year0-01-01 00:00:00
      year0 = -1
      yeari = nc_year0
      timei = 0
      do while (year0.eq.-1)
          nyrdays = (365 + merge(1,0,(mod(yeari,4).eq.0)))
          if ((time0-timei).ge.0.0_rk.and.(time0-timei).lt.(nyrdays*86400)) then !Found initial year in file (originally nyrdays)
              year0 = yeari
              time00 = timei
          elseif (time0.gt.0.0_rk) then !Advance candidate year0 = yeari forwards by one year
              yeari = yeari + 1
              timei = timei + nyrdays*86400
          elseif (time0.lt.0.0_rk) then !Advance candidate year0 = yeari backwards by one year
              nyrdaysm1 = (365 + merge(1,0,(mod(yeari-1,4).eq.0)))
              yeari = yeari - 1
              timei = timei - nyrdaysm1*86400
          end if
      end do
      !Find the starting and finishing time indices (istart, iend) for the selected year "year"
      istart = 1
      iend = 1
      yeari = year0
      timei = time00
      do i=1,i_time
          nyrdays = (365 + merge(1,0,(mod(yeari,4).eq.0)))
          if ((time_temp(i)-timei).ge.(nyrdays*86400)) then !Advance the year and time counters
              yeari = yeari + 1
              timei = timei + nyrdays*86400
          end if
          if (yeari.eq.year.and.istart.eq.-1) istart = i    !istart = first index where (yeari==year)
      end do
      iend=istart+steps_in_yr-1
      ni = iend-istart+1
        if (ni.lt.steps_in_yr) then
          write(*,*) "Could not find steps_in_yr time inputs starting from 1st day of selected year (ni = ", ni, ", stopping"
          stop
        end if
  
      !Check - these results can be validated by importing ocean_time into Matlab and converting to human dates using http://www.epochconverter.com/
      write(*,*) "First year in netCDF input: year0 = ", year0
      write(*,*) "Reference year for netCDF time in seconds: nc_year0 = ", nc_year0
      write(*,*) "NetCDF time in seconds at beginning of first recorded year: time00 = ", time00
      write(*,*) "Target year for BROM model input: year = ", year
      write(*,*) "Calculated starting index: istart = ", istart
      write(*,*) "Calculated finishing index: iend = ", iend
      write(*,*) "Number of input time points: ni = ", ni
      write(*,*) "Initial input time in netCDF units for chosen year: time_temp(istart) = ", time_temp(istart)
      write(*,*) "Final input time in netCDF units for chosen year: time_temp(iend) = ", time_temp(iend)
  
  !------fill in model' arrays:
  
      !depths of layer midpoints (z_w)
      z_w(:) = z_temp(:)
      !Spacing between layer midpoints (dz_w)
      dz_w(1:k_wat_bbl-1) = z_w(2:k_wat_bbl) - z_w(1:k_wat_bbl-1)
      dz_w(k_wat_bbl) = 0.0_rk
      !Layer thicknesses (hz_w)
      hz_w(1) = z_w(1) + 0.5_rk*dz_w(1)
      do j=2,k_wat_bbl
          hz_w(j) = 0.5_rk*(dz_w(j-1)+dz_w(j))
      end do

      !Set the water temperature (t_w), salinity (s_w), vertical diffusivity (kz_w),
      !and (if required) the surface irradiance (Eair) and ice thickness (hice)
    !  do istep=1,steps_in_yr-1 !Loop over hours
      do istep=0,steps_in_yr-1 !Loop over hours
          t_w(i_max,1:k_wat_bbl,istep+1) = t_temp(istart+istep,1:k_wat_bbl)
          s_w(i_max,1:k_wat_bbl,istep+1) = s_temp(istart+istep,1:k_wat_bbl)
          Kz_w(i_max,1:k_wat_bbl,istep+1) = kz_temp(istart+istep,1:k_wat_bbl)
          u_x_w(i_max,1:k_wat_bbl,istep+1) = u_temp(istart+istep,1:k_wat_bbl)
  !        u_x_w(i_max,1:k_wat_bbl,istep) = v_temp(istart+istep,1:k_wat_bbl)
          if (use_hice.eq.1) hice(istep+1) = hice_temp(istart+istep)
          if (use_hice.eq.1) aice(istep+1) = aice_temp(istart+istep)
          if (use_swradWm2.eq.1) swradWm2(istep+1) = swradWm2_temp(istart+istep)
  
  !        t_w(i_max,1:k_wat_bbl,istep) = t_temp(istart+istep+1,1:k_wat_bbl)
  !        s_w(i_max,1:k_wat_bbl,istep) = s_temp(istart+istep+1,1:k_wat_bbl)
  !        Kz_w(i_max,1:k_wat_bbl,istep) = kz_temp(istart+istep+1,1:k_wat_bbl)
  !        u_x_w(i_max,1:k_wat_bbl,istep) = u_temp(istart+istep+1,1:k_wat_bbl)
  !!        u_x_w(i_max,1:k_wat_bbl,istep) = v_temp(istart+istep+1,1:k_wat_bbl)
  !        if (use_hice.eq.1) hice(istep) = hice_temp(istart+istep+1)
  !        if (use_hice.eq.1) aice(istep) = aice_temp(istart+istep+1)
  !        if (use_swradWm2.eq.1) swradWm2(istep) = swradWm2_temp(istart+istep)
  
      enddo
  
      if (use_gargett.eq.1) then
      ! Calculate Kz using Gargett formula if needed:
        do istep=1,steps_in_yr
            do j=1, k_wat_bbl
                call svan(s_w(i_max,j,istep),t_w(i_max,j,istep),z_w(j),dens(j)) !calculate density as f(p,t,s)
            end do
            do j=1, k_wat_bbl-1
                Kz_w(i_max,j,istep)=gargett_a0 &
                    /((9.81/(1000.+(dens(j)+dens(j+1))/2.)&
                    *(abs(dens(j+1)-dens(j))/dz_w(j)) &
                    )**gargett_q)
            end do
                Kz_w(i_max,k_wat_bbl,istep)=Kz_w(i_max,k_wat_bbl-1,istep)
        end do
      endif
      ! fill all the horizontal columns
      do i = 1, i_max
          t_w(i,:,:)  =  t_w(i_max,:,:)
          s_w(i,:,:)  =  s_w(i_max,:,:)
          kz_w(i,:,:) = min(0.10,kz_w(i_max,:,:)) !max(0.00000001,min(0.10,kz_w(i_max,:,:)))
          u_x_w(i,:,:)= u_x_w(i_max,:,:)
      enddo
  
      ! deallocate(t_temp)
      ! deallocate(s_temp)
      ! deallocate(kz_temp)
      ! deallocate(u_temp)
      ! deallocate(v_temp)
      ! if (use_hice.eq.1) deallocate(hice_temp)
      ! if (use_hice.eq.1) deallocate(aice_temp)
      ! if (use_swradWm2.eq.1) deallocate(swradWm2_temp)
      ! deallocate(z_temp)
      ! deallocate(time_temp)
      end subroutine
  
  
  
  
  
  
  
  !=======================================================================================================================
      subroutine init_netcdf(fn, i_max, k_max, model, use_swradWm2, use_hice, year)
  
      !Input variables
      character(len=*), intent(in)     :: fn
      integer, intent(in)              :: k_max, i_max, use_swradWm2, use_hice, year
      class (type_fabm_model), pointer :: model
  
      !Local variables
      integer                          :: i_dim_id, z_dim_id, z2_dim_id, time_dim_id
      integer                          :: ip, iret, ilast
      integer, parameter               :: time_len = NF90_UNLIMITED
      character(len=4)                 :: yearstr
      integer                          :: dim1d
      integer                          :: dim_ids(3), dim_ids2(3), dim_ids0(2)
  
    write(*,*) "k_max = ", k_max
  
      first = .true.
      print *, 'NetCDF version: ', trim(nf90_inq_libvers())
      nc_id = -1
      call check_err(nf90_create(fn, NF90_CLOBBER, nc_id))
  
      !Define the dimensions
      call check_err(nf90_def_dim(nc_id, "i", i_max, i_dim_id))
      call check_err(nf90_def_dim(nc_id, "z", k_max, z_dim_id))
      call check_err(nf90_def_dim(nc_id, "z2", k_max+1, z2_dim_id))
      call check_err(nf90_def_dim(nc_id, "time", time_len, time_dim_id))
  
      !Define coordinates
      dim1d = z_dim_id
      call check_err(nf90_def_var(nc_id, "z", NF90_REAL, dim1d, z_id))
      call check_err(nf90_put_att(nc_id, z_id, "positive", "down"))
      call check_err(nf90_put_att(nc_id, z_id, "long_name", "depth at layer midpoints"))
      call check_err(nf90_put_att(nc_id, z_id, "units", "metres"))
      call check_err(nf90_put_att(nc_id, z_id, "axis", "Z"))
      dim1d = z2_dim_id
      call check_err(nf90_def_var(nc_id, "z2", NF90_REAL, dim1d, z2_id))
      call check_err(nf90_put_att(nc_id, z2_id, "positive", "down"))
      call check_err(nf90_put_att(nc_id, z2_id, "long_name", "depth at layer interfaces"))
      call check_err(nf90_put_att(nc_id, z2_id, "units", "metres"))
      call check_err(nf90_put_att(nc_id, z2_id, "axis", "Z"))
      dim1d = i_dim_id
      call check_err(nf90_def_var(nc_id, "i", NF90_REAL, dim1d, i_id))
      call check_err(nf90_put_att(nc_id, i_id, "positive", "right"))
      call check_err(nf90_put_att(nc_id, i_id, "long_name", "distance at horizontal axis"))
      call check_err(nf90_put_att(nc_id, i_id, "units", "metres"))
      call check_err(nf90_put_att(nc_id, i_id, "axis", "X"))
      dim1d = time_dim_id
      call check_err(nf90_def_var(nc_id, "time", NF90_REAL, dim1d, time_id))
      write(yearstr,'(i4)') year
      call check_err(nf90_put_att(nc_id, time_id, "long_name", "time"))
      call check_err(nf90_put_att(nc_id, time_id, "units", "days since "//trim(yearstr)//"-01-01 00:00:00"))
      call check_err(nf90_put_att(nc_id, time_id, "axis", "T"))
  
      write(*,*) "Init params"
      write(*,*) "i_dim_id: ", i_id
      write(*,*) "z_dim_id: ", z_id
      write(*,*) "z2_dim_id: ", z2_dim_id
  
      !Define state variables
      dim_ids = (/i_dim_id, z_dim_id, time_dim_id/)
      dim_ids2 = (/i_dim_id, z2_dim_id, time_dim_id/)
      dim_ids0 = (/i_dim_id, time_dim_id/)
      allocate(parameter_id(size(model%interior_state_variables)))
      allocate(parameter_fick_id(size(model%interior_state_variables)))
      allocate(parameter_sink_id(size(model%interior_state_variables)))
      do ip=1,size(model%interior_state_variables)
          ilast = index(model%interior_state_variables(ip)%path,'/',.true.)
          call check_err(nf90_def_var(nc_id, model%interior_state_variables(ip)%path(ilast+1:), NF90_REAL, dim_ids, parameter_id(ip)))  ! was ilast+1:
          call check_err(nf90_def_var(nc_id, 'fick:'//model%interior_state_variables(ip)%path(ilast+1:), NF90_REAL, dim_ids2, parameter_fick_id(ip)))
          call check_err(nf90_def_var(nc_id, 'sink:'//model%interior_state_variables(ip)%path(ilast+1:), NF90_REAL, dim_ids2, parameter_sink_id(ip)))
          call check_err(set_attributes(ncid=nc_id, id=parameter_id(ip), units=model%interior_state_variables(ip)%units, &
              long_name=model%interior_state_variables(ip)%long_name, missing_value=model%interior_state_variables(ip)%missing_value))
          call check_err(set_attributes(ncid=nc_id, id=parameter_fick_id(ip), units='mmol/m^2/day', &
              long_name='fick:'//model%interior_state_variables(ip)%long_name,missing_value=model%interior_state_variables(ip)%missing_value))
          call check_err(set_attributes(ncid=nc_id, id=parameter_sink_id(ip), units='mmol/m^2/day', &
              long_name='sink:'//model%interior_state_variables(ip)%long_name,missing_value=model%interior_state_variables(ip)%missing_value))
          call check_err(nf90_put_att(nc_id, parameter_fick_id(ip), "positive", "down"))
          call check_err(nf90_put_att(nc_id, parameter_sink_id(ip), "positive", "down"))
      end do
  
      !Define diagnostic variables
      allocate(parameter_id_diag(size(model%interior_diagnostic_variables)))
      !do ip=1,size(model%interior_diagnostic_variables)
      !    write(*,*) model%interior_diagnostic_variables(ip)%name
      !enddo
      do ip=1,size(model%interior_diagnostic_variables)
          if (model%interior_diagnostic_variables(ip)%save) then
              ilast = index(model%interior_diagnostic_variables(ip)%path,'/',.true.)
              call check_err(nf90_def_var(nc_id, model%interior_diagnostic_variables(ip)%path(ilast+1:), NF90_REAL, dim_ids, parameter_id_diag(ip)))
              call check_err(set_attributes(ncid=nc_id, id=parameter_id_diag(ip), units=model%interior_diagnostic_variables(ip)%units, &
                  long_name=model%interior_diagnostic_variables(ip)%long_name,missing_value=model%interior_diagnostic_variables(ip)%missing_value))
          end if
      end do
  
      !Define forcing variables used in the run
      call check_err(nf90_def_var(nc_id, "T", NF90_REAL, dim_ids, T_id))
      call check_err(nf90_put_att(nc_id, T_id, "long_name", "temperature"))
      call check_err(nf90_put_att(nc_id, T_id, "units", "degC"))
      call check_err(nf90_def_var(nc_id, "S", NF90_REAL, dim_ids, S_id))
      call check_err(nf90_put_att(nc_id, S_id, "long_name", "salinity"))
      call check_err(nf90_def_var(nc_id, "Kz", NF90_REAL, dim_ids2, Kz_id))
      call check_err(nf90_put_att(nc_id, Kz_id, "long_name", "vertical eddy diffusivity"))
      call check_err(nf90_put_att(nc_id, Kz_id, "units", "m2/s"))
      call check_err(nf90_def_var(nc_id, "Ux", NF90_REAL, dim_ids, u_x_id))
      call check_err(nf90_put_att(nc_id, u_x_id, "long_name", "horizontal advection"))
      call check_err(nf90_put_att(nc_id, u_x_id, "units", "m/s"))
      call check_err(nf90_def_var(nc_id, "Kz_sol", NF90_REAL, dim_ids2, Kz_sol_id))
      call check_err(nf90_put_att(nc_id, Kz_sol_id, "long_name", "total vertical diffusivity of a solute"))
      call check_err(nf90_put_att(nc_id, Kz_sol_id, "units", "m2/s"))
      call check_err(nf90_def_var(nc_id, "Kz_par", NF90_REAL, dim_ids2, Kz_par_id))
      call check_err(nf90_put_att(nc_id, Kz_par_id, "long_name", "total vertical diffusivity of a particulate"))
      call check_err(nf90_put_att(nc_id, Kz_par_id, "units", "m2/s"))
      call check_err(nf90_def_var(nc_id, "w_sol", NF90_REAL, dim_ids2, w_sol_id))
      call check_err(nf90_put_att(nc_id, w_sol_id, "long_name", "total advective velocity of a solute"))
      call check_err(nf90_put_att(nc_id, w_sol_id, "units", "m/s"))
      call check_err(nf90_def_var(nc_id, "w_par", NF90_REAL, dim_ids2, w_par_id))
      call check_err(nf90_put_att(nc_id, w_par_id, "long_name", "total advective velocity of a particulate"))
      call check_err(nf90_put_att(nc_id, w_par_id, "units", "m/s"))
      call check_err(nf90_def_var(nc_id, "gas_air_sea", NF90_REAL, dim_ids, gas_air_sea_id))  ! DIC air-sea flux
      call check_err(nf90_put_att(nc_id, gas_air_sea_id, "long_name", "gas_air_sea"))
      call check_err(nf90_put_att(nc_id, gas_air_sea_id, "units", "mmol/m2/d"))
      if (use_swradWm2.eq.1) then
          call check_err(nf90_def_var(nc_id, "swradWm2", NF90_REAL, time_dim_id, swradWm2_id))
          call check_err(nf90_put_att(nc_id, swradWm2_id, "units", "W/m2"))
      end if
      if (use_hice.eq.1) then
          call check_err(nf90_def_var(nc_id, "hice", NF90_REAL, time_dim_id, hice_id))
          call check_err(nf90_put_att(nc_id, hice_id, "units", "m"))
      end if
  
      call check_err(nf90_enddef(nc_id))
  
      end subroutine init_netcdf
  !=======================================================================================================================
  
  
  
  
  
  
  
  !=======================================================================================================================
      subroutine save_netcdf(i_max, k_max, julianday, cc, t, s, kz, kzti, wti, &
          model, z, hz, swradWm2, use_swradWm2, hice, use_hice, fick_per_day, sink_per_day, &
          ip_sol, ip_par, x, u_x, i_day, id, idt, output_step, input_step, gas_air_sea, multiyears_physics ) !i_sec_pr) ! i_day here = i_day + 1
  
      !Input variables
      integer, intent(in)                    :: i_max, k_max, julianday, use_swradWm2, input_step
      integer, intent(in)                    :: use_hice, ip_sol, ip_par, i_day, id, idt, output_step, multiyears_physics
      real(rk), dimension(:,:,:), intent(in) :: cc, t, s, kz, kzti, wti, fick_per_day, sink_per_day, u_x, gas_air_sea
      class (type_fabm_model), pointer :: model
      real(rk), dimension(:), intent(in)     :: z, hz, swradWm2, hice, x
  
      !Local variables
      integer, dimension(1)                  :: start_z, count_z, start_z2, count_z2, start_time, count_time, start_x, count_x
      integer, dimension(3)                  :: start_cc, count_cc, start_flux, count_flux
      real(rk)                               :: temp_matrix(i_max,k_max), dum(1), z2(k_max+1), day_part
      !Note: The input arguments to nf90_put_var MUST be vectors, even if the length is 1
      !      Removing the dimension(1) or (1) from dum above triggers a spurious error "not finding nf90_put_var"
      integer                                :: ip, i, i_sec, istep_out, i_day_share
        !integer,parameter         :: timestepkind = selected_int_kind(12)   !how to make int(8)
        !integer(kind=timestepkind):: i_sec  !declared as int(8) it will not work with netcdf functions 
  !    day_part=real(id)/real(idt)
  !    i_sec=((i_day-1)*86400 + int((86400*(id/100))/(idt/100)))/output_step ! as in Horten
  !    istep_out = int((julianday)*86400/input_step)
  
      day_part=real(id)/real(idt)
      i_day_share = 86400/output_step  !a multipier allowing to decrease max int number in i_sec
!      i_sec=int(((i_day-1)*86400 + int(86400*id/idt))/output_step) ! time count for saving  (in array numbers)
      i_sec=int(((i_day-1)*i_day_share + int(i_day_share*id/idt))) ! time count for saving  (in array numbers)      
      istep_out = max(1, int(((julianday-1)*86400 + int(86400*id/idt))/input_step)) ! time count to select data from arrays, i.e. temp, salt
  
   !Define nf90_put_var arguments "start" and "count" for z, z2, time, (cc,t,s) and (fick,kz)
      start_z = 1
      count_z = k_max
      start_z2 = 1
      count_z2 = k_max+1
      start_time = i_sec
      count_time = 1
      start_x = 1
      count_x = i_max !number of columns
      start_cc(1:3) = (/1, 1, i_sec/) ! i_sec
      count_cc(1:3) = (/i_max, k_max, 1/)
      start_flux(1:3) = (/1, 1, i_sec/) ! i_sec
      count_flux(1:3) = (/i_max, k_max+1, 1/)
      !At first call only, output depth variable mz = -1*z
      if (first) then
          call check_err(nf90_put_var(nc_id, z_id, z, start_z, count_z))
          z2(1:k_max) = z(1:k_max) - 0.5_rk*hz(1:k_max)
          z2(k_max+1) = z(k_max) + 0.5_rk*hz(k_max)
          call check_err(nf90_put_var(nc_id, z2_id, z2, start_z2, count_z2))
          call check_err(nf90_put_var(nc_id, i_id, x, start_x, count_x))
          !Note: nc_id, z_id and z2_id are available to the entire module and defined in init_netcdf
          first = .false.
      end if
      dum(1) = (real(i_day)+day_part)
      !For all calls output cc, fick, diagnostics and forcings (t,s,kz)
      if (nc_id.ne.-1) then
          call check_err(nf90_put_var(nc_id, time_id, dum, start_time, count_time))
          do ip=1,size(model%interior_state_variables)
              call check_err(nf90_put_var(nc_id, parameter_id(ip), cc(:,:,ip), start_cc, count_cc))
              call check_err(nf90_put_var(nc_id, gas_air_sea_id, gas_air_sea(:,:,ip), start_cc, count_cc))
              call check_err(nf90_put_var(nc_id, parameter_fick_id(ip), fick_per_day(:,:,ip), start_flux, count_flux))
              call check_err(nf90_put_var(nc_id, parameter_sink_id(ip), sink_per_day(:,:,ip), start_flux, count_flux))
          end do
          do ip=1,size(model%interior_diagnostic_variables)
              if (model%interior_diagnostic_variables(ip)%save) then
                  temp_matrix = model%get_interior_diagnostic_data(ip)
                  if (maxval(abs(temp_matrix)).lt.1.0E37) then
                      !This is to avoid "NetCDF numeric conversion error" for some NetCDF configurations
                      call check_err(nf90_put_var(nc_id, parameter_id_diag(ip), temp_matrix, start_cc, count_cc))
                  end if
              end if
          end do
          if(multiyears_physics.gt.0) then
            call check_err(nf90_put_var(nc_id, T_id, t(:,:,i_sec), start_cc, count_cc))
            call check_err(nf90_put_var(nc_id, S_id, s(:,:,i_sec), start_cc, count_cc))
            call check_err(nf90_put_var(nc_id, Kz_id, kz(:,:,i_sec), start_flux, count_flux))
          else
            call check_err(nf90_put_var(nc_id, T_id, t(:,:,istep_out), start_cc, count_cc))
            call check_err(nf90_put_var(nc_id, S_id, s(:,:,istep_out), start_cc, count_cc))
            call check_err(nf90_put_var(nc_id, Kz_id, kz(:,:,istep_out), start_flux, count_flux))
          endif

          if (i_max.gt.1) call check_err(nf90_put_var(nc_id, u_x_id, u_x(:,:,istep_out), start_cc, count_cc))
  
          call check_err(nf90_put_var(nc_id, Kz_sol_id, kzti(:,:,ip_sol), start_flux, count_flux))
          call check_err(nf90_put_var(nc_id, Kz_par_id, kzti(:,:,ip_par), start_flux, count_flux))
          call check_err(nf90_put_var(nc_id, w_sol_id, wti(:,:,ip_sol), start_flux, count_flux))
          call check_err(nf90_put_var(nc_id, w_par_id, wti(:,:,ip_par), start_flux, count_flux))
  
          if (use_swradWm2.eq.1) then
              dum(1) = swradWm2(i_day) !julianday
              call check_err(nf90_put_var(nc_id, swradWm2_id, dum, start_time, count_time))
          end if
  
          if (use_hice.eq.1) then
              dum(1) = hice(i_day) !julianday
              call check_err(nf90_put_var(nc_id, hice_id, dum, start_time, count_time))
          end if
          call check_err(nf90_sync(nc_id))
      end if
  
      end subroutine save_netcdf
  !=======================================================================================================================
  
  
  
  
  
  
  
  
  !=======================================================================================================================
      subroutine close_netcdf()
  
      implicit none
      if (nc_id.ne.-1) then
          call check_err(nf90_close(nc_id))
          deallocate(parameter_id)
          deallocate(parameter_fick_id)
          deallocate(parameter_sink_id)
          deallocate(parameter_id_diag)
          write (*,'(a)') "finished"
      end if
      nc_id = -1
  
      end subroutine close_netcdf
  !=======================================================================================================================
  
  
  
  
  
  !=======================================================================================================================
      integer function set_attributes(ncid,id,                         &
                                      units,long_name,                 &
                                      valid_min,valid_max,valid_range, &
                                      scale_factor,add_offset,         &
                                      FillValue,missing_value,         &
                                      C_format,FORTRAN_format)
      !
      ! !DESCRIPTION:
      !  This routine is used to set a number of attributes for
      !  variables. The routine makes heavy use of the {\tt optional} keyword.
      !  The list of recognized keywords is very easy to extend. We have
      !  included a sub-set of the COARDS conventions.
      !
      ! !USES:
      !  IMPLICIT NONE
      !
      ! !INPUT PARAMETERS:
      integer, intent(in)                     :: ncid,id
      character(len=*), optional              :: units,long_name
      real, optional                          :: valid_min,valid_max
      real, optional                          :: valid_range(2)
      real, optional                          :: scale_factor,add_offset
      double precision, optional              :: FillValue,missing_value
      character(len=*), optional              :: C_format,FORTRAN_format
      !
      ! !REVISION HISTORY:
      !  Original author(s): Karsten Bolding & Hans Burchard
      !
      ! !LOCAL VARIABLES:
      integer                                 :: iret
      real                                    :: vals(2)
      !
      !
      !-----------------------------------------------------------------------
      !
      if (present(units)) then
          iret = nf90_put_att(ncid,id,'units',trim(units))
      end if
  
      if (present(long_name)) then
          iret = nf90_put_att(ncid,id,'long_name',trim(long_name))
      end if
  
      if (present(C_format)) then
          iret = nf90_put_att(ncid,id,'C_format',trim(C_format))
      end if
  
      if (present(FORTRAN_format)) then
          iret = nf90_put_att(ncid,id,'FORTRAN_format',trim(FORTRAN_format))
      end if
  
      if (present(valid_min)) then
          vals(1) = valid_min
          iret = nf90_put_att(ncid,id,'valid_min',vals(1:1))
      end if
  
      if (present(valid_max)) then
          vals(1) = valid_max
          iret = nf90_put_att(ncid,id,'valid_max',vals(1:1))
      end if
  
      if (present(valid_range)) then
          vals(1) = valid_range(1)
          vals(2) = valid_range(2)
          iret = nf90_put_att(ncid,id,'valid_range',vals(1:2))
      end if
  
      if (present(scale_factor)) then
          vals(1) = scale_factor
          iret = nf90_put_att(ncid,id,'scale_factor',vals(1:1))
      end if
  
      if (present(add_offset)) then
          vals(1) = add_offset
          iret = nf90_put_att(ncid,id,'add_offset',vals(1:1))
      end if
  
      if (present(FillValue)) then
          vals(1) = FillValue
          iret = nf90_put_att(ncid,id,'_FillValue',vals(1:1))
      end if
  
      if (present(missing_value)) then
          vals(1) = missing_value
          iret = nf90_put_att(ncid,id,'missing_value',vals(1:1))
      end if
  
      set_attributes = 0
  
      return
  
      end function set_attributes
  !=======================================================================================================================
  
  
  
  
  
  
  !=======================================================================================================================
      subroutine check_err(status)
  
      integer, intent (in) :: status
  
      if (status .ne. NF90_NOERR) then
          print *, trim(nf90_strerror(status))
          stop
      endif
  
      end subroutine check_err
  !=======================================================================================================================
  
  
  
  
  !=======================================================================================================================
      subroutine svan(s, t, po, sigma)
  !/*
  !c------specific volume anomaly based on 1980 equation of state for
  !c      seawater and 1978 practical salinity scale
  !c      pressure          PO     decibars
  !c      temperature        T     degree celsius (IPTS-68)
  !c      salinitty          S     (PSS-78)
  !c      spec.vol.anom.  SVAN     1.0E-8 m**3/Kg
  !c      density anom.   SIGMA    Kg/m**3
  !*/
  
      real(rk)  p,sig,sr,r1,r2,r3,s,t,po,sigma
      real(rk)  a,b,c,d,e,a1,b1,aw,bw,k,ko,kw,k35,v350p,sva,dk,gam,pk,dr35p,dvan
  
      real(rk) r3500, r4 ,dr350
      data  r3500 /1028.1063/, r4/4.8314E-4/,dr350/28.106331/
  
      p=po/10.
      sr=sqrt(abs(s))
  
      r1= ((((6.536332E-9*t-1.120083E-6)*t+1.001685E-4)*t &
             -9.095290E-3)*t+6.793952E-2)*t-28.263737
      r2= (((5.3875E-9*t-8.2467E-7)*t+7.6438E-5)*t-4.0899E-3)*t &
            +8.24493E-1
      r3= (-1.6546E-6*t+1.0227E-4)*t-5.72466E-3
  
      sig=(r4*s + r3*sr + r2)*s +r1
  
      v350p=1.0/r3500
      sva=-sig*v350p/(r3500+sig)
      sigma= sig + dr350
  
      if (p.eq.0.0) return
  
      e = (9.1697E-10*t+2.0816E-8)*t-9.9348E-7
         bw = (5.2787E-8*t-6.12293E-6)*t+3.47718E-5
      b = bw + e*s
  
      d= 1.91075E-4
      c = (-1.6078E-6*t-1.0981E-5)*t+2.2838E-3
      aw = ((-5.77905E-7*t+1.16092E-4)*t+1.43713E-3)*t-0.1194975
      a = (d*sr + c)*s + aw
  
      b1 = (-5.3009E-4*t+1.6483E-2)*t+7.944E-2
      a1 = ((-6.1670E-5*t+1.09987E-2)*t-0.603459)*t+54.6746
      kw = (((-5.155288E-5*t+1.360477E-2)*t-2.327105)*t &
             +148.4206)*t-1930.06
      ko = (b1*sr + a1)*s + kw
  
      dk = (b*p+a)*p+ko
      k35 = (5.03217E-5*p+3.359406)*p+21582.27
      gam=p/k35
      pk=1.0-gam
      sva = sva * pk + (v350p+sva)*p*dk/(k35*(k35+dk))
  
      v350p= v350p*pk
      dr35p=gam/v350p
      dvan= sva/(v350p*(v350p+sva))
      sigma = dr350 + dr35p -dvan
      return
      end subroutine svan
  !=======================================================================================================================
  
  
      end module io_netcdf
  