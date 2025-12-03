! This file is part of Bottom RedOx Model (BROM, v.1.1).
! BROM is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the BROM distribution.
!-----------------------------------------------------------------------
! Original author(s): Evgeniy Yakushev, Elizaveta Protsenko, Phil Wallhead,
!                     Anfisa Berezina  
!-----------------------------------------------------------------------




    module calculate


        use fabm_types, only: rk
        use io_ascii, only: get_brom_par
    
        use fabm_omp, only: type_fabm_model => type_fabm_omp_model
    
        implicit none
        private
        public calculate_phys, calculate_sed, calculate_sed_eya, calculate_bubble
    
    
        contains
    
    !=======================================================================================================================
        subroutine calculate_phys(i, k_max, par_max, model, cc, kzti, fick, dcc, bctype_top, bctype_bottom, bc_top, bc_bottom, &
            surf_flux, bott_flux, bott_source, k_bbl_sed, dz, hz, kz, kz_mol, kz_bio, istep, julianday, id_O2, K_O2s, dt, freq_turb, &
            diff_method, cnpar, surf_flux_with_diff, bott_flux_with_diff, bioturb_across_SWI, pF1, pF2, phi_inv, is_solid, cc0)
    
        !Calculate vertical diffusion in the water column and sediments
    
        implicit none
    
        !Input variables
        integer, intent(in)                         :: k_max, par_max, k_bbl_sed, istep, julianday, id_O2, freq_turb
        integer, intent(in)                         :: surf_flux_with_diff, bott_flux_with_diff, bioturb_across_SWI
        integer, dimension(:,:), intent(in)         :: bctype_top, bctype_bottom
        real(rk), dimension(:,:), intent(in)        :: bc_top, bc_bottom, kz_bio, phi_inv
        real(rk), dimension(:,:,:), intent(in)      :: kz, kz_mol, pF1, pF2
        real(rk), dimension(:), intent(in)          :: dz, hz       !Grid spacings and layer thicknesses
        real(rk), intent(in)                        :: K_O2s        !Half-saturation constant for O2 effect on kz_bio
        real(rk), intent(in)                        :: dt
        integer, intent(in)                         :: diff_method  !Numerical method to treat vertical diffusion (set in brom.yaml):
        real(rk), intent(in)                        :: cnpar        !"Implicitness" parameter for GOTM vertical diffusion (set in brom.yaml)
        integer, dimension(:), intent(in)           :: is_solid     !Length par_max vector with 1s where variable is particulate, 0s elsewhere
        real(rk), intent(in)                        :: cc0          !Resilient concentration (same for all variables)
    
        !Output variables
        real(rk), dimension(:,:), intent(out)       :: surf_flux, bott_flux, bott_source
        real(rk), dimension(:,:,:), intent(inout)     :: kzti, fick, dcc    !Note: fick and dcc have units [mmol/m2/s] and [mmol/m3/s] respectively
    
        !Input/output variables
        class (type_fabm_model), pointer :: model
        real(rk), dimension(:,:,:), intent(inout)   :: cc
    
        !Local variables
        real(rk) :: h(0:k_max), nuY(0:k_max), Af(0:k_max), Vc(0:k_max), Y(0:k_max), ccnew1(1:k_max), pF(0:k_max) !GOTM routine inputs (Y also output)
        real(rk) :: Lsour(0:k_max), Qsour(0:k_max), Taur(0:k_max), Yobs(0:k_max) !Currently unused GOTM inputs
        real(rk) :: O2stat, dtt, Diffccup, Diffccdw, pFSWIup, pFSWIdw
        integer  :: i, k, ip, posconc, DiffBcup, DiffBcdw, idf, i_sed_top
    
    
    
        dtt = dt/freq_turb !Turbulence model time step [seconds]
    
        do idf=1,freq_turb !Do freq_turb time steps of size dt/freq_turb to perform turbulent mixing
    
            surf_flux = 0.0_rk
            if (surf_flux_with_diff.eq.1) call model%get_surface_sources(i, i, surf_flux(i:i,:))
            !Note: We MUST pass the range "i:i" to model%get_surface_sources( -- a single value "i" will produce compiler error!
            !Note: surf_flux has units [mass/m2/s] and is positive into the ocean (air->sea), see e.g. brom_carb.F90 and brom_bio.F90
            !!Calculate total vertical diffusivity [m2/s] and diffusive Courant number
            if(k_bbl_sed.ne.k_max) then
                O2stat = cc(i,k_bbl_sed,id_O2) / (cc(i,k_bbl_sed,id_O2) + K_O2s)   !Oxygen status of sediments set by O2 level just above sediment surface
            else
                O2stat = 0.0_rk
            endif
            do ip=1,par_max
                kzti(i,:,ip) = kz(i,:,istep) + kz_mol(i,:,ip) + kz_bio(i,:)*O2stat
                !Total diffusivity = turbulent (zero in sediments) + molecular (zero in water column) + bioturbation (zero in water column)
                !Note: molecular diffusivity is zero for variables that are solid in the sediments (see io_ascii.f90/make_physics_bbl_sed)
            end do
    
    
            !!Calculate Fickian diffusive fluxes (fick) if necessary
            if (diff_method.eq.0.or.idf.eq.freq_turb) then
                !Note this is done using centred differences at time t, irrespective of the numerical scheme for vertical diffusion
                !Units of fick are strictly [mass/unit total area/second]
    
                !Air-sea interface
                fick(i,1,:) = surf_flux(i,:)
    
                !Water column layer interfaces, not including the SWI
                do k=2,k_bbl_sed
    !                fick(i,k,:) = -1.0_rk * kzti(i,k,:) * (cc(i,k,:)-cc(i,k-1,:)) / dz(k-1)
                    fick(i,k,:) = -1.0_rk * 100000.0_rk*kzti(i,k,:) * (cc(i,k,:)-cc(i,k-1,:)) / dz(k-1) /100000.0_rk
                end do
    
                !Sediment-water interface (SWI)
                
                k = k_bbl_sed+1
                if (bioturb_across_SWI.eq.0) fick(i,k,:) = -1.0_rk * pF2(i,k,:) * 100000.0_rk*kzti(i,k,:) * (pF1(i,k,:)*cc(i,k,:)-cc(i,k-1,:)) / dz(k-1)/100000.0_rk
                !If not including bioturbation across the SWI, kzti will be zero for particles so we do not need to worry
                !about the Fickian gradient. Here fick will only represent molecular diffusion.
             if (k_bbl_sed.lt.k_max) then            
                if (bioturb_across_SWI.eq.1) then
                    fick(i,k,:) = -1.0_rk * pF2(i,k,:) * 100000.0_rk*kz_mol(i,k,:) * (pF1(i,k,:)*cc(i,k,:)-cc(i,k-1,:)) / dz(k-1)/100000.0_rk
                    !Molecular part is only applied to solutes (kz_mol = 0 for particulates); this mixes [mass/unit volume water] (intraphase mixing)
                    fick(i,k,:) = fick(i,k,:) - kz_bio(i,k) * O2stat * (cc(i,k,:)-cc(i,k-1,:)) / dz(k-1)
                    !Bioturbation part applies to both and mixes the concentrations in [mass/unit total volume] either side of SWI (interphase mixing)
                end if
             endif   
    
                !Sediment layer interfaces, not including the SWI
                do k=k_bbl_sed+2,k_max
                    fick(i,k,:) = -1.0_rk * pF2(i,k,:) * 100000.0_rk*kzti(i,k,:) * (pF1(i,k,:)*cc(i,k,:)-pF1(i,k-1,:)*cc(i,k-1,:)) / dz(k-1)/100000.0_rk
                    !Note: pF1 accounts for a porosity factor in the Fickian gradient for solutes and solids in the sediments
                    !      Since we are modelling as [mass per unit total volume], this is needed to convert to [mass per unit volume pore water] or [mass per unit volume solid]
                    !      in order to calculate the Fickian gradient for solutes and solids respectively (Berner 1980; Boudreau 1997; Soetaert et al. 1996).
                    !      pF1 = 1/phi for solutes, pF1 = 1/(1-phi) for solids
                    !Note: pF2 accounts for area restriction due to porosity and converts [mass/unit area pore water or solids/second] to [mass/unit total area/second]
                    !      pF2 = phi for solutes, pF2 = (1-phi) for solids
                end do
    
                !Bottom boundary
                if (bott_flux_with_diff.eq.1) then
                    bott_flux = 0.0_rk
                    bott_source = 0.0_rk                
    !                model%get_bottom_sources(i, i, bott_flux(i:i,:),bott_source(i:i,:))
                    fick(i,k_max+1,:) = bott_flux(i,:)
                else
                    fick(i,k_max+1,:) = 0.0_rk
                endif
            end if
    
    
            !!Perform turbulent diffusion update----------------------------------------------------------------------------
    
            if (diff_method.eq.0) then ! Use  FTCS------------------------------------------------------
                !     This is a FTCS (forward-time, central-space) algorithm
                !     It uses the fluxes fick in a consistent manner and therefore conserves mass
                !     Note that fluxes fick are defined on the layer interfaces:
                !     ============ kzti(1) (unused), fick(1) = surf_flux ( = 0 if surf_flux_with_diff = 0)
                !           o      cc(1), hz(1)
                !     ------------ kzti(2), fick(2)
                !           :
                !     ------------ kzti(k), fick(k)
                !           o      cc(k), hz(k)
                !     ------------ kzti(k+1), fick(k+1)
                !           :
                !     ============ kzti(k_max+1) (unused), fick(k_max+1) = bott_flux ( = 0 if bott_flux_with_diff = 0)
    
                !!Calculate tendencies dcc = dcc/dt = -dF/dz on layer midpoints (top/bottom not used when Dirichlet bc imposed)
                do k=1,k_max
                    dcc(i,k,:) = -1.0_rk * (fick(i,k+1,:)-fick(i,k,:)) / hz(k)
                end do
    
                !!Time integration
                !cc surface layer
                do ip=1,par_max
                    if (bctype_top(i,ip).gt.0) then
                        cc(i,1,ip) = bc_top(i,ip)
                    else
                        cc(i,1,ip) = cc(i,1,ip) + dtt*dcc(i,1,ip) !Simple Euler time step of size dtt = 86400.*dt/freq_turb [s]
                    end if
                end do
                !cc intermediate layers
                do k=2,(k_max-1)
                    cc(i,k,:) = cc(i,k,:) + dtt*dcc(i,k,:)
                end do
                !cc bottom layer
                do ip=1,par_max
                    if (bctype_bottom(i,ip).gt.0) then
                        cc(i,k_max,ip) = bc_bottom(i,ip)
                    else
                        cc(i,k_max,ip) = cc(i,k_max,ip) + dtt*dcc(i,k_max,ip)
                    end if
                end do
            end if!---------------------------------------------------------------------------------------------------------
    
    
    
            if (diff_method.gt.0) then !Use GOTM approach-------------------------------------------------------------------
                !The GOTM routines solve a semi-implicit scheme which is not limited by diffusion Courant number
                !In order to use the GOTM routines we must reorder vectors to conform with the GOTM grid,
                !which has concentrations (Y) and layer thicknesses (h) defined at layer centres and
                !diffusivity (nuY) defined at layer interfaces, with both indexed from bottom to top.
                !          GOTM             BROM
                !========= nuY(N)         = kzti(1)            (air-sea interface)
                !    o     Y(N), h(N)     = cc(1), hz(1)
                !--------- nuY(N-1)       = kzti(2)
                !    o     Y(N-1), h(N-1) = cc(2), hz(2)
                !    :                    :
                !    :                    :
                !========= nuY(i_sed_top) = kzti(k_bbl_sed+1)  (sediment-water interface)
                !    o     Y(i_sed_top)   = cc(k_bbl_sed+1)
                !    :                    :
                !    :                    :
                !--------- nuY(1)         = kzti(k_max)
                !    o     Y(1), h(1)     = cc(k_max), hz(k_max)
                !========= nuY(0)         = kzti(k_max+1)      (bottom in sediments)
                !
                !Note: N = k_max, so i_sed_top = k_max - k_bbl_sed
    
                h(1:k_max) = hz(k_max:1:-1) !Note: h(0) is expected by diff_center but is not used
                Lsour(0:k_max) = 0.0_rk
                Qsour(0:k_max) = 0.0_rk
                Taur(0:k_max) = 1.d20
                posconc = 1
                do ip=1,par_max
                    if (bctype_top(i,ip).gt.0) then
                        DiffBcup = 0 !Dirichlet (see util.f90)
                        Diffccup = bc_top(i,ip)
                    else
                        DiffBcup = 1 !Neumann
                        Diffccup = surf_flux(i,ip) !Note units [mass/m2/s] are correct since dtt [s] is passed to diff_center
                        !Note: In GOTM, fluxes *entering* a boundary cell are always positive (as in FABM), so no need for a minus sign
                    end if
                    if (bctype_bottom(i,ip).gt.0) then
                        DiffBcdw = 0 !Dirichlet (see util.f90)
                        Diffccdw = bc_bottom(i,ip)  !"
                    else
                        DiffBcdw = 1 !Neumann
                        Diffccdw = 0.0_rk
                    end if
    
                    if (bioturb_across_SWI.eq.0.and.diff_method.eq.1) then
                        !Use diff_center from GOTM lake, converting units of input and output
                        nuY(0:k_max) = kzti(i,k_max+1:1:-1,ip)
                        Af(0:k_max) = pF2(i,k_max+1:1:-1,ip) !pF2 = phi for solutes, (1-phi) for particulates, on layer interfaces
                        Vc(1:k_max) = hz(k_max:1:-1) / pF1(i,k_max:1:-1,ip)       !pF1 = 1/phi for solutes, 1/(1-phi) for particulates, on layer midpoints
                                                                                  !Note: pF(0) is expected by diff_center but is not used
                        Y(1:k_max) =  cc(i,k_max:1:-1,ip) * pF1(i,k_max:1:-1,ip)  !Note: Y(0) is expected by diff_center but is not used
                        !Input concentrations in units C' = C/phi [mass per unit volume pore water] or C' = C/(1-phi) [mass per unit volume sediments]
    
                        call diff_center(k_max,dtt,cnpar,posconc,h,Vc,Af,DiffBcup,DiffBcdw,          &
                                    Diffccup,Diffccdw,nuY,Lsour,Qsour,Taur,Yobs,Y)
    
                        ccnew1(1:k_max) = Y(k_max:1:-1) / pF1(i,1:k_max,ip) !Updated concentrations in [mass per unit total volume]
                        dcc(i,1:k_max,ip) = (ccnew1(1:k_max)-cc(i,1:k_max,ip))/dtt
                        cc(i,1:k_max,ip) = ccnew1(1:k_max)
                    end if
    
                    if (bioturb_across_SWI.eq.0.and.diff_method.eq.2) then
                        !Use diff_center1 derived from original GOTM, solving dC/dt directly without units conversion
                        nuY(0:k_max) = pF2(i,k_max+1:1:-1,ip) * kzti(i,k_max+1:1:-1,ip) !pF2 = phi for solutes, (1-phi) for particulates, on layer interfaces
                        pF(1:k_max) = pF1(i,k_max:1:-1,ip) !pF1 = 1/phi for solutes, 1/(1-phi) for particulates, on layer midpoints
                                                           !Note: pF(0) is expected by diff_center but is not used
                        Y(1:k_max) = cc(i,k_max:1:-1,ip)   !Note: Y(0) is expected by diff_center but is not used
    
                        call diff_center1(k_max,dtt,cnpar,posconc,h,DiffBcup,DiffBcdw,          &
                                    Diffccup,Diffccdw,nuY,Lsour,Qsour,Taur,Yobs,pF,Y)
    
                        dcc(i,1:k_max,ip) = (Y(k_max:1:-1)-cc(i,1:k_max,ip))/dtt
                        cc(i,1:k_max,ip) = Y(k_max:1:-1) !Updated concentrations
                    end if
    
                    if (bioturb_across_SWI.eq.1) then
                        !Here we must use a modified diff_center (diff_center2), solving dC/dt directly without units conversion
                        i_sed_top = k_max - k_bbl_sed
                        nuY(0:k_max) = pF2(i,k_max+1:1:-1,ip) * kzti(i,k_max+1:1:-1,ip) !pF2 = phi for solutes, (1-phi) for particulates, on layer interfaces
                        pF(1:k_max) = pF1(i,k_max:1:-1,ip)  !pF1 = 1/phi for solutes, 1/(1-phi) for particulates, on layer midpoints
                                                            !Note: pF(0) is expected by diff_center but is not used
                        Y(1:k_max) = cc(i,k_max:1:-1,ip)    !Note: Y(0) is expected by diff_center but is not used
                    if (k_bbl_sed.lt.k_max) then
                        if (is_solid(ip).eq.0) then
                            !fick_SWI = -phi_0*Km/dz*(C_1/phi_1 - C_-1) - Kb/dz*(C_1 - C_-1)   [intraphase molecular + interphase bioturb]
                            !         = -phi_0*(Km+Kb)/dz * (p_1*C_1 - p_-1*C_-1)
                            !where: p_1  = (phi_0/phi_1*Km + Kb) / (phi_0*(Km+Kb))
                            !       p_-1 = (phi_0*Km + Kb) / (phi_0*(Km+Kb))
                            !and subscripts refer to SWI (0), below (1), and above (-1)
                            pFSWIup = (pF2(i,k_bbl_sed+1,ip)*kz_mol(i,k_bbl_sed+1,ip) + kz_bio(i,k_bbl_sed+1)*O2stat) / &
                                      (pF2(i,k_bbl_sed+1,ip)*(kz_mol(i,k_bbl_sed+1,ip) + kz_bio(i,k_bbl_sed+1)*O2stat))
                            pFSWIdw = (pF1(i,k_bbl_sed+1,ip)*pF2(i,k_bbl_sed+1,ip)*kz_mol(i,k_bbl_sed+1,ip) + kz_bio(i,k_bbl_sed+1)*O2stat) / &
                                      (pF2(i,k_bbl_sed+1,ip)*(kz_mol(i,k_bbl_sed+1,ip) + kz_bio(i,k_bbl_sed+1)*O2stat))
                        else
                            !fick_SWI = -Kb/dz*(C_1 - C_-1)   [interphase bioturb]
                            !         = -(1-phi_0)*Kb/dz * (p_1*C_1 - p_-1*C_-1)
                            !where: p_1  = 1 / (1 - phi_0)
                            !       p_-1 = 1 / (1 - phi_0)
                            !and subscripts refer to SWI (0), below (1), and above (-1)
                            pFSWIup = 1.0_rk / (1.0_rk - pF2(i,k_bbl_sed+1,ip))
                            pFSWIdw = pFSWIup
                        end if
                    endif
    
                        call diff_center2(k_max,dtt,cnpar,posconc,h,DiffBcup,DiffBcdw,          &
                                    Diffccup,Diffccdw,nuY,Lsour,Qsour,Taur,Yobs,pF,             &
                                    pFSWIup, pFSWIdw, i_sed_top, Y)
    
                        dcc(i,1:k_max,ip) = (Y(k_max:1:-1)-cc(i,1:k_max,ip))/dtt
                        cc(i,1:k_max,ip) = Y(k_max:1:-1) !Updated concentrations
                    end if
                end do
            end if!---------------------------------------------------------------------------------------------------------
    
            cc(i,:,:) = max(cc0, cc(i,:,:)) !Impose resilient concentration
        end do
    
            end subroutine calculate_phys
    !=======================================================================================================================
    
    
    
    
    
    
    
    
    
    !=======================================================================================================================
        subroutine calculate_sed(i, k_max, par_max, model, cc, wti, sink, dcc, dcc_R, bctype_top, bctype_bottom, &
            bc_top, bc_bottom, hz, dz, k_bbl_sed, wbio, w_b, u_b, julianday, dt, freq_sed, dynamic_w_sed, is_solid, &
            rho, phi1, fick, k_sed1, K_O2s, kz_bio, id_O2, dphidz_SWI, cc0, bott_flux, bott_source)
    
        !Calculates vertical advection (sedimentation) in the water column and sediments
    
        implicit none
    
        !Input variables
        integer, intent(in)                         :: k_max, par_max, id_O2
        integer, intent(in)                         :: julianday, freq_sed, k_bbl_sed, dynamic_w_sed
        integer, dimension(:,:), intent(in)         :: bctype_top, bctype_bottom
        integer, dimension(:), intent(in)           :: is_solid, k_sed1
        real(rk), dimension(:,:,:), intent(in)      :: dcc_R, fick
        real(rk), dimension(:,:), intent(in)        :: bc_top, bc_bottom, phi1, w_b, u_b, kz_bio
        real(rk), dimension(:), intent(in)          :: hz, dz, rho
        real(rk), intent(in)                        :: dt, K_O2s, dphidz_SWI, cc0
    
        !Output variables
        real(rk), dimension(:,:,:), intent(in)      :: wbio
        real(rk), dimension(:,:,:), intent(inout)   :: wti
        real(rk), dimension(:,:,:), intent(inout)   :: sink, dcc
        real(rk), dimension(:,:), intent(out)       :: bott_flux, bott_source
     
        !Input/output variables
        class (type_fabm_model), pointer :: model
        real(rk), dimension(:,:,:), intent(inout)   :: cc
    
        !Local variables
        real(rk) :: dtt, O2stat
        real(rk) :: w_1m(k_max+1,par_max), w_1(k_max+1), u_1(k_max+1), w_1c(k_max+1), u_1c(k_max+1)
        integer  :: i, k, ip, idf
    
    
        dtt = dt/freq_sed !Sedimentation model time step [seconds]
    !    sink = 0.0_rk
        w_1 = 0.0_rk
        u_1 = 0.0_rk
        w_1c = 0.0_rk
        u_1c = 0.0_rk
        sink(i,:,:) = 0.0_rk
    
        !Compute vertical velocity components due to modelled (reactive) particles, if required
        !
        !Assuming no externally impressed porewater flow, the equations for liquid and solid volume fractions are:
        !
        !dphi/dt = -d/dz(phi*u - Dbip*dphi/dz) - sum_i Rp_i/rhop_i                (1)
        !d(1-phi)/dt = -d/dz((1-phi)*w - Dbip*d/dz(1-phi)) + sum_i Rp_i/rhop_i    (2)
        !
        !where u is the solute advective velocity, w is the particulate advective velocity,
        !Dbip is the interphase component of bioturbation diffusivity ( = Db at SWI, 0 elsewhere)
        !Rp_i is the net reaction term for the i^th particulate substance [mmol/m3/s]
        !rhop_i is the molar density of the i^th particulate substance [mmol/m3]
        !
        !(1)+(2) => phi*u + (1-phi)*w = const                                     (3)
        !
        !Given u = w at some (possibly infinite) depth where compaction ceases and phi = phi_inf, w = w_inf:
        !phi*u + (1-phi)*w = phi_inf*w_inf + (1-phi_inf)*w_inf
        !        => phi*u = w_inf - (1-phi)*w                                     (4)
        !
        !We calculate w(z) by assuming steady state compaction (dphi/dt = 0)
        !and matching the solid volume flux across the SWI with the (approximated) sinking flux of suspended particles in the fluff layer:
        !
        !(2) -> (1-phi)*w + Dbip*dphi/dz = (1-phi_0)*w_0 + Dbip*dphi/dz_0 + sum_i (1/rhop_i)*int_z0^z Rp_i(z') dz'      (5)
        !
        !Matching the volume flux across the SWI gives:
        !       (1-phi_0)*w_0 + Dbip*dphi/dz_0 = F_b0 + sum_i (1/rhop_i)*wbio_f(i)*Cp_sf(i)
        !
        !where F_b0 is the background volume flux [m/s]
        !      wbio_f(i) is the sinking speed in the fluff layer [m/s]
        !      Cp_sf(i) is the suspended particle concentration in the fluff layer [mmol/m3]
        !      (approximated by the minimum of the particle concentration in the fluff layer and layer above)
        !
        !
        !For dynamic_w_sed = 0, we set the reactions and modelled volume fluxes in (5) to zero:
        !
        !(1-phi)*w + Dbip*dphi/dz = F_b0 = (1-phi_inf)*w_binf
        !
        !Then using (4):
        !
        !phi*u = w_inf - (1-phi)*w = phi_inf*w_binf + Dbip*dphi/dz
        !
        !Decomposing the velocities as (w,u) = (w,u)_b + (w,u)_1, we get:
        !
        !(1-phi)*w_1 = -Dbip*dphi/dz
        !phi*u_1     = Dbip*dphi/dz
        !
        !where (1-phi)*w_b = (1-phi_inf)*w_binf and phi*u_b = phi_inf*w_binf
        !
        !
        !For dynamic_w_sed = 1, we have further corrections (w,u)_1c, where:
        !
        !(1-phi)*w_1c = sum_i [ (1/rhop_i)*(wbio_f(i)*Cp_sf(i) + int_zTF^z Rp_i(z') dz') ]
        !
        !and using (4):
        !
        !phi*u_1c     = w_1cinf - (1-phi)*w_1c
        !
        !where w_1cinf can be approximated by the deepest value of w_1c
        !
        if(k_bbl_sed.ne.k_max) then  
            O2stat = cc(i,k_bbl_sed,id_O2) / (cc(i,k_bbl_sed,id_O2) + K_O2s)   !Oxygen status of sediments set by O2 level just above sediment surface
        else
            O2stat = 0.0_rk
        endif
        w_1(k_bbl_sed+1) = -1.0_rk*O2stat*kz_bio(i,k_bbl_sed+1)*dphidz_SWI / (1.0_rk-phi1(i,k_bbl_sed+1))
        u_1(k_bbl_sed+1) = O2stat*kz_bio(i,k_bbl_sed+1)*dphidz_SWI / phi1(i,k_bbl_sed+1)
        if (dynamic_w_sed.eq.1) then
            w_1m = 0.0_rk
            do ip=1,par_max !Sum over contributions from each particulate variable
                if (is_solid(ip).eq.1) then
                    !First set rhop_i*(1-phi)*w_1i at the SWI
                    w_1m(k_bbl_sed+1,ip) = wbio(i,k_bbl_sed,ip)*min(cc(i,k_bbl_sed,ip),cc(i,k_bbl_sed-1,ip))
                    !Now set rhop_i*(1-phi)*w_1i in the sediments by integrating the reaction terms
                    do k=k_bbl_sed+2,k_max+1
                        w_1m(k,ip) = w_1m(k-1,ip) + dcc_R(i,k-1,ip)*hz(k-1)
                    end do
                    w_1m(k_sed1,ip) = w_1m(k_sed1,ip) / (rho(ip)*(1.0_rk-phi1(i,k_sed1))) !Divide by rhop_i*(1-phi) to get w_1c(i)
                    w_1c(k_sed1) = w_1c(k_sed1) + w_1m(k_sed1,ip)                         !Add to total w_1c
                end if
            end do
            !Now calculate u from w using (4) above
            u_1c(k_sed1) = (w_1c(k_max+1) - (1.0_rk-phi1(i,k_sed1))*w_1(k_sed1)) / phi1(i,k_sed1)
        end if
    
    
        !Interpolate wbio from FABM (defined on layer midpoints, as for concentrations) to wti on the layer interfaces
        wti(i,1,:) = 0.0_rk
        do ip=1,par_max
            !Air-sea interface (unused)
            wti(i,1,ip) = wbio(i,1,ip)
            !Water column layer interfaces, not including SWI
            wti(i,2:k_bbl_sed,ip) = wbio(i,1:k_bbl_sed-1,ip) + 0.5_rk*hz(1:k_bbl_sed-1)*(wbio(i,2:k_bbl_sed,ip) - wbio(i,1:k_bbl_sed-1,ip)) / dz(1:k_bbl_sed-1)
            !Sediment layer interfaces, including SWI
            if (is_solid(ip).eq.0) then
                wti(i,k_sed1,ip) = u_b(i,k_sed1) + u_1(k_sed1) + u_1c(k_sed1)
            else
                wti(i,k_sed1,ip) = w_b(i,k_sed1) + w_1(k_sed1) + w_1c(k_sed1)
            end if
            wti(i,k_max+1,ip) = wti(i,k_max,ip)
        end do
    
        !! Perform advective flux calculation and cc update
        !This uses a simple first order upwind differencing scheme (FUDM)
        !It uses the fluxes sink in a consistent manner and therefore conserves mass
        do idf=1,freq_sed
            !!Calculate sinking fluxes at layer interfaces (sink, units strictly [mass/unit total area/second])
            !Air-sea interface
            sink(i,1,:) = 0.0_rk
            !Water column and sediment layer interfaces
            do k=2,k_max+1
                sink(i,k,:) = wti(i,k,:)*cc(i,k-1,:)
                !This is an upwind differencing approx., hence the use of cc(k-1)
                !Note: factors phi, (1-phi) are not needed in the sediments because cc is in units [mass per unit total volume]
            end do
    
    
            !!Calculate tendencies dcc = dcc/dt = -dF/dz on layer midpoints (top/bottom not used where Dirichlet bc imposed)
            do k=1,k_max
                dcc(i,k,:) = -1.0_rk * (sink(i,k+1,:)-sink(i,k,:)) / hz(k)
            end do
    
    
            !!Time integration
            !cc surface layer
            do ip=1,par_max
                if (bctype_top(i,ip).gt.0) then
                    cc(i,1,ip) = bc_top(i,ip)
                else
                    cc(i,1,ip) = cc(i,1,ip) + dtt*dcc(i,1,ip) !Simple Euler time step of size dtt = 86400.*dt/freq_sed [s]
                end if
            end do
            !cc intermediate layers
            do k=2,(k_max-1)
                cc(i,k,:) = cc(i,k,:) + dtt*dcc(i,k,:)
            end do
            !cc bottom layer
            bott_flux = 0.0_rk
            bott_source = 0.0_rk       
            call model%get_bottom_sources(i, i, bott_flux(i:i,:),bott_source(i:i,:))
            do ip=1,par_max
                if (is_solid(ip).eq.1) then            
                    dcc(i,k_max,ip) = dcc(i,k_max,ip) + bott_flux(i,ip) / hz(k_max)
                    sink(i,k_max+1,:) = bott_flux(i,:)
                endif
            end do
    !            cc(i,k_max,ip) = cc(i,k_max,ip) + dtt*dcc(i,k_max,ip)               
            do ip=1,par_max
                if (bctype_bottom(i,ip).gt.0) then
                    cc(i,k_max,ip) = bc_bottom(i,ip)
                else
                    cc(i,k_max,ip) = cc(i,k_max,ip) + dtt*dcc(i,k_max,ip)
                end if
            end do
            cc(i,:,:) = max(cc0, cc(i,:,:)) !Impose resilient concentration
        end do
    
            end subroutine calculate_sed
    !=======================================================================================================================
    
    
    
    
    
    !=======================================================================================================================
        subroutine calculate_sed_eya(i, k_max, par_max, model, cc, wti, sink, &
            dcc, dVV, bctype_top, bctype_bottom, bc_top, bc_bottom, &
            hz, dz, k_bbl_sed, wbio, w_b, u_b, julianday, dt, freq_sed, &
            dynamic_w_sed, constant_w_sed, is_solid, rho, phi1, fick, &
            k_sed1, K_O2s, kz_bio, id_O2, dphidz_SWI, &
            cc0, bott_flux, bott_source, w_binf, bu_co, is_gas)
    
        !Calculates vertical advection (sedimentation) in the water column and sediments
    
        implicit none
    
        !Input variables
        integer, intent(in)                         :: k_max, par_max, id_O2, julianday, freq_sed
        integer, intent(in)                         :: k_bbl_sed, dynamic_w_sed, constant_w_sed
        integer, dimension(:,:), intent(in)         :: bctype_top, bctype_bottom
        integer, dimension(:), intent(in)           :: is_solid, is_gas, k_sed1
        real(rk), dimension(:,:,:), intent(in)      :: fick
        real(rk), dimension(:,:), intent(in)        :: bc_top, bc_bottom, phi1, w_b, u_b, kz_bio
        real(rk), dimension(:), intent(in)          :: hz, dz, rho
        real(rk), intent(in)                        :: dt, K_O2s, dphidz_SWI, cc0
        real(rk), intent(in)                        :: w_binf, bu_co
    
        !Output variables
        real(rk), dimension(:,:,:), intent(in)      :: wbio
        real(rk), dimension(:,:,:), intent(inout)   :: wti
        real(rk), dimension(:,:,:), intent(inout)   :: sink, dcc
        real(rk), dimension(:,:), intent(out)       :: bott_flux, bott_source
     
        !Input/output variables
        class (type_fabm_model), pointer :: model
        real(rk), dimension(:,:,:), intent(inout)   :: cc, dVV
    
        !Local variables
        real(rk) :: dtt, O2stat
        real(rk) :: w_1m(k_max+1,par_max), w_1(k_max+1), u_1(k_max+1), w_1c(k_max+1), u_1c(k_max+1)
        integer  :: i, k, ip, idf
    
    
        dtt = dt/freq_sed !Sedimentation model time step [seconds]
    !    sink = 0.0_rk
        w_1 = 0.0_rk
        u_1 = 0.0_rk
        w_1c = 0.0_rk
        u_1c = 0.0_rk
        sink(i,:,:) = 0.0_rk
        wti(i,:,:) = 0.0_rk 
!--------------------------------------------------------------------   
        ! wti() at water column layer midpoints including Air-sea interface (unused) 
        !    and not including SWI:  (as wbio in FABM) 
        !    wti(i,1:k_bbl_sed,:) = max(0.0_rk,wbio(i,1:k_bbl_sed,:))
            wti(i,1:(k_bbl_sed-1),:) = wbio(i,1:(k_bbl_sed-1),:)
            wti(i,1:(k_bbl_sed-1),:) = wbio(i,1:(k_bbl_sed-1),:)
            wti(i,1,:) = 0.0_rk !as boundary condition 
        ! check for light microplast and/or gas
        ! check for light microplast and/or gas
            do ip=1,par_max
               if (wti(i,2,ip)<0.0_rk)  then
                    wti(i,2,ip) = 0.0_rk  !check for is_gas??
               endif
            enddo
!--------------------------------------------------------------------   
        if (constant_w_sed.eq.1) then  ! case constant burial velosity
            ! wti() at sediment layer midpoints: w_b, backgound burying velocity
            do k=k_bbl_sed,k_max+1
                do ip=1,par_max
                   wti(i,k,ip) = w_binf
                enddo
            enddo
        endif
!--------------------------------------------------------------------        
        ! wti() increases on burying rate, that depends of change of particles volume in the layer above the SWI 
        ! and can not be higher than defined by CFL criteria)
        if (dynamic_w_sed.eq.1) then ! case  burial velosity depending on particles accumulation above SWI
        ! we accelerate burying rate due to an increase of particles volume dVV()[m3/sec] in water layer just above SWI
            do ip=1,par_max
                do k=k_bbl_sed,k_max+1       
                    wti(i,k,ip) = wti(i,k,ip) + dVV(i,k_bbl_sed,1)
                    ! dVV is in m/s
                end do
            enddo
        endif
!--------------------------------------------------------------------           
            ! Apply "Burial coeficient" to increase wti () exactly to the SWI at proportional to the
            !    settling velocity in the water column (0<bu_co<1) or (if bu_co <1) calulate it as a ratio of dVV/hz)
            !    settling velocity in the water column (0<bu_co<1) or (if bu_co <1) calulate it as a ratio of dVV/hz)
        do ip=1,par_max
            if(bu_co.gt.0) then
                wti(i,k_bbl_sed,ip) = wti(i,k_bbl_sed,ip) + bu_co*wti(i,k_bbl_sed-1,ip)
            else
                wti(i,k_bbl_sed,ip) = wti(i,k_bbl_sed,ip) +  dVV(i,k_bbl_sed,1)/hz(k_bbl_sed)*wti(i,k_bbl_sed-1,ip) 
                ! dVV/hz is the share of particles present in the bottom layer that should be burried  
            endif    
        end do        
!--------------------------------------------------------------------             
        do ip=1,par_max
          if (is_gas(ip).eq.1) then
              wti(:,:,ip)=0.0_rk
          endif
        enddo  
    
        ! Perform sinking advective flux calculation and cc update
        ! This uses a simple first order upwind differencing scheme (FUDM)
        ! It uses the fluxes sink in a consistent manner and therefore conserves mass
        do idf=1,freq_sed
            !Calculate sinking fluxes at layer interfaces (sink, units strictly [mass/unit total area/second])
            sink(i,1,:) = 0.0_rk
            !Water column and sediment layer interfaces
            do ip=1,par_max
                do k=1,k_max
                    sink(i,k,ip) = ((100000.0_rk*wti(i,k,ip))*cc(i,k,ip))/100000.0_rk
                end do
            enddo
    
            !!Calculate tendencies dcc = dcc/dt = -dF/dz on layer midpoints (top/bottom not used where Dirichlet bc imposed)
            do k=2,k_max
                dcc(i,k,:) = (sink(i,k-1,:)-sink(i,k,:)) / hz(k) !-1)
            end do
                dcc(i,1,:) = - sink(i,1,:)/ hz(1)
                dcc(i,k_max,:) = -sink(i,k_max,:) / hz(k_max) !0.0
            !!Time integration
            !cc water surface   
            do ip=1,par_max
                if (bctype_top(i,ip).gt.0) then
                    cc(i,1,ip) = bc_top(i,ip)
                else
                    cc(i,1,ip) = cc(i,1,ip) + dtt*dcc(i,1,ip) !Simple Euler time step of size dtt = dt/freq_sed [s]
                end if
            end do
            !cc intermediate layers in the water and  BBL+sediments (if considered)
            do k=2,(k_max-1)
                cc(i,k,:) = cc(i,k,:) + dtt*dcc(i,k,:)
            end do
            !cc bottom layer
            bott_flux = 0.0_rk
            bott_source = 0.0_rk       
            call model%get_bottom_sources(i, i, bott_flux(i:i,:),bott_source(i:i,:))
            do ip=1,par_max
                if (is_solid(ip).eq.1) then   
                    dcc(i,k_max,ip) = dcc(i,k_max,ip) + bott_flux(i,ip) / hz(k_max)
     !               dcc(i,k_max,ip) = dcc(i,k_max-1,ip) + bott_flux(i,ip) / hz(k_max) ! with  Matvey
                    sink(i,k_max+1,:) = bott_flux(i,:)
                endif
            end do
    !                   
            do ip=1,par_max
                if (bctype_bottom(i,ip).gt.0) then
                    cc(i,k_max,ip) = bc_bottom(i,ip)
                else
                    cc(i,k_max,ip) = cc(i,k_max,ip) + dtt*dcc(i,k_max,ip)
                end if
            end do
            cc(i,:,:) = max(cc0, cc(i,:,:)) !Impose resilient concentration
        end do
    
        end subroutine calculate_sed_eya
    !=======================================================================================================================
            
            
            
            
    !=======================================================================================================================
        subroutine calculate_bubble(i_min, i_max, k_max, par_max, model, cc, sink, N_bubbles, &
            dcc, hz, dz, z, t, k_bbl_sed, julianday, dt, freq_float, is_gas, wbio, cc0, use_hice, aice)
    
        !Calculates floating of bubbles in the water column and sediments
    
        implicit none
    
        !Input variables
        integer, intent(in)                         :: i_min, i_max, k_max, par_max, julianday, freq_float
        integer, intent(in)                         :: k_bbl_sed, use_hice
        integer, dimension(:), intent(in)           :: is_gas 
        real(rk), dimension(:), intent(in)          :: hz, dz, z, aice
        real(rk), dimension(:,:,:), intent(in)      :: t
        real(rk), intent(in)                        :: dt, cc0, N_bubbles
    
        !Output variables
        real(rk), dimension(:,:,:), intent(in)      :: wbio
        real(rk), dimension(:,:,:), intent(out)     :: sink, dcc
     !   real(rk), dimension(:,:), intent(out)       :: bott_flux, bott_source
     
        !Input/output variables
        class (type_fabm_model), pointer :: model
        real(rk), dimension(:,:,:), intent(inout)   :: cc
    
        !Local variables
        real(rk) :: dtt, rb, wbub(k_max+1)
        integer  :: i, k, ip, idf
    
    !$OMP PARALLEL DEFAULT(SHARED)
    
        dtt = dt/freq_float !Sedimentation model time step [seconds]
    !    sink(i,:,:) = 0.0_rk
    
        ! Perform rise advective flux calculation and cc update
        ! This uses a simple first order upwind differencing scheme (FUDM)
        ! It uses the fluxes sink in a consistent manner and therefore conserves mass
        wbub= 0.0_rk
    
        do ip=1,par_max
          if (is_gas(ip).eq.1) then
            do idf=1,freq_float
        !Calculate floating fluxes at layer interfaces (sink, units strictly [mass/unit total area/second])
    !$OMP DO PRIVATE(i,rb,wbub)
              do k=1,k_max-1
                 do i=i_min, i_max
                rb = 100._rk * &     !bubble radius in [cm] (details in brom_bubble), we convert to [m]
                    (cc(i,k,ip)*0.001_rk*8.314_rk*(273.15_rk+t(i,k,julianday))/(101325.0_rk &
                    +z(k)*101325.0_rk/10.3_rk) &                        ! Volume of gas
                    /N_bubbles*3.0_rk/4.0_rk/3.14159_rk)**(1.0_rk/3.0_rk)
    !            if(rb.gt.0.0000001) then
    !              if (rb.lt.0.4) then
    !                wbub(k)=  -(22.16_rk+0.733_rk*(rb-0.0000000584_rk)**(-0.0849_rk)) &   ! rate of rise
    !                         *exp(4.792e-4_rk*t(i,k,julianday)*(rb-0.0000000584_rk)**(-0.815_rk)) &
    !                         /100._rk !into m/s          
    !              else
    !                wbub(k)= -19.16_rk/100._rk !into m/s
    !              endif
    !!                if(wbub(k).lt.wbio(i,k,ip)) wbub(k)=wbio(i,k,ip)
    !!                if(wbub(k).gt.0.0) wbub(k)=0.0_rk
    !            endif
                  if(rb.gt.0.01) then ! in cm
                    if(rb.lt.1.08) then !with 0.01_rk* we convert to [m/s]
                      wbub(k) = -0.01_rk*rb*(276._rk+ rb*(-1648._rk+rb*(4882._rk &   !Clift et al., 1978
                                 +rb*(-7429._rk+rb*(5618._rk+rb*(-1670._rk))))))     ! cm/s
                    else
                      wbub(k) = -0.01_rk*23._rk ! 23 cm/s, rate typical for r>4 mm (Leifer et al., 2002)
                    !wbub(k) = (22.16_rk+0.773_rk*((rb-0.0584)**(-0.849)))*exp(4.792e-4*20._rk*(rb-0.0584)**(-0.815))
                    !   ! for r<4mm Leifer,Patro, 2002, wbub[cm/s], rb[cm]  !  for ToC=20. 
                    !wbub(k) = (19.16_rk+11.05_rk*(rb-0.0584)**(-0.0))   !for rb=4-40 mm, Leifer,Patro, 2002
                    endif                          
                    !if(wbub(k).lt.wbio(i,k,ip)) then
                    !    wbub(k)=wbio(i,k,ip)
                    !endif
                  endif
                ! if(rb.gt.0.01) then
                !    write (*,*) "k=", k,"i=", i,"rb[cm]=",rb,"wbub[m/s]=", wbub(k)
                ! endif
                  wbub(k)=min(wbub(k),0.0_rk)
                sink(i,k,ip) = -100000.0_rk*wbub(k)*cc(i,k,ip)/100000.0_rk
                end do
              end do
    
              if (use_hice.eq.1.and.aice(julianday).eq.0.0_rk)   sink(i,1,ip) =0.0_rk      
        !Calculate tendencies dcc = dcc/dt = -dF/dz on layer midpoints (top/bottom not used where Dirichlet bc imposed)
    !$OMP DO PRIVATE(i)
              do k=1,k_max-1
                 do i=i_min, i_max
                dcc(i,k,ip) = (sink(i,k+1,ip)-sink(i,k,ip))/hz(k) !-1)
        !Integrate cc() on layer midpoints 
                cc(i,k,ip) = cc(i,k,ip) + dtt*dcc(i,k,ip)
                end do
              end do
            enddo
          endif
          enddo
    
        !Impose resilient concentration
    !$OMP DO
        do ip=1,par_max
          cc(:,:,ip) = max(cc0, cc(:,:,ip)) 
        end do
    
    !$OMP END PARALLEL
    
            end subroutine calculate_bubble
    !===================================================================================
    
    
        end module calculate
    