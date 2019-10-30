!-----------------------------------------------------------------------
! BROM2 is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the BROM2 distribution.
!-----------------------------------------------------------------------

#include "fabm_driver.h"

module fabm_niva_brom_carbon
  use fabm_types

  implicit none
  private
  type,extends(type_base_model),public:: type_niva_brom_carbon
    !all descriptions are in the initialize subroutine
    type(type_state_variable_id):: id_DIC,id_Alk

    type(type_diagnostic_variable_id):: id_pCO2
    type(type_diagnostic_variable_id):: id_CO3
    !diagnostic dependencies
    type(type_dependency_id):: id_Hplus,id_Kc0,id_Kc1,id_Kc2
    !for do_surface
    type(type_dependency_id):: id_temp,id_salt,id_pCO2w
    type(type_horizontal_dependency_id):: id_windspeed,id_pCO2a
  contains
    procedure :: initialize
    procedure :: do_surface
    procedure :: do
  end type
contains
  !
  !
  !
  subroutine initialize(self,configunit)
    class(type_niva_brom_carbon),intent(inout),target :: self
    integer,                     intent(in)           :: configunit

    call self%register_state_variable(&
         self%id_DIC,'DIC','mmol/m**3','DIC',minimum=0.0_rk)
    call self%add_to_aggregate_variable(&
         standard_variables%total_carbon,self%id_DIC)
    call self%register_state_variable(&
         self%id_Alk,'Alk','umol/kg','Alk',2300._rk,minimum=1.e-4_rk,&
         standard_variable=&
         standard_variables%alkalinity_expressed_as_mole_equivalent)

    !register diagnostic variables
    call self%register_diagnostic_variable(&
         self%id_pCO2,'pCO2','ppm','partial pressure of CO2')
    call self%register_diagnostic_variable(&
         self%id_CO3,'CO3','mmol/m**3','CO3--',standard_variable=&
         standard_variables%&
         mole_concentration_of_carbonate_expressed_as_carbon)

    !dependencies
    !for do_surface
    call self%register_dependency(&
         self%id_temp,standard_variables%temperature)
    call self%register_dependency(&
         self%id_salt,standard_variables%practical_salinity)
    call self%register_dependency(&
         self%id_pco2a,&
         standard_variables%mole_fraction_of_carbon_dioxide_in_air)
    call self%register_dependency(&
         self%id_windspeed,standard_variables%wind_speed)
    call self%register_dependency(&
         self%id_pCO2w,'pCO2','ppm','partial pressure of CO2')

    !diagnostic variables
    call self%register_dependency(self%id_Hplus,'Hplus','mmol/m**3','H+')
    call self%register_dependency(&
         self%id_Kc0,'Kc0','-','Henry''s constant')
    call self%register_dependency(&
         self%id_Kc1,'Kc1','-','[H+][HCO3-]/[H2CO3]')
    call self%register_dependency(&
         self%id_Kc2,'Kc2','-','[H+][CO3--]/[HCO3-]')
  end subroutine initialize
  !
  !Sea water CO2 exchange, adapted from PML's ERSEM code
  !
  subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
    class (type_niva_brom_carbon),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_SURFACE_
    real(rk):: Q_pCO2, Q_DIC
    real(rk):: temp, salt, Kc0
    real(rk):: Sc, fwind !PML
    real(rk):: pCO2w, pCO2a
    real(rk):: windspeed

    _HORIZONTAL_LOOP_BEGIN_
      _GET_(self%id_temp,temp) !temperature
      _GET_(self%id_salt,salt) !salinity
      !previouss pCO2 which calculates here in the subroutine do
      _GET_(self%id_pCO2w,pCO2w) 
      !Equilibrium constants
      _GET_(self%id_Kc0,Kc0)
      _GET_HORIZONTAL_(self%id_windspeed,windspeed)
      _GET_HORIZONTAL_(self%id_pCO2a,pCO2a) !from atmosphere

      !PML
      !calculate the scmidt number and unit conversions
      Sc = 2073.1_rk-125.62_rk*temp+3.6276_rk*temp**2._rk-0.043219_rk*&
           temp**3.0_rk
      fwind = (0.222_rk*windspeed**2_rk+0.333_rk*windspeed)*&
              (Sc/660._rk)**(-0.5_rk)
      fwind=fwind*24._rk/100._rk !convert to m/day
      !flux depends on the difference in partial pressures, wind and henry
      !here it is rescaled to mmol/m2/d
      !flux = fwind * HENRY * ( PCO2A - PCO2W ) * dcf
      Q_pCO2= fwind*(pCO2a-max(0e0,pCO2w))
      Q_DIC = Q_pCO2*Kc0/86400._rk

      _SET_SURFACE_EXCHANGE_(self%id_DIC,Q_DIC)
    _HORIZONTAL_LOOP_END_
  end subroutine do_surface
  !
  !concentrations in uM, but constants in moles/kg
  !
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_carbon),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_
    !state variables
    real(rk):: DIC
    !diagnostic variables
    real(rk):: Kc0,Kc1,Kc2
    real(rk):: H_
    real(rk):: hco3,co3,co2,pCO2

    _LOOP_BEGIN_
      !state variables
      _GET_(self%id_DIC,DIC)
      !diagnostic variables
      _GET_(self%id_Hplus,H_)
      !Equilibrium constants
      _GET_(self%id_Kc0,Kc0)
      _GET_(self%id_Kc1,Kc1)
      _GET_(self%id_Kc2,Kc2)

      !Calculate all the others as a function of
      ![H+], DIC and constants
      hco3 = DIC/(1._rk+H_/Kc1+Kc2/H_)!these are in [uM]
      co3  = DIC/(1._rk+H_/Kc2+H_*H_/Kc1/Kc2)
      co2  = DIC/(1._rk+Kc1/H_+Kc1*Kc2/H_/H_)
      pco2 = co2/Kc0 ![uatm]

      _SET_DIAGNOSTIC_(self%id_pCO2,pCO2)
      _SET_DIAGNOSTIC_(self%id_CO3,co3)
    _LOOP_END_
  end subroutine do
end module
