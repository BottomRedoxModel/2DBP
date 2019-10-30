!-----------------------------------------------------------------------
! brom_injection is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the FABM distribution.
!-----------------------------------------------------------------------
! Original author(s): Evgeniy Yakushev, Anfisa Berezina
!-----------------------------------------------------------------------
#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE:
!
! !INTERFACE:
   module fabm_niva_brom_injection
!
! !DESCRIPTION:
!
! brom_injection parameterizes chemical/ biogeochemical transfomration of an injected matter.
! brom_injection consists of 1 state variables ( in nitrogen units):
! - waste - is a compound that can be can be decayed and /or oxidized to inorganic nutrients 
!   with a strong oxidizing agent under acidic condition (oxygen).
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !REVISION HISTORY:!
!  Original author(s): Evgeniy Yakushev, Anfisa Berezina
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_injection
!     Variable identifiers
      type (type_state_variable_id)        :: id_oxy, id_nut, id_pom,  id_dom, id_waste
      type (type_dependency_id)            :: id_par, id_temp, id_salt
!      type (type_horizontal_dependency_id) :: id_windspeed
      type (type_diagnostic_variable_id)   :: id_waste_decay,id_waste_miner
!     Model parameters
      !---Organic matter mineralization---- !
       real(rk) :: r_waste_decay,r_waste_miner !,r_waste_dom_oxy, r_waste_dom_nut, 
       real(rk) :: Wwaste, Bu
   contains
      procedure :: initialize
      procedure :: do
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the OXYDEP-COD model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
!
! !INPUT PARAMETERS:
   class (type_niva_brom_injection), intent(inout), target :: self
   integer,                          intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
       real(rk),parameter :: d_per_s = 1.0_rk/86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
  ! waste
   call self%get_parameter(self%r_waste_decay, 'r_waste_decay', '1/d', 'Specific rate of waste decay', default=0.10_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_waste_miner, 'r_waste_miner', '1/d', 'Specific rate of waste mineralization',default=0.10_rk,scale_factor=d_per_s)
!   call self%get_parameter(self%beta_da,        'beta_da',         'nd', 'Coefficient for dependence of mineralization on t ', default=20._rk)
!   call self%get_parameter(self%Tda,            'Tda',             'nd', 'Coefficient for dependence of mineralization on t ', default=13._rk)
   call self%get_parameter(self%Wwaste,         'Wwaste',          'm/s', 'vertical velocity of waste (<0 for sinking)', default=-1.0_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Bu,             'Bu',             'nd',  'Burial coeficient for lower boundary',      default=0.25_rk)
   ! Register state variables
   call self%register_state_variable(self%id_waste,'waste','mmol/m**3','waste   compound ', 0.0_rk, minimum=0.0_rk, vertical_movement=self%Wwaste)
   ! Register link to external variables
   call self%register_state_dependency(self%id_oxy,'Oxy','mmol/m**3','OXY')
!   call self%register_state_dependency(self%id_oxy,'NUT','mmol/m**3','NUT')
   call self%register_state_dependency(self%id_pom,'POM','mmol/m**3','POM')
 !  call self%register_state_dependency(self%id_oxy,'DOM','mmol/m**3','DOM')
   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_waste_decay,'waste_decay','mmol/m**3/d',  'waste_decay,  Decay of waste',&
                    output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_waste_miner,'waste_miner','mmol/m**3/d',  'waste_miner,  Mineralization of waste with oxygen',&
                    output=output_time_step_integrated)

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
!   call self%register_dependency(self%id_windspeed,standard_variables%wind_speed)
!   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
!
! !INPUT PARAMETERS:
   class (type_niva_brom_injection),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
   real(rk)                   :: waste,oxy,nut,dom,pom,t,iopt

   real(rk) :: doxy,dwaste,dpom
 ! Rates of biogeochemical processes
   real(rk) :: waste_decay            ! decay mineralization (1/d)
   real(rk) :: waste_miner            ! oxic mineralization of waste (1/d)
!   real(rk) :: waste_decay_denitr        ! suboxic mineralization of waste (denitrification) (1/d)
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_waste,waste)
!   _GET_(self%id_nut,nut)
   _GET_(self%id_pom,pom)
!   _GET_(self%id_dom,dom)

   ! Retrieve current environmental conditions.
!   _GET_(self%id_par,Iz)              ! local photosynthetically active radiation
   _GET_(self%id_temp,t)              ! temperature

!--------------------------------------------------------------
! Oxic mineralization of waste depends on T
   waste_decay = self%r_waste_decay*waste
!   waste_decay = self%r_waste_nut_oxy*(1.+self%beta_da*yy(self%tda,t))*waste
! Suboxic mineralization depends on T,O2,NO3/NO2
   !waste_decay_denitr = self%r_pom_nut_nut*(1.+self%beta_da*yy(self%tda,t)) &
   !                        * (0.5-0.5*tanh(self%O2LimC-oxy)) &
   !                        * (1-tanh(1.-nut))*pom
! Mineralization of OM, ammonification and growth of NUT


! Now we can summarize processes and write state variables sink/sources:
!--------------------------------------------------------------
! OXY
!--------------------------------------------------------------
! Changes of OXY due to OM production and decay!
!   doxy  = -2.0_rk*waste_decay_ox
   dwaste = -waste_decay
!   dnut = 
   dpom   =  waste_decay
!--------------------------------------------------------------



!derivatives for FABM
   _SET_ODE_(self%id_waste,dwaste)
   _SET_ODE_(self%id_pom,dpom)
!   _SET_ODE_(self%id_nut,dnut)
!   _SET_ODE_(self%id_dom,ddom)
!   _SET_ODE_(self%id_oxy, doxy)

   
   ! Export diagnostic variables
_SET_DIAGNOSTIC_(self%id_waste_decay,waste_decay)
   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !DESCRIPTION:
!

! !INPUT PARAMETERS:
   class (type_niva_brom_injection),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk)                   :: waste

   _HORIZONTAL_LOOP_BEGIN_
   _GET_(self%id_waste,waste)

   ! BURYING into the sediments, mmol/m2/s (sinking rates "Wxxx" are in m/s and positive upward)
   _SET_BOTTOM_EXCHANGE_(self%id_waste,self%Bu*self%Wwaste*waste)

   _HORIZONTAL_LOOP_END_

   end subroutine
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Saturation function squared
!
! !INTERFACE:
   real(rk) function yy(a,x)
!
! !DESCRIPTION:
! This is a squared Michaelis-Menten type of limiter:
! \begin{equation}\label{Y}
! Y(x_w,x) = \frac{x^2}{x_w^2+x^2}.
! \end{equation}
!
! !IN2PUT PARAMETERS:
   real(rk), intent(in)                :: a,x
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   yy=x**2/(a**2+x**2)

   end function yy
!EOC

   end module fabm_niva_brom_injection

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
