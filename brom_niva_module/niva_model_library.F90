module niva_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: niva_model_factory

contains

   subroutine create(self,name,model)

      use fabm_niva_brom_acetate
      use fabm_niva_brom_bact
      use fabm_niva_brom_bio
      use fabm_niva_brom_calcium
      use fabm_niva_brom_carbon
      use fabm_niva_brom_cod
      use fabm_niva_brom_eq_constants
      use fabm_niva_brom_fe
      use fabm_niva_brom_hg
      use fabm_niva_brom_main_nutrients
      use fabm_niva_brom_manganese
      use fabm_niva_brom_methane
      use fabm_niva_brom_ni
      use fabm_niva_brom_nitrogen
      use fabm_niva_brom_partitioning
      use fabm_niva_brom_pH
      use fabm_niva_brom_salt
      use fabm_niva_brom_silicon
      use fabm_niva_brom_sulfur
      use fabm_niva_brom_ba
      use fabm_niva_oxydep
      use fabm_niva_brom_injection
      ! Add new NIVA models here

      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('brom_acetate');       allocate(type_niva_brom_acetate::model)
         case ('brom_bact');          allocate(type_niva_brom_bact::model)
         case ('brom_bio');           allocate(type_niva_brom_bio::model)
         case ('brom_calcium');       allocate(type_niva_brom_calcium::model)
         case ('brom_carbon');        allocate(type_niva_brom_carbon::model)
         case ('brom_cod');           allocate(type_niva_brom_cod::model)
         case ('brom_eq_constants');  allocate(type_niva_brom_eq_constants::model)
         case ('brom_fe');            allocate(type_niva_brom_fe::model)
         case ('brom_hg');            allocate(type_niva_brom_hg::model)
         case ('brom_main_nutrients');allocate(type_niva_brom_main_nutrients::model)
         case ('brom_manganese');     allocate(type_niva_brom_manganese::model)
         case ('brom_methane');       allocate(type_niva_brom_methane::model)
         case ('brom_ni');            allocate(type_niva_brom_ni::model)
         case ('brom_nitrogen');      allocate(type_niva_brom_nitrogen::model)
         case ('brom_partitioning');  allocate(type_niva_brom_partitioning::model)
         case ('brom_pH');            allocate(type_niva_brom_pH::model)
         case ('brom_salt');          allocate(type_niva_brom_salt::model)
         case ('brom_silicon');       allocate(type_niva_brom_silicon::model)
         case ('brom_sulfur');        allocate(type_niva_brom_sulfur::model)
         case ('brom_ba');            allocate(type_niva_brom_ba::model)
         case ('oxydep');             allocate(type_niva_oxydep::model)
         case ('brom_injection');     allocate(type_niva_brom_injection::model)
         ! Add new NIVA models here
      end select

   end subroutine

end module
