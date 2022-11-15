subroutine source_adm_mass_WL_EMF(soug)
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg
  use def_metric, only     : psi
  use def_SEM_tensor_EMF, only : rhoH_EMF
  use interface_interpo_linear_type0
  implicit none
  real(long), pointer :: soug(:,:,:)
  integer     ::   irg, itg, ipg
  real(long)  ::   psiw, rhoHw
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        call interpo_linear_type0(psiw,psi,irg,itg,ipg)
        rhoHw = rhoH_EMF(irg,itg,ipg)
        soug(irg,itg,ipg) = 2.0d0*pi*psiw**5*rhoHw
      end do
    end do
  end do
!
end subroutine source_adm_mass_WL_EMF
