subroutine source_komar_mass_peos_EMF(soug)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrg, ntg, npg
  use def_metric, only  :  psi, alph
  use def_SEM_tensor_EMF, only : rhoH_EMF, trsm_EMF
  use interface_interpo_linear_type0
  implicit none
  real(long),pointer ::   soug(:,:,:)
  real(long)  ::   alphw, psiw, rhoHw, esseS
  integer     ::   irg,itg,ipg
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        call interpo_linear_type0(psiw,psi,irg,itg,ipg)
        call interpo_linear_type0(alphw,alph,irg,itg,ipg)
!
        rhoHw = rhoH_EMF(irg,itg,ipg)
        esseS = trsm_EMF(irg,itg,ipg)
! 
        soug(irg,itg,ipg) = 4.0d0*pi*alphw*psiw**6*(esseS+rhoHw)
!
      end do
    end do
  end do
!
end subroutine source_komar_mass_peos_EMF
