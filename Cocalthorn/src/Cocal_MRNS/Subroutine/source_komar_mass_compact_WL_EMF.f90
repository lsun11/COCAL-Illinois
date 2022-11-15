subroutine source_komar_mass_compact_WL_EMF(soug)
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg
  use def_metric, only : psi, alph, bvxu, bvyu, bvzu
  use def_SEM_tensor_EMF, only : rhoH_EMF, jmd_EMF, trsm_EMF
  use interface_interpo_linear_type0
  implicit none
  real(long), pointer :: soug(:,:,:)
  real(long) :: alphw, psiw, rhoHw, esseS
  real(long) :: rjjx, rjjy, rjjz, rjjbeta
  real(long) :: bvxuw, bvyuw, bvzuw
  integer    :: irg, itg, ipg
!
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        call interpo_linear_type0(psiw,psi,irg,itg,ipg)
        call interpo_linear_type0(alphw,alph,irg,itg,ipg)
        call interpo_linear_type0(bvxuw,bvxu,irg,itg,ipg)
        call interpo_linear_type0(bvyuw,bvyu,irg,itg,ipg)
        call interpo_linear_type0(bvzuw,bvzu,irg,itg,ipg)
!
        rhoHw = rhoH_EMF(irg,itg,ipg)
        esseS = trsm_EMF(irg,itg,ipg)
        rjjx  =  jmd_EMF(irg,itg,ipg,1)
        rjjy  =  jmd_EMF(irg,itg,ipg,2)
        rjjz  =  jmd_EMF(irg,itg,ipg,3)
        rjjbeta = rjjx*bvxuw + rjjy*bvyuw + rjjz*bvzuw
!
        soug(irg,itg,ipg) = (alphw*(esseS+rhoHw) - 2.0d0*rjjbeta)*psiw**6
      end do
    end do
  end do
!
end subroutine source_komar_mass_compact_WL_EMF
