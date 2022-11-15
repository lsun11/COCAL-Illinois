subroutine source_ang_mom_EMF(soug)
  use phys_constant, only :  long, pi
  use grid_parameter, only :  nrg, ntg, npg
  use def_SEM_tensor_EMF, only : jmd_EMF
  use def_metric, only : psi
  use interface_interpo_linear_type0
  use def_vector_phi, only : hvec_phig
  implicit none
  real(long), pointer ::  soug(:,:,:)
  real(long)  ::   rjjx, rjjy, rjjz, rjjphi
  real(long)  ::   vphig(1:3), psigc, fac8pi=8.0d0*pi
  integer     ::   ipg, itg, irg, ia, ib
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        call interpo_linear_type0(psigc,psi,irg,itg,ipg)
        vphig(1) = hvec_phig(irg,itg,ipg,1)
        vphig(2) = hvec_phig(irg,itg,ipg,2)
        vphig(3) = hvec_phig(irg,itg,ipg,3)
        rjjx  =  jmd_EMF(irg,itg,ipg,1)
        rjjy  =  jmd_EMF(irg,itg,ipg,2)
        rjjz  =  jmd_EMF(irg,itg,ipg,3)
        rjjphi = rjjx*vphig(1) + rjjy*vphig(2) + rjjz*vphig(3)
!
        soug(irg,itg,ipg) = fac8pi*rjjphi*psigc**6
!
      end do
    end do
  end do
!
end subroutine source_ang_mom_EMF
