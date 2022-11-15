subroutine calc_4velocity_ut_irrot
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf, ntfeq, npfyzp
  use def_matter, only : emd, utf, utg, lambda 
  use def_velocity_potential, only : grad_vpot
  use def_vector_phi, only : vec_phif
  use def_metric, only :   alph, bvxu, bvyu, bvzu
  use make_array_3d
  use interface_interpo_gr2fl
  use interface_interpo_fl2gr
  implicit none
  real(long) :: emdfc, rhofc, prefc, hhfc, ene
  real(long) :: bvxfc, bvyfc, bvzfc, vphix, vphiy, vphiz
  real(long) :: dxvp, dyvp, dzvp, lam, lamfc
  integer :: ir, it, ip
!
  real(long), pointer :: alphf(:,:,:)
  real(long), pointer :: bvxuf(:,:,:), bvyuf(:,:,:), bvzuf(:,:,:)
!
  call alloc_array3d(alphf, 0, nrf, 0, ntf, 0, npf) 
  call alloc_array3d(bvxuf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(bvyuf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(bvzuf, 0, nrf, 0, ntf, 0, npf)
!
  call interpo_gr2fl(alph,alphf)
  call interpo_gr2fl(bvxu,bvxuf)
  call interpo_gr2fl(bvyu,bvyuf)
  call interpo_gr2fl(bvzu,bvzuf)
!
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
!
        if (ir.eq.0) then 
          bvyfc = bvyuf(ir,ntfeq,npfyzp)
          vphiy = vec_phif(ir,ntfeq,npfyzp,2)
          dyvp  = grad_vpot(ir,ntfeq,npfyzp,2)
          lam   = ber + (bvyfc + ome*vphiy)*dyvp
        else 
          bvxfc = bvxuf(ir,it,ip)
          bvyfc = bvyuf(ir,it,ip)
          bvzfc = bvzuf(ir,it,ip)
          vphix = vec_phif(ir,it,ip,1)
          vphiy = vec_phif(ir,it,ip,2)
          vphiz = vec_phif(ir,it,ip,3)
          dxvp  = grad_vpot(ir,it,ip,1)
          dyvp  = grad_vpot(ir,it,ip,2)
          dzvp  = grad_vpot(ir,it,ip,3)
          lam   = ber + (bvxfc + ome*vphix)*dxvp &
          &           + (bvyfc + ome*vphiy)*dyvp &
          &           + (bvzfc + ome*vphiz)*dzvp
        end if
!
        lamfc = lam
        alpfc = alphf(ir,it,ip)
        emdfc = emd(ir,it,ip)
        if (emdfc <= 1.0d-15) emdfc = 1.0d-15
        call peos_q2hprho(emdfc, hhfc, prefc, rhofc, ene)
!
        lambda(ir,it,ip) = lamfc
        utf(ir,it,ip) = lamfc/(alpfc**2*hhfc)
!
      end do
    end do
  end do
  call interpo_fl2gr(utf,utg)
!
  deallocate(alphf)
  deallocate(bvxuf)
  deallocate(bvyuf)
  deallocate(bvzuf)
end subroutine calc_4velocity_ut_irrot
