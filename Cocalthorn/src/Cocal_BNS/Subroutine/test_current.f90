subroutine test_current
  use phys_constant, only : long, pi
  use grid_parameter, only : nrg, ntg, npg
  use def_emfield, only : jtu, jxu, jyu, jzu
  use def_vector_x, only : vec_xg
  use def_vector_phi, only : vec_phig
  implicit none
  real(long) :: ome_cur, dis_cur, sigma, charge
  real(long) :: xx, yy, zz, phix, phiy, phiz, rminusR, factor
  integer :: irg, itg, ipg
!
! --- Source for Maxwell eq normal component
! --  current term
!
  ome_cur = 0.3d0
  dis_cur = 0.5d0
  sigma   = 0.3d0
  charge  = 0.003d0
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
!
        xx = vec_xg(irg,itg,ipg,1)
        yy = vec_xg(irg,itg,ipg,2)
        zz = vec_xg(irg,itg,ipg,3)
        phix = vec_phig(irg,itg,ipg,1)
        phiy = vec_phig(irg,itg,ipg,2)
        phiz = vec_phig(irg,itg,ipg,3)
!
        rminusR  = (sqrt(xx**2+yy**2) - dis_cur)**2 + zz**2
        factor = charge/(sqrt(2.0d0*pi)*sigma)**2 &
        &      * exp(-rminusR /(2.0d0*sigma**2))
! 
        jtu(irg,itg,ipg) = factor
        jxu(irg,itg,ipg) = factor*ome_cur*phix
        jyu(irg,itg,ipg) = factor*ome_cur*phiy
        jzu(irg,itg,ipg) = 0.0d0
!
      end do
    end do
  end do
!
end subroutine test_current

