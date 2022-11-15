subroutine excurve_TrpBH
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_pBH,  only : aij_trpBH, aijaij_trpBH
  use def_bh_parameter, only : mass_pBH, mom_pBH, spin_pBH
  use coordinate_grav_r, only : hrg
  use def_vector_x, only : hvec_xg
  use make_array_4d
  implicit none
  integer :: info
  integer :: irg, itg, ipg, i, j, k, l, m
  real(long) :: gm(3,3), fac, fac_p, fac_s, hnk_pk, hn(3), p(3), s(3)
  real(long) :: epsi(3,3,3), eps_s(3)
!
! --- Compute components of extringic curvature 
! --- whose values are assigned on the mid points. 
!
  info = 0
!
  fac = 3.0d0*dsqrt(3.0d0)*mass_pBH**2.0d0/4.0d0
  fac_p = 3.0d0/2.0d0
  fac_s = 6.0d0
  epsi(1:3,1:3,1:3) = 0.0d0
  epsi(1,2,3) =  1.0d0 ; epsi(2,3,1) =  1.0d0 ; epsi(3,1,2) =  1.0d0
  epsi(2,1,3) = -1.0d0 ; epsi(1,3,2) = -1.0d0 ; epsi(3,2,1) = -1.0d0
  gm(1:3,1:3) = 0.0d0
  gm(1,1) = 1.0d0 ; gm(2,2) = 1.0d0 ; gm(3,3) = 1.0d0
  p(1) = mom_pBH(1) ; s(1) = spin_pBH(1)
  p(2) = mom_pBH(2) ; s(2) = spin_pBH(2)
  p(3) = mom_pBH(3) ; s(3) = spin_pBH(3)
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        aijaij_trpBH(irg,itg,ipg) = 0.0d0
        hn(1) =  hvec_xg(irg,itg,ipg,1)/hrg(irg)
        hn(2) =  hvec_xg(irg,itg,ipg,2)/hrg(irg)
        hn(3) =  hvec_xg(irg,itg,ipg,3)/hrg(irg)
        hnk_pk = hn(1)*p(1) + hn(2)*p(2) + hn(3)*p(3)
        do k = 1, 3
          eps_s(k) = 0.0d0
          do l = 1, 3
            do m = 1, 3
              eps_s(k) = eps_s(k) + epsi(k,l,m)*s(l)*hn(m)
            end do
          end do
        end do
!
        do i = 1, 3
          do j = 1, 3
!
            aij_trpBH(irg,itg,ipg,i,j) = fac/hrg(irg)**3.0d0 &
            &                          * (gm(i,j) - 3.0d0*hn(i)*hn(j))&
            &                          + fac_p/hrg(irg)**2.0d0 &
            &                          *(p(i)*hn(j) + p(j)*hn(i) &
            &                          - hnk_pk*(gm(i,j) - hn(i)*hn(j))) &
            &                          + fac_s/hrg(irg)**3.0d0 &
            &                          *0.5d0*(hn(i)*eps_s(j) &
            &                                + hn(j)*eps_s(i))
!
            aijaij_trpBH(irg,itg,ipg) = aijaij_trpBH(irg,itg,ipg) &
            &                         + aij_trpBH(irg,itg,ipg,i,j)**2.0d0
          end do
        end do
        if (aijaij_trpBH(irg,itg,ipg) /= 0.0d0) info = 1
      end do
    end do
  end do
!
  if (info /= 1) write(6,*) ' ### Warning Aij_trpBH = 0 *** '
!
end subroutine excurve_TrpBH 
