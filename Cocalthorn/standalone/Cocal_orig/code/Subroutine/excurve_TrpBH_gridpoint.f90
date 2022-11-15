subroutine excurve_TrpBH_gridpoint
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_pBH,  only : aij_trpBH_grid, aijaij_trpBH_grid
  use def_bh_parameter, only : mass_pBH, mom_pBH, spin_pBH
  use coordinate_grav_r, only : rg
  use def_vector_x, only : vec_xg
  use make_array_4d
  implicit none
  integer :: info
  integer :: irg, itg, ipg, i, j, k, l, m
  real(long) :: gm(3,3), fac, fac_p, fac_s, nk_pk, n(3), p(3), s(3)
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
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 1, nrg
        aijaij_trpBH_grid(irg,itg,ipg) = 0.0d0
        n(1) = vec_xg(irg,itg,ipg,1)/rg(irg)
        n(2) = vec_xg(irg,itg,ipg,2)/rg(irg)
        n(3) = vec_xg(irg,itg,ipg,3)/rg(irg)
        nk_pk = n(1)*p(1) + n(2)*p(2) + n(3)*p(3)
        do k = 1, 3
          eps_s(k) = 0.0d0
          do l = 1, 3
            do m = 1, 3
              eps_s(k) = eps_s(k) + epsi(k,l,m)*s(l)*n(m)
            end do
          end do
        end do
!

        do i = 1, 3
          do j = 1, 3
!
            aij_trpBH_grid(irg,itg,ipg,i,j) = fac/rg(irg)**3.0d0 &
            &                               * (gm(i,j) - 3.0d0*n(i)*n(j))&
            &                               + fac_p/rg(irg)**2.0d0 &
            &                               *(p(i)*n(j) + p(j)*n(i) &
            &                               - nk_pk*(gm(i,j) - n(i)*n(j))) &
            &                               + fac_s/rg(irg)**3.0d0 &
            &                               *0.5d0*(n(i)*eps_s(j) &
            &                                     + n(j)*eps_s(i))
!
            aijaij_trpBH_grid(irg,itg,ipg) = aijaij_trpBH_grid(irg,itg,ipg) &
            &                      + aij_trpBH_grid(irg,itg,ipg,i,j)**2.0d0
          end do
        end do
        if (aijaij_trpBH_grid(irg,itg,ipg) /= 0.0d0) info = 1
      end do
    end do
  end do
!
  aij_trpBH_grid(0,0:ntg,0:npg,1:3,1:3) = 1.0d+40
  aijaij_trpBH_grid(0,0:ntg,0:npg) = 1.0d+40
!
  if (info /= 1) write(6,*) ' ### Warning Aij_trpBH = 0 *** '
!
end subroutine excurve_TrpBH_gridpoint
