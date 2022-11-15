subroutine test_source_helical
  use phys_constant, only     : long, pi
  use coordinate_grav_r, only : rg
  use grid_parameter, only    : nrg, ntg, npg
  use def_matter, only        : emdg, rs
  use def_matter_parameter, only : ome
  use def_vector_x, only      :vec_xg
  implicit none
  integer     ::   irg, itg, ipg
  real(long)  ::   zfac, small = 1.0d-15
  real(long), parameter ::   a = 1.0d0, q = 1.0d0, sigma = 0.5d0
  real(long)  ::   xx, yy, zz, rplusR, rminusR
!
!
  ome = 0.3d0
  call calc_vector_x_grav(1)
  do ipg = 0, npg
    do itg = 0, ntg
    rs(itg,ipg) = 1.0d0
      do irg = 0, nrg
        xx = vec_xg(irg,itg,ipg,1)
        yy = vec_xg(irg,itg,ipg,2)
        zz = vec_xg(irg,itg,ipg,3)
        rplusR  = (xx + a)**2 + yy**2 + zz**2
        rminusR = (xx - a)**2 + yy**2 + zz**2
        emdg(irg,itg,ipg) = q/((sqrt(2.0d0*pi)*sigma)**3) &
        &                 * ( exp(-rplusR /(2.0d0*sigma**2)) &
        &                   + exp(-rminusR/(2.0d0*sigma**2)) )
      end do
    end do
  end do
!
!
end subroutine test_source_helical
