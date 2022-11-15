subroutine test_source_helical_binary_mpt(impt)
  use phys_constant, only     : long, pi
  use coordinate_grav_r, only : rg
  use grid_parameter, only    : nrg, ntg, npg
  use def_binary_parameter, only    : sepa, mass_ratio
  use def_matter, only        : emdg, rs
  use def_matter_parameter, only : ome
  use def_bh_parameter, only : ome_bh
  use def_vector_x, only      :vec_xg
  implicit none
  integer     ::   irg, itg, ipg, impt
  real(long)  ::   zfac, small = 1.0d-15
  real(long), parameter ::   a = 1.0d0, q = 1.0d0, sigma = 0.5d0
  real(long)  ::   xx, yy, zz, rplusR, rminusR
!
!
!  ome = 0.1d0
!  ome_bh = ome
  ome = ome_bh
  call calc_vector_x_grav(1)
  if(impt.eq.1) then
    do ipg = 0, npg
      do itg = 0, ntg
      rs(itg,ipg) = 1.0d0
        do irg = 0, nrg
          xx = vec_xg(irg,itg,ipg,1)
          yy = vec_xg(irg,itg,ipg,2)
          zz = vec_xg(irg,itg,ipg,3)
          rplusR  = (xx - sepa)**2 + yy**2 + zz**2
          rminusR = xx**2 + yy**2 + zz**2
          emdg(irg,itg,ipg) = q/((sqrt(2.0d0*pi)*sigma)**3) &
          &                 * ( 1.0d0     *exp(-rplusR /(2.0d0*sigma**2)) &
          &                   + mass_ratio*exp(-rminusR/(2.0d0*sigma**2)) )
        end do
      end do
    end do
  end if
  if(impt.eq.2) then
    do ipg = 0, npg
      do itg = 0, ntg
      rs(itg,ipg) = 1.0d0
        do irg = 0, nrg
          xx = vec_xg(irg,itg,ipg,1)
          yy = vec_xg(irg,itg,ipg,2)
          zz = vec_xg(irg,itg,ipg,3)
          rplusR  =  xx**2 + yy**2 + zz**2
          rminusR = (xx - sepa)**2 + yy**2 + zz**2
          emdg(irg,itg,ipg) = q/((sqrt(2.0d0*pi)*sigma)**3) &
          &                 * ( 1.0d0     *exp(-rplusR /(2.0d0*sigma**2)) &
          &                   + mass_ratio*exp(-rminusR/(2.0d0*sigma**2)) )
        end do
      end do
    end do
  end if
!
!
end subroutine test_source_helical_binary_mpt
