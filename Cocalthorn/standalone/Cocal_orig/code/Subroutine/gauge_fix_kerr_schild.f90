subroutine gauge_fix_kerr_schild
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric, only : trk
  use def_metric_excurve_grid, only : trk_grid
  use def_transverse_part, only : Ftvx,      Ftvy,      Ftvz, &
  &                               Ftvx_grid, Ftvy_grid, Ftvz_grid
  use def_vector_x, only : vec_xg, hvec_xg
  implicit none
  real(long) :: traceK, dgamma(3)
  real(long) :: x, y, z
  integer :: irg, itg, ipg
!
! --- Impose Kerr-Schild gauge condition.
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        x = hvec_xg(irg,itg,ipg,1)
        y = hvec_xg(irg,itg,ipg,2)
        z = hvec_xg(irg,itg,ipg,3)
        call kerr_schild_traceK(x,y,z,traceK)
        call kerr_schild_transverse_part(x,y,z,dgamma)
        trk( irg,itg,ipg) = traceK
        Ftvx(irg,itg,ipg) = dgamma(1)
        Ftvy(irg,itg,ipg) = dgamma(2)
        Ftvz(irg,itg,ipg) = dgamma(3)
!
      end do
    end do
  end do
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
!
        x = vec_xg(irg,itg,ipg,1)
        y = vec_xg(irg,itg,ipg,2)
        z = vec_xg(irg,itg,ipg,3)
        call kerr_schild_traceK(x,y,z,traceK)
        call kerr_schild_transverse_part(x,y,z,dgamma)
        trk_grid( irg,itg,ipg) = traceK
        Ftvx_grid(irg,itg,ipg) = dgamma(1)
        Ftvy_grid(irg,itg,ipg) = dgamma(2)
        Ftvz_grid(irg,itg,ipg) = dgamma(3)
!
      end do
    end do
  end do
!
end subroutine gauge_fix_kerr_schild
