subroutine gauge_Coulomb
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r, only : rg
  use def_matter,  only : rs
  use def_emfield, only : vaxd, vayd, vazd
  use interface_poisson_solver
  use make_array_3d
  use interface_grgrad_4th_gridpoint
  use interface_grgrad_midpoint_type0
  implicit none
!
  real(8), pointer :: sou(:,:,:), gauge(:,:,:)
  real(8) :: dvaxdx, dvaxdy, dvaxdz, dvaydx, dvaydy, dvaydz, &
  &          dvazdx, dvazdy, dvazdz
  real(8) :: dgadx, dgady, dgadz
  integer :: irg, itg, ipg
!
  call alloc_array3d(sou,0,nrg,0,ntg,0,npg)
  call alloc_array3d(gauge,0,nrg,0,ntg,0,npg)
!
! --- Impose extended Dirac gauge condition.
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        call grgrad_midpoint_type0(vaxd,dvaxdx,dvaxdy,dvaxdz,irg,itg,ipg)
        call grgrad_midpoint_type0(vayd,dvaydx,dvaydy,dvaydz,irg,itg,ipg)
        call grgrad_midpoint_type0(vazd,dvazdx,dvazdy,dvazdz,irg,itg,ipg)
        sou(irg,itg,ipg) = -(dvaxdx + dvaydy + dvazdz)
!
      end do
    end do
  end do
!
  call poisson_solver(sou,gauge)
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
!
        call grgrad_4th_gridpoint(gauge,dgadx,dgady,dgadz,irg,itg,ipg)
!
!test
!!        if (rg(irg).gt.rs(itg,ipg)) then
!test
        vaxd(irg,itg,ipg) = vaxd(irg,itg,ipg) + dgadx
        vayd(irg,itg,ipg) = vayd(irg,itg,ipg) + dgady
        vazd(irg,itg,ipg) = vazd(irg,itg,ipg) + dgadz
!test
!!        end if
!test
!
      end do
    end do
  end do
!
  deallocate(sou)
  deallocate(gauge)
end subroutine gauge_Coulomb
