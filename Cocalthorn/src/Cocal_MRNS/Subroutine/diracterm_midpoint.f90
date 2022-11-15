subroutine diracterm_midpoint
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_metric_hij_dirac, only : dagmabu
  use interface_grgrad1g_midpoint
  implicit none
!
  real(8) :: grad1(3), dhd(3,3,3), dhu(3,3,3)
  integer :: ipg, itg, irg
!
! --- Compute Dirac term partial_a gamma^ab
!     whose value is assigned on the mid points. 
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
! --- Derivatives of h_ab.
!
        call grgrad1g_midpoint(hxxu,grad1,irg,itg,ipg)
        dhu(1,1,1:3) = grad1(1:3)
        call grgrad1g_midpoint(hxyu,grad1,irg,itg,ipg)
        dhu(1,2,1:3) = grad1(1:3)
        dhu(2,1,1:3) = grad1(1:3)
        call grgrad1g_midpoint(hxzu,grad1,irg,itg,ipg)
        dhu(1,3,1:3) = grad1(1:3)
        dhu(3,1,1:3) = grad1(1:3)
        call grgrad1g_midpoint(hyyu,grad1,irg,itg,ipg)
        dhu(2,2,1:3) = grad1(1:3)
        call grgrad1g_midpoint(hyzu,grad1,irg,itg,ipg)
        dhu(2,3,1:3) = grad1(1:3)
        dhu(3,2,1:3) = grad1(1:3)
        call grgrad1g_midpoint(hzzu,grad1,irg,itg,ipg)
        dhu(3,3,1:3) = grad1(1:3)
!
        dagmabu(irg,itg,ipg,1) = dhu(1,1,1) + dhu(1,2,2) + dhu(1,3,3)
        dagmabu(irg,itg,ipg,2) = dhu(2,1,1) + dhu(2,2,2) + dhu(2,3,3)
        dagmabu(irg,itg,ipg,3) = dhu(3,1,1) + dhu(3,2,2) + dhu(3,3,3)
!
      end do
    end do
  end do
!
end subroutine diracterm_midpoint
