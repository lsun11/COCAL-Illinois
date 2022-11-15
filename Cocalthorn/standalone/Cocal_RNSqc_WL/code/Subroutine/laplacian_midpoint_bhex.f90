subroutine laplacian_midpoint_bhex(fnc,lapfnc)
  use phys_constant, only :  long
  use grid_parameter, only : nrg, ntg, npg
  use interface_interpo_linear_type0
  use interface_grgrad_midpoint
  use interface_grgrad1g_midpoint
!  use interface_dadbscalar_type0
  use interface_dadbscalar_type3_bhex
  use make_array_2d
  implicit none
  real(long), pointer :: fnc(:,:,:), lapfnc(:,:,:)
  real(long), pointer :: d2fnc(:,:)
  integer :: ipg, itg, irg
!
  call alloc_array2d(d2fnc,1,3,1,3)

  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
!        call dadbscalar_type0(fnc,d2fnc,irg,itg,ipg)
        call dadbscalar_type3_bhex(fnc,d2fnc,irg,itg,ipg)
!
        lapfnc(irg,itg,ipg) = d2fnc(1,1) + d2fnc(2,2) + d2fnc(3,3)
      end do
    end do
  end do
!
  deallocate(d2fnc)
!
end subroutine laplacian_midpoint_bhex
