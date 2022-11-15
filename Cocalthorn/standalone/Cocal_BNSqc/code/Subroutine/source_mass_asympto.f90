subroutine source_mass_asympto(fnc,sousf,irg)
  use phys_constant, only  : long
  use grid_parameter, only : ntg, npg
  use make_array_2d
  use interface_interpo_linear_type0_2Dsurf
  use interface_grdr_gridpoint_type0_nosym
  implicit none
  real(long), pointer :: fnc(:,:,:),dfnc(:,:), sousf(:,:)
  integer    :: irg, itg, ipg
  real(long) :: deriv, val 
!
  call alloc_array2d(dfnc, 0,ntg, 0,npg)
!
  do ipg = 0, npg
    do itg = 0, ntg
      call grdr_gridpoint_type0_nosym(fnc,deriv,irg,itg,ipg)
      dfnc(itg,ipg) = deriv
    end do
  end do
!
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(val,dfnc,itg,ipg)
      sousf(itg,ipg) = val
    end do
  end do
!
  deallocate(dfnc)
end subroutine source_mass_asympto
