subroutine lentim_conversion(i,inval,outval)
  use phys_constant, only : long,nmpt, g, c, solmas
  use def_matter_parameter_mpt
  implicit none
  real(long), intent(in)  :: inval
  real(long), intent(out) :: outval
  real(long) :: R0
  integer    :: i
  real(long) :: MM = solmas, LL = g*solmas/c**2, TT = g*solmas/c**3

  R0 = def_matter_param_real_(5,1)

! cocal units to km
  if(i .eq. 1) outval = inval*LL*R0*1.0d-5

! km to cocal units
  if(i .eq. 2) outval = inval/(LL*R0*1.0d-5)

! cocal units to rad/sec
  if(i .eq. 3) outval = inval/R0/TT

! rad/sec to cocal units
  if(i .eq. 4) outval = inval*R0*TT
!
end subroutine lentim_conversion
