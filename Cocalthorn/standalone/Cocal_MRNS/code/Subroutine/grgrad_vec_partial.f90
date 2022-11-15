subroutine grgrad_vec_partial(fnx,fny,fnz,pdfnx,pdfny,pdfnz)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use interface_grgrad_midpoint
  use make_array_3d
  implicit none
!
  real(long), pointer :: fnx(:,:,:), fny(:,:,:), fnz(:,:,:)
  real(long), pointer :: pdfnx(:,:,:,:), pdfny(:,:,:,:), pdfnz(:,:,:,:)
  real(long), pointer :: dfdx(:,:,:), dfdy(:,:,:), dfdz(:,:,:)
!
  call alloc_array3d(dfdx,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfdy,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfdz,1,nrg,1,ntg,1,npg)
!
  call grgrad_midpoint(fnx,dfdx,dfdy,dfdz)
  pdfnx(1:nrg,1:ntg,1:npg,1) = dfdx(1:nrg,1:ntg,1:npg)
  pdfnx(1:nrg,1:ntg,1:npg,2) = dfdy(1:nrg,1:ntg,1:npg)
  pdfnx(1:nrg,1:ntg,1:npg,3) = dfdz(1:nrg,1:ntg,1:npg)
!
  call grgrad_midpoint(fny,dfdx,dfdy,dfdz)
  pdfny(1:nrg,1:ntg,1:npg,1) = dfdx(1:nrg,1:ntg,1:npg)
  pdfny(1:nrg,1:ntg,1:npg,2) = dfdy(1:nrg,1:ntg,1:npg)
  pdfny(1:nrg,1:ntg,1:npg,3) = dfdz(1:nrg,1:ntg,1:npg)
!
  call grgrad_midpoint(fnz,dfdx,dfdy,dfdz)
  pdfnz(1:nrg,1:ntg,1:npg,1) = dfdx(1:nrg,1:ntg,1:npg)
  pdfnz(1:nrg,1:ntg,1:npg,2) = dfdy(1:nrg,1:ntg,1:npg)
  pdfnz(1:nrg,1:ntg,1:npg,3) = dfdz(1:nrg,1:ntg,1:npg)
!
  deallocate(dfdx)
  deallocate(dfdy)
  deallocate(dfdz)
end subroutine grgrad_vec_partial
