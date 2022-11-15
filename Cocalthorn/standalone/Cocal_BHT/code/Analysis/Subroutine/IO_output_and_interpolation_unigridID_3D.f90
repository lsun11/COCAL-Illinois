subroutine IO_output_and_interpolation_unigridID_3D
  use phys_constant, only : long
  use grid_parameter_cartesian, only : nx, ny, nz, nxhalf
  use coordinate_grav_xyz
  use def_metric, only : psi, alph, bvxd, bvyd, bvzd
  use def_matter, only : emd, rhof, utf, uxf, uyf, uzf
  use interface_modules_cartesian
  use make_array_3d
  implicit none
  real(long), pointer :: fncca(:,:,:)
  integer :: nvar = 11, ix, iy, iz, igrid, ii
  character(len=5) :: char_file(11) = (/ 'Psi__', 'alpha', &
  &  'betax', 'betay', 'betaz', 'rho__', 'ut___', 'ux___', 'uy___', 'uz___', 'emden' /)
!
  call  alloc_array3d(fncca,1,nx,1,ny,1,nz)
!
  igrid = mod(nx,2)
  do ii = 1, nvar
    if (ii.eq.1) call interpolation_metric(psi,fncca)
    if (ii.eq.2) call interpolation_metric(alph,fncca)
    if (ii.eq.3) call interpolation_metric(bvxd,fncca)
    if (ii.eq.4) call interpolation_metric(bvyd,fncca)
    if (ii.eq.5) call interpolation_metric(bvzd,fncca)
!
    if (ii.eq.6) call interpolation_matter(rhof,fncca)
    if (ii.eq.7) call interpolation_matter(utf,fncca)
    if (ii.eq.8) call interpolation_matter(uxf,fncca)
    if (ii.eq.9) call interpolation_matter(uyf,fncca)
    if (ii.eq.10)call interpolation_matter(uzf,fncca)
    if (ii.eq.11) call interpolation_matter(emd,fncca)
!
    open(15,file='./ID/'//char_file(ii)//'.dat',status='unknown')
!    do iz = 1, nxhalf + 3 + igrid
    do iz = nxhalf - 3 + igrid, nx
      do iy = 1, nx
        do ix = 1, nx
          write(15,'(1p,e20.12)') fncca(ix,iy,iz)
        end do
      end do
    end do
    close(15)
!
  end do
!
  open(15,file='CartesianCoordinates.dat',status='unknown')
  igrid = mod(nx,2)
!  do iz = 1, nxhalf + 3 + igrid
  do iz = nxhalf - 3 + igrid, nx
    do iy = 1, nx
      do ix = 1, nx
        write(15,'(1p,3e20.12)') x(ix),y(iy),z(iz)
      end do
    end do
  end do
  close(15)
!
  deallocate(fncca)
end subroutine IO_output_and_interpolation_unigridID_3D
