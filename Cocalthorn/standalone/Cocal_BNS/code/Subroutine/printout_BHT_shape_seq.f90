subroutine printout_BHT_shape_seq(iseq)
  use phys_constant, only : long, pi
  use grid_parameter, only : nrgin, nrg, ntg, npg, ntgeq, ntgxy, &
  &                          npgxzp, npgxzm, npgyzp, npgyzm
  use def_matter, only  :   emdg
  use def_matter_parameter
  use def_formulation
  use def_bht_parameter
  use coordinate_grav_r, only : rg
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi, only : sinphig, cosphig
  implicit none
  integer :: iseq
  real(long) :: r1,x1,y1,z1,rin,rout,x2,y2,z2
  real(long) :: emd0, emdgc, ray_width
  integer    :: ipg, itg, irg, iflagin

!
  if (iseq.eq.1) then
    open(20,file='bhtinterior_xz.dat',status='unknown')
    open(21,file='bhtshape_xz.dat',status='unknown')
!    open(23,file='bhtshape_seq_yz.dat',status='unknown')
!    open(24,file='bhtshape_NS_gnuplot.dat',status='unknown')
  else
    open(20,file='bhtinterior_xz.dat',status='old', position="append")
    open(21,file='bhtshape_xz.dat',status='old', position="append")
!    open(23,file='bhtshape_seq_yz.dat',status='old', position="append")
!    open(24,file='bhtshape_NS_gnuplot.dat',status='old', position="append")
  end if
!
  write(20,'(1p,2e23.15)')  rg(nrgin),0.0d0
  write(21,'(1p,2e23.15)')  rg(nrgin),0.0d0
  ipg=npgxzp
  do irg = nrgin+1, nrg
    emd0 = emdg(irg,0,ipg)*10.0d0
    iflagin=0

    do itg = 0, ntgeq
      emdgc = emdg(irg,itg,ipg)
      x1 = rg(irg)*sinthg(itg)*cosphig(ipg)
      y1 = rg(irg)*sinthg(itg)*sinphig(ipg)
      z1 = rg(irg)*costhg(itg)
!      write(6,'(a3,2i5,1p,3e23.15)')  "",irg,itg, x1,z1, emdgc
      if (emdgc > emd0) then
!        x2 = rg(irg-1)*sinthg(itg)*cosphig(ipg)
!        y2 = rg(irg-1)*sinthg(itg)*sinphig(ipg)
!        z2 = rg(irg-1)*costhg(itg)

!       Output all points inside the torus
        write(20,'(1p,2e23.15)')  x1, z1

!       Output the inner part of the torus envelope
        if(emdg(irg-1,itg,ipg)<emd0 .and. iflagin==0) write(21,'(1p,2e23.15)')  x1,z1
        if(irg==nrgin+1)  iflagin=1 

!       Output the outer part of the torus envelope
        if(emdg(irg+1,itg,ipg)<emd0 ) write(21,'(1p,2e23.15)')  x1,z1  
      end if
    end do
  end do
  write(20,'(1x)') 
  write(21,'(1x)') 
!
  close(20)
  close(21)
!  close(23)
!  close(24)
!
end subroutine printout_BHT_shape_seq
