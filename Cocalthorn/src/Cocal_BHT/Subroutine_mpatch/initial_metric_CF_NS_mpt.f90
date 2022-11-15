subroutine initial_metric_CF_NS_mpt(impt)
  use phys_constant, only  : long, pi
  use grid_parameter
  use coordinate_grav_r
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi,   only : sinphig, cosphig
  use def_binary_parameter,    only : sepa, dis
  use def_metric
  implicit none
  real(long) :: st, ct, sp, cp, xa,ya,za, rcm2, xycm2, rcm, xycm, tcm, pcm, xcm,ycm,zcm
  real(long) :: rr, work_shift, pari
  integer    :: irg, itg, ipg, impt, npg_l, npg_r
  character(30) :: char1, char2, char3, char4, char5
!
  write( 6,*)  sepa, dis

  if(impt==1)  then
    work_shift = -dis;     pari=  1.0d0;
  else if(impt==2) then
    work_shift = -dis;     pari=  1.0d0;
  else if(impt==3) then
    work_shift = 0.0d0;    pari=  1.0d0;
  end if


  do irg = 0, nrg
    do ipg = 0, npg
      do itg = 0, ntg
        st = sinthg(itg)
        ct = costhg(itg)
        sp = sinphig(ipg)
        cp = cosphig(ipg)
        rr = rg(irg)
! Coordinates wrt the center of black hole
        xa = rr*st*cp
        ya = rr*st*sp
        za = rr*ct
! Coordinates wrt the CM
        rcm2 = (xa-0.5d0*sepa)**2 + ya**2 + za**2
        xycm2= (xa-0.5d0*sepa)**2 + ya**2
        rcm = dsqrt(rcm2)
        xycm= dsqrt(xycm2)
!       tcm = atan2(sqrt(xycm2),za)
!       pcm = dmod(2.0d0*pi+datan2(ya,xa-0.5d0*sepa),2.0d0*pi)
!        xcm = xa - 0.5d0*sepa
!        ycm = ya
!        zcm = za
        psi(irg,itg,ipg)  = 1.0d0
        alph(irg,itg,ipg) = 1.0d0
        alps(irg,itg,ipg) = 1.0d0
!        bvxd(irg,itg,ipg) = 0.0d0
!        bvyd(irg,itg,ipg) = 0.0d0
!        bvzd(irg,itg,ipg) = 0.0d0
        bvxd(irg,itg,ipg) = 0.0d0
        bvyd(irg,itg,ipg) = 0.01d0*pari*( -(xycm-work_shift)*((xycm-work_shift)**4 + 0.5d0**4)**(-0.5) )
        bvzd(irg,itg,ipg) = 0.0d0
      end do
    end do
  end do
!
        if (impt.eq.1) then
          npg_l = npgxzm;  npg_r = npgxzp;  work_shift = - dis;  pari = 1.0;
        else if(impt.eq.2) then
          npg_l = npgxzp;  npg_r = npgxzm;  work_shift = dis;   pari = -1.0;
        else if(impt.eq.3) then
          npg_l = npgxzm;  npg_r = npgxzp;  work_shift = 0.0d0;  pari = 1.0;
        end if
        write(char4, '(i5)') impt
        char5 = adjustl(char4)
        char3 = 'bv_mpt' // trim(char5) // '.txt'
!         open(15,file='../../'//char3,status='unknown')
        open(15,file='./'//char3,status='unknown')
          do irg = nrg, 0, -1
            write(15,'(1p,9e20.12)') -rg(irg)+work_shift                  &
            &                                , bvxd(irg,ntgeq,npg_l)      &
            &                                , bvyd(irg,ntgeq,npg_l) 
          end do
          write(15,*) ' '
          do irg = 0, nrg
            write(15,'(1p,9e20.12)')  rg(irg)+work_shift                  &
            &                                , bvxd(irg,ntgeq,npg_r)      &
            &                                , bvyd(irg,ntgeq,npg_r) 
          end do
        close(15)

end subroutine initial_metric_CF_NS_mpt
