subroutine calc_surface(rsnew,emd)
  use phys_constant, only  : long
  use grid_parameter, only : nrf, ntf, npf
  use coordinate_grav_r, only : drg, rg
  use def_matter, only : rs
!  use def_matter_parameter, only : emdc
  implicit none
  real(long), pointer :: emd(:,:,:)
  real(long), pointer :: rsnew(:,:)
  real(long) :: delta
  integer    :: itf, ipf
  real(long) :: sgg(1:4), gg(1:4,1:4)
  real(long) :: r1,r2,r3,q1,q2,q3, aa,bb,cc,det,rs1,rs2,rs2a,rs2b
!
  do ipf = 0, npf
    do itf = 0, ntf
      delta = drg(nrf)*(0.0d0            - emd(nrf  ,itf,ipf)) &
    &                 /(emd(nrf,itf,ipf) - emd(nrf-1,itf,ipf))
!robust    delta = (drg(nrf)+drg(nrf-1)) &
!robust  &        *(0.0d0            - emd(nrf  ,itf,ipf)) &
!robust  &        /(emd(nrf,itf,ipf) - emd(nrf-2,itf,ipf))
      rsnew(itf,ipf) = rs(itf,ipf)*(1.0d0 + delta)
!
! =============== 3 point interpolation ============================
      q1 = emd(nrf  ,itf,ipf);    r1 = rs(itf,ipf)*rg(nrf  )
      q2 = emd(nrf-1,itf,ipf);    r2 = rs(itf,ipf)*rg(nrf-1)
      q3 = emd(nrf-2,itf,ipf);    r3 = rs(itf,ipf)*rg(nrf-2)

      sgg(1) = q1
      sgg(2) = q2
      sgg(3) = q3

      gg(1,1) = r1*r1
      gg(1,2) = r1
      gg(1,3) = 1.0d0

      gg(2,1) = r2*r2
      gg(2,2) = r2
      gg(2,3) = 1.0d0

      gg(3,1) = r3*r3
      gg(3,2) = r3
      gg(3,3) = 1.0d0

      call minv(gg,sgg,3,4)
      aa  = sgg(1)
      bb  = sgg(2)
      cc  = sgg(3)
      det = bb*bb-4.0d0*aa*cc     

      rs1 = rs(itf,ipf)*(1.0d0 + delta)

      rs2 = 0.0d0
      !rs2a = (-bb+sqrt(det))/2.0d0/aa
!uryu      if( det >= 0 )  then
      if( det >= 0 .and. aa > 0 )  then
        rs2 = (-bb-sqrt(det))/2.0d0/aa     ! Correct solution 
!uryu        if ((itf==(ntf/2).or.itf==(ntf/4)).and.ipf<=2) &
!uryu          write(6,'(2i4,1p,10e15.5)') &
!uryu          itf,ipf,rs1,rs2,r1,r2,r3,aa,bb,cc,det
      else
!uryu        rs2 = -bb/2.0d0/aa       ! Approximate solution. Connect to minimum.
        rs2 = rsnew(itf,ipf)
!uryu        if ((itf==(ntf/2).or.itf==(ntf/4)).and.ipf<=2)  &
!uryu               write(6,'(2i4,1p,e15.5,1p,10e15.5)')  &
!uryu               itf,ipf,  rs1, rs2, r1,r2,r3, aa,bb,cc,det
      end if
      rsnew(itf,ipf) = rs2
! ==================================================================
!
    end do
  end do
  rsnew(0,1:npf) = rsnew(0,0)
end subroutine calc_surface
