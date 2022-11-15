subroutine calc_quad_pole_qeos
  use phys_constant, only  :   long, pi
  use grid_parameter, only :   nrf, ntf, npf
  use def_matter_parameter, only : radi, ome
  use make_array_3d
  use make_array_5d
  use def_quantities, only : restmass, Iij, Itf, dt1Itf, dt2Itf, dt3Itf, &
  &                          LGW, dJdt, hplus, hcross
  use interface_source_quad_pole_qeos
  use interface_vol_int_fluid
  implicit none
  real(long)          :: fac2pi, fac4pi
  real(long)          :: volf
  real(long), pointer :: souf(:,:,:), souf5d(:,:,:,:,:)
  real(long)          :: theta, phi, rr
  integer :: ir, it, ip, i, j, k, l, iquad
!
  call alloc_array3d(souf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array5d(souf5d, 0, nrf, 0, ntf, 0, npf, 1, 3, 1, 3)
!
  do iquad = 0, 4
    call source_quad_pole_qeos(souf5d, iquad)
    do i = 1, 3
      do j = 1, 3
        souf(0:nrf,0:ntf,0:npf) = souf5d(0:nrf,0:ntf,0:npf,i,j)
        call vol_int_fluid(souf,volf)
        if (iquad == 0) Iij(i,j) = radi**5*volf
        if (iquad == 1) Itf(i,j) = radi**5*volf
        if (iquad == 2) dt1Itf(i,j) = ome*radi**4*volf
        if (iquad == 3) dt2Itf(i,j) = -ome**2*radi**3*volf
        if (iquad == 4) dt3Itf(i,j) = -ome**3*radi**2*volf
      end do
    end do
  end do
!
  LGW = 0.0d0
  dJdt(1:3) = 0.0d0
  do i = 1, 3
    do j = 1, 3
      LGW = LGW + dt3Itf(i,j)*dt3Itf(i,j)/5.0d0
    end do
    if (i.eq.1) then; j = 2; k = 3; end if
    if (i.eq.2) then; j = 3; k = 1; end if
    if (i.eq.3) then; j = 1; k = 2; end if
    do l = 1, 3
      dJdt(i) = dJdt(i) &
    & + (dt2Itf(j,l)*dt3Itf(l,k) - dt2Itf(k,l)*dt3Itf(l,j))*2.0d0/5.0d0
    end do
  end do
!
  theta = 0.0d0
  phi = 0.0d0
  rr = 1.0d0
  hplus = ((dt2Itf(1,1)-dt2Itf(2,2))*(cos(theta)**2+1.0d0)*cos(2*phi)/4.0d0 &
     &  -  (dt2Itf(1,1)+dt2Itf(2,2)-2.0d0*dt2Itf(3,3))*sin(theta)**2/4.0d0 &
     &  +   dt2Itf(1,2)*((cos(theta)**2+1.0d0)/2.0d0)*sin(2*phi) &
     &  -   dt2Itf(1,3)*sin(theta)*cos(theta)*cos(phi) &
     &  -   dt2Itf(2,3)*sin(theta)*cos(theta)*sin(phi) &
     &    )*2.0d0/rr
  hcross =( - (dt2Itf(1,1)-dt2Itf(2,2))*cos(theta)*sin(2*phi)/2.0d0 &
     &   +     dt2Itf(1,2)*cos(theta)*cos(2*phi) &
     &   +     dt2Itf(1,3)*sin(theta)*sin(phi) &
     &   -     dt2Itf(2,3)*sin(theta)*cos(phi) &
     &    )*2.0d0/rr
!
!  write (6,'(a20,1p,e14.6)') ' Rest  mass =       ',  restmass
!
  deallocate(souf)
  deallocate(souf5d)
end subroutine calc_quad_pole_qeos
