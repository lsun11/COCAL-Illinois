subroutine printout_quad_pole(iseq)
  use phys_constant, only : c, g, msol, pi
  use def_quantities, only : gravmass_sph
  use def_matter_parameter, only : ome, radi
  use def_quantities, only : Iij, Itf, dt1Itf, dt2Itf, dt3Itf, &
  &                          LGW, dJdt, hplus, hcross
  implicit none
  real(8) :: ff, rr, hplus30k, hcross30k, hplus30M, hcross30M
  integer :: i, j
  integer, intent(in) :: iseq
!
  if (iseq.eq.1) then
    open(300,file='rnsquadpole.dat',status='unknown')
  else
    open(300,file='rnsquadpole.dat',status='old',position='append')
  end if
!
  do i=1,3
    do j=1,3
      write (300, '(a36, i1, a1, i1, a4, 1es23.15)') &
      &  'Mass quadrupole moment I(', i, ',', j, ') = ', Iij(i, j)
    end do
  end do
  write (300, *)
  do i=1,3
    do j=1,3
      write (300, '(a36, i1, a1, i1, a4, 1es23.15)') &
      &  'Reduced quadrupole moment Itf(',i,',',j, ') = ', Itf(i,j)
    end do
  end do
  write (300, *)
  do i=1,3
    do j=1,3
      write (300, '(a36, i1, a1, i1, a4, 1es23.15)') &
      &  'First time derivative of QP dt1Itf(',i,',',j, ') = ', dt1Itf(i,j)
    end do
  end do
  write (300, *)
  do i=1,3
    do j=1,3
      write (300, '(a36, i1, a1, i1, a4, 1es23.15)') &
      &  'Second time derivative of QP dt2Itf(',i,',',j, ') = ', dt2Itf(i,j)
    end do
  end do
  write (300, *)
  do i=1,3
    do j=1,3
      write (300, '(a36, i1, a1, i1, a4, 1es23.15)') &
      &  'Third time derivative of QP dt3Itf(',i,',',j, ') = ', dt3Itf(i,j)
    end do
  end do
!
  write (300, *)
  write(300, '(a10, 1es23.15)') 'LGW = ', LGW
  do i = 1, 3
    write(300, '(a5, i1, a4, 1es23.15)') 'dJdt(', i, ') = ', dJdt(i)
  end do
!
  write(300, *)
!
!3.24d-14pc=1km, 1=147670cm
!1=3.24d-14*137670/100000pc
  rr = 30.0d+8/(3.24d-14*147670.0d0)
  hplus30k  = hplus/rr
  hcross30k = hcross/rr
  write(300, '(a16, 1es23.15)') 'h+ (r=30kpc) = ', hplus30k
  write(300, '(a16, 1es23.15)') 'hx (r=30kpc) = ', hcross30k
!
  rr = 30.0d+11/(3.24d-14*147670.0d0)
  hplus30M  = hplus/rr
  hcross30M = hcross/rr
  write(300, '(a16, 1es23.15)') 'h+ (r=30Mpc) = ', hplus30M
  write(300, '(a16, 1es23.15)') 'hx (r=30Mpc) = ', hcross30M
!
  ff = (ome/radi) *c**3 / (pi*g*msol)
  write(300, '(a16, 1es23.15)') 'frequency[Hz] = ', ff
!
  close(300)
!
end subroutine printout_quad_pole
