subroutine printout_emfield_strength(iseq)
  use phys_constant,  only : long, Gauss, statV
  use grid_parameter, only : nrg, ntg, npg, ntgeq, nrf
  use def_matter_parameter, only : radi
  use def_faraday_tensor, only : fxd_grid, fyd_grid, fzd_grid, fijd_grid
  implicit none
  integer :: irg, itg, ipg, iseq, ii, jj
  integer :: E_max_gp(3,3), B_max_gp(3,3)
  integer :: Ep_max_gp(3), Bp_max_gp(3)
  integer :: Et_max_gp(3), Bt_max_gp(3)
  real(long) :: E_max(3),B_max(3)
  real(long) :: Ep_max, Et_max, Bp_max, Bt_max  ! Poloidal and Toroidal
  real(long) :: E_cens(4,3), B_cens(4,3)  ! 1 = Center,     2 = Equator
  real(long) :: Ep_cens(4), Et_cens(4)    ! 3 = North Pole, 4 = South Pole
  real(long) :: Bp_cens(4), Bt_cens(4)
  real(long) :: tmpE(3), tmpB(3), tmpEt, tmpEp, tmpBt, tmpBp
  character(len=31) char_31
!
! --- Find maximum value
!
  E_max(1:3) = 0.0d0 ; Ep_max = 0.0d0 ; Et_max = 0.0d0
  B_max(1:3) = 0.0d0 ; Bp_max = 0.0d0 ; Bt_max = 0.0d0
  ipg = 0
  do itg = 0, ntg
    do irg = 0, nrg
!
      tmpE(1) = fxd_grid(irg,itg,ipg)/radi
      tmpE(2) = fyd_grid(irg,itg,ipg)/radi
      tmpE(3) = fzd_grid(irg,itg,ipg)/radi
      tmpB(1) =   fijd_grid(irg,itg,ipg,3)/radi
      tmpB(2) = - fijd_grid(irg,itg,ipg,2)/radi
      tmpB(3) =   fijd_grid(irg,itg,ipg,1)/radi
!
      tmpEp = sqrt(tmpE(1)**2+tmpE(3)**2)
      tmpBp = sqrt(tmpB(1)**2+tmpB(3)**2)
      tmpEt = dabs(tmpE(2))
      tmpBt = dabs(tmpB(2))
!
      do ii = 1, 3
        if (dabs(tmpE(ii)).ge.dabs(E_max(ii))) then 
          E_max(ii) = tmpE(ii)
          E_max_gp(ii,1) = irg ; E_max_gp(ii,2) = itg ; E_max_gp(ii,3) = ipg
        end if
        if (dabs(tmpB(ii)).ge.dabs(B_max(ii))) then
          B_max(ii) = tmpB(ii)
          B_max_gp(ii,1) = irg ; B_max_gp(ii,2) = itg ; B_max_gp(ii,3) = ipg
        end if
      end do
      if (tmpEp.ge.Ep_max) then
        Ep_max = tmpEp
        Ep_max_gp(1) = irg ; Ep_max_gp(2) = itg ; Ep_max_gp(3) = ipg
      end if
      if (tmpEt.ge.Et_max) then
        Et_max = tmpEt
        Et_max_gp(1) = irg ; Et_max_gp(2) = itg ; Et_max_gp(3) = ipg
      end if
      if (tmpBp.ge.Bp_max) then
        Bp_max = tmpBp
        Bp_max_gp(1) = irg ; Bp_max_gp(2) = itg ; Bp_max_gp(3) = ipg
      end if
      if (tmpBt.ge.Bt_max) then
        Bt_max = tmpBt
        Bt_max_gp(1) = irg ; Bt_max_gp(2) = itg ; Bt_max_gp(3) = ipg
      end if
    end do
  end do
!
  do ii = 1, 4
    if (ii.eq.1) then ; irg = 0   ; itg = 0     ; ipg = 0 ; end if
    if (ii.eq.2) then ; irg = nrf ; itg = ntgeq ; ipg = 0 ; end if
    if (ii.eq.3) then ; irg = nrf ; itg = 0     ; ipg = 0 ; end if
    if (ii.eq.4) then ; irg = nrf ; itg = ntg   ; ipg = 0 ; end if
    E_cens(ii,1) = fxd_grid(irg,itg,ipg)/radi
    E_cens(ii,2) = fyd_grid(irg,itg,ipg)/radi
    E_cens(ii,3) = fzd_grid(irg,itg,ipg)/radi
    B_cens(ii,1) =   fijd_grid(irg,itg,ipg,3)/radi
    B_cens(ii,2) = - fijd_grid(irg,itg,ipg,2)/radi
    B_cens(ii,3) =   fijd_grid(irg,itg,ipg,1)/radi
    Ep_cens(ii) = sqrt(E_cens(ii,1)**2 + E_cens(ii,3)**2)
    Bp_cens(ii) = sqrt(B_cens(ii,1)**2 + B_cens(ii,3)**2)
    Et_cens(ii) = dabs(E_cens(ii,2))
    Bt_cens(ii) = dabs(B_cens(ii,2))
  end do
!
  if (iseq.eq.1) then 
    open(190,file='rnsEMFseq.dat',status='unknown')
  end if
  write(190,*) '== Sequence number == ', iseq
  write(190,*) '## EM fields in G = c = M = 1/4pi epsi = 1 unit ##'
!
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' Ex   max  =', E_max(1), &
  &             '(',E_max_gp(1,1),',',E_max_gp(1,2),',',E_max_gp(1,3),')'
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' Ey   max  =', E_max(2), &
  &             '(',E_max_gp(2,1),',',E_max_gp(2,2),',',E_max_gp(2,3),')'
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' Ez   max  =', E_max(3), &
  &             '(',E_max_gp(3,1),',',E_max_gp(3,2),',',E_max_gp(3,3),')'
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' Epol max  =', Ep_max, &
  &                   '(',Ep_max_gp(1),',',Ep_max_gp(2),',',Ep_max_gp(3),')'
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' Etor max  =', Et_max, &
  &                   '(',Et_max_gp(1),',',Et_max_gp(2),',',Et_max_gp(3),')'
  do ii = 1, 4
    if (ii.eq.1) then ; irg = 0   ; itg = 0     ; ipg = 0 ; end if
    if (ii.eq.2) then ; irg = nrf ; itg = ntgeq ; ipg = 0 ; end if
    if (ii.eq.3) then ; irg = nrf ; itg = 0     ; ipg = 0 ; end if
    if (ii.eq.4) then ; irg = nrf ; itg = ntg   ; ipg = 0 ; end if
    if (ii.eq.1) char_31 = ' -- Value at the center     -- '
    if (ii.eq.2) char_31 = ' -- Value at the equator    -- '
    if (ii.eq.3) char_31 = ' -- Value at the north pole -- '
    if (ii.eq.4) char_31 = ' -- Value at the south pole -- '
    write(190,'((a31,3x,3(a1,i3),a1))') char_31,'(',irg,',',itg,',',ipg,')'
    write(190,'(a12,1p,3e22.14)') '(Ex,Ey,Ez) =', (E_cens(ii,jj),jj=1,3)
    write(190,'(a12,1p,2e22.14)') '(Epol,Etor)=', Ep_cens(ii), Et_cens(ii)
  end do
!
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' Bx   max  =', B_max(1), &
  &             '(',B_max_gp(1,1),',',B_max_gp(1,2),',',B_max_gp(1,3),')'
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' By   max  =', B_max(2), &
  &             '(',B_max_gp(2,1),',',B_max_gp(2,2),',',B_max_gp(2,3),')'
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' Bz   max  =', B_max(3), &
  &             '(',B_max_gp(3,1),',',B_max_gp(3,2),',',B_max_gp(3,3),')'
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' Bpol max  =', Bp_max, &
  &                   '(',Bp_max_gp(1),',',Bp_max_gp(2),',',Bp_max_gp(3),')'
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' Btor max  =', Bt_max, &
  &                   '(',Bt_max_gp(1),',',Bt_max_gp(2),',',Bt_max_gp(3),')'
  do ii = 1, 4
    if (ii.eq.1) then ; irg = 0   ; itg = 0     ; ipg = 0 ; end if
    if (ii.eq.2) then ; irg = nrf ; itg = ntgeq ; ipg = 0 ; end if
    if (ii.eq.3) then ; irg = nrf ; itg = 0     ; ipg = 0 ; end if
    if (ii.eq.4) then ; irg = nrf ; itg = ntg   ; ipg = 0 ; end if
    if (ii.eq.1) char_31 = ' -- Value at the center     -- '
    if (ii.eq.2) char_31 = ' -- Value at the equator    -- '
    if (ii.eq.3) char_31 = ' -- Value at the north pole -- '
    if (ii.eq.4) char_31 = ' -- Value at the south pole -- '
    write(190,'((a31,3x,3(a1,i3),a1))') char_31,'(',irg,',',itg,',',ipg,')'
    write(190,'(a12,1p,3e22.14)') '(Bx,By,Bz) =', (B_cens(ii,jj),jj=1,3)
    write(190,'(a12,1p,2e22.14)') '(Bpol,Btor)=', Bp_cens(ii), Bt_cens(ii)
  end do
!
  write(190,*) '## EM fields in Gauss ##' 
!
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' Ex   max  =',E_max(1)*statV,&
  &             '(',E_max_gp(1,1),',',E_max_gp(1,2),',',E_max_gp(1,3),')'
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' Ey   max  =',E_max(2)*statV,&
  &             '(',E_max_gp(2,1),',',E_max_gp(2,2),',',E_max_gp(2,3),')'
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' Ez   max  =',E_max(3)*statV,&
  &             '(',E_max_gp(3,1),',',E_max_gp(3,2),',',E_max_gp(3,3),')'
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' Epol max  =',Ep_max*statV,&
  &                   '(',Ep_max_gp(1),',',Ep_max_gp(2),',',Ep_max_gp(3),')'
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' Etor max  =',Et_max*statV,&
  &                   '(',Et_max_gp(1),',',Et_max_gp(2),',',Et_max_gp(3),')'
  do ii = 1, 4
    if (ii.eq.1) then ; irg = 0   ; itg = 0     ; ipg = 0 ; end if
    if (ii.eq.2) then ; irg = nrf ; itg = ntgeq ; ipg = 0 ; end if
    if (ii.eq.3) then ; irg = nrf ; itg = 0     ; ipg = 0 ; end if
    if (ii.eq.4) then ; irg = nrf ; itg = ntg   ; ipg = 0 ; end if
    if (ii.eq.1) char_31 = ' -- Value at the center     -- '
    if (ii.eq.2) char_31 = ' -- Value at the equator    -- '
    if (ii.eq.3) char_31 = ' -- Value at the north pole -- '
    if (ii.eq.4) char_31 = ' -- Value at the south pole -- '
    write(190,'((a31,3x,3(a1,i3),a1))') char_31,'(',irg,',',itg,',',ipg,')'
    write(190,'(a12,1p,3e22.14)') '(Ex,Ey,Ez) =', (E_cens(ii,jj)*statV,jj=1,3)
    write(190,'(a12,1p,2e22.14)') '(Epol,Etor)=', Ep_cens(ii)*statV, &
    &                                             Et_cens(ii)*statV
  end do
!
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' Bx   max  =',B_max(1)*Gauss,&
  &             '(',B_max_gp(1,1),',',B_max_gp(1,2),',',B_max_gp(1,3),')'
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' By   max  =',B_max(2)*Gauss,&
  &             '(',B_max_gp(2,1),',',B_max_gp(2,2),',',B_max_gp(2,3),')'
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' Bz   max  =',B_max(3)*Gauss,&
  &             '(',B_max_gp(3,1),',',B_max_gp(3,2),',',B_max_gp(3,3),')'
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' Bpol max  =',Bp_max*Gauss,&
  &                   '(',Bp_max_gp(1),',',Bp_max_gp(2),',',Bp_max_gp(3),')'
  write(190,'(a12,1p,1e22.14,3x,3(a1,i3),a1)') ' Btor max  =',Bt_max*Gauss,&
  &                   '(',Bt_max_gp(1),',',Bt_max_gp(2),',',Bt_max_gp(3),')'
  do ii = 1, 4
    if (ii.eq.1) then ; irg = 0   ; itg = 0     ; ipg = 0 ; end if
    if (ii.eq.2) then ; irg = nrf ; itg = ntgeq ; ipg = 0 ; end if
    if (ii.eq.3) then ; irg = nrf ; itg = 0     ; ipg = 0 ; end if
    if (ii.eq.4) then ; irg = nrf ; itg = ntg   ; ipg = 0 ; end if
    if (ii.eq.1) char_31 = ' -- Value at the center     -- '
    if (ii.eq.2) char_31 = ' -- Value at the equator    -- '
    if (ii.eq.3) char_31 = ' -- Value at the north pole -- '
    if (ii.eq.4) char_31 = ' -- Value at the south pole -- '
    write(190,'((a31,3x,3(a1,i3),a1))') char_31,'(',irg,',',itg,',',ipg,')'
    write(190,'(a12,1p,3e22.14)') '(Bx,By,Bz) =', (B_cens(ii,jj)*Gauss,jj=1,3)
    write(190,'(a12,1p,2e22.14)') '(Bpol,Btor)=', Bp_cens(ii)*Gauss, &
    &                                             Bt_cens(ii)*Gauss
  end do
!
end subroutine printout_emfield_strength
