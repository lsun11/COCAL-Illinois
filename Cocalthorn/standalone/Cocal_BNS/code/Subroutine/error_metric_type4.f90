subroutine error_metric_type4(pot,pot_bak,error,flag,ctype,cpot)
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg, eps, npgyzp, npgyzm, &
  &                          ntgpolp, ntgpolm, ntgeq, npgxzp, npgxzm
  use coordinate_grav_r, only : rg
  implicit none
  real(long), pointer :: pot(:,:,:), pot_bak(:,:,:)
  real(long), intent(out) :: error
  integer,    intent(out) :: flag
  character(len=2), intent(in) :: ctype
  character(len=4), intent(in) :: cpot
  real(long) :: error_pot = 0.0d0, small = 1.0d-14
  integer    :: irg, itg, ipg, nrgmin, nrgmax, ire, ite, ipe
!
  error = 0.0d0
  ire = 0
  ite = 0
  ipe = 0
  flag = 0
  nrgmin = 0; nrgmax = nrg-1
  if (ctype.eq.'bh') nrgmin = 1
  do irg = nrgmin, nrgmax
    do itg = 1, ntg-1
      do ipg = 1, npg-1
!        if ( dabs(pot(irg,itg,ipg)) < 1.0d-12 ) cycle
        if ( dabs(pot(irg,itg,ipg)) < 1.0d-8 ) cycle
!        if (itg.eq.ntgpolp.or.itg.eq.ntgpolm.or.itg.eq.ntgeq.or. &
!        &   ipg.eq.npgxzp .or.ipg.eq.npgxzm.or.ipg.eq.npgyzp.or.ipg.eq.npgyzm) cycle
        if (itg.eq.ntgeq.or.ipg.eq.npgxzm.or.ipg.eq.npgyzp.or.ipg.eq.npgyzm) cycle

        error_pot = 2.0d0*dabs(pot(irg,itg,ipg) -     pot_bak(irg,itg,ipg)) &
      &                 /(dabs(pot(irg,itg,ipg))+dabs(pot_bak(irg,itg,ipg)) &
      &                 + small)
        if (error_pot > eps) flag = 1
        if(error_pot > error) then
          error = error_pot
          ire = irg
          ite = itg
          ipe = ipg
        end if
      end do
    end do
  end do
  write(6,'(1a19,3i5,a10,1p,e14.6,a2,a4,a1,1p,e20.12,a8,1p,e14.6)') '-- max error grid =',ire,ite,ipe, &
    &            "  r coord=", rg(ire), "  ", cpot, "=", pot(ire,ite,ipe), "  Error=", error
end subroutine error_metric_type4
