subroutine error_metric_type0(pot,pot_bak,error,flag,ctype)
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg, eps_coc, &
  &                          ntgpolp, ntgpolm, ntgeq, npgxzp, npgxzm
  implicit none
  real(long), pointer :: pot(:,:,:), pot_bak(:,:,:)
  real(long), intent(out) :: error
  integer,    intent(out) :: flag
  real(long) :: error_pot = 0.0d0, small = 1.0d-14
  integer    :: irg, itg, ipg, nrgmin, nrgmax, ire, ite, ipe
  character(len=2), intent(in) :: ctype
!
  error = 0.0d0
  flag = 0
  nrgmin = 0; nrgmax = nrg-1
  if (ctype.eq.'bh') nrgmin = 1
  do irg = nrgmin, nrgmax
    do itg = 0, ntg
      do ipg = 0, npg
        if (itg.eq.ntgpolp.or.itg.eq.ntgpolm.or.itg.eq.ntgeq.or. &
        &   ipg.eq.npgxzp .or.ipg.eq.npgxzm) cycle
        error_pot = 2.0d0*dabs(pot(irg,itg,ipg) -     pot_bak(irg,itg,ipg)) &
      &                 /(dabs(pot(irg,itg,ipg))+dabs(pot_bak(irg,itg,ipg)) &
      &                 + small)
        if (error_pot > eps_coc) flag = 1
        if(error_pot > error) then
          error = error_pot
          ire = irg
          ite = itg
          ipe = ipg
        end if
      end do
    end do
  end do
write(6,'(1a19,3i5,1p,2e15.6)') '-- max error grid =',ire,ite,ipe, &
&                    pot(ire,ite,ipe), error
end subroutine error_metric_type0
