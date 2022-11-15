subroutine error_metric_type2_mpt(pot,pot_bak,error,flag,ctype,impt)
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg, eps_coc, &
  &                          ntgpolp, ntgpolm, ntgeq, npgxzp, npgxzm
  use coordinate_grav_r, only : rg
  implicit none
  real(long), pointer :: pot(:,:,:), pot_bak(:,:,:)
  real(long), intent(out) :: error
  integer,    intent(out) :: flag
  integer,    intent(in)  :: impt
  real(long) :: error_pot = 0.0d0, small = 1.0d-14
  integer    :: irg, itg, ipg, nrgmin, nrgmax, ire, ite, ipe
  character(len=2), intent(in) :: ctype
!
  error = 0.0d0
  ire = 0
  ite = 0
  ipe = 0
  flag = 0
  nrgmin = 0; nrgmax = nrg-1
  if (ctype.eq.'bh') nrgmin = 1

  if (impt>2) nrgmin=2
  if (impt==1 .or. impt==2)  nrgmax=nrg-2

  do irg = nrgmin, nrgmax
    do itg = 0, ntg
      do ipg = 0, npg
!orig        if ( dabs(pot(irg,itg,ipg)) < 1.0e-12 ) cycle
        if ( dabs(pot(irg,itg,ipg)) < 1.0d-8 ) cycle
        if (itg.eq.ntgpolp.or.itg.eq.ntgpolm.or.itg.eq.ntgeq.or. &
        &   ipg.eq.npgxzp .or.ipg.eq.npgxzm .or. ipg.eq.npg) cycle
        error_pot = 2.0d0*dabs(pot(irg,itg,ipg) -     pot_bak(irg,itg,ipg)) &
      &                 /(dabs(pot(irg,itg,ipg))+dabs(pot_bak(irg,itg,ipg)) &
      &                 + small)
        if (error_pot > eps_coc)     flag = 1   ! eps_coc = 1.0e-06
!        error = dmax1(error,error_pot)
        if(error_pot > error) then
          error = error_pot
          ire = irg
          ite = itg
          ipe = ipg
        end if        
      end do
    end do
  end do
!  write(6,*) ''
  write(6,'(1a19,3i5,a10,1p,e14.6,a12,1p,e20.12,a8,1p,e14.6)') '-- max error grid =',ire,ite,ipe, &
    &            "  r coord=", rg(ire), "  Potential=", pot(ire,ite,ipe), "  Error=", error
!
end subroutine error_metric_type2_mpt
