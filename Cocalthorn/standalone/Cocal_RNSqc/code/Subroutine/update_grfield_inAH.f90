subroutine update_grfield_inAH(pot,grfield,convf)
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg, rgin
  use coordinate_grav_r, only : rg
  use def_kerr_schild
  use def_horizon, only : ahz
  use make_array_3d
  implicit none
  real(long), pointer    :: pot(:,:,:)
  real(long), pointer    :: grfield(:,:,:)
  real(long), intent(in) :: convf
  real(long) :: cf
  real(long) :: conv_rgin,conv_ah,slope,rcenter,f_rgin,f_ah,aa,bb,ff
  integer    :: irg,itg,ipg
!
! Convergence factors at AH and at rgin. slope parameter determines
! the steepness of tanh() between these two values.
  conv_rgin = 0.05d0
  conv_ah   = convf
  slope     = 1.0d+03 

! rcenter is the center of the tanh function. 0.5 corresponds
! exactly at the middle, 0.7 makes the center of tanh to be closer to 
! to the AH thus most of the domain will have cf~conv_rgin
  rcenter = rgin + 0.8*(reh_ks - rgin)
  f_rgin  = 1.0d0+tanh(slope*(rgin-rcenter))
  f_ah    = 1.0d0+tanh(slope*(reh_ks-rcenter))

  bb = (conv_ah - conv_rgin)/(f_ah - f_rgin)
  aa = conv_ah - bb*f_ah 

  do irg=0,nrg
    do itg=0,ntg
      do ipg=0,npg
        if (rg(irg)> ahz(itg,ipg))  then
          cf = conv_ah
        else
          ff = 1.0d0+tanh(slope*(rg(irg)-rcenter)) 
          cf = aa+bb*ff
        end if
!        if(itg==0.and.ipg==0.and.irg<20)  write(6,*) irg,rg(irg),cf

        grfield(irg,itg,ipg) = cf*pot(irg,itg,ipg) + (1.0d0-cf)*grfield(irg,itg,ipg)   
      end do
    end do
  end do

end subroutine update_grfield_inAH
