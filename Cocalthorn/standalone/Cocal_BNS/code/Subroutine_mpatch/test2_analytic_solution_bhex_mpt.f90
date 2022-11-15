subroutine test2_analytic_solution_bhex_mpt(impt)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrg, ntg, npg, nrf, rgin
  use def_metric, only  :   bvxd, bvyd
  use coordinate_grav_r, only : rg
  use grid_points_binary_excision, only : rb
  use grid_parameter_mpt, only : grid_param_real_
  use grid_points_asymptotic_patch_mpt, only : ra_
  use make_array_3d
  use copy_array_4dto3d_mpt
  implicit none
  integer     ::   irg,itg,ipg, impt, impt_ce, impt_ex
  real(long)  ::   zfac, small = 1.0d-15, alps_tmp
  real(long)  ::   mass_ce, mass_ex, rgin_ce, rgin_ex, rad_1, rad_2
  real(long), pointer :: ra_1(:,:,:), ra_2(:,:,:)
!
  call alloc_array3d(ra_1,-2,nrg+2,0,ntg,0,npg)
  call alloc_array3d(ra_2,-2,nrg+2,0,ntg,0,npg)
!
  if (impt.eq.1) then ; impt_ce = impt ; impt_ex = 2 ; end if
  if (impt.eq.2) then ; impt_ce = impt ; impt_ex = 1 ; end if
  if (impt.eq.3) then ; impt_ce = 1    ; impt_ex = 2 ; end if
!
  rgin_ce = grid_param_real_(1,impt_ce)
  rgin_ex = grid_param_real_(1,impt_ex)
  mass_ce = 2.0d0*rgin_ce*0.8d0
  mass_ex = 2.0d0*rgin_ex*0.8d0
!
  if (impt.eq.3) then 
    call copy_array4dto3d_mpt(impt_ce,ra_,ra_1,-2,nrg+2,0,ntg,0,npg)
    call copy_array4dto3d_mpt(impt_ex,ra_,ra_2,-2,nrg+2,0,ntg,0,npg)
  end if
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        if (impt.ne.3) then 
          rad_1 = rg(irg)
          rad_2 = rb(irg,itg,ipg)
        end if 
        if (impt.eq.3) then 
          rad_1 = ra_1(irg,itg,ipg) 
          rad_2 = ra_2(irg,itg,ipg)
        end if
        bvxd(irg,itg,ipg) = 1.0d0 + mass_ce/(2.0d0*rad_1) &
        &                         + mass_ex/(2.0d0*rad_2)
        alps_tmp          = 1.0d0 - mass_ce/(2.0d0*rad_1) &
        &                         - mass_ex/(2.0d0*rad_2)
        bvyd(irg,itg,ipg) = alps_tmp/bvxd(irg,itg,ipg)
      end do
    end do
  end do
!
  deallocate(ra_1)
  deallocate(ra_2)
end subroutine test2_analytic_solution_bhex_mpt
