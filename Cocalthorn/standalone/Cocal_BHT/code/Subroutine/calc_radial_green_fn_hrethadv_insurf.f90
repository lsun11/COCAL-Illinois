subroutine calc_radial_green_fn_hrethadv_insurf(omega)
  use coordinate_grav_r, only : rg, hrg
  use phys_constant, only : long, pi
  use grid_parameter, only : nrg, npg, ntg, nlg
  use radial_green_fn_hrethadv
  implicit none
  integer :: il, im, ir, irr
  real(8) :: omega, omegam, omegam_rg, omegam_hrg, omegam_rgin, om2lp1
  real(8) :: rtpi2 = 1.253314154753822d0, rgoutsq, small = 1.0d-13
  real(8) :: x, xx, xnu, xrj, xrjp, xry, xryp
  real(8) :: sbj,  sj,  sby,  sy
  real(8) :: sbjp, sjp, sbyp, syp
!
!-------------------------------------------------------------------
!--- Computation of function sbsjy at the inner surface boundary ---
!--- This over writes the sbsjy_hrha sbsjyp_hrha
!-------------------------------------------------------------------
!
!     bsjy_hrha(irr,il,im,ir) computes, for l=il, m nonzero 
!     -(2l+1)*m*omega*j_l(m*omega r_<)*y_l(m*omega r_>)
!    The omega -> 0 limit of this is r_<^l/r_>^(l+1),
!    and bsjy_hrha replaces the function hfsn in grins
!    For m=0, bsjy_hrha is still r_<^l/r_>^(l+1) 
!
! For the surface integral, construct, for m nonzero
!  sbsjy_hrha(l,m,ir) :=  j_l(m omega R) y_l(m omega r) (m omega)(2l+1)
!  sbsjyp_hrha(l,m,ir) =  j_l'(m omega R) y_l(m omega R)(m omega)^2(2l+1)
!
! For m = 0, use 
!  sbsjy_hrha(l,0,ir)  = R^l/r^(l+1)
!  sbsjyp_hrha(l,0,ir) = l R^(l-1)/r^(l+1)
!
! Note: Minus sign for inward normal is included in radial green's functions.
! Note: This routine is used only for BH, rg(0) > 0.
!
  sbsjy_hrha(:,:,:)  = 0.0d0
  sbsjyp_hrha(:,:,:) = 0.0d0
  do il = 0, nlg
    do im = 0, il
      do ir = 0, nrg
        if(im.eq.0)then
          sbsjy_hrha(il,im,ir) = - rg(0)**il/rg(irg)**(il+1)
          if (il.ne.0) sbsjyp_hrha(il,im,ir) & 
          &       = - dble(il)*rg(0)**(il-1)/rg(irg)**(il+1)
        else
          omegam = dble(im)*omega
          omegam_rgin = omegam*rg(0)
          om2lp1 = omegam*(2.d0*dble(il)+1.d0)
          call sphbess_and_dx(omegam_rgin,il,sj,sy,sjp,syp)
          sbj = sj
          sbjp = omegam*sjp
          omegam_rg = omegam*rg(ir)
          call sphbess(omegam_rg,il,sj,sy)
          sby = sy
          sbsjy_hrha(il,im,ir)  = om2lp1*sbj*sby
          sbsjyp_hrha(il,im,ir) = om2lp1*sbjp*sby
        endif
      enddo !ir
    enddo ! im 
  enddo ! il
  rgoutsq = rg(0)*rg(0)
  sbsjy_hrha(:,:,:)  =  rgoutsq*sbsjy_hrha(:,:,:)
  sbsjyp_hrha(:,:,:) =  rgoutsq*sbsjyp_hrha(:,:,:) 
!
end subroutine calc_radial_green_fn_hrethadv_insurf
