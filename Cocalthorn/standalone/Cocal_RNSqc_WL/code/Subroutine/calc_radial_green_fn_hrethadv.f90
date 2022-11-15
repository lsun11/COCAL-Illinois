subroutine calc_radial_green_fn_hrethadv(omega)
  use coordinate_grav_r, only : rg, hrg
  use phys_constant, only : long, pi
  use grid_parameter, only : nrg, npg, ntg, nlg
  use radial_green_fn_hrethadv
  implicit none
  integer :: il, im, ir, irr, irini
  real(8) :: omega, omegam, omegam_rg, omegam_hrg, omegam_rgout, om2lp1
  real(8) :: rtpi2 = 1.253314154753822d0, rgoutsq, small = 1.0d-13
  real(8) :: x, xx, xnu, xrj, xrjp, xry, xryp
  real(8) :: sbj,  sj,  sby,  sy
  real(8) :: sbjp, sjp, sbyp, syp
!
!-------------------------------------
!--- Computation of function bsjy ---
!-------------------------------------
!
!     bsjy_hrha(irr,il,im,ir) computes, for l=il, m nonzero 
!     -(2l+1)*m*omega*j_l(m*omega r_<)*y_l(m*omega r_>)
!    The omega -> 0 limit of this is r_<^l/r_>^(l+1),
!    and bsjy_hrha replaces the function hfsn in grins
!    For m=0, bsjy_hrha is still r_<^l/r_>^(l+1) 
!                        
!
  bsjy_hrha(:,:,:,:) = 0.0d0
!
  irini = 0
  if (rg(0).le.small) then
    irini = 1
    ir = 0
    il = 0
    do irr = 1, nrg
      bsjy_hrha(irr,il,0,ir) = 1.d0/hrg(irr)
    end do
  end if
! 
  do ir  = irini, nrg
    do irr = 1, nrg
      if (hrg(irr).lt.rg(ir)) then
        do il  = 0, nlg
          do im  = 0, il
            if(im.eq.0)then
              bsjy_hrha(irr,il,im,ir) = hrg(irr)**il/rg(ir)**(il+1)
            else
              omegam = dble(im)*omega
              omegam_hrg = omegam*hrg(irr)
              call sphbess(omegam_hrg,il,sj,sy)
              sbj=sj
              omegam_rg = omegam*rg(ir)
              call sphbess(omegam_rg,il,sj,sy)
              sby=sy
              bsjy_hrha(irr,il,im,ir) = -omegam*(2.d0*dble(il)+1.d0)*sbj*sby
            endif
          end do
        end do
      endif
      if (hrg(irr).ge.rg(ir)) then
        do il  = 0, nlg
          do im  = 0, il
            if(im.eq.0)then
               bsjy_hrha(irr,il,im,ir) = rg(ir)**il/hrg(irr)**(il+1)
            else
              omegam = dble(im)*omega
              omegam_rg = omegam*rg(ir)
              call sphbess(omegam_rg,il,sj,sy)
              sbj=sj
              omegam_hrg = omegam*hrg(irr)
              call sphbess(omegam_hrg,il,sj,sy)
              sby=sy
              bsjy_hrha(irr,il,im,ir) = -omegam*(2.d0*dble(il)+1.d0)*sbj*sby
            endif
          end do
        end do
      endif
    end do
  end do
!
! For the surface integral, construct, for m nonzero
!  sbsjy_hrha(l,m,ir) :=  j_l(m omega r) y_l(m omega R) (m omega)(2l+1)
!  sbsjyp_hrha(l,m,ir) =  j_l(m omega r) y_l'(m omega R)(m omega)^2(2l+1)
!
! For m = 0, use 
!  sbsjy_hrha(l,0,ir)  = r^l/R^(l+1)
!  sbsjyp_hrha(l,0,ir) = - (l+1) r^l/R^(l+2)
!
  sbsjy_hrha(:,:,:)  = 0.0d0
  sbsjyp_hrha(:,:,:) = 0.0d0
!
  irini = 0
! ir = 0 case: nonzero contribution only when l=m=0 for rg(0) = 0.0
  if (rg(0).le.small) then
    irini = 1
    sbsjy_hrha(0,0,0) =  1.0d0/rg(nrg)
    sbsjyp_hrha(0,0,0)= -1.0d0/(rg(nrg)*rg(nrg))
  end if
! 
  do il = 0, nlg
    do im = 0, il
      do ir = irini, nrg
        if(im.eq.0)then
          sbsjy_hrha(il,im,ir) = rg(ir)**il/rg(nrg)**(il+1)
          sbsjyp_hrha(il,im,ir) = -(dble(il)+1.d0)*rg(ir)**il/rg(nrg)**(il+2)
        else
          omegam = dble(im)*omega
          omegam_rgout = omegam*rg(nrg)
          om2lp1 = omegam*(2.d0*dble(il)+1.d0)
          call sphbess_and_dx(omegam_rgout,il,sj,sy,sjp,syp)
          sby = sy
          sbyp = omegam*syp
          omegam_rg = omegam*rg(ir)
          call sphbess(omegam_rg,il,sj,sy)
          sbj = sj
          sbsjy_hrha(il,im,ir)  = - om2lp1*sbj*sby
          sbsjyp_hrha(il,im,ir) = - om2lp1*sbj*sbyp
        endif
      enddo !ir
    enddo ! im 
  enddo ! il
!
  rgoutsq = rg(nrg)*rg(nrg)
  sbsjy_hrha(:,:,:)  =  rgoutsq*sbsjy_hrha(:,:,:)
  sbsjyp_hrha(:,:,:) =  rgoutsq*sbsjyp_hrha(:,:,:) 
!
end subroutine calc_radial_green_fn_hrethadv
