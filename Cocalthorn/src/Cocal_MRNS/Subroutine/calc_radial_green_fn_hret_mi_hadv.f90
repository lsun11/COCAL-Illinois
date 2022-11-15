subroutine calc_radial_green_fn_hret_mi_hadv(omega)
  use coordinate_grav_r, only : rg, hrg
  use phys_constant, only : long, pi
  use grid_parameter, only : nrg, npg, ntg, nlg
  use radial_green_fn_hret_mi_hadv
  implicit none
  integer :: il, im, ir, irr
  real(8) :: omega, omegam, omegam_rg, omegam_hrg, om2np1
  real(8) :: rtpi2 = 1.253314154753822d0
  real(8) :: x, xx, xnu, xrj, xrjp, xry, xryp
  real(8) :: sbj, sj
  real(8) :: sby, sy, syp
!
!-------------------------------------
!--- Computation of function bsjy ---
!-------------------------------------
!
!     bsjy(irr,il,im,ir) computes, for l=il, m nonzero 
!     -(2l+1)*m*omega*j_l(m*omega r_<)*j_l(m*omega r_>)
!    The omega -> 0 limit of this is r_<^l/r_>^(l+1),
!    and bsjy replaces the function hfsn in grins
!    For m=0, bsjy is still r_<^l/r_>^(l+1) 
!                        
!
  bsjy_hrmiha(:,:,:,:) = 0.0d0
!
  do ir  = 1, nrg
    do irr = 1, nrg
      if (hrg(irr).lt.rg(ir)) then
        do il  = 0, nlg
          do im  = 0, il
            if(im.eq.0)then
!               bsjy_hrmiha(irr,il,im,ir) = hrg(irr)**il/rg(ir)**(il+1)
               bsjy_hrmiha(irr,il,im,ir) = 0.0d0  ! 1/2(Gout - Gin)
            else
              omegam = dble(im)*omega
              omegam_hrg = omegam*hrg(irr)
              call sphbess(omegam_hrg,il,sj,sy)
              sbj=sj
              omegam_rg = omegam*rg(ir)
              call sphbess(omegam_rg,il,sj,sy)
              sby=sj  ! 1/2(Gout - Gin)
              bsjy_hrmiha(irr,il,im,ir) = -omegam*(2.d0*dble(il)+1.d0)*sbj*sby
            endif
          end do
        end do
      endif
      if (hrg(irr).ge.rg(ir)) then
        do il  = 0, nlg
          do im  = 0, il
            if(im.eq.0)then
!              bsjy_hrmiha(irr,il,im,ir) = rg(ir)**il/hrg(irr)**(il+1)
              bsjy_hrmiha(irr,il,im,ir) = 0.0d0  ! 1/2(Gout - Gin)
            else
              omegam = dble(im)*omega
              omegam_rg = omegam*rg(ir)
              call sphbess(omegam_rg,il,sj,sy)
              sbj=sj
              omegam_hrg = omegam*hrg(irr)
              call sphbess(omegam_hrg,il,sj,sy)
              sby=sj  ! 1/2(Gout - Gin)
              bsjy_hrmiha(irr,il,im,ir) = -omegam*(2.d0*dble(il)+1.d0)*sbj*sby
            endif
          end do
        end do
      endif
    end do
  end do
!
!
  ir  = 0
  il  = 0
  do irr = 1, nrg
    bsjy_hrmiha(irr,il,0,ir) = 1.d0/hrg(irr)
  end do
! 
! For the surface integral, construct, for m nonzero
!  sbsjy(l,m,ir) :=  j_l(m omega r) y_l(m omega R) (m omega)(2l+1)
!  sbsjyp(l,m,ir) =  j_l(m omega r) y_l'(m omega R)(m omega)^2(2l+1)
!
! For m = 0, use 
!  sbsjy(l,0,ir)  = r^l/R^(l+1)
!  sbsjyp(l,0,ir) = - (l+1) r^l/R^(l+2)
!
  do il = 0,nlg
    do im = 1,il
      omegam=dble(im)*omega
      x = omegam*rg(nrg)
      xnu = dble(il)+0.5d0
      om2np1=omegam*(2.d0*dble(il)+1.d0)
      call bessjy(x,xnu,xrj,xry,xrjp,xryp)
      sy = xrj*rtpi2/sqrt(x)  ! 1/2(Gout - Gin)
      syp = rtpi2/sqrt(x)*(xrjp - xrj*0.5d0/x)   ! 1/2(Gout - Gin)
      do ir = 1,nrg
        if(im.eq.0)then
          sbsjy_hrmiha(il,im,ir) = rg(ir)**il/rg(nrg)**(il+1)
          sbsjyp_hrmiha(il,im,ir) = -(dble(il)+1.d0)*rg(ir)**il/rg(nrg)**(il+2)
        else
          xx=omegam*rg(ir)
          call bessjy(xx,xnu,xrj,xry,xrjp,xryp)
          sj = xrj*rtpi2/sqrt(xx)
          sbsjy_hrmiha(il,im,ir) =  om2np1*sj*sy
          sbsjyp_hrmiha(il,im,ir) =  omega*om2np1*sj*syp
        endif
      enddo !ir
!   ir = 0 case: nonzero contribution only when l=m=0
!      
      sbsjy_hrmiha(0,0,0) = 1./rg(nrg)
      sbsjyp_hrmiha(0,0,0)= -1.d0/( rg(nrg)*rg(nrg) )
    enddo ! im 
  enddo ! il
!
end subroutine calc_radial_green_fn_hret_mi_hadv

