subroutine calc_radial_green_fn_hrethadv_homosol(omega,char_dn)
  use coordinate_grav_r, only : rg, hrg
  use phys_constant, only : long, pi
  use grid_parameter, only : nrg, npg, ntg, nlg
  use radial_green_fn_hrethadv_homosol
  implicit none
  integer :: il, im, ir, irr
  real(8) :: omega, omegam, omegam_rg, omegam_rgin, l2p1, l2p1dlp1, W, Wp
  real(8) :: rtpi2 = 1.253314154753822d0
  real(8) :: sbj, sbj_in, sbjp_in, sj, sjp
  real(8) :: sby, sby_in, sbyp_in, sy, syp
  character(len=2) :: char_dn
!
!-------------------------------------
!--- Computation of radial green function for homogeneous sol.
!    char_dn = 'dh', dirichlet needs minus sign
!    char_dn = 'nh', neumann
!-------------------------------------
!
!
! For the homogeneous solution, for m nonzero
!  sbsjy_hrha_ho(l,m,ir) :=  j_l(m omega r)/W
!  sbsjyp_hrha_ho(l,m,ir) =  j'_l(m omega r)/W'
!  W  = sqrt(j_l(m omega R)^2  + y_l(m omega R)^2)
!  W' = sqrt(j'_l(m omega R)^2 + y'_l(m omega R)^2)
!
! For m = 0, use 
!  sbsjy_hrha(l,0,ir)  = r^l/R^(l+1)
!  sbsjyp_hrha(l,0,ir) = - (l+1) r^l/R^(l+2)
!
  sbsjy_hrha_ho  = 0.0d0
  sbsjyp_hrha_ho = 0.0d0
  do il = 0, nlg
    do im = 0, il
      do ir = 0, nrg
        if(im.eq.0)then
          l2p1 = dble(2*il+1)
          l2p1dlp1 = l2p1/dble(il+1)
          if (char_dn.eq.'nh') then 
            sbsjy_hrha_ho(il,im,ir) = - l2p1dlp1*rg(0)**(il+2)/rg(ir)**(il+1)
          end if  
          if (char_dn.eq.'dh') then
            sbsjyp_hrha_ho(il,im,ir)= - l2p1*(rg(0)/rg(ir))**(il+1)
          end if  
        else
          omegam = dble(im)*omega
          omegam_rgin = omegam*rg(0)
          l2p1 = dble(2*il+1)
          call sphbess_and_dx(omegam_rgin,il,sj,sy,sjp,syp)
          sby_in = sy
          sbj_in = sj
          sbyp_in = omegam*syp
          sbjp_in = omegam*sjp
          omegam_rg = omegam*rg(ir)
          W  = sqrt(sbj_in**2 + sby_in**2)
          Wp = sqrt(sbjp_in**2 + sbyp_in**2)
          call sphbess(omegam_rg,il,sj,sy)
          sby = sy
          if (char_dn.eq.'nh') then 
!            sbsjy_hrha_ho(il,im,ir)  = l2p1*sby/sbyp_in ! original
            sbsjy_hrha_ho(il,im,ir)  = l2p1*sby/Wp
          end if
          if (char_dn.eq.'dh') then
!            sbsjyp_hrha_ho(il,im,ir) = - l2p1*sby/sby_in  ! original
            sbsjyp_hrha_ho(il,im,ir) = l2p1*sby/W
          end if
        endif
      enddo !ir
    enddo ! im 
  enddo ! il
!
end subroutine calc_radial_green_fn_hrethadv_homosol
