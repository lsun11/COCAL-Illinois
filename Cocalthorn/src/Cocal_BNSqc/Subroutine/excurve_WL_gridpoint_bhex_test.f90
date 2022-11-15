subroutine excurve_WL_gridpoint_bhex
  use phys_constant, only : nnrg, nntg, nnpg, pi
  use grid_parameter, only : nrg, ntg, npg
  use def_matter_parameter, only : ome, ber, radi
  use def_metric, only : alph, bvxd, bvyd, bvzd, bvxu, bvyu, bvzu
  use def_metric_hij, only : hxxu, hxyu, hxzu, hyyu, hyzu, hzzu  
  use def_kerr_schild, only: Fxu_grid_ks, Fyu_grid_ks, Fzu_grid_ks
  use def_metric_excurve_grid, only : tfkij_grid, tfkijkij_grid, Lij_grid
  use coordinate_grav_r, only : rg
!
  use def_cristoffel_grid, only : cri_grid
  use def_Lie_derivatives_grid, only : elpxx_grid, elpxy_grid, elpxz_grid, &
  &                                    elpyy_grid, elpyz_grid, elpzz_grid, &
  &                                    rlbxx_grid, rlbxy_grid, rlbxz_grid, &
  &                                    rlbyy_grid, rlbyz_grid, rlbzz_grid
  use def_shift_derivatives_grid, only : cdbvxd_grid, cdbvyd_grid, &
  &                                      cdbvzd_grid, cdivbv_grid, &
  &                                      pdbvxd_grid, pdbvyd_grid, pdbvzd_grid
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_cutsw, only : cutfac
  use def_formulation, only : chgra
!
  use interface_grgrad_4th_gridpoint_bhex
  use make_array_3d
  use make_array_4d
  implicit none
!
  real(8) :: dfdx, dfdy, dfdz
  real(8) :: gamu(3,3)
  real(8) :: ainvh, bvxha, bvyha, bvzha, cutoff, diver, fa23, &
  &          dbvxdx, dbvxdy, dbvxdz, &
  &          dbvydx, dbvydy, dbvydz, &
  &          dbvzdx, dbvzdy, dbvzdz, &
  &          gmxxd, gmxyd, gmxzd, gmxxu, gmxyu, gmxzu, &
  &          gmyxd, gmyyd, gmyzd, gmyxu, gmyyu, gmyzu, &
  &          gmzxd, gmzyd, gmzzd, gmzxu, gmzyu, gmzzu, &
  &          oelpxx, oelpxy, oelpxz, &
  &          oelpyy, oelpyz, oelpzz, &
  &          pdbvxdx, pdbvxdy, pdbvxdz, &
  &          pdbvydx, pdbvydy, pdbvydz, &
  &          pdbvzdx, pdbvzdy, pdbvzdz
  real(8) :: pdhxxdx,pdhxxdy,pdhxxdz, pdhxydx,pdhxydy,pdhxydz, &
     &       pdhxzdx,pdhxzdy,pdhxzdz, pdhyydx,pdhyydy,pdhyydz, &
     &       pdhyzdx,pdhyzdy,pdhyzdz, pdhzzdx,pdhzzdy,pdhzzdz, &
     &       bxu,byu,bzu,txx,tyy,tzz,txy,txz,tyz,bdhxx,bdhxy,bdhxz, &
     &       bdhyy, bdhyz, bdhzz
  integer :: ia, ib, ic, id, info, ipg, irg, itg
!
! --- Compute extringic curvature.  
! --- Whose value is assigned on the grid points. 
!
!testtesttest
  info = 0
!testtesttest
!
  fa23 = 2.0d0/3.0d0
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
!
        call grgrad_4th_gridpoint_bhex(bvxd,dfdx,dfdy,dfdz,irg,itg,ipg)
        pdbvxd_grid(irg,itg,ipg,1) = dfdx
        pdbvxd_grid(irg,itg,ipg,2) = dfdy
        pdbvxd_grid(irg,itg,ipg,3) = dfdz
!
        call grgrad_4th_gridpoint_bhex(bvyd,dfdx,dfdy,dfdz,irg,itg,ipg)
        pdbvyd_grid(irg,itg,ipg,1) = dfdx
        pdbvyd_grid(irg,itg,ipg,2) = dfdy
        pdbvyd_grid(irg,itg,ipg,3) = dfdz
!
        call grgrad_4th_gridpoint_bhex(bvzd,dfdx,dfdy,dfdz,irg,itg,ipg)
        pdbvzd_grid(irg,itg,ipg,1) = dfdx
        pdbvzd_grid(irg,itg,ipg,2) = dfdy
        pdbvzd_grid(irg,itg,ipg,3) = dfdz
!
        call grgrad_4th_gridpoint_bhex(bvxu,pdbvxdx,pdbvxdy,pdbvxdz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(bvyu,pdbvydx,pdbvydy,pdbvydz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(bvzu,pdbvzdx,pdbvzdy,pdbvzdz,irg,itg,ipg)
!
        cutoff = 1.0d0
        if (chgra == 'c'.or.chgra == 'C'.or.chgra == 'W') then
          if (rg(irg) > cutfac*pi/ome) cutoff = 0.0d0
        end if
!
        bvxha = bvxd(irg,itg,ipg)
        bvyha = bvyd(irg,itg,ipg)
        bvzha = bvzd(irg,itg,ipg)
        ainvh = alph(irg,itg,ipg)
        ainvh = 0.5d0/ainvh
!
        cdbvxd_grid(irg,itg,ipg,1) = pdbvxd_grid(irg,itg,ipg,1) &
        &                      - cri_grid(irg,itg,ipg,1,1)*bvxha &
        &                      - cri_grid(irg,itg,ipg,2,1)*bvyha &
        &                      - cri_grid(irg,itg,ipg,3,1)*bvzha
        cdbvyd_grid(irg,itg,ipg,1) = pdbvyd_grid(irg,itg,ipg,1) &
        &                      - cri_grid(irg,itg,ipg,1,2)*bvxha &
        &                      - cri_grid(irg,itg,ipg,2,2)*bvyha &
        &                      - cri_grid(irg,itg,ipg,3,2)*bvzha
        cdbvzd_grid(irg,itg,ipg,1) = pdbvzd_grid(irg,itg,ipg,1) &
        &                      - cri_grid(irg,itg,ipg,1,3)*bvxha &
        &                      - cri_grid(irg,itg,ipg,2,3)*bvyha &
        &                      - cri_grid(irg,itg,ipg,3,3)*bvzha
!
        cdbvxd_grid(irg,itg,ipg,2) = pdbvxd_grid(irg,itg,ipg,2) &
        &                      - cri_grid(irg,itg,ipg,1,2)*bvxha &
        &                      - cri_grid(irg,itg,ipg,2,2)*bvyha &
        &                      - cri_grid(irg,itg,ipg,3,2)*bvzha
        cdbvyd_grid(irg,itg,ipg,2) = pdbvyd_grid(irg,itg,ipg,2) &
        &                      - cri_grid(irg,itg,ipg,1,4)*bvxha &
        &                      - cri_grid(irg,itg,ipg,2,4)*bvyha &
        &                      - cri_grid(irg,itg,ipg,3,4)*bvzha
        cdbvzd_grid(irg,itg,ipg,2) = pdbvzd_grid(irg,itg,ipg,2) &
        &                      - cri_grid(irg,itg,ipg,1,5)*bvxha &
        &                      - cri_grid(irg,itg,ipg,2,5)*bvyha &
        &                      - cri_grid(irg,itg,ipg,3,5)*bvzha
!
        cdbvxd_grid(irg,itg,ipg,3) = pdbvxd_grid(irg,itg,ipg,3) &
        &                      - cri_grid(irg,itg,ipg,1,3)*bvxha &
        &                      - cri_grid(irg,itg,ipg,2,3)*bvyha &
        &                      - cri_grid(irg,itg,ipg,3,3)*bvzha
        cdbvyd_grid(irg,itg,ipg,3) = pdbvyd_grid(irg,itg,ipg,3) &
        &                      - cri_grid(irg,itg,ipg,1,5)*bvxha &
        &                      - cri_grid(irg,itg,ipg,2,5)*bvyha &
        &                      - cri_grid(irg,itg,ipg,3,5)*bvzha
        cdbvzd_grid(irg,itg,ipg,3) = pdbvzd_grid(irg,itg,ipg,3) &
        &                      - cri_grid(irg,itg,ipg,1,6)*bvxha &
        &                      - cri_grid(irg,itg,ipg,2,6)*bvyha &
        &                      - cri_grid(irg,itg,ipg,3,6)*bvzha
!
        dbvxdx = cdbvxd_grid(irg,itg,ipg,1)
        dbvydx = cdbvyd_grid(irg,itg,ipg,1)
        dbvzdx = cdbvzd_grid(irg,itg,ipg,1)
!
        dbvxdy = cdbvxd_grid(irg,itg,ipg,2)
        dbvydy = cdbvyd_grid(irg,itg,ipg,2)
        dbvzdy = cdbvzd_grid(irg,itg,ipg,2)
!
        dbvxdz = cdbvxd_grid(irg,itg,ipg,3)
        dbvydz = cdbvyd_grid(irg,itg,ipg,3)
        dbvzdz = cdbvzd_grid(irg,itg,ipg,3)
!
        gmxxd = hxxd(irg,itg,ipg)
        gmxyd = hxyd(irg,itg,ipg)
        gmxzd = hxzd(irg,itg,ipg)
        gmyyd = hyyd(irg,itg,ipg)
        gmyzd = hyzd(irg,itg,ipg)
        gmzzd = hzzd(irg,itg,ipg)
        gmxxd = gmxxd + 1.0d0
        gmyyd = gmyyd + 1.0d0
        gmzzd = gmzzd + 1.0d0
        gmyxd = gmxyd
        gmzxd = gmxzd
        gmzyd = gmyzd
        gmxxu = hxxu(irg,itg,ipg)
        gmxyu = hxyu(irg,itg,ipg)
        gmxzu = hxzu(irg,itg,ipg)
        gmyyu = hyyu(irg,itg,ipg)
        gmyzu = hyzu(irg,itg,ipg)
        gmzzu = hzzu(irg,itg,ipg)
        gmyyu = gmyyu + 1.0d0
        gmxxu = gmxxu + 1.0d0
        gmzzu = gmzzu + 1.0d0
        gmyxu = gmxyu
        gmzxu = gmxzu
        gmzyu = gmyzu
!
!ccp      cdivbv(irg,itg,ipg) = gmxxu*dbvxdx + gmxyu*dbvydx + gmxzu*dbvzdx
!ccp     &                    + gmyxu*dbvxdy + gmyyu*dbvydy + gmyzu*dbvzdy
!ccp     &                    + gmzxu*dbvxdz + gmzyu*dbvydz + gmzzu*dbvzdz
        cdivbv_grid(irg,itg,ipg) = pdbvxdx + pdbvydy + pdbvzdz
        diver = fa23*cdivbv_grid(irg,itg,ipg)
!
! --  For rotating shift
!      
        oelpxx = 0.0d0
        oelpxy = 0.0d0
        oelpxz = 0.0d0
        oelpyy = 0.0d0
        oelpyz = 0.0d0
        oelpzz = 0.0d0
        if (chgra == 'h'.or.chgra == 'c'.or.chgra == 'C' &
           &.or.chgra == 'H'.or.chgra == 'W') then
          oelpxx = ainvh*ome*elpxx_grid(irg,itg,ipg)*cutoff
          oelpxy = ainvh*ome*elpxy_grid(irg,itg,ipg)*cutoff
          oelpxz = ainvh*ome*elpxz_grid(irg,itg,ipg)*cutoff
          oelpyy = ainvh*ome*elpyy_grid(irg,itg,ipg)*cutoff
          oelpyz = ainvh*ome*elpyz_grid(irg,itg,ipg)*cutoff
          oelpzz = ainvh*ome*elpzz_grid(irg,itg,ipg)*cutoff
        end if
!
        tfkij_grid(irg,itg,ipg,1,1) = ainvh*(2.0d0*dbvxdx-gmxxd*diver) +oelpxx
        tfkij_grid(irg,itg,ipg,2,2) = ainvh*(2.0d0*dbvydy-gmyyd*diver) +oelpyy
        tfkij_grid(irg,itg,ipg,3,3) = ainvh*(2.0d0*dbvzdz-gmzzd*diver) +oelpzz
        tfkij_grid(irg,itg,ipg,1,2) = ainvh*(dbvydx+dbvxdy-gmxyd*diver)+oelpxy
        tfkij_grid(irg,itg,ipg,1,3) = ainvh*(dbvzdx+dbvxdz-gmxzd*diver)+oelpxz
        tfkij_grid(irg,itg,ipg,2,3) = ainvh*(dbvzdy+dbvydz-gmyzd*diver)+oelpyz
        tfkij_grid(irg,itg,ipg,2,1) = tfkij_grid(irg,itg,ipg,1,2)
        tfkij_grid(irg,itg,ipg,3,1) = tfkij_grid(irg,itg,ipg,1,3)
        tfkij_grid(irg,itg,ipg,3,2) = tfkij_grid(irg,itg,ipg,2,3)
!


        call grgrad_4th_gridpoint_bhex(hxxu,pdhxxdx,pdhxxdy,pdhxxdz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(hxyu,pdhxydx,pdhxydy,pdhxydz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(hxzu,pdhxzdx,pdhxzdy,pdhxzdz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(hyyu,pdhyydx,pdhyydy,pdhyydz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(hyzu,pdhyzdx,pdhyzdy,pdhyzdz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(hzzu,pdhzzdx,pdhzzdy,pdhzzdz,irg,itg,ipg)

!        pdhxxdx = Fxu_grid_ks(irg,itg,ipg) - pdhxydy - pdhxzdz
!        pdhyydy = Fyu_grid_ks(irg,itg,ipg) - pdhxydx - pdhyzdz
!        pdhzzdz = Fzu_grid_ks(irg,itg,ipg) - pdhxzdx - pdhyzdy

        bxu = bvxu(irg,itg,ipg)
        byu = bvyu(irg,itg,ipg)
        bzu = bvzu(irg,itg,ipg)

        bdhxx = bxu*pdhxxdx + byu*pdhxxdy + bzu*pdhxxdz
        bdhyy = bxu*pdhyydx + byu*pdhyydy + bzu*pdhyydz 
        bdhzz = bxu*pdhzzdx + byu*pdhzzdy + bzu*pdhzzdz
        bdhxy = bxu*pdhxydx + byu*pdhxydy + bzu*pdhxydz
        bdhxz = bxu*pdhxzdx + byu*pdhxzdy + bzu*pdhxzdz
        bdhyz = bxu*pdhyzdx + byu*pdhyzdy + bzu*pdhyzdz


        txx = gmxxd**2*bdhxx + gmxyd**2*bdhyy + gmxzd**2*bdhzz   &
         &  + 2.0d0*gmxxd*gmxyd*bdhxy    &
         &  + 2.0d0*gmxxd*gmxzd*bdhxz    &
         &  + 2.0d0*gmxyd*gmxzd*bdhyz    

        tyy = gmyxd**2*bdhxx + gmyyd**2*bdhyy + gmyzd**2*bdhzz     &
         &  + 2.0d0*gmyyd*gmxyd*bdhxy    &
         &  + 2.0d0*gmyxd*gmyzd*bdhxz    &
         &  + 2.0d0*gmyyd*gmyzd*bdhyz     

        tzz = gmzxd**2*bdhxx + gmzyd**2*bdhyy + gmzzd**2*bdhzz     &
         &  + 2.0d0*gmzxd*gmzyd*bdhxy    &
         &  + 2.0d0*gmzxd*gmzzd*bdhxz    &
         &  + 2.0d0*gmzyd*gmzzd*bdhyz      

        txy = gmxxd*gmxyd*bdhxx + gmxyd*gmyyd*bdhyy + gmxzd*gmyzd*bdhzz    &
         &  + (gmxyd*gmyxd+gmxxd*gmyyd)*bdhxy    &
         &  + (gmxzd*gmyxd+gmxxd*gmyzd)*bdhxz    &
         &  + (gmxzd*gmyyd+gmxyd*gmyzd)*bdhyz     

        txz = gmxxd*gmxzd*bdhxx + gmxyd*gmzyd*bdhyy + gmxzd*gmzzd*bdhzz    &
         &  + (gmxyd*gmzxd+gmxxd*gmzyd)*bdhxy    &
         &  + (gmxzd*gmzxd+gmxxd*gmzzd)*bdhxz    &
         &  + (gmxzd*gmzyd+gmxyd*gmzzd)*bdhyz    

        tyz = gmyxd*gmzxd*bdhxx + gmyyd*gmzyd*bdhyy + gmyzd*gmzzd*bdhzz  &
         &  + (gmyyd*gmzxd+gmyxd*gmzyd)*bdhxy     &
         &  + (gmyzd*gmzxd+gmyxd*gmzzd)*bdhxz     &
         &  + (gmyzd*gmzyd+gmyyd*gmzzd)*bdhyz      



        Lij_grid(irg,itg,ipg,1,1) = ainvh*(gmxxd*pdbvxdx + gmxyd*pdbvydx + gmxzd*pdbvzdx &
                                  &      + gmxxd*pdbvxdx + gmxyd*pdbvydx + gmxzd*pdbvzdx - gmxxd*diver - txx) 

        Lij_grid(irg,itg,ipg,2,2) = ainvh*(gmyxd*pdbvxdy + gmyyd*pdbvydy + gmyzd*pdbvzdy &
                                  &      + gmyxd*pdbvxdy + gmyyd*pdbvydy + gmyzd*pdbvzdy - gmyyd*diver - tyy)

        Lij_grid(irg,itg,ipg,3,3) = ainvh*(gmzxd*pdbvxdz + gmzyd*pdbvydz + gmzzd*pdbvzdz &
                                  &      + gmzxd*pdbvxdz + gmzyd*pdbvydz + gmzzd*pdbvzdz - gmzzd*diver - tzz)

        Lij_grid(irg,itg,ipg,1,2) = ainvh*(gmxxd*pdbvxdy + gmxyd*pdbvydy + gmxzd*pdbvzdy &
                                  &      + gmyxd*pdbvxdx + gmyyd*pdbvydx + gmyzd*pdbvzdx - gmxyd*diver - txy)

        Lij_grid(irg,itg,ipg,1,3) = ainvh*(gmxxd*pdbvxdz + gmxyd*pdbvydz + gmxzd*pdbvzdz &
                                  &      + gmzxd*pdbvxdx + gmzyd*pdbvydx + gmzzd*pdbvzdx - gmxzd*diver - txz)

        Lij_grid(irg,itg,ipg,2,3) = ainvh*(gmyxd*pdbvxdz + gmyyd*pdbvydz + gmyzd*pdbvzdz &
                                  &      + gmzxd*pdbvxdy + gmzyd*pdbvydy + gmzzd*pdbvzdy - gmyzd*diver - tyz)

        Lij_grid(irg,itg,ipg,2,1) = Lij_grid(irg,itg,ipg,1,2)
        Lij_grid(irg,itg,ipg,3,1) = Lij_grid(irg,itg,ipg,1,3)
        Lij_grid(irg,itg,ipg,3,2) = Lij_grid(irg,itg,ipg,2,3)
!


        gamu(1,1) = gmxxu 
        gamu(1,2) = gmxyu 
        gamu(1,3) = gmxzu 
        gamu(2,1) = gmyxu 
        gamu(2,2) = gmyyu 
        gamu(2,3) = gmyzu 
        gamu(3,1) = gmzxu 
        gamu(3,2) = gmzyu 
        gamu(3,3) = gmzzu
! 
        tfkijkij_grid(irg,itg,ipg) = 0.0d0
        do id = 1, 3
          do ic = 1, 3
            do ib = 1, 3
              do ia = 1, 3
                tfkijkij_grid(irg,itg,ipg) = tfkijkij_grid(irg,itg,ipg) &
                   &       + gamu(ia,ic)*gamu(ib,id) &
                   &        *tfkij_grid(irg,itg,ipg,ia,ib)*tfkij_grid(irg,itg,ipg,ic,id)
              end do
            end do
          end do
        end do
!
        if (tfkijkij_grid(irg,itg,ipg) /= 0.) info = 1
!
! --  Lie_beta tgamma
!
        rlbxx_grid(irg,itg,ipg) = 2.0d0*dbvxdx
        rlbxy_grid(irg,itg,ipg) = dbvydx + dbvxdy
        rlbxz_grid(irg,itg,ipg) = dbvzdx + dbvxdz
        rlbyy_grid(irg,itg,ipg) = 2.0d0*dbvydy
        rlbyz_grid(irg,itg,ipg) = dbvzdy + dbvydz
        rlbzz_grid(irg,itg,ipg) = 2.0d0*dbvzdz
!
      end do
    end do
  end do
  if (info /= 1) write(6,*) ' ### Warning K_ij = 0 *** '
!
end subroutine excurve_WL_gridpoint_bhex
