subroutine excurve_WL_gridpoint
  use phys_constant, only : nnrg, nntg, nnpg, pi
  use grid_parameter, only : nrg, ntg, npg
  use def_matter_parameter, only : ome, ber, radi
  use def_metric, only : alph, bvxd, bvyd, bvzd, bvxu, bvyu, bvzu
  use def_metric_excurve_grid, only : tfkij_grid, tfkijkij_grid, trk_grid
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
  use def_vector_x, only : vec_xg
!
  use interface_grgrad_4th_gridpoint
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
        call grgrad_4th_gridpoint(bvxd,dfdx,dfdy,dfdz,irg,itg,ipg)
        pdbvxd_grid(irg,itg,ipg,1) = dfdx
        pdbvxd_grid(irg,itg,ipg,2) = dfdy
        pdbvxd_grid(irg,itg,ipg,3) = dfdz
!
        call grgrad_4th_gridpoint(bvyd,dfdx,dfdy,dfdz,irg,itg,ipg)
        pdbvyd_grid(irg,itg,ipg,1) = dfdx
        pdbvyd_grid(irg,itg,ipg,2) = dfdy
        pdbvyd_grid(irg,itg,ipg,3) = dfdz
!
        call grgrad_4th_gridpoint(bvzd,dfdx,dfdy,dfdz,irg,itg,ipg)
        pdbvzd_grid(irg,itg,ipg,1) = dfdx
        pdbvzd_grid(irg,itg,ipg,2) = dfdy
        pdbvzd_grid(irg,itg,ipg,3) = dfdz
!
        call grgrad_4th_gridpoint(bvxu,pdbvxdx,pdbvxdy,pdbvxdz,irg,itg,ipg)
        call grgrad_4th_gridpoint(bvyu,pdbvydx,pdbvydy,pdbvydz,irg,itg,ipg)
        call grgrad_4th_gridpoint(bvzu,pdbvzdx,pdbvzdy,pdbvzdz,irg,itg,ipg)
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
        &                             - cri_grid(irg,itg,ipg,1,1)*bvxha &
        &                             - cri_grid(irg,itg,ipg,2,1)*bvyha &
        &                             - cri_grid(irg,itg,ipg,3,1)*bvzha
        cdbvyd_grid(irg,itg,ipg,1) = pdbvyd_grid(irg,itg,ipg,1) &
        &                             - cri_grid(irg,itg,ipg,1,2)*bvxha &
        &                             - cri_grid(irg,itg,ipg,2,2)*bvyha &
        &                             - cri_grid(irg,itg,ipg,3,2)*bvzha
        cdbvzd_grid(irg,itg,ipg,1) = pdbvzd_grid(irg,itg,ipg,1) &
        &                             - cri_grid(irg,itg,ipg,1,3)*bvxha &
        &                             - cri_grid(irg,itg,ipg,2,3)*bvyha &
        &                             - cri_grid(irg,itg,ipg,3,3)*bvzha
!
        cdbvxd_grid(irg,itg,ipg,2) = pdbvxd_grid(irg,itg,ipg,2) &
        &                             - cri_grid(irg,itg,ipg,1,2)*bvxha &
        &                             - cri_grid(irg,itg,ipg,2,2)*bvyha &
        &                             - cri_grid(irg,itg,ipg,3,2)*bvzha
        cdbvyd_grid(irg,itg,ipg,2) = pdbvyd_grid(irg,itg,ipg,2) &
        &                             - cri_grid(irg,itg,ipg,1,4)*bvxha &
        &                             - cri_grid(irg,itg,ipg,2,4)*bvyha &
        &                             - cri_grid(irg,itg,ipg,3,4)*bvzha
        cdbvzd_grid(irg,itg,ipg,2) = pdbvzd_grid(irg,itg,ipg,2) &
        &                             - cri_grid(irg,itg,ipg,1,5)*bvxha &
        &                             - cri_grid(irg,itg,ipg,2,5)*bvyha &
        &                             - cri_grid(irg,itg,ipg,3,5)*bvzha
!
        cdbvxd_grid(irg,itg,ipg,3) = pdbvxd_grid(irg,itg,ipg,3) &
        &                             - cri_grid(irg,itg,ipg,1,3)*bvxha &
        &                             - cri_grid(irg,itg,ipg,2,3)*bvyha &
        &                             - cri_grid(irg,itg,ipg,3,3)*bvzha
        cdbvyd_grid(irg,itg,ipg,3) = pdbvyd_grid(irg,itg,ipg,3) &
        &                             - cri_grid(irg,itg,ipg,1,5)*bvxha &
        &                             - cri_grid(irg,itg,ipg,2,5)*bvyha &
        &                             - cri_grid(irg,itg,ipg,3,5)*bvzha
        cdbvzd_grid(irg,itg,ipg,3) = pdbvzd_grid(irg,itg,ipg,3) &
        &                             - cri_grid(irg,itg,ipg,1,6)*bvxha &
        &                             - cri_grid(irg,itg,ipg,2,6)*bvyha &
        &                             - cri_grid(irg,itg,ipg,3,6)*bvzha
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
!       better for Aij
        cdivbv_grid(irg,itg,ipg) = gmxxu*dbvxdx + gmxyu*dbvydx + gmxzu*dbvzdx &
        &                        + gmyxu*dbvxdy + gmyyu*dbvydy + gmyzu*dbvzdy &
        &                        + gmzxu*dbvxdz + gmzyu*dbvydz + gmzzu*dbvzdz
        diver = fa23*cdivbv_grid(irg,itg,ipg)
!       better for traceK
        cdivbv_grid(irg,itg,ipg) = pdbvxdx + pdbvydy + pdbvzdz
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
        trk_grid(irg,itg,ipg) = 0.0d0
        do id = 1, 3
          do ic = 1, 3
            do ib = 1, 3
              do ia = 1, 3
                tfkijkij_grid(irg,itg,ipg) = tfkijkij_grid(irg,itg,ipg) &
                &                          + gamu(ia,ic)*gamu(ib,id) &
                &                          * tfkij_grid(irg,itg,ipg,ia,ib) &
                &                          * tfkij_grid(irg,itg,ipg,ic,id)
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
!!irg = 2; itg = ntg/2-2; ipg = 2
!!write(6,'(1p,3e16.8)') vec_xg(irg,itg,ipg,1), vec_xg(irg,itg,ipg,2), &
!!&                      vec_xg(irg,itg,ipg,3)
!!write(6,'(1p,3e16.8)') tfkij_grid(irg,itg,ipg,1,1), &
!!&                      tfkij_grid(irg,itg,ipg,1,2), &
!!&                      tfkij_grid(irg,itg,ipg,1,3)
!!write(6,'(1p,3e16.8)') tfkij_grid(irg,itg,ipg,2,1), &
!!&                      tfkij_grid(irg,itg,ipg,2,2), &
!!&                      tfkij_grid(irg,itg,ipg,2,3)
!!write(6,'(1p,3e16.8)') tfkij_grid(irg,itg,ipg,3,1), &
!!&                      tfkij_grid(irg,itg,ipg,3,2), &
!!&                      tfkij_grid(irg,itg,ipg,3,3)
!!stop
end subroutine excurve_WL_gridpoint
