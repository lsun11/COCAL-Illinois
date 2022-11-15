subroutine excurve_WL(chgra)
  use phys_constant, only : nnrg, nntg, nnpg, pi
  use grid_parameter, only : nrg, ntg, npg
  use def_matter_parameter, only : ome, ber, radi
  use def_metric, only : alph, tfkij, tfkijkij, trk, &
  &                      bvxd, bvyd, bvzd, bvxu, bvyu, bvzu
  use coordinate_grav_r, only : rg
!
  use def_cristoffel, only : cri
  use def_Lie_derivatives, only : elpxx, elpxy, elpxz, elpyy, elpyz, elpzz, &
  &                               rlbxx, rlbxy, rlbxz, rlbyy, rlbyz, rlbzz
  use def_shift_derivatives, only : cdbvxd, cdbvyd, cdbvzd, cdivbv, &
  &                                 pdbvxd, pdbvyd, pdbvzd
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_cutsw, only : cutfac
!
  use interface_interpo_linear_type0
  use interface_grgrad_midpoint
  use make_array_3d
  use make_array_4d
  implicit none
!
  real(8), pointer :: dbvxu(:,:,:,:), dbvyu(:,:,:,:), dbvzu(:,:,:,:)
  real(8), pointer :: dfdx(:,:,:), dfdy(:,:,:), dfdz(:,:,:)
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
  character(len=1), intent(in) :: chgra
!
  call alloc_array3d(dfdx,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfdy,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfdz,1,nrg,1,ntg,1,npg)
  call alloc_array4d(dbvxu,1,nrg,1,ntg,1,npg,1,3)
  call alloc_array4d(dbvyu,1,nrg,1,ntg,1,npg,1,3)
  call alloc_array4d(dbvzu,1,nrg,1,ntg,1,npg,1,3)
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
  call grgrad_midpoint(bvxd,dfdx,dfdy,dfdz)
  pdbvxd(1:nrg,1:ntg,1:npg,1) = dfdx(1:nrg,1:ntg,1:npg)
  pdbvxd(1:nrg,1:ntg,1:npg,2) = dfdy(1:nrg,1:ntg,1:npg)
  pdbvxd(1:nrg,1:ntg,1:npg,3) = dfdz(1:nrg,1:ntg,1:npg)
!
  call grgrad_midpoint(bvyd,dfdx,dfdy,dfdz)
  pdbvyd(1:nrg,1:ntg,1:npg,1) = dfdx(1:nrg,1:ntg,1:npg)
  pdbvyd(1:nrg,1:ntg,1:npg,2) = dfdy(1:nrg,1:ntg,1:npg)
  pdbvyd(1:nrg,1:ntg,1:npg,3) = dfdz(1:nrg,1:ntg,1:npg)
!
  call grgrad_midpoint(bvzd,dfdx,dfdy,dfdz)
  pdbvzd(1:nrg,1:ntg,1:npg,1) = dfdx(1:nrg,1:ntg,1:npg)
  pdbvzd(1:nrg,1:ntg,1:npg,2) = dfdy(1:nrg,1:ntg,1:npg)
  pdbvzd(1:nrg,1:ntg,1:npg,3) = dfdz(1:nrg,1:ntg,1:npg)
!
  call grgrad_midpoint(bvxu,dfdx,dfdy,dfdz)
  dbvxu(1:nrg,1:ntg,1:npg,1) = dfdx(1:nrg,1:ntg,1:npg)
  dbvxu(1:nrg,1:ntg,1:npg,2) = dfdy(1:nrg,1:ntg,1:npg)
  dbvxu(1:nrg,1:ntg,1:npg,3) = dfdz(1:nrg,1:ntg,1:npg)
!
  call grgrad_midpoint(bvyu,dfdx,dfdy,dfdz)
  dbvyu(1:nrg,1:ntg,1:npg,1) = dfdx(1:nrg,1:ntg,1:npg)
  dbvyu(1:nrg,1:ntg,1:npg,2) = dfdy(1:nrg,1:ntg,1:npg)
  dbvyu(1:nrg,1:ntg,1:npg,3) = dfdz(1:nrg,1:ntg,1:npg)
!
  call grgrad_midpoint(bvzu,dfdx,dfdy,dfdz)
  dbvzu(1:nrg,1:ntg,1:npg,1) = dfdx(1:nrg,1:ntg,1:npg)
  dbvzu(1:nrg,1:ntg,1:npg,2) = dfdy(1:nrg,1:ntg,1:npg)
  dbvzu(1:nrg,1:ntg,1:npg,3) = dfdz(1:nrg,1:ntg,1:npg)
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        cutoff = 1.0d0
        if (chgra == 'c'.or.chgra == 'C'.or.chgra == 'W') then
          if (rg(irg) > cutfac*pi/ome) cutoff = 0.0d0
        end if
!
        call interpo_linear_type0(bvxha,bvxd,irg,itg,ipg)
        call interpo_linear_type0(bvyha,bvyd,irg,itg,ipg)
        call interpo_linear_type0(bvzha,bvzd,irg,itg,ipg)
        call interpo_linear_type0(ainvh,alph,irg,itg,ipg)
        ainvh = 0.5d0/ainvh
! 
        pdbvxdx = dbvxu(irg,itg,ipg,1)
        pdbvxdy = dbvxu(irg,itg,ipg,2)
        pdbvxdz = dbvxu(irg,itg,ipg,3)
        pdbvydx = dbvyu(irg,itg,ipg,1)
        pdbvydy = dbvyu(irg,itg,ipg,2)
        pdbvydz = dbvyu(irg,itg,ipg,3)
        pdbvzdx = dbvzu(irg,itg,ipg,1)
        pdbvzdy = dbvzu(irg,itg,ipg,2)
        pdbvzdz = dbvzu(irg,itg,ipg,3)
!
        cdbvxd(irg,itg,ipg,1) = pdbvxd(irg,itg,ipg,1) &
           &                      - cri(irg,itg,ipg,1,1)*bvxha &
           &                      - cri(irg,itg,ipg,2,1)*bvyha &
           &                      - cri(irg,itg,ipg,3,1)*bvzha
        cdbvyd(irg,itg,ipg,1) = pdbvyd(irg,itg,ipg,1) &
           &                      - cri(irg,itg,ipg,1,2)*bvxha &
           &                      - cri(irg,itg,ipg,2,2)*bvyha &
           &                      - cri(irg,itg,ipg,3,2)*bvzha
        cdbvzd(irg,itg,ipg,1) = pdbvzd(irg,itg,ipg,1) &
           &                      - cri(irg,itg,ipg,1,3)*bvxha &
           &                      - cri(irg,itg,ipg,2,3)*bvyha &
           &                      - cri(irg,itg,ipg,3,3)*bvzha
!
        cdbvxd(irg,itg,ipg,2) = pdbvxd(irg,itg,ipg,2) &
           &                      - cri(irg,itg,ipg,1,2)*bvxha &
           &                      - cri(irg,itg,ipg,2,2)*bvyha &
           &                      - cri(irg,itg,ipg,3,2)*bvzha
        cdbvyd(irg,itg,ipg,2) = pdbvyd(irg,itg,ipg,2) &
           &                      - cri(irg,itg,ipg,1,4)*bvxha &
           &                      - cri(irg,itg,ipg,2,4)*bvyha &
           &                      - cri(irg,itg,ipg,3,4)*bvzha
        cdbvzd(irg,itg,ipg,2) = pdbvzd(irg,itg,ipg,2) &
           &                      - cri(irg,itg,ipg,1,5)*bvxha &
           &                      - cri(irg,itg,ipg,2,5)*bvyha &
           &                      - cri(irg,itg,ipg,3,5)*bvzha
!
        cdbvxd(irg,itg,ipg,3) = pdbvxd(irg,itg,ipg,3) &
           &                      - cri(irg,itg,ipg,1,3)*bvxha &
           &                      - cri(irg,itg,ipg,2,3)*bvyha &
           &                      - cri(irg,itg,ipg,3,3)*bvzha
        cdbvyd(irg,itg,ipg,3) = pdbvyd(irg,itg,ipg,3) &
           &                      - cri(irg,itg,ipg,1,5)*bvxha &
           &                      - cri(irg,itg,ipg,2,5)*bvyha &
           &                      - cri(irg,itg,ipg,3,5)*bvzha
        cdbvzd(irg,itg,ipg,3) = pdbvzd(irg,itg,ipg,3) &
           &                      - cri(irg,itg,ipg,1,6)*bvxha &
           &                      - cri(irg,itg,ipg,2,6)*bvyha &
           &                      - cri(irg,itg,ipg,3,6)*bvzha
!
        dbvxdx = cdbvxd(irg,itg,ipg,1)
        dbvydx = cdbvyd(irg,itg,ipg,1)
        dbvzdx = cdbvzd(irg,itg,ipg,1)
!
        dbvxdy = cdbvxd(irg,itg,ipg,2)
        dbvydy = cdbvyd(irg,itg,ipg,2)
        dbvzdy = cdbvzd(irg,itg,ipg,2)
!
        dbvxdz = cdbvxd(irg,itg,ipg,3)
        dbvydz = cdbvyd(irg,itg,ipg,3)
        dbvzdz = cdbvzd(irg,itg,ipg,3)
!
        call interpo_linear_type0(gmxxd,hxxd,irg,itg,ipg)
        call interpo_linear_type0(gmxyd,hxyd,irg,itg,ipg)
        call interpo_linear_type0(gmxzd,hxzd,irg,itg,ipg)
        call interpo_linear_type0(gmyyd,hyyd,irg,itg,ipg)
        call interpo_linear_type0(gmyzd,hyzd,irg,itg,ipg)
        call interpo_linear_type0(gmzzd,hzzd,irg,itg,ipg)
        gmxxd = gmxxd + 1.0d0
        gmyyd = gmyyd + 1.0d0
        gmzzd = gmzzd + 1.0d0
        gmyxd = gmxyd
        gmzxd = gmxzd
        gmzyd = gmyzd
        call interpo_linear_type0(gmxxu,hxxu,irg,itg,ipg)
        call interpo_linear_type0(gmxyu,hxyu,irg,itg,ipg)
        call interpo_linear_type0(gmxzu,hxzu,irg,itg,ipg)
        call interpo_linear_type0(gmyyu,hyyu,irg,itg,ipg)
        call interpo_linear_type0(gmyzu,hyzu,irg,itg,ipg)
        call interpo_linear_type0(gmzzu,hzzu,irg,itg,ipg)
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
        cdivbv(irg,itg,ipg) = pdbvxdx + pdbvydy + pdbvzdz
        diver = fa23*cdivbv(irg,itg,ipg)
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
          oelpxx = ainvh*ome*elpxx(irg,itg,ipg)*cutoff
          oelpxy = ainvh*ome*elpxy(irg,itg,ipg)*cutoff
          oelpxz = ainvh*ome*elpxz(irg,itg,ipg)*cutoff
          oelpyy = ainvh*ome*elpyy(irg,itg,ipg)*cutoff
          oelpyz = ainvh*ome*elpyz(irg,itg,ipg)*cutoff
          oelpzz = ainvh*ome*elpzz(irg,itg,ipg)*cutoff
        end if
!
        tfkij(irg,itg,ipg,1,1) = ainvh*(2.0d0*dbvxdx-gmxxd*diver) +oelpxx
        tfkij(irg,itg,ipg,2,2) = ainvh*(2.0d0*dbvydy-gmyyd*diver) +oelpyy
        tfkij(irg,itg,ipg,3,3) = ainvh*(2.0d0*dbvzdz-gmzzd*diver) +oelpzz
        tfkij(irg,itg,ipg,1,2) = ainvh*(dbvydx+dbvxdy-gmxyd*diver)+oelpxy
        tfkij(irg,itg,ipg,1,3) = ainvh*(dbvzdx+dbvxdz-gmxzd*diver)+oelpxz
        tfkij(irg,itg,ipg,2,3) = ainvh*(dbvzdy+dbvydz-gmyzd*diver)+oelpyz
        tfkij(irg,itg,ipg,2,1) = tfkij(irg,itg,ipg,1,2)
        tfkij(irg,itg,ipg,3,1) = tfkij(irg,itg,ipg,1,3)
        tfkij(irg,itg,ipg,3,2) = tfkij(irg,itg,ipg,2,3)
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
        tfkijkij(irg,itg,ipg) = 0.0d0
        trk(irg,itg,ipg) = 0.0d0
        do id = 1, 3
          do ic = 1, 3
            do ib = 1, 3
              do ia = 1, 3
                tfkijkij(irg,itg,ipg) = tfkijkij(irg,itg,ipg) &
                   &       + gamu(ia,ic)*gamu(ib,id) &
                   &        *tfkij(irg,itg,ipg,ia,ib)*tfkij(irg,itg,ipg,ic,id)
              end do
            end do
          end do
        end do
!
        if (tfkijkij(irg,itg,ipg) /= 0.) info = 1
!
! --  Lie_beta tgamma
!
        rlbxx(irg,itg,ipg) = 2.0d0*dbvxdx
        rlbxy(irg,itg,ipg) = dbvydx + dbvxdy
        rlbxz(irg,itg,ipg) = dbvzdx + dbvxdz
        rlbyy(irg,itg,ipg) = 2.0d0*dbvydy
        rlbyz(irg,itg,ipg) = dbvzdy + dbvydz
        rlbzz(irg,itg,ipg) = 2.0d0*dbvzdz
!
      end do
    end do
  end do
  if (info /= 1) write(6,*) ' ### Warning K_ij = 0 *** '
!
  deallocate(dbvxu)
  deallocate(dbvyu)
  deallocate(dbvzu)
  deallocate(dfdx)
  deallocate(dfdy)
  deallocate(dfdz)
!
end subroutine excurve_WL
