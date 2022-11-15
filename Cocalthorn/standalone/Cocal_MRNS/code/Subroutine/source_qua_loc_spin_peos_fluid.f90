subroutine source_qua_loc_spin_peos_fluid(souf,irs)
  use phys_constant, only :  long, nmpt
  use grid_parameter, only : nrf, ntf, npf
  use coordinate_grav_r, only : hrg
  use def_metric, only : psi, tfkij
  use def_metric_on_SFC_CF, only : psif
  use def_excurve_on_SFC_CF
  use def_matter, only : rs
  use def_vector_phi, only : hvec_phif
  use def_vector_x, only : hvec_xf
  use interface_interpo_gr2fl
  use interface_interpo_linear_surface_type0
  use interface_calc_surface_normal_tangent_midpoint
  use make_array_1d
  use make_array_3d
  use make_array_4d
  implicit none
  real(long), pointer :: souf(:,:), nv(:), tv(:)
  real(long), pointer :: tfkij_surf(:,:,:,:), tempf(:,:,:)
  real(long) :: Aij_surf, psif6, psifc
  integer    :: irf, itf, ipf, irs, ia, ib, ii
!
  call alloc_array4d(tfkij_surf,0,ntf,0,npf,1,3,1,3)
  call alloc_array3d(tempf,0,nrf,0,ntf,0,npf)
  call alloc_array1d(nv, 1, 3)
  call alloc_array1d(tv, 1, 3)
  call interpo_gr2fl(psi, psif)

  irf = irs
  do ib = 1, 3
    do ia = 1, 3
      tempf(0:nrf,0:ntf,0:npf) = tfkij_grid_fluid(0:nrf,0:ntf,0:npf,ia,ib) 
      do itf = 1, ntf
        do ipf = 1, npf
          call interpo_linear_surface_type0(Aij_surf,tempf,irf,itf,ipf)
          tfkij_surf(itf,ipf,ia,ib) = Aij_surf
        end do
      end do
    end do
  end do

  call calc_vector_x_matter(1)
  call calc_vector_phi_matter(1)

! irf is a gridpoint, and itf,ipf are midpoints
  souf = 0.0d0
  do ipf = 1, npf
    do itf = 1, ntf
      nv = 0.0d0;   tv=0.0d0
      call interpo_linear_surface_type0(psifc,psif,irf,itf,ipf)
      psif6 = psifc**6

      if (irf==nrf) then
        call calc_surface_normal_tangent_midpoint(rs,nv,tv,itf,ipf)
!        tv(1:3) = hvec_phif(irf,itf,ipf,1:3)
      else
        nv(1:3) = hvec_xf(irf,itf,ipf,1:3)/hrg(irf)
        tv(1:3) = hvec_phif(irf,itf,ipf,1:3)
      end if

      do ib = 1, 3
        do ia = 1, 3
          Aij_surf = tfkij_surf(itf,ipf,ia,ib)
          souf(itf,ipf) = souf(itf,ipf) + Aij_surf*tv(ia)*nv(ib)
        end do
      end do
!      if (irf==nrf.and.itf==ntf.and.ipf==1) then
!        write(6,'(a7,1p,7e14.6)') "sou,psi", souf(itf,ipf), psif6
!        write(6,'(a7,1p,7e14.6)') "tv,nv= ", tv(1), tv(2), tv(3), nv(1), nv(2), nv(3)
!        write(6,'(a7,1p,7e14.6)') "tfkij= ", tfkij_surf(itf,ipf,1,1), tfkij_surf(itf,ipf,1,2), tfkij_surf(itf,ipf,1,3), &
!         &                     tfkij_surf(itf,ipf,2,2), tfkij_surf(itf,ipf,2,3), tfkij_surf(itf,ipf,3,3)
!      end if

      souf(itf,ipf) = souf(itf,ipf)*psif6
    end do
  end do
!
  deallocate(tfkij_surf)
  deallocate(tempf)
  deallocate(nv)
  deallocate(tv)
end subroutine source_qua_loc_spin_peos_fluid
