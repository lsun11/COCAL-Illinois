subroutine source_quad_pole_peos(souf5d, iquad)
  use phys_constant, only  :   long
  use grid_parameter, only  :   nrf, ntf, npf
  use def_matter, only  :   emd, rs, utf
  use def_matter_parameter, only  :   ber
  use def_metric, only  :   psi, alph
  use coordinate_grav_r, only : rg
  use trigonometry_grav_phi, only : sinphig, cosphig
  use trigonometry_grav_theta, only : sinthg, costhg
  use def_vector_x,   only : vec_xf
  use def_vector_phi, only : vec_phif
  use make_array_3d
  use interface_interpo_gr2fl
  implicit none
  real(long), pointer ::  souf5d(:,:,:,:,:)
  real(long), pointer ::  alphf(:,:,:), psif(:,:,:)
  real(long)  ::   emdw, alphw, psiw, rhow, hhw, utw
  real(long)  ::   small = 1.0d-15
  real(long)  ::   prew, ene
  real(long)  ::   xx(1:3), rr, ox(1:3), px(1:3)
  integer     ::   ir, it, ip, i, j, iquad
!
  call alloc_array3d(psif, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(alphf, 0, nrf, 0, ntf, 0, npf)
  call interpo_gr2fl(alph, alphf)
  call interpo_gr2fl(psi, psif)
!
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
        emdw = emd(ir,it,ip)
        if (emdw <= small) emdw = small
        psiw = psif(ir,it,ip)
        alphw = alphf(ir,it,ip)
        call peos_q2hprho(emdw, hhw, prew, rhow, ene)
!        utw = hhw/ber
        utw = utf(ir,it,ip)
!
        xx(1) = vec_xf(ir,it,ip,1)
        xx(2) = vec_xf(ir,it,ip,2)
        xx(3) = vec_xf(ir,it,ip,3)
        rr = xx(1)**2+xx(2)**2+xx(3)**2
        ox(1:2) = xx(1:2)
        ox(3) = 0.0d0
        px(1) = vec_phif(ir,it,ip,1)
        px(2) = vec_phif(ir,it,ip,2)
        px(3) = vec_phif(ir,it,ip,3)
!
        do i=1, 3
          do j=1, 3
            if (iquad == 0 .or. iquad == 1) &
            & souf5d(ir,it,ip,i,j) = xx(i)*xx(j)
            if (iquad == 1 .and. i == j) &
            & souf5d(ir,it,ip,i,j) = xx(i)*xx(j) - rr/3.0d0
            if (iquad == 2) &
            & souf5d(ir,it,ip,i,j) = px(i)*xx(j) + xx(i)*px(j)
            if (iquad == 3) &
            & souf5d(ir,it,ip,i,j) = &
            & ox(i)*xx(j) + xx(i)*ox(j) - 2.0d0*px(i)*px(j)
            if (iquad == 4) &
            & souf5d(ir,it,ip,i,j) = &
            & px(i)*xx(j) + 3.0d0*ox(i)*px(j) + 3.0d0*px(i)*ox(j) + xx(i)*px(j)
            souf5d(ir,it,ip,i,j) = souf5d(ir,it,ip,i,j)*rhow*alphw*utw*psiw**6
          end do
        end do
!
      end do
    end do
  end do
!
  deallocate(alphf)
  deallocate(psif)
end subroutine source_quad_pole_peos
