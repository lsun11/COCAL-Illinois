subroutine source_vep_CF_peos_spin(souv)
  use phys_constant, only  : long
  use grid_parameter, only : nrf, ntf, npf
  use def_metric_on_SFC_CF
  use def_matter
  use def_matter_parameter, only : ome, ber, confpow
  use def_vector_phi, only : hvec_phif, vec_phif
  use def_velocity_rot
  use def_velocity_potential
  use interface_interpo_linear_type0
  use interface_flgrad_midpoint_type0
  use interface_flgrad_4th_gridpoint
  use interface_flgrad_2nd_gridpoint
  use make_array_3d
  implicit none
  real(long), pointer :: souv(:,:,:)
  real(long), pointer :: hut(:,:,:), aloh(:,:,:)
  real(long), pointer :: hutwx(:,:,:), hutwy(:,:,:), hutwz(:,:,:)
  real(long) :: psifc, alpfc, bvxdfc, bvydfc, bvzdfc
  real(long) :: rotshx, rotshy, rotshz
  real(long) :: emdfc, utfc, hhfc, prefc, rhofc, enefc, abin, abct
  real(long) :: hhut, pfinv, pf2inv, pf4
  real(long) :: dxemd, dyemd, dzemd, dxpsi, dypsi, dzpsi
  real(long) :: dxvpd, dyvpd, dzvpd, hwx, hwy, hwz
  real(long) :: dxhut, dyhut, dzhut, wdpsi, wdaloh, divw
  real(long) :: dxaloh, dyaloh, dzaloh
  real(long) :: dxwx, dywx, dzwx, dxwy, dywy, dzwy, dxwz, dywz, dzwz
  real(long) :: dxbvxd, dybvxd, dzbvxd, dxbvyd, dybvyd, dzbvyd, &
  &             dxbvzd, dybvzd, dzbvzd, divbvf
  real(long) :: vphif(3)
  real(long) :: ovdfc(3), ovdfc2
  real(long) :: dxvep, dyvep, dzvep, lam, alpfc2
  real(long) :: dvep2, psifc4, psifcp, uih2,w2,wdvep,wterm,wx,wy,wz
!  real(long) :: gcx, gcy, gcz
!  real(long) :: dabvep(3,3)
  integer :: ir, it, ip, irf,itf,ipf
!
  call alloc_array3d(hut , 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(aloh, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(hutwx, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(hutwy, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(hutwz, 0, nrf, 0, ntf, 0, npf)
!
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
        psifc = psif(ir,it,ip)
        alpfc= alphf(ir,it,ip)
        emdfc = emd(ir,it,ip)
!
        call peos_q2hprho(emdfc, hhfc, prefc, rhofc, enefc)
        vphif(1) = vec_phif(ir,it,ip,1)
        vphif(2) = vec_phif(ir,it,ip,2)
        vphif(3) = vec_phif(ir,it,ip,3)
        wx       = wxspf(ir,it,ip)
        wy       = wyspf(ir,it,ip)
        wz       = wzspf(ir,it,ip)
 
        ovdfc(1) = bvxdf(ir,it,ip) + ome*vphif(1) - wx
        ovdfc(2) = bvydf(ir,it,ip) + ome*vphif(2) - wy
        ovdfc(3) = bvzdf(ir,it,ip) + ome*vphif(3) - wz
        call flgrad_2nd_gridpoint(vep,dxvep,dyvep,dzvep,ir,it,ip)
  
        psifc4   = psifc**4
!        psifcp   = psifc**confpow
        alpfc2   = alpfc**2
        lam      = ber + ovdfc(1)*dxvep + ovdfc(2)*dyvep + ovdfc(3)*dzvep

        dvep2    = (dxvep**2 + dyvep**2 + dzvep**2)/psifc4
        wdvep    = (wx*dxvep + wy*dyvep + wz*dzvep)
        w2       = psifc4*(wx*wx + wy*wy + wz*wz)

!        wterm    = wdvep + w2
!        uih2     = dvep2 + 2.0d0*wdvep + w2

        hut(ir,it,ip) = (lam + 2.0d0*wdvep)/(alpfc2 - w2)
        if ( hut(ir,it,ip)<0.0d0 ) then
          write(6,*)  "hut negative....exiting"
          stop
        end if

        hutwx(ir,it,ip) = wx*hut(ir,it,ip)
        hutwy(ir,it,ip) = wy*hut(ir,it,ip)
        hutwz(ir,it,ip) = wz*hut(ir,it,ip)

        utf(ir,it,ip)  = hut(ir,it,ip)/hhfc

        aloh(ir,it,ip) = dlog(alpfc*rhofc/hhfc)
      end do
    end do
  end do
!
! =================================
!     Equation of continuity.
!     Homogeneous solution.
! =================================
! ----
!     Source term.
! ----
  souv(1:nrf,1:ntf,1:npf) = 0.0d0
!
  do ip = 1, npf
    do it = 1, ntf
      do ir = 1, nrf-1
        call interpo_linear_type0(bvxdfc,bvxdf,ir,it,ip)
        call interpo_linear_type0(bvydfc,bvydf,ir,it,ip)
        call interpo_linear_type0(bvzdfc,bvzdf,ir,it,ip)

        call interpo_linear_type0(hwx,wxspf,ir,it,ip)
        call interpo_linear_type0(hwy,wyspf,ir,it,ip)
        call interpo_linear_type0(hwz,wzspf,ir,it,ip)

        rotshx = bvxdfc + ome*hvec_phif(ir,it,ip,1) 
        rotshy = bvydfc + ome*hvec_phif(ir,it,ip,2) 
        rotshz = bvzdfc + ome*hvec_phif(ir,it,ip,3) 
!
        call interpo_linear_type0(hhut,hut,ir,it,ip)
        call interpo_linear_type0(psifc,psif,ir,it,ip)
        pfinv = 1.0d0/psifc
        pf2inv = pfinv**2
        pf4 = psifc**4
!
        call flgrad_midpoint_type0(psif,dxpsi, dypsi, dzpsi, ir,it,ip)
        call flgrad_midpoint_type0( vep,dxvpd, dyvpd, dzvpd, ir,it,ip)
        call flgrad_midpoint_type0( hut,dxhut, dyhut, dzhut, ir,it,ip)
        call flgrad_midpoint_type0(aloh,dxaloh,dyaloh,dzaloh,ir,it,ip)
!
        wdpsi = hhut*(hwx*dxpsi + hwy*dypsi + hwz*dzpsi)
        wdaloh= hhut*(hwx*dxaloh + hwy*dyaloh + hwz*dzaloh)

        call flgrad_midpoint_type0(hutwx,dxwx,dywx,dzwx,ir,it,ip)
        call flgrad_midpoint_type0(hutwx,dxwy,dywy,dzwy,ir,it,ip)
        call flgrad_midpoint_type0(hutwx,dxwz,dywz,dzwz,ir,it,ip)
        divw = dxwx + dywy + dzwz
!
!        if(ir==10.and.it==10)   write(6,*)  divw, wdaloh,wdpsi, hhut, rotshx, dxvpd, dxaloh

        souv(ir,it,ip) = - 2.0d0*pfinv*(dxvpd*dxpsi + dyvpd*dypsi + dzvpd*dzpsi)      &
          &              + pf4*(rotshx*dxhut + rotshy*dyhut + rotshz*dzhut)           &
          &              + (pf4*hhut*rotshx - dxvpd)*dxaloh                           &
          &              + (pf4*hhut*rotshy - dyvpd)*dyaloh                           &
          &              + (pf4*hhut*rotshz - dzvpd)*dzaloh                           &
          &              - pf4*( 6.0d0*wdpsi/psifc + divw + wdaloh )
      end do
    end do
  end do
!
  do ip = 1, npf
    do it = 1, ntf
      souv(nrf,it,ip) = 2.0d0*souv(nrf-1,it,ip) - souv(nrf-2,it,ip)
    end do
  end do

  deallocate(hut)
  deallocate(aloh)
  deallocate(hutwx)
  deallocate(hutwy)
  deallocate(hutwz)

end subroutine source_vep_CF_peos_spin
