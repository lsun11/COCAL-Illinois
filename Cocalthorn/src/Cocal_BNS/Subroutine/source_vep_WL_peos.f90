subroutine source_vep_WL_peos(souv)
  use phys_constant, only  : long
  use grid_parameter, only : nrf, ntf, npf
  use def_metric_on_SFC, only : alphf, psif, bvxuf, bvyuf, bvzuf, &
  &                             hxxuf, hxyuf, hxzuf, hyyuf, hyzuf, hzzuf
  use def_gamma_crist
  use def_vector_phi, only : hvec_phif
  use interface_interpo_linear_type0
  use interface_flgrad_midpoint_type0
  use make_array_3d
  implicit none
  real(long), pointer :: souv(:,:,:)
  real(long), pointer :: hut(:,:,:), hutp6(:,:,:), aloh(:,:,:)
  real(long) :: psifc, alphfc, bvxufc, bvyufc, bvzufc
  real(long) :: rotshx, rotshy, rotshz
  real(long) :: hhxxu, hhxyu, hhxzu, hhyxu, hhyyu, hhyzu, &
  &             hhzxu, hhzyu, hhzzu
  real(long) :: gamxxu, gamxyu, gamxzu, gamyxu, gamyyu, gamyzu, &
  &             gamzxu, gamzyu, gamzzu
  real(long) :: emdfc, utfc, hhfc, prefc, rhofc, enefc, abin, abct
  real(long) :: hhut, pfinv, pf2inv, pf4
  real(long) :: dxemd, dyemd, dzemd, dxpsi, dypsi, dzpsi
  real(long) :: dxvpd, dyvpd, dzvpd, dxvpu, dyvpu, dzvpu
  real(long) :: dxhutp6, dyhutp6, dzhutp6
  real(long) :: dxaloh, dyaloh, dzaloh, dxlnarh, dylnarh, dzlnarh
  real(long) :: dxbvxu, dybvxu, dzbvxu, dxbvyu, dybvyu, dzbvyu, &
  &             dxbvzu, dybvzu, dzbvzu, divbvf
  real(long) :: gcx, gcy, gcz
  real(long) :: dabvep(3,3)
  integer :: ir, it, ip
!
  call alloc_array3d(hut, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(hutp6, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(aloh, 0, nrf, 0, ntf, 0, npf)
!
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
        psifc = psif(ir,it,ip)
        alphfc= alphf(ir,it,ip)
        emdfc = emd(ir,it,ip)
        utfc  = utf(ir,it,ip)
!
        call peos_q2hprhoab(emdfc, hhfc, prefc, rhofc, enefc, abin, abct)
!
        hut(ir,it,ip) = hhfc*utfc
        hutp6(ir,it,ip) = hut(ir,it,ip)*psifc**6
        aloh(ir,it,ip) = dlog(alphfc/hhfc)
!!!        pinx_peos(ir,it,ip) = 1.0d0/(abin-1.0d0)
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
  souv(:,:,:) = 0.0d0
!
  do ip = 1, npf
    do it = 1, ntf
      do ir = 1, nrf
        call interpo_linear_type0(hhxxu,hxxuf,ir,it,ip)
        call interpo_linear_type0(hhxyu,hxyuf,ir,it,ip)
        call interpo_linear_type0(hhxzu,hxzuf,ir,it,ip)
        call interpo_linear_type0(hhyyu,hyyuf,ir,it,ip)
        call interpo_linear_type0(hhyzu,hyzuf,ir,it,ip)
        call interpo_linear_type0(hhzzu,hzzuf,ir,it,ip)
        hhyxu = hhxyu
        hhzxu = hhxzu
        hhzyu = hhyzu
        gamxxu = hhxxu + 1.0d0
        gamxyu = hhxyu
        gamxzu = hhxzu
        gamyyu = hhyyu + 1.0d0
        gamyzu = hhyzu
        gamzzu = hhzzu + 1.0d0
        gamyxu = gamxyu
        gamzxu = gamxzu
        gamzyu = gamyzu
        call interpo_linear_type0(bvxufc,bvxuf,ir,it,ip)
        call interpo_linear_type0(bvyufc,bvyuf,ir,it,ip)
        call interpo_linear_type0(bvzufc,bvzuf,ir,it,ip)
        rotshx = bvxufc + ome*hvec_phif(ir,it,ip,1)
        rotshy = bvyufc + ome*hvec_phif(ir,it,ip,2)
        rotshz = bvzufc + ome*hvec_phif(ir,it,ip,3)
!
        call interpo_gr2fl_midpoint_type0(gcx,gmcrix,ir,it,ip)
        call interpo_gr2fl_midpoint_type0(gcy,gmcriy,ir,it,ip)
        call interpo_gr2fl_midpoint_type0(gcz,gmcriz,ir,it,ip)
!
        call interpo_linear_type0(psifc,psif,ir,it,ip)
        pfinv = 1.0d0/psifc
        pf2inv = pfinv**2
        pf4 = psifc**4
!
        call flgrad_midpoint_type0(psif,dxpsi,dypsi,dzpsi,ir,it,ip)
        call flgrad_midpoint_type0(vep,dxvpd,dyvpd,dzvpd,ir,it,ip)
        dxvpu = gamxxu*dxvpd + gamxyu*dyvpd + gamxzu*dzvpd
        dyvpu = gamyxu*dxvpd + gamyyu*dyvpd + gamyzu*dzvpd
        dzvpu = gamzxu*dxvpd + gamzyu*dyvpd + gamzzu*dzvpd
!
        call interpo_linear_type0(hhut,hut,ir,it,ip)
        call flgrad_midpoint_type0(hutp6,dxhutp6,dyhutp6,dzhutp6,ir,it,ip)
!
        call interpo_linear_type0(emdfc,emd,ir,it,ip)
        call peos_q2hprhoab(emdfc, hhfc, prefc, rhofc, enefc, abin, abct)
        pind = 1.0d0/(abin-1.0d0)
!!!        pinx_peos(ir,it,ip) = 1.0d0/(abin-1.0d0)
        call flgrad_midpoint_type0(emd,dxemd,dyemd,dzemd,ir,it,ip)
        call flgrad_midpoint_type0(aloh,dxaloh,dyaloh,dzaloh,ir,it,ip)
        dxlnarh = dxaloh + pind*dxemd/emdfc
        dylnarh = dyaloh + pind*dyemd/emdfc
        dzlnarh = dzaloh + pind*dzemd/emdfc
!
        call flgrad_midpoint_type0(bvxuf,dxbvxu,dybvxu,dzbvxu,ir,it,ip)
        call flgrad_midpoint_type0(bvyuf,dxbvyu,dybvyu,dzbvyu,ir,it,ip)
        call flgrad_midpoint_type0(bvzuf,dxbvzu,dybvzu,dzbvzu,ir,it,ip)
        divbvf = dxbvxu + dybvyu + dzbvzu
!
        call dadbscalar_fluid_type0(vep,dabvep,ir,it,ip)
!
        souv(ir,it,ip) = &
     &  -(hhxxu*dabvep(1,1) + hhxyu*dabvep(1,2) + hhxzu*dabvep(1,3)  &
     &  + hhyxu*dabvep(2,1) + hhyyu*dabvep(2,2) + hhyzu*dabvep(2,3)  &
     &  + hhzxu*dabvep(3,1) + hhzyu*dabvep(3,2) + hhzzu*dabvep(3,3)) &
     &  + gcx*dxvpd + gcy*dyvpd + gcz*dzvpd                          &
     &  - 2.0d0*pfinv*(dxvpu*dxpsi + dyvpu*dypsi + dzvpu*dzpsi)      &
     &  + pf2inv*(rotshx*dxhutp6 + rotshy*dyhutp6 + rotshz*dzhutp6)  &
     &  + pf4*hhut*divbvf                                            &
     &  -(dxvpu - pf4*hhut*rotshx)*dxlnarh                           &
     &  -(dyvpu - pf4*hhut*rotshy)*dylnarh                           &
     &  -(dzvpu - pf4*hhut*rotshz)*dzlnarh
!
      end do
    end do
  end do
!
  deallocate()
!
end subroutine source_vep_WL_peos
