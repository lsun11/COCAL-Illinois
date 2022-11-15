subroutine test_source_mpt(impt)
  use phys_constant, only     : long, pi
  use coordinate_grav_r, only : rg
  use grid_parameter, only    : nrf, ntf, npf
  use def_matter, only        : emd, rs
  implicit none
  integer     ::   irf, itf, ipf, impt
  real(long)  ::   zfac, small = 1.0d-15
!
!
  if(impt.eq.1) then
    do ipf = 0, npf
      do itf = 0, ntf
        rs(itf,ipf) = 1.0d0
          emd(0,itf,ipf) = 4.0d0*pi
!!          emd(0,itf,ipf) = 4.0d0*pi
        do irf = 1, nrf
          emd(irf,itf,ipf) = 4.0d0*pi
!!          emd(irf,itf,ipf) = 4.0d0*pi*sin(pi*rg(irf))/(pi*rg(irf))
        end do
      end do
    end do
  end if
!
  if(impt.eq.2) then
    do ipf = 0, npf
      do itf = 0, ntf
        rs(itf,ipf) = 1.0d0
        !rs(itf,ipf) = 1.0d0
          emd(0,itf,ipf) = 8.0d0*pi
          !emd(0,itf,ipf) = 8.0d0*pi
!!          emd(0,itf,ipf) = 4.0d0*pi
        do irf = 1, nrf
          emd(irf,itf,ipf) = 8.0d0*pi
          !emd(irf,itf,ipf) = 8.0d0*pi
!!          emd(irf,itf,ipf) = 4.0d0*pi*sin(pi*rg(irf))/(pi*rg(irf))
        end do
      end do
    end do
  end if
end subroutine test_source_mpt
