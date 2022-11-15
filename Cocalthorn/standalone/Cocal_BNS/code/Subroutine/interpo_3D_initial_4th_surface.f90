subroutine interpo_3D_initial_4th_surface(fnc)
  use phys_constant,  only : long
  use grid_parameter, only : ntg, npg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi,   only : phig
  use grid_temporary, only : ntftmp, npftmp, thgtmp, phigtmp
  use make_array_2d
  implicit none
  real(long), pointer  :: fnc(:,:)
  real(long), pointer  :: fnctmp(:,:)
  real(long) ::  thc, phic
  real(long) ::  th4(4), phi4(4), ft4(4), fp4(4)
  integer :: itg, ipg, it, ip
  integer :: it0, ip0, itg0, ipg0, ii, jj, kk
  real(long), external :: lagint_4th
!
! --- Interpolation to a grid point (thc,phic)
! --- using 4th order Lagrange formula.
!
  call alloc_array2d(fnctmp,0,ntftmp,0,npftmp)
!
  fnctmp(0:ntftmp,0:npftmp) = 0.0d0
  fnctmp(0:ntftmp,0:npftmp) = fnc(0:ntftmp,0:npftmp)
!
  do ipg = 0, npg
    phic = phig(ipg)
    do ip = 0, npftmp
      if (phic.le.phigtmp(ip)) then
        ip0 = min0(max0(0,ip-2),npftmp-3)
        exit
      end if
    end do
    do itg = 0, ntg
      thc = thg(itg)
      do it = 0, ntftmp
        if (thc.le.thgtmp(it)) then 
          it0 = min0(max0(0,it-2),ntftmp-3)
          exit
        end if
      end do 
!
      do ii = 1, 4
        itg0 = it0 + ii - 1
        ipg0 = ip0 + ii - 1
        th4(ii) = thgtmp(itg0)
        phi4(ii) = phigtmp(ipg0)
      end do
!
      do kk = 1, 4
        ipg0 = ip0 + kk - 1
        do jj = 1, 4
          itg0 = it0 + jj - 1
          ft4(jj) = fnctmp(itg0,ipg0)
        end do
        fp4(kk) = lagint_4th(th4,ft4,thc)
      end do
      fnc(itg,ipg) = lagint_4th(phi4,fp4,phic)
!
    end do
  end do
!
  deallocate(fnctmp)
!
end subroutine interpo_3D_initial_4th_surface
