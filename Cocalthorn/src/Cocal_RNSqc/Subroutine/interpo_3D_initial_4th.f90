subroutine interpo_3D_initial_4th(fnc,char_co)
  use phys_constant, only : long, pi, nnrg, nntg, nnpg
  use grid_parameter, only : nrf, nrg, ntg, npg
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig
  use grid_temporary, only : nrftmp, nrgtmp, ntgtmp, npgtmp, &
  &                                  rgtmp, thgtmp, phigtmp
  use make_array_3d
  implicit none
  real(long), pointer  :: fnc(:,:,:)
  real(long), pointer  :: fnctmp(:,:,:)
  real(long) ::  rc, thc, phic
  real(long) ::  r4(4), th4(4), phi4(4), fr4(4), ft4(4), fp4(4)
  integer :: irg, itg, ipg, ir, it, ip
  integer :: ir0, it0, ip0, irg0 , itg0 , ipg0, ii, jj, kk
  integer :: nrfnrg, nrfnrgtmp
  real(long), external :: lagint_4th
  character(len=4) :: char_co 
!
! --- Interpolation to a grid point (rc,thc,phic)
! --- using 4th order Lagrange formula.
!
  call alloc_array3d(fnctmp,0,nrgtmp,0,ntgtmp,0,npgtmp)
!
  fnctmp(0:nrgtmp,0:ntgtmp,0:npgtmp) = 0.0d0
  fnctmp(0:nrgtmp,0:ntgtmp,0:npgtmp) = fnc(0:nrgtmp,0:ntgtmp,0:npgtmp)
!
  if (char_co.eq.'sfco') then 
    nrfnrg    = nrf
    nrfnrgtmp = nrftmp
  else if (char_co.eq.'grco') then 
    nrfnrg    = nrg
    nrfnrgtmp = nrgtmp
  end if
!
  do ipg = 0, npg
    phic = phig(ipg)
    do ip = 0, npgtmp
      if (phic.le.phigtmp(ip)) then
        ip0 = min0(max0(0,ip-2),npgtmp-3)
        exit
      end if
    end do
    do itg = 0, ntg
      thc = thg(itg)
      do it = 0, ntgtmp
        if (thc.le.thgtmp(it)) then 
          it0 = min0(max0(0,it-2),ntgtmp-3)
          exit
        end if
      end do 
      do irg = 0, nrfnrg
        rc = rg(irg)
        do ir = 0, nrgtmp
          if (rc.le.rgtmp(ir)) then 
            ir0 = min0(max0(0,ir-2),nrfnrgtmp-3)
            exit
          end if
        end do
!
        do ii = 1, 4
          irg0 = ir0 + ii - 1
          itg0 = it0 + ii - 1
          ipg0 = ip0 + ii - 1
          r4(ii) = rgtmp(irg0)
          th4(ii) = thgtmp(itg0)
          phi4(ii) = phigtmp(ipg0)
        end do
!
        do kk = 1, 4
          ipg0 = ip0 + kk - 1
          do jj = 1, 4
            itg0 = it0 + jj - 1
            do ii = 1, 4
              irg0 = ir0 + ii - 1
              fr4(ii) = fnctmp(irg0,itg0,ipg0)
            end do
            ft4(jj) = lagint_4th(r4,fr4,rc)
          end do
          fp4(kk) = lagint_4th(th4,ft4,thc)
        end do
        fnc(irg,itg,ipg) = lagint_4th(phi4,fp4,phic)
!
      end do
    end do
  end do
!
  deallocate(fnctmp)
!
end subroutine interpo_3D_initial_4th
