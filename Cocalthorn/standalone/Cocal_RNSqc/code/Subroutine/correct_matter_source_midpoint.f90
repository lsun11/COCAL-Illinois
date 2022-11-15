subroutine correct_matter_source_midpoint(sou)
  use phys_constant, only : long, pi
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r, only  :   drg, rg, hrg
  use def_matter, only : rs
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long), pointer :: sou(:,:,:)
  real(long) :: hsurf
  integer    :: irg, itg, ipg
  real(long) :: x(2),f(2), v, sougc
  real(long), external :: lagint_2nd
!
! -- Correct the weight for a matter source in a volume integral
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(hsurf,rs,itg,ipg)
      do irg = 1, nrg
        if (hrg(irg) <= hsurf.and.hsurf <= rg(irg)) then 
!!        sou(irg,itg,ipg) = sou(irg,itg,ipg)*(hsurf-rg(irg-1))/drg(irg)
          x(1) = hrg(irg-1) ; f(1) = sou(irg-1,itg,ipg)
          x(2) = hrg(irg)   ; f(2) = sou(irg,itg,ipg)
          v = 0.5d0*(hsurf + rg(irg-1))
          sougc = lagint_2nd(x,f,v)
          sou(irg,itg,ipg) = sougc*(v/hrg(irg))**2 &
          &                *(hsurf-rg(irg-1))/drg(irg)
          exit
        end if
        if (rg(irg) < hsurf.and.hsurf < hrg(irg+1)) then 
!!!        sou(irg,itg,ipg) = sou(irg,itg,ipg)
!!!                         *(1.d0+(hsurf-rg(irg))/drg(irg+1))
          x(1) = hrg(irg-1) ; f(1) = sou(irg-1,itg,ipg)
          x(2) = hrg(irg)   ; f(2) = sou(irg,itg,ipg)
!bak          v = 0.5d0*(hsurf + rg(irg-1))
          v = 0.5d0*(hsurf + rg(irg))
          sougc = lagint_2nd(x,f,v)
!bak          sou(irg,itg,ipg) = sougc*(v/hrg(irg))**2 &
!bak          &                *(1.d0+(hsurf-rg(irg))/drg(irg+1))
          sou(irg+1,itg,ipg) = sougc*(v/hrg(irg+1))**2 &
          &                  *(hsurf-rg(irg))/drg(irg+1)
          exit
        end if
      end do
    end do
  end do
end subroutine correct_matter_source_midpoint
