subroutine source_bioa_sph(sou_bioa)
  use grid_parameter, only : long, nrg, ntg, npg
  use phys_constant, only : pi
  use coordinate_grav_r, only : hrg
  use def_metric, only : psi, alph, bvxu, bvyu, bvzu
  use trigonometry_grav_theta, only : hsinthg, hcosthg
  use trigonometry_grav_phi, only : hsinphig, hcosphig
  use interface_interpo_linear_type0
  implicit none
  real(long), pointer :: sou_bioa(:,:,:,:)
  real(long) :: psigc,alpgc, bvxgc, bvygc, bvzgc, xx,yy,zz,rr,rxy
  real(long) :: br,bth,bph, psi4, rr2, rxy2
  integer :: ipg, irg, itg
!
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        xx   = hrg(irg)*hsinthg(itg)*hcosphig(ipg)
        yy   = hrg(irg)*hsinthg(itg)*hsinphig(ipg)
        zz   = hrg(irg)*hcosthg(itg)
        rr   = hrg(irg)
        rr2  = rr*rr
        rxy  = dsqrt(xx*xx+yy*yy)
        rxy2 = xx*xx+yy*yy

        call interpo_linear_type0(psigc, psi,irg,itg,ipg)
        call interpo_linear_type0(alpgc,alph,irg,itg,ipg)
        call interpo_linear_type0(bvxgc,bvxu,irg,itg,ipg)
        call interpo_linear_type0(bvygc,bvyu,irg,itg,ipg)
        call interpo_linear_type0(bvzgc,bvzu,irg,itg,ipg)
        psi4 = psigc*psigc*psigc*psigc
!
        br  = (bvxgc*xx + bvygc*yy + bvzgc*zz)/rr
        bth = (bvxgc*xx*zz/rxy + bvygc*yy*zz/rxy - bvzgc*rxy)/rr2 
        bph = (-bvxgc*yy + bvygc*xx)/rxy2

        sou_bioa(irg,itg,ipg,1) = -2.0d0*psi4*br/alpgc 
        sou_bioa(irg,itg,ipg,2) = -2.0d0*psi4*bth/alpgc
        sou_bioa(irg,itg,ipg,3) = -2.0d0*psi4*bph/alpgc
      end do
    end do
  end do
!
end subroutine source_bioa_sph
