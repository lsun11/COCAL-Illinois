subroutine sourceterm_helmholtz_solver_test(soug,char_sou)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrg, ntg, npg
  use coordinate_grav_r, only  : rg
  use def_metric, only  :   psi
  use def_matter, only  :   emdg, rs
  use def_matter_parameter, only  :  radi, pinx, ber, ome
  use def_vector_phi, only : vec_phig, hvec_phig
  use interface_interpo_linear_type0
  use interface_grgrad_midpoint
  use interface_grgrad_gridpoint
  use interface_interpo_linear_type0
  use make_array_3d
  implicit none
  real(long), pointer :: soug(:,:,:)
  real(long), pointer :: fncdp(:,:,:)
  real(long), pointer :: dfdx(:,:,:), dfdy(:,:,:), dfdz(:,:,:)
  integer      :: irg,itg,ipg
  real(long)   :: emdw,  dxpsi, dypsi, dzpsi
  real(long)   :: d2psi, vphix, vphiy, vphiz
  character(2) :: char_sou
!
  if (char_sou.eq.'he') then
    do ipg = 1, npg
      do itg = 1, ntg
        do irg = 1, nrg
          call interpo_linear_type0(emdw,emdg,irg,itg,ipg)
          soug(irg,itg,ipg) = emdw
        end do
      end do
    end do
  else if (char_sou.eq.'po') then
!
    call alloc_array3d(dfdx,0,nrg,0,ntg,0,npg)
    call alloc_array3d(dfdy,0,nrg,0,ntg,0,npg)
    call alloc_array3d(dfdz,0,nrg,0,ntg,0,npg)
    call alloc_array3d(fncdp,0,nrg,0,ntg,0,npg)
    call grgrad_gridpoint(psi,dfdx,dfdy,dfdz)
    do ipg = 0, npg
      do itg = 0, ntg
        do irg = 0, nrg
          dxpsi = dfdx(irg,itg,ipg)
          dypsi = dfdy(irg,itg,ipg)
          dzpsi = dfdz(irg,itg,ipg)
          vphix = vec_phig(irg,itg,ipg,1)
          vphiy = vec_phig(irg,itg,ipg,2)
          vphiz = vec_phig(irg,itg,ipg,3)
          fncdp(irg,itg,ipg) = vphix*dxpsi + vphiy*dypsi + vphiz*dzpsi
        end do
      end do
    end do
    call grgrad_midpoint(fncdp,dfdx,dfdy,dfdz)
!
    do ipg = 1, npg
      do itg = 1, ntg
        do irg = 1, nrg
          dxpsi = dfdx(irg,itg,ipg)
          dypsi = dfdy(irg,itg,ipg)
          dzpsi = dfdz(irg,itg,ipg)
          vphix = hvec_phig(irg,itg,ipg,1)
          vphiy = hvec_phig(irg,itg,ipg,2)
          vphiz = hvec_phig(irg,itg,ipg,3)
          d2psi = vphix*dxpsi + vphiy*dypsi + vphiz*dzpsi
          call interpo_linear_type0(emdw,emdg,irg,itg,ipg)
          soug(irg,itg,ipg) = emdw + ome**2*d2psi
        end do
      end do
    end do
!
    deallocate(dfdx,dfdy,dfdz,fncdp)
  end if
!
!
end subroutine sourceterm_helmholtz_solver_test
