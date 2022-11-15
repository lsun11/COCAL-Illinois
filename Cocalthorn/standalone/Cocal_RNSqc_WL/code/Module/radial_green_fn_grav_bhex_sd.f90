! --- Computation of the radial green's function hgfn. ---
! ________________________________________________________
module radial_green_fn_grav_bhex_sd
  use phys_constant, only : long
  use grid_parameter, only : nrg, nlg, rgin, rgout
  use coordinate_grav_r, only : rg, rginv, hrg, hrginv
  use make_array_3d
  implicit none
  real(long), pointer  ::  hgfn_sd(:,:,:), gfnsf_sd(:,:,:)
contains
subroutine allocate_hgfn_bhex_sd
  implicit none
  call alloc_array3d(hgfn_sd,1,nrg,0,nlg,0,nrg)
  call alloc_array3d(gfnsf_sd,0,nlg,0,nrg,1,4)
end subroutine allocate_hgfn_bhex_sd
! ----------------------------------------------
! Subroutine
! hgfn(irr,il,ir) --> for r' < r
!                     for r  < r'
! _sd --> NS (no boundary), Asymptotics Dirichlet
! N.B.  rg(0) = rgin should be zero
! ----------------------------------------------
subroutine calc_hgfn_bhex_sd
  implicit none
  integer     ::  ir, irr, il
  real(long)  ::  radext, radextinv, radsq
  real(long)  ::  fac1, fac2, fac3, dil1, dil2
!
  do ir  = 0, nrg
    do il = 0, nlg
      do irr = 1, nrg
        hgfn_sd(irr,il,ir) = 0.d0
      end do
    end do
  end do
! 
! --  3. Drichlet at rgin and rgout
!
  radext    = rgout
  radextinv = 1.0d0/rgout
!
  do ir  = 0, nrg
    do irr = 1, nrg
      if (hrg(irr).lt.rg(ir)) then
        do il  = 0, nlg
          fac1 = 1.0d0
          fac2 = radextinv**(il+1)
          hgfn_sd(irr,il,ir) = fac1*fac2 &
          &   *hrg(irr)**il &
          &   *((radext*rginv(ir))**(il+1) - (rg(ir)*radextinv)**il)
        end do
      end if
      if (hrg(irr).ge.rg(ir)) then
        do il  = 0, nlg
          fac1 = 1.0d0
          fac2 = radextinv**(il+1)
          hgfn_sd(irr,il,ir) = fac1*fac2 &
          &   *rg(ir)**il &
          &   *((radext*hrginv(irr))**(il+1) - (hrg(irr)*radextinv)**il)
        end do
      end if
    end do
  end do
!
! --- Green's functions for the surface term (TIMES R^2)
!
!     ia --> 1 for d phi/dr condition at rgin (inward normal). 
!            2 for   phi(r) condition at rgin (inward normal).
!            3 for d phi/dr at rgout (outward normal). 
!            4 for   phi(r) at rgout (outward normal).
!     gfnsf_(l,r,ia)
!
  radext    = rgout
  radextinv = 1.0d0/rgout
  radsq     = radext**2
!
  do ir = 0, nrg
    il = 0
!
    fac1 = 1.0d0
    fac2 = 0.0d0
    fac3 = dble(2*il+1)*radextinv**(il+2)
    gfnsf_sd(il,ir,1) = - 0.0d0
    gfnsf_sd(il,ir,2) = - 0.0d0
    gfnsf_sd(il,ir,3) = 0.0d0
    gfnsf_sd(il,ir,4) = - fac1*fac3*radsq
!
    do il = 1, nlg
!
      fac1 = 1.0d0
      fac2 = 0.0d0
      fac3 = dble(2*il+1)*radextinv**(il+2)
      gfnsf_sd(il,ir,1) = - 0.0d0
      gfnsf_sd(il,ir,2) = - 0.0d0
      gfnsf_sd(il,ir,3) = 0.0d0
      gfnsf_sd(il,ir,4) = - fac1*fac3 &
      &                 *rg(ir)**il*radsq
!
    end do
  end do
!
end subroutine calc_hgfn_bhex_sd
end module radial_green_fn_grav_bhex_sd
