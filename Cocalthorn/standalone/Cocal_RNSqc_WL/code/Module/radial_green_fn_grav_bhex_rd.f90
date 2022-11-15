! --- Computation of the radial green's function hgfn. ---
! ________________________________________________________
module radial_green_fn_grav_bhex_rd
  use phys_constant, only : long
  use grid_parameter, only : nrg, nlg, rgin, rgout
  use coordinate_grav_r, only : rg, rginv, hrg, hrginv
  use make_array_3d
  implicit none
  real(long), pointer  ::  hgfn_rd(:,:,:), gfnsf_rd(:,:,:)
contains
subroutine allocate_hgfn_bhex_rd
  implicit none
  call alloc_array3d(hgfn_rd,1,nrg,0,nlg,0,nrg)
  call alloc_array3d(gfnsf_rd,0,nlg,0,nrg,1,4)
end subroutine allocate_hgfn_bhex_rd
! ----------------------------------------------
! Subroutine
! hgfn(irr,il,ir) --> for r' < r
!                     for r  < r'
! _nb --> no boundary
! _di --> BH dirichlet, Asymptotics no boundary
! _dd --> BH dirichlet, Asymptotics Dirichlet
! _nd --> BH Neumann,   Asymptotics Dirichlet
! _rd --> BH Robin,     Asymptotics Dirichlet
! N.B.  rg(0) = rgin should not equal to zero
! ----------------------------------------------
subroutine calc_hgfn_bhex_rd
  implicit none
  integer     ::  ir, irr, il
  real(long)  ::  bhrad,  bhradinv,  bhsq
  real(long)  ::  radext, radextinv, radsq
  real(long)  ::  fac1, fac2, fac3, dil1, dil2
!
  do ir  = 0, nrg
    do il = 0, nlg
      do irr = 1, nrg
        hgfn_rd(irr,il,ir) = 0.d0
      end do
    end do
  end do
! 
! --  3. Robin at rgin and Dirichlet rgout
!
  bhrad     = rgin
  bhradinv  = 1.0d0/rgin
  radext    = rgout
  radextinv = 1.0d0/rgout
!
  do ir  = 0, nrg
    do irr = 1, nrg
      if (hrg(irr).lt.rg(ir)) then
        do il  = 0, nlg
          fac1 = 1.0d0/(1.0d0 + (bhrad*radextinv)**(2*il+1))
          fac2 = (bhrad*radextinv)**il*radextinv
          hgfn_rd(irr,il,ir) = fac1*fac2 &
          &   *((hrg(irr)*bhradinv)**il + (bhrad*hrginv(irr))**(il+1)) &
          &   *((radext*rginv(ir))**(il+1) - (rg(ir)*radextinv)**il)
        end do
      end if
      if (hrg(irr).ge.rg(ir)) then
        do il  = 0, nlg
          fac1 = 1.0d0/(1.0d0 + (bhrad*radextinv)**(2*il+1))
          fac2 = (bhrad*radextinv)**il*radextinv
          hgfn_rd(irr,il,ir) = fac1*fac2 &
          &   *((rg(ir)*bhradinv)**il + (bhrad*rginv(ir))**(il+1)) &
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
  bhrad     = rgin
  bhradinv  = 1.0d0/rgin
  radext    = rgout
  radextinv = 1.0d0/rgout
  bhsq      = bhrad**2
  radsq     = radext**2
!
  do ir = 0, nrg
    il = 0
!
    fac1 = 1.0d0/(1.0d0 + (bhrad*radextinv)**(2*il+1))
    fac2 = 2.0d0*bhrad**il*radextinv**(il+1)
    fac3 = dble(2*il+1)*bhrad**il*radextinv**(il+2)
    gfnsf_rd(il,ir,1) = - fac1*fac2    &
      &    *(radext*rginv(ir)-1.0d0)*bhsq
    gfnsf_rd(il,ir,2) = 0.0d0
    gfnsf_rd(il,ir,3) = 0.0d0
    gfnsf_rd(il,ir,4) = - fac1*fac3*(1.0d0 + bhrad*rginv(ir))*radsq
!
    do il = 1, nlg
!
      fac1 = 1.0d0/(1.0d0 + (bhrad*radextinv)**(2*il+1))
      fac2 = 2.0d0*bhrad**il*radextinv**(il+1)
      fac3 = dble(2*il+1)*bhrad**il*radextinv**(il+2)
      gfnsf_rd(il,ir,1) = - fac1*fac2    &
      &    *((radext*rginv(ir))**(il+1)-(rg(ir)*radextinv)**il)*bhsq
      gfnsf_rd(il,ir,2) = 0.0d0
      gfnsf_rd(il,ir,3) = 0.0d0
      gfnsf_rd(il,ir,4) = - fac1*fac3 &
      &    *((rg(ir)*bhradinv)**il+(bhrad*rginv(ir))**(il+1))*radsq
!
    end do
  end do
!
end subroutine calc_hgfn_bhex_rd
end module radial_green_fn_grav_bhex_rd
