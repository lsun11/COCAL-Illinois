! --- Computation of the radial green's function hgfn. ---
! ________________________________________________________
module radial_green_fn_grav_bhex_di
  use phys_constant, only : long
  use grid_parameter, only : nrg, nlg, rgin, rgout
  use coordinate_grav_r, only : rg, rginv, hrg, hrginv
  use make_array_3d
  implicit none
  real(long), pointer  ::  hgfn_di(:,:,:), gfnsf_di(:,:,:)
contains
subroutine allocate_hgfn_bhex_di
  implicit none
  call alloc_array3d(hgfn_di,1,nrg,0,nlg,0,nrg)
  call alloc_array3d(gfnsf_di,0,nlg,0,nrg,1,4)
end subroutine allocate_hgfn_bhex_di
! ----------------------------------------------
! Subroutine
! hgfn(irr,il,ir) --> for r' < r
!                     for r  < r'
! _nb --> no boundary
! _di --> BH dirichlet, Asymptotics no boundary
! _dd --> BH dirichlet, Asymptotics Dirichlet
! _nd --> BH Neumann,   Asymptotics Dirichlet
! N.B.  rg(0) = rgin should not equal to zero
! ----------------------------------------------
subroutine calc_hgfn_bhex_di
  implicit none
  integer     ::  ir, irr, il
  real(long)  ::  bhrad,  bhradinv,  bhsq
  real(long)  ::  radext, radextinv, radsq
  real(long)  ::  fac1, fac2, fac3, dil1, dil2
!
  do ir  = 0, nrg
    do il = 0, nlg
      do irr = 1, nrg
        hgfn_di(irr,il,ir) = 0.d0
      end do
    end do
  end do
! 
! --  2. Dirichlet at rgin
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
          hgfn_di(irr,il,ir) =(hrg(irr)**il/rg(ir)**(il+1) &
          &                  - bhrad**(2*il+1)/(hrg(irr)*rg(ir))**(il+1))
        end do
      end if
      if (hrg(irr).ge.rg(ir)) then
        do il  = 0, nlg
          hgfn_di(irr,il,ir) =(rg(ir)**il/hrg(irr)**(il+1) &
          &                  - bhrad**(2*il+1)/(hrg(irr)*rg(ir))**(il+1))
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
    gfnsf_di(il,ir,1) = - 0.0d0
    gfnsf_di(il,ir,2) = - bhrad*rginv(ir)
    gfnsf_di(il,ir,3) = radext*(1.0d0-bhrad*rginv(ir))
    gfnsf_di(il,ir,4) = - (1.0d0-bhrad*rginv(ir))
!
    do il = 1, nlg
!
      gfnsf_di(il,ir,1) = - 0.0d0
      gfnsf_di(il,ir,2) = &
      &    - (2.0d0*dble(il)+1.0d0)*(bhrad*rginv(ir))**(il+1)
      gfnsf_di(il,ir,3) = radext*(bhrad/radext)**il &
      &    *((rg(ir)*bhradinv)**il-(bhrad*rginv(ir))**(il+1))
      gfnsf_di(il,ir,4) = - dble(il+1)*(bhrad/radext)**il &
      &    *((rg(ir)*bhradinv)**il-(bhrad*rginv(ir))**(il+1))
!
    end do
  end do
!
end subroutine calc_hgfn_bhex_di
end module radial_green_fn_grav_bhex_di
