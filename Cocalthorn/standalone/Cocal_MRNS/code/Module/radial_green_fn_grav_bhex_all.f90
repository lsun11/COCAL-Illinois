! --- Computation of the radial green's function hgfn. ---
! ________________________________________________________
module radial_green_fn_grav_bhex_all
  use phys_constant, only : long
  use grid_parameter, only : nrg, nlg, rgin, rgout
  use coordinate_grav_r, only : rg, rginv, hrg, hrginv
  use make_array_3d
  implicit none
  real(long), pointer  ::  hgfn_nb(:,:,:), gfnsf_nb(:,:,:)
  real(long), pointer  ::  hgfn_di(:,:,:), gfnsf_di(:,:,:)
  real(long), pointer  ::  hgfn_dd(:,:,:), gfnsf_dd(:,:,:)
  real(long), pointer  ::  hgfn_nd(:,:,:), gfnsf_nd(:,:,:)
contains
subroutine allocate_hgfn_bhex_all
  implicit none
!
  call alloc_array3d(hgfn_nb,1,nrg,0,nlg,0,nrg)
  call alloc_array3d(hgfn_di,1,nrg,0,nlg,0,nrg)
  call alloc_array3d(hgfn_dd,1,nrg,0,nlg,0,nrg)
  call alloc_array3d(hgfn_nd,1,nrg,0,nlg,0,nrg)
  call alloc_array3d(gfnsf_nb,0,nlg,0,nrg,1,4)
  call alloc_array3d(gfnsf_di,0,nlg,0,nrg,1,4)
  call alloc_array3d(gfnsf_dd,0,nlg,0,nrg,1,4)
  call alloc_array3d(gfnsf_nd,0,nlg,0,nrg,1,4)
!
end subroutine allocate_hgfn_bhex_all
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
subroutine calc_hgfn_bhex_all
  implicit none
  integer     ::  ir, irr, il
  real(long)  ::  bhrad,  bhradinv,  bhsq
  real(long)  ::  radext, radextinv, radsq
  real(long)  ::  fac1, fac2, fac3, dil1, dil2
!
  do ir  = 0, nrg
    do il = 0, nlg
      do irr = 1, nrg
        hgfn_nb(irr,il,ir) = 0.d0
        hgfn_di(irr,il,ir) = 0.d0
        hgfn_dd(irr,il,ir) = 0.d0
        hgfn_nd(irr,il,ir) = 0.d0
      end do
    end do
  end do
! 
! --  1. No boundary Green's fn.
! --  2. Dirichlet at rgin
! --  3. Drichlet at rgin and rgout
! --  4. Neumann at rgin and Dirichlet at rgout
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
! -- 1.
          hgfn_nb(irr,il,ir) = hrg(irr)**il/rg(ir)**(il+1)
! -- 2.
          hgfn_di(irr,il,ir) =(hrg(irr)**il/rg(ir)**(il+1) &
          &                  - bhrad**(2*il+1)/(hrg(irr)*rg(ir))**(il+1))
! -- 3.
          fac1 = 1.0d0/(1.0d0 - (bhrad*radextinv)**(2*il+1))
          fac2 = (bhrad*radextinv)**il*radextinv
          hgfn_dd(irr,il,ir) = fac1*fac2 &
          &   *((hrg(irr)*bhradinv)**il - (bhrad*hrginv(irr))**(il+1)) &
          &   *((radext*rginv(ir))**(il+1) - (rg(ir)*radextinv)**il)
! -- 4.
          dil1 = dble(il)/dble(il+1)
          fac1 = 1.0d0/(1.0d0 + dil1*(bhrad*radextinv)**(2*il+1))
          fac2 = (bhrad*radextinv)**il*radextinv
          hgfn_nd(irr,il,ir) = fac1*fac2 &
          &   *((hrg(irr)*bhradinv)**il + dil1*(bhrad*hrginv(irr))**(il+1)) &
          &   *((radext*rginv(ir))**(il+1) - (rg(ir)*radextinv)**il)
        end do
      end if
      if (hrg(irr).ge.rg(ir)) then
        do il  = 0, nlg
! -- 1.
          hgfn_nb(irr,il,ir) = rg(ir)**il/hrg(irr)**(il+1)
! -- 2.
          hgfn_di(irr,il,ir) =(rg(ir)**il/hrg(irr)**(il+1) &
          &                  - bhrad**(2*il+1)/(hrg(irr)*rg(ir))**(il+1))
! -- 3.
          fac1 = 1.0d0/(1.0d0 - (bhrad*radextinv)**(2*il+1))
          fac2 = (bhrad*radextinv)**il*radextinv
          hgfn_dd(irr,il,ir) = fac1*fac2 &
          &   *((rg(ir)*bhradinv)**il - (bhrad*rginv(ir))**(il+1)) &
          &   *((radext*hrginv(irr))**(il+1) - (hrg(irr)*radextinv)**il)
! -- 4.
          dil1 = dble(il)/dble(il+1)
          fac1 = 1.0d0/(1.0d0 + dil1*(bhrad*radextinv)**(2*il+1))
          fac2 = (bhrad*radextinv)**il*radextinv
          hgfn_nd(irr,il,ir) = fac1*fac2 &
          &   *((rg(ir)*bhradinv)**il + dil1*(bhrad*rginv(ir))**(il+1)) &
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
! -- 1.
    gfnsf_nb(il,ir,1) = - bhrad*(bhrad*rginv(ir))
    gfnsf_nb(il,ir,2) = 0.0d0
    gfnsf_nb(il,ir,3) = radext
    gfnsf_nb(il,ir,4) = - dble(il+1)
! -- 2.
    gfnsf_di(il,ir,1) = - 0.0d0
    gfnsf_di(il,ir,2) = - bhrad*rginv(ir)
    gfnsf_di(il,ir,3) = radext*(1.0d0-bhrad*rginv(ir))
    gfnsf_di(il,ir,4) = - (1.0d0-bhrad*rginv(ir))
! -- 3.
    fac1 = 1.0d0/(1.0d0 - (bhrad*radextinv)**(2*il+1))
    fac2 = dble(2*il+1)*bhrad**(il-1)*radextinv**(il+1)
    fac3 = dble(2*il+1)*bhrad**il*radextinv**(il+2)
    gfnsf_dd(il,ir,1) = - 0.0d0
    gfnsf_dd(il,ir,2) = - fac1*fac2*(radext*rginv(ir) - 1.0d0)*bhsq
    gfnsf_dd(il,ir,3) = 0.0d0
    gfnsf_dd(il,ir,4) = - fac1*fac3*(1.0d0 - bhrad*rginv(ir))*radsq
! -- 4.
    dil1 = dble(il)/dble(il+1)
    dil2 = dble(2*il+1)/dble(il+1)
    fac1 = 1.0d0/(1.0d0 + dil1*(bhrad*radextinv)**(2*il+1))
    fac2 = dil2*(bhrad*radextinv)**il*radextinv
    fac3 = dble(2*il+1)*bhrad**il*radextinv**(il+2)
    gfnsf_nd(il,ir,1) = - fac1*fac2*(radext*rginv(ir) - 1.0d0)*bhsq
    gfnsf_nd(il,ir,2) = - 0.0d0
    gfnsf_nd(il,ir,3) = 0.0d0
    gfnsf_nd(il,ir,4) = - fac1*fac3 &
    &                 *(1.0d0 + dil1*bhrad*rginv(ir))*radsq
!
    do il = 1, nlg
! -- 1.
      gfnsf_nb(il,ir,1) = - bhrad*(bhrad*rginv(ir))**(il+1)
      gfnsf_nb(il,ir,2) = - dble(il)*(bhrad*rginv(ir))**(il+1)
      gfnsf_nb(il,ir,3) = radext*(rg(ir)*radextinv)**il
      gfnsf_nb(il,ir,4) = - dble(il+1)*(rg(ir)*radextinv)**il
! -- 2.
      gfnsf_di(il,ir,1) = - 0.0d0
      gfnsf_di(il,ir,2) = &
      &    - (2.0d0*dble(il)+1.0d0)*(bhrad*rginv(ir))**(il+1)
      gfnsf_di(il,ir,3) = radext*(bhrad/radext)**il &
      &    *((rg(ir)*bhradinv)**il-(bhrad*rginv(ir))**(il+1))
      gfnsf_di(il,ir,4) = - dble(il+1)*(bhrad/radext)**il &
      &    *((rg(ir)*bhradinv)**il-(bhrad*rginv(ir))**(il+1))
! -- 3.
      fac1 = 1.0d0/(1.0d0 - (bhrad*radextinv)**(2*il+1))
      fac2 = dble(2*il+1)*bhrad**(il-1)*radextinv**(il+1)
      fac3 = dble(2*il+1)*bhrad**il*radextinv**(il+2)
      gfnsf_dd(il,ir,1) = - 0.0d0
      gfnsf_dd(il,ir,2) = - fac1*fac2 &
      &    *((radext*rginv(ir))**(il+1)-(rg(ir)*radextinv)**il)*bhsq
      gfnsf_dd(il,ir,3) = 0.0d0
      gfnsf_dd(il,ir,4) = - fac1*fac3 &
      &    *((rg(ir)*bhradinv)**il-(bhrad*rginv(ir))**(il+1))*radsq
! -- 4.
      dil1 = dble(il)/dble(il+1)
      dil2 = dble(2*il+1)/dble(il+1)
      fac1 = 1.0d0/(1.0d0 + dil1*(bhrad*radextinv)**(2*il+1))
      fac2 = dil2*(bhrad*radextinv)**il*radextinv
      fac3 = dble(2*il+1)*bhrad**il*radextinv**(il+2)
      gfnsf_nd(il,ir,1) = - fac1*fac2 &
      &    *((radext*rginv(ir))**(il+1)-(rg(ir)*radextinv)**il)*bhsq
      gfnsf_nd(il,ir,2) = - 0.0d0
      gfnsf_nd(il,ir,3) = 0.0d0
      gfnsf_nd(il,ir,4) = - fac1*fac3 &
      &    *((rg(ir)*bhradinv)**il+dil1*(bhrad*rginv(ir))**(il+1))*radsq
!
    end do
  end do
!
end subroutine calc_hgfn_bhex_all
end module radial_green_fn_grav_bhex_all
