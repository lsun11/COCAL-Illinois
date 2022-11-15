! --- Computation of the radial green's function hgfn. ---
! ________________________________________________________
module radial_green_fn_grav
  use phys_constant, only : long
  use grid_parameter, only : nrg, nlg, rgmid, rgout
  use coordinate_grav_r, only : rg, rginv, hrg, hrginv
  use make_array_3d
  implicit none
  real(long), pointer  ::  hgfn(:,:,:), gfnsf(:,:,:)
contains
subroutine allocate_hgfn
  implicit none
!
  call alloc_array3d(hgfn,1,nrg,0,nlg,0,nrg)
!
end subroutine allocate_hgfn
subroutine allocate_hgfn_bhex
  implicit none
!
  call alloc_array3d(hgfn,1,nrg,0,nlg,0,nrg)
  call alloc_array3d(gfnsf,0,nlg,0,nrg,1,4)
!
end subroutine allocate_hgfn_bhex
! -------------------------------------
!! Subroutine
!    hgfn(irr,ir) --> for r' < r
!                     for r < r'
subroutine calc_hgfn
  implicit none
  integer  ::  ir, irr, nn, status
!
  do ir  = 0, nrg
    do nn = 0, nlg
      do irr = 1, nrg
        hgfn(irr,nn,ir) = 0.e0
      end do
    end do
  end do
! 
  do ir  = 1, nrg
    do irr = 1, nrg
    IF (hrg(irr)<rg(ir)) THEN
      do nn = 0, nlg
!        hgfn(irr,nn,ir) = hrg(irr)**nn/rg(ir)**(nn+1)
        hgfn(irr,nn,ir) = (hrg(irr)/rg(ir))**nn/rg(ir)
      end do
    end IF
    IF (hrg(irr)>=rg(ir)) THEN
      do nn = 0, nlg
!        hgfn(irr,nn,ir) = rg(ir)**nn/hrg(irr)**(nn+1)
        hgfn(irr,nn,ir) = (rg(ir)/hrg(irr))**nn/hrg(irr)
      end do
    end IF
    end do
  end do
! 
  ir  = 0
  nn  = 0
  do irr = 1, nrg
    hgfn(irr,nn,ir) = 1.0e0/hrg(irr)
  end do
end subroutine calc_hgfn
end module radial_green_fn_grav
