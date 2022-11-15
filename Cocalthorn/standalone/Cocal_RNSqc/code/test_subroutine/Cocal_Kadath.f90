include '../Include_file/include_modulefiles_BBH_CF_3mpt.f90'
include '../Include_file/include_interface_modulefiles_BBH_CF_3mpt.f90'
include '../Include_file/include_subroutines_BBH_CF_3mpt.f90'
include '../Include_file/include_functions.f90'
!--
PROGRAM Cocal_Kadath
  use phys_constant, only : long, nnrg, nntg, nnpg
  implicit none
  integer :: ir, it, ip, nrtmp, nttmp, nptmp, nrg, ntg, npg
  real(long) :: psi(0:nnrg,0:nntg,0:nnpg),  psi_K(0:nnrg),  psi_C
  real(long) :: alph(0:nnrg,0:nntg,0:nnpg), alph_K(0:nnrg), alph_C
  real(long) :: bvxd(0:nnrg,0:nntg,0:nnpg), bvxd_K(0:nnrg), bvxd_C
  real(long) :: bvyd(0:nnrg,0:nntg,0:nnpg), bvyd_K(0:nnrg), bvyd_C
  real(long) :: bvzd(0:nnrg,0:nntg,0:nnpg), bvzd_K(0:nnrg), bvzd_C
  real(long) :: rg(0:nnrg), thg(0:nntg), phig(0:nnpg)
  real(long) :: rg_K(0:nnrg), rv
!
! --- Metric potentials.
  open(13,file='input_Cocal.dat',status='old')
  read(13,'(5i5)') nrtmp, nttmp, nptmp
  do ip = 0, nptmp
    do it = 0, nttmp
      do ir = 0, nrtmp
        read(13,'(1p,6e20.12)')  psi(ir,it,ip), &
    &                           alph(ir,it,ip), &
    &                           bvxd(ir,it,ip), &
    &                           bvyd(ir,it,ip), &
    &                           bvzd(ir,it,ip)
      end do
    end do
  end do
! --- Coordinate grids.
  open(14,file='input_Cocal_grid.dat',status='old')
  read(14,'(5i5)') nrtmp, nttmp, nptmp
  do ir = 0, nrtmp
    read(14,'(1p,6e20.12)') rg(ir)
  end do
  do it = 0, nttmp
    read(14,'(1p,6e20.12)') thg(it)
  end do
  do ip = 0, nptmp
    read(14,'(1p,6e20.12)') phig(ip)
  end do
  close(14)
! --- Kadath input
  nrg = nrtmp; ntg = nttmp; npg = nptmp
  open(13,file='input_Kadath.dat',status='old')
  it = ntg/2; ip = 0
  do ir = 0, nrtmp
    read(13,*)  rg_K(ir), psi_K(ir), &
    &                    alph_K(ir), &
    &                    bvxd_K(ir), &
    &                    bvyd_K(ir), &
    &                    bvzd_K(ir)
  end do
  close(13)
!
! --- Interpolation
!
  open(13,file='output_Cocal_Kadath.dat',status='unknown')
  rg(0:nrg) = rg(0:nrg)/rg(0)
  it = ntg/2; ip = 0
  do ir = 0, nrg
    rv = rg_K(ir)
    if (rv.gt.rg(nrg)) exit
    call interpo_radial(psi,psi_C,rv,rg,nrg,it,ip)
    call interpo_radial(alph,alph_C,rv,rg,nrg,it,ip)
    call interpo_radial(bvxd,bvxd_C,rv,rg,nrg,it,ip)
    call interpo_radial(bvyd,bvyd_C,rv,rg,nrg,it,ip)
    call interpo_radial(bvzd,bvzd_C,rv,rg,nrg,it,ip)
    write(13,'(1p,11e20.12)') rg_K(ir), psi_C,  psi_K(ir), &
    &                                  alph_C, alph_K(ir), &
    &                                  bvxd_C, bvxd_K(ir), &
    &                                  bvyd_C, bvyd_K(ir), &
    &                                  bvzd_C, bvzd_K(ir)
  end do
  close(13) 
!
end program Cocal_Kadath
!
subroutine interpo_radial(grv,val,rv,rg,nrg,it,ip)
  use phys_constant, only : long, nnrg, nntg, nnpg
  implicit none
  real(long), external :: lagint_4th
  real(long) :: rg(0:nnrg)
  real(long) :: grv(0:nnrg,0:nntg,0:nnpg)
  real(long), intent(out) :: val
  real(long), intent(in)  :: rv
  integer, intent(in)     :: it, ip
  real(long) :: x(4), f(4)
  integer :: nrg, irg, ir0
!
  do irg = 0, nrg-1
    if (rv.le.rg(irg)) then 
      ir0 = min0(max0(0,irg-2),nrg-3)
      exit
    end if
  end do
  x(1:4) = rg(ir0:ir0+3)
  f(1:4) = grv(ir0:ir0+3,it,ip)
  val = lagint_4th(x,f,rv)
!
end subroutine interpo_radial
