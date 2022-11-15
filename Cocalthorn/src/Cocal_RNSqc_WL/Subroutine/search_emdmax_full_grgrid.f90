subroutine search_emdmax_full_grgrid(iremax,ithemax,iphemax)
  use phys_constant,  only  : long
  use grid_parameter, only  : nrg,ntg,npg
  use def_matter, only  :   emdg
  use def_quantities, only : rho_max, pre_max, epsi_max, q_max
  use coordinate_grav_r, only : rg  
  implicit none
  integer, intent(out) :: iremax, ithemax,iphemax
  real(long) :: emdmax, hmax
  integer    :: irg,itg,ipg
!
  emdmax = -1.0d0
  iremax = 0
  ithemax= 0
  iphemax= 0

  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        if (emdg(irg,itg,ipg).gt.emdmax) then
          emdmax = emdg(irg,itg,ipg)
          iremax = irg
          ithemax= itg
          iphemax= ipg
        end if
      end do
    end do
  end do

  q_max = emdmax
  call peos_q2hprho(q_max, hmax, pre_max, rho_max, epsi_max)
!
end subroutine search_emdmax_full_grgrid
