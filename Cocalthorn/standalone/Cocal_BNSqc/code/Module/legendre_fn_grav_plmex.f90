!  Associated Legendre function and factorials
!______________________________________________
module legendre_fn_grav_plmex
  use phys_constant, only : long
  implicit none
  real(long), pointer ::    plm_ex(:,:,:,:,:,:)
contains

subroutine allocate_legendre_plmex
  use phys_constant,  only : long
  use grid_parameter, only : nrg,ntg,npg, nlg     ! nlg : maximum value of multipole
  use make_array_6d
  implicit none
  real(long)  :: leg_mem
!
  leg_mem = 2.0d0*dble(nrg)*dble(ntg)*dble(npg)*dble(nlg)*dble(nlg)*8.0d0
  leg_mem = leg_mem/1024.0d0/1024.0d0/1024.0d0
  write(6,'(a50)') "++++++++++++++++++++++++++++++++++++++++++++++++++"
  write(6,'(a28,e15.3,a3)') "Memory allocation of plm_ex:", leg_mem, " GB"  
  write(6,'(a50)') "++++++++++++++++++++++++++++++++++++++++++++++++++"

  call alloc_array6d(plm_ex, 1,2, 0,nrg, 0,ntg, 0,npg, 0,nlg, 0,nlg)
!
end subroutine allocate_legendre_plmex

SUBROUTINE legendre_plmex(impt)
  use phys_constant,  only : nnlg, long, pi
  use grid_parameter, only : nrg, ntg, npg, nlg     ! nlg : maximum value of multipole
!  use coordinate_grav_theta, only : thg, hthg
  use grid_parameter_binary_excision, only : ex_radius
  use grid_points_binary_excision, only : rb, thb, phib, irg_exin, irg_exout
  use make_array_1d
  implicit none
  integer    :: irg,itg,ipg, mm, nn, il, kk, impt
  real(long) :: fmm, fnn
  real(long) :: q1, q2, cc1, ss1, fkk
  real(long) :: rra1, rra1l
  real(long), pointer  ::  fac(:)
!
! --- computation of associated Legendre polynomials wrt the excised origin.
!
  call alloc_array1d(fac, 0, nlg)
  fac(0) = 1.0d0
  do mm = 1, nlg
    fmm = dble(mm)
    fac(mm) = (2.0d0*fmm-1.0d0) * fac(mm-1)
  end do

  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        if (irg.gt.irg_exin(itg,ipg).and. irg.lt.irg_exout(itg,ipg)) then
          plm_ex(impt,irg,itg,ipg,0:nlg,0:nlg) = 0.0d0
        else
          plm_ex(impt,irg,itg,ipg,0:nlg,0:nlg) = 0.0d0
          cc1 = dcos(thb(irg,itg,ipg))
          ss1 = dsin(thb(irg,itg,ipg))
!
          plm_ex(impt,irg,itg,ipg,0,0) = 1.0d0
          do mm = 1, nlg
            plm_ex(impt,irg,itg,ipg,mm,mm) = fac(mm) * (-ss1)**mm
          end do
!
          do mm = 0, nlg-1
            fmm = dble(mm)
            plm_ex(impt,irg,itg,ipg,mm+1,mm) = (2.d0*fmm + 1.d0)*cc1*plm_ex(impt,irg,itg,ipg,mm,mm)
          end do
!
          do mm = 0, nlg-2
            fmm = dble(mm)
            do kk = 2, nlg-mm
              fkk = dble(kk)
              q1 = ( 2.0d0 * fmm + 2.0d0 * fkk - 1.0d0 ) / fkk
              q2 = ( 2.0d0 * fmm + fkk - 1.0d0 ) / fkk
              plm_ex(impt,irg,itg,ipg,mm+kk,mm) = q1 * cc1 * plm_ex(impt,irg,itg,ipg,mm+kk-1,mm) &
             &                                  - q2       * plm_ex(impt,irg,itg,ipg,mm+kk-2,mm)
            end do
          end do

          rra1 = ex_radius/rb(irg,itg,ipg)
          rra1l = 1.0d0
          do il = 0, nlg
            rra1l = rra1l*rra1
            plm_ex(impt,irg,itg,ipg, il,0:nlg) = plm_ex(impt,irg,itg,ipg, il,0:nlg)*rra1l
          end do
!
!          if(ipg==0 .and. itg==0 .and. irg<=5)   write(6,'(a10,1p,2e20.12)') "plm_ex=", plm_ex(impt,irg,itg,ipg,0,0)

        end if
      end do
    end do
  end do
!
  deallocate(fac)
!
end subroutine legendre_plmex

end module legendre_fn_grav_plmex
