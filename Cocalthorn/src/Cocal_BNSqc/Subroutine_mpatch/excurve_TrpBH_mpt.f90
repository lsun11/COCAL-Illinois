subroutine excurve_TrpBH_mpt(impt)
  use phys_constant, only : long, nmpt
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_pBH, only : aij_trpBH, aijaij_trpBH
  use def_binary_parameter, only : dis, sepa
  use def_bh_parameter, only : mom_pBH, spin_pBH, mass_pBH
  use def_vector_x, only : hvec_xg
  implicit none
  integer :: impt, impt_bh1, impt_bh2, info = 0
  integer :: irg, itg, ipg, i, j, k, l, m
  real(long) :: dis_bh1, dis_bh2
  real(long) :: hvx_bh1(3), hvx_bh2(3), hn_bh1(3), hn_bh2(3), hr_bh1, hr_bh2
  real(long) :: p_bh1(3), p_bh2(3), s_bh1(3), s_bh2(3)
  real(long) :: fac_bh1, fac_bh2, fac_p, fac_s, hnk_pk_bh1, hnk_pk_bh2, gamma
  real(long) :: epsi(3,3,3), eps_s_n_bh1(3), eps_s_n_bh2(3)
!
  fac_p = 3.0d0/2.0d0
  fac_s = 6.0d0
  epsi(1:3,1:3,1:3) = 0.0d0
  epsi(1,2,3) =  1.0d0 ; epsi(2,3,1) =  1.0d0 ; epsi(3,1,2) =  1.0d0
  epsi(2,1,3) = -1.0d0 ; epsi(1,3,2) = -1.0d0 ; epsi(3,2,1) = -1.0d0
!
  if (impt.eq.1)    then ; impt_bh1 = 1 ; impt_bh2 = 2 ; end if
  if (impt.eq.2)    then ; impt_bh1 = 2 ; impt_bh2 = 1 ; end if
  if (impt.eq.nmpt) then ; impt_bh1 = 1 ; impt_bh2 = 2 ; end if
!
  call copy_def_binary_parameter_from_mpt(impt_bh1)
  call copy_def_bh_parameter_from_mpt(impt_bh1)
  dis_bh1 = 0.0d0
  if (impt.eq.nmpt) dis_bh1 = - dis
  p_bh1(1) = mom_pBH(1) ; s_bh1(1) = spin_pBH(1)
  p_bh1(2) = mom_pBH(2) ; s_bh1(2) = spin_pBH(2)
  p_bh1(3) = mom_pBH(3) ; s_bh1(3) = spin_pBH(3)
  fac_bh1 = 3.0d0*dsqrt(3.0d0)*mass_pBH**2.0d0/4.0d0
  call copy_def_binary_parameter_from_mpt(impt_bh2)
  call copy_def_bh_parameter_from_mpt(impt_bh2)
  dis_bh2 = sepa
  if (impt.eq.nmpt) dis_bh2 = dis
  p_bh2(1) = - mom_pBH(1) ; s_bh2(1) = - spin_pBH(1)
  p_bh2(2) = - mom_pBH(2) ; s_bh2(2) = - spin_pBH(2)
  p_bh2(3) =   mom_pBH(3) ; s_bh2(3) =   spin_pBH(3)
  fac_bh2 = 3.0d0*dsqrt(3.0d0)*mass_pBH**2.0d0/4.0d0
!
  call copy_def_binary_parameter_from_mpt(impt)  ! back to active patch
  call copy_def_bh_parameter_from_mpt(impt)      ! back to active patch
!
  aijaij_trpBH(0:nrg,0:ntg,0:npg) = 0.0d0
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        hvx_bh1(1) = hvec_xg(irg,itg,ipg,1) - dis_bh1
        hvx_bh1(2) = hvec_xg(irg,itg,ipg,2)
        hvx_bh1(3) = hvec_xg(irg,itg,ipg,3)
        hr_bh1     = dsqrt(hvx_bh1(1)**2 + hvx_bh1(2)**2 + hvx_bh1(3)**2)
        hn_bh1(1)  = hvx_bh1(1)/hr_bh1
        hn_bh1(2)  = hvx_bh1(2)/hr_bh1
        hn_bh1(3)  = hvx_bh1(3)/hr_bh1
        hnk_pk_bh1 = hn_bh1(1)*p_bh1(1)+hn_bh1(2)*p_bh1(2)+hn_bh1(3)*p_bh1(3)
!
        hvx_bh2(1) = hvec_xg(irg,itg,ipg,1) - dis_bh2
        hvx_bh2(2) = hvec_xg(irg,itg,ipg,2)
        hvx_bh2(3) = hvec_xg(irg,itg,ipg,3)
        hr_bh2     = dsqrt(hvx_bh2(1)**2 + hvx_bh2(2)**2 + hvx_bh2(3)**2)
        hn_bh2(1) =  hvx_bh2(1)/hr_bh2
        hn_bh2(2) =  hvx_bh2(2)/hr_bh2
        hn_bh2(3) =  hvx_bh2(3)/hr_bh2
        hnk_pk_bh2 = hn_bh2(1)*p_bh2(1)+hn_bh2(2)*p_bh2(2)+hn_bh2(3)*p_bh2(3)
        do k = 1, 3
          eps_s_n_bh1(k) = 0.0d0
          eps_s_n_bh2(k) = 0.0d0
          do l = 1, 3
            do m = 1, 3
              eps_s_n_bh1(k) = eps_s_n_bh1(k)+epsi(k,l,m)*s_bh1(l)*hn_bh1(m)
              eps_s_n_bh2(k) = eps_s_n_bh2(k)+epsi(k,l,m)*s_bh2(l)*hn_bh2(m)
            end do
          end do
        end do
!
        do i = 1, 3
          do j = 1, 3
            if (i == j) then 
              gamma = 1.0d0
            else 
              gamma = 0.0d0
            end if
!
            aij_trpBH(irg,itg,ipg,i,j) = &
            &                     fac_bh1/hr_bh1**3.0d0 &
            &                   *(gamma - 3.0d0*hn_bh1(i)*hn_bh1(j))&
            &                   + fac_p/hr_bh1**2.0d0 &
            &                   *(p_bh1(i)*hn_bh1(j) + p_bh1(j)*hn_bh1(i) &
            &                   - hnk_pk_bh1*(gamma - hn_bh1(i)*hn_bh1(j))) & 
            &                   + fac_s/hr_bh1**3.0d0 &
            &                   *0.5d0*(hn_bh1(i)*eps_s_n_bh1(j)  &
            &                         + hn_bh1(j)*eps_s_n_bh1(i)) &
!
            &                   + fac_bh2/hr_bh2**3.0d0 &
            &                   *(gamma - 3.0d0*hn_bh2(i)*hn_bh2(j)) &
            &                   + fac_p/hr_bh2**2.0d0 &
            &                   *(p_bh2(i)*hn_bh2(j) + p_bh2(j)*hn_bh2(i) &
            &                   - hnk_pk_bh2*(gamma - hn_bh2(i)*hn_bh2(j))) &
            &                   + fac_s/hr_bh2**3.0d0 &
            &                   *0.5d0*(hn_bh2(i)*eps_s_n_bh2(j) &
            &                         + hn_bh2(j)*eps_s_n_bh2(i)) 
!
            aijaij_trpBH(irg,itg,ipg) = aijaij_trpBH(irg,itg,ipg) &
            &                         + aij_trpBH(irg,itg,ipg,i,j)**2.0d0
          end do
        end do
        if (aijaij_trpBH(irg,itg,ipg) /= 0.0d0) info = 1
      end do
    end do
  end do
!
  if (info /= 1) write(6,*) ' ### Warning Aij_trpBH = 0 *** '
!
end subroutine excurve_TrpBH_mpt
