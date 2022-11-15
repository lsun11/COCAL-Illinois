subroutine faraday_raise_CF
  use phys_constant,  only : long
  use make_array_3d
  use grid_parameter, only : nrg, ntg, npg
  use def_faraday_tensor,  only : fxd, fyd, fzd, fxu, fyu, fzu, &
  &                               fijd, fiju, fijdu, fidfiu fijfij
  use interface_interpo_linear_type0
  implicit none
  real(8) :: fijdc(3,3), fijuc(3,3), fid(3), fiu(3)
  integer :: irg, itg, ipg, ia, ib, ic
!
! --- Compute components of faraday tensor
! --- whose values are assigned on the mid points. 
!
  fxu(1:nrg,1:ntg,1:npg) = fxd(1:nrg,1:ntg,1:npg)
  fyu(1:nrg,1:ntg,1:npg) = fyd(1:nrg,1:ntg,1:npg)
  fzu(1:nrg,1:ntg,1:npg) = fzd(1:nrg,1:ntg,1:npg)
!
  fiju(1:nrg,1:ntg,1:npg,1) = fijd(1:nrg,1:ntg,1:npg,1)
  fiju(1:nrg,1:ntg,1:npg,2) = fijd(1:nrg,1:ntg,1:npg,2)
  fiju(1:nrg,1:ntg,1:npg,3) = fijd(1:nrg,1:ntg,1:npg,3)
!
  fijdu(1:nrg,1:ntg,1:npg,1,1) = 0.0d0
  fijdu(1:nrg,1:ntg,1:npg,1,2) = fijd(1:nrg,1:ntg,1:npg,1)
  fijdu(1:nrg,1:ntg,1:npg,1,3) = fijd(1:nrg,1:ntg,1:npg,2)
  fijdu(1:nrg,1:ntg,1:npg,2,1) = - fijd(1:nrg,1:ntg,1:npg,1)
  fijdu(1:nrg,1:ntg,1:npg,2,2) = 0.0d0
  fijdu(1:nrg,1:ntg,1:npg,2,3) = fijd(1:nrg,1:ntg,1:npg,3)
  fijdu(1:nrg,1:ntg,1:npg,3,1) = - fijd(1:nrg,1:ntg,1:npg,2)
  fijdu(1:nrg,1:ntg,1:npg,3,2) = - fijd(1:nrg,1:ntg,1:npg,3)
  fijdu(1:nrg,1:ntg,1:npg,3,3) = 0.0d0
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        fijdc(1,1) = 0.0d0
        fijdc(1,2) = fijd(irg,itg,ipg,1)
        fijdc(1,3) = fijd(irg,itg,ipg,2)
        fijdc(2,1) = - fijdc(1,2)
        fijdc(2,2) = 0.0d0
        fijdc(2,3) = fijd(irg,itg,ipg,3)
        fijdc(3,1) = - fijdc(1,3)
        fijdc(3,2) = - fijdc(2,3)
        fijdc(3,3) = 0.0d0
!
        fijuc(1,1) = 0.0d0
        fijuc(1,2) = fiju(irg,itg,ipg,1)
        fijuc(1,3) = fiju(irg,itg,ipg,2)
        fijuc(2,1) = - fijuc(1,2)
        fijuc(2,2) = 0.0d0
        fijuc(2,3) = fiju(irg,itg,ipg,3)
        fijuc(3,1) = - fijuc(1,3)
        fijuc(3,2) = - fijuc(2,3)
        fijuc(3,3) = 0.0d0
!
        fid(1) = fxd(irg,itg,ipg) ; fiu(1) = fxu(irg,itg,ipg)
        fid(2) = fyd(irg,itg,ipg) ; fiu(2) = fyu(irg,itg,ipg)
        fid(3) = fzd(irg,itg,ipg) ; fiu(3) = fzu(irg,itg,ipg)
!
        fidfiu(irg,itg,ipg) = fid(1)*fiu(1) + fid(2)*fiu(2) + fid(3)*fiu(3)
        fijfij(irg,itg,ipg) = 0.0d0
        do ib = 1, 3
          do ia = 1, 3
            fijfij(irg,itg,ipg) = fijfij(irg,itg,ipg) &
            &                   + fijdc(ia,ib)*fijuc(ia,ib)
          end do
        end do
!
      end do
    end do
  end do
!
end subroutine faraday_raise_CF
