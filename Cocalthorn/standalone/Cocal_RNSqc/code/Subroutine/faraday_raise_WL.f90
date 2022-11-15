subroutine faraday_raise_WL
  use phys_constant,  only : long
  use make_array_3d
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_faraday_tensor,  only : fxd, fyd, fzd, fxu, fyu, fzu, &
  &                               fijd, fiju, fijdu, fidfiu, fijfij
  use interface_interpo_linear_type0
  use interface_index_vec_down2up_midpoint
  implicit none
  real(8) :: hhxxu, hhxyu, hhxzu, hhyyu, hhyzu, hhzzu
  real(8) :: gmiju(3,3), fijdc(3,3), fijuc(3,3), fid(3), fiu(3)
  integer :: irg, itg, ipg, ia, ib, ic
!
! --- Compute components of faraday tensor
! --- whose values are assigned on the mid points. 
!
  call index_vec_down2up_midpoint(fxu,fyu,fzu,fxd,fyd,fzd)
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        call interpo_linear_type0(hhxxu,hxxu,irg,itg,ipg)
        call interpo_linear_type0(hhxyu,hxyu,irg,itg,ipg)
        call interpo_linear_type0(hhxzu,hxzu,irg,itg,ipg)
        call interpo_linear_type0(hhyyu,hyyu,irg,itg,ipg)
        call interpo_linear_type0(hhyzu,hyzu,irg,itg,ipg)
        call interpo_linear_type0(hhzzu,hzzu,irg,itg,ipg)
        gmiju(1,1) = 1.0d0 + hhxxu
        gmiju(1,2) =         hhxyu
        gmiju(1,3) =         hhxzu
        gmiju(2,1) =         hhxyu
        gmiju(2,2) = 1.0d0 + hhyyu
        gmiju(2,3) =         hhyzu
        gmiju(3,1) =         hhxzu
        gmiju(3,2) =         hhyzu
        gmiju(3,3) = 1.0d0 + hhzzu
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
        fiju(irg,itg,ipg,1) = 0.0d0
        fiju(irg,itg,ipg,2) = 0.0d0
        fiju(irg,itg,ipg,3) = 0.0d0
        do ib = 1, 3
          do ia = 1, 3
            fiju(irg,itg,ipg,1) = fiju(irg,itg,ipg,1) &
            &                   + gmiju(1,ia)*gmiju(2,ib)*fijdc(ia,ib)
            fiju(irg,itg,ipg,2) = fiju(irg,itg,ipg,2) & 
            &                   + gmiju(1,ia)*gmiju(3,ib)*fijdc(ia,ib)
            fiju(irg,itg,ipg,3) = fiju(irg,itg,ipg,3) & 
            &                   + gmiju(2,ia)*gmiju(3,ib)*fijdc(ia,ib)
            fijdu(irg,itg,ipg,ia,ib) = 0.0d0
            do ic = 1, 3
              fijdu(irg,itg,ipg,ia,ib) = fijdu(irg,itg,ipg,ia,ib) &
              &                        + gmiju(ib,ic)*fijdc(ia,ic)
            end do
          end do
        end do        
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
end subroutine faraday_raise_WL
