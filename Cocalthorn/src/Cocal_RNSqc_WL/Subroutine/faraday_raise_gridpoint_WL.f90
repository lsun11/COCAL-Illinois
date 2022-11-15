subroutine faraday_raise_gridpoint_WL
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_faraday_tensor,  only : fxd_grid, fyd_grid, fzd_grid, &
  &                               fxu_grid, fyu_grid, fzu_grid, &
  &                               fijd_grid, fiju_grid, fijdu_grid, &
  &                               fidfiu_grid, fijfij_grid
  use interface_index_vec_down2up
  implicit none
  real(8) :: gmiju(3,3), fijdc(3,3), fijuc(3,3), fid(3), fiu(3)
  integer :: irg, itg, ipg, ia, ib, ic
!
! --- Compute components of faraday tensor
! --- whose values are assigned on the mid points. 
!
  call index_vec_down2up(fxu_grid,fyu_grid,fzu_grid,fxd_grid,fyd_grid,fzd_grid)
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        gmiju(1,1) = 1.0d0 + hxxu(irg,itg,ipg)
        gmiju(1,2) =         hxyu(irg,itg,ipg)
        gmiju(1,3) =         hxzu(irg,itg,ipg)
        gmiju(2,1) =         hxyu(irg,itg,ipg)
        gmiju(2,2) = 1.0d0 + hyyu(irg,itg,ipg)
        gmiju(2,3) =         hyzu(irg,itg,ipg)
        gmiju(3,1) =         hxzu(irg,itg,ipg)
        gmiju(3,2) =         hyzu(irg,itg,ipg)
        gmiju(3,3) = 1.0d0 + hzzu(irg,itg,ipg)
!
        fijdc(1,1) = 0.0d0
        fijdc(1,2) = fijd_grid(irg,itg,ipg,1)
        fijdc(1,3) = fijd_grid(irg,itg,ipg,2)
        fijdc(2,1) = - fijdc(1,2)
        fijdc(2,2) = 0.0d0
        fijdc(2,3) = fijd_grid(irg,itg,ipg,3)
        fijdc(3,1) = - fijdc(1,3)
        fijdc(3,2) = - fijdc(2,3)
        fijdc(3,3) = 0.0d0
!
        fiju_grid(irg,itg,ipg,1) = 0.0d0
        fiju_grid(irg,itg,ipg,2) = 0.0d0
        fiju_grid(irg,itg,ipg,3) = 0.0d0
        do ib = 1, 3
          do ia = 1, 3
            fiju_grid(irg,itg,ipg,1) = fiju_grid(irg,itg,ipg,1) &
            &                        + gmiju(1,ia)*gmiju(2,ib)*fijdc(ia,ib)
            fiju_grid(irg,itg,ipg,2) = fiju_grid(irg,itg,ipg,2) & 
            &                        + gmiju(1,ia)*gmiju(3,ib)*fijdc(ia,ib)
            fiju_grid(irg,itg,ipg,3) = fiju_grid(irg,itg,ipg,3) & 
            &                        + gmiju(2,ia)*gmiju(3,ib)*fijdc(ia,ib)
            fijdu_grid(irg,itg,ipg,ia,ib) = 0.0d0
            do ic = 1, 3
              fijdu_grid(irg,itg,ipg,ia,ib) = fijdu_grid(irg,itg,ipg,ia,ib) &
              &                             + gmiju(ib,ic)*fijdc(ia,ic)
            end do
          end do
        end do        
!
        fijuc(1,1) = 0.0d0
        fijuc(1,2) = fiju_grid(irg,itg,ipg,1)
        fijuc(1,3) = fiju_grid(irg,itg,ipg,2)
        fijuc(2,1) = - fijuc(1,2)
        fijuc(2,2) = 0.0d0
        fijuc(2,3) = fiju_grid(irg,itg,ipg,3)
        fijuc(3,1) = - fijuc(1,3)
        fijuc(3,2) = - fijuc(2,3)
        fijuc(3,3) = 0.0d0
!
        fid(1) = fxd_grid(irg,itg,ipg) ; fiu(1) = fxu_grid(irg,itg,ipg)
        fid(2) = fyd_grid(irg,itg,ipg) ; fiu(2) = fyu_grid(irg,itg,ipg)
        fid(3) = fzd_grid(irg,itg,ipg) ; fiu(3) = fzu_grid(irg,itg,ipg)
!
        fidfiu_grid(irg,itg,ipg) = fid(1)*fiu(1)+fid(2)*fiu(2)+fid(3)*fiu(3)
        fijfij_grid(irg,itg,ipg) = 0.0d0
        do ib = 1, 3
          do ia = 1, 3
            fijfij_grid(irg,itg,ipg) = fijfij_grid(irg,itg,ipg) &
            &                        + fijdc(ia,ib)*fijuc(ia,ib)
          end do
        end do
!
      end do
    end do
  end do
!
end subroutine faraday_raise_gridpoint_WL
