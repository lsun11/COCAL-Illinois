subroutine allocate_emfield
  use phys_constant, only : long
  use grid_parameter
  use def_emfield
  use def_emfield_derivatives
  use def_faraday_tensor
  use make_array_3d
  use make_array_4d
  use make_array_5d
  implicit none
!
  call alloc_array3d(va , 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(alva,0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(vaxd,0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(vayd,0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(vazd,0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(vaxu,0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(vayu,0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(vazu,0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(jtd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(jxd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(jyd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(jzd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(jtu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(jxu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(jyu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(jzu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(jtuf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(jxuf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(jyuf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(jzuf, 0, nrf, 0, ntf, 0, npf)
!
  call alloc_array3d(fxd, 1, nrg, 1, ntg, 1, npg)
  call alloc_array3d(fyd, 1, nrg, 1, ntg, 1, npg)
  call alloc_array3d(fzd, 1, nrg, 1, ntg, 1, npg)
  call alloc_array3d(fxu, 1, nrg, 1, ntg, 1, npg)
  call alloc_array3d(fyu, 1, nrg, 1, ntg, 1, npg)
  call alloc_array3d(fzu, 1, nrg, 1, ntg, 1, npg)
  call alloc_array4d(fijd,1, nrg, 1, ntg, 1, npg, 1, 3) ! In order of
  call alloc_array4d(fiju,1, nrg, 1, ntg, 1, npg, 1, 3) ! f12, f13, f23
  call alloc_array5d(fijdu, 1, nrg, 1, ntg, 1, npg, 1, 3, 1, 3) 
  call alloc_array3d(fidfiu,1, nrg, 1, ntg, 1, npg)
  call alloc_array3d(fijfij,1, nrg, 1, ntg, 1, npg)
  call alloc_array3d(fxd_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(fyd_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(fzd_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(fxu_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(fyu_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(fzu_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array4d(fijd_grid,0, nrg, 0, ntg, 0, npg, 1, 3) ! In order of
  call alloc_array4d(fiju_grid,0, nrg, 0, ntg, 0, npg, 1, 3) ! f12, f13, f23
  call alloc_array5d(fijdu_grid, 0, nrg, 0, ntg, 0, npg, 1, 3, 1, 3) 
  call alloc_array3d(fidfiu_grid,0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(fijfij_grid,0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(Lie_bFxd, 1, nrg, 1, ntg, 1, npg)
  call alloc_array3d(Lie_bFyd, 1, nrg, 1, ntg, 1, npg)
  call alloc_array3d(Lie_bFzd, 1, nrg, 1, ntg, 1, npg)
!
  call alloc_array4d(pdvaxd, 1, nrg, 1, ntg, 1, npg, 1, 3) 
  call alloc_array4d(pdvayd, 1, nrg, 1, ntg, 1, npg, 1, 3) 
  call alloc_array4d(pdvazd, 1, nrg, 1, ntg, 1, npg, 1, 3) 
  call alloc_array4d(cdvaxd, 1, nrg, 1, ntg, 1, npg, 1, 3) 
  call alloc_array4d(cdvayd, 1, nrg, 1, ntg, 1, npg, 1, 3) 
  call alloc_array4d(cdvazd, 1, nrg, 1, ntg, 1, npg, 1, 3) 
  call alloc_array3d(Lie_bAxd, 1, nrg, 1, ntg, 1, npg)
  call alloc_array3d(Lie_bAyd, 1, nrg, 1, ntg, 1, npg)
  call alloc_array3d(Lie_bAzd, 1, nrg, 1, ntg, 1, npg)
  call alloc_array3d(Lie_bAxd_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(Lie_bAyd_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(Lie_bAzd_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(gLie_bAxu_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(gLie_bAyu_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(gLie_bAzu_grid, 0, nrg, 0, ntg, 0, npg)
end subroutine allocate_emfield
