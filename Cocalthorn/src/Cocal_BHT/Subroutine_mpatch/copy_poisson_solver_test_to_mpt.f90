subroutine copy_poisson_solver_test_to_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
  use def_metric
  use def_metric_mpt
  use def_matter
  use def_matter_mpt
  use copy_array_3dto4d_mpt
  use copy_array_2dto3d_mpt
  implicit none
  integer :: impt
!
  call copy_array3dto4d_mpt(impt, psi, psi_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, alph, alph_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, emd, emd_, 0, nrf, 0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, emdg, emdg_, 0, nrg, 0, ntg, 0, npg)
  call copy_array2dto3d_mpt(impt, rs, rs_, 0, ntf, 0, npf)
!
end subroutine copy_poisson_solver_test_to_mpt
