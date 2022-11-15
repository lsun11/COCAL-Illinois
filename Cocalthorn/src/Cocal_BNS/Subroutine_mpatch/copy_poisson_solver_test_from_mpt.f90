subroutine copy_poisson_solver_test_from_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
  use def_metric
  use def_metric_mpt
  use def_matter
  use def_matter_mpt
  use copy_array_4dto3d_mpt
  use copy_array_3dto2d_mpt
  implicit none
  integer :: impt
!
  call copy_array4dto3d_mpt(impt, psi_, psi, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, alph_, alph, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, emd_, emd, 0, nrf, 0, ntf, 0, npf)
  call copy_array4dto3d_mpt(impt, emdg_, emdg, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto2d_mpt(impt, rs_, rs, 0, ntf, 0, npf)
!
end subroutine copy_poisson_solver_test_from_mpt
