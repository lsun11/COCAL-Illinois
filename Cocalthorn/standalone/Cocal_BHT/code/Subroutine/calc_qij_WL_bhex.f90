subroutine calc_qij_WL_bhex(iter_count,sou_qij)
  use phys_constant, only  : long, pi
  use grid_parameter
  use def_bh_parameter
  use def_metric_hij
  use def_kerr_schild
  use interface_IO_output_1D_general
  use interface_poisson_solver_1bh
  use make_array_3d
  use make_array_2d
  implicit none
  real(long), pointer :: sou_qij(:,:,:,:)
  real(long), pointer :: sou1(:,:,:), sou2(:,:,:), sou3(:,:,:)
  real(long), pointer :: sou4(:,:,:), sou5(:,:,:), sou6(:,:,:)
  real(long), pointer :: sou_bhsurf(:,:), dsou_bhsurf(:,:)
  real(long), pointer :: sou_outsurf(:,:), dsou_outsurf(:,:)
  character(len=2) :: chgreen, chpa, chpB
  character(30) :: char1, char2, char3, char4, char5
  integer :: its, ips, irg, itg, ipg, iter_count

  its = 1
  ips = 1

  call alloc_array3d(sou1,0,nrg,0,ntg,0,npg)
  call alloc_array3d(sou2,0,nrg,0,ntg,0,npg)
  call alloc_array3d(sou3,0,nrg,0,ntg,0,npg)
  call alloc_array3d(sou4,0,nrg,0,ntg,0,npg)
  call alloc_array3d(sou5,0,nrg,0,ntg,0,npg)
  call alloc_array3d(sou6,0,nrg,0,ntg,0,npg)
  call alloc_array2d(sou_bhsurf,0,ntg,0,npg)
  call alloc_array2d(dsou_bhsurf,0,ntg,0,npg)
  call alloc_array2d(sou_outsurf,0,ntg,0,npg)
  call alloc_array2d(dsou_outsurf,0,ntg,0,npg)


  sou1(1:nrg,1:ntg,1:npg) = sou_qij(1:nrg,1:ntg,1:npg,1)
  sou2(1:nrg,1:ntg,1:npg) = sou_qij(1:nrg,1:ntg,1:npg,2)
  sou3(1:nrg,1:ntg,1:npg) = sou_qij(1:nrg,1:ntg,1:npg,3)
  sou4(1:nrg,1:ntg,1:npg) = sou_qij(1:nrg,1:ntg,1:npg,4)
  sou5(1:nrg,1:ntg,1:npg) = sou_qij(1:nrg,1:ntg,1:npg,5)
  sou6(1:nrg,1:ntg,1:npg) = sou_qij(1:nrg,1:ntg,1:npg,6)


  sou_bhsurf(0:ntg,0:npg) = 0.0d0 ; dsou_bhsurf(0:ntg,0:npg) = 0.0d0
  sou_outsurf(0:ntg,0:npg)= 0.0d0 ; dsou_outsurf(0:ntg,0:npg)= 0.0d0
  chgreen = 'dd'

  call poisson_solver_1bh(chgreen,sou1, &
  &                       sou_bhsurf, dsou_bhsurf, &
  &                       sou_outsurf,dsou_outsurf,qxxd)

  char3 = 'qxxd_xa.txt'
  call IO_output_1D_general(char3, 'g', 'g', qxxd, -1, ntg/2,0)
  char3 = 'qxxd_ya.txt'
  call IO_output_1D_general(char3, 'g', 'g', qxxd, -1, ntg/2,npg/4)
  char3 = 'qxxd_za.txt'
  call IO_output_1D_general(char3, 'g', 'g', qxxd, -1, 0,0)
  char3 = 'qxxd_ar.txt'
  call IO_output_1D_general(char3, 'g', 'g', qxxd, -1, its,ips)


  sou_bhsurf(0:ntg,0:npg) = 0.0d0 ; dsou_bhsurf(0:ntg,0:npg) = 0.0d0
  sou_outsurf(0:ntg,0:npg)= 0.0d0 ; dsou_outsurf(0:ntg,0:npg)= 0.0d0
  chgreen = 'dd'

  call poisson_solver_1bh(chgreen,sou2, &
  &                       sou_bhsurf, dsou_bhsurf, &
  &                       sou_outsurf,dsou_outsurf,qxyd)

  char3 = 'qxyd_xa.txt'
  call IO_output_1D_general(char3, 'g', 'g', qxyd, -1, ntg/2,0)
  char3 = 'qxyd_ya.txt'
  call IO_output_1D_general(char3, 'g', 'g', qxyd, -1, ntg/2,npg/4)
  char3 = 'qxyd_za.txt'
  call IO_output_1D_general(char3, 'g', 'g', qxyd, -1, 0,0)
  char3 = 'qxyd_ar.txt'
  call IO_output_1D_general(char3, 'g', 'g', qxyd, -1, its,ips)


  sou_bhsurf(0:ntg,0:npg) = 0.0d0 ; dsou_bhsurf(0:ntg,0:npg) = 0.0d0
  sou_outsurf(0:ntg,0:npg)= 0.0d0 ; dsou_outsurf(0:ntg,0:npg)= 0.0d0
  chgreen = 'dd'

  call poisson_solver_1bh(chgreen,sou3, &
  &                       sou_bhsurf, dsou_bhsurf, &
  &                       sou_outsurf,dsou_outsurf,qxzd)

  char3 = 'qxzd_xa.txt'
  call IO_output_1D_general(char3, 'g', 'g', qxzd, -1, ntg/2,0)
  char3 = 'qxzd_ya.txt'
  call IO_output_1D_general(char3, 'g', 'g', qxzd, -1, ntg/2,npg/4)
  char3 = 'qxzd_za.txt'
  call IO_output_1D_general(char3, 'g', 'g', qxzd, -1, 0,0)
  char3 = 'qxzd_ar.txt'
  call IO_output_1D_general(char3, 'g', 'g', qxzd, -1, its,ips)


  sou_bhsurf(0:ntg,0:npg) = 0.0d0 ; dsou_bhsurf(0:ntg,0:npg) = 0.0d0
  sou_outsurf(0:ntg,0:npg)= 0.0d0 ; dsou_outsurf(0:ntg,0:npg)= 0.0d0
  chgreen = 'dd'

  call poisson_solver_1bh(chgreen,sou4, &
  &                       sou_bhsurf, dsou_bhsurf, &
  &                       sou_outsurf,dsou_outsurf,qyyd)

  char3 = 'qyyd_xa.txt'
  call IO_output_1D_general(char3, 'g', 'g', qyyd, -1, ntg/2,0)
  char3 = 'qyyd_ya.txt'
  call IO_output_1D_general(char3, 'g', 'g', qyyd, -1, ntg/2,npg/4)
  char3 = 'qyyd_za.txt'
  call IO_output_1D_general(char3, 'g', 'g', qyyd, -1, 0,0)
  char3 = 'qyyd_ar.txt'
  call IO_output_1D_general(char3, 'g', 'g', qyyd, -1, its,ips)


  sou_bhsurf(0:ntg,0:npg) = 0.0d0 ; dsou_bhsurf(0:ntg,0:npg) = 0.0d0
  sou_outsurf(0:ntg,0:npg)= 0.0d0 ; dsou_outsurf(0:ntg,0:npg)= 0.0d0
  chgreen = 'dd'

  call poisson_solver_1bh(chgreen,sou5, &
  &                       sou_bhsurf, dsou_bhsurf, &
  &                       sou_outsurf,dsou_outsurf,qyzd)

  char3 = 'qyzd_xa.txt'
  call IO_output_1D_general(char3, 'g', 'g', qyzd, -1, ntg/2,0)
  char3 = 'qyzd_ya.txt'
  call IO_output_1D_general(char3, 'g', 'g', qyzd, -1, ntg/2,npg/4)
  char3 = 'qyzd_za.txt'
  call IO_output_1D_general(char3, 'g', 'g', qyzd, -1, 0,0)
  char3 = 'qyzd_ar.txt'
  call IO_output_1D_general(char3, 'g', 'g', qyzd, -1, its,ips)


  sou_bhsurf(0:ntg,0:npg) = 0.0d0 ; dsou_bhsurf(0:ntg,0:npg) = 0.0d0
  sou_outsurf(0:ntg,0:npg)= 0.0d0 ; dsou_outsurf(0:ntg,0:npg)= 0.0d0
  chgreen = 'dd'

  call poisson_solver_1bh(chgreen,sou6, &
  &                       sou_bhsurf, dsou_bhsurf, &
  &                       sou_outsurf,dsou_outsurf,qzzd)

  char3 = 'qzzd_xa.txt'
  call IO_output_1D_general(char3, 'g', 'g', qzzd, -1, ntg/2,0)
  char3 = 'qzzd_ya.txt'
  call IO_output_1D_general(char3, 'g', 'g', qzzd, -1, ntg/2,npg/4)
  char3 = 'qzzd_za.txt'
  call IO_output_1D_general(char3, 'g', 'g', qzzd, -1, 0,0)
  char3 = 'qzzd_ar.txt'
  call IO_output_1D_general(char3, 'g', 'g', qzzd, -1, its,ips)

!
  deallocate(sou1)
  deallocate(sou2)
  deallocate(sou3)
  deallocate(sou4)
  deallocate(sou5)
  deallocate(sou6)
  deallocate(sou_bhsurf)
  deallocate(dsou_bhsurf)
  deallocate(sou_outsurf)
  deallocate(dsou_outsurf)

end subroutine calc_qij_WL_bhex
