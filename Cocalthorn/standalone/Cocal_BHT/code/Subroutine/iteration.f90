subroutine iteration(iter_count)
  use phys_constant, only :  long, nnrg, nntg, nnpg
  use grid_parameter
  use coordinate_grav_r
  use coordinate_grav_phi
  use coordinate_grav_theta
  use weight_midpoint_grav
  use make_array_2d
  use make_array_3d
  use make_array_4d
  use def_metric
  use def_matter
  use interface_interpo_fl2gr
  use interface_poisson_solver
  use interface_sourceterm_HaC
  use interface_sourceterm_trG
  use interface_sourceterm_MoC
  use interface_compute_shift
  use interface_compute_alps2alph
  use interface_hydrostatic_eq
  use interface_calc_surface
  use interface_update_grfield
  use interface_update_matter
  use interface_update_parameter
  use interface_update_surface
  use interface_error_matter
  implicit none
  real(long), pointer :: sou(:,:,:), pot(:,:,:), sou2(:,:,:)
  real(long), pointer :: potf(:,:,:), emd_bak(:,:,:)
  real(long), pointer :: potrs(:,:)
  real(long), pointer :: potx(:,:,:), poty(:,:,:), potz(:,:,:)
  real(long), pointer :: souvec(:,:,:,:)
  real(long) :: work(0:nnrg,0:nntg,0:nnpg)
  real(long) :: error_emd, count
  integer    :: iter_count, flag = 0, hydro_iter = 4
  integer    :: irf, itf, ipf, irg, itg, ipg, ihy
!
  call alloc_array3d(sou,0,nrg,0,ntg,0,npg)
  call alloc_array3d(sou2,0,nrg,0,ntg,0,npg)
  call alloc_array3d(potx,0,nrg,0,ntg,0,npg)
  call alloc_array3d(poty,0,nrg,0,ntg,0,npg)
  call alloc_array3d(potz,0,nrg,0,ntg,0,npg)
  call alloc_array3d(pot,0,nrg,0,ntg,0,npg)
  call alloc_array3d(potf,0,nrf,0,ntf,0,npf)
  call alloc_array3d(emd_bak,0,nrf,0,ntf,0,npf)
  call alloc_array2d(potrs,0,ntf,0,npf)
  call alloc_array4d(souvec,0,nrg,0,ntg,0,npg,1,3)
!
  iter_count = 0
  do
    iter_count = iter_count + 1      
    count = dble(iter_count) 
    conv_gra = dmin1(conv0_gra,conv_ini*count)
    conv_den = dmin1(conv0_den,conv_ini*count)
!
    emdg = 0.0d0
    call interpo_fl2gr(emd, emdg)
    call calc_vector_x_grav(1)
    call calc_vector_x_matter(1)
    call calc_vector_phi_grav(1)
    call calc_vector_phi_matter(1)
    call excurve
! --
    call sourceterm_HaC(sou)
    call poisson_solver(sou,pot)
    pot = pot + 1.0d0
    call update_grfield(pot,psi,conv_gra)
    call update_parameter(conv_den)
! --
!
    call sourceterm_trG(sou)
    call poisson_solver(sou,pot)
    pot = pot + 1.0d0
    call compute_alps2alph(pot,psi)
    call update_grfield(pot,alph,conv_gra)
    call update_parameter(conv_den)
! --
    call sourceterm_MoC(souvec,sou)
    work(1:nrg, 1:ntg, 1:npg) = souvec(1:nrg, 1:ntg, 1:npg, 1)
    sou2(1:nrg, 1:ntg, 1:npg) = work(1:nrg, 1:ntg, 1:npg)
    call poisson_solver(sou2,potx)
    work(1:nrg, 1:ntg, 1:npg) = souvec(1:nrg, 1:ntg, 1:npg, 2)
    sou2(1:nrg, 1:ntg, 1:npg) = work(1:nrg, 1:ntg, 1:npg)
    call poisson_solver(sou2,poty)
    work(1:nrg, 1:ntg, 1:npg) = souvec(1:nrg, 1:ntg, 1:npg, 3)
    sou2(1:nrg, 1:ntg, 1:npg) = work(1:nrg, 1:ntg, 1:npg)
    call poisson_solver(sou2,potz)
    call poisson_solver(sou,pot)
    work(  0:nrg, 0:ntg, 0:npg)    = potx(0:nrg, 0:ntg, 0:npg)
    souvec(0:nrg, 0:ntg, 0:npg, 1) = work(0:nrg, 0:ntg, 0:npg)
    work(  0:nrg, 0:ntg, 0:npg)    = poty(0:nrg, 0:ntg, 0:npg)
    souvec(0:nrg, 0:ntg, 0:npg, 2) = work(0:nrg, 0:ntg, 0:npg)
    work(  0:nrg, 0:ntg, 0:npg)    = potz(0:nrg, 0:ntg, 0:npg)
    souvec(0:nrg, 0:ntg, 0:npg, 3) = work(0:nrg, 0:ntg, 0:npg)
    call compute_shift(potx, poty, potz, souvec, pot)
    call update_grfield(potx,bvxd,conv_gra)
    call update_grfield(poty,bvyd,conv_gra)
    call update_grfield(potz,bvzd,conv_gra)
    call update_parameter(conv_den)
! -- Hydro equations.
do ihy = 1, hydro_iter
    call hydrostatic_eq(potf)
    call calc_surface(potrs, potf)
    emd_bak = emd 
    call update_matter(potf,emd,conv_den)
    call update_surface(potrs,rs,conv_den)
    call reset_fluid
    call update_parameter(conv_den)
    call calc_vector_x_matter(1)
    call calc_vector_phi_matter(1)
! -- calculate error and print it out.
if (ihy.eq.1) then
    call error_matter(emd,emd_bak,error_emd,flag)
    call printout_error(iter_count,error_emd)
end if
end do
!
    if (flag == 0) exit
    if (iter_count >= iter_max) exit
    flag = 0
!
  end do
!
  deallocate(sou)
  deallocate(sou2)
  deallocate(potx)
  deallocate(poty)
  deallocate(potz)
  deallocate(pot)
  deallocate(potf)
  deallocate(emd_bak)
  deallocate(potrs)
  deallocate(souvec)
end subroutine iteration
