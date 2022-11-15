subroutine sourceterm_qij_WL_bhex(iter_count,sou_qij)
  use phys_constant, only : pi, long
  use def_metric_hij
  use def_metric
  use def_metric_excurve_grid
  use def_kerr_schild
  use grid_parameter, only : nrg, ntg, npg
  use interface_grgrad_midpoint
  use interface_interpo_linear_type0
  use interface_dadbscalar_type3_bhex
  use make_array_2d
  use make_array_3d
  use make_array_5d
  use interface_IO_output_1D_general

  implicit none
  real(long), pointer :: sou_qij(:,:,:,:)
  real(long), pointer :: work(:,:,:), fnc0(:,:,:)
  real(8), pointer :: dfdx(:,:,:), dfdy(:,:,:), dfdz(:,:,:)
  real(8), pointer :: gradk(:,:,:,:,:), gradb(:,:,:,:,:)
  real(8) :: dabalps2(1:3,1:3), daalps2(3), liebaij(3,3)
  real(8) :: dbv(1:3,1:3), daij(1:6,1:3), aij(1:3,1:3)
  real(8) :: gamu(1:3,1:3), gamd(1:3,1:3)
  real(long) :: d2hxx(3,3), d2hxy(3,3), d2hxz(3,3), &
     &          d2hyy(3,3), d2hyz(3,3), d2hzz(3,3)
  character(30) :: char1, char2, char3, char4, char5
  real(long) :: psim,alphm,bvxgc,bvygc,bvzgc,ps4oal, trliebaij, &
  &          hhxxu, hhxyu, hhxzu, hhyxu, hhyyu, hhyzu, &
  &          hhzxu, hhzyu, hhzzu, &
  &          hhxxd, hhxyd, hhxzd, hhyxd, hhyyd, hhyzd, &
  &          hhzxd, hhzyd, hhzzd, sgam, san, tfliebaij
  integer :: iter_count
  integer :: ipg, irg, itg, i,j,k,m, ia,ib,ic, ii
!
! 
  call alloc_array3d(dfdx,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfdy,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfdz,1,nrg,1,ntg,1,npg)
  call alloc_array3d(work,0,nrg,0,ntg,0,npg)
  call alloc_array3d(fnc0,0,nrg,0,ntg,0,npg)
  call alloc_array5d(gradk,1,nrg,1,ntg,1,npg,1,6,1,3)
  call alloc_array5d(gradb,1,nrg,1,ntg,1,npg,1,3,1,3)


  san = 1.0d0/3.0d0

  do ii = 1, 3
! ### bvx, bvy, bvz are non-rotating shift!!
    if (ii == 1) call grgrad_midpoint(bvxu,dfdx,dfdy,dfdz)
    if (ii == 2) call grgrad_midpoint(bvyu,dfdx,dfdy,dfdz)
    if (ii == 3) call grgrad_midpoint(bvzu,dfdx,dfdy,dfdz)
    gradb(1:nrg,1:ntg,1:npg,ii,1) = dfdx(1:nrg,1:ntg,1:npg)
    gradb(1:nrg,1:ntg,1:npg,ii,2) = dfdy(1:nrg,1:ntg,1:npg)
    gradb(1:nrg,1:ntg,1:npg,ii,3) = dfdz(1:nrg,1:ntg,1:npg)
  end do

  do ic = 1, 6
    ia = 1 + ic/4 + ic/6
    ib = ic - (ic/4)*2 - ic/6
    fnc0(0:nrg,0:ntg,0:npg) = tfkij_grid(0:nrg,0:ntg,0:npg,ia,ib)

    call grgrad_midpoint(fnc0,dfdx,dfdy,dfdz)
    gradk(1:nrg,1:ntg,1:npg,ic,1) = dfdx(1:nrg,1:ntg,1:npg)
    gradk(1:nrg,1:ntg,1:npg,ic,2) = dfdy(1:nrg,1:ntg,1:npg)
    gradk(1:nrg,1:ntg,1:npg,ic,3) = dfdz(1:nrg,1:ntg,1:npg)
  end do

  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        call interpo_linear_type0(psim, psi, irg,itg,ipg)
        call interpo_linear_type0(alphm,alph,irg,itg,ipg)
        call interpo_linear_type0(bvxgc,bvxu,irg,itg,ipg)
        call interpo_linear_type0(bvygc,bvyu,irg,itg,ipg)
        call interpo_linear_type0(bvzgc,bvzu,irg,itg,ipg)
        ps4oal = psim**4/alphm 
    
        do ia = 1, 3
          aij(ia,1:3) = tfkij(irg,itg,ipg,ia,1:3)
          dbv(ia,1:3) = gradb(irg,itg,ipg,ia,1:3)
        end do
        daij(1:6,1:3) = gradk(irg,itg,ipg,1:6,1:3)

        do ic = 1, 6
          ia = 1 + ic/4 + ic/6
          ib = ic - (ic/4)*2 - ic/6
!
          liebaij(ia,ib) = &
          &     bvxgc*daij(ic,1) + bvygc*daij(ic,2) + bvzgc*daij(ic,3)    &
          &   + aij(ia,1)*dbv(1,ib)+aij(ia,2)*dbv(2,ib)+aij(ia,3)*dbv(3,ib) &
          &   + aij(1,ib)*dbv(1,ia)+aij(2,ib)*dbv(2,ia)+aij(3,ib)*dbv(3,ia)
        end do
        liebaij(2,1) = liebaij(1,2)
        liebaij(3,1) = liebaij(1,3)
        liebaij(3,2) = liebaij(2,3)

        call interpo_linear_type0(hhxxd,hxxd,irg,itg,ipg)
        call interpo_linear_type0(hhxyd,hxyd,irg,itg,ipg)
        call interpo_linear_type0(hhxzd,hxzd,irg,itg,ipg)
        call interpo_linear_type0(hhyyd,hyyd,irg,itg,ipg)
        call interpo_linear_type0(hhyzd,hyzd,irg,itg,ipg)
        call interpo_linear_type0(hhzzd,hzzd,irg,itg,ipg)
        hhyxd = hhxyd
        hhzxd = hhxzd
        hhzyd = hhyzd
        gamd(1,1) = hhxxd + 1.0d0
        gamd(1,2) = hhxyd
        gamd(1,3) = hhxzd
        gamd(2,2) = hhyyd + 1.0d0
        gamd(2,3) = hhyzd
        gamd(3,3) = hhzzd + 1.0d0
        gamd(2,1) = gamd(1,2)
        gamd(3,1) = gamd(1,3)
        gamd(3,2) = gamd(2,3)

        call interpo_linear_type0(hhxxu,hxxu,irg,itg,ipg)
        call interpo_linear_type0(hhxyu,hxyu,irg,itg,ipg)
        call interpo_linear_type0(hhxzu,hxzu,irg,itg,ipg)
        call interpo_linear_type0(hhyyu,hyyu,irg,itg,ipg)
        call interpo_linear_type0(hhyzu,hyzu,irg,itg,ipg)
        call interpo_linear_type0(hhzzu,hzzu,irg,itg,ipg)
        hhyxu = hhxyu
        hhzxu = hhxzu
        hhzyu = hhyzu
        gamu(1,1) = hhxxu + 1.0d0
        gamu(1,2) = hhxyu
        gamu(1,3) = hhxzu
        gamu(2,2) = hhyyu + 1.0d0
        gamu(2,3) = hhyzu
        gamu(3,3) = hhzzu + 1.0d0
        gamu(2,1) = gamu(1,2)
        gamu(3,1) = gamu(1,3)
        gamu(3,2) = gamu(2,3)


        trliebaij = 0.0d0
        do ib = 1, 3
          do ia = 1, 3
            trliebaij = trliebaij + gamu(ia,ib)*liebaij(ia,ib)
          end do
        end do

        do ic = 1, 6
          ia = 1 + ic/4 + ic/6
          ib = ic - (ic/4)*2 - ic/6
!
          sgam = san*gamd(ia,ib)
          tfliebaij = liebaij(ia,ib) - sgam*trliebaij

          sou_qij(irg,itg,ipg, ic) = ps4oal*tfliebaij

!          sou_qij(irg,itg,ipg, ic) = gamd(ia,ib)
        end do
 
!        if (irg==-3 .and. itg==4 .and. ipg==5)  then
!          write(6,'(a19,1p,3e23.15)') "At (x,y,z)        =", xyz(1), xyz(2), xyz(3)
!        end if

      end do
    end do
  end do
!

  write(char1, '(i5)') iter_count
  char2 = adjustl(char1)
  work(1:nrg,1:ntg,1:npg) = sou_qij(1:nrg,1:ntg,1:npg, 1)
  char3 = 'souxx_xa' // trim(char2) // '.txt'
  call IO_output_1D_general(char3,'g','m',work,-1, ntg/2, 1)
  char3 = 'souxx_ya' // trim(char2) // '.txt'
  call IO_output_1D_general(char3,'g','m',work,-1, ntg/2,npg/4)
  char3 = 'souxx_za' // trim(char2) // '.txt'
  call IO_output_1D_general(char3,'g','m',work,-1, 1, 1)
!
  work(1:nrg,1:ntg,1:npg) = sou_qij(1:nrg,1:ntg,1:npg, 2)
  char3 = 'souxy_xa' // trim(char2) // '.txt'
  call IO_output_1D_general(char3,'g','m',work,-1, ntg/2, 1)
  char3 = 'souxy_ya' // trim(char2) // '.txt'
  call IO_output_1D_general(char3,'g','m',work,-1, ntg/2,npg/4)
  char3 = 'souxy_za' // trim(char2) // '.txt'
  call IO_output_1D_general(char3,'g','m',work,-1, 1, 1)
!
  work(1:nrg,1:ntg,1:npg) = sou_qij(1:nrg,1:ntg,1:npg, 4)
  char3 = 'souyy_xa' // trim(char2) // '.txt'
  call IO_output_1D_general(char3,'g','m',work,-1, ntg/2, 1)
  char3 = 'souyy_ya' // trim(char2) // '.txt'
  call IO_output_1D_general(char3,'g','m',work,-1, ntg/2,npg/4)
  char3 = 'souyy_za' // trim(char2) // '.txt'
  call IO_output_1D_general(char3,'g','m',work,-1, 1, 1)
!
  work(1:nrg,1:ntg,1:npg) = sou_qij(1:nrg,1:ntg,1:npg, 6)
  char3 = 'souzz_xa' // trim(char2) // '.txt'
  call IO_output_1D_general(char3,'g','m',work,-1, ntg/2, 1)
  char3 = 'souzz_ya' // trim(char2) // '.txt'
  call IO_output_1D_general(char3,'g','m',work,-1, ntg/2,npg/4)
  char3 = 'souzz_za' // trim(char2) // '.txt'
  call IO_output_1D_general(char3,'g','m',work,-1, 1, 1)
!
  deallocate(work)
  deallocate(fnc0)
  deallocate(dfdx)
  deallocate(dfdy)
  deallocate(dfdz)
  deallocate(gradk)
  deallocate(gradb)

end subroutine sourceterm_qij_WL_bhex
