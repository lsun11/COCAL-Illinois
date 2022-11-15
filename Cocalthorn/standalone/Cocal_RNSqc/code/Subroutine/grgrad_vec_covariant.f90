subroutine grgrad_vec_covariant(index,fnx,fny,fnz,pdfnx,pdfny,pdfnz, &
           &                                      cdfnx,cdfny,cdfnz)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_cristoffel, only : cri
  use make_array_3d
  implicit none
  real(long), pointer :: fnx(:,:,:), fny(:,:,:), fnz(:,:,:)
  real(long), pointer :: pdfnx(:,:,:,:), pdfny(:,:,:,:), pdfnz(:,:,:,:)
  real(long), pointer :: cdfnx(:,:,:,:), cdfny(:,:,:,:), cdfnz(:,:,:,:)
  real(long) :: pdfn(3,3), cdfn(3,3), c(3,3,3)
  real(long) :: vecx, vecy, vecz
  integer :: irg, itg, ipg, ia, ib
  character(len=1) :: index
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        vecx = fnx(irg,itg,ipg)
        vecy = fny(irg,itg,ipg)
        vecz = fnz(irg,itg,ipg)
        do ia = 1, 3
          pdfn(1,ia) = pdfnx(irg,itg,ipg,ia)
          pdfn(2,ia) = pdfny(irg,itg,ipg,ia)
          pdfn(3,ia) = pdfnz(irg,itg,ipg,ia)
          c(ia,1,1) = cri(irg,itg,ipg,ia,1)
          c(ia,1,2) = cri(irg,itg,ipg,ia,2)
          c(ia,1,3) = cri(irg,itg,ipg,ia,3)
          c(ia,2,1) = c(ia,1,2)
          c(ia,2,2) = cri(irg,itg,ipg,ia,4)
          c(ia,2,3) = cri(irg,itg,ipg,ia,5)
          c(ia,3,1) = c(ia,1,3)
          c(ia,3,2) = c(ia,2,3)
          c(ia,3,3) = cri(irg,itg,ipg,ia,6)
        end do
!
        if (index.eq.'d') then 
          do ia = 1, 3
            do ib = 1, 3
              cdfn(ia,ib) = pdfn(ia,ib) - c(1,ib,ia)*vecx &
              &                         - c(2,ib,ia)*vecy &
              &                         - c(3,ib,ia)*vecz
            end do
          end do
        end if
        if (index.eq.'u') then 
          do ia = 1, 3
            do ib = 1, 3
              cdfn(ia,ib) = pdfn(ia,ib) + c(ia,ib,1)*vecx &
              &                         + c(ia,ib,2)*vecy &
              &                         + c(ia,ib,3)*vecz
            end do
          end do
        end if
!
        do ia = 1, 3
          cdfnx(irg,itg,ipg,ia) = cdfn(1,ia)
          cdfny(irg,itg,ipg,ia) = cdfn(2,ia)
          cdfnz(irg,itg,ipg,ia) = cdfn(3,ia)
        end do
!
      end do
    end do
  end do
!
end subroutine grgrad_vec_covariant
