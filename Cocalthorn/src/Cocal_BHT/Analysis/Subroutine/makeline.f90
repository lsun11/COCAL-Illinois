subroutine makeline
  use phys_constant, only : long
  use def_metric
  use grid_parameter, only  : nrg, ntg, npg
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig 
  implicit none
  integer :: ok, ir0,it0,ip0, irg, itg, ipg
  character(10) :: char1, char2, char3, char11, char22, char33
  character(30) :: char_file
!
  open(11,file='line.dat',status='unknown')
  ok=0
  do while (ok==0)
    read(11,'(5i5)',iostat=ok) ir0,it0,ip0
!   close(11)
  
    write(char1, '(i5)') ir0
    write(char2, '(i5)') it0
    write(char3, '(i5)') ip0
    char11 = adjustl(char1)
    char22 = adjustl(char2)
    char33 = adjustl(char3)

    if      (ir0.eq.(-1)) then
      char_file = 'line_it' // trim(char22) // '_ip' // trim(char33) // '.txt'
    else if (it0.eq.(-1)) then
      char_file = 'line_ir' // trim(char11) // '_ip' // trim(char33) // '.txt'
    else if (ip0.eq.(-1)) then
      char_file = 'line_ir' // trim(char11) // '_it' // trim(char22) // '.txt'
    endif
!
!   If ir0=-1 prints a radial line (constant theta,phi)
!   If it0=-1 prints a theta line  (constant r,phi)
!   If ip0=-1 prints a phi line    (constant r,theta)
!
    open(12, file=char_file, status='unknown')
    if      (ir0.eq.(-1).and.it0.le.ntg.and.ip0.le.npg) then
      do irg = 0, nrg
        write(12,'(1p,6e20.12)') rg(irg),  psi(irg,it0,ip0)  &
        &                               , alph(irg,it0,ip0)  &
        &                               , bvxd(irg,it0,ip0)  &
        &                               , bvyd(irg,it0,ip0)  &
        &                               , bvzd(irg,it0,ip0)
      end do
    else if (it0.eq.(-1).and.ir0.le.nrg.and.ip0.le.npg) then
      do itg = 0, ntg
        write(12,'(1p,6e20.12)') thg(itg),  psi(ir0,itg,ip0)  &
        &                                , alph(ir0,itg,ip0)  &
        &                                , bvxd(ir0,itg,ip0)  &
        &                                , bvyd(ir0,itg,ip0)  &
        &                                , bvzd(ir0,itg,ip0)
      end do    
    else if (ip0.eq.(-1).and.it0.le.ntg.and.ir0.le.nrg) then
      do ipg = 0, ntg
        write(12,'(1p,6e20.12)') phig(ipg),  psi(ir0,it0,ipg)  &
        &                                 , alph(ir0,it0,ipg)  &
        &                                 , bvxd(ir0,it0,ipg)  &
        &                                 , bvyd(ir0,it0,ipg)  &
        &                                 , bvzd(ir0,it0,ipg)
      end do        
    endif
    close(12)
  end do
  close(11)
end subroutine makeline
