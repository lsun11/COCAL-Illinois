PROGRAM xy2rp
  IMPLICIT NONE
  real(8)  :: pi = 3.141592653589793d+0
  real(8)  :: g = 6.67428d-08
  real(8)  :: c = 2.99792458d+10
  real(8)  :: msol   = 1.98892d+33
  real(8)  :: solmas = 1.98892d+33

  character(LEN=100) :: inputfile, cmd, tempchar1
!
  real(8), allocatable :: alldata(:,:)
  real(8) :: aa,bb,cc,vio, rr,ph,th, ph_prev
  integer :: stat, i, j, k, nlines, narg, iofile, rflag
!
!  narg=command_argument_count()
!  if (narg.ne.1) then
!    write(6,*) "The arguments must be 1: filename "
!    stop
!  endif
!
!  aa=1.0d0;  bb=0.0d0
!  ph = dmod(2.0d0*pi+datan2(bb,aa),2.0d0*pi)
!  write(6,'(1p,3e23.15)')  aa,bb,180.0d0*ph/pi
!  write(6,*)'---------------------------------------------------------'
!  aa=1.0d0;  bb=1.0d0
!  ph = dmod(2.0d0*pi+datan2(bb,aa),2.0d0*pi)
!  write(6,'(1p,3e23.15)')  aa,bb,180.0d0*ph/pi
!  write(6,*)'---------------------------------------------------------'
!  aa=0.0d0;  bb=1.0d0
!  ph = dmod(2.0d0*pi+datan2(bb,aa),2.0d0*pi)
!  write(6,'(1p,3e23.15)')  aa,bb,180.0d0*ph/pi
!  write(6,*)'---------------------------------------------------------'
!  aa=-1.0d0;  bb=1.0d0
!  ph = dmod(2.0d0*pi+datan2(bb,aa),2.0d0*pi)
!  write(6,'(1p,3e23.15)')  aa,bb,180.0d0*ph/pi
!  write(6,*)'---------------------------------------------------------'
!  aa=-1.0d0;  bb=0.0d0
!  ph = dmod(2.0d0*pi+datan2(bb,aa),2.0d0*pi)
!  write(6,'(1p,3e23.15)')  aa,bb,180.0d0*ph/pi
!  write(6,*)'---------------------------------------------------------'
!  aa=-1.0d0;  bb=-1.0d0
!  ph = dmod(2.0d0*pi+datan2(bb,aa),2.0d0*pi)
!  write(6,'(1p,3e23.15)')  aa,bb,180.0d0*ph/pi
!  write(6,*)'---------------------------------------------------------'
!  aa=0.0d0;  bb=-1.0d0
!  ph = dmod(2.0d0*pi+datan2(bb,aa),2.0d0*pi)
!  write(6,'(1p,3e23.15)')  aa,bb,180.0d0*ph/pi
!  write(6,*)'---------------------------------------------------------'
!  aa=1.0d0;  bb=-1.0d0
!  ph = dmod(2.0d0*pi+datan2(bb,aa),2.0d0*pi)
!  write(6,'(1p,3e23.15)')  aa,bb,180.0d0*ph/pi
!  write(6,*)'---------------------------------------------------------'

!  call get_command_argument(1,inputfile)
  open(3, file="HaC_xy_mpt1.txt", status='unknown')
  open(4, file="HaC_rp_mpt1.txt", status='unknown')
  rflag=0
  do 
    read(3,'(1p,3d23.15)',iostat=iofile) aa, bb, vio
    if (iofile > 0)  then
      write(6,*) "Problem reading file teos_parameter.dat...exiting"
      stop
    else if (iofile < 0) THEN
      exit
    else
      rr = sqrt(aa**2 + bb**2)
      if (rr==0.0d0) then
        if (rflag==0 )  then
          rflag=1
          write(4,'(1p,3e23.10)') 0.0d0, 0.0d0, vio
        end if
        cycle
      end if
      ph = dmod(2.0d0*pi+datan2(bb,aa),2.0d0*pi)
!     th = = atan2(rr,cc)

      if (ph_prev>5 .and. ph==0.0d0)  then
        ph = 2.0d0*pi
      end if

      if (abs(ph-ph_prev)<1.0d-10) then
        if(rr.ne.0.0d0)  write(4,'(1p,4e23.10)') rr, ph, vio
      else
        if(rr.ne.0.0d0)  then
          write(4,'(a1)') '' 
          write(4,'(1p,4e23.10)') rr, ph, vio
        end if
      end if
      ph_prev = ph
    end if
  enddo
  close(4)
  close(3)


END PROGRAM xy2rp
