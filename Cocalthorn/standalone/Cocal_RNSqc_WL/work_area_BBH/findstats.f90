PROGRAM findstats
  IMPLICIT NONE
interface
  subroutine read_datafile(inputfile, nlines, colnum, alldata)
    character(LEN=100), intent(in)  :: inputfile
    integer, intent(out)  :: colnum
    integer, intent(out) :: nlines
    real(8), pointer, intent(out) :: alldata(:,:)
  end subroutine read_datafile

  subroutine find_maxmin(nlines, colnum, alldata, maxmin)
    integer, intent(in)  :: colnum
    integer, intent(in) :: nlines
    real(8), pointer, intent(in) :: alldata(:,:)
    real(8), pointer, intent(out) :: maxmin(:)
  end subroutine find_maxmin

  subroutine interpo(nlines, colnum, alldata, interpoint)
    integer, intent(in)  :: colnum
    integer, intent(in) :: nlines
    real(8), pointer, intent(in) :: alldata(:,:)
    real(8), pointer, intent(out) :: interpoint(:)
  end subroutine interpo
end interface

  real(8)  :: pi = 3.141592653589793d+0
  real(8)  :: g = 6.67428d-08
  real(8)  :: c = 2.99792458d+10
  real(8)  :: msol   = 1.98892d+33
  real(8)  :: solmas = 1.98892d+33

  character(LEN=100) :: inputfile, cmd, tempchar1, tempchar2
!
!  real(8), allocatable :: alldata(:,:)
  real(8), pointer :: alldata(:,:), maxmin(:), interpoint(:)
  integer :: i, j, k, nlines, colnum, narg, iact
!
  narg=command_argument_count()
  if (narg.ne.1) then
    write(6,*) "Provide only the data filename that will be processed. "
    stop
  else 
    call get_command_argument(1,inputfile)
!    write(6,*) "Note that the format for TOV file is different than the RNS file. Modify findstats accordingly."
    call read_datafile(inputfile, nlines, colnum, alldata)
  endif
!
  write(6,*)  "1) Find max, min of two columns."
  write(6,*)  "2) Interpolate."
  write(6,'(a21)',advance='no') "Choose action 1 or 2:"
  read (*,'(i5)') iact

  if(iact==1)  then
    call find_maxmin(nlines, colnum, alldata, maxmin)
  else if(iact==2) then
    call interpo(nlines, colnum, alldata, interpoint)
  else
    stop
  end if
!
  deallocate(alldata)
  if(iact==1)  then
    deallocate(maxmin)
  else if (iact==2) then
    deallocate(interpoint)
  end if
END PROGRAM findstats

subroutine read_datafile(inputfile, nlines, colnum, alldata)
  implicit none
  character(LEN=100), intent(in) :: inputfile
  integer, intent(out)  :: colnum
  integer, intent(out) :: nlines
  real(8), pointer, intent(out) :: alldata(:,:)
  integer :: i,j
  character(LEN=100) :: cmd, tempchar1, tempchar2

  nlines = 0
  cmd = "cat " // trim(inputfile) // " | grep '[^ ]' | wc -l > nlines.txt"
  call system(cmd)
  open(1,file='nlines.txt')
  read(1,*) nlines
  write(6,'(a21,i3)') "#Number of lines are ", nlines
  cmd = 'rm nlines.txt'
  call system(cmd)
!
  colnum = 0
  cmd = "awk '{print NF}' " // trim(inputfile) // " | sort -nu > colnum.txt"
  call system(cmd)
  open(1,file='colnum.txt')
  read(1,*) colnum
  write(6,'(a23,i3)') "#Number of columns are ", colnum
  cmd = 'rm colnum.txt'
  call system(cmd)

  if(nlines==0 .or. colnum==0)  then
    write(6,*) "Number of lines, columns:", nlines, colnum, "...exiting"
    stop
  end if
!
  allocate (alldata(nlines,colnum))
!
  open(3, file=inputfile, status='unknown')
  do i=1,nlines
!   This is the format of TOV isotropic output
!    read(3,'(100es14.6)')  (alldata(i,j), j=1,colnum)

!   This is the format of RNS output in printout_physq_plot.f90
!    read(3,'(2e5.0,1p,48e23.15)')  (alldata(i,j), j=1,colnum)
    read(3,*)  (alldata(i,j), j=1,colnum)
  enddo
  close(3)

end subroutine read_datafile

subroutine find_maxmin(nlines, colnum, alldata, maxmin)
  implicit none
  integer, intent(in)  :: colnum
  integer, intent(in) :: nlines
  real(8), pointer, intent(in) :: alldata(:,:)
  real(8), pointer, intent(out) :: maxmin(:)
  character(LEN=100) :: cmd, tempchar1, tempchar2
!
  real(8) :: fmax, xfmax, fmin, xfmin,max_madm, epsiloncgs, precgs
  integer :: stat, i, j, k, xcol, ycol, imax
!
  fmax=-1.0e+20
  fmin= 1.0e+20

  allocate (maxmin(10))

  write(6,'(a39)',advance='no') "Which column for f_max, f_min (y-axis):"
  read (*,'(i5)') ycol
  write(6,'(a39)',advance='no') "Which column for x_max, x_min (x-axis):"
  read (*,'(i5)') xcol

! find the maximum gravitational mass of the sequence and print main quantities
!  max_madm=0.0d0
!  do i=1,nlines
!    if (alldata(i,11).gt.max_madm) then
!      max_madm=alldata(i,11)
!      imax=i
!    endif
!  enddo
!  precgs     = alldata(imax,5)*c**8/(g**3*solmas**2)   ! pre
!  epsiloncgs = alldata(imax,6)*c**6/(g**3*solmas**2)   ! ene
!  write(6,'(8a12)') " adm ", " restmass ", " propermass ", " radius[km] ", " compa ",  &
!     &          " rhocgs ", " precgs ", " epsiloncgs "
!  write(6,*) "$\ M_{ADM}[M_\odot]\ $ & ", "$\ M_0[M_\odot]\ $ & ", "$\ M_p[M_\odot]\ $ & ", &
!     & "$\ R[km]\ $ & ", "$\quad\ M_{ADM}/R\quad $ & " , "$\quad \log(\rho_0)_c\quad $ & ", &
!     & "$\quad \log P_c\quad $ & ", "$\quad \log\rho_c\quad$"

!  write(6,'(f6.3,a3, f6.3,a3, f6.3,a3, f7.3,a3, f7.4,a3, f7.3,a3, f7.3,a3, f7.3)') &
!     &   alldata(imax,11)," & ", alldata(imax,9), " & ", &
!     &   alldata(imax,10), " & ", alldata(imax,15), " & ", alldata(imax,1), " &",  &
!     &  log10(alldata(imax,4)), " & ", log10(precgs), " & ", log10(epsiloncgs)
!
! find the max of ycol (fmax) and its position in xcol (xfmax)
  do i=1,nlines
    if (alldata(i,ycol).gt.fmax) then
      fmax=alldata(i,ycol)
      xfmax=alldata(i,xcol)
    endif
    if (alldata(i,ycol).lt.fmin) then
      fmin=alldata(i,ycol)
      xfmin=alldata(i,xcol)
    endif
  enddo
  write(6,*) '#---------------------------------------------------------'
  write(6,'(a15,i3,a4,es14.6,a6,es14.6)') '#max of column ', ycol,' is ', fmax, ' at x=', xfmax
  write(6,'(a15,i3,a4,es14.6,a6,es14.6)') '#min of column ', ycol,' is ', fmin, ' at x=', xfmin
  write(6,'(4es14.6)') xfmax, fmax, xfmin, fmin
  maxmin(1) = xfmax
  maxmin(2) =  fmax
  maxmin(3) = xfmin
  maxmin(4) =  fmin
!
end subroutine find_maxmin

subroutine interpo(nlines, colnum, alldata, interpoint)
  implicit none
  integer, intent(in)  :: colnum
  integer, intent(in) :: nlines
  real(8), pointer, intent(in) :: alldata(:,:)
  real(8), pointer, intent(inout) :: interpoint(:)
  real(8), external :: fn_lagint, fn_lagint_2nd
  real(8) :: xval, fval, quad, epsiloncgs, precgs
  integer :: i0, i, j, k, xcol, ycol, iphase
  real(8) :: x4(4), f4(4), x2(2), f2(2)

  allocate (interpoint(colnum))

  write(6,'(a50)',advance='no') "Give the value to be searched for (ex. rest mass):"
  read (*,*) xval
  write(6,'(a42)',advance='no') "Which column corresponds to this quantity:"
  read (*,'(i5)') xcol
  if(xcol > colnum) then
    write(6,*) "There are only ", colnum, " columns."
    stop
  end if
!  write(6,'(a38)',advance='no') "Which column to search for the x-axis:"
!  read (*,'(i5)') xcol

  iphase=0
  do i = 2, nlines
    quad = (xval - alldata(i-1,xcol))*(xval - alldata(i,xcol))
    if( quad .le. 0) then
      iphase = i
      exit
    end if
  end do
  write(6,*)   iphase
  if(iphase==0) then
    write(6,'(a43,i3,a1)') "The value entered is out of range of column",xcol,"."
    stop
  end if

  i0 = min0(max0(iphase-2,0),nlines-3)
  x4(1:4) = alldata(i0:i0+3,xcol)
  do i=1,colnum
    if (i==xcol) then
      interpoint(i) = xval
      continue
    end if
    f4(1:4) = alldata(i0:i0+3, i)
    interpoint(i) = fn_lagint(x4,f4,xval)
  end do
!  write(6,'(a3,1p,4e23.15)')  "x4:", (x4(i),i=1,4)
!  write(6,'(a3,1p,4e23.15)')  "f4:", (f4(i),i=1,4)
!  write(6,'(a10,1p,4e23.15)')  "xval,fval:", xval,fval
!  write(6,*)  "-----------------------------------------------------------------------------------------------"  

  write(6,'(1p,100e23.15)')  (interpoint(i),i=1,colnum)

end subroutine interpo

function fn_lagint(x,y,v)
  implicit none
  real(8), intent(in) :: x(4),y(4), v
  real(8) :: dx12, dx13, dx14, dx21, dx23, dx24, dx31, dx32, dx34, &
  &          dx41, dx42, dx43, wex1, wex2, wex3, wex4, &
  &          xv1, xv2, xv3, xv4, fn_lagint
!
  dx12 = x(1) - x(2)
  dx13 = x(1) - x(3)
  dx14 = x(1) - x(4)
  dx23 = x(2) - x(3)
  dx24 = x(2) - x(4)
  dx34 = x(3) - x(4)
  dx21 = - dx12
  dx31 = - dx13
  dx32 = - dx23
  dx41 = - dx14
  dx42 = - dx24
  dx43 = - dx34
  xv1 = v - x(1)
  xv2 = v - x(2)
  xv3 = v - x(3)
  xv4 = v - x(4)
  wex1 = xv2*xv3*xv4/(dx12*dx13*dx14)
  wex2 = xv1*xv3*xv4/(dx21*dx23*dx24)
  wex3 = xv1*xv2*xv4/(dx31*dx32*dx34)
  wex4 = xv1*xv2*xv3/(dx41*dx42*dx43)
!
  fn_lagint = wex1*y(1) + wex2*y(2) + wex3*y(3) + wex4*y(4)
!
end function fn_lagint

function fn_lagint_2nd(x,y,v)
  implicit none
  real(8), intent(in) :: x(2),y(2), v
  real(8) :: dx12, dx21, wex1, wex2, xv1, xv2, fn_lagint_2nd
!
  dx12 = x(1) - x(2)
  dx21 = - dx12
  xv1 = v - x(1)
  xv2 = v - x(2)
  wex1 = xv2/dx12
  wex2 = xv1/dx21
!
  fn_lagint_2nd = wex1*y(1) + wex2*y(2)
!
end function fn_lagint_2nd

