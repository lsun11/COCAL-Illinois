!______________________________________________
!
!      READ TYPE OF COCAL ID
!______________________________________________
PROGRAM coc2pri_ini
  implicit none

INTERFACE 
   integer function read_id_type(filename,id_type)
     character(400) :: filename
     character(2) :: id_type
   end function read_id_type

   SUBROUTINE coc2pri_ir(dir_path,Xb,Yb,Zb,Nx,Ny,Nz,&
        gbtt,gbtx,gbty,gbtz,gbxx,gbxy,gbxz,gbyy,gbyz,gbzz,&
        gbxx_t,gbxy_t,gbxz_t,gbyy_t,gbyz_t,gbzz_t,&
        Pb,rhob,vbx,vby,vbz)
     character*400 :: dir_path
     integer :: Nx,Ny,Nz
     real(8) :: gbtt(Nx,Ny,Nz),gbtx(Nx,Ny,Nz),gbty(Nx,Ny,Nz),gbtz(Nx,Ny,Nz)
     real(8) :: gbxx(Nx,Ny,Nz),gbxy(Nx,Ny,Nz),gbxz(Nx,Ny,Nz),gbyy(Nx,Ny,Nz),gbyz(Nx,Ny,Nz),gbzz(Nx,Ny,Nz)
     real(8) :: gbxx_t(Nx,Ny,Nz),gbxy_t(Nx,Ny,Nz),gbxz_t(Nx,Ny,Nz),gbyy_t(Nx,Ny,Nz),gbyz_t(Nx,Ny,Nz),gbzz_t(Nx,Ny,Nz)
     real(8) :: Pb(Nx,Ny,Nz),rhob(Nx,Ny,Nz),vbx(Nx,Ny,Nz),vby(Nx,Ny,Nz),vbz(Nx,Ny,Nz),Xb(Nx),Yb(Ny),Zb(Nz)
   END SUBROUTINE coc2pri_ir
   
   SUBROUTINE coc2pri_co(dir_path,Xb,Yb,Zb,Nx,Ny,Nz,&
        gbtt,gbtx,gbty,gbtz,gbxx,gbxy,gbxz,gbyy,gbyz,gbzz,&
        gbxx_t,gbxy_t,gbxz_t,gbyy_t,gbyz_t,gbzz_t,&
        Pb,rhob,vbx,vby,vbz)
     character*400 :: dir_path
     integer :: Nx,Ny,Nz
     real(8) :: gbtt(Nx,Ny,Nz),gbtx(Nx,Ny,Nz),gbty(Nx,Ny,Nz),gbtz(Nx,Ny,Nz)
     real(8) :: gbxx(Nx,Ny,Nz),gbxy(Nx,Ny,Nz),gbxz(Nx,Ny,Nz),gbyy(Nx,Ny,Nz),gbyz(Nx,Ny,Nz),gbzz(Nx,Ny,Nz)
     real(8) :: gbxx_t(Nx,Ny,Nz),gbxy_t(Nx,Ny,Nz),gbxz_t(Nx,Ny,Nz),gbyy_t(Nx,Ny,Nz),gbyz_t(Nx,Ny,Nz),gbzz_t(Nx,Ny,Nz)
     real(8) :: Pb(Nx,Ny,Nz),rhob(Nx,Ny,Nz),vbx(Nx,Ny,Nz),vby(Nx,Ny,Nz),vbz(Nx,Ny,Nz),Xb(Nx),Yb(Ny),Zb(Nz)
   END SUBROUTINE coc2pri_co
   
   SUBROUTINE coc2pri_sp(dir_path,Xb,Yb,Zb,Nx,Ny,Nz,&
        gbtt,gbtx,gbty,gbtz,gbxx,gbxy,gbxz,gbyy,gbyz,gbzz,&
        gbxx_t,gbxy_t,gbxz_t,gbyy_t,gbyz_t,gbzz_t,&
        Pb,rhob,vbx,vby,vbz)
     character*400 :: dir_path
     integer :: Nx,Ny,Nz
     real(8) :: gbtt(Nx,Ny,Nz),gbtx(Nx,Ny,Nz),gbty(Nx,Ny,Nz),gbtz(Nx,Ny,Nz)
     real(8) :: gbxx(Nx,Ny,Nz),gbxy(Nx,Ny,Nz),gbxz(Nx,Ny,Nz),gbyy(Nx,Ny,Nz),gbyz(Nx,Ny,Nz),gbzz(Nx,Ny,Nz)
     real(8) :: gbxx_t(Nx,Ny,Nz),gbxy_t(Nx,Ny,Nz),gbxz_t(Nx,Ny,Nz),gbyy_t(Nx,Ny,Nz),gbyz_t(Nx,Ny,Nz),gbzz_t(Nx,Ny,Nz)
     real(8) :: Pb(Nx,Ny,Nz),rhob(Nx,Ny,Nz),vbx(Nx,Ny,Nz),vby(Nx,Ny,Nz),vbz(Nx,Ny,Nz),Xb(Nx),Yb(Ny),Zb(Nz)
   END SUBROUTINE coc2pri_sp
   
END INTERFACE

character(400) :: dir_path,dir_path_file
integer        :: dir_path_len,i,j,k
character(2)   :: id_type
integer        :: ierr=0
integer :: Nx=10000,Ny=1,Nz=1
real(8), dimension(10000,1,1) :: gbtt,gbtx,gbty,gbtz
real(8), dimension(10000,1,1) :: gbxx,gbxy,gbxz,gbyy,gbyz,gbzz
real(8), dimension(10000,1,1) :: gbxx_t,gbxy_t,gbxz_t,gbyy_t,gbyz_t,gbzz_t
real(8), dimension(10000,1,1) :: Pb,rhob,vbx,vby,vbz
real(8), dimension(10000) :: Xb
real(8), dimension(1) :: Yb
real(8), dimension(1) :: Zb
real(8) :: xmin=0.0d0, xmax=100.0d0


! real(8) :: gbtt(Nx,Ny,Nz),gbtx(Nx,Ny,Nz),gbty(Nx,Ny,Nz),gbtz(Nx,Ny,Nz)
! real(8) :: gbxx(Nx,Ny,Nz),gbxy(Nx,Ny,Nz),gbxz(Nx,Ny,Nz),gbyy(Nx,Ny,Nz),gbyz(Nx,Ny,Nz),gbzz(Nx,Ny,Nz)
! real(8) :: gbxx_t(Nx,Ny,Nz),gbxy_t(Nx,Ny,Nz),gbxz_t(Nx,Ny,Nz),gbyy_t(Nx,Ny,Nz),gbyz_t(Nx,Ny,Nz),gbzz_t(Nx,Ny,Nz)
! real(8) :: Pb(Nx,Ny,Nz),rhob(Nx,Ny,Nz),vbx(Nx,Ny,Nz),vby(Nx,Ny,Nz),vbz(Nx,Ny,Nz),Xb(Nx),Yb(Ny),Zb(Nz)

do i=1,Nx
   Xb(i)=xmin+(xmax-xmin)*i/Nx
end do
Yb=0.d0
Zb=0.d0

!
!TODO remove this
  !dir_path="/home/astro/mundim/tmp/ET_2014_05_wheeler/Cactus/repos/Cocal/standalone/Cocal/ID_BNS"
  !dir_path="../../standalone/Cocal/ID_BNS"
  !dir_path="/home/ruizm/Code/COCAL/standalone/"
  dir_path="/home/ruizm/Code/COCAL/standalone/CID_irHs3d_adm1.5_K123.6_45km"
  dir_path_len = len(trim(dir_path))
  dir_path_file=dir_path(1:dir_path_len)//"/rnspar_mpt1.dat"
  dir_path_len = len(trim(dir_path_file))

! -- Read ID type
!  ierr = read_id_type(dir_path(1:dir_path_len)//"/rnspar_mpt1.dat",id_type)
  ierr = read_id_type(dir_path_file(1:dir_path_len),id_type)
  if (ierr.ne.0) then
     write(6,'(a50)') "COCAL_ID:: problem reading file rnspar_mpt1.dat."
  else
     write(*,*) "COCAL_ID:: success reading rnspar_mpt1.dat."
  end if

  print *, "******"
  print *, "dir_path", dir_path
  print *, "dir_path_len",dir_path_len
  print *, "dir_path_file",dir_path_file
  print *, "dir_path_len",dir_path_len
  print *, "dir_path(1:dir_path_len)=", dir_path(1:dir_path_len)
  print *, "******"

stop
  select case (id_type)
    case ("CO")
      write(6,*) "COCAL_ID:: Reading corotating BNS ID"
      call coc2pri_co(dir_path(1:dir_path_len),Xb,Yb,Zb,Nx,Ny,Nz,&
        gbtt,gbtx,gbty,gbtz,gbxx,gbxy,gbxz,gbyy,gbyz,gbzz,&
        gbxx_t,gbxy_t,gbxz_t,gbyy_t,gbyz_t,gbzz_t,&
        Pb,rhob,vbx,vby,vbz)
    case ("IR")
      write(6,*) "COCAL_ID:: Reading irrotational BNS ID"
      call coc2pri_ir(dir_path(1:dir_path_len),Xb,Yb,Zb,Nx,Ny,Nz,&
        gbtt,gbtx,gbty,gbtz,gbxx,gbxy,gbxz,gbyy,gbyz,gbzz,&
        gbxx_t,gbxy_t,gbxz_t,gbyy_t,gbyz_t,gbzz_t,&
        Pb,rhob,vbx,vby,vbz)
    case ("SP")
      write(6,*) "COCAL_ID:: Reading spinning BNS ID"
      call coc2pri_sp(dir_path(1:dir_path_len),Xb,Yb,Zb,Nx,Ny,Nz,&
        gbtt,gbtx,gbty,gbtz,gbxx,gbxy,gbxz,gbyy,gbyz,gbzz,&
        gbxx_t,gbxy_t,gbxz_t,gbyy_t,gbyz_t,gbzz_t,&
        Pb,rhob,vbx,vby,vbz)
  end select

  Do k=1,Nz
     Do j=1,Ny
        Do i=1,Nx         
           write(*,*) Xb(i),gbtt(i,j,k),rhob(i,j,k)
        end do
     end do
  end do
  
  !
END PROGRAM coc2pri_ini
