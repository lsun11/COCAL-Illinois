subroutine disk_dimensions
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrg, ntg, npg
  use def_matter, only  :   emdg
  use def_matter_parameter
  use def_formulation
  use def_bht_parameter
  use def_metric
  use coordinate_grav_r, only : rg
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi, only : sinphig, cosphig
  use make_array_2d
  implicit none
  real(long) :: r1,x1,y1,z1,rin,rout
  real(long) :: emd0, emdgc, ray_width
  integer    :: ipg, itg, irg, iflagin

!
  disk_height=0.0d0
  do ipg = 0,0  ! npg
    do irg = 0, nrg
      emd0 = emdg(irg,0,ipg)*10.0d0

      do itg = 0, ntg/2
        emdgc = emdg(irg,itg,ipg)
        x1 = rg(irg)*sinthg(itg)*cosphig(ipg)
        y1 = rg(irg)*sinthg(itg)*sinphig(ipg)
        z1 = rg(irg)*costhg(itg)
!        write(6,'(a3,2i5,1p,3e23.15)')  "",irg,itg, x1,z1, emdgc
        if (emdgc > emd0 .and. z1>disk_height) then
          disk_xh = rg(irg)*sinthg(itg-1)*cosphig(ipg) 
          disk_yh = rg(irg)*sinthg(itg-1)*sinphig(ipg)
          disk_zh = rg(irg)*costhg(itg-1)
          disk_height = disk_zh
!          write(6,'(a3,2i5,1p,3e23.15)')  "***", irg,itg-1, disk_xh, disk_zh, emdg(irg,itg-1,ipg)
          exit
        end if

      end do
    end do
  end do
!
  disk_width=0.0d0
  do ipg = 0,0  ! npg
    do itg = 0, ntg/2
      rin=0.0d0;   rout=0.0d0;   ray_width=0.0d0;
      iflagin=0
      emd0 = emdg(0,itg,ipg)*10.0d0

      do irg = 1, nrg
        emdgc = emdg(irg,itg,ipg)
        x1 = rg(irg)*sinthg(itg)*cosphig(ipg)
        y1 = rg(irg)*sinthg(itg)*sinphig(ipg)
        z1 = rg(irg)*costhg(itg)
        r1 = rg(irg)
!        write(6,'(a3,2i5,1p,4e23.15)')  "",itg,irg, x1,z1,r1,emdgc
        if (emdgc<emd0 .and. iflagin==1) then
          rout=r1;  disk_xou=x1; disk_you=y1; disk_zou=z1;
          ray_width = rout-rin
          if (ray_width>disk_width)   disk_width = ray_width
!          write(6,'(a3,2i5,1p,5e23.15)')  "",itg,irg, x1,z1,r1,emdgc,ray_width
          exit
        end if
        if (emdgc>emd0 .and. iflagin==0) then
          iflagin=1
          rin = rg(irg-1)
          disk_xin = rg(irg-1)*sinthg(itg)*cosphig(ipg)
          disk_yin = rg(irg-1)*sinthg(itg)*sinphig(ipg)
          disk_zin = rg(irg-1)*costhg(itg)
        end if

      end do
    end do
  end do
!

end subroutine disk_dimensions
