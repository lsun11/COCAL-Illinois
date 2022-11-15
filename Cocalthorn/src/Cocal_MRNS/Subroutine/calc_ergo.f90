subroutine calc_ergo
  use phys_constant, only : long
  use def_metric
  use grid_parameter, only : nrg, ntg, npg
  use def_matter, only : ergoin, ergoout
  implicit none
  integer   :: itg, ipg, irg
  real(long) :: gtt
  do ipg = 0, npg
    do itg = 0, ntg
       ergoin(itg,ipg)=0
       ergoout(itg,ipg)=0
       do irg = 0, nrg
          gtt = psi(irg,itg,ipg)**4*(bvxd(irg,itg,ipg)*bvxd(irg,itg,ipg)+ &
              & bvyd(irg,itg,ipg)*bvyd(irg,itg,ipg)                     + &
              & bvzd(irg,itg,ipg)*bvzd(irg,itg,ipg))                    - &
              & alph(irg,itg,ipg)**2
          if (gtt.gt.0.and.ergoin(itg,ipg).eq.0) then
             ergoin(itg,ipg)=irg
             ergoout(itg,ipg)=irg
          end if
          if (gtt.gt.0.and.ergoin(itg,ipg).gt.0) then
             ergoout(itg,ipg)=irg
          end if
       end do
     end do
   end do
end subroutine calc_ergo
