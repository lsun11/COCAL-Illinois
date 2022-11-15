subroutine interpolation_fl2gr_IDall
  use phys_constant,  only : long
  use def_matter,     only : emd,  rhof, utf, uxf, uyf, uzf, &
  &                          emdg, rhog, utg, uxg, uyg, uzg
  use interface_interpo_fl2gr
!
  call interpo_fl2gr(emd, emdg)
  call interpo_fl2gr(rhof,rhog)
  call interpo_fl2gr(utf, utg)
  call interpo_fl2gr(uxf, uxg)
  call interpo_fl2gr(uyf, uyg)
  call interpo_fl2gr(uzf, uzg)
!
end subroutine interpolation_fl2gr_IDall
