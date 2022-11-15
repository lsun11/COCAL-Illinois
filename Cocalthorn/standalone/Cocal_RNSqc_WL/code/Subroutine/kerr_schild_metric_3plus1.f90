subroutine kerr_schild_metric_3plus1(x,y,z,psi,alph,bvu,bvd,hijd,hiju)
  use phys_constant,  only : long
  implicit none
  real(long) :: psi, alph, bvu(3), bvd(3), hijd(3,3), hiju(3,3)
  real(long) :: x, y, z, r, lvu(0:3), lvd(0:3), H
  real(long) :: Hp1, H1pH, pp4, pm4, fij(3,3)
!
  call kerr_schild_r_lv_home(x,y,z,r,lvu,lvd,H)
  Hp1  = 1.0d0 + H
  H1pH = H/Hp1
!
  psi  = Hp1**(1.0d0/12.0d0)
  alph = 1.0d0/sqrt(Hp1) 
!
  pp4 = Hp1**(1.0d0/3.0d0)
  pm4 = 1.0d0/pp4
!
  bvu(1:3) = H1pH *lvu(1:3)
  bvd(1:3) = pm4*H*lvd(1:3)
!
  fij(1:3,1:3) = 0.0d0
  fij(1,1) = 1.0d0 ; fij(2,2) = 1.0d0 ; fij(3,3) = 1.0d0
  hijd(1:3,1) = pm4*(fij(1:3,1) +    H*lvd(1:3)*lvd(1)) - fij(1:3,1)
  hijd(1:3,2) = pm4*(fij(1:3,2) +    H*lvd(1:3)*lvd(2)) - fij(1:3,2)
  hijd(1:3,3) = pm4*(fij(1:3,3) +    H*lvd(1:3)*lvd(3)) - fij(1:3,3)
  hiju(1:3,1) = pp4*(fij(1:3,1) - H1pH*lvu(1:3)*lvu(1)) - fij(1:3,1)
  hiju(1:3,2) = pp4*(fij(1:3,2) - H1pH*lvu(1:3)*lvu(2)) - fij(1:3,2)
  hiju(1:3,3) = pp4*(fij(1:3,3) - H1pH*lvu(1:3)*lvu(3)) - fij(1:3,3)
!
end subroutine kerr_schild_metric_3plus1
