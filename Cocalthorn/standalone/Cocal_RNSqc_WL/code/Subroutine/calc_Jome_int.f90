subroutine calc_Jome_int(typeR,omega,Jome,Jome_int)
  use phys_constant, only  : long
  use def_matter_parameter
  implicit none
  real(long), intent(inout) :: omega, Jome
  real(long), intent(out):: Jome_int
! -- Jome: functional of Omega (specific angular momentum in Newtonian)
! -- Jome_int: integration of Jome w.r.t. omega
  real(long)             :: pmfac
  character(len=5)       :: typeR
!-------------------------------------------------------------
! Assume axisymmetry.  
! Compute only on xz plane then copy to the other phi=const planes.
!
select case (typeR)
  case('type0')
    if (index_DR.ne.2.0d0) then
      Jome_int = A2DR*omega**2*((ome/omega)**index_DR/(2.0d0-index_DR)-0.5d0) &
      &        - index_DR*A2DR*ome**2/(4.0d0-2.0d0*index_DR)
    else ! v-constant law
      Jome_int = - A2DR*ome**2*dlog(ome/omega) + 0.5d0*A2DR*(ome**2-omega**2)
    end if
!
  case('type1','type2')
    if (typeR.eq.'type1') pmfac = -1.0d0
    if (typeR.eq.'type2') pmfac =  1.0d0
    Jome_int = A2DR*(omega-ome)*((1.0d0+index_DR)*omega+ome)  &
    &        *(pmfac*(omega/ome - 1.0d0))**index_DR  &
    &        /((1.0+index_DR)*(2.0+index_DR))
!
  case('OJ2nd')
    Jome_int = &
    &     A2DR*ome**2*((-0.7853981633974483*A2DR)/B2DR +   &
    &    (1.*A2DR**3*Jome*ome**3)/(Jome**4 + A2DR**4*ome**4) +   &
    &    (1.*A2DR**3*Jome**2*ome**2)/  &
    &     (B2DR*Jome**4 + A2DR**4*B2DR*ome**4) +   &
    &    (0.3535533905932738 + (0.5*A2DR)/B2DR)*  &
    &     atan(1. - (1.4142135623730951*Jome)/(A2DR*ome)) +   &
    &    (-0.3535533905932738 + (0.5*A2DR)/B2DR)*  &
    &     atan(1. + (1.4142135623730951*Jome)/(A2DR*ome)) +   &
    &    0.1767766952966369*Log(Jome**2 -   &
    &       1.4142135623730951*A2DR*Jome*ome + A2DR**2*ome**2) -   &
    &    0.1767766952966369*Log(Jome**2 +   &
    &       1.4142135623730951*A2DR*Jome*ome + A2DR**2*ome**2))
!
  case('OJ4th')
    Jome_int = &
    &     (-0.1414213562373095*A2DR*  &
    &     (3.0575518754255118*A2DR**2 + 0.7217900873323894*B2DR**2)*  &
    &     ome**2)/B2DR**2 + 0.05*A2DR*ome**2*  &
    &   ((20.*A2DR**4*Jome**3*ome**2)/  &
    &      (B2DR**2*(Jome**5 + A2DR**5*ome**5)) +   &
    &     (20.*A2DR**4*Jome*ome**4)/(Jome**5 + A2DR**5*ome**5) +   &
    &     (-7.608452130361228 + (4.702282018339785*A2DR**2)/B2DR**2)*  &
    &      atan((0.2628655560595668*  &
    &          (4.*Jome + 1.2360679774997898*A2DR*ome))/(A2DR*ome)) +   &
    &     (4.702282018339785 + (7.608452130361228*A2DR**2)/B2DR**2)*  &
    &      atan((0.42532540417601994*  &
    &          (-4.*Jome + 3.23606797749979*A2DR*ome))/(A2DR*ome)) -   &
    &     4.*Log(Jome + A2DR*ome) -   &
    &     (4.*A2DR**2*Log(Jome + A2DR*ome))/B2DR**2 +   &
    &     3.23606797749979*Log(Jome**2 -   &
    &        1.618033988749895*A2DR*Jome*ome + A2DR**2*ome**2) -   &
    &     (1.2360679774997898*A2DR**2*  &
    &        Log(Jome**2 - 1.618033988749895*A2DR*Jome*ome +   &
    &          A2DR**2*ome**2))/B2DR**2 -   &
    &     1.2360679774997898*  &
    &      Log(Jome**2 + 0.6180339887498949*A2DR*Jome*ome +   &
    &        A2DR**2*ome**2) +   &
    &     (3.23606797749979*A2DR**2*  &
    &        Log(Jome**2 + 0.6180339887498949*A2DR*Jome*ome +   &
    &          A2DR**2*ome**2))/B2DR**2)
!
  case('OJjco')
    Jome = max(Jome,1.0d-30)
    Jome_int = &
    &   (-0.5*Jome*(2.*(1. + index_DRp)**2*Jome**(1. + index_DRp) -   &
    &      2.*A2DR*index_DRp*(2. + index_DRp)*Jome**index_DRp*ome +   &
    &      B2DR**index_DRp*(2. + 3.*index_DRp + index_DRp**2)*Jome*  &
    &       ome**index_DRp))/  &
    &  (A2DR*B2DR**(1.*index_DRp)*(1. + index_DRp)*(2. + index_DRp)*  &
    &    ome**(1.*index_DRp))
end select
!
end subroutine calc_Jome_int
