subroutine printout_error_all_metric(iter_count,error_psi,ire_psi,ite_psi,ipe_psi, &
&                                               error_alph,ire_alph,ite_alph,ipe_alph )
  implicit none
  real(8) :: error_psi,error_alph
  integer :: iter_count, ire_psi,ite_psi,ipe_psi, ire_alph,ite_alph,ipe_alph
  write(6,'(a1,i2,a15,1p,e14.6,a1,i3,a1,i3,a1,i3,a1,a15,1p,e14.6,a1,i3,a1,i3,a1,i3,a1)') '#', iter_count, &
  &       '     Error psi=', error_psi, '(',ire_psi,',',ite_psi,',',ipe_psi,')', &
  &       '    Error alph=', error_alph, '(',ire_alph,',',ite_alph,',',ipe_alph,')' 
end subroutine printout_error_all_metric
