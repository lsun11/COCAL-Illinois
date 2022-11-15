module def_TOV_quantities
  implicit none
  integer, parameter :: neq = 6, neq_wbc=7, neq_iso=10, neq_iti=12
  real(8)  :: y0(neq), yn(neq), work(neq, 2), &
           &  adm, compa, dqdq, dr, dr_back, &
           &  erer, ermx, ermx_MorC, radi_error, h, &
           &  compab, qini, hinib, radi, radib, radiini, &
           &  x0, xe, xn, yr, hini, ene, pre, rho0, &
           &  rhocgs, precgs, epsiloncgs, radicgs, hh, xecgs, &
           &  y0_wbc(neq_wbc), yn_wbc(neq_wbc), work_wbc(neq_wbc, 2), &
           &  y0_iso(neq_iso), yn_iso(neq_iso), work_iso(neq_iso, 2), &
           &  rho0cgs, rhoscgs, enes, ene0cgs, enescgs, eneini, &
           &  rhoini, dpdrho, enecgs,                           &
           &  y0_iti(neq_iti), yn_iti(neq_iti), work_iti(neq_iti, 2), &
           &  dpde, k2

end module def_TOV_quantities
