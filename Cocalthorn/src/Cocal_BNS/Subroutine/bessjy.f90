subroutine bessjy(x,xnu,rj,ry,rjp,ryp)
  use phys_constant, only : pi
  implicit none
  integer :: i, isign, l, nl
  integer, parameter :: maxit = 10000
  real(8) :: rj, rjp, ry, ryp, x, xnu
  real(8), parameter :: eps = 1.0d-16, fpmin = 1.0d-30, xmin = 2.0d0
  real(8) ::  a, b, br, bi, c, cr, ci, d, del, del1, den, di, dlr, dli, dr, e, &
    & f, fact, fact2, fact3, ff, gam, gam1, gam2, gammi, gampl, h, p, pimu, pimu2, q, &
    & r, rjl, rjl1, rjmu, rjp1, rjpl, rjtemp, ry1, rymu, rymup, rytemp, sum, sum1, &
    & temp, w, x2, xi, xi2, xmu, xmu2
!
  if (x.le.0.0d0.or.xnu.lt.0.0d0) stop 'bad arguments in bessjy'
  if (x.lt.xmin) then
    nl = int(xnu + 0.5d0)
  else
    nl = max(0,int(xnu-x+1.5d0))
  endif
  xmu = xnu-nl
  xmu2 = xmu*xmu
  xi = 1.0d0/x
  xi2 = 2.0d0*xi
  w = xi2/pi
  isign = 1
  h = xnu*xi
  if (h.lt.fpmin) h = fpmin
  b = xi2*xnu
  d = 0.0d0
  c = h
  do i = 1,maxit
    b = b+xi2
    d = b-d
    if (abs(d).lt.fpmin) d = fpmin
    c = b-1.0d0/c
    if (abs(c).lt.fpmin) c = fpmin
    d = 1.0d0/d
    del = c*d
    h = del*h
    if (d.lt.0.0d0) isign = -isign
    if (abs(del-1.0d0).lt.eps) exit
    if (i.eq.maxit) stop 'x too large in bessjy; try asymptotic expansion'
  end do
  rjl = isign*fpmin
  rjpl = h*rjl
  rjl1 = rjl
  rjp1 = rjpl
  fact = xnu*xi
  do l = nl,1,-1
    rjtemp = fact*rjl+rjpl
    fact = fact-xi
    rjpl = fact*rjtemp-rjl
    rjl = rjtemp
  end do
  if (rjl.eq.0.d0) rjl = eps
  f = rjpl/rjl
  if (x.lt.xmin) then
    x2 = .5d0*x
    pimu = pi*xmu
    if (abs(pimu).lt.eps) then
      fact = 1.0d0
    else
      fact = pimu/sin(pimu)
    endif
    d = -log(x2)
    e = xmu*d
    if (abs(e).lt.eps) then
      fact2 = 1.0d0
    else
      fact2 = sinh(e)/e
    endif
    call beschb(xmu, gam1, gam2, gampl, gammi)
    ff = 2.0d0/pi*fact*(gam1*cosh(e)+gam2*fact2*d)
    e = exp(e)
    p = e/(gampl*pi)
    q = 1.0d0/(e*pi*gammi)
    pimu2 = 0.50d0*pimu
    if (abs(pimu2).lt.eps) then
      fact3 = 1.0d0
    else
      fact3 = sin(pimu2)/pimu2
    endif
    r = pi*pimu2*fact3*fact3
    c = 1.0d0
    d = -x2*x2
    sum = ff+r*q
    sum1 = p
    do i = 1,maxit
      ff = (i*ff+p+q)/(i*i-xmu2)
      c = c*d/i
      p = p/(i-xmu)
      q = q/(i+xmu)
      del = c*(ff+r*q)
      sum = sum+del
      del1 = c*p-i*del
      sum1 = sum1+del1
      if (abs(del).lt.(1.0d0+abs(sum))*eps) exit
      if (i.eq.maxit) stop 'bessy series failed to converge'
    end do
    rymu = -sum
    ry1 = -sum1*xi2
    rymup = xmu*xi*rymu-ry1
    rjmu = w/(rymup-f*rymu)
  else
    a = .25d0-xmu2
    p = -.5d0*xi
    q = 1.0d0
    br = 2.0d0*x
    bi = 2.0d0
    fact = a*xi/(p*p+q*q)
    cr = br+q*fact
    ci = bi+p*fact
    den = br*br+bi*bi
    dr = br/den
    di = -bi/den
    dlr = cr*dr-ci*di
    dli = cr*di+ci*dr
    temp = p*dlr-q*dli
    q = p*dli+q*dlr
    p = temp
    do i = 2,maxit
      a = a+2*(i-1)
      bi = bi+2.d0
      dr = a*dr+br
      di = a*di+bi
      if (abs(dr)+abs(di).lt.fpmin) dr = fpmin
      fact = a/(cr*cr+ci*ci)
      cr = br+cr*fact
      ci = bi-ci*fact
      if (abs(cr)+abs(ci).lt.fpmin) cr = fpmin
      den = dr*dr+di*di
      dr = dr/den
      di = -di/den
      dlr = cr*dr-ci*di
      dli = cr*di+ci*dr
      temp = p*dlr-q*dli
      q = p*dli+q*dlr
      p = temp
      if (abs(dlr-1.d0)+abs(dli).lt.eps) exit
      if (i.eq.maxit) stop 'cf2 failed in bessjy'
    end do
    gam = (p-f)/q
    rjmu = sqrt(w/((p-f)*gam+q))
    rjmu = sign(rjmu,rjl)
    rymu = rjmu*gam 
    rymup = rymu*(p+q/gam)
    ry1 = xmu*xi*rymu-rymup
  endif
  fact = rjmu/rjl
  rj = rjl1*fact
  rjp = rjp1*fact
  do i=1,nl
    rytemp = (xmu+i)*xi2*ry1-rymu
    rymu = ry1
    ry1 = rytemp
  end do
  ry = rymu
  ryp = xnu*xi*rymu-ry1
  return
end subroutine bessjy  
