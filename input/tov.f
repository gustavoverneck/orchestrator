      Program tov
c     By Marcelo Chiapparini, DFT-IF/UERJ, october 1996.
c     First revision: 8/3/1998
c     Second revision: 3/7/2001
c
c     This program solves the Tolman-Oppenheimer-Volkoff equations
c     for the inner structure of a spherical star.
c     The units for the energy density E and pressure P are fm^-4,
c     the unit for the radial coordinate r is km and the unit for
c     the star mass M is the sun mass Msun.
c     Two files are need: tov.inp and tov.dat. In tov.inp the user must
c     give the interval (E0inf,E0sup) of central energy density E0 in
c     which the curves M(E0) and R(E0) will be draw, where M(E0) and
c     R(E0) are the mass and the radius of the neutron star with central
c     energy density E0. Npoints is the number of points (stars) in this
c     interval. If Npoints is set to an integer greater than one, then
c     in one of the output files of the program, tov.out, are written E0,
c     M(E0) and R(E0). The values for the maximum mass Mmax
c     and radius Rmax are written in the other output file, star.out,
c     together with the corresponding central energy density. If Npoints
c     is set to one, then only one star is calculated with central energy
c     density E0inf, its mass M and radius R are written in star.out and
c     the curve (r,E(r),P(r)), for r between 0 and R,  is written in tov.out.
c     The file tov.dat is the user provided file with the equation of
c     state (E,P(E)).  This file must have no more than "dime" lines and
c     MUST be ended with the point (-1.,-1.). This last point is used by the
c     program only for determinate how many lines has the tov.dat file (I
c     have experimented problems with some compilers when the eof is
c     reached in tov.dat, so that I use (-1.,-1.) to end the file. This is
c     a workaround which works :) ).
c     Note that the TOV equations are integrated until a value of the
c     radius (r=R) where the pressure obtained is less than the lower
c     pressure value (pmin=p(1)) in tov.dat, this point defines the boundary
c     of the star.
c     All the internal calculations are performed in double precision,
c     the outputs are given in single precision.
c     This program uses a fifth-order Runge-Kutta algorithm with an adaptative
c     stepsize. It need the subroutines SPLINE, SPLINT, ODEINT, RKQS and
c     RKCK from Numerical Recipes (included). I have added some statements to
c     ODEINT in order to fit the requirements of this program.
c
c     variables declaration
      Implicit None
      Integer dime,i,ndat,npoints,flag,flag1,flag2,j
      Integer nvar,nok,nbad,kmax,kmax1,kount
      Parameter (dime=9000,kmax=3000)
      Real*8 e(dime),p(dime),e2(dime),p2(dime),ep1,ep2,pp1,pp2
      Real*8 rho(dime),b2(dime),rp1,rp2
      Real*8 e0inf,e0sup,e0,de0,ener
      Real*8 pmin,r,m,rmax,mmax,e0rmax,e0mmax,rmmax,mrmax
      Real*8 mrbmax,mb,mbmax
      Real*8 dxsav,dxsav1,xp(kmax),yp(3,kmax)
      Real*8 y(3),r1,r2,eps,h1,hmin
      Common/Cf/e,p,rho,e2,p2,b2,ndat
      Common/Cpmin/pmin
      COMMON /path/ kmax1,kount,dxsav1,xp,yp
      External derivs,rkqs
      Data pp1,pp2,ep1,ep2,dxsav/5.d30,5.d30,5.d30,5.d30,1.d-2/
      Data y(2),y(3),nvar,r1,r2,eps,h1,hmin/0.d0,0.d0,3,0.d0,1.d30,
     ^                                 3.d-16,1.d-4,0.d0/

c      write(*,*) ""
c      Write(*,*) "Calculando eq's de Tolman-Oppenheimer-Volkoff..."
c      Write(*,*) " "

c     data input
c      Open(Unit=1,File='tov.inp')
c            Read(1,*)e0inf        !lower central energy density (fm^-3)
c            Read(1,*)e0sup        !higher central energy density (fm^-3)
c            Read(1,*)npoints      !numbers of stars
c      Close(1)
      Open(Unit=2,File='eos.dat')
      flag=0
      i=1
      Do While (flag.eq.0)
            If(i.gt.dime)Then
                  Write(6,*)'Warning!, too much points in TOV.dat.
     ^It will be used only the firts',dime
                  ndat=i-1
                  flag=1
            Else
                  Read(2,*)rho(i),e(i),p(i)
                  If(e(i).lt.0.d0)Then
                        rho(i)=0.d0
                        e(i)=0.d0
                        p(i)=0.d0
                        ndat=i-1
                        flag=1
                  End If
                  i=i+1
            End If
       End Do
       Close(2)
c	INPUT FROM TOV.INP
       e0inf = minval(e)
       e0sup = maxval(e)
       !write(*,*) "Insira o número de pontos (recomendado 100~1000): "
       !READ(*,*) npoints
       !npoints = 500
       npoints = 800
       pmin=p(1)
       If(pmin.le.0.d0)Then
                     Write(6,*)'Fatal!, negative or zero pmin will fail
     ^the convergence of the program'
                     Stop
       End If
c      data interpolation
       Call spline(e,p,ndat,pp1,pp2,p2)
       Call spline(p,e,ndat,ep1,ep2,e2)
       Call spline(rho,e,ndat,ep1,ep2,b2)
c      open output files
       Open(Unit=3,File='tov.out')
       Open(UNit=4,File='star.out')
c      pass kmax and dxsav to odeint
       kmax1=kmax
       dxsav1=dxsav
c      initialize variables
       rmax=0.d0
       mmax=0.d0
       de0=(e0sup-e0inf)/npoints
c      loop over the total number of stars
       Do i=1,npoints
            e0=e0inf+(i-1)*de0
c           p0 and e0
            Call sfp(e0,y(1),flag1)
            If(flag1.eq.0)Then
c                 integrate the TOV equations
                  Call odeint(y,nvar,r1,r2,eps,h1,hmin,nok,nbad,
     &                        derivs,rkqs)
c                 m and r
                  r=xp(kount)
                  m=yp(2,kount)
                  mb=yp(3,kount)
c                 output statements
                  If(npoints.eq.1)Then
c                       output statement for a single star
                        Write(4,*)' --------------'
                        Write(4,*)' - E0 [fm^-4] =',Sngl(e0)
                        Write(4,*)' - M  [Msun]  =',Sngl(m)
                        Write(4,*)' - R  [Km]    =',Sngl(r)
                        Write(4,*)' - Mb  [Msun]  =',Sngl(mb)
                        Write(4,*)' --------------'
                        Do j=1,kount
                              Call SFE(yp(1,j),ener,flag2)
                              Write(3,10)Sngl(xp(j)),Sngl(ener),
     &                                   Sngl(yp(1,j))
                        End do
                  Else
c                       output statement for a set of stars
                        Write(3,10)Sngl(e0),Sngl(m),Sngl(r),Sngl(mb)
c                       choose the maximum radius and mass
                        If(r.gt.rmax)Then
                              rmax=r
                              mrmax=m
                              mrbmax=mb
                              e0rmax=e0
                        End If
                        If(m.gt.mmax)Then
                              mmax=m
                              mbmax=mb
                              rmmax=r
                              e0mmax=e0
                        End If
                  End If
            End If
      End Do
c     write maximum mass and maximum radius
      If(npoints.ne.1)Then
             Write(4,*)' ----------------'
             Write(4,*)' - Mmax [Msun]  =',Sngl(mmax)
             Write(4,*)' - Mbmax [Msun]  =',Sngl(mbmax)
             Write(4,*)'    R   [km]    =',Sngl(rmmax)
             Write(4,*)'    E0  [fm^-4] =',Sngl(e0mmax)
             Write(4,*)' - Rmax [km]    =',Sngl(rmax)
             Write(4,*)'    M   [Msun]  =',Sngl(mrmax)
             Write(4,*)'    Mb   [Msun]  =',Sngl(mrbmax)
             Write(4,*)'    E0  [fm^-4] =',Sngl(e0rmax)
             Write(4,*)' ----------------'
      End If
c     format statement, close output files and end
10    format(11(1x,e25.15))
      Close(3)
      Close(4)
      End
c----------------------------------------------------------------------
c     my subroutines
c----------------------------------------------------------------------
      Subroutine derivs(r,y,dydr)
      Implicit None
      Integer flag
      Real*8 r,y(3),dydr(3)
      Real*8 p,e,m,dpdr,dmdr,dadr,baryon,pi
      Real*8 ms,gms
      Data ms,gms/5660.57d0,1.47556d0/
      Data pi/3.1415926535897d0/

      If(r.eq.0.d0)Then
            dydr(1)=0.d0
            dydr(2)=0.d0
            dydr(3)=0.d0
      Else
            p=y(1)
            Call sfe(p,e,flag)
            m=y(2)
            Call sfb(e,baryon,flag)
            baryon=baryon*939.d0/197.33d0
            dpdr=-(e+p)*(m+4.d0*pi*r**3*p/ms)/(r**2/gms-2.d0*r*m)
            dmdr=4.d0*pi*r**2*e/ms
            dadr=4.d0*pi*(1.d0-2.d0*m/(r/gms))**(-0.5d0)*r*r*baryon/ms
            dydr(1)=dpdr
            dydr(2)=dmdr
            dydr(3)=dadr
c            write(6,*)y(1),y(2),y(3)
      End If

      Return
      End
c----------------------------------------------------------------------
      Subroutine sfb(ex,fb,flag)
      Implicit None
      Integer dime
      Parameter (dime=9000)
      Integer ndat,flag
      Real*8 e(dime),p(dime),e2(dime),p2(dime)
      Real*8 rho(dime),b2(dime)
      Real*8 ex,fb
      Common/Cf/e,p,rho,e2,p2,b2,ndat
      flag=0
      If(ex.lt.e(1).or.ex.gt.e(ndat))Then
c            Write(6,*)'Fatal: attemp to use a central energy out of
c     ^the bounds of tov.inp!'
c            Write(6,*)'E0=',ex,'Emin=',e(1),'Emax=',e(ndat)
            flag=1
      Else
            Call splint(e,rho,b2,ndat,ex,fb)
      End If

      Return
      End
c----------------------------------------------------------------------
      Subroutine sfp(ex,fp,flag)
      Implicit None
      Integer dime
      Parameter (dime=9000)
      Integer ndat,flag
      Real*8 e(dime),p(dime),e2(dime),p2(dime)
      Real*8 rho(dime),b2(dime)
      Real*8 ex,fp
      Common/Cf/e,p,rho,e2,p2,b2,ndat

      flag=0
      If(ex.lt.e(1).or.ex.gt.e(ndat))Then
c            Write(6,*)'Fatal: attemp to use a central energy out of
c     ^the bounds of tov.inp!'
c            Write(6,*)'E0=',ex,'Emin=',e(1),'Emax=',e(ndat)
            flag=1
      Else
            Call splint(e,p,p2,ndat,ex,fp)
      End If

      Return
      End
c----------------------------------------------------------------------
      Subroutine sfe(px,fe,flag)
      Implicit None
      Integer dime
      Parameter (dime=9000)
      Real*8 e(dime),p(dime),e2(dime),p2(dime)
      Real*8 rho(dime),b2(dime)
      Integer ndat,flag
      Real*8 px,fe
      Common/Cf/e,p,rho,e2,p2,b2,ndat

      flag=0
      If(px.lt.p(1))Then
            fe=0.d0
            Write(6,*)'Warning: attemp to use a pressure less than
     ^Pmin in tov.inp!'
            Write(6,*)'P=',px,'Pmin=',p(1)
            flag=1
      Else
            Call splint(p,e,e2,ndat,px,fe)
      End If

      Return
      End
c----------------------------------------------------------------------
c     Numerical Recipes subroutines
c----------------------------------------------------------------------
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=9000)
      INTEGER i,k
      REAL*8 p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
c----------------------------------------------------------------------
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL*8 x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL*8 a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) stop ! 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.d0
      return
      END
c----------------------------------------------------------------------
      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,
     *rkqs)
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL*8 eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=90000,NMAX=3,KMAXX=3000,TINY=1.d-30)
      INTEGER i,kmax,kount,nstp
      REAL*8 dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
     *yp(NMAX,KMAXX),yscal(NMAX)
      Real*8 Pmin,Pnew                      !MC 10/96
      Common/CPmin/Pmin                     !MC 10/96
      COMMON /path/ kmax,kount,dxsav,xp,yp
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.d0*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
12      continue
        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.d0) h=x2-x
        Pnew=y(1)+h*dydx(1)                 !MC 10/96
        If(Pnew.le.Pmin)Return              !MC 10/96
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.d0)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
        if(abs(hnext).lt.hmin) pause
     *'stepsize smaller than minimum in odeint'
        h=hnext
16    continue
      stop ! 'too many steps in odeint'
      return
      END
c----------------------------------------------------------------------
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER n,NMAX
      DOUBLE PRECISION eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
CU    USES derivs,rkck
      INTEGER i
      DOUBLE PRECISION errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY
     *,PGROW,
     *PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,ERRCON=1.89d-4)
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.d0
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.d0)then
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp),0.1d0*abs(h)),h)
        xnew=x+h
        if(xnew.eq.x)stop
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.d0*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END
c----------------------------------------------------------------------
      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
      INTEGER n,NMAX
      DOUBLE PRECISION h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
CU    USES derivs
      INTEGER i
      DOUBLE PRECISION ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX)
     *,
     *ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,
     *B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,B21=.2d0,B31
     *=3.d0/40.d0,
     *B32=9.d0/40.d0,B41=.3d0,B42=-.9d0,B43=1.2d0,B51=-11.d0/54.d0,B52
     *=2.5d0,
     *B53=-70.d0/27.d0,B54=35.d0/27.d0,B61=1631.d0/55296.d0,B62=175.d0
     */512.d0,
     *B63=575.d0/13824.d0,B64=44275.d0/110592.d0,B65=253.d0/4096.d0,C1
     *=37.d0/378.d0,
     *C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0,DC1=C1-2825.d0
     */27648.d0,
     *DC3=C3-18575.d0/48384.d0,DC4=C4-13525.d0/55296.d0,DC5=-277.d0
     */14336.d0,
     *DC6=C6-.25d0)
      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
11    continue
      call derivs(x+A2*h,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
      call derivs(x+A3*h,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      call derivs(x+A4*h,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      call derivs(x+A5*h,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     *B65*ak5(i))
15    continue
      call derivs(x+A6*h,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     *ak6(i))
17    continue
      return
      END
c----------------------------------------------------------------------
