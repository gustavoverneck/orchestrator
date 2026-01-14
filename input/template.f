      program eosfuncionando

      implicit double precision(a-h,k-z)

      logical check
      integer n, store_cnt, jj, kk, max_store, npoint
      parameter (max_store=9000)
      double precision store_data(max_store, 13)
      dimension x(4),fvec(4)
      dimension ml(2),mb(8)
      dimension nb(8),nl(2),fpu(9999),fpd(9999),fe(9999)
      dimension fmu(9999), fsm(9999),fsp(9999),fxm(9999)
      common/cdata/pi,pi2,ml,rb,rc,rxi,gs,gv,gr,qe,b,ammp,ammn,mb
      common/cnbt/nbt,nb,nl
      common/cfvec/fvec
      common/out/mns,efn,efp,efe,efmu,mup,mun,fpu,fpd,fe,fmu
     &,npu,npd,ne,nu,rkfl0,efl0,ml0s, mul0
     &,msms,efsm,musm,fsm,nsm,rkfs0,efs0,ms0s,mus0
     &,msps,efsp,musp,fsp,nsp,rkfx0,efx0,mx0s,mux0
     &,mxms,efxm,muxm,fxm,nxm
      COMMON/lvs/csi
      
      character(len=32) :: csi_string    ! Variável para capturar o argumento
      character(len=32) :: B_string
      character(len=32) :: model_string
      integer :: iostat                  ! Variável para verificar se a leitura foi bem-sucedida

      ! Obtém o argumento da linha de comando
      call get_command_argument(1, csi_string)
      call get_command_argument(2, B_string)
      call get_command_argument(3, model_string)

      ! Converte o argumento de string para double precision
      read(csi_string, *, iostat=iostat) csi
      
      if (iostat /= 0) then
            print *, "Erro ao converter o argumento para número real."
            stop
      endif
      
      read(B_string, *, iostat=iostat) bg

      pi=dacos(-1.d0)
      pi2=pi*pi
      hc=197.326d0		      ! Conversion 197.326 mev.fm = 3.16151093e-17 j.m


      m = 938.9187137765d0          ! Nucleon mass
      denss=0.153d0		      ! Nuclear saturation density n0 = 0.153 fm^−3


c	Lepton's Masses
      ml(1)=0.511d0/m		      ! Normalizesd Electron mass
      ml(2)=105.66d0/m		      ! Normalized Muon mass

c	massas dos barions
      mb(1) = 939.56534623d0/m	! Normalized Neutron mass
      mb(2) = 938.272081323d0/m	! Normalized Proton mass     
      mb(3) = 1116.d0/m 	      ! Normalized Lambda0 mass  
      mb(4) = 1193.d0/m 	      ! Normalized Sigma- mass   
      mb(5) = 1193.d0/m 	      ! Normalized Sigma0 mass     
      mb(6) = 1193.d0/m 	      ! Normalized Sigma+ mass   
      mb(7) = 1318.d0/m 	      ! Normalized Xi- mass   
      mb(8) = 1318.d0/m 	      ! Normalized Xi0 mass

c	Meson's masses
      ms=400.d0/m		            ! Normalized Scalar Meson (sigma) mass
      mv=783.d0/m		            ! Normalized Vector Meson (omega) mass
      mrho=770.d0/m		      ! Normalized Iso-Vector Meson (rho) mass

c---------------------------------------------------------------------
c  PARAMETRIZATION
c---------------------------------------------------------------------
      if (model_string .eq. "GM1") then
            gs = dsqrt(11.785d0) / hc * m
            gv = dsqrt(7.148d0) / hc * m
            gr = dsqrt(4.41d0) / hc * m
            rb = 0.002948d0
            rc = -0.001071d0
            rxi = 0.d0
            xs = 0.7d0
            xv = 0.783d0
      else if (model_string .eq. "GM3") then
            gs = dsqrt(9.927d0) / hc * m
            gv = dsqrt(4.820d0) / hc * m
            gr = dsqrt(4.791d0) / hc * m
            rb = 0.008659d0
            rc = -0.002421d0
            rxi = 0.d0
            xs = 0.7d0
            xv = 0.783d0
      else
            print *, "Erro: Modelo não reconhecido. Use GM1 ou GM3."
            stop
      endif

c---------------------------------------------------------------------
c      B(n) = B_surf + B0{1 - exp[-beta(n/n0)^alpha]}
c---------------------------------------------------------------------

      qe= dsqrt(4.0d0*pi/137.0d0)         ! Electron charge
      bce=ml(1)**2/qe    		      ! Critical magnetic field
      rncm=qe/2.0d0*0.0d0          	      ! Nuclear magneton

      ammp=rncm*(5.5857d0/2.0d0-1.0d0)	! Proton g = 5.5857,  the magnetic moment is μ = 2.7928 nuclear magnetons.
      ammn=-rncm*3.8261d0/2.0d0		! Neutron g = 3.8261,  the magnetic moment is μ = 2.7928 nuclear magnetons.

      be=ml(1)**2          
      b0=bg/4.41d13             	      ! Magnetic field inside the star
      b=b0*bce
      btsl = bg*1.0d-4
      bsurf = 1.0d11		            ! Magnetic field in the surface

      betaa = 1.0d-2		            ! eq  6.32 parameter
      alphaa = 3.0d0		            ! eq  6.32 parameter


      muninf=0.92d0		            ! Inferior limit of Neutron Chemical potential (k^2+ms^2)^1/2   [fm^-3]
      munsup=1.80d0		            ! Superior limit of Neutron Chemical potential (k^2+ms^2)^1/2   [fm^-3]
      npoint=1201		      	      ! Number of points


      dmub=(munsup-muninf)/(npoint-1)	! Finite element of Chemical Potential


      open(unit=4,file='eos.dat')         ! output file
      store_cnt = 0


c     Nonlinear system variables
      x(1)=0.8d0                   	      ! Electronic chemical potential 
      x(2)=0.01d0            		      ! vsigma field
      x(3)=0.8d0                          ! vomega field
      x(4)=0.02d0                         ! vrho field

c---------------------------------------------------------------------
      ! Main loop
      do i=1, npoint
            mun=munsup-(i-1)*dmub

            call broydn(x,4,check)

            if(check) exit

            call mapping(x,mue,vsigma,vomega,vrho)
            mns=1.0d0-vsigma                   	      ! neutron effective mass M*/m = m/m - gs*sigma/m 
            efn=mun-vomega+vrho/2.0d0		      ! energy fermi neutron eq(5.7)  mun = efn + vomega + (-1/2)vrho
            kfnu2=efn**2-(mns-ammn*b)**2		      ! nucleon fermi momenta^2 for spin up

            ! Keep fermi momenta positive
            if(kfnu2.gt.0.0d0) then
                  kfnu=dsqrt(kfnu2)
            else
                  kfnu=0.0d0
            endif
            kfnd2=efn**2-(mns+ammn*b)**2	            !nucleon fermi momenta^2 for spin down

            if(kfnd2.gt.0.0d0) then
                  kfnd=dsqrt(kfnd2)
            else
                  kfnd=0.0d0
            endif

            call eos(mue,kfnu,kfnd,vsigma,vomega,vrho,ener,press)       ! Calculate Equation of State

            nbtd=nbt*(m/197.32d0)**3

            ff1 = 0.333d0*(nb(3)+nb(4)+nb(5)+nb(6)+2*(nb(7)+nb(8)))     ! Strangeness fraction

            emn=(mns*((nb(1)+nb(2))/nbt))*939.d0
            eml=(ml0s*((nb(3))/nbt))*939.d0
            ems=(ms0s*((nb(4)+nb(5)+nb(6))/nbt))*939.d0
            emx=(mx0s*((nb(7)+nb(8))/nbt))*939.d0
            embay=emn+eml+ems+emx

            bdd = bsurf + btsl*(1 - exp(-betaa*(nbtd/denss)**alphaa))   ! magnetic fild density-dependent model
            ebsi = ((bdd)**2)/(8*pi*1.0d-7)                             ! energy of magnetic fild in si
            ebsd = ebsi/(3.161d34)                                      ! energy of magnetic fild in crazy units
            dener=ener/nbt-1.0d0
            
            press=press*m*(m/197.32d0)**3
            ener=ener*m*(m/197.32d0)**3
            ener00 = ener/197.32d0
            press00 = press/197.32d0
            ener = ener/197.326d0 + ebsd		      ! ener_T = ener + B^2/(2*4*pi)
            press = press/197.326d0  + ebsd		! press_T = press + B^2/(2*4*pi)
            nbtd = nbtd*0.153d0

            if (ener .GE. 0.0d0 .and. press .GE. 0.0d0) then
                  store_cnt = store_cnt + 1
                  if (store_cnt .le. max_store) then
                        store_data(store_cnt, 1) = nbtd/denss
                        store_data(store_cnt, 2) = ener
                        store_data(store_cnt, 3) = press
                        store_data(store_cnt, 4) = nl(1)
                        store_data(store_cnt, 5) = nl(2)
                        store_data(store_cnt, 6) = nb(1)
                        store_data(store_cnt, 7) = nb(2)
                        store_data(store_cnt, 8) = nb(3)
                        store_data(store_cnt, 9) = nb(4)
                        store_data(store_cnt, 10) = nb(5)
                        store_data(store_cnt, 11) = nb(6)
                        store_data(store_cnt, 12) = nb(7)
                        store_data(store_cnt, 13) = nb(8)
                  else
                        print *, "Storage limit exceeded"
                  endif
            endif

      end do
      
      do jj = store_cnt, 1, -1
         write(4,10) (store_data(jj, kk), kk=1,13)
      end do
      
      write(4, 10) -1., -1., -1., -1.

      close(2)
      close(3)
      close(4)
10    format(13(1x,e12.5))

      stop
      end
c--------------------------------------------------------------------------
      subroutine mapping(x,mue,vsigma,vomega,vrho)
c     this subroutine sets the potential fields between the physical values
      implicit double precision(a-h,k-z)  
      integer n
      dimension x(4)
      COMMON/lvs/csi

      mue=x(1)                        	!to assure kfe greater than zero
      vsigma=dsin(x(2))**2   		      !mapping vsigma in (0,mb(1))
      vomega=x(3)                         !to assure vomega greater than zero
      vrho=x(4)                       	!no mapping on vrho

      return
      end
c-----------------------------------------------------------------------
      subroutine funcv(n,x,fvec)
c     this subroutine is used by broydn and gives the functions to be zeroed 
c     in the vector fvec
      implicit double precision(a-h,k-z)
      integer n

      dimension x(4),fvec(4),fvec1(4)
      dimension ml(2),mb(8)
      dimension nb(8),nl(2),fpu(9999),fpd(9999),fe(9999)
      dimension fmu(9999), fsm(9999),fsp(9999),fxm(9999)
      common/cdata/pi,pi2,ml,rb,rc,rxi,gs,gv,gr,qe,b,ammp,ammn,mb
      common/cnbt/nbt,nb,nl
      common/cfvec/fvec1
      common/out/mns,efn,efp,efe,efmu,mup,mun,fpu,fpd,fe,fmu
     &,npu,npd,ne,nu,rkfl0,efl0,ml0s, mul0
     &,msms,efsm,musm,fsm,nsm,rkfs0,efs0,ms0s,mus0
     &,msps,efsp,musp,fsp,nsp,rkfx0,efx0,mx0s,mux0
     &,mxms,efxm,muxm,fxm,nxm
      COMMON/lvs/csi

      call mapping(x,mue,vsigma,vomega,vrho)

      fsigma =0.0d0
      fomega =0.0d0
      frho   =0.0d0
      charge =0.0d0                           !electric charge density
      rnumber=0.0d0                           !baryonic number density 

      xs = 0.7d0
      xv = 0.783d0

      mns=1.0d0-vsigma                          ! neutron e proton effective mass
      ml0s = (1116.d0/938.99d0) - xs*vsigma     ! lambda0 efetive mass
      msms = (1193.d0/938.99d0) - xs*vsigma     ! sigma- effective mass
      ms0s = (1193.d0/938.99d0) - xs*vsigma     ! sigma0 effective mass
      msps = (1193.d0/938.99d0) - xs*vsigma     ! sigma+ effective mass
      mxms = (1318.d0/938.99d0) - xs*vsigma     ! xi- effective mass
      mx0s = (1318.d0/938.99d0) - xs*vsigma     ! xi0 effective mass

c	===========================================
c	tratamento energi de fermi neutron
c	spin up
      efn=mun-vomega+vrho/2.0d0                 ! neutron fermi energy
      kfnu2=efn**2-(mns)**2-csi**2
      
      if(kfnu2.gt.0.0d0)then
         kfnu=dsqrt(kfnu2)             ! neutron fermi momenta for spin up
      else
         kfnu=0.0d0
      endif
c	spin down       
      kfnd2=efn**2-(mns)**2
      if(kfnd2.gt.0.0d0)then
         kfnd=dsqrt(kfnd2)                     ! neutron fermi momenta for spin down
      else
         kfnd=0.0d0
      endif

c	definição dos momentos de fermi a partir dos numeros de landau  
      do i=1,9999               !numero relacionado aos niveis de Landau, nu_max
         fe(i)=0.d0
         fpu(i)=0.d0
         fpd(i)=0.d0
         fmu(i)=0.d0
         fsm(i)=0.d0
         fsp(i)=0.d0
         fxm(i)=0.d0
      enddo

c	===========================================
c	tratamento energi de fermi eletron
      ne=0
      do i=1,9999
         fe2=mue**2-(ml(1)**2+2*qe*b*(i-1)+csi**2)
         if(fe2.gt.0.0d0)then
            rkfe=dsqrt(fe2)
         else
            rkfe=0.0d0
         endif
         fe(i)=rkfe
         if(rkfe.eq.0.d0)goto 3
         ne=i
      enddo
 3    continue
      efe=dsqrt(fe(1)**2+ml(1)**2+csi**2)             !eletron fermi energy

c decaimentos de potenciais quimicos
      mup=mun-mue		!decaimento beta-  potenciais quimicos
      mul0 = mun		!decaimento lambda_0 no neutron
      musm = mun+mue		!decaimento sigma- neutron + eletron
      mus0 = mun		!decaimento sigma0 no neutron
      musp = mun-mue		!decaimento sigma+ em neutron e eletron
      muxm = mun+mue		!decaimento xi- em neutron e eletron
      mux0 = mun		!decaimento xi0 no neutron

c	===========================================
c	tratamento energi de fermi proton
c	spin up
      npu=0
      do i=1,9999
         rkfpu2=(mup-vomega-vrho/2.0d0)**2           	!quadrado proton fermi momenta for spin up
     &-(dsqrt(mns**2+2*qe*b*(i-1)+csi**2))**2 			!proton fermi mom squared
         if(rkfpu2.gt.0.0d0)then
            rkfpu=dsqrt(rkfpu2)
         else
            rkfpu=0.0d0
         endif
         fpu(i)=rkfpu
         if(rkfpu.eq.0.d0)goto 4
         npu=i                	!numero de landau para proton up
      enddo
 4    continue                  !proton chemical potential

c      efp=mup-vomega-vrho/2.0d0 	!proton fermi energy
      efp=dsqrt(fpu(1)**2+(mns)**2+csi**2) 	!proton fermi energy

c	===========================================
c	tratamento energi de fermi sigma-
c	spin up xxxxx acho q não incluiu o magneton de bohr aqui
      nsm=0
      do i=1,9999
         rkfsm2=(musm-(xv*vomega)+(xv*vrho))**2
     &-(dsqrt(msms**2+2*qe*b*(i-1)+csi**2))**2 !sigma- fermi mom squared
         if(rkfsm2.gt.0.0d0)then
            rkfsm=dsqrt(rkfsm2)
         else
            rkfsm=0.0d0
         endif
         fsm(i)=rkfsm
         if(rkfsm.eq.0.d0)goto 61
         nsm=i
      enddo
 61    continue                  

      efsm=dsqrt(fsm(1)**2+(msms)**2+csi**2) !sigma- fermi energy

c	===========================================
c	tratamento energi de fermi sigma+
      nsp=0
      do i=1,9999
         rkfsp2=(musp-(xv*vomega)-(xv*vrho))**2
     &-(dsqrt(msps**2+2*qe*b*(i-1)+csi**2))**2 !sigma+ fermi mom squared
         if(rkfsp2.gt.0.0d0)then
            rkfsp=dsqrt(rkfsp2)
         else
            rkfsp=0.0d0
         endif
         fsp(i)=rkfsp
         if(rkfsp.eq.0.d0)goto 62
         nsp=i
      enddo
 62    continue                  

      efsp=dsqrt(fsp(1)**2+(msps)**2+csi**2) !sigma+ fermi energy

c	===========================================
c	tratamento energi de fermi xi-
      nxm=0
      do i=1,9999
         rkfxm2=(muxm-(xv*vomega)+(xv*vrho/2.d0))**2
     &-(dsqrt(mxms**2+2*qe*b*(i-1)+csi**2))**2 !xi- fermi mom squared
         if(rkfxm2.gt.0.0d0)then
            rkfxm=dsqrt(rkfxm2)
         else
            rkfxm=0.0d0
         endif
         fxm(i)=rkfxm
         if(rkfxm.eq.0.d0)goto 63
         nxm=i
      enddo
 63    continue                  

      efxm=dsqrt(fxm(1)**2+(mxms)**2+csi**2) !sigma- fermi energy

c	===========================================
c	tratamento energi de fermi proton down
      npd=0
      do i=1,9999
         rkfpd2=(mup-vomega-vrho/2.0d0)**2
     &-(dsqrt(mns**2+2*qe*b*i+csi**2))**2 !proton fermi mom squared
         if(rkfpd2.gt.0.0d0)then
            rkfpd=dsqrt(rkfpd2)
         else
            rkfpd=0.0d0
         endif
         fpd(i)=rkfpd
         if(rkfpd.eq.0.d0) goto 5
         npd=i
      enddo
 5    continue                  !proton chemical potential

c	===========================================
c	tratamento energi de fermi lambda0
       rkfl02 = (mul0 - xv*vomega)**2 - ml0s**2-csi**2             !lambda0 momento de fermi ao quadrado
       if (rkfl02 .gt. 0.d0) then                             ! condiçaõ de existência
           rkfl0 = dsqrt(rkfl02)           
           else
           rkfl0 = 0.0d0           
       endif
       efl0 = dsqrt(ml0s**2 + rkfl0**2+csi**2)           !lambda0 fermi energy

c	===========================================
c	tratamento energi de fermi sigma0
       rkfs02 = (mus0 - xv*vomega)**2 - ms0s**2-csi**2              !sigma0 momento de fermi ao quadrado
       if (rkfs02 .gt. 0.d0) then                             ! condiçaõ de existência
           rkfs0 = dsqrt(rkfs02)           
           else
           rkfs0 = 0.0d0           
       endif
       efs0 = dsqrt(ms0s**2 + rkfs0**2+csi**2)           !sigma0 fermi energy

c	===========================================
c	tratamento energi de fermi xi0
       rkfx02 = (mux0 - xv*vomega -(xv*vrho/2.d0))**2 - mb(8)**2
     &-csi**2              !xi0 momento de fermi ao quadrado
       if (rkfx02 .gt. 0.d0) then                             ! condiçaõ de existência
           rkfx0 = dsqrt(rkfx02)           
           else
           rkfx0 = 0.0d0           
       endif
       efx0 = dsqrt((mx0s**2 + rkfx0**2+csi**2))           !xi0 fermi energy
          
c	===========================================
c	tratamento energi de fermi muon
      rkfmu2=0.0d0
      nu=0
      do i=1,9999
         rkfmu2=mue**2-(ml(2)**2+2*qe*b*(i-1)+csi**2)
         if(rkfmu2.gt.0.0d0)then
            rkfmu=dsqrt(rkfmu2)
         else
            rkfmu=0.0d0
         endif
         fmu(i)=rkfmu
         if(rkfmu.eq.0.d0)goto 7
         nu=i
      enddo
 7    continue
      efmu=dsqrt(fmu(1)**2+ml(2)**2+csi**2)	!muon fermi energy

c	===========================================
c	Calculos densidades
c	===========================================

c	densidade do neutron

      rhosn=mns*(efn*kfnu
     &-(mns)**2*dlog(dabs((kfnu+efn)/(mns))))*0.5d0

c	densidade do proton up
      rhosp=0.0d0
      do i=1,npu
      mnsu=dsqrt(mns**2+2.0d0*qe*b*(i-1)+csi**2)
      rhosp=rhosp+qe*b*mns*((mnsu)/mnsu
     &*dlog(dabs((fpu(i)+efp)/(mnsu)))) ! proton scalar density
      enddo

c	densidade do proton down
      do i=1,npd
      mnsd=dsqrt(mns**2+2.0d0*qe*b*i+csi**2)
      rhosp=rhosp+qe*b*mns*((mnsd)/mnsd
     &*dlog(dabs((fpd(i)+efp)/(mnsd))))
      enddo
       rhosp = rhosp*0.5d0

c	densidade sigma-
      rhossm=0.0d0
      do i=1,nsm
      msmsu=dsqrt(msms**2+2.0d0*qe*b*(i-1)+csi**2)
      rhossm=rhossm+dg(i)*qe*b*msms*((msmsu)/msmsu
     &*dlog(dabs((fsm(i)+efsm)/(msmsu)))) ! sigma- scalar density
      enddo
       rhossm = rhossm*0.5d0

c	densidade sigma+
      rhossp=0.0d0
      do i=1,nsp
      mspsu=dsqrt(msps**2+2.0d0*qe*b*(i-1)+csi**2)
      rhossp=rhossp+dg(i)*qe*b*msps*((mspsu)/mspsu
     &*dlog(dabs((fsp(i)+efsp)/(mspsu)))) ! sigma+ scalar density
      enddo
       rhossp = rhossp*0.5d0

c	densidade xi-
      rhosxm=0.0d0
      do i=1,nxm
      mxmsu=dsqrt(mxms**2+2.0d0*qe*b*(i-1)+csi**2)
      rhosxm=rhosxm+dg(i)*qe*b*mxms*((mxmsu)/mxmsu
     &*dlog(dabs((fxm(i)+efxm)/(mxmsu)))) ! xi- scalar density
      enddo
       rhosxm = rhosxm*0.5d0

c	densidade escalar do lambda0
      rhosl01 = ml0s*(efl0*rkfl0)  
      rhosl02 = ml0s*(ml0s**2)*dlog(abs((rkfl0+efl0)/ml0s))
      rhosl0 = (rhosl01 - rhosl02)/2.d0                               !densidade escalar do lambda0

c	densidade escalar do sigma0
      rhoss01 = ms0s*(efs0*rkfs0)  
      rhoss02 = ms0s*(ms0s**2)*dlog(abs((rkfs0+efs0)/ms0s))
      rhoss0 = (rhoss01 - rhoss02)/2.d0                               !densidade escalar do sigma0

c	densidade escalar do xi0
      rhosx01 = mx0s*(efx0*rkfx0)  
      rhosx02 = mx0s*(mx0s**2)*dlog(abs((rkfx0+efx0)/mx0s))
      rhosx0 = (rhosx01 - rhosx02)/2.d0                               !densidade escalar do xi0
   
c	soma das densidade escalar   
      rhosh=rhosl0+rhossm+rhoss0+rhossp+rhosxm+rhosx0

      cahs = xs*rhosh                  ! constante de acoplamento aplicado a densidade escalar dos hiperons

c	soma total das densidades escalar
      rhosb=(rhosn+rhosp+cahs)/pi2 ! total scalar density

c	equação da funcao de movimento nula do boson sigma
      fsigma=gs**2*(rhosb-rb*vsigma**2-rc*vsigma**3)-vsigma       ! sigma equation of motion

c	densidade do neutron up	
      densn=(kfnu**3/3.0d0)/pi2

c	densidade do proton up
      densp=0.0d0
      do i=1,npu
         densp=densp+qe*b*fpu(i)/2.d0/pi2
      enddo

c	densidade do proton down
      do i=1,npd
         densp=densp+qe*b*fpd(i)/2.d0/pi2
      enddo

c	densidade do lambda0
      densl0 =  (rkfl0**3)/(3.0d0*pi2)                            !lambda zero density

c	densidade do sigma0
      denss0 =  (rkfs0**3)/(3.0d0*pi2)                            !sigma zero density

c	densidade do xi0
      densx0 =  (rkfx0**3)/(3.0d0*pi2)                            !xi zero density

c	densidade do sigma-      
      denssm=0.d0
      do i=1,nsm
         denssm=denssm+dg(i)*qe*b*fsm(i)/2.0d0/pi2
      enddo

c	densidade do sigma+
      denssp=0.d0
      do i=1,nsp
         denssp=denssp+dg(i)*qe*b*fsp(i)/2.0d0/pi2
      enddo

c	densidade do xi-
      densxm=0.d0
      do i=1,nxm
         densxm=densxm+dg(i)*qe*b*fxm(i)/2.0d0/pi2
      enddo

c	Soma das densidades
      densh=densl0+denssm+denss0+denssp+densxm+densx0
      cahv = xv*densh                     !constante de acoplamento aplicada a densidade numerica dos hiperons

c	Soma total das dens barionicas	
      densb=densn+densp+densh                                           ! total baryon density

c	Definindo funcoes para alimentar as funcoes 
c	fomega, frho, charge, number
      ddd = densn+densp+cahv

      frhnuc=(densp-densn)/2.0d0
      fhh= -denssm+denssp-(densxm/2.d0)+(densx0/2.d0)
      frhh=xv*fhh
  
c	Funcoes
      fomega=gv**2*(ddd)-vomega               		! omega equation of motion
      frho=gr**2*(frhnuc+frhh)-vrho           		! rho equation of motion
    
      chargeh= -denssm+denssp-densxm
c
      charge =charge+densp+chargeh
      rnumber=rnumber+densb

c	densidades leptonicas
c	eletron
      dense=0.d0
      do i=1,ne
         dense=dense+dg(i)*qe*b*fe(i)/2.0d0/pi2
      enddo
      charge=charge-dense

c	muon
      densmu=0.d0
      do i=1,nu
         densmu=densmu+dg(i)*qe*b*fmu(i)/2.0d0/pi2
      enddo

c	densidade barionica total
      nbt=rnumber

c	funcoes de cargas
      charge=charge-densmu

c	redefinindo as densidades barionicas
      nb(1)=densn                                          	!neutron density
      nb(2)=densp                                          	!proton density
      nb(3)=densl0                                          	!lambda0 density
      nb(4)=denssm                                          	!sigma- density
      nb(5)=denss0                                            !sigma0 density
      nb(6)=denssp                                            !sigma+ density
      nb(7)=densxm                                            !xi - density
      nb(8)=densx0                                            !x0 density

c	redefinindo as densidades leptonicas
      nl(1)=dense                                       ! electron density
      nl(2)=densmu                                      ! muon density

c	Funcoes para broyden
      fvec(1)=fsigma
      fvec(2)=fomega
      fvec(3)=frho
      fvec(4)=charge

      do i=1,4
         fvec1(i)=fvec(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine eos(mue,kfnu,kfnd,vsigma,vomega,vrho,ener,press)
c     this subroutine calculates the eos of this stellar matter
      implicit double precision(a-h,k-z)   
      integer n

      dimension ml(2),mb(8)
      dimension nb(8),nl(2),fpu(9999),fpd(9999),fe(9999)
      dimension fmu(9999), fsm(9999), fsp(9999),fxm(9999)
      common/cdata/pi,pi2,ml,rb,rc,rxi,gs,gv,gr,qe,b,ammp,ammn,mb
      common/cnbt/nbt,nb,nl
      common/out/mns,efn,efp,efe,efmu,mup,mun,fpu,fpd,fe,fmu
     &,npu,npd,ne,nu,rkfl0,efl0, ml0s, mul0
     &,msms,efsm,musm,fsm,nsm,rkfs0,efs0,ms0s,mus0
     &,msps,efsp,musp,fsp,nsp,rkfx0,efx0,mx0s,mux0
     &,mxms,efxm,muxm,fxm,nxm
     common/lvs/csi
c  ============================================
c    Termos de dens de energia

c  energia dos termos dos mesons
      enerf  =0.0d0		

      enerf=(vsigma/gs)**2/2.0d0+(vomega/gv)**2/2.0d0
     &+(vrho/gr)**2/2.0d0+rb*vsigma**3/3.0d0+rc*vsigma**4/4.0d0
     &+rxi*vomega**4/4.0d0

c  energia neutron
      enerbar=0.d0
      re3n=0.d0
      re3n=(efn**3*kfnu/2.0d0
     &-(mns/4.0d0)*(mns*kfnu*efn
     &+(mns)**3*dlog(dabs((kfnu+efn)/mns)))
     &)/2.0d0/pi2
     
      enerbar=enerbar+re3n

c  energia lambda0
      re3l0=0.d0
      re3l0=(efl0**3*rkfl0/2.0d0
     &-(ml0s/4.0d0)*(ml0s*rkfl0*efl0
     &+(ml0s)**3*dlog(dabs((rkfl0+efl0)/ml0s)))
     &)/2.0d0/pi2
     
      enerbar=enerbar+re3l0

c  energia sigma0
      re3s0=0.d0
      re3s0=(efs0**3*rkfs0/2.0d0
     &-(ms0s/4.0d0)*(ms0s*rkfs0*efs0
     &+(ms0s)**3*dlog(dabs((rkfs0+efs0)/ms0s)))
     &)/2.0d0/pi2
    
      enerbar=enerbar+re3s0

c  energia xi0
      re3x0=0.d0
      re3x0=(efx0**3*rkfx0/2.0d0
     &-(mx0s/4.0d0)*(mx0s*rkfx0*efx0
     &+(mx0s)**3*dlog(dabs((rkfx0+efx0)/mx0s)))
     &)/2.0d0/pi2
     
      enerbar=enerbar+re3x0

c  energia proton
c	up
      re3p=0.0d0
      do i=1,npu
      mnsu=dsqrt(mns**2+2*qe*b*(i-1)+csi**2)
      re3p=re3p+qe*b*(efp*fpu(i)
     &+mnsu**2*dlog(dabs((fpu(i)+efp)/mnsu)))/4.0d0/pi2
      enddo

c	down
      do i=1,npd
      mnsd=dsqrt(mns**2+2*qe*b*i+csi**2)+ammp*b
      re3p=re3p+qe*b*(efp*fpd(i)
     &+mnsd**2*dlog(dabs((fpd(i)+efp)/mnsd)))/4.0d0/pi2
      enddo

      enerbar=enerbar+re3p

c  energia sigma-
      re3sm=0.0d0
      do i=1,nsm
         msmsu=dsqrt(msms**2+2*qe*b*(i-1)+csi**2)
         re3sm=re3sm+dg(i)*qe*b*(efsm*fsm(i)
     &  +msmsu**2*dlog(dabs((fsm(i)+efsm)/msmsu)))/4.0d0/pi2
      enddo

      enerbar=enerbar+re3sm

c  energia sigma+
      re3sp=0.0d0
      do i=1,nsp
         mspsu=dsqrt(msps**2+2*qe*b*(i-1)+csi**2)
         re3sp=re3sp+dg(i)*qe*b*(efsp*fsp(i)
     &  +mspsu**2*dlog(dabs((fsp(i)+efsp)/mspsu)))/4.0d0/pi2
      enddo

      enerbar=enerbar+re3sp

c  energia xi-
      re3xm=0.0d0
      do i=1,nxm
         mxmsu=dsqrt(mxms**2+2*qe*b*(i-1)+csi**2)
         re3xm=re3xm+dg(i)*qe*b*(efxm*fxm(i)
     &  +mxmsu**2*dlog(dabs((fxm(i)+efxm)/mxmsu)))/4.0d0/pi2
      enddo

      enerbar=enerbar+re3xm

c  energia eletron
      enerlep=0.d0
      re3le=0.0d0
      do i=1,ne
         ml1=dsqrt(ml(1)**2+2*qe*b*(i-1)+csi**2)
         re3le=re3le+dg(i)*qe*b*(efe*fe(i)
     &   +ml1**2*dlog(dabs((fe(i)+efe)/ml1)))/4.0d0/pi2
      enddo
         enerlep=enerlep+re3le

c  energia muon
      re3lmu=0.0d0
      do i=1,nu
         ml2=dsqrt(ml(2)**2+2*qe*b*(i-1)+csi**2)
         re3lmu=re3lmu+dg(i)*qe*b*(efmu*fmu(i)
     &   +ml2**2*dlog(dabs((fmu(i)+efmu)/ml2)))/4.0d0/pi2
      enddo

c   soma energia leptons
      enerlep=enerlep+re3lmu
      ener=0.0d0
      ener=enerf+enerbar+enerlep
c     &+b**2/8.0d0/pi2

c  ============================================
c    Termos de Pressao
      press=0.0d0
      press1=mun*(nb(1)+nb(3)+nb(5)+nb(8))        	!barions sem carga
      press2=mup*(nb(2)+nb(6))                     	!barions positivos
      press3=mue*(nl(1)+nl(2))               		!leptons
      press4=musm*(nb(4)+nb(7))                     	!barions negativos
      press=(press1+press2+press3+press4)-ener

c      press = press*1.0084d0 ! nomalization factor
      c1 = nb(2)+nb(6)
      c2 = nb(4)+nl(1)+nl(2)+nb(7)
      c3 = mul0-mun
      c4 = mup-(mun-mue)
      c5=  musm-(mun+mue)

      return
      end
c-----------------------------------------------------------------------
      function dg(i)
      implicit double precision(a-h,k-z)  
      if(i.eq.1)dg=1.0d0
      if(i.gt.1)dg=2.0d0
      return
      end
c-----------------------------------------------------------------------      
c=======================================================================
c     other subroutines from numerical recipes
c=======================================================================

      subroutine broydn(x,n,check)
      integer n,nn,np,maxits
      double precision x(n),fvec,eps,tolf,tolmin,tolx,stpmx
      logical check
      parameter (np=40,maxits=40000,eps=1.d-19,tolf=1.d-12,
     &tolmin=1.d-12,tolx=eps,stpmx=100.d0)
      common /newtv/ fvec(np),nn
      save/newtv/
cu    uses fdjac,fmin,lnsrch,qrdcmp,qrupdt,rsolv
      integer i,its,j,k
      double precision den,f,fold,stpmax,sum,temp,test,c(np),d(np)
     &,fvcold(np),g(np),p(np),qt(np,np),r(np,np),s(np),t(np),w(np)
     &,xold(np),fmin
      logical restrt,sing,skip
      external fmin
      nn=n
      f=fmin(x)
      test=0.d0

      do 11 i=1,n
        if(dabs(fvec(i)).gt.test)test=dabs(fvec(i))
11    continue
      if(test.lt..01d0*tolf)then
        check=.false.
        return
      endif
      sum=0.d0
      do 12 i=1,n
        sum=sum+x(i)**2
12    continue
      stpmax=stpmx*max(dsqrt(sum), dfloat(n))
      restrt=.true.
      do 44 its=1,maxits
        if(restrt)then
          call fdjac(n,x,fvec,np,r)
          call qrdcmp(r,n,np,c,d,sing)
          if(sing) stop
          do 14 i=1,n
            do 13 j=1,n
              qt(i,j)=0.d0
13          continue
            qt(i,i)=1.d0
14        continue
          do 18 k=1,n-1
            if(c(k).ne.0.d0)then
              do 17 j=1,n
                sum=0.d0
                do 15 i=k,n
                  sum=sum+r(i,k)*qt(i,j)
15              continue
                sum=sum/c(k)
                do 16 i=k,n
                  qt(i,j)=qt(i,j)-sum*r(i,k)
16              continue
17            continue
            endif
18        continue
          do 21 i=1,n
            r(i,i)=d(i)
            do 19 j=1,i-1
              r(i,j)=0.d0
19          continue
21        continue
        else
          do 22 i=1,n
            s(i)=x(i)-xold(i)
22        continue
          do 24 i=1,n
            sum=0.d0
            do 23 j=i,n
              sum=sum+r(i,j)*s(j)
23          continue
            t(i)=sum
24        continue
          skip=.true.
          do 26 i=1,n
            sum=0.d0
            do 25 j=1,n
              sum=sum+qt(j,i)*t(j)
25          continue
            w(i)=fvec(i)-fvcold(i)-sum
            if(dabs(w(i)).ge.eps*(dabs(fvec(i))+dabs(fvcold(i))))then
              skip=.false.
            else
              w(i)=0.d0
            endif
26        continue
          if(.not.skip)then
            do 28 i=1,n
              sum=0.d0
              do 27 j=1,n
                sum=sum+qt(i,j)*w(j)
27            continue
              t(i)=sum
28          continue
            den=0.d0
            do 29 i=1,n
              den=den+s(i)**2
29          continue
            do 31 i=1,n
              s(i)=s(i)/den
31          continue
            call qrupdt(r,qt,n,np,t,s)
            do 32 i=1,n
              if(r(i,i).eq.0.d0) stop
              d(i)=r(i,i)
32          continue
          endif
        endif
        do 34 i=1,n
          sum=0.d0
          do 33 j=1,n
            sum=sum+qt(i,j)*fvec(j)
33        continue
          g(i)=sum
34      continue
        do 36 i=n,1,-1
          sum=0.d0
          do 35 j=1,i
            sum=sum+r(j,i)*g(j)
35        continue
          g(i)=sum
36      continue
        do 37 i=1,n
          xold(i)=x(i)
          fvcold(i)=fvec(i)
37      continue
        fold=f
        do 39 i=1,n
          sum=0.d0
          do 38 j=1,n
            sum=sum+qt(i,j)*fvec(j)
38        continue
          p(i)=-sum
39      continue
        call rsolv(r,n,np,d,p)
        call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin)
        test=0.d0
        do 41 i=1,n
          if(dabs(fvec(i)).gt.test)test=dabs(fvec(i))
41      continue
        if(test.lt.tolf)then
          check=.false.
          return
        endif
        if(check)then
          if(restrt)then
            return
          else
            test=0.d0
            den=max(f,.5d0*n)
            do 42 i=1,n
              temp=dabs(g(i))*max(dabs(x(i)),1.d0)/den
              if(temp.gt.test)test=temp
42          continue
            if(test.lt.tolmin)then
              return
            else
              restrt=.true.
            endif
          endif
        else
          restrt=.false.
          test=0.d0
          do 43 i=1,n
            temp=(dabs(x(i)-xold(i)))/max(dabs(x(i)),1.d0)
            if(temp.gt.test)test=temp
43        continue
          if(test.lt.tolx)return
        endif
44    continue
c      pause 'maxits exceeded in broydn'
      write(*,*)'maxits exceeded in broydn'
      stop
      end
c-----------------------------------------------------------------------
      subroutine fdjac(n,x,fvec,np,df)
      integer n,np,nmax
      double precision df(np,np),fvec(n),x(n),eps
      parameter (nmax=40,eps=1.d-8)
cu    uses funcv
      integer i,j
      double precision h,temp,f(nmax)
      do 12 j=1,n
        temp=x(j)
        h=eps*dabs(temp)
        if(h.eq.0.d0)h=eps
        x(j)=temp+h
        h=x(j)-temp
        call funcv(n,x,f)
        x(j)=temp
        do 11 i=1,n
          df(i,j)=(f(i)-fvec(i))/h
11      continue
12    continue
      return
      end
c-----------------------------------------------------------------------
      function fmin(x)
      integer n,np
      double precision fmin,x(*),fvec
      parameter (np=40)
      common /newtv/ fvec(np),n
      save /newtv/
cu    uses funcv
      integer i
      double precision sum
      call funcv(n,x,fvec)
      sum=0.d0
      do 11 i=1,n
      sum=sum+fvec(i)**2
11    continue
      fmin=0.5d0*sum
      return
      end
c-----------------------------------------------------------------------
      subroutine lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)
      integer n
      logical check
      double precision f,fold,stpmax,g(n),p(n),x(n),xold(n),func,alf
     &,tolx
      parameter (alf=1.d-4,tolx=1.d-16)
      external func
cu    uses func
      integer i
      double precision a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2
     &,slope,sum,temp,test,tmplam

      check=.false.
      sum=0.d0
      do 11 i=1,n
        sum=sum+p(i)*p(i)
11    continue
      sum=dsqrt(sum)
      if(sum.gt.stpmax)then
        do 12 i=1,n
          p(i)=p(i)*stpmax/sum
12      continue
      endif
      slope=0.d0
      do 13 i=1,n
        slope=slope+g(i)*p(i)
13    continue
      test=0.d0
      do 14 i=1,n
        temp=dabs(p(i))/max(dabs(xold(i)),1.d0)
        if(temp.gt.test)test=temp
14    continue
      alamin=tolx/test
      alam=1.d0
1     continue
        do 15 i=1,n
          x(i)=xold(i)+alam*p(i)
15      continue
        f=func(x)
        if(alam.lt.alamin)then
          do 16 i=1,n
            x(i)=xold(i)
16        continue
          check=.true.
          return
        else if(f.le.fold+alf*alam*slope)then
          return
        else
          if(alam.eq.1.d0)then
            tmplam=-slope/(2.d0*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold2-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.d0)then
              tmplam=-slope/(2.d0*b)
            else
              disc=b*b-3.d0*a*slope
              if(disc.lt.0.d0) stop 'roundoff problem in lnsrch'
              tmplam=(-b+dsqrt(disc))/(3.d0*a)
            endif
            if(tmplam.gt..5d0*alam)tmplam=.5d0*alam
          endif
        endif
        alam2=alam
        f2=f
        fold2=fold
        alam=max(tmplam,.1d0*alam)
      goto 1
      end
c-----------------------------------------------------------------------
      subroutine qrdcmp(a,n,np,c,d,sing)
      integer n,np
      double precision a(np,np),c(n),d(n)
      logical sing
      integer i,j,k
      double precision scale,sigma,sum,tau
      sing=.false.
      do 17 k=1,n-1
        scale=0.d0
        do 11 i=k,n
          scale=max(scale,dabs(a(i,k)))
11      continue
        if(scale.eq.0.d0)then
          sing=.true.
          c(k)=0.d0
          d(k)=0.d0
        else
          do 12 i=k,n
            a(i,k)=a(i,k)/scale
12        continue
          sum=0.d0
          do 13 i=k,n
            sum=sum+a(i,k)**2
13        continue
          sigma=sign(dsqrt(sum),a(k,k))
          a(k,k)=a(k,k)+sigma
          c(k)=sigma*a(k,k)
          d(k)=-scale*sigma
          do 16 j=k+1,n
            sum=0.d0
            do 14 i=k,n
              sum=sum+a(i,k)*a(i,j)
14          continue
            tau=sum/c(k)
            do 15 i=k,n
              a(i,j)=a(i,j)-tau*a(i,k)
15          continue
16        continue
        endif
17    continue
      d(n)=a(n,n)
      if(d(n).eq.0.d0)sing=.true.
      return
      end
c-----------------------------------------------------------------------
      subroutine qrupdt(r,qt,n,np,u,v)
      integer n,np
      double precision r(np,np),qt(np,np),u(np),v(np)
cu    uses rotate
      integer i,j,k
      do 11 k=n,1,-1
        if(u(k).ne.0.d0)goto 1
11    continue
      k=1
1     do 12 i=k-1,1,-1
        call rotate(r,qt,n,np,i,u(i),-u(i+1))
        if(u(i).eq.0.d0)then
          u(i)=dabs(u(i+1))
        else if(dabs(u(i)).gt.dabs(u(i+1)))then
          u(i)=dabs(u(i))*dsqrt(1.d0+(u(i+1)/u(i))**2)
        else
          u(i)=dabs(u(i+1))*dsqrt(1.d0+(u(i)/u(i+1))**2)
        endif
12    continue
      do 13 j=1,n
        r(1,j)=r(1,j)+u(1)*v(j)
13    continue
      do 14 i=1,k-1
        call rotate(r,qt,n,np,i,r(i,i),-r(i+1,i))
14    continue
      return
      end
c-----------------------------------------------------------------------
      subroutine rsolv(a,n,np,d,b)
      integer n,np
      double precision a(np,np),b(n),d(n)
      integer i,j
      double precision sum
      b(n)=b(n)/d(n)
      do 12 i=n-1,1,-1
        sum=0.d0
        do 11 j=i+1,n
          sum=sum+a(i,j)*b(j)
11      continue
        b(i)=(b(i)-sum)/d(i)
12    continue
      return
      end
c-----------------------------------------------------------------------
      subroutine rotate(r,qt,n,np,i,a,b)
      integer n,np,i
      double precision a,b,r(np,np),qt(np,np)
      integer j
      double precision c,fact,s,w,y
      if(a.eq.0.d0)then
        c=0.d0
        s=sign(1.d0,b)
      else if(dabs(a).gt.dabs(b))then
        fact=b/a
        c=sign(1.d0/dsqrt(1.d0+fact**2),a)
        s=fact*c
      else
        fact=a/b
        s=sign(1.d0/dsqrt(1.d0+fact**2),b)
        c=fact*s
      endif
      do 11 j=i,n
        y=r(i,j)
        w=r(i+1,j)
        r(i,j)=c*y-s*w
        r(i+1,j)=s*y+c*w
11    continue
      do 12 j=1,n
        y=qt(i,j)
        w=qt(i+1,j)
        qt(i,j)=c*y-s*w
        qt(i+1,j)=s*y+c*w
12    continue
      return
      end
c========================================================================
