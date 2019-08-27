      Program CORAL1DC  
C    variable grid 1-D gasdynamics
C     flow with photoionization from source at origin
C=======================================================================
C> @file coral1d.f
C     > @author Ary Rodriguez, Gazyna Stasinska, Alex Raga & Zakaria Meliani
C     (photoionization G. Mellema)
C     > @date 3/Abr/2019 Version 1.0     
C     Copyright (c) 2019 Ary Rodriguez & A. Raga
C     Coral-1D is free software; you can redistribute it and/or modify
C     it under the terms of the GNU General Public License as published by
C     the Free Software Foundation; either version 3 of the License, or
C     (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C     GNU General Public License for more details.
C     You should have received a copy of the GNU General Public License
C     along with this program.  If not, see http://www.gnu.org/licenses/
C=======================================================================
C     > @brief Coral-1D Main Program
C     > @n The code itegrates Euler equations in one dimensions including the
C     source terms of chemistry solved in Mellema (1997),
C     @n The flow (conserved) variables are taken to be:
C
C     Variables in code
c      prim(1) pressure  [erg/cm^3] 
c      prim(2) velocity  [cm/s ]
c      prim(3) mass density gas  [gr/cm^3] 
c      prim(4) numerical density gas  [1/cm^3] 
c      prim(5:NEQ) numerical density of each specie [1/cm^3]
c      xh, xhe, xc …   fraccion de especie por numero
c      denelec     [1/cm^3]
c      AbH_H ...  [es numeric abundance of H with respect to H 
c      ABHe-H es numeric abundance of He with respect to H  etc...
c      AMOL   
c      AMASS    [gr] proton mass
c      AMU         =   AMOL* AMASS    (masa molecular media) [gr]
c      den             [1/cm^3] cual es la definition exacta? density de que?
c      rho             [gr/cm^3]  cual es  la definition exacta?
c      XJ       fraccion de especie por numero
C
C     1 : energy, 
C     2 : rho u,
C     3 : rho
C     4 : n
C     5 : n(HI)
C     6 : n(HII)
C     7 : n(HeI)
C     8 : n(HeII)
C     9 : n(CIII)
C    10 : n(CIV)
C    11 : n(CV)
C    12 : n(OII)
C    13 : n(OIII)
C    14 : n(OIV)
C    15 : n(OV)
C    16 : n(OVI)
C    17 : n(NeII)
C    18 : n(NeIII)
C    19 : n(NeIV)
C    20 : n(NeV)
C    21 : n(NeVI)
C    22 : n(SIII)
C    23 : n(SIV)
C    24 : n(SV)
C    25 : n(SVI)
C    26 : n(NII)
C    27 : n(NIII)
C    28 : n(NIV)
C    29 : n(NV)
C    30 : n(NVI)
C
      implicit real*8 (a-h,o-z)
      PARAMETER(NR=15000,NEQ1=5,NEQ=30)
      PARAMETER(YR=3.15E7)
      REAL*8 LSTAR
      DIMENSION U(NR,NEQ),UP(NR,NEQ),RM(NR),UU(NEQ),PRIM(NEQ)
      DIMENSION COLHEAT(2,NR)
      COMMON/IT/IT
      COMMON/RADIUS/DR,RMIN
      COMMON/GAMMA/GAMMA
      COMMON/COLHEATF/COLHEAT
      common/starprop/Teff,alst,qh      
      save
C
      call cpu_time(tini)
      TIME=0.
      CALL INITFLOW(U,TMAX,dtprint,tprint,itprint)
C
      call output(U,dtprint,itprint)
C
      DO WHILE (TIME.LE.TMAX)
C
         WRITE(*,*) 'TIME='
     &        ,TIME/YR,'[YR]  '
C
         call TIMESTEP(U,DT)
C
         CALL TSTEP(U,DT,TIME)
C
         if(time.ge.tprint) then                                                
            write(*,*) 'output------->',itprint                             
            write(*,*) '----> Tiempo, dt, tmax', time, dt, tmax
            call output(U,dtprint,itprint)
            tprint=tprint+dtprint                                           
            itprint=itprint+1                                               
         end if
C         
         TIME=TIME+DT
C

      END DO
      call cpu_time(tfin)
      write(*,*)'Simulation Time=',(tfin-tini)/60,' min'

C     
      STOP
      END
C
C     
      SUBROUTINE INITFLOW(U,TMAX,dtprint,tprint,itprint)
C
C    reads initial parameters, etc.
C
      implicit real*8 (a-h,o-z)
      PARAMETER(NR=15000,NEQ1=5,NEQ=30)
      PARAMETER(AMASS=1.674D-24,AK=1.38D-16,PC=3.08E18)
      PARAMETER(AMSUN=1.989E33,YR=3.15E7,pi=3.141592)
      CHARACTER*40 FILE1,FILE2
c
c     Note: abundances also defined in inirad !!!
c
      PARAMETER(Nmax=1000)
      DIMENSION U(NR,NEQ)
      DIMENSION RM(NR),XJ(NEQ)
      DIMENSION anum(Nmax),tmod(Nmax),alogl(Nmax),tstar(Nmax)
     &     ,amstar(Nmax),alogdmdt(Nmax),tempmod(Nmax),vterm(Nmax)
     &     ,alstars(Nmax),dmdt(Nmax),qhmod(Nmax),qhemod(Nmax)
     &     ,tmod1(Nmax),tmod2(Nmax),tmod3(Nmax)
      COMMON/Modelo/tmod,alogl,tstar,amstar
     &     ,alogdmdt,tempmod,vterm,qhmod,qhemod
     &     ,tmod1,tmod2,tmod3
      COMMON/nmaxx/nmax2,nmax3,nmax1
      character*20 modelfile
      character*1 modeloption
      COMMON/PAR/XJ,DEN
      COMMON/RADIUS/DR,RMIN
      common/hhphot/hhphot
      Common/amoln/amol
      COMMON/GAMMA/GAMMA
      COMMON/HDSTUFF/Co,eta
      COMMON/ABUND/XH,XHE,XC,XN,XO,XS,XAr,XNe
      common/modelf/modelfile
      integer date_time(8)
      character*10 b(3)
      common/celd/NRW
      common/wohd/DEN0,T0,V0
      common/nummethod/imet,ichem
      common/ifmodel/ifmod
      common/soundspeed0/cspeed0,twindagb,rho0
      common/percent/perc_fot,perc_temp
      common/iterat/iterat_limit,iflqh,ifltemp
C
      call date_and_time(b(1), b(2), b(3), date_time)
      print *,'Day :  ',b(1)(7:8),'/',b(1)(5:6),'/',b(1)(1:4)
      print *,'Time : ',b(2)(1:2),'h ',b(2)(3:4),'m ',b(2)(5:6),'s'
C
      GAMMA=5./3.
C
      OPEN(UNIT=9,FILE='Coral-1D.in',STATUS='UNKNOWN')
C
      read(9,'(A20)') modelfile
      READ(9,*)
      READ(9,*)
      read(9,'(A20)') modeloption
C     
 70   CONTINUE
      READ(9,*,ERR=70) Co
 71   CONTINUE
      READ(9,*,ERR=71) imet
 72   CONTINUE
      READ(9,*,ERR=72) ichem
 73   CONTINUE
      READ(9,*,ERR=73) ETA
 74   CONTINUE
      READ(9,*,ERR=74) perc_fot
 75   CONTINUE
      READ(9,*,ERR=75) perc_temp
 76   CONTINUE
      READ(9,*,ERR=76) iterat_limit
 77   CONTINUE
      READ(9,*,ERR=77) RMAX, TMAX, dtprint
 78   CONTINUE
      READ(9,*,ERR=78) amagb1,amagb2,amagb3
 79   CONTINUE
      READ(9,*,ERR=79) tagb1,tagb2,tagb3
 80   CONTINUE
      READ(9,*,ERR=80) vwindagb,twindagb
 81   CONTINUE
      READ(9,*,ERR=81) AbH_H,AbHE_H,AbC_H, AbN_H, AbO_H, AbS_H, AbAr_H
     & , AbNe_H
C
      close(9)
      if (modeloption(1:1) .eq. 'M') ifmod=0
      if (modeloption(1:1) .eq. 'S') ifmod=1
      if (modeloption(1:1) .eq. 'P') ifmod=2
      print*,modeloption(1:1),ifmod
C     
      AMOL=(1+4.*AbHE_H+12.*AbC_H+14*AbN_H+16.*AbO_H+32.*AbS_H
     &     +20.*AbNe_H+40.*AbAr_H)
      AMU=AMOL*AMASS
C
      Ab_Total=AbH_H+AbHe_H+AbC_H+AbN_H+AbO_H+AbS_H+AbAr_H+Ab_Ne
      XH =AbH_H/Ab_Total
      XHe=AbHe_H/Ab_Total
      XC =AbC_H/Ab_Total
      XN =AbN_H/Ab_Total
      XO =AbO_H/Ab_Total
      XS =AbS_H/Ab_Total
      XAr=AbAr_H/Ab_Total
      XNe=AbNe_H/Ab_Total
C
      perc_fot=perc_fot/100.
      perc_temp=perc_temp/100.
C
      RMIN=5e14!RMAX/float(NR)
C     Hydrogen I,II
C
c
      XJ(5)=XH*1.
      XJ(6)=XH*1.E-5
C
C     Helium I,II
C
      XJ(7)=XHE*1.
      XJ(8)=XHE*1.E-5
C
C     Carbon III,IV,V
C
      XJ(9)=XC*1.E-4
      XJ(10)=XC*1.E-4
      XJ(11)=XC*1.E-4
C
C     Oxygen II,III,IV,V,VI
C
      XJ(12)=XO*1.E-4
      XJ(13)=XO*1.E-4
      XJ(14)=XO*1.E-4
      XJ(15)=XO*1.E-4
      XJ(16)=XO*1.E-4
C
C     Neon II,III,IV,V,VI
C
      XJ(17)=XNE*1.E-4
      XJ(18)=XNE*1.E-4
      XJ(19)=XNE*1.E-4
      XJ(20)=XNE*1.E-4
      XJ(21)=XNE*1.E-4
C
C     Sulphur III,IV,V,VI
C
      XJ(22)=XS*1.E-4
      XJ(23)=XS*1.E-4
      XJ(24)=XS*1.E-4
      XJ(25)=XS*1.E-4
C
C     Nitrogen II,III,IV,V,VI
C
      XJ(26)=XO*1.E-4
      XJ(27)=XO*1.E-4
      XJ(28)=XO*1.E-4
      XJ(29)=XO*1.E-4
      XJ(30)=XO*1.E-4
C
C
      vwindagb=vwindagb*1.e5
      TMAX=TMAX*YR
      dtprint=dtprint*YR
      RMAX=RMAX
      DR=(RMAX)/FLOAT(NR-1)
C
      amagb1=amagb1*AMSUN/YR
      amagb2=amagb2*AMSUN/YR
      amagb3=amagb3*AMSUN/YR
      tagb1=tagb1*YR+1.
      tagb2=tagb2*YR
      tagb3=tagb3*YR

      if (ifmod .eq. 0) call inputm
      if (ifmod .eq. 1) call inputp
C      
      CALL INIRAD
C     sound speed in the interstellar medium
      cspeed0=dsqrt(ak*twindagb/AMASS/AMOL)
      if (ichem .eq. 1) Then
      DEN0=1.e4
      T0=1.e4
      V0=0.
      endif

C
      NRW=2
C
      slM1=(dlog10(amagb2)-dlog10(amagb1))/(dlog10(tagb2)-dlog10(tagb1))
      slM2=(dlog10(amagb3)-dlog10(amagb2))
     &     /(dlog10(tagb3)-dlog10(tagb2))
C     
      DO I=1,NR
         R=(FLOAT(I)-1.)*DR+RMIN
         RM(I)=1./R
         tdyn=R/vwindagb
         tm=tagb3-tdyn
         if(tm .lt. tagb1) amdotagb=amagb1
         if(tm .gt. tagb1 .and. tm .lt. tagb2) then
            amdotagb=10**(slM1*(dlog10(tm)-dlog10(tagb1))
     &           +dlog10(amagb1))
         elseif(tm .gt. time2) then
            amdotagb=10**(slM2*(dlog10(tm)-dlog10(tagb2))
     &           +dlog10(amagb2))
         endif
            rho=amdotagb/(4.*pi*vwindagb*(R)**2)
C
            DEN=rho/AMASS/AMOL
            if (ichem .eq. 1) then
               DEN=DEN0
               RHO=DEN*AMASS*AMOL
               Twindagb=T0
               vwindagb=V0
            endif
C
             U(I,3)=rho
             U(I,4)=DEN
            
            U(I,2)=rho*vwindagb
            U(I,1)=twindagb*den*ak/(GAMMA-1.)+0.5*rho*vwindagb**2
            DO IEQ=NEQ1,NEQ
               U(I,IEQ)=DEN*XJ(IEQ)
            END DO
         END DO
C     C
      WRITE(FILE1,1000)modelfile(1:12)                               
C
      open(unit=30,file=FILE1)
      write(30,'(A38,A2,A1,A2,A1,A2)')'Day                               
     &  :  ',b(1)(7:8),'/',b(1)(5:6),'/',b(1)(1:4)
      write(30,'(A38,A2,A1,A2,A1,A2,A1)')'Time                  
     &              :  ',b(2)(1:2),'h',b(2)(3:4),'m ',b(2)(5:6),'s'
      write(30,'(A38,A12)')'********   Stellar Model  *********:',
     &     modelfile
      write(30,'(A38,F7.1)')'Slow Wind Vwind [km/s]              :',
     &     vwindagb/1.e5
      write(30,'(A38,F8.1)')'Slow Wind Tempwind [K]              :',
     &     twindagb
      write(30,'(A38,ES13.4)')'Minimal Mdot [Msun/yr]              :',
     &     amagb1*YR/AMSUN
      write(30,'(A38,ES13.4)')'Interm  Mdot [Msun/yr]              :',
     &     amagb2*YR/AMSUN
      write(30,'(A38,ES13.4)')'Maximal Mdot [Msun/yr]              :',
     &     amagb3*YR/AMSUN
      write(30,'(A38,ES13.4)')'Initial Time [yr]                   :',
     &     (tagb1-1.)/YR
      write(30,'(A38,ES13.4)')'Interm  Time [yr]                   :',
     &     tagb2/YR
      write(30,'(A38,ES13.4)')'Final   Time [yr]                   :',
     &     tagb3/YR
      write(30,'(A38)')' ************************************'
      write(30,'(A38,ES13.4)')'AbH_H                               :',
     &     AbH_H
      write(30,'(A38,ES13.4)')'AbHe_H                              :',
     &     AbHe_H
      write(30,'(A38,ES13.4)')'AbC_H                               :',
     &     AbC_H
      write(30,'(A38,ES13.4)')'AbN_H                               :',
     &     AbN_H
      write(30,'(A38,ES13.4)')'AbO_H                               :',
     &     AbO_H
      write(30,'(A38,ES13.4)')'AbS_H                               :',
     &     AbS_H
      write(30,'(A38,ES13.4)')'AbAr_H                              :',
     &     AbAr_H
      write(30,'(A38,ES13.4)')'AbNe_H                              :',
     &     AbNe_H
      write(30,'(A38)')'----------Simulation stuff- --------'
      if (imet .eq. 0) then
         write(30,'(A38)')'Solver:        2o. order Godunov'
      endif
      if (imet .eq. 1) then
         write(30,'(A38)')'Solver:        1er order Godunov'
      endif
      if (imet .eq. 2) then
         write(30,'(A38)')'Solver:         Predictor-Corrector'
      endif
      if (imet .eq. 2) then
         write(30,'(A38)')'Solver:         Fred-Lax'
      endif
      write(30,'(A38,ES13.4)')'Maximum Radius [cm]                 :',
     &     RMAX
      write(30,'(A38,ES13.4)')'Maximum Time   [yr]                 :',
     &     TMAX/3.15e7
      write(30,'(A38,ES13.4)')'Output Time    [yr]                 :',
     &     dtprint/3.15e7
      write(30,'(A38,I7)')'Number of cells                     :',NR
      write(30,'(A38,I6)')'Number of cells in the wind         :',NRW
      write(30,'(A38,F6.1)')'Courant Number                      :',Co
      if (imet .ge. 2) then
         write(30,'(A38,F4.2)')'Artificial Viscosity          :',eta
      endif
      write(30,'(A38,F6.1)')' % flux phot for call photoion      :',
     &     perc_fot*100
      write(30,'(A38,F6.1)')' % flux temp for -> photoion/atomic :',
     &     perc_temp*100
      write(30,'(A38,I6)')' Max iter before -> photoion        :',
     &     iterat_limit 
      write(30,'(A38,I6)')'Number of variables                 :',NEQ
      write(30,'(A34)')' 1 : energy                             '
      write(30,'(A34)')' 2 : rho u                              ' 
      write(30,'(A34)')' 3 : n_gas                              '
      write(30,'(A34)')' 4 : passive scalar                     '
      write(30,'(A34)')' 5 : n(HI)                              '
      write(30,'(A34)')' 6 : n(HII)                             '
      write(30,'(A34)')' 7 : n(HeI)                             '
      write(30,'(A34)')' 8 : n(HeII)                            '
      write(30,'(A34)')' 9 : n(CIII)                            '
      write(30,'(A34)')' 10 : n(CIV)                            '
      write(30,'(A34)')' 11 : n(CV)                             '
      write(30,'(A34)')' 12 : n(OII)                            '
      write(30,'(A34)')' 13 : n(OIII)                           '
      write(30,'(A34)')' 14 : n(OIV)                            '
      write(30,'(A34)')' 15 : n(OV)                             '
      write(30,'(A34)')' 16 : n(OVI)                            '
      write(30,'(A34)')' 17 : n(NeII)                           '
      write(30,'(A34)')' 18 : n(NeIII)                          '
      write(30,'(A34)')' 19 : n(NeIV)                           '
      write(30,'(A34)')' 20 : n(NeV)                            '
      write(30,'(A34)')' 21 : n(NeVI)                           '
      write(30,'(A34)')' 22 : n(SIII)                           '
      write(30,'(A34)')' 23 : n(SIV)                            '
      write(30,'(A34)')' 24 : n(SV)                             '
      write(30,'(A34)')' 25 : n(SVI)                            '
      write(30,'(A34)')' 26 : n(NII)                            '
      write(30,'(A34)')' 27 : n(NIII)                           '
      write(30,'(A34)')' 28 : n(NIV)                            '
      write(30,'(A34)')' 29 : n(NV)                             '
      write(30,'(A34)')' 30 : n(NVI)                            '
      close(30)
C
      call uflow(U,TIME,DT,ALST,TEFF)
      time=0.                                                                       
      tprint=0.                                                                     
      itprint=0
 1000 FORMAT('./',A,'/used_parameters.txt')
      RETURN
      END
C
      subroutine inputp
      implicit real*8 (a-h,o-z)
      PARAMETER(Nmax=1000)
      PARAMETER(AMASS=1.674D-24,AK=1.38D-16,PC=3.08E18)
      PARAMETER(AMSUN=1.989E33,YR=3.15E7,pi=3.141592)
      DIMENSION anum(Nmax),tmod(Nmax),alogl(Nmax),tstar(Nmax)
     &     ,amstar(Nmax),alogdmdt(Nmax),tempmod(Nmax),vterm(Nmax)
     &     ,alstars(Nmax),dmdt(Nmax),qhmod(Nmax),qhemod(Nmax)
     &     ,tmod1(Nmax),tmod2(Nmax),tmod3(Nmax)
      COMMON/Modelo/tmod,alogl,tstar,amstar
     &     ,alogdmdt,tempmod,vterm,qhmod,qhemod
     &     ,tmod1,tmod2,tmod3
      COMMON/nmaxx/nmax2,nmax3,nmax1
C
      open(1,file='sh0605_sb.pagb',STATUS='UNKNOWN')
      open(2,file='perinotto_dmdt.txt',STATUS='UNKNOWN')
      open(3,file='perinotto_v.txt',STATUS='UNKNOWN')
C
      Nmax1=31
      Nmax2=31
      Nmax3=25
C
      Do k=1,Nmax1
C
         READ(1,*,ERR=101)tm,tst,alo
         tmod1(k)=tm*YR
         alogl(k)=alo
         tstar(k)=tst
C
      End do
C
      close(1)
C
 101  CONTINUE
      Do k=1,Nmax2 
         READ(2,*,ERR=103)tm,alogd
         tmod2(k)=tm*YR
         alogdmdt(k)=alogd
      End do
 103  CONTINUE
c
      close(2)
      Do k=1,Nmax3
         READ(3,*,ERR=105)tm,vtermw
         tmod3(k)=tm*YR
         vterm(k)=vtermw
      End do
 105  CONTINUE

      close(3)
C
      Return
      End
C
      SUBROUTINE INPUTM
      implicit real*8 (a-h,o-z)
      PARAMETER(Nmax=1000)
      PARAMETER(AMASS=1.674D-24,AK=1.38D-16,PC=3.08E18)
      PARAMETER(AMSUN=1.989E33,YR=3.15E7,pi=3.141592)
      DIMENSION anum(Nmax),tmod(Nmax),alogl(Nmax),tstar(Nmax)
     &     ,amstar(Nmax),alogdmdt(Nmax),tempmod(Nmax),vterm(Nmax)
     &     ,alstars(Nmax),dmdt(Nmax),qhmod(Nmax),qhemod(Nmax)
     &     ,tmod1(Nmax),tmod2(Nmax),tmod3(Nmax)
      COMMON/Modelo/tmod,alogl,tstar,amstar
     &     ,alogdmdt,tempmod,vterm,qhmod,qhemod
     &     ,tmod1,tmod2,tmod3
      COMMON/nmaxx/nmax2,nmax3,nmax1
      common/modelf/modelfile
            character*20 modelfile
      open(1,file=modelfile,STATUS='UNKNOWN')
      kk=1
      Nmax2=0
C
      iflag=1
C
 100     CONTINUE 
         READ(1,*,ERR=100)tmod(1),alogl(1),tstar(1),amstar(1)
     &        ,alogdmdt(1),tempmod(1),vterm(1),qhmod(1),qhemod(1)

      Do k=1,Nmax
C
         READ(1,*,ERR=101) tm,alo,tst,ams,alogd,tempmd,vtermw,qh,qhe
c
         if (iflag .eq. 0) then
            tmod(kk)=tm
            alogl(kk)=alo
            tstar(kk)=tst
            amstar(kk)=ams
            alogdmdt(kk)=alogd
            tempmod(kk)=tempmd
            vterm(kk)=vtermw
            qhmod(kk)=qh
            qhemod(kk)=qhe
C
            kk=kk+1
         endif
C
         if (iflag .eq. 1 .and. tm .ge. 0) then
            tmod(kk)=tm
            alogl(kk)=alo
            tstar(kk)=tst
            amstar(kk)=ams
            alogdmdt(kk)=alogd
            tempmod(kk)=tempmd
            vterm(kk)=vtermw
            qhmod(kk)=qh
            qhemod(kk)=qhe
C
            kk=kk+1
         endif
         
         if (iflag .gt. 1) write(*,*)'The table does not contain
     & valid times'
C
      End do
 101  CONTINUE

      Nmax2=kk

      close(1)
      return
      end

C      
      SUBROUTINE UFLOW(U,TIME,DT,ALST,TEFF)
      implicit real*8 (a-h,o-z)
      PARAMETER(NR=15000,NEQ1=5,NEQ=30)
      PARAMETER(Nmax=1000)
      DIMENSION U(NR,NEQ),XJ(NEQ)
      PARAMETER(AMASS=1.674D-24,AK=1.38D-16,PC=3.08E18
     &     ,pi=3.141592)
      DIMENSION anum(Nmax),tmod(Nmax),alogl(Nmax),tstar(Nmax)
     &     ,amstar(Nmax),alogdmdt(Nmax),tempmod(Nmax),vterm(Nmax)
     &     ,alstars(Nmax),dmdt(Nmax),qhmod(Nmax),qhemod(Nmax)
     &     ,tmod1(Nmax),tmod2(Nmax),tmod3(Nmax)
      COMMON/Modelo/tmod,alogl,tstar,amstar
     &     ,alogdmdt,tempmod,vterm,qhmod,qhemod
     &     ,tmod1,tmod2,tmod3
      COMMON/nmaxx/nmax2,nmax3,nmax1
      real*8 ALST
      common/nummethod/imet,ichem
      COMMON/RADIUS/DR,RMIN
      COMMON/GAMMA/GAMMA
      Common/amoln/AMOL
      COMMON/ABUND/XH,XHE,XC,XN,XO,XS,XAr,XNe
      common/celd/NRW
      common/wohd/DEN0,T0,V0
      common/ifmodel/ifmod
      common/soundspeed0/cspeed0,twindagb,rho0
      common/percent/perc_fot,perc_temp
      common/iterat/iterat_limit,iflqh,ifltemp
      common/qhbb/qhbb,qhaprox  
      save it
      save qh00
C
C
      if (ifmod .eq. 0) then
         tmod(1:Nmax2)=tmod(1:Nmax2)-tmod(1)
         qh=ainterpol(time,tmod,qhmod,Nmax2)
         timep=TIME-DT
         teff=ainterpol(time,tmod,tempmod,Nmax2)
C     
         alst=ainterpol(time,tmod,alogl,Nmax2)
         amdot=ainterpol(time,tmod,alogdmdt,Nmax2)
         twind=ainterpol(time,tmod,tempmod,Nmax2)   
C     
         vwind=ainterpol(time,tmod,vterm,Nmax2)
         alst=10**(alst)*3.89e33
         amdot=10**(amdot)*1.98e33/3.15e7
!
      elseif (ifmod .eq. 1) then
C
         alst=ainterpol(time,tmod1,alogl,Nmax1)
         amdot=ainterpol(time,tmod2,alogdmdt,Nmax2)
         teff=ainterpol(time,tmod1,tstar,Nmax1)
         vwind=ainterpol(time,tmod3,vterm,Nmax3)
         if(time .lt. tmod1(1))then
            teff=3.8
            alst=3.8
         endif
         alst=10**(alst)*3.89e33
         amdot=amdot*1.98e33/3.15e7
         vwind=vwind*1.e5
         teff=10**(teff)
         twind=max(5000.,teff)
      endif
!
      if (ichem .eq. 1) then
         TEFF=50000.
         ALST=1.e38
      endif
!     Number of ionizing photons
      qhbb=aionizphotons(TEFF,alst)
!     h nu/kT
      anuh=3.28d15
      x0=ah*anuh/(AK*TEFF)
! Flux = sigma*Teff
      totflux=5.670d-05*teff**4
!
      qhaprox=(2.*pi*alst)/(totflux*clight**2)*(AK*TEFF/ah)**3.
      qhaprox=qhaprox*(x0**2+x0+2.)*exp(-x0)
!
      if (ifmod .eq. 1) qh=qhbb
      if (time .lt. 1.d-100) then
         qh0=0.
         teff0=0.
         it=0
      else
         if(ifmod .eq. 0) then
            qh0=ainterpol(timep,tmod,qhmod,Nmax2)
            teff0=ainterpol(timep,tmod,tempmod,Nmax2)
         else
            qh0=qh00
         endif
      endif

      if (abs(qh-qh0)/qh .gt. perc_fot.or. it .eq. iterat_limit) then
         write(*,*)'GenPhot=',abs(qh-qh0)/qh,qh,qh0,Teff
      qh00=qh
C 
c
         CALL GENPHOTO(RSTAR,TEFF,ALST)
c
         it=0
         iflqh=1
      else
         iflq=0
      endif
c
      if (abs(teff-teff0)/teff .gt. (perc_temp)) then
         ifltemp=1
      else
         ifltemp=0
      endif
      it=it+1
c
      XJ(5)=XH*1.
      XJ(6)=XH*1.E-5
C
C     Helium I,II
C
      XJ(7)=XHE*1.
      XJ(8)=XHE*1.E-5
C
C     Carbon III,IV,V
C
      XJ(9)=XC*1.E-4
      XJ(10)=XC*1.E-4
      XJ(11)=XC*1.E-4
C
C     Oxygen II,III,IV,V,VI
C
      XJ(12)=XO*1.E-4
      XJ(13)=XO*1.E-4
      XJ(14)=XO*1.E-4
      XJ(15)=XO*1.E-4
      XJ(16)=XO*1.E-4
C
C     Neon II,III,IV,V,VI
C
      XJ(17)=XNE*1.E-4
      XJ(18)=XNE*1.E-4
      XJ(19)=XNE*1.E-4
      XJ(20)=XNE*1.E-4
      XJ(21)=XNE*1.E-4
C
C     Sulphur III,IV,V,VI
C
      XJ(22)=XS*1.E-4
      XJ(23)=XS*1.E-4
      XJ(24)=XS*1.E-4
      XJ(25)=XS*1.E-4
C
C     Nitrogen II,III,IV,V,VI
C
      XJ(26)=XO*1.E-4
      XJ(27)=XO*1.E-4
      XJ(28)=XO*1.E-4
      XJ(29)=XO*1.E-4
      XJ(30)=XO*1.E-4
C
c    C
      CV=1.5
C
      Do i=1,NRW
         R=(float(i)-1.)*DR+RMIN
         rho=amdot/(4.*pi*vwind*R**2.)
         DEN=rho/AMASS/AMOL
         if (ichem .eq. 1) then
            DEN=DEN0
            RHO=DEN*AMASS*AMOL
            Twind=T0
            vwind=V0
         endif
         U(i,3)=RHO
         U(i,2)=vwind*rho
         U(i,1)=DEN*AK*twind/(GAMMA-1.)+0.5*rho*vwind**2
         U(i,4)=DEN
         do IEQ=NEQ1,NEQ
            U(i,IEQ)=DEN*XJ(ieq)
         end do
c
         Enddo
C     
      RETURN
      END
C
      function aionizphotons(TEFF,alst)
      implicit real*8 (a-h,o-z)
      PARAMETER(AMASS=1.674D-24,AK=1.38D-16,PC=3.08E18
     &     ,pi=3.141592, ah=6.62d-27, clight=3d10)
c
      vh=13.6
      vf=200.
c
      vh=vh*1.6e-12/ah
      vf=vf*1.6e-12/ah  
      dv=(vf-vh)/1000
      v=vh
      aioni=0.
      do while (v .lt. vf)
         fluxblack=fluxbb(TEFF,v)
         rstar2=alst/(4.*pi*5.670d-05*teff**4.)
c
         alv=4.*pi*rstar2*fluxblack
         aioni=aioni+alv*dv/(ah*v)
         v=v+dv
      enddo
      aionizphotons=aioni
c
      return
      end
c
      function fluxbb(Teff,anu)
      implicit real*8 (a-h,o-z)
      PARAMETER(AMASS=1.674D-24,AK=1.38D-16,PC=3.08E18
     &     ,pi=3.141592, ah=6.62d-27, clight=3d10)
c
      fluxb=2.*pi*ah*anu**3./clight**2.
      fluxbb=fluxb/(exp(ah*anu/(aK*Teff))-1.)
c
      return
      end 
c      
      subroutine output(U,dtprint,itprint)
      implicit real*8 (a-h,o-z)
      PARAMETER(NR=15000,NEQ=30)
      DIMENSION U(NR,NEQ),PRIM(NEQ)
      Dimension al4861(NR),al3727(NR),al4686(NR),al5007(NR),al5876(NR)
     &     ,al6716(NR),al6730(NR),amassion(NR),COLHEAT(2,NR)
      Dimension primp(NEQ)
      Dimension emhb(NR)   
      CHARACTER*40 FILE1,FILE2
      character*20 modelfile
      PARAMETER(AMASS=1.674D-24,AK=1.38D-16,PC=3.08E18,PI=3.141592)
      PARAMETER(ALSUN=4e33,AMSUN=2e33)
      COMMON/RADIUS/DR,RMIN
      COMMON/ABUND/XH,XHE,XC,XN,XO,XS,XAr,XNe
      Common/amoln/AMOL
      common/modelf/modelfile
      COMMON/COLHEATF/COLHEAT
      common/wohd/DEN0,T0,V0
      common/emhbe/emhb
      common/qhbb/qhbb,qhaprox 
      save      
C     
      ITP=ITPRINT+1
C     
      IF(ITP.LT.100) THEN
c         file1='/DAT/'modelfile(*)'_'ITP'.dat'
         WRITE(FILE2,1001)modelfile(1:12),ITP
         WRITE(FILE1,1000)modelfile(1:12),ITP

 1000    FORMAT('./',A,'/points',I3.3,'.dat')
 1001    FORMAT('./',A,'/lumin',I3.3,'.dat')
      END IF
      OPEN(UNIT=11,FILE=FILE1,STATUS='UNKNOWN')
      OPEN(UNIT=8,FILE=FILE2,STATUS='UNKNOWN')
C
      TIME=ITPRINT*dtprint
c      WRITE(11,1101)TIME,Teff,LSTAR,qh
      call alumsandmass(U,al4861,al3727,al4686,
     &     al5007,al5876,al6716,al6730,amassion)
      WRITE(8,*) 
     & '       Radii    ',
     & '    gas density ',
     & '      velocity  ',
     & '  temperature   ',
     & '   n_hi         ',
     & '   n_hii        ',
     & '   n_hei        ',
     & '   n_heii       ', 
     & '   n_oii        ',   
     & '   n_oiii       ',          
     & '   n_sii        ', 
     & '   l4861        ', 
     & '   l3727        ',        
     & '   l4686        ',   
     & '   l5007        ',  
     & '   l5876        ',  
     & '   l6716        ',  
     & '   l6730        ',  
     & '   M_HII        ',
     & '   COOLING      ',  
     & '   HEATING      ',
     & '   M_ion_obs    ',  
     & '   n_e_sii      ',
     & '   denelec      ',
     & '   em4861       ',
     & '    QH_bb       ',
     & '    Rst         '    
      WRITE(11,*) 
     & '     Radii      ',
     & '     pressure   ',
     & '      velocity  ',
     & '  rho           ',
     & '  temperature   ',
     & '   n_gas        ',
     & '   n_hi         ',
     & '   n_hii        ',
     & '   n_hei        ', 
     & '   n_heii       ',   
     & '   n_heiii      ',          
     & '   n_cii        ', 
     & '   n_ciii       ', 
     & '   n_civ        ',        
     & '   n_cv         ',   
     & '   n_oii        ', 
     & '   n_oiii       ', 
     & '   n_oiv        ',        
     & '   n_ov         ',
     & '   n_ovi        ',
     & '   n_neii       ', 
     & '   n_neiii      ', 
     & '   n_neiv       ',        
     & '   n_nev        ',
     & '   n_nevi       ',
     & '   n_sii        ',
     & '   n_siii       ',
     & '   n_siv        ',
     & '   n_sv         ',
     & '   n_svi        ',
     & '   n_nii        ',
     & '   n_niii       ',
     & '   n_niv        ',
     & '   n_nv         ',
     & '   n_nvi        ',
     & '   n_e          '
      
C
      Do i=1,NR
c   
         R=i*DR
         
         call UPRIM(PRIM(1:NEQ),T,U(i,1:NEQ))
C     SII density calculation                                                                            
      DENS=PRIM(4)*XS
      DENSII=DENS-PRIM(22)-PRIM(23)-PRIM(24)-PRIM(25)
      DENHeIII=PRIM(4)*XHe-PRIM(7)-PRIM(8)
      DENCII=PRIM(4)*XC-PRIM(9)-PRIM(10)-PRIM(11)
      DENSII=max(DENSII,0.)
      DENHeIII=max(DENHeIII,0.)
      DENCII=max(DENCII,0.)
      cociente_sii=(al6730(i)+1.0e-5)/(al6716(i)+1e-15)
c     calculation of electron density of SII
c    Fitting target of lowest sum of squared absolute error = 1.3167229088376683E-01
      amplitude=4.4287962750101889E-01
      center=1.5044187041344719E+00
      width=1.8449798641798925E+00
      Offset= 3.3226999462955078E+00
      elecden_sii_log=amplitude*tan(pi*(cociente_sii-center)/width)
     &     +Offset
      elecden_sii=10**(elecden_sii_log)
      am_ion_obs=37.5*(al4861(i)/ALSUN)/elecden_sii
      denelec=delec(PRIM)

C
      alfaH=2.56e-13
      rst=(3.*qhbb*(1.-exp(-DEN0*alfaH*time))/(4.*pi*alfaH*DEN0**2.))
      rst=rst**0.3333333
      WRITE(8,1100)R,PRIM(4),PRIM(2),T
     &     ,PRIM(5),PRIM(6),PRIM(7),PRIM(8),PRIM(12),PRIM(13),DENSII
     &         ,al4861(i),al3727(i),al4686(i),al5007(i),al5876(i)
     &     ,al6716(i),al6730(i),amassion(i),COLHEAT(1,i)
     &     ,COLHEAT(2,i),am_ion_obs*AMSUN,elecden_sii,denelec,emhb(i)
     &     ,qhbb,rst
      WRITE(11,1101)R,PRIM(1),PRIM(2),PRIM(3),T
     &     ,PRIM(4),PRIM(5),PRIM(6),PRIM(7),PRIM(8),DENHeIII,DENCII
     &     ,PRIM(9),PRIM(10)
     &     ,PRIM(11),PRIM(12),PRIM(13),PRIM(14),PRIM(15),PRIM(16)
     &     ,PRIM(17),PRIM(18),PRIM(19),PRIM(20),PRIM(21),DENSII
     &     ,PRIM(22),PRIM(23),PRIM(24),PRIM(25),PRIM(26),PRIM(27)
     &     ,PRIM(28),PRIM(29),PRIM(30),denelec
      Enddo
 1100 FORMAT(27(' ',ES15.4))
 1101 FORMAT(36(' ',ES15.4))
c
      close(8)
      close(11)
C                                                                                                       
      Return
      End
C     
      subroutine alumsandmass(U,al4861,al3727,al4686,
     &     al5007,al5876,al6716,al6730,amassion)
      implicit real*8 (a-h,o-z)
      PARAMETER(NR=15000,NEQ=30)
      PARAMETER(AMASS=1.674D-24,AK=1.38D-16,pi=3.141592)
      PARAMETER(NCOL=30,NLINE=230)
      Dimension dens(NCOL),temp(NLINE),coef4861(NCOL,NLINE)
     &     ,coef3727(NCOL,NLINE),coef4686(NCOL,NLINE)
     &     ,coef5007(NCOL,NLINE),coef5876(NCOL,NLINE)
     &     ,coef6716(NCOL,NLINE),coef6730(NCOL,NLINE)
      Dimension U(NR,NEQ),PRIM(NEQ)
      Dimension al4861(NR),al3727(NR),al4686(NR),al5007(NR),al5876(NR)
     &     ,al6716(NR),al6730(NR),amassion(NR)
      Dimension emhb(NR)   
      COMMON/RADIUS/DR,RMIN
      Common/coef/coef4861,coef3727,coef4686,coef5007,coef5876
     &     ,coef6716,coef6730
      Common/dentemp/dens,temp
      Common/amoln/AMOL
      COMMON/ABUND/XH,XHE,XC,XN,XO,XS,XAr,XNe
      common/celd/NRW
      Common/emhbe/emhb 
      save      
      
C
      call opentablesiones
C
      emhb(1:NR)=0.
      al4861(1:NR)=0.
      al3727(1:NR)=0.
      al4686(1:NR)=0.
      al5007(1:NR)=0.
      al5876(1:NR)=0.
      al6716(1:NR)=0.
      al6730(1:NR)=0.
      amassion(1:NR)=0.
C                                                                                                                 
      Do i=NRW+1,NR
         R=i*DR
         call UPRIM(PRIM(1:NEQ),T,U(i,1:NEQ))
C     ELECTRONIC DENSITY CALCULATION
         DENHeIII=PRIM(4)*XHe-PRIM(7)-PRIM(8)
         DENCII=PRIM(4)*XC-PRIM(9)-PRIM(10)-PRIM(11)
C
         DENSII=PRIM(4)*XS-PRIM(22)-PRIM(23)-PRIM(24)-PRIM(25)
         DENHeIII=max(DENHeIII,0.)
         DENCII=max(DENCII,0.)
         DENSII=max(DENSII,0.)
C
         denelec=delec(PRIM)
C
         Fact=4.*pi*R**2.*denelec*DR
C
         If (T .lt. 1.e6) Then
            T1=T
            if(T .lt. 200) T1=200.
C     Hbeta
            em4861=bidiminterpol(T1,denelec,temp,dens,coef4861)
            emhb(i)=em4861
            if (em4861 .lt. 0) em4861=0.
            al4861(i)=Fact*PRIM(6)*em4861+al4861(i-1)
C     [OII]3727
            em3727=bidiminterpol(T1,denelec,temp,dens,coef3727)
            if (em3727 .lt. 0) em3727=0.
            al3727(i)=Fact*PRIM(12)*em3727+al3727(i-1)
c
C     [HeII]4686     
C
            em4686=bidiminterpol(T1,denelec,temp,dens,coef4686)
            if (em4686 .lt. 0) em4686=0.
            al4686(i)=Fact*DENHeIII*em4686+al4686(i-1)
c
C     [OIII]5007
            em5007=bidiminterpol(T1,denelec,temp,dens,coef5007)
            if (em5007 .lt. 0) em5007=0.
            al5007(i)=Fact*PRIM(13)*em5007+al5007(i-1)
c
C     [HeI]5876
            em5876=bidiminterpol(T1,denelec,temp,dens,coef5876)
            if (em5876 .lt. 0) em5876=0.
           al5876(i)=Fact*PRIM(8)*em5876+al5876(i-1)
c
C     [SII]6716
C     
           em6716=bidiminterpol(T1,denelec,temp,dens,coef6716)
           if (em6716 .lt. 0) em6716=0.
           al6716(i)=Fact*DENSII*em6716+al6716(i-1)
c
C     [SII]6730
           em6730=bidiminterpol(T1,denelec,temp,dens,coef6730)
           if (em6730 .lt. 0) em6730=0.
            al6730(i)=Fact*DENSII*em6730+al6730(i-1)

C
         Else
            al4861(i)=al4861(i-1)
            al3727(i)=al3727(i-1)   
            al4686(i)=al4686(i-1)
            al5007(i)=al5007(i-1)
            al5876(i)=al5876(i-1)
            al6716(i)=al6716(i-1)
            al6730(i)=al6730(i-1)
         Endif
C     Ionized Mass
         amassion(i)=Fact*AMASS*AMOL+amassion(i-1)      
      Enddo
C         
      close(19)
      return
      end
C
      subroutine opentablesiones
      implicit real*8 (a-h,o-z)
      PARAMETER(NCOL=30,NLINE=230)
      Dimension dens(NCOL),temp(NLINE),coef4861(NCOL,NLINE)
     &     ,coef3727(NCOL,NLINE),coef4686(NCOL,NLINE)
     &     ,coef5007(NCOL,NLINE),coef5876(NCOL,NLINE)
     &     ,coef6716(NCOL,NLINE),coef6730(NCOL,NLINE)
C                                                                                                                 
      Common/coef/coef4861,coef3727,coef4686,coef5007,coef5876
     &     ,coef6716,coef6730
      Common/dentemp/dens,temp
      
      save
       
      OPEN(UNIT=19,FILE='ABELION.EMISSIVITIES',STATUS='UNKNOWN')
C
 96   CONTINUE
      READ(19,*,ERR=96)(dens(j),j=1,NCOL)
      Do i=1,NLINE
         READ(19,*) temp(i),(coef4861(j,i),j=1,NCOL)
      EndDo
 97   CONTINUE
      READ(19,*,ERR=97)(dens(j),j=1,NCOL)
      Do i=1,NLINE
         READ(19,*)temp(i),(coef3727(j,i),j=1,NCOL)
      EndDo
 98   CONTINUE
      READ(19,*,ERR=98)(dens(j),j=1,NCOL)
      Do i=1,NLINE
         READ(19,*)temp(i),(coef4686(j,i),j=1,NCOL)
      EndDo
 99   CONTINUE
      READ(19,*,ERR=99)(dens(j),j=1,NCOL)
      Do i=1,NLINE
         READ(19,*)temp(i),(coef5007(j,i),j=1,NCOL)
      EndDo
 100  CONTINUE
      READ(19,*,ERR=100)(dens(j),j=1,NCOL)

      Do i=1,NLINE
         READ(19,*) temp(i),(coef5876(j,i),j=1,NCOL)
      EndDo
 101  CONTINUE
      READ(19,*,ERR=101)(dens(j),j=1,NCOL)
      Do i=1,NLINE
         READ(19,*) temp(i),(coef6716(j,i),j=1,NCOL)
      EndDo
 102  CONTINUE
      READ(19,*,ERR=102)(dens(j),j=1,NCOL)
      Do i=1,NLINE
         READ(19,*)temp(i),(coef6730(j,i),j=1,NCOL)
      EndDo
C             
C                                                                                                               
      return
      end
C                                 
      function bidiminterpol(Temp,den,ta,da,coef)
      implicit real*8 (a-h,o-z)
      PARAMETER (NCOL=30,NLINE=230)
      Dimension ta(NLINE),da(NCOL),coef(NCOL,NLINE)
      
      save      

      i=1
      k=2
      kd=2
      Do While (ta(i) .le. Temp .and. i .lt. NLINE)
          i=i+1
          k=i
      EndDo
      t0=ta(k-1)
      t1=ta(k)
      i=1
      Do While (da(i) .le. den .and. i .lt. NCOL)
         i=i+1
         kd=i
      EndDo
      d0=da(kd-1)
      d1=da(kd)
C                                              
      tc=(Temp-t0)/(t1-t0)
      dc=(den-d0)/(d1-d0)
C
      emissivity=(1.-tc)*(1.-dc)*coef(kd-1,k-1)+
     &     tc*(1.-dc)*coef(kd-1,k)+
     &     tc*dc*coef(kd,k)+
     &     (1.-tc)*dc*coef(kd,k-1)
C     
      bidiminterpol=emissivity
C                                                                                                                 
      Return
      End
C      
      Function ainterpol(X,FX,FY,Nmax)
C         
      implicit real*8 (a-h,o-z)
      DIMENSION FX(Nmax),FY(Nmax)
      
      save

      If(abs(X) .lt. 1.d-100) Then 
      ainterpol=FY(1)                                                   
      Else 
         Do i=2,Nmax   
            If (FX(i) .ge. X) Then                                             
               pend=(FY(i)-FY(i-1))/(FX(i)-FX(i-1))
               ainterpol=pend*(X-FX(i-1))+FY(i-1)
               return 
            EndIF        
         Enddo 
      Endif                                                                     
      Return                                                                 
      End 
C
      subroutine hllcfluxes(U,f)
      implicit real*8 (a-h,o-z)
      PARAMETER(NR=15000,NEQ=30)
       DIMENSION U(NR,NEQ),f(NR,NEQ),primL(NEQ),primR(NEQ),ff(NEQ)
      DIMENSION PRIM(NEQ),uu(NEQ)
      save      
      
      do i=1,NR-1
C
            call UPRIM(PRIML,T,U(i,1:NEQ))
C
            call UPRIM(PRIMR,T,U(i+1,1:NEQ))
C
            call prim2hllc(primL, primR, ff)
C
            f(i,1:NEQ)=ff(1:NEQ)
C
         end do
      return
      end
C
      subroutine prim2hllc(primL,primR,ff)
      implicit real*8 (a-h,o-z)
      PARAMETER (NEQ=30)
      DIMENSION fl(NEQ), fr(NEQ), ust(NEQ), uL(NEQ), uR(NEQ), uu(NEQ)
      DIMENSION primL(NEQ),primR(NEQ),ff(NEQ)
      PARAMETER(AMASS=1.674D-24,AK=1.38D-16,PC=3.08E18)
      COMMON/GAMMA/GAMMA
      Common/amoln/AMOL
      
      save      
C
      csl= dsqrt(GAMMA*primL(1)/primL(3))
      csr= dsqrt(GAMMA*primR(1)/primR(3))

      sl = dmin1(primL(2)-csl, primR(2)-csr )
      sr = dmax1(primL(2)+csl, primR(2)+csr )
C
      sst=(primR(1)-primL(1)+primL(3)*primL(2)*(sl-primL(2))
     &     -primR(3)*primR(2)*(sr-primR(2)))
     &     /(primL(3)*(sl-primL(2))-primR(3)*(sr-primR(2)))
C

      if(sl .gt. 0. ) then
         
         call eulerfluxes(primL,ff)
         
      else if (sr .lt. 0.) then
c         
         call eulerfluxes(primR,ff)
C
      else if (sst .ge. 0.) then
         ust(3)=primL(3)*(sl-primL(2) )/(sl-sst)
         ust(2)=ust(3)*sst
C
         ek=(0.5*primL(3)*primL(2)**2)+primL(1)/(gamma-1.)
         ust(1)=ust(3)*(ek/primL(3)+ (sst-primL(2))*
     &        (sst+primL(1)/( primL(3)*(sl-primL(2) ) ) ) )

         call eulerfluxes(primL(:),fl)
         call primu     (primL(:),uu)

         ff(1:NEQ)= fl(1:NEQ)+ sl * (ust(1:NEQ)-uu(1:NEQ) )

      else if (sst <= 0.) then

         ust(3)=primR(3)*(sr-primR(2)) /(sr-sst)
         ust(2)=ust(3)*sst
!
         ek=(0.5*primR(3)*primR(2)**2)+primR(1)/(gamma-1.)
         ust(1)=ust(3)*(ek/primR(3)+ (sst-primR(2))*
     &        (sst+primR(1)/( primR(3)*(sr-primR(2)) ) ) )

         call eulerfluxes(primR(:),fr)
         call primu     (primR(:),uu)
         ff(:)= fr(:)+sr*(ust(:)-uu(:) )

      else
         write(*,*)'Error in HLLC', sl, sr, sst
         stop
         
      endif
      return
      end
      !
      subroutine sources(i,pp,ss)
      implicit real*8 (a-h,o-z)
      PARAMETER (NEQ=30)
      DIMENSION ss(NEQ), pp(NEQ)
      parameter (alpha=2.)
      COMMON/GAMMA/GAMMA
      COMMON/RADIUS/DR,RMIN
      save      
      
C
      R=(float(i)-1.)*DR
      Etot=0.5*pp(3)*pp(2)*pp(2)+pp(1)/(GAMMA-1.)
      term=+alpha/R
      ss(1)=term*pp(2)*(Etot+pp(1))
      ss(2)=term*pp(3)*pp(2)*pp(2)
      ss(3)=term*pp(3)*pp(2)
C
      do ieq=4,NEQ
         ss(ieq)=term*pp(ieq)*pp(2)
      enddo
c      
      return
      end
C
      subroutine eulerfluxes(pp,ff)
      implicit real*8 (a-h,o-z)
      PARAMETER (NEQ=30)
      DIMENSION ff(NEQ), pp(NEQ)
      PARAMETER(AMASS=1.674D-24,AK=1.38D-16)
      COMMON/GAMMA/GAMMA
      Common/amoln/AMOL
      save      
C
      Etot=0.5*pp(3)*pp(2)*pp(2)+pp(1)/(GAMMA-1.)
      ff(1)=pp(2)*(Etot+pp(1))
      ff(2)=pp(3)*pp(2)*pp(2)+pp(1)
      ff(3)=pp(3)*pp(2)
C     
      do ieq=4,NEQ
         ff(ieq)=pp(ieq)*pp(2)
      enddo
      return
      end 
C
      subroutine boundaries(up)
      implicit real*8 (a-h,o-z)
      PARAMETER(NR=15000,NEQ=30)
      DIMENSION up(NR,NEQ)
      save      
C
      up(1,1:NEQ)=up(2,1:NEQ)
      up(NR,1:NEQ)=up(NR-1,1:NEQ)
C
      return
      end
C
      SUBROUTINE TSTEP(U,DT,TIME)
C
      implicit real*8 (a-h,o-z)
      PARAMETER(NR=15000,NEQ1=5,NEQ=30)
      DIMENSION U(NR,NEQ),UP(NR,NEQ),HICOL(NR),HEICOL(NR),HEIICOL(NR)
      DIMENSION F(NR,NEQ),PRIML(NEQ),UUL(NEQ),FFL(NEQ),COLHEAT(2,NR)
      DIMENSION PRIMR(NEQ),UUR(NEQ),FFR(NEQ),UPP(NR,NEQ)
      DIMENSION SS(NEQ),UU(NEQ),UPREV(NR,NEQ)
      COMMON/RADIUS/DR,RMIN
      PARAMETER(AMASS=1.674D-24,AK=1.38D-16,YR=3.15E7)
      PARAMETER(pi=3.141592)
      real*8 ALST
      COMMON/HDSTUFF/Co,eta
      Common/amoln/AMOL
      COMMON/COLHEATF/COLHEAT
      common/celd/NRW
      common/nummethod/imet,ichem
      common/iterat/iterat_limit,iflqh,ifltemp
      common/qhbb/qhbb,qhaprox 
      save
c   Numerical integrator
c     Simple Lax solver = 3
c     Maccormack corrector-predictor method  = 2
c     Godunov Methodo 1er orden= 1
c     Godunov Methodo 2o. orden= 0
      dtr=DT/(DR)
      UPREV(1:NR,1:NEQ)=U(1:NR,1:NEQ)
C Lax Solver
      if (imet .eq. 3) then
         do i=2,NR-1
            UUL(1:NEQ)=U(i-1,1:NEQ)
            UUR(1:NEQ)=U(i+1,1:NEQ)
            call UPRIM(PRIML,T,UUL)
            call UPRIM(PRIMR,T,UUR)
            call eulerfluxes(PRIML,FFL)
            call eulerfluxes(PRIMR,FFR)
           call sources(i,PRIML,SS)
            up(i,1:NEQ)=u(I,1:NEQ)-dtr*(FFR(1:NEQ)-FFL(1:NEQ))
     &           -dt*ss(1:NEQ)
     &        +eta*(U(i+1,1:NEQ)+U(i-1,1:NEQ)-2.*U(i,1:NEQ))
 
         end do
         call boundaries(UP)
         UPP(1:NR,1:NEQ)=UP(1:NR,1:NEQ)
      endif
c
      if (imet .eq. 2) then
C
      do i=2,NR-1
         UUL(1:NEQ)=U(i,1:NEQ)
         UUR(1:NEQ)=U(i+1,1:NEQ)
         call UPRIM(PRIML,T,UUL)
         call UPRIM(PRIMR,T,UUR)
         call eulerfluxes(PRIML,FFL)
         call eulerfluxes(PRIMR,FFR)

         up(i,1:NEQ)=u(I,1:NEQ)-dtr*(FFR(1:NEQ)-FFL(1:NEQ))
     &        +eta*(U(i+1,1:NEQ)+U(i-1,1:NEQ)-2.*U(i,1:NEQ))
         call sources(i,PRIML,SS)
         up(i,1:NEQ)=up(i,1:NEQ)-dt*ss(1:NEQ)
      end do
C
      call boundaries(UP)
C
      do i=2,NR-1
         UUL(1:NEQ)=UP(i-1,1:NEQ)
         UUR(1:NEQ)=UP(i,1:NEQ)
         call UPRIM(PRIML,T,UUL)
         call UPRIM(PRIMR,T,UUR)
         call eulerfluxes(PRIML,FFL)
         call eulerfluxes(PRIMR,FFR)
         upp(i,1:NEQ)=0.5*(u(I,1:NEQ)+up(I,1:NEQ))-dtr*
     &        (FFR(1:NEQ)-FFL(1:NEQ))
     &        +eta*(up(i+1,1:NEQ)+up(i-1,1:NEQ)-2.*up(i,1:NEQ))
         call sources(i,PRIMR,SS)
         upp(i,1:NEQ)=upp(i,1:NEQ)-0.5*dt*ss(1:NEQ)

      end do
c
      call boundaries(UPP)
      endif
C
      if (imet .eq. 1) then
C
         call hllcfluxes(U,f)
C
         do i=2,nR-1
            upp(i,1:NEQ)=u(i,1:NEQ)-DTR*(f(i,1:NEQ)-f(i-1,1:NEQ))
            call UPRIM(PRIMR,T,UPP(i,1:NEQ))
            call sources(i,PRIMR,SS)
            upp(i,1:NEQ)=upp(i,1:NEQ)-dt*ss(1:NEQ)
         end do
         call boundaries(UPP)
C
      endif
      if (imet .eq. 0) then
c         
         call hllcfluxes(U,f)
c     
         do i=2,nR-1
C     
            DTR2=0.5*DTR
C     
            up(i,1:NEQ)=u(i,1:NEQ)-DTR2*(f(i,1:NEQ)-f(i-1,1:NEQ))
C     
         end do
C     
         call boundaries(UP)
C     
         call hllcfluxes2(UP,f)
C     
         do i=2,nR-1
            upp(i,1:NEQ)=u(i,1:NEQ)-DTR*(f(i,1:NEQ)-f(i-1,1:NEQ))
            call UPRIM(PRIMR,T,UP(i,1:NEQ))
            call sources(i,PRIMR,SS)
            upp(i,1:NEQ)=upp(i,1:NEQ)-dt*ss(1:NEQ)
         end do
         call boundaries(UPP)
C     
      endif
c
      dtss=1e30
      sanhiialf=0.
      sanhii_r=0.
      sanhii_l=0.
C
      Do I=NRW+1,NR
C
         call UPRIM(PRIMR,TR,UPREV(i,1:NEQ))
         call UPRIM(PRIML,TL,UPP(i,1:NEQ))
c
      anhiirho_r=max(primr(6)/primr(4),0.01)
      anhiirho_l=max(priml(6)/priml(4),0.01)
      dnhiirho=abs(anhiirho_r-anhiirho_l)
      dnhiirhoDT=dnhiirho/DT
      dtp=(0.1*anhiirho_l)/dnhiirhoDT
cC
      R=i*DR
      anhii_r=primr(6)
      anhii_l=priml(6)
      denelec_l=delec(priml)
cc
cc     Recombination rate case B, Pequignot et al. 1991,A&A,251, 680
cc     
      apqui=4.309
      bpqui=-0.6166
      cpqui=0.6703
      dpqui=0.5300
      tpqui=tl/1e4
cc      
      alfhb=1e-13*apqui*tpqui**bpqui/(1.+cpqui*tpqui**dpqui)
cc
cc
      sanhiialf=sanhiialf+
     & 4.*pi*r**2*dr*(denelec_l*anhii_l)*alfhb
      sanhii_r=sanhii_r+
     & 4.*pi*r**2*dr*anhii_r
      sanhii_l=sanhii_l+
     & 4.*pi*r**2*dr*anhii_l
      if(dtp.lt.dtss) then
        dtss=dtp
        anhiir=anhiirho_r
        anhiil=anhiirho_l
        dn=dnhiirhoDT
        ii=I
      end if
c
C
c
      EndDo
      dsanhiidt=(sanhii_l-sanhii_r)/dt
C
C    For the ATOMIC DT we solve dt=(ne/nh)/d(ne/nh)/dt
C
      Nsteps=int(DT/dtss)
      Nsteps=max(1,Nsteps)
c      Nsteps=1
      write(*,*) 'TestAlex->',qhbb,sanhiialf,dsanhiidt,
     &     qhbb-(sanhiialf+dsanhiidt),Nsteps
c     
      DTS=DT/float(Nsteps)
      Do isteps=1,Nsteps
         CALL COLDENS(UPP,HICOL,HEICOL,HEIICOL)
         DO I=NRW+2,NR
            CALL ATOMIC(UPP,I,DTS,HICOL,HEICOL,HEIICOL,AL,HEAT)
C
            COLHEAT(1,I)=AL
            COLHEAT(2,I)=HEAT

         END DO
      Enddo
C
      call uflow(UPP,TIME,DT,ALST,TEFF)
C      
      U(1:NR,1:NEQ)=UPP(1:NR,1:NEQ)
C
      RETURN
      END
C
      subroutine hllcfluxes2(U,f) 
      implicit real*8 (a-h,o-z)
      PARAMETER(NR=15000,NEQ=30)
      DIMENSION U(NR,NEQ),f(NR,NEQ),primL(NEQ),primR(NEQ),prim(neq),
     &     primc(nr,neq),grad(nr,neq),ff(NEQ)
      DIMENSION uu(NEQ)
      do i=1,NR
         call UPRIM(PRIM,T,U(i,1:NEQ))
         DO ieq=1,neq
            primc(i,ieq)=prim(ieq)
         end do
      end do
      do ieq=1,neq
         grad(1,ieq)=0.
         grad(nr,ieq)=0.
         do i=2,nr-1
            grad(i,ieq)=0.5*aver(primc(i,ieq)-primc(i-1,ieq),
     &           primc(i+1,ieq)-primc(i,ieq))
         end do
      end do
c
      do i=1,nr-1
         do ieq=1,neq
            priml(ieq)=primc(i,ieq)+grad(i,ieq)
            primr(ieq)=primc(i+1,ieq)-grad(i+1,ieq)
         end do
C
         call prim2hllc(primL, primR, ff)
C
         f(i,1:NEQ)=ff(1:NEQ)
C
      end do
      return
      end
C
      FUNCTION AVER(A,B)
C
C    Averaging function
C
      implicit real*8 (a-h,o-z)
      AB=A*B
      IF(AB.LE.0) THEN
         AVER=0.
      ELSE
         AVER=(AB*B+AB*A)/(A*A+B*B)
      END IF
C
      RETURN
      END  
c
      SUBROUTINE TIMESTEP(U,DT)
C
      implicit real*8 (a-h,o-z)
      PARAMETER(NR=15000,NEQ=30)
      DIMENSION U(NR,NEQ),PRIM(NEQ)
      PARAMETER(AMASS=1.674D-24,AK=1.38D-16)
      COMMON/RADIUS/DR,RMIN
      COMMON/GAMMA/GAMMA
      Common/amoln/AMOL
      COMMON/HDSTUFF/Co,eta
      save      
C
      del=1E30
      csmax=1e-30
      vmax=1e-30
C
      do i=1,NR
         call UPRIM(PRIM,T,U(i,1:NEQ))
         cs=dsqrt(GAMMA*prim(1)/(prim(3)))
         del=dmin1(del,DR/(abs(prim(2))+cs))
         vmax=max(prim(2),vmax)
         csmax=max(cs,csmax)
      end do
      DT=Co*del

      return
      END SUBROUTINE
C
      SUBROUTINE UPRIM(PRIM,T,UU)
C
      implicit real*8 (a-h,o-z)
      PARAMETER(NEQ=30,NEQ1=5)
      PARAMETER(AMASS=1.674D-24,AK=1.38D-16)
      DIMENSION PRIM(NEQ),UU(NEQ)
      COMMON/GAMMA/GAMMA
      Common/amoln/AMOL
            common/nummethod/imet,ichem
      save      
      
C
      DO IEQ=3,NEQ
         PRIM(IEQ)=UU(IEQ)
      END DO
      denelec=delec(prim)
c      
      DEN=PRIM(3)/AMASS/AMOL+denelec
      PRIM(2)=UU(2)/PRIM(3)
      PRIM(1)=(UU(1)-0.5*PRIM(3)*PRIM(2)*PRIM(2))*(GAMMA-1.)
c
      PRIM(4)=PRIM(3)/AMASS/AMOL
      PRIM(1)=max(PRIM(1),1.*AK*DEN)
      T=dmax1(PRIM(1)/(AK*DEN),200.)

C     
      RETURN
      END
C
      Function delec(prim)
      implicit real*8 (a-h,o-z)
      PARAMETER(NEQ=30)
      DIMENSION PRIM(NEQ)
      COMMON/ABUND/XH,XHE,XC,XN,XO,XS,XAr,XNe
      save
C   ELECTRONIC DENSITY CALCULATION
         DENHeIII=PRIM(4)*XHe-PRIM(7)-PRIM(8)
         DENCII=PRIM(4)*XC-PRIM(9)-PRIM(10)-PRIM(11)
         DENSII=PRIM(4)*XS-PRIM(22)-PRIM(23)-PRIM(24)-PRIM(25)
         DENHeIII=max(DENHeIII,0.)
         DENCII=max(DENCII,0.)
         DENSII=max(DENSII,0.)
C
         delec=0.
         delec=PRIM(6)+PRIM(8)+2.*DENHeIII
C    CARBON
c     &        +2.*PRIM(9)+3.*PRIM(10)+4.*PRIM(11)
C    OXYGEN
c     &        +PRIM(12)+2.*PRIM(13)+3.*PRIM(14)+4.*PRIM(15)+5.*PRIM(16)
C    NEON
c     &        +PRIM(17)+2.*PRIM(18)+3.*PRIM(19)+4.*PRIM(20)
c     &        +5.*PRIM(21)
C    SULFUR
c     &        +2.*PRIM(22)+3.*PRIM(23)+4.*PRIM(24)
c     &        +5.*PRIM(25)
C    NITROGEN
c     &        +PRIM(26)+2.*PRIM(27)+3.*PRIM(28)+4.*PRIM(29)+5.*PRIM(30)
c         if ((prim(6)+prim(8))/(delec) .lt. 0.7 ) then
C    OXYGEN
C    NEON
C    SULFUR
C    NITROGEN
c            stop
c         endif
            Return
            
         End

C
      SUBROUTINE PRIMU(PRIM,UU)
C
      implicit real*8 (a-h,o-z)
      PARAMETER(NEQ=30,NEQ1=5)
      DIMENSION PRIM(NEQ),UU(NEQ)
      COMMON/GAMMA/GAMMA
       save     
      
C     
      DO IEQ=3,NEQ
         UU(IEQ)=PRIM(IEQ)
      END DO
      UU(1)=PRIM(1)/(GAMMA-1.)+0.5*PRIM(3)*PRIM(2)*PRIM(2)
      UU(2)=PRIM(2)*PRIM(3)
C
      RETURN
      END
C
C
      SUBROUTINE ATOMIC(U,I,DT,HICOL,HEICOL,HEIICOL,AL,HEAT)
C
C    calculates the atomic source terms
C
      implicit real*8 (a-h,o-z)
      PARAMETER(NEQ=30,NEQ1=5,LEVMAX=5,NR=15000)
      DIMENSION PRIM(NEQ),PRIM1(NEQ),XE(NEQ),XJ(NEQ),
     * UU(NEQ),YH(0:1),YHE(0:2),
     &     YC(0:LEVMAX),YN(0:LEVMAX),YO(0:LEVMAX),YNE(0:LEVMAX),
     &     YS(0:LEVMAX),RM(NR)
      DIMENSION U(NR,NEQ),HICOL(NR),HEICOL(NR),HEIICOL(NR)
      PARAMETER(AMASS=1.674D-24,AK=1.38D-16)
      
      COMMON/RADIUS/DR,RMIN
      Common/amoln/AMOL
      COMMON/ABUND/XH,XHE,XC,XN,XO,XS,XAr,XNe
      save      
            
C
      R=(FLOAT(I)-1)*DR+RMIN                      
      RM(I)=1./R
      DIL=RM(I)**2
      DO IEQ=1,NEQ
         UU(IEQ)=U(I,IEQ)
      END DO
      CALL UPRIM(PRIM,T,UU)
C      
      DO IEQ=1,NEQ
         PRIM1(IEQ)=PRIM(IEQ)
      END DO
      T1=T
      P0=PRIM(1)
c
c     Set the photo-ionization rates (outside of the iteration loop)
c
      CALL PHOTOION(HICOL(I),HEICOL(I),HEIICOL(I),DIL)
c
c     Iterate to find a good temperature
c
         DO IT=1,3
            T1=DMAX1(T1,10.)
c           TAV=AVER(T1,T)
            TAV=0.5*(T1+T)
C     
            tav=t
            YH(1)=PRIM(6)/(XH*PRIM(4))
            YH(0)=PRIM(5)/(XH*PRIM(4))
            TOT=0.
            DO IS=0,1
               YH(IS)=DMAX1(YH(IS),0.)
               YH(IS)=DMIN1(YH(IS),1.)
               TOT=TOT+YH(IS)
            END DO
            TOT=1./TOT
            DO IS=0,1
               YH(IS)=YH(IS)*TOT
            END DO
C     
            YHE(1)=PRIM(8)/(XHE*PRIM(4))
            YHE(0)=PRIM(7)/(XHE*PRIM(4))
            YHE(2)=1.-YHE(1)-YHE(0)
            TOT=0.
            DO IS=0,2
               YHE(IS)=DMAX1(YHE(IS),0.)
               YHE(IS)=DMIN1(YHE(IS),1.)
               TOT=TOT+YHE(IS)
            END DO
            TOT=1./TOT
            DO IS=0,2
               YHE(IS)=YHE(IS)*TOT
            END DO
C     
            YC(2)=PRIM(9)/(XC*PRIM(4))
            YC(3)=PRIM(10)/(XC*PRIM(4))
            YC(4)=PRIM(11)/(XC*PRIM(4))
            YC(1)=1.-YC(2)-YC(3)-YC(4)
            TOT=0.
            DO IS=1,4
               YC(IS)=DMAX1(YC(IS),0.)
               YC(IS)=DMIN1(YC(IS),1.)
               TOT=TOT+YC(IS)
            END DO
            TOT=1./TOT
            DO IS=1,4
               YC(IS)=YC(IS)*TOT
            END DO
C     
            YO(1)=PRIM(12)/(XO*PRIM(4))
            YO(2)=PRIM(13)/(XO*PRIM(4))
            YO(3)=PRIM(14)/(XO*PRIM(4))
            YO(4)=PRIM(15)/(XO*PRIM(4))
            YO(5)=PRIM(16)/(XO*PRIM(4))
            YO(0)=1.-YO(1)-YO(2)-YO(3)-YO(4)-YO(5)
            TOT=0.
            DO IS=0,5
               YO(IS)=DMAX1(YO(IS),0.)
               YO(IS)=DMIN1(YO(IS),1.)
               TOT=TOT+YO(IS)
            END DO
            TOT=1./TOT
            DO IS=0,5
               YO(IS)=YO(IS)*TOT
            END DO
C     
            YNE(1)=PRIM(17)/(XNE*PRIM(4))
            YNE(2)=PRIM(18)/(XNE*PRIM(4))
            YNE(3)=PRIM(19)/(XNE*PRIM(4))
            YNE(4)=PRIM(20)/(XNE*PRIM(4))
            YNE(5)=PRIM(21)/(XNE*PRIM(4))
            YNE(0)=1.-YNE(1)-YNE(2)-YNE(3)-YNE(4)-YNE(5)
            TOT=0.
            DO IS=0,5
               YNE(IS)=DMAX1(YNE(IS),0.)
               YNE(IS)=DMIN1(YNE(IS),1.)
               TOT=TOT+YNE(IS)
            END DO
            TOT=1./TOT
            DO IS=0,5
               YNE(IS)=YNE(IS)*TOT
            END DO
C     
            YS(2)=PRIM(22)/(XS*PRIM(4))
            YS(3)=PRIM(23)/(XS*PRIM(4))
            YS(4)=PRIM(24)/(XS*PRIM(4))
            YS(5)=PRIM(25)/(XS*PRIM(4))
            YS(1)=1.-YS(2)-YS(3)-YS(4)-YS(5)
            TOT=0.
            DO IS=1,5
               YS(IS)=DMAX1(YS(IS),0.)
               YS(IS)=DMIN1(YS(IS),1.)
               TOT=TOT+YS(IS)
            END DO
            TOT=1./TOT
            DO IS=1,5
               YS(IS)=YS(IS)*TOT
            END DO
C     
            YN(1)=PRIM(26)/(XN*PRIM(4))
            YN(2)=PRIM(27)/(XN*PRIM(4))
            YN(3)=PRIM(28)/(XN*PRIM(4))
            YN(4)=PRIM(29)/(XN*PRIM(4))
            YN(5)=PRIM(30)/(XN*PRIM(4))
            YN(0)=1.-YN(1)-YN(2)-YN(3)-YN(4)-YN(5)
            TOT=0.
            DO IS=0,5
               YN(IS)=DMAX1(YN(IS),0.)
               YN(IS)=DMIN1(YN(IS),1.)
               TOT=TOT+YN(IS)
            END DO
            TOT=1./TOT
            DO IS=0,5
               YN(IS)=YN(IS)*TOT
            END DO
C     
            DE=PRIM(6)+PRIM(8)+2.*ABS(XHE*PRIM(4)-PRIM(7)-PRIM(8))
C
            CALL DORIC(DT,TAV,DE,PRIM(4),YH,YHE
     &           ,YC,YN,YO,YNE,YS,AL,HEAT)
c
c     heating by photoionization of H
c
            CV=1.5
C
            IF(TAV .LT. 400.) al=0.
            IF(TAV.GT.5.e4 .OR.TAV.LT.300.) THEN
C
               DAL=AL-HEAT
               PRIM1(1)=PRIM(1)-(DT*DAL)/CV
C
            ELSE
               AA=AL/(CV*P0)
               PRIM1(1)=(HEAT/AA+(CV*P0-HEAT/AA)*DEXP(-AA*DT))/CV
C 
            END IF
c            
            DDEN=PRIM(4)*(1.+XH*YH(1)+XHE*(YHE(1)+2.*YHE(2)))
            T1=PRIM1(1)/(AK*DDEN)
c            prim1(1)=max(prim1(1),1.*ak*DDEN)
         END DO
C     
         PRIM1(5)=(XH*PRIM(4))*YH(0)
         PRIM1(6)=(XH*PRIM(4))*YH(1)
C
         PRIM1(7)=(XHE*PRIM(4))*YHE(0)
         PRIM1(8)=(XHE*PRIM(4))*YHE(1)
C     
         PRIM1(9)=(XC*PRIM(4))*YC(2)
         PRIM1(10)=(XC*PRIM(4))*YC(3)
         PRIM1(11)=(XC*PRIM(4))*YC(4)
C     
         PRIM1(12)=(XO*PRIM(4))*YO(1)
         PRIM1(13)=(XO*PRIM(4))*YO(2)
         PRIM1(14)=(XO*PRIM(4))*YO(3)
         PRIM1(15)=(XO*PRIM(4))*YO(4)
         PRIM1(16)=(XO*PRIM(4))*YO(5)
C     
         PRIM1(17)=(XNE*PRIM(4))*YNE(1)
         PRIM1(18)=(XNE*PRIM(4))*YNE(2)
         PRIM1(19)=(XNE*PRIM(4))*YNE(3)
         PRIM1(20)=(XNE*PRIM(4))*YNE(4)
         PRIM1(21)=(XNE*PRIM(4))*YNE(5)
C     
         PRIM1(22)=(XS*PRIM(4))*YS(2)
         PRIM1(23)=(XS*PRIM(4))*YS(3)
         PRIM1(24)=(XS*PRIM(4))*YS(4)
         PRIM1(25)=(XS*PRIM(4))*YS(5)
C     
         PRIM1(26)=(XN*PRIM(4))*YN(1)
         PRIM1(27)=(XN*PRIM(4))*YN(2)
         PRIM1(28)=(XN*PRIM(4))*YN(3)
         PRIM1(29)=(XN*PRIM(4))*YN(4)
         PRIM1(30)=(XN*PRIM(4))*YN(5)
C
c
         CALL PRIMU(PRIM1,UU)
         U(I,1)=UU(1)
C
         DO IEQ=4,NEQ
            U(I,IEQ)=UU(IEQ)
         END DO
c     END IF
C     
      RETURN
      END
C
      SUBROUTINE COLDENS(U,HICOL,HEICOL,HEIICOL)
C
      implicit real*8 (a-h,o-z)
      PARAMETER(NEQ=30,NR=15000)
      DIMENSION U(NR,NEQ),HICOL(NR),HEICOL(NR),HEIICOL(NR),RM(NR)
      COMMON/RADIUS/DR,RMIN
      common/celd/NRW       
      save      
C
      HICOL(1:NRW)=0.
      HEICOL(1:NRW)=0.
      HEIICOL(1:NRW)=0.
c integration with trapezes

      DO I=NRW+2,NR
         HICOL(I)=HICOL(I-1)+0.5*DR*(U(I,5)+U(I-1,5))
         HEICOL(I)=HEICOL(I-1)+0.5*DR*(U(I,7)+U(I-1,7))
         HEIICOL(I)=HEIICOL(I-1)+0.5*DR*(U(I,8)+U(I-1,8))
      END DO
C
      RETURN
      END
C
C

C
C
C
c=======================================================================
c
c
      subroutine photoion (hcolum,he0colum,he1colum,dilu)
c
c
c=======================================================================
c
c             calculates photo-ionization rates and heating
c  
c                           garrelt mellema
c 
c                            04-mar-1997
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c                          introduction
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c-----------------------------------------------------------------------
c     goal
c-----------------------------------------------------------------------
c
c
c      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c          formats, parameters, declarations, common blocks, data
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c-----------------------------------------------------------------------
c     implicits
c-----------------------------------------------------------------------
c
c     
c-----------------------------------------------------------------------
c     the formats
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     parameters
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      parameter(numspec=5,levmax=5)
c
c-----------------------------------------------------------------------
c     declarations
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     common blocks
c-----------------------------------------------------------------------
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /abunda/ transports heavy element abundances 
c
      common /abunda/ abu_he,abu_c,abu_n,abu_o,abu_ne,abu_s
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /hydrogen/ hydrogen atomic data
c
      common /hydrogen/ bh00,albpow,colh0,temph0,hionen,
     &     hphot(0:100,3),hheat(0:100,3)
c
c     /helium/ helium atomic data
c
      common /helium/ bhe00,bhe10,alcpow,colhe(0:1),temphe(0:1),
     &     hephot(0:100,0:1,3),heheat(0:100,0:1,3)
c
c     /atomtab/ atomic data for metals
c
      common /atomtab/ rec(41,numspec,0:levmax),
     &     col(41,numspec,0:levmax),h0chex(41,numspec,0:levmax),
     &     h1chex(41,numspec,0:levmax),
     &     phot(0:100,numspec,0:levmax,3)
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /optdep/ constants needed to calculate optical depths
c
      common /optdep/ sigh,sighe0,sighe1,tf2h,tf3h,tf3he0
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /photrates/ photo-ionization rates
c     
      common /photrates/ phih,phihe(0:1),aphi(numspec,0:levmax),
     &     hvphih,hvphhe(0:1)
     
      save     

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      calculate ionization states and energy source terms
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c-----------------------------------------------------------------------
c     find the optical depths
c-----------------------------------------------------------------------
c
      tauh=sigh*hcolum
      tauhe0=sighe0*he0colum
      tauhe1=sighe1*he1colum
c
c-----------------------------------------------------------------------
c     find the table positions for the optical depth
c     from tenorio-tagle et al. (1983)
c-----------------------------------------------------------------------
c
c      print*,'tauh',tauh,sigh,hcolum
      tau1=dlog10(dmax1(1.0e-20,tauh))
      odpos1=dmin1(100.0,dmax1(0.0,1.0+(tau1+4.0)*12.5))
      iodpo1=int(odpos1)
      dodpo1=odpos1-float(iodpo1)
      iodp11=min(100,iodpo1+1)
      
      tau2=dlog10(dmax1(1.0e-20,tauh*tf2h+tauhe0))
      odpos2=dmin1(100.0,dmax1(0.0,1.0+(tau2+4.0)*12.5))
      iodpo2=int(odpos2)
      dodpo2=odpos2-float(iodpo2)
      iodp21=min(100,iodpo2+1)
      
      tau3=dlog10(dmax1(1.0e-20,tauh*tf3h+tauhe0*tf3he0+tauhe1))
      odpos3=dmin1(100.0,dmax1(0.0,1.0+(tau3+4.0)*12.5))
      iodpo3=int(odpos3)
      dodpo3=odpos3-float(iodpo3)
      iodp31=min(100,iodpo3+1)
c
c-----------------------------------------------------------------------
c     find the hydrogen photo-ionization rate
c-----------------------------------------------------------------------
c
      phih=(hphot(iodpo1,1)+(hphot(iodp11,1)-hphot(iodpo1,1))*dodpo1+
     &     hphot(iodpo2,2)+(hphot(iodp21,2)-hphot(iodpo2,2))*dodpo2+
     &     hphot(iodpo3,3)+(hphot(iodp31,3)-hphot(iodpo3,3))*dodpo3)*
     &     dilu
c
c
c-----------------------------------------------------------------------
c     find the helium photo-ionization rate
c-----------------------------------------------------------------------
c
      do n=0,1
         phihe(n)=(hephot(iodpo2,n,2)+(hephot(iodp21,n,2)-
     &        hephot(iodpo2,n,2))*dodpo2+
     &        hephot(iodpo3,n,3)+(hephot(iodp31,n,3)-
     &        hephot(iodpo3,n,3))*dodpo3)*dilu
      enddo
c
c-----------------------------------------------------------------------
c     find the photo-ionization rate for the metals
c-----------------------------------------------------------------------
c
      do l=0,levmax-1
         do isp=1,numspec
            aphi(isp,l)=(phot(iodpo1,isp,l,1)+(phot(iodp11,isp,l,1)-
     &           phot(iodpo1,isp,l,1))*dodpo1+
     &           phot(iodpo2,isp,l,2)+(phot(iodp21,isp,l,2)-
     &           phot(iodpo2,isp,l,2))*dodpo2+
     &           phot(iodpo3,isp,l,3)+(phot(iodp31,isp,l,3)-
     &           phot(iodpo3,isp,l,3))*dodpo3)*dilu
         enddo
      enddo

      do isp=1,numspec
c**         write(*,*) (aphi(isp,l),l=0,levmax-1)
      enddo
c
c-----------------------------------------------------------------------
c     find the hydrogen photo-ionization heating rate
c-----------------------------------------------------------------------
c
      hvphih=(hheat(iodpo1,1)+(hheat(iodp11,1)-
     &     hheat(iodpo1,1))*dodpo1+
     &     hheat(iodpo2,2)+(hheat(iodp21,2)-
     &     hheat(iodpo2,2))*dodpo2+
     &     hheat(iodpo3,3)+(hheat(iodp31,3)-
     &     hheat(iodpo3,3))*dodpo3)*dilu
c
c
c-----------------------------------------------------------------------
c     find the helium photo-ionization rate
c-----------------------------------------------------------------------
c
      do n=0,1
         hvphhe(n)=(heheat(iodpo2,n,2)+(heheat(iodp21,n,2)-
     &        heheat(iodpo2,n,2))*dodpo2+
     &        heheat(iodpo3,n,3)+(heheat(iodp31,n,3)-
     &        heheat(iodpo3,n,3))*dodpo3)*dilu
      enddo
c      if (tauh.gt.1.0) write(6,*) tauh,phih,tauhe0,phihe(0),
c     &     tauhe1,phihe(1)

      return
      end
c=======================================================================
c
c
      subroutine doric
     &     (dt,temp0,rhe,rhh,xh,xhe,xc,xn,xo,xne,xs,emin,eplus)
c
c
c=======================================================================
c
c             calculates time dependent ionization state and 
c
c                   radiative energy source terms
c  
c                           garrelt mellema
c 
c                             06-mar-1997
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c                          introduction
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c-----------------------------------------------------------------------
c     goal
c-----------------------------------------------------------------------
c
c     this subroutine calculates the ionization state of the gas and
c     at the same time the energy source terms (cooling by recombination 
c     and lines, heating by photo-ionization). The ionization state is 
c     calculated time dependently for all species.
c
c      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c          formats, parameters, declarations, common blocks, data
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c-----------------------------------------------------------------------
c     implicits
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     the formats
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     parameters
c       levmax: maximum number of ionization states (0 to levmax)
c       numspec: number of metals
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      parameter(numspec=5,levmax=5)
c
c-----------------------------------------------------------------------
c     declarations
c-----------------------------------------------------------------------
c
      dimension xh(0:1),xhe(0:2),xc(0:levmax),xn(0:levmax),
     &     xo(0:levmax),xne(0:levmax),xs(0:levmax)
c
c-----------------------------------------------------------------------
c     common blocks
c-----------------------------------------------------------------------
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /abunda/ transports heavy element abundances 
c
      common /abunda/ abu_he,abu_c,abu_n,abu_o,abu_ne,abu_s
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /hydrogen/ hydrogen atomic data
c
      common /hydrogen/ bh00,albpow,colh0,temph0,hionen,
     &     hphot(0:100,3),hheat(0:100,3)
c
c     /helium/ helium atomic data
c
      common /helium/ bhe00,bhe10,alcpow,colhe(0:1),temphe(0:1),
     &     hephot(0:100,0:1,3),heheat(0:100,0:1,3)
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /tabpos/ temperature and density positions in the tables
c
      common /tabpos/ tpos,dtpos,dpos,ddpos,itpos,itpos1,idpos,idpos1
      common /tabposgz/ tposgz,dtposgz,dposgz,ddposgz,itposgz
     &     ,itpos1gz,idposgz,idpos1gz

c     
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /photrates/ local photo-ionization rates
c     
      common /photrates/ phih,phihe(0:1),aphi(numspec,0:levmax),
     &     hvphih,hvphhe(0:1)
c

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     COOLING TABLES STUFF
      parameter(ntetable=49,ndetable=8,ncoolgz=ntetable*ndetable)
      
      save

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      calculate ionization states and energy source terms
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c-----------------------------------------------------------------------
c     find the hydrogen recombination rate at the local temperature
c-----------------------------------------------------------------------
c
c      brech0=2.07e-11/sqrt(temp0)*(0.4288+0.5*log(temph0/temp0)+
c     &    0.469*(temph0/temp0)**(-0.33333))
c      brech0=0.
      t2bpow=temp0**albpow
      brech0=bh00*t2bpow
c
c      find the hydrogen recombination rate at the local temperature
c     (Ferland 1980, PASP, 92,596)
c      if (temp0 .le. 2.6e4) then
c         coefbh00=2.53e-22
c         expbh00=-0.833
c      endif
c      if (temp0 .gt. 2.6e4) then
c         coefbh00=1.12e-20
c         expbh00=-1.20
c      endif
c     c
c     Recombination rate case B, Pequignot et al. 1991,A&A,251, 680
c     
      apqui=4.309
      bpqui=-0.6166
      cpqui=0.6703
      dpqui=0.5300
      tpqui=temp0/1e4
c      
      brech0=1e-13*apqui*tpqui**bpqui/(1.+cpqui*tpqui**dpqui)
c
c     
c-----------------------------------------------------------------------
c     find the helium recombination rate at the local temperature
c-----------------------------------------------------------------------
c
      breche0=bhe00*temp0**alcpow
      breche1=bhe10*t2bpow
c
c-----------------------------------------------------------------------
c     find the hydrogen collisional ionization rate at the local 
c     temperature
c-----------------------------------------------------------------------
c
      sqrtt0=dsqrt(temp0)
      acolh0=colh0*sqrtt0*dexp(-temph0/temp0)
c
c-----------------------------------------------------------------------
c     find the helium collisional ionization rate at the local 
c     temperature
c-----------------------------------------------------------------------
c
      acolhe0=colhe(0)*sqrtt0*dexp(-temphe(0)/temp0)
      acolhe1=colhe(1)*sqrtt0*dexp(-temphe(1)/temp0)
c
c-----------------------------------------------------------------------
c     determine the hydrogen and helium ionization states and 
c     electron density using an iterarive method 
c     (schmidt-voigt & koeppen 1987)
c-----------------------------------------------------------------------
c
      rhe0=rhe
      dt5=dt
      xh1old=xh(1)
      xh0old=xh(0)
      xhe1old=xhe(1)
      xhe2old=xhe(2)
c
c**      write(*,*) acolh0,brech0,phih
      do itt=1,5
         ah0=phih+rhe*acolh0
         delth=ah0+rhe*brech0
         eqxh1=ah0/delth
         eqxh0=rhe*brech0/delth
         deltht=delth*dt5
         ee=dexp(-deltht)
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
         xh(1)=(xh1old-eqxh1)*ee+eqxh1
         xh(0)=(xh0old-eqxh0)*ee+eqxh0
c         xh(1)=eqxh1
c         xh(0)=eqxh0
c         rhe=xh(1)*rhh
         rhe=((1.-abu_he)*xh(1)+abu_he*xhe(1)+
     &        2.0*abu_he*xhe(2)+abu_c)*rhh
         ahe0=phihe(0)+rhe*acolhe0
         ahe1=phihe(1)+rhe*acolhe1
         rbhe0=rhe*breche0
         rbhe1=rhe*breche1
         hd0d1=(ahe0+ahe1+rbhe0+rbhe1)*0.5
         d01=rbhe0*rbhe1+ahe0*(ahe1+rbhe1)
         delthe=hd0d1+dsqrt(hd0d1*hd0d1-d01)
         eqxhe1=ahe0*rbhe1/d01
         eqxhe2=eqxhe1*ahe1/rbhe1
c**   eqxhe2=ahe0*ahe1/d01
         xhe(1)=(xhe1old-eqxhe1)*dexp(-delthe*dt)+eqxhe1
         xhe(2)=(xhe2old-eqxhe2)*dexp(-delthe*dt)+eqxhe2
         rhe=((1.-abu_he)*xh(1)+abu_he*xhe(1)
     &        +2.0*abu_he*xhe(2)+abu_c)*rhh
      enddo
c     
c-----------------------------------------------------------------------
c     determine neutral densities (take care of precision fluctuations)
c-----------------------------------------------------------------------
c
      if (xh(0).lt.0.0.and.abs(xh(0)).lt.1.0e-8)
     &         xh(0)=0.0
      xhe(0)=1.0-(xhe(1)+xhe(2))
      if (xhe(0).lt.0.0.and.abs(xhe(0)).lt.1.0e-8)
     &         xhe(0)=0.0
c
c-----------------------------------------------------------------------
c     find the temperature and density positions in the tables
c   ORIGINAL COOLING STUFF
c-----------------------------------------------------------------------
c
C
      temp00=dmax1(temp0,200.)
      tpos=dmin1(41.0,dmax1(1.0,(dlog10(temp00)-2.0)*10.0))
      itpos=int(tpos)
      dtpos=tpos-float(itpos)
      itpos1=min(41,itpos+1)
      dpos=dmin1(13.0,dmax1(1.0,dlog10(rhe)*2.0))
      idpos=int(dpos)
      ddpos=dpos-float(idpos)
      idpos1=min(13,idpos+1)
C
c-----------------------------------------------------------------------
c     find the temperature and density positions in the tables
C     COOLING STUFF Grazyna
c-----------------------------------------------------------------------
c
      denmin=10.
      aldenm=dlog10(denmin)   
      alogdeld=1.
      tempmin=100.             
      altempm=dlog10(tempmin)
      alogdelt=0.1                                                                                  
C     Temperature position                                                                          
      temp00gz=dmax1(temp0,200.)           
      tposgz=dmin1(1.e0*ntetable,dmax1(1.0,
     &     (dlog10(temp00gz)-altempm)/alogdelt)+1.)                                                  
      itposgz=min(int(tposgz),ntetable)                                                             
      dtposgz=10**(alogdelt)                                                                        
      itpos1gz=min(ntetable,itposgz+1)                                                              
C     Density position                                                                              
      dposgz=dmin1(1.e0*ndetable,                                                                   
     &     dmax1(1.0,(dlog10(rhe)-aldenm)/alogdeld)+1.)                                              
      idposgz=min(int(dposgz),ndetable)                                                             
      ddposgz=10**alogdeld                                                                          
      idpos1gz=min(ndetable,idposgz+1)  
c
c-----------------------------------------------------------------------
c     find time-dependent ionization state of the other elements
c-----------------------------------------------------------------------
c
      call timeion(dt,temp00gz,rhe,rhh*xh(0),rhh*xh(1),xc,1,1,4)
      call timeion(dt,temp00gz,rhe,rhh*xh(0),rhh*xh(1),xn,2,0,5)
      call timeion(dt,temp00gz,rhe,rhh*xh(0),rhh*xh(1),xo,3,0,5)
      call timeion(dt,temp00gz,rhe,rhh*xh(0),rhh*xh(1),xne,4,0,5)
      call timeion(dt,temp00gz,rhe,rhh*xh(0),rhh*xh(1),xs,5,1,5)
c
c-----------------------------------------------------------------------
c     find the value of the due to collisional ionizations of hydrogen
c-----------------------------------------------------------------------
c
      fracol=(ah0-phih)/ah0
      sph=rhh*((xh(1)-xh1old)*fracol/dt-xh(1)*rhe*brech0*(1.-fracol))
     * *hionen
c
c-----------------------------------------------------------------------
c     find the value of the total radiative cooling.
c-----------------------------------------------------------------------
c
C
      emin=coolin(rhh,rhe,xh,xhe,xc,xn,xo,xs,xne,sph,temp00gz)
C
c
c-----------------------------------------------------------------------
c     find the value of the photoionization heating
c-----------------------------------------------------------------------
c
c      
      eplus=rhh*(xh(0)*(1.-abu_he)*hvphih+abu_he*(xhe(0)*hvphhe(0)+
     &     xhe(1)*hvphhe(1)))

      return
      end
c
c=======================================================================
c
c
      subroutine timeion (dt,temp0,rhoe,rh0,rh1,x,isp,lmin,lmax)
c
c
c=======================================================================
c
c           solves for the time dependent ionization structure
c  
c                           garrelt mellema
c 
c                            06-mar-1997
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c                          introduction
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c-----------------------------------------------------------------------
c     goal
c-----------------------------------------------------------------------
c
c     this subroutine calculates the ionization state of a species
c     in a time dependent way by iterating a simple two-level
c     calculation.
c
c      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c          formats, parameters, declarations, common blocks, data
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c-----------------------------------------------------------------------
c     implicits
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     the formats
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     parameters
c       levmax: maximum number of ionization states (0 to levmax)
c       numspec: number of metals
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      parameter (levmax=5,numspec=5)
c
c-----------------------------------------------------------------------
c     declarations
c-----------------------------------------------------------------------
c
      dimension x(0:levmax),s(0:levmax+1),a(0:levmax),
     &     delta(0:levmax),equi(levmax)
c
c-----------------------------------------------------------------------
c     common blocks
c-----------------------------------------------------------------------
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /atomtab/ atomic data for metals
c
      common /atomtab/ rec(41,numspec,0:levmax),
     &     col(41,numspec,0:levmax),h0chex(41,numspec,0:levmax),
     &     h1chex(41,numspec,0:levmax),
     &     phot(0:100,numspec,0:levmax,3)
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /tabpos/ temperature and density positions in the tables
c
      common /tabpos/ tpos,dtpos,dpos,ddpos,
     &     itpos,itpos1,idpos,idpos1
      common /tabposgz/ tposgz,dtposgz,dposgz,ddposgz,itposgz
     &     ,itpos1gz,idposgz,idpos1gz

c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /photrates/ local photo-ionization rates
c     
      common /photrates/ phih,phihe(0:1),aphi(numspec,0:levmax),
     &     hvphih,hvphhe(0:1)
     
      save     

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     time dependent ionization calculation
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c-----------------------------------------------------------------------
c     set the 'boundary values' for the ionization changes
c-----------------------------------------------------------------------
c
      s(lmin)=0.0
      s(lmax+1)=0.0
c
c-----------------------------------------------------------------------
c     calculate the ionization and recombination rates
c     acol  - collisional ionization
c     achex - charge exchange ionization
c     aphi  - photo ionization (calculated in photoion)
c     a     - total ionization rate
c-----------------------------------------------------------------------
c
      do l=lmin+1,lmax
C
         acol=col(itpos,isp,l-1)+(col(itpos1,isp,l-1)-
     &        col(itpos,isp,l-1))*dtpos
         achex=h1chex(itpos,isp,l-1)+(h1chex(itpos1,isp,l-1)-
C
     &       h1chex(itpos,isp,l-1))*dtpos
         a(l-1)=rhoe*acol+rh1*achex+aphi(isp,l-1)
c
c-----------------------------------------------------------------------
c     calculate the recombination rates
c     arec  - radiative and dielectronic recombination
c     achex - charge exchange recombination
c-----------------------------------------------------------------------
c
C
         brec=rec(itpos,isp,l)+(rec(itpos1,isp,l)-rec(itpos,isp,l))*
     &        dtpos 
         bchex=h0chex(itpos,isp,l)+(h0chex(itpos1,isp,l)-
     &       h0chex(itpos,isp,l))*dtpos
C
c
c-----------------------------------------------------------------------
c     calculate the sum of the total ionization and total 
c     recombination rate. This are used below in the iteration.
c-----------------------------------------------------------------------
c
         delta(l)=a(l-1)+rhoe*brec+rh0*bchex
      enddo
c
c-----------------------------------------------------------------------
c     iterate the two-level solutions to obtain convergence
c-----------------------------------------------------------------------
c
      nitt=lmax*2
      dtt=dt/float(nitt)
      do nt=1,nitt
         do l=lmin+1,lmax
            equi(l)=a(l-1)*(x(l-1)+x(l))/delta(l)
            s(l)=(equi(l)-x(l))*(1.0-dexp(-dtt*delta(l)))
         enddo
         do l=lmin,lmax
            x(l)=x(l)+(s(l)-s(l+1))
         enddo
      enddo

      return
      end

c=======================================================================
c
c
      double precision function coolin (rhh,rhe,xh,xhe,xc,xn,xo,xs,xne
     & ,sph,temp0)
c
c
c=======================================================================
c
c                  calculates radiative cooling terms
c  
c                           garrelt mellema
c 
c                            06-Mar-1997
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c                          introduction
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c-----------------------------------------------------------------------
c     goal
c-----------------------------------------------------------------------
c
c      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c          formats, parameters, declarations, common blocks, data
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c-----------------------------------------------------------------------
c     implicits
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c     the formats
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     parameters
c       levmax: maximum number of ionization states (0 to levmax)
c       numspec: number of metals
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      parameter(numspec=5,levmax=5)
c
c-----------------------------------------------------------------------
c     declarations
c-----------------------------------------------------------------------
c
      dimension xh(0:1),xc(0:levmax),xn(0:levmax),xo(0:levmax),
     &     xne(0:levmax),xhe(0:2),xs(0:levmax)
c      ,xs(0:levmax)
      dimension emh(0:1),emc(1:3),emn(0:4),emo(0:5),emne(1:5),emhe(0:2)
      dimension ems(1:5)
c
c-----------------------------------------------------------------------
c     common blocks
c-----------------------------------------------------------------------
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /abunda/ transports heavy element abundances 
c
      common /abunda/ abu_he,abu_c,abu_n,abu_o,abu_ne,abu_s
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /cooltab/ cooling tables
c
c      common /cooltab/ coolh(41,0:1),coolc(13,41,0:levmax),
c     &     cooln(13,41,0:levmax),coolo(13,41,0:levmax),
c     &     coolne(13,41,0:levmax)
      parameter(ntetable=49,ndetable=8,ncoolgz=ntetable*ndetable
     &,levmax2=19)
C
      common /cooltabgz/coolhgz(ndetable,ntetable,0:levmax2)
     &     ,coolhegz(ndetable,ntetable,0:levmax2)
     &     ,coolcgz(ndetable,ntetable,0:levmax2)
     &     ,coolangz(ndetable,ntetable,0:levmax2)
     &     ,coologz(ndetable,ntetable,0:levmax2)
     &     ,coolanegz(ndetable,ntetable,0:levmax2)
     &     ,coolamggz(ndetable,ntetable,0:levmax2)
     &     ,coolsigz(ndetable,ntetable,0:levmax2)
     &     ,coolsgz(ndetable,ntetable,0:levmax2)
     &     ,coolclgz(ndetable,ntetable,0:levmax2)
     &     ,coolargz(ndetable,ntetable,0:levmax2) 
c
      common/tmaxim/tmaximo 
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /tabpos/ temperature and density positions in the tables
c
      common /tabpos/ tpos,dtpos,dpos,ddpos,itpos,itpos1,idpos,idpos1
      common /tabposgz/ tposgz,dtposgz,dposgz,ddposgz,itposgz
     &     ,itpos1gz,idposgz,idpos1gz
            
      save
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c               calculate cooling
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     Hydrogen
c
      do l=0,1
         emh(l)=rhe*rhh*xh(l)
     &       *((1.0-ddposgz)*(1.0-dtposgz)*coolhgz(idposgz,itposgz,l)+
     &        ddposgz*(1.0-dtposgz)*coolhgz(idposgz+1,itposgz,l)+
     &        ddposgz*dtposgz*coolhgz(idposgz+1,itpos1gz,l)+
     &        (1.0-ddposgz)*dtposgz*coolhgz(idposgz,itpos1gz,l))
      enddo
c
c     Helio
      do l=0,2
         emhe(l)=rhe*abu_he*xhe(l)
     &       *((1.0-ddposgz)*(1.0-dtposgz)*coolhegz(idposgz,itposgz,l)+
     &        ddposgz*(1.0-dtposgz)*coolhegz(idposgz+1,itposgz,l)+
     &        ddposgz*dtposgz*coolhegz(idposgz+1,itpos1gz,l)+
     &        (1.0-ddposgz)*dtposgz*coolhegz(idposgz,itpos1gz,l))
      enddo
c      
c     Carbon

      do l=1,3
         emc(l)=rhe*abu_c*rhh*xc(l)*
     &        ((1.0-ddposgz)*(1.0-dtposgz)*coolcgz(idposgz,itposgz,l)+
     &        ddposgz*(1.0-dtposgz)*coolcgz(idposgz+1,itposgz,l)+
     &        ddposgz*dtposgz*coolcgz(idposgz+1,itpos1gz,l)+
     &        (1.0-ddposgz)*dtposgz*coolcgz(idposgz,itpos1gz,l))
      enddo
c
c     Nitrogen
c
      do l=0,4
         emn(l)=rhe*abu_n*rhh*xn(l)*
     &        ((1.0-ddposgz)*(1.0-dtposgz)*coolangz(idposgz,itposgz,l)+
     &        ddposgz*(1.0-dtposgz)*coolangz(idposgz+1,itposgz,l)+
     &        ddposgz*dtposgz*coolangz(idposgz+1,itpos1gz,l)+
     &        (1.0-ddposgz)*dtposgz*coolangz(idposgz,itpos1gz,l))
      enddo
c
c     Oxygen
C
      do l=0,5
         emo(l)=rhe*abu_o*rhh*xo(l)*
     &        ((1.0-ddposgz)*(1.0-dtposgz)*coologz(idposgz,itposgz,l)+
     &        ddposgz*(1.0-dtposgz)*coologz(idposgz+1,itposgz,l)+
     &        ddposgz*dtposgz*coologz(idposgz+1,itpos1gz,l)+
     &        (1.0-ddposgz)*dtposgz*coologz(idposgz,itpos1gz,l))
      enddo
c
c     Neon
c
      do l=1,5
         emne(l)=rhe*abu_ne*rhh*xne(l)*
     &        ((1.0-ddposgz)*(1.0-dtposgz)*coolanegz(idposgz,itposgz,l)+
     &        ddposgz*(1.0-dtposgz)*coolanegz(idposgz+1,itposgz,l)+
     &        ddposgz*dtposgz*coolanegz(idposgz+1,itpos1gz,l)+
     &        (1.0-ddposgz)*dtposgz*coolanegz(idposgz,itpos1gz,l))
      enddo
c
      do l=1,5
         ems(l)=rhe*abu_s*rhh*xs(l)*
     &        ((1.0-ddposgz)*(1.0-dtposgz)*coolsgz(idposgz,itposgz,l)+
     &        ddposgz*(1.0-dtposgz)*coolsgz(idposgz+1,itposgz,l)+
     &        ddposgz*dtposgz*coolsgz(idposgz+1,itpos1gz,l)+
     &        (1.0-ddposgz)*dtposgz*coolsgz(idposgz,itpos1gz,l))
      enddo
c
c
C      Total cooling
c
      coolin=emh(0)+emh(1)+emhe(0)+emhe(1)+emhe(2)+emc(1)+emc(2)+emc(3)+
     &     emn(0)+emn(1)+emn(2)+emn(3)+emn(4)+
     &     emo(0)+emo(1)+emo(2)+emo(3)+emo(4)+emo(5)+
     &     emne(1)+emne(2)+emne(3)+emne(4)+emne(5)+
     &     ems(1)+ems(2)+ems(3)+ems(4)+ems(5)
c
c COOLING OLD STUFF
c     Hydrogen
c
c      emh(0)=rhe*rhh*xh(0)*(coolh(itpos,0)+
c     &     (coolh(itpos1,0)-coolh(itpos,0))*dtpos)
c      emh(1)=rhe*rhh*xh(1)*(coolh(itpos,1)+
c     &     (coolh(itpos1,1)-coolh(itpos,1))*dtpos)+sph
c
c     Carbon
c
c      do l=1,3
c         emc(l)=rhe*abu_c*rhh*xc(l)*
c     &        ((1.0-ddpos)*(1.0-dtpos)*coolc(idpos,itpos,l)+
c     &        ddpos*(1.0-dtpos)*coolc(idpos+1,itpos,l)+
c     &        ddpos*dtpos*coolc(idpos+1,itpos1,l)+
c     &        (1.0-ddpos)*dtpos*coolc(idpos,itpos1,l))
c      enddo
c
c     Nitrogen
c
c      do l=0,4
c         emn(l)=rhe*abu_n*rhh*xn(l)*
c     &        ((1.0-ddpos)*(1.0-dtpos)*cooln(idpos,itpos,l)+
c     &        ddpos*(1.0-dtpos)*cooln(idpos+1,itpos,l)+
c     &        ddpos*dtpos*cooln(idpos+1,itpos1,l)+
c     &        (1.0-ddpos)*dtpos*cooln(idpos,itpos1,l))
c      enddo
c
c     Oxygen
c
c      do l=0,4
c         emo(l)=rhe*abu_o*rhh*xo(l)*
c     &        ((1.0-ddpos)*(1.0-dtpos)*coolo(idpos,itpos,l)+
c     &        ddpos*(1.0-dtpos)*coolo(idpos+1,itpos,l)+
c     &        ddpos*dtpos*coolo(idpos+1,itpos1,l)+
c     &        (1.0-ddpos)*dtpos*coolo(idpos,itpos1,l))
c      enddo
c
c     Neon
c
c      do l=2,4
c         emne(l)=rhe*abu_ne*rhh*xne(l)*
c     &        ((1.0-ddpos)*(1.0-dtpos)*coolne(idpos,itpos,l)+
c     &        ddpos*(1.0-dtpos)*coolne(idpos+1,itpos,l)+
c     &        ddpos*dtpos*coolne(idpos+1,itpos1,l)+
c     &        (1.0-ddpos)*dtpos*coolne(idpos,itpos1,l))
c      enddo
c
c     Total cooling
c
c      coolin=emh(0)+emh(1)+emc(1)+emc(2)+emc(3)+
c**     &     emn(0)+emn(1)+emn(2)+emn(3)+emn(4)+
c     &     emo(0)+emo(1)+emo(2)+emo(4)+
c     &     emne(2)+emne(3)+emne(4)

      return
      end

c=======================================================================
c
c
                           subroutine inirad
c
c
c=======================================================================
c            
c            initializes constants and tables for radiation processes
c                    (heating, cooling and ionization)
c
c                           garrelt mellema
c
c                             06-Mar-1997
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c                          introduction
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c-----------------------------------------------------------------------
c     goal
c-----------------------------------------------------------------------
c
c     This subroutine sets constants and tables having to do 
c     with radiation effects. These are recombination and collisional
c     ionization rates, cooling rates and heating rates. Some of the 
c     data is read out of data files.
c 
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c          formats, parameters, declarations, common blocks, data
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c-----------------------------------------------------------------------
c     implicits
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     parameters
c       levmax: maximum number of ionization states (0 to levmax)
c       numspec: number of metals
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      parameter(numspec=5,levmax=5)
c
c-----------------------------------------------------------------------
c     declarations
c-----------------------------------------------------------------------
c
      integer nspec(numspec)
      real*8 m_p
      dimension ethe(0:1),xihe(0:1),fhe(0:1)
c
c-----------------------------------------------------------------------
c     common blocks
c-----------------------------------------------------------------------
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /abunda/ transports heavy element abundances 
c
      common /abunda/ abu_he,abu_c,abu_n,abu_o,abu_ne,abu_s
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /hydrogen/ hydrogen atomic data
c
      common /hydrogen/ bh00,albpow,colh0,temph0,hionen,
     &     hphot(0:100,3),hheat(0:100,3)
c
c     /helium/ helium atomic data
c
      common /helium/ bhe00,bhe10,alcpow,colhe(0:1),temphe(0:1),
     &     hephot(0:100,0:1,3),heheat(0:100,0:1,3)
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /atomtab/ metals atomic data tables
c
      common /atomtab/ rec(41,numspec,0:levmax),
     &     col(41,numspec,0:levmax),h0chex(41,numspec,0:levmax),
     &     h1chex(41,numspec,0:levmax),
     &     phot(0:100,numspec,0:levmax,3)
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /cooltab/ cooling tables
c
      common /cooltab/ coolh(41,0:1),coolc(13,41,0:levmax),
     &     cooln(13,41,0:levmax),coolo(13,41,0:levmax),
     &     coolne(13,41,0:levmax)
      parameter(ntetable=49,ndetable=8,ncoolgz=ntetable*ndetable
     &     ,levmax2=19)
C
      common /cooltabgz/coolhgz(ndetable,ntetable,0:levmax2)
     &     ,coolhegz(ndetable,ntetable,0:levmax2)
     &     ,coolcgz(ndetable,ntetable,0:levmax2)
     &     ,coolangz(ndetable,ntetable,0:levmax2)
     &     ,coologz(ndetable,ntetable,0:levmax2)
     &     ,coolanegz(ndetable,ntetable,0:levmax2)
     &     ,coolamggz(ndetable,ntetable,0:levmax2)
     &     ,coolsigz(ndetable,ntetable,0:levmax2)
     &     ,coolsgz(ndetable,ntetable,0:levmax2)
     &     ,coolclgz(ndetable,ntetable,0:levmax2)
     &     ,coolargz(ndetable,ntetable,0:levmax2)
C
      Dimension tgz(ncoolgz),degz(ncoolgz)
      Dimension hgz(0:1,ncoolgz),hegz(0:2,ncoolgz),cgz(0:6,ncoolgz)
      Dimension angz(0:7,ncoolgz),ogz(0:8,ncoolgz),anegz(0:10,ncoolgz)
      Dimension amggz(0:12,ncoolgz),sigz(0:14,ncoolgz)
      Dimension sgz(0:16,ncoolgz),clgz(0:17,ncoolgz),argz(0:18,ncoolgz)
      integer num,i
      character*100 fila1,fila2,fila3,fila4,fila5
      character*50 temptgz
c      Save tgz,degz,hgz,hegz,cgz,angz,ogz,anegz,amggz,sigz,sggz,clgz
c     &     , argz
      Common/tmaxim/tmaximo 
c      save rec,col,h0chex,h1chex,phot
c     
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /optdep/ constants needed to calculate optical depths
c
      common /optdep/ sigh,sighe0,sighe1,tf2h,tf3h,tf3he0
c
      common/hhphot/hhphot
      COMMON/ABUND/XH,XHE,XC,XN,XO,XS,XAr,XNe
      
      save      
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c            set up the radiative constants and tables 
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c-----------------------------------------------------------------------
c     set heavy element abundances
c     from Osterbrock (1989), Table 5.13
c-----------------------------------------------------------------------

C
      abu_he=XHE
      abu_c=XC
      abu_n=XN
      abu_o=XO
      abu_ne=XNe
      abu_s=XS
C
c**      abu_c=3.3e-4
c**      abu_n=0.0
c**      abu_o=6.6e-4
c**      abu_ne=8.3e-5
c**      abu_s=1.6e-5
c
c-----------------------------------------------------------------------
c     set the values of physical constants
c     m_p    - proton mass (in g)
c     ev2k   - conversion factor between evs and kelvins
c     ev2erg   - conversion factor between evs and ergs
c-----------------------------------------------------------------------
c
      m_p=1.672661e-24
      ev2k=1.0/8.617e-05
      ev2erg=1.602e-12
      ev2fr=0.241838d15
c
c-----------------------------------------------------------------------
c     initialize hydrogen recombination coefficient (case B)
c     atomic data from Osterbrock (1989), tables 2.1
c-----------------------------------------------------------------------
c
c      albpow=-0.82
      bh00 =2.59e-13
C      bh00 =4.18e-13
      bh00=bh00*(1.0e-4)**albpow
c
c----------------------------------------------------------------------      
c-----------------------------------------------------------------------
c     initialize helium recombination coefficients 
c     (case A for He0 and He1)
c     atomic data from Osterbrock (1989), tables 2.4 and 2.8
c-----------------------------------------------------------------------
c
      alcpow=-0.672
      bhe00 =2.73e-13+1.59e-13
      bhe00=bhe00*(1.0e-4)**alcpow
      bhe10 =5.91e-12
      bhe10=bhe10*(1.0e-4)**albpow
c
c-----------------------------------------------------------------------
c     set the collisional ionization data for hydrogen
c     (from Cox 1970)
c-----------------------------------------------------------------------
c
      eth0=13.598
      xih0=1.0
      fh0=0.83
      colh0=1.3e-8*fh0*xih0/(eth0*eth0)
      temph0=eth0*ev2k
      hionen=eth0*ev2erg
c
c-----------------------------------------------------------------------
c     set the collisional ionization data for helium
c     (from Cox 1970)
c-----------------------------------------------------------------------
c
      ethe(0)=24.587
      xihe(0)=2.0
      fhe(0)=0.63
      ethe(1)=54.416
      xihe(1)=1.0
      fhe(1)=1.30
      do k=0,1
         colhe(k)=1.3d-8*fhe(k)*xihe(k)/(ethe(k)*ethe(k))
         temphe(k)=ethe(k)*ev2k
      enddo
c
c-----------------------------------------------------------------------
c     photoionization constants (for h, he0, he1)
c       sig - cross section at threshold
c       frt - threshold frequency
c       s   - power of frequency dependence of cross section
c       tf  - correction factor to add H, He0, He1 optical depths
c-----------------------------------------------------------------------
c
      sigh=6.30e-18
      frth0=ev2fr*eth0
      sh0=2.8
      sighe0=7.83e-18
      frthe0=ev2fr*ethe(0)
      she0=1.7
      sighe1=1.58e-18
      frthe1=ev2fr*ethe(1)
      she1=2.8      
      tf2h=(0.63*frth0/frthe0)**she0
      tf3h=(frth0/frthe1)**she1
      tf3he0=(1.51*frthe0/frthe1)**she1
c
c-----------------------------------------------------------------------
c     read in atomic data for other elements (C,N,O,Ne,S)
c       rec - total recombination rate
c       col - collisional ionization rate
c       h0chex - charge exchange with h0 (recombination)
c       h1chex - charge exchange with h1 (ionization)
c-----------------------------------------------------------------------
c
      open(unit=21,file='coeff.tab',status='old')
      do isp=1,numspec
         do l=1,levmax
            read(21,*) nspec(isp),nion
            do i=1,41
               read(21,*) temp,rec(i,isp,nion),col(i,isp,nion-1),
     &              h0chex(i,isp,nion),h1chex(i,isp,nion-1)
            enddo
c
         enddo
      enddo
      close(21)
c
c-----------------------------------------------------------------------
c     Convert atomic data from log to linear 
c-----------------------------------------------------------------------
c
      do l=1,levmax
         do isp=1,numspec
            do it=1,41
               rec(it,isp,l)=10.0**(rec(it,isp,l))
               col(it,isp,l-1)=10.0**(col(it,isp,l-1))
               h0chex(it,isp,l)=10.0**(h0chex(it,isp,l))
               h1chex(it,isp,l-1)=10.0**(h1chex(it,isp,l-1))
            enddo
         enddo
      enddo
c**      write(*,*) rec(20,3,2),col(20,3,2),h0chex(20,3,2),h1chex(20,3,2)
c
c-----------------------------------------------------------------------
c     read in log(cooling) for cooling species 
c-----------------------------------------------------------------------
      open(unit=22, file='EMISLOSS.DAT')
      read(22,1010)fila1
      read(22,1010)fila2
      read(22,1010)fila3
      read(22,1010)fila4
      read(22,1010)fila5
C
      Do i=1,ncoolgz
         read(22,1010)temptgz
         read(temptgz(13:20),*)tgz(i)
         read(temptgz(27:36),*)degz(i)
         read(22,*)num,hgz(0,i),hgz(1,i)
         read(22,*)num,hegz(0,i),hegz(1,i),hegz(2,i)
         read(22,*)num,cgz(0,i),cgz(1,i),cgz(2,i),cgz(3,i)
     &        ,cgz(4,i),cgz(5,i),cgz(6,i)
         read(22,*)num,angz(0,i),angz(1,i),angz(2,i),angz(3,i)
     &        ,angz(4,i),angz(5,i),angz(6,i),angz(7,i)
         read(22,*)num,ogz(0,i),ogz(1,i),ogz(2,i),ogz(3,i),ogz(4,i)
     &        ,ogz(5,i),ogz(6,i),ogz(7,i),ogz(8,i)
         read(22,*)num,anegz(0,i),anegz(1,i),anegz(2,i),anegz(3,i)
     &        ,anegz(4,i),anegz(5,i),anegz(6,i),anegz(7,i),anegz(8,i)
     &        ,anegz(9,i),anegz(10,i)
         read(22,*)num,amggz(0,i),amggz(1,i),amggz(2,i),amggz(3,i)
     &        ,amggz(4,i),amggz(5,i),amggz(6,i),amggz(7,i),amggz(8,i)
     &        ,amggz(9,i),amggz(10,i),amggz(11,i),amggz(12,i)
     &        
         read(22,*)num,sigz(0,i),sigz(1,i),sigz(2,i),sigz(3,i),sigz(4,i)
     &        ,sigz(5,i),sigz(6,i),sigz(7,i),sigz(8,i),sigz(9,i)
     &        ,sigz(10,i),sigz(11,i),sigz(12,i),sigz(13,i),sigz(14,i)
         read(22,*)num,sgz(0,i),sgz(1,i),sgz(2,i),sgz(3,i),sgz(4,i)
     &        ,sgz(5,i),sgz(6,i),sgz(7,i),sgz(8,i),sgz(9,i),sgz(10,i)
     &        ,sgz(11,i),sgz(12,i),sgz(13,i),sgz(14,i),sgz(15,i)
     &        ,sgz(16,i)
         read(22,*)num,clgz(0,i),clgz(1,i),clgz(2,i),clgz(3,i),clgz(4,i)
     &        ,clgz(5,i),clgz(6,i),clgz(7,i),clgz(8,i),clgz(9,i)
     &        ,clgz(10,i),clgz(11,i),clgz(12,i),clgz(13,i),clgz(14,i)
     &        ,clgz(15,i),clgz(16,i),clgz(17,i)
         read(22,*)num,argz(0,i),argz(1,i),argz(2,i),argz(3,i)
     &        ,argz(4,i),argz(5,i)
     &        ,argz(6,i),argz(7,i),argz(8,i),argz(9,i),argz(10,i)
     &        ,argz(11,i),argz(12,i),argz(13,i),argz(14,i),argz(15,i)
     &        ,argz(16,i),argz(17,i),argz(18,i)
C       
      enddo
      close(22)
c
      tmaximo=tgz(ncoolgz) 
c
      ii=1
      Do ittt=1,ntetable
         do iddd=1,ndetable
            coolhgz(iddd,ittt,0:1)=hgz(0:1,ii)
            coolhegz(iddd,ittt,0:2)=hegz(0:2,ii)
            coolcgz(iddd,ittt,0:6)=cgz(0:6,ii)
            coolangz(iddd,ittt,0:7)=angz(0:7,ii)
            coologz(iddd,ittt,0:8)=ogz(0:8,ii)
            coolanegz(iddd,ittt,0:10)=anegz(0:10,ii)
            coolamggz(iddd,ittt,0:12)=amggz(0:12,ii)
            coolsigz(iddd,ittt,0:14)=sigz(0:14,ii)
            coolsgz(iddd,ittt,0:16)=sgz(0:16,ii)
            coolclgz(iddd,ittt,0:17)=clgz(0:17,ii)
            coolargz(iddd,ittt,0:18)=argz(0:18,ii)
            ii=1+ii
         enddo
      enddo
     
 1010 Format(A)


c      open(unit=22,file='hcool.tab',status='old')
c      open(unit=23,file='cool.tab',status='old')
c
c     H
c
c      do l=1,2
c         read(22,*) nsp,nion,nchck
c         do it=1,41
c            read(22,*) temp,coolh(it,nion-1)
c         enddo
c         write(6,*) ' cool (H) : ',coolh(1,nion-1)
c      enddo
c
c     C
c
c      do l=1,4
c         read(23,*) nsp,nion,nchck
c         if (nchck.eq.1) then
c            do it=1,41
c               read(23,*) temp,(coolc(id,it,nion-1),id=1,13)
c            enddo
c         else
c            do it=1,41
c               read(23,*) temp,coolc(1,it,nion-1)
c               do id=2,13
c                  coolc(id,it,nion-1)=coolc(1,it,nion-1)
c               enddo
c            enddo
c         endif
c         write(6,*) ' cool (C) : ',coolc(1,1,nion-1)
c      enddo
c     
c     N
c
c      do l=1,5
c         read(23,*) nsp,nion,nchck
c        if (nchck.eq.1) then
c            do it=1,41
c               read(23,*) temp,(cooln(id,it,nion-1),id=1,13)
c            enddo
c         else
c            do it=1,41
c               read(23,*) temp,cooln(1,it,nion-1)
c               do id=2,13
c                  cooln(id,it,nion-1)=cooln(1,it,nion-1)
c               enddo
c           enddo
c         endif
c         write(6,*) ' cool (N) : ',cooln(1,1,nion-1)
c      enddo
c     
c     O
c
c      do l=1,5
c         read(23,*) nsp,nion,nchck
c         if (nchck.eq.1) then
c            do it=1,41
c               read(23,*) temp,(coolo(id,it,nion-1),id=1,13)
c            enddo
c         else
c            do it=1,41
c               read(23,*) temp,coolo(1,it,nion-1)
c               do id=2,13
c                  coolo(id,it,nion-1)=coolo(1,it,nion-1)
c               enddo
c            enddo
c         endif
c         write(6,*) ' cool (O) : ',coolo(1,1,nion-1)
c      enddo
c
c
c     Ne
c
c      do l=3,5
c         read(23,*) nsp,nion,nchck
c         do it=1,41
c            read(23,*) temp,(coolne(id,it,nion-1),id=1,13)
c         enddo
c      enddo

c      close(22)
c      close(23)
c
c-----------------------------------------------------------------------
c     Convert cooling to from log to linear 
c-----------------------------------------------------------------------
c
c      do it=1,41
c         do l=0,1
c            coolh(it,l)=10.0**(coolh(it,l))
c         enddo
c      enddo
c      do it=1,41
c         do id=1,13
c            do l=0,3
c               coolc(id,it,l)=10.0**(coolc(id,it,l))
c            enddo
c            do l=0,4
c               cooln(id,it,l)=10.0**(cooln(id,it,l))
c            enddo
c            do l=0,4
c               coolo(id,it,l)=10.0**(coolo(id,it,l))
c            enddo
c            do l=2,4
c               coolne(id,it,l)=10.0**(coolne(id,it,l))
c            enddo
c         enddo
c      enddo
      RETURN
      end
c
      SUBROUTINE GENPHOTO(RSTAR,TEFF,LSTAR)
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     parameters
c       levmax: maximum number of ionization states (0 to levmax)
c       numspec: number of metals
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      parameter(numspec=5,levmax=5)
c
c-----------------------------------------------------------------------
c     declarations
c-----------------------------------------------------------------------
c
      integer nspec(numspec)
      real*8 m_p
      real*8 lstar
      dimension ethe(0:1),xihe(0:1),fhe(0:1)
c
c-----------------------------------------------------------------------
c     common blocks
c-----------------------------------------------------------------------
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /abunda/ transports heavy element abundances 
c
      common /abunda/ abu_he,abu_c,abu_n,abu_o,abu_ne,abu_s
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /hydrogen/ hydrogen atomic data
c
      common /hydrogen/ bh00,albpow,colh0,temph0,hionen,
     &     hphot(0:100,3),hheat(0:100,3)
c
c     /helium/ helium atomic data
c
      common /helium/ bhe00,bhe10,alcpow,colhe(0:1),temphe(0:1),
     &     hephot(0:100,0:1,3),heheat(0:100,0:1,3)
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /atomtab/ metals atomic data tables
c
      common /atomtab/ rec(41,numspec,0:levmax),
     &     col(41,numspec,0:levmax),h0chex(41,numspec,0:levmax),
     &     h1chex(41,numspec,0:levmax),
     &     phot(0:100,numspec,0:levmax,3)
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /cooltab/ cooling tables
c
c      common /cooltab/ coolh(41,0:1),coolc(13,41,0:levmax),
c     &     cooln(13,41,0:levmax),coolo(13,41,0:levmax),
c     &     coolne(13,41,0:levmax)
cC
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /optdep/ constants needed to calculate optical depths
c
      common /optdep/ sigh,sighe0,sighe1,tf2h,tf3h,tf3he0
c
      common/hhphot/hhphot
      
      save
c
      CALL PHOTOTABLES(RSTAR,TEFF,LSTAR)
c
c     H
c
c      hhphot=hphot(0,1)+hphot(0,2)+hphot(0,3)
c
c     He
c
c      do ion=0,1
c         do lfr=1,3
c            read(24,*) (hephot(n,ion,lfr),n=0,100)
c         enddo
c         write(6,*) ' He : ',hephot(0,ion,1)+hephot(0,ion,2)+
c     *    hephot(0,ion,3)
c      enddo
c      do ion=0,1
c         do lfr=1,3
c            read(24,*) (heheat(n,ion,lfr),n=0,100)
c         enddo
c      enddo
c
c     C
c
c      do ion=0,4
c         do lfr=1,3
c            read(24,*) (phot(n,1,ion,lfr),n=0,100)
c         enddo
c         write(6,*) ' C : ',phot(0,1,ion,1)+phot(0,1,ion,2)+
c     *    phot(0,1,ion,3)
c      enddo
c
c     N
c
c     do ion=0,4
c     do lfr=1,3
c           read(24,*) (phot(n,2,ion,lfr),n=0,100)
c        enddo
c         write(6,*) ' N : ',phot(0,2,ion,1)+phot(0,2,ion,2)+
c     *    phot(0,2,ion,3)
c      enddo
c     
c     O
c
c      do ion=0,4
c         do lfr=1,3
c            read(24,*) (phot(n,3,ion,lfr),n=0,100)
c         enddo
c         write(6,*) ' O : ',phot(0,3,ion,1)+phot(0,3,ion,2)+
c     *    phot(0,3,ion,3)
c      enddo
c
c     Ne
c
c      do ion=0,4
c         do lfr=1,3
c            read(24,*) (phot(n,4,ion,lfr),n=0,100)
c         enddo
c         write(6,*) ' Ne : ',phot(0,4,ion,1)+phot(0,4,ion,2)+
c     *    phot(0,4,ion,3)
c      enddo
c
c     S
c
c      do ion=0,4
c         do lfr=1,3
c            read(24,*) (phot(n,5,ion,lfr),n=0,100)
c         enddo
c         write(6,*) ' S : ',phot(0,5,ion,1)+phot(0,5,ion,2)+
c     *    phot(0,5,ion,3)
c      enddo

      close(24)

      return
      end
c
c
c
c=======================================================================
c
c
      subroutine phototables(rstar,teff,lstar)
c
c
c=======================================================================
c
c                Calculation of photo-ionization tables
c
c                     to be used by other programs
c
c                        Author: Garrelt Mellema
c
c                         Update: 03-Mar-1997
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     Introduction
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     This program calculates the photoionization tables to be used
c     with a ionization and heating calculation.
c     The tables are organised as follows:
c
c      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c          formats, parameters, declarations, common blocks, data
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c-----------------------------------------------------------------------
c     implicits
c-----------------------------------------------------------------------
c
c      implicit doubleprecision (a-h,o-z), integer*4 (i-n)
c
c-----------------------------------------------------------------------
c     the formats
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c     parameters
c
c     'gridsize' contains the values of the size of the grid: meshr, meshrp.
c     lint - the number of integration points in one of the three frequency
c           interval.
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      parameter(pi=3.141592654,lint=32)
      parameter(numspec=5,levmax=5)
c
c-----------------------------------------------------------------------
c     declarations
c-----------------------------------------------------------------------
c
c      doubleprecision nphot,nephot,lstar
      real*8 nphot,nephot
      real*8 lstar,llstar
c
c-----------------------------------------------------------------------
c     common blocks
c-----------------------------------------------------------------------
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /scalin/ transmits the scaling parameters
c
      common /scalin/ scleng,scmass,sctime,scvelo,scdens,
     &                scmome,scener,sccool
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /radtab/ transmits radiative tables
c
c      common /radtab/ hphot(0:100,3),hheat(0:100,3),
c     &     hephot(0:100,0:1,3),heheat(0:100,0:1,3),
c     &     cphot(0:100,0:5,3),nphot(0:100,0:5,3),
c     &     ophot(0:100,0:5,3),nephot(0:100,0:5,3),
c     $     sphot(0:100,0:5,3)
      common /radtab/ ahphot(0:100,3),ahheat(0:100,3),
     &     ahephot(0:100,0:1,3),aheheat(0:100,0:1,3),
     &     cphot(0:100,0:5,3),nphot(0:100,0:5,3),
     &     ophot(0:100,0:5,3),nephot(0:100,0:5,3),
     $     sphot(0:100,0:5,3)
      common /hydrogen/ bh00,albpow,colh0,temph0,hionen,
     &     hphot(0:100,3),hheat(0:100,3)
c
c     /helium/ helium atomic data
c
      common /helium/ bhe00,bhe10,alcpow,colhe(0:1),temphe(0:1),
     &     hephot(0:100,0:1,3),heheat(0:100,0:1,3)
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /atomtab/ metals atomic data tables
c
      common /atomtab/ rec(41,numspec,0:levmax),
     &     col(41,numspec,0:levmax),h0chex(41,numspec,0:levmax),
     &     h1chex(41,numspec,0:levmax),
     &     phot(0:100,numspec,0:levmax,3)
c
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /star/ transports time dependent stellar properties
c
      common /star/ teff1,rrstar,llstar
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c               calculate and write the tables
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      ON REAL*4 DIV 0    ABORT
c      ON REAL*4 OVERFLOW ABORT
c      ON REAL*4 UNDERFLOW IGNORE
c      ON REAL*4 INEXACT IGNORE
c      ON REAL*4 ILLEGAL  ABORT

c      ON REAL*8 DIV 0    ABORT
c      ON REAL*8 OVERFLOW ABORT
c      ON REAL*8 UNDERFLOW IGNORE
c      ON REAL*8 INEXACT IGNORE
c      ON REAL*8 ILLEGAL  ABORT
c      ON INTEGER*4 DIV 0    ABORT
c      ON INTEGER*4 OVERFLOW ABORT
c
c-----------------------------------------------------------------------
c     scalings
c-----------------------------------------------------------------------
c

      save
      scleng=1.0
      scmass=1.0
      sctime=1.0
      scvelo=1.0
      scdens=1.0
      scmome=1.0
      scener=1.0
      sccoolc=1.0
c
         totflux=5.670d-05*teff**4
c         rstar=dsqrt(lstar*3.862d33/(4.0d0*pi*totflux))
         rstar=dsqrt(lstar/(4.0d0*pi*totflux))
c       write(*,*) 'rstar: ',rstar,teff,lstar
      teff1=teff
      llstar=lstar
      rrstar=rstar
c
c-----------------------------------------------------------------------
c     Initialize some variables
c-----------------------------------------------------------------------
c
      call inirad2
c
c-----------------------------------------------------------------------
c     Do the integrals
c-----------------------------------------------------------------------
c
      call bbintegr
c
c-----------------------------------------------------------------------
c     write the tables to files
c-----------------------------------------------------------------------
c
c copiamos tasas de photoionizacion de CII-CVI a su lugar
      do i=0,4
         do n=0,100
            do l=1,3
               phot(n,1,i,l)=cphot(n,i,l)
c               write(*,*) 'c=',n,i,l,cphot(n,i,l),hphot(n,i)
            enddo
         enddo
      enddo
c tasas de N
      do i=0,4
         do n=0,100
            do l=1,3
               phot(n,2,i,l)=nphot(n,i,l)
C               write(*,*) 'n=',n,i,l,nphot(n,i,l),ahephot(n,i,l)
            enddo
         enddo
      enddo
c tasas de O
      do i=0,4
         do n=0,100
            do l=1,3
               phot(n,3,i,l)=ophot(n,i,l)
            enddo
         enddo
      enddo
c tasas de Ne
      do i=0,4
         do n=0,100
            do l=1,3
               phot(n,4,i,l)=nephot(n,i,l)
            enddo
         enddo
      enddo
c tasas de S
      do i=0,4
         do n=0,100
            do l=1,3
               phot(n,5,i,l)=sphot(n,i,l)
            enddo
         enddo
      enddo
      do i=1,3
         do n=0,100
            hphot(n,i)=ahphot(n,i)
            hheat(n,i)=ahheat(n,i)
            do l=0,1
               hephot(n,l,i)=ahephot(n,l,i)
               heheat(n,l,i)=aheheat(n,l,i)
            end do
         end do
      end do
c-----------------------------------------------------------------------
c
      return
      end
c
c
c=======================================================================
c
c
                           subroutine bbintegr
c
c
c=======================================================================
c
c                 calculates photo ionization integrals
c  
c                         for use with phototables
c  
c                           garrelt mellema
c 
c                              3-Mar-1997
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c                          introduction
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c-----------------------------------------------------------------------
c     goals
c-----------------------------------------------------------------------
c
c     for new blackbody parameters (effective temperature and radius),
c     this subroutine calcultes the photo-ionisation tables (these tables
c     are used in ionic to calculate the ionisation and heating rates).
c     the photo-ionisation tables contain the value of the photo-ionisation
c     integrals for a range of optical depths at the threshold -by assuming
c     a fixed frequency dependence for the optical depth (a powerlaw) the
c     integration over frequency can be carried out-.
c     the photon absorbing species are h0, he0, and he1, so there are three
c     thresholds and thus three tables for every ionisation process. see
c     ionic to learn how _one_ optical depth is constructed at each 
c     threshold.
c     the photo-ionisation tables (hphot, hephot, nphot, and ophot) supply
c     the photo-ionisation rate, the tables hheat and heheat, the heating 
c     rate due to photo-ionisation.
c     each integral contains a functional part that is independent of
c     stellar properties and thus static. to save execution time this
c     part (the integral core), is alculated in inirad and stored as hintc,
c     he0intc, he1intc,nintc,ointc, for the photo-ionisation rate integrals,
c     and hhint, hhe0int, hhe1int for the heating integrals. the first index
c     of these cores indicates frequency according to:
c     fr(i)=frs + step * i.
c     frs is the start frequency (frsh, frshe, frsn, frso), step the stepsize
c     (steph, stephe, stepn, stepo).
c     the second index (n) denotes optical depth, the last index indicates
c     the frequency range (1,2,3). in ointc and nintc the third index 
c     determines the ionisation species.
c     the final tables have the same index pattern, but without any
c     frequency dependence.
c
c      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c          formats, parameters, declarations, common blocks, data
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c-----------------------------------------------------------------------
c     implicits
c-----------------------------------------------------------------------
c
c      implicit doubleprecision (a-h,o-z), integer*4 (i-n)
c
c-----------------------------------------------------------------------
c     the formats
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c     parameters
c
c     lint - the number of integration points in one of the three frequency
c           interval.
c-----------------------------------------------------------------------
c
      parameter(lint=32)
c
c-----------------------------------------------------------------------
c     declarations
c-----------------------------------------------------------------------
c
c      doubleprecision nintc,neintc,nphot,nephot,lstar
c      doubleprecision fr(0:lint),func1(0:lint,0:100),
c     &                func2(0:lint,0:100),
c     &                weight(0:lint,0:100),phot(0:100)
c
      implicit real*8 (a-h,o-z)
      real*8 nintc,neintc,nphot,nephot
      real*8 lstar
      real*8 fr(0:lint),func1(0:lint,0:100),
     &                func2(0:lint,0:100),
     &                weight(0:lint,0:100),phot(0:100)
c
c-----------------------------------------------------------------------
c     common blocks
c-----------------------------------------------------------------------
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /scalin/ transmits the scaling parameters
c
      common /scalin/ scleng,scmass,sctime,scvelo,scdens,
     &                scmome,scener,sccool
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /radtab/ transmits radiative tables
c
      common /radtab/ ahphot(0:100,3),ahheat(0:100,3),
     &     ahephot(0:100,0:1,3),aheheat(0:100,0:1,3),
     &     cphot(0:100,0:5,3),nphot(0:100,0:5,3),
     &     ophot(0:100,0:5,3),nephot(0:100,0:5,3),
     ^     sphot(0:100,0:5,3)
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /radint/ transmits radiation integral cores
c
      common /radint/ frth0,frthe0,frthe1,steph0(3),
     &                h0int(0:lint,0:100,3),hh0int(0:lint,0:100,3),
     &                he0int(0:lint,0:100,3),hhe0int(0:lint,0:100,3),
     &                he1int(0:lint,0:100,3),hhe1int(0:lint,0:100,3),
     &                frsc(0:5,3),stepc(0:5,3),
     &                cintc(0:lint,0:100,0:5,3),
     &                frsn(0:5,3),stepn(0:5,3),
     &                nintc(0:lint,0:100,0:5,3),
     &                frso(0:5,3),stepo(0:5,3),
     &                ointc(0:lint,0:100,0:5,3),
     &                frsne(0:5,3),stepne(0:5,3),
     &                neintc(0:lint,0:100,0:5,3),
     &                frss(0:5,3),steps(0:5,3),
     &                sintc(0:lint,0:100,0:5,3)
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /star/ transports time dependent stellar properties
c
      common /star/ teff,rstar,lstar

      save
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c               calculate the photoionization integrals
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c-----------------------------------------------------------------------
c     if so, redo them and multiply by rstar2 later.
c     if not test whether rstar has changed and if so, determine the
c     factor needed to correct for that (stored in rstar2).
c-----------------------------------------------------------------------
c
      rstar2=rstar*rstar
      rfr=47979.72484/teff
c
c-----------------------------------------------------------------------
c     calculate the photoionisation integrals for hydrogen and helium
c-----------------------------------------------------------------------
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     frequency interval 1
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do i=0,lint
         fr(i)=frth0+steph0(1)*dfloat(i)
         do n=0,100
            weight(i,n)=steph0(1)
            func1(i,n)=h0int(i,n,1)/(dexp(fr(i)*rfr)-1.0)
            func2(i,n)=hh0int(i,n,1)/(dexp(fr(i)*rfr)-1.0)
         enddo
      enddo

      call vector_romberg (func1,weight,lint,lint,100,phot)
      do n=0,100
c     hphot(n,1)=phot(n)
         ahphot(n,1)=phot(n)
      enddo

      call vector_romberg (func2,weight,lint,lint,100,phot)
      do n=0,100
c         hheat(n,1)=phot(n)
         ahheat(n,1)=phot(n)
      enddo
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     frequency interval 2
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do i=0,lint
         fr(i)=frthe0+steph0(2)*dfloat(i)
         do n=0,100
            weight(i,n)=steph0(2)
            func1(i,n)=h0int(i,n,2)/(dexp(fr(i)*rfr)-1.0)
            func2(i,n)=hh0int(i,n,2)/(dexp(fr(i)*rfr)-1.0)
         enddo
      enddo


      call vector_romberg (func1,weight,lint,lint,100,phot)
      do n=0,100
c     hphot(n,2)=phot(n)
      ahphot(n,2)=phot(n)
      enddo
          
      call vector_romberg (func2,weight,lint,lint,100,phot)
      do n=0,100
c     hheat(n,2)=phot(n)
         ahheat(n,2)=phot(n)
      enddo

      do i=0,lint
         do n=0,100
            func1(i,n)=he0int(i,n,2)/(dexp(fr(i)*rfr)-1.0)
            func2(i,n)=hhe0int(i,n,2)/(dexp(fr(i)*rfr)-1.0)
         enddo
      enddo


      call vector_romberg (func1,weight,lint,lint,100,phot)
      do n=0,100
c         hephot(n,0,2)=phot(n)
         ahephot(n,0,2)=phot(n)
      enddo

      call vector_romberg (func2,weight,lint,lint,100,phot)
      do n=0,100
c         heheat(n,0,2)=phot(n)
         aheheat(n,0,2)=phot(n)
      enddo
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     frequency interval 3
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do i=0,lint
         fr(i)=frthe1+steph0(3)*dfloat(i)
         do n=0,100
            weight(i,n)=steph0(3)
            func1(i,n)=h0int(i,n,3)/(dexp(fr(i)*rfr)-1.0)
            func2(i,n)=hh0int(i,n,3)/(dexp(fr(i)*rfr)-1.0)
         enddo
      enddo

      call vector_romberg (func1,weight,lint,lint,100,phot)
      do n=0,100
c         hphot(n,3)=phot(n)
         ahphot(n,3)=phot(n)
      enddo

      call vector_romberg (func2,weight,lint,lint,100,phot)
      do n=0,100
c         hheat(n,3)=phot(n)
         ahheat(n,3)=phot(n)
      enddo

      do i=0,lint
         do n=0,100
            func1(i,n)=he0int(i,n,3)/(dexp(fr(i)*rfr)-1.0)
            func2(i,n)=hhe0int(i,n,3)/(dexp(fr(i)*rfr)-1.0)
         enddo
      enddo

      call vector_romberg (func1,weight,lint,lint,100,phot)
      do n=0,100
c     hephot(n,0,3)=phot(n)
         ahephot(n,0,3)=phot(n)
      enddo

      call vector_romberg (func2,weight,lint,lint,100,phot)
      do n=0,100
c     heheat(n,0,3)=phot(n)
         aheheat(n,0,3)=phot(n)
      enddo

      do i=0,lint
         do n=0,100
            func1(i,n)=he1int(i,n,3)/(dexp(fr(i)*rfr)-1.0)
            func2(i,n)=hhe1int(i,n,3)/(dexp(fr(i)*rfr)-1.0)
         enddo
      enddo

      call vector_romberg (func1,weight,lint,lint,100,phot)
      do n=0,100
c         hephot(n,1,3)=phot(n)
         ahephot(n,1,3)=phot(n)
      enddo

      call vector_romberg (func2,weight,lint,lint,100,phot)
      do n=0,100
c         heheat(n,1,3)=phot(n)
         aheheat(n,1,3)=phot(n)
      enddo
c
c-----------------------------------------------------------------------
c     calculate the photoionisation integrals for carbon
c-----------------------------------------------------------------------
c
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     frequency interval 1
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,1
         do i=0,lint
            fr(i)=frsc(k,1)+stepc(k,1)*dfloat(i)
            do n=0,100
               weight(i,n)=stepc(k,1)
               func1(i,n)=cintc(i,n,k,1)/(dexp(fr(i)*rfr)-1.0)
            enddo
         enddo
         call vector_romberg (func1,weight,lint,lint,100,phot)
         do n=0,100
            cphot(n,k,1)=phot(n)
         enddo
      enddo
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     frequency interval 2
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,2
         do i=0,lint
            fr(i)=frsc(k,2)+stepc(k,2)*dfloat(i)
            do n=0,100
               weight(i,n)=stepc(k,2)
               func1(i,n)=cintc(i,n,k,2)/(dexp(fr(i)*rfr)-1.0)
            enddo
         enddo
         call vector_romberg (func1,weight,lint,lint,100,phot)
         do n=0,100
            cphot(n,k,2)=phot(n)
         enddo
      enddo
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     frequency interval 3
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,4
         do i=0,lint
            fr(i)=frsc(k,3)+stepc(k,3)*dfloat(i)
            do n=0,100
               weight(i,n)=stepc(k,3)
               func1(i,n)=cintc(i,n,k,3)/(dexp(fr(i)*rfr)-1.0)
            enddo
         enddo
         call vector_romberg (func1,weight,lint,lint,100,phot)
         do n=0,100
            cphot(n,k,3)=phot(n)
         enddo
      enddo
c
c-----------------------------------------------------------------------
c     calculate the photoionisation integrals for nitrogen
c-----------------------------------------------------------------------
c
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     frequency interval 1
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do i=0,lint
         fr(i)=frsn(0,1)+stepn(0,1)*dfloat(i)
         do n=0,100
            weight(i,n)=stepn(0,1)
            func1(i,n)=nintc(i,n,0,1)/(dexp(fr(i)*rfr)-1.0)
         enddo
      enddo
      call vector_romberg (func1,weight,lint,lint,100,phot)
      do n=0,100
         nphot(n,0,1)=phot(n)
      enddo
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     frequency interval 2
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,2
         do i=0,lint
            fr(i)=frsn(k,2)+stepn(k,2)*dfloat(i)
            do n=0,100
               weight(i,n)=stepn(k,2)
               func1(i,n)=nintc(i,n,k,2)/(dexp(fr(i)*rfr)-1.0)
            enddo
         enddo
         call vector_romberg (func1,weight,lint,lint,100,phot)
         do n=0,100
            nphot(n,k,2)=phot(n)
         enddo
      enddo
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     frequency interval 3
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,4
         do i=0,lint
            fr(i)=frsn(k,3)+stepn(k,3)*dfloat(i)
            do n=0,100
               weight(i,n)=stepn(k,3)
               func1(i,n)=nintc(i,n,k,3)/(dexp(fr(i)*rfr)-1.0)
            enddo
         enddo
         call vector_romberg (func1,weight,lint,lint,100,phot)
         do n=0,100
            nphot(n,k,3)=phot(n)
         enddo
      enddo
c
c-----------------------------------------------------------------------
c     calculate the photoionisation integrals for oxygen
c-----------------------------------------------------------------------
c
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     frequency interval 1
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do i=0,lint
         fr(i)=frso(0,1)+stepo(0,1)*dfloat(i)
         do n=0,100
            weight(i,n)=stepo(0,1)
            func1(i,n)=ointc(i,n,0,1)/(dexp(fr(i)*rfr)-1.0)
         enddo
      enddo
      call vector_romberg (func1,weight,lint,lint,100,phot)
      do n=0,100
         ophot(n,0,1)=phot(n)
      enddo
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     frequency interval 2
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,1
         do i=0,lint
            fr(i)=frso(k,2)+stepo(k,2)*dfloat(i)
            do n=0,100
               weight(i,n)=stepo(k,2)
               func1(i,n)=ointc(i,n,k,2)/(dexp(fr(i)*rfr)-1.0)
            enddo
         enddo
         call vector_romberg (func1,weight,lint,lint,100,phot)
         do n=0,100
            ophot(n,k,2)=phot(n)
         enddo
      enddo
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     frequency interval 3
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,4
         do i=0,lint
            fr(i)=frso(k,3)+stepo(k,3)*dfloat(i)
            do n=0,100
               weight(i,n)=stepo(k,3)
               func1(i,n)=ointc(i,n,k,3)/(dexp(fr(i)*rfr)-1.0)
            enddo
         enddo
         call vector_romberg (func1,weight,lint,lint,100,phot)
         do n=0,100
            ophot(n,k,3)=phot(n)
         enddo
      enddo
c
c-----------------------------------------------------------------------
c     calculate the photoionisation integrals for neon
c-----------------------------------------------------------------------
c
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     frequency interval 1
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do i=0,lint
         fr(i)=frsne(0,1)+stepne(0,1)*dfloat(i)
         do n=0,100
            weight(i,n)=stepne(0,1)
            func1(i,n)=neintc(i,n,0,1)/(dexp(fr(i)*rfr)-1.0)
         enddo
      enddo
      call vector_romberg (func1,weight,lint,lint,100,phot)
      do n=0,100
         nephot(n,0,1)=phot(n)
      enddo
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     frequency interval 2
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,1
         do i=0,lint
            fr(i)=frsne(k,2)+stepne(k,2)*dfloat(i)
            do n=0,100
               weight(i,n)=stepne(k,2)
               func1(i,n)=neintc(i,n,k,2)/(dexp(fr(i)*rfr)-1.0)
            enddo
         enddo
         call vector_romberg (func1,weight,lint,lint,100,phot)
         do n=0,100
            nephot(n,k,2)=phot(n)
         enddo
      enddo
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     frequency interval 3
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,4
         do i=0,lint
            fr(i)=frsne(k,3)+stepne(k,3)*dfloat(i)
            do n=0,100
               weight(i,n)=stepne(k,3)
               func1(i,n)=neintc(i,n,k,3)/(dexp(fr(i)*rfr)-1.0)
            enddo
         enddo
         call vector_romberg (func1,weight,lint,lint,100,phot)
         do n=0,100
            nephot(n,k,3)=phot(n)
         enddo
      enddo
c
c-----------------------------------------------------------------------
c     calculate the photoionisation integrals for sulfur
c-----------------------------------------------------------------------
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     frequency interval 1
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,1
         do i=0,lint
            fr(i)=frss(k,1)+steps(k,1)*dfloat(i)
            do n=0,100
               weight(i,n)=steps(k,1)
               func1(i,n)=sintc(i,n,k,1)/(dexp(fr(i)*rfr)-1.0)
            enddo
         enddo
         call vector_romberg (func1,weight,lint,lint,100,phot)
         do n=0,100
            sphot(n,k,1)=phot(n)
         enddo
      enddo
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     frequency interval 2
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,3
         do i=0,lint
            fr(i)=frss(k,2)+steps(k,2)*dfloat(i)
            do n=0,100
               weight(i,n)=steps(k,2)
               func1(i,n)=sintc(i,n,k,2)/(dexp(fr(i)*rfr)-1.0)
            enddo
         enddo
         call vector_romberg (func1,weight,lint,lint,100,phot)
         do n=0,100
            sphot(n,k,2)=phot(n)
         enddo
      enddo
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     frequency interval 3
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,4
         do i=0,lint
            fr(i)=frss(k,3)+steps(k,3)*dfloat(i)
            do n=0,100
               weight(i,n)=steps(k,3)
               func1(i,n)=sintc(i,n,k,3)/(dexp(fr(i)*rfr)-1.0)
            enddo
         enddo
         call vector_romberg (func1,weight,lint,lint,100,phot)
         do n=0,100
            sphot(n,k,3)=phot(n)
         enddo
      enddo
c
c-----------------------------------------------------------------------
c     multiply the integrals with stellar radius squared
c-----------------------------------------------------------------------
c
      do n=0,100
         ahphot(n,1)=ahphot(n,1)*rstar2
         ahphot(n,2)=ahphot(n,2)*rstar2
         ahphot(n,3)=ahphot(n,3)*rstar2
         ahheat(n,1)=ahheat(n,1)*rstar2
         ahheat(n,2)=ahheat(n,2)*rstar2
         ahheat(n,3)=ahheat(n,3)*rstar2
C
         ahephot(n,0,2)=ahephot(n,0,2)*rstar2
         ahephot(n,0,3)=ahephot(n,0,3)*rstar2
         aheheat(n,0,2)=aheheat(n,0,2)*rstar2
         aheheat(n,0,3)=aheheat(n,0,3)*rstar2
         ahephot(n,1,3)=ahephot(n,1,3)*rstar2
         aheheat(n,1,3)=aheheat(n,1,3)*rstar2
         cphot(n,0,1)=cphot(n,0,1)*rstar2
         cphot(n,0,2)=cphot(n,0,2)*rstar2
         cphot(n,0,3)=cphot(n,0,3)*rstar2
         cphot(n,1,2)=cphot(n,1,2)*rstar2
         cphot(n,1,3)=cphot(n,1,3)*rstar2
         cphot(n,2,2)=cphot(n,2,2)*rstar2
         cphot(n,2,3)=cphot(n,2,3)*rstar2
         cphot(n,3,3)=cphot(n,3,3)*rstar2
         cphot(n,4,3)=cphot(n,4,3)*rstar2
         nphot(n,0,1)=nphot(n,0,1)*rstar2
         nphot(n,0,2)=nphot(n,0,2)*rstar2
         nphot(n,0,3)=nphot(n,0,3)*rstar2
         nphot(n,1,2)=nphot(n,1,2)*rstar2
         nphot(n,1,3)=nphot(n,1,3)*rstar2
         nphot(n,2,2)=nphot(n,2,2)*rstar2
         nphot(n,2,3)=nphot(n,2,3)*rstar2
         nphot(n,3,3)=nphot(n,3,3)*rstar2
         nphot(n,4,3)=nphot(n,4,3)*rstar2
         ophot(n,0,1)=ophot(n,0,1)*rstar2
         ophot(n,0,2)=ophot(n,0,2)*rstar2
         ophot(n,0,3)=ophot(n,0,3)*rstar2
         ophot(n,1,2)=ophot(n,1,2)*rstar2
         ophot(n,1,3)=ophot(n,1,3)*rstar2
         ophot(n,2,3)=ophot(n,2,3)*rstar2
         ophot(n,3,3)=ophot(n,3,3)*rstar2
         ophot(n,4,3)=ophot(n,4,3)*rstar2
         ophot(n,5,3)=ophot(n,5,3)*rstar2
         nephot(n,0,1)=nephot(n,0,1)*rstar2
         nephot(n,0,2)=nephot(n,0,2)*rstar2
         nephot(n,0,3)=nephot(n,0,3)*rstar2
         nephot(n,1,2)=nephot(n,1,2)*rstar2
         nephot(n,1,3)=nephot(n,1,3)*rstar2
         nephot(n,2,3)=nephot(n,2,3)*rstar2
         nephot(n,3,3)=nephot(n,3,3)*rstar2
         nephot(n,4,3)=nephot(n,4,3)*rstar2
         nephot(n,5,3)=nephot(n,5,3)*rstar2
         sphot(n,0,1)=sphot(n,0,1)*rstar2
         sphot(n,0,2)=sphot(n,0,2)*rstar2
         sphot(n,0,3)=sphot(n,0,3)*rstar2
         sphot(n,1,2)=sphot(n,1,2)*rstar2
         sphot(n,1,3)=sphot(n,1,3)*rstar2
         sphot(n,2,3)=sphot(n,2,3)*rstar2
         sphot(n,3,3)=sphot(n,3,3)*rstar2
         sphot(n,4,3)=sphot(n,4,3)*rstar2
         sphot(n,5,3)=sphot(n,5,3)*rstar2
      enddo
      
      return
      end
c=======================================================================
c
c
                           subroutine inirad2
c
c
c=======================================================================
c            
c            initializes constants and arrays for radiation effects
c
c                          for use with phototables
c
c                           garrelt mellema
c
c                             3-Mar-1997
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c                          introduction
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c-----------------------------------------------------------------------
c     goal
c-----------------------------------------------------------------------
c
c     this subroutine sets constants and constant arrays having to do 
c     with radiation effects. by calculating constant quantities here,
c     we save execution time, not having to calculate them every time
c     step.
c     the quantities set are:
c     1) photoionisation integral cores
c
c-----------------------------------------------------------------------
c     photoionisation integral cores
c-----------------------------------------------------------------------
c
c     the photoionisation integrals are calculated as a function of
c     optical depth at every time step. this is done to handle an
c     evolving star. a large portion of the function under the integral,
c     however, is just a function of frequency and optical 
c     depth. this part is calculated here once. they are 3d and 4d arrays.
c     the first index (i) controls frequency according to
c     fr(i) = frt + step*i
c     with frt the threshold freqency. the maximum frequency is 
c     determined in a cunning way.
c     the second index (n) controls the optical depth in 100 logarithmic
c     steps from 10**-4 to 10**4.
c     these arrays are calculated for h0 (photoionisation and heating), 
c     n0-4, and o0-5. for o0 the ionisation cross section is rather 
c     complex.
c
c-----------------------------------------------------------------------
c     references
c-----------------------------------------------------------------------
c
c     osterbrock d.e. 1989, astrophysics of gaseous nebulae and active 
c                 galactic nuclei, univ. sc. books, mill valley, ca, 
c                 usa.
c     lang, k.r. 1980, astrophysical formulae, springer-verlag.
c 
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c          formats, parameters, declarations, common blocks, data
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c-----------------------------------------------------------------------
c     implicits
c-----------------------------------------------------------------------
c
c      implicit doubleprecision (a-h,o-z), integer*4 (i-n)
c
c-----------------------------------------------------------------------
c     parameters
c
c     pi -   just boring old pi.
c     lint - the number of integration points in one of the three frequency
c           intervals.
c
c-----------------------------------------------------------------------
c
       implicit real*8 (a-h,o-z)
       parameter(pi=3.141592654,lint=32,econv=3.5742435d18)
c
c-----------------------------------------------------------------------
c     declarations
c-----------------------------------------------------------------------
c
c      doubleprecision tau(0:100),fr(0:lint),ao0fr(0:lint),
c     $     ane1fr(0:lint),ane2fr(0:lint),
c     &     h0ffr(0:lint),he0ffr(0:lint),he1ffr(0:lint),
c     $     ethe(0:1),betahe(0:1),
c     &     etc(0:5),frtc(0:5),ac(0:5),betac(0:5),sc(0:5),
c     &     etn(0:5),frtn(0:5),an(0:5),betan(0:5),sn(0:5),
c     &     eto(0:5),frto(0:5),ao(0:5),betao(0:5),so(0:5),
c     &     etne(0:5),frtne(0:5),ane(0:5),betane(0:5),sne(0:5),
c     &     ets(0:5),frts(0:5),as(0:5),betas(0:5),ss(0:5),
c     &     bb(0:lint),weight(0:lint)
c      doubleprecision nintc,neintc,m_p,lstar
      real*8 tau(0:100),fr(0:lint),ao0fr(0:lint),
     $     ane1fr(0:lint),ane2fr(0:lint),
     &     h0ffr(0:lint),he0ffr(0:lint),he1ffr(0:lint),
     $     ethe(0:1),betahe(0:1),
     &     etc(0:5),frtc(0:5),ac(0:5),betac(0:5),sc(0:5),
     &     etn(0:5),frtn(0:5),an(0:5),betan(0:5),sn(0:5),
     &     eto(0:5),frto(0:5),ao(0:5),betao(0:5),so(0:5),
     &     etne(0:5),frtne(0:5),ane(0:5),betane(0:5),sne(0:5),
     &     ets(0:5),frts(0:5),as(0:5),betas(0:5),ss(0:5),
     &     bb(0:lint),weight(0:lint)
      real*8 nintc,neintc,m_p
      real*8 lstar
      integer nk(0:3)
c
c-----------------------------------------------------------------------
c     common blocks
c-----------------------------------------------------------------------
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /scalin/ transmits the scaling parameters
c
      common /scalin/ scleng,scmass,sctime,scvelo,scdens,
     &                scmome,scener,sccool
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /radint/ transmits radiation integral cores
c
      common /radint/ frth0,frthe0,frthe1,steph0(3),
     &                h0int(0:lint,0:100,3),hh0int(0:lint,0:100,3),
     &                he0int(0:lint,0:100,3),hhe0int(0:lint,0:100,3),
     &                he1int(0:lint,0:100,3),hhe1int(0:lint,0:100,3),
     &                frsc(0:5,3),stepc(0:5,3),
     &                cintc(0:lint,0:100,0:5,3),
     &                frsn(0:5,3),stepn(0:5,3),
     &                nintc(0:lint,0:100,0:5,3),
     &                frso(0:5,3),stepo(0:5,3),
     &                ointc(0:lint,0:100,0:5,3),
     &                frsne(0:5,3),stepne(0:5,3),
     &                neintc(0:lint,0:100,0:5,3),
     &                frss(0:5,3),steps(0:5,3),
     &                sintc(0:lint,0:100,0:5,3)
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     /star/ transports time dependent stellar properties
c
      common /star/ teff,rstar,lstar
c
      save
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c                   set up the radiative tables 
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c-----------------------------------------------------------------------
c     initialize romberg routine
c     used in integr to calculate the photo-ionisation integrals
c-----------------------------------------------------------------------
c
      call romberg_initialisation(lint)
c
c-----------------------------------------------------------------------
c     set frequency scaling
c     this scaling parameter is independent of the main program scaling
c     (see main), and is only used in this subroutine.
c-----------------------------------------------------------------------
c
      sclfre=1.0d15
c
c-----------------------------------------------------------------------
c     set the values of physical quantities
c     m_p    - proton mass (in kg)
c     c      - velocity of light (in m/s)
c     hconst - (planck's constant)/(proton mass)
c              used to determine heating integral.
c     ev2fr  - conversion factor between evs and sclfre hertz
c     ev2k   - conversion factor between evs and kelvins
c     tpic2  - 2*pi/c^2 times scaling factors needed for the integral cores
c-----------------------------------------------------------------------
c
      m_p=1.672661e-24
      c=2.997925e+10
      hconst =6.6260755d-27*sclfre
c**      hconst =3.961d-03
      ev2fr=0.241838d15/sclfre
      ev2k=1.0/8.617d-05
      tpic2=2.0*pi/(c*c)*sclfre**3*sctime
c     
c-----------------------------------------------------------------------
c     find upper limits for integrals_
c     frtop1: this is the upper limit due to arithmetic precision:
c     dexp(700) exceeds double precision limit
c     frtop2: this is the upper limit due to the form of the planck
c     curve: take 10 times the frequency of maximum intensity
c-----------------------------------------------------------------------
c     
      thigh=200000.0
      tlow=2000.0
      
      frtop1=7.0d2*tlow/47979.72484
      frtop2=5.88e-05*thigh
c     
c-----------------------------------------------------------------------
c     fill the optica depth array used to fill the tables 
c     if is filled in 100 logarithmic steps from 10^-4 to 10^4.
c-----------------------------------------------------------------------
c     
      do n=1,100
         tau(n)=10.0**(-4.0+0.08*dfloat(n-1))
      enddo
      tau(0)=0.0
c
c-----------------------------------------------------------------------
c     tabulate the photo ionization integral cores for hydrogen and helium
c     one for photoionisation, one for heating.
c     atomic data from osterbrock (1989), tab.2.7 
c     and tenorio-tagle at al. (1987)
c     ionisation potentials from lang (1980), tab.34
c-----------------------------------------------------------------------
c
      eth0=13.598
      frth0=ev2fr*eth0
      ah0=6.30d-18
c**      betah0=1.34
      betah0=1.0
      sh0=2.8
      sigma1 =ah0/(m_p)*scmass/(scleng*scleng)

      ethe(0)=24.587
      frthe0=ev2fr*ethe(0)
      ahe0=7.83d-18
      betahe(0)=1.0
      she0=1.7
      sigma2=ahe0/(m_p)*scmass/(scleng*scleng)

      ethe(1)=54.416
      frthe1=ev2fr*ethe(1)
      ahe1=1.58d-18
      betahe(1)=1.0
      she1=2.8
      sigma3 =ahe1/(m_p)*scmass/(scleng*scleng)
c
c-----------------------------------------------------------------------
c     ionizing flux
c-----------------------------------------------------------------------
c
      rfr=47979.72484/teff
      frmax=dmin1(frtop1,10.0*frtop2)
      stepfl=(frmax-frth0)/dfloat(lint)
      do i=0,lint
         fr(i)=frth0+stepfl*dfloat(i)
         weight(i)=stepfl
         bb(i)=tpic2*fr(i)*fr(i)/(dexp(fr(i)*rfr)-1.0)
      enddo
      flux=romberg(bb,weight,lint,lint,0)
      flux=4.0*pi*rstar*rstar*flux
c**      call romberg_initialisation(lint)
c
c-----------------------------------------------------------------------
c     calculate the correction factors for tau2 and tau3 (see ionic)
c-----------------------------------------------------------------------
c
c**      tau2h=(0.63*frth0/frthe0)**she0
c**      tau3h=(frth0/frthe1)**she1
c**      tau3he0=(1.51*frthe0/frthe1)**she1
c
c-----------------------------------------------------------------------
c     h0, he0-1
c-----------------------------------------------------------------------
c
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 1
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      if (frth0.lt.frtop1) then

         frmax=dmin1(frthe0,frtop1)
         steph0(1)=(frmax-frth0)/dfloat(lint)
         do i=0,lint
            fr(i)=frth0+steph0(1)*dfloat(i)
            h0ffr(i)=(betah0*(fr(i)/frth0)**(-sh0)+
     &           (1.0-betah0)*(fr(i)/frth0)**(-sh0-1.0))
            do n=0,100
               h0int(i,n,1)=tpic2*ah0*h0ffr(i)*
     &              fr(i)*fr(i)*dexp(-tau(n)*h0ffr(i))
               hh0int(i,n,1)=hconst*(fr(i)-frth0)*h0int(i,n,1)
            enddo
         enddo
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 2
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
         if (frthe0.lt.frtop1) then

            frmax=dmin1(frthe1,frtop1)
            steph0(2)=(frmax-frthe0)/dfloat(lint)
            do i=0,lint
               fr(i)=frthe0+steph0(2)*dfloat(i)
               h0ffr(i)=(betah0*(fr(i)/frth0)**(-sh0)+
     &              (1.0-betah0)*(fr(i)/frth0)**(-sh0-1.0))
               he0ffr(i)=(betahe(0)*(fr(i)/frthe0)**(-she0)+
     &              (1.0-betahe(0))*(fr(i)/frthe0)**(-she0-1.0))
               do n=0,100
                  h0int(i,n,2)=tpic2*ah0*h0ffr(i)*
     &                 fr(i)*fr(i)*dexp(-tau(n)*he0ffr(i))
                  hh0int(i,n,2)=hconst*(fr(i)-frth0)*h0int(i,n,2)
                  he0int(i,n,2)=tpic2*ahe0*he0ffr(i)*
     &                 fr(i)*fr(i)*dexp(-tau(n)*he0ffr(i))
                  hhe0int(i,n,2)=hconst*(fr(i)-frthe0)*he0int(i,n,2)
               enddo
            enddo
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 3
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
            if (frthe1.lt.frtop1) then

               frmax=dmin1(frtop1,10.0*dmax1(frthe1,frtop2))
               steph0(3)=(frmax-frthe1)/dfloat(lint)
               do i=0,lint
                  fr(i)=frthe1+steph0(3)*dfloat(i)
                  h0ffr(i)=(betah0*(fr(i)/frth0)**(-sh0)+
     &                 (1.0-betah0)*(fr(i)/frth0)**(-sh0-1.0))
                  he0ffr(i)=(betahe(0)*(fr(i)/frthe0)**(-she0)+
     &                 (1.0-betahe(0))*(fr(i)/frthe0)**(-she0-1.0))
                  he1ffr(i)=(betahe(1)*(fr(i)/frthe1)**(-she1)+
     &                 (1.0-betahe(1))*(fr(i)/frthe1)**(-she1-1.0))
                  do n=0,100
                     h0int(i,n,3)=tpic2*ah0*h0ffr(i)*
     &                    fr(i)*fr(i)*dexp(-tau(n)*he1ffr(i))
                     hh0int(i,n,3)=hconst*(fr(i)-frth0)*h0int(i,n,3)
                     he0int(i,n,3)=tpic2*ahe0*he0ffr(i)*
     &                    fr(i)*fr(i)*dexp(-tau(n)*he1ffr(i))
                     hhe0int(i,n,3)=hconst*(fr(i)-frthe0)*he0int(i,n,3)
                     he1int(i,n,3)=tpic2*ahe1*he1ffr(i)*
     &                    fr(i)*fr(i)*dexp(-tau(n)*he1ffr(i))
                     hhe1int(i,n,3)=hconst*(fr(i)-frthe1)*he1int(i,n,3)
                  enddo
               enddo
            endif
         endif
      endif
c
c-----------------------------------------------------------------------
c     tabulate the photo ionization integral cores for carbon.
c     atomic data from osterbrock (1989), tab.2.7 and Raymond (1976)
c     ionisation potentials from Cox (1970)
c-----------------------------------------------------------------------
c
      etc(0)=11.3
      frtc(0)=ev2fr*etc(0)
      ac(0)=12.2d-18
      betac(0)=3.32
      sc(0)=2.0

      etc(1)=24.4
      frtc(1)=ev2fr*etc(1)
      ac(1)=4.60d-18
      betac(1)=1.95
      sc(1)=3.0

      etc(2)=47.9
      frtc(2)=ev2fr*etc(2)
      ac(2)=1.60d-18
      betac(2)=2.60
      sc(2)=3.0

      etc(3)=64.5
      frtc(3)=ev2fr*etc(3)
      ac(3)=0.68d-18
      betac(3)=1.0
      sc(3)=2.0

      etc(4)=392.0
      frtc(4)=ev2fr*etc(4)
      ac(4)=0.472d-18
      betac(4)=1.42
      sc(4)=2.0
c
c-----------------------------------------------------------------------
c     c0-4
c-----------------------------------------------------------------------
c
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 1
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c

      do k=0,1
         frsc(k,1)=dmax1(frtc(k),frth0)
         if (frsc(k,1).lt.frtop1) then
            frmax=dmin1(frthe0,frtop1)
            stepc(k,1)=(frmax-frsc(k,1))/dfloat(lint)
            do i=0,lint
               fr(i)=frsc(k,1)+stepc(k,1)*dfloat(i)
               h0ffr(i)=(betah0*(fr(i)/frth0)**(-sh0)+
     &              (1.0-betah0)*(fr(i)/frth0)**(-sh0-1.0))
               do n=0,100
                  cintc(i,n,k,1)=tpic2*ac(k)*
     &                 (betac(k)*(fr(i)/frtc(k))**(-sc(k))+
     &                 (1.0-betac(k))*(fr(i)/frtc(k))**(-sc(k)-1.0))*
     &              fr(i)*fr(i)*dexp(-tau(n)*h0ffr(i))
               enddo
            enddo
         endif
      enddo
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 2
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,2
         frsc(k,2)=dmax1(frtc(k),frthe0)
         if (frsc(k,2).lt.frtop1) then
            frmax=dmin1(frthe1,frtop1)
            stepc(k,2)=(frmax-frsc(k,2))/dfloat(lint)
            do i=0,lint
               fr(i)=frsc(k,2)+stepc(k,2)*dfloat(i)
               he0ffr(i)=(betahe(0)*(fr(i)/frthe0)**(-she0)+
     &              (1.0-betahe(0))*(fr(i)/frthe0)**(-she0-1.0))
               do n=0,100
                  cintc(i,n,k,2)=tpic2*ac(k)*
     &                 (betac(k)*(fr(i)/frtc(k))**
     &                 (-sc(k))+
     &                 (1.0-betac(k))*(fr(i)/frtc(k))**
     &                 (-sc(k)-1.0))*
     &                 fr(i)*fr(i)*dexp(-tau(n)*he0ffr(i))
               enddo
            enddo
         endif
      enddo
c     
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 3
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,4
         frsc(k,3)=dmax1(frtc(k),frthe1)
         if (frsc(k,3).lt.frtop1) then
            frmax=dmin1(frtop1,10.0*dmax1(frsc(k,3),frtop2))
            stepc(k,3)=(frmax-frsc(k,3))/dfloat(lint)
            do i=0,lint
               fr(i)=frsc(k,3)+stepc(k,3)*dfloat(i)
               he1ffr(i)=(betahe(1)*(fr(i)/frthe1)**(-she1)+
     &              (1.0-betahe(1))*(fr(i)/frthe1)**(-she1-1.0))
               do n=0,100
                  cintc(i,n,k,3)=tpic2*ac(k)*
     &                 (betac(k)*(fr(i)/frtc(k))**
     &                 (-sc(k))+
     &                 (1.0-betac(k))*(fr(i)/frtc(k))**
     &                 (-sc(k)-1.0))*
     &                 fr(i)*fr(i)*dexp(-tau(n)*he1ffr(i))
               enddo
            enddo
         endif
      enddo
c
c-----------------------------------------------------------------------
c     tabulate the photo ionization integral cores for nitrogen.
c     atomic data from osterbrock (1989), tab.2.7
c     ionisation potentials from lang (1980), tab.34
c-----------------------------------------------------------------------
c
      etn(0)=14.354
      frtn(0)=ev2fr*etn(0)
      an(0)=11.4d-18
      betan(0)=4.29
      sn(0)=2.0

      etn(1)=29.601
      frtn(1)=ev2fr*etn(1)
      an(1)=6.65d-18
      betan(1)=2.86
      sn(1)=3.0

      etn(2)=47.448
      frtn(2)=ev2fr*etn(2)
      an(2)=2.06d-18
      betan(2)=1.63
      sn(2)=3.0

      etn(3)=77.472
      frtn(3)=ev2fr*etn(3)
      an(3)=1.08d-18
      betan(3)=2.6
      sn(3)=3.0

      etn(4)=97.888
      frtn(4)=ev2fr*etn(4)
      an(4)=0.48d-18
      betan(4)=1.0
      sn(4)=2.0
c
c-----------------------------------------------------------------------
c     n0-4
c-----------------------------------------------------------------------
c
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 1
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c

      k=0
      frsn(k,1)=dmax1(frtn(k),frth0)
      if (frsn(k,1).lt.frtop1) then
         frmax=dmin1(frthe0,frtop1)
         stepn(k,1)=(frmax-frsn(k,1))/dfloat(lint)
         do i=0,lint
            fr(i)=frsn(k,1)+stepn(k,1)*dfloat(i)
            h0ffr(i)=(betah0*(fr(i)/frth0)**(-sh0)+
     &           (1.0-betah0)*(fr(i)/frth0)**(-sh0-1.0))
            do n=0,100
               nintc(i,n,k,1)=tpic2*an(k)*
     &              (betan(k)*(fr(i)/frtn(k))**(-sn(k))+
     &              (1.0-betan(k))*(fr(i)/frtn(k))**(-sn(k)-1.0))*
     &              fr(i)*fr(i)*dexp(-tau(n)*h0ffr(i))
            enddo
         enddo
      endif
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 2
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,2
         frsn(k,2)=dmax1(frtn(k),frthe0)
         if (frsn(k,2).lt.frtop1) then
            frmax=dmin1(frthe1,frtop1)
            stepn(k,2)=(frmax-frsn(k,2))/dfloat(lint)
            do i=0,lint
               fr(i)=frsn(k,2)+stepn(k,2)*dfloat(i)
               he0ffr(i)=(betahe(0)*(fr(i)/frthe0)**(-she0)+
     &              (1.0-betahe(0))*(fr(i)/frthe0)**(-she0-1.0))
               do n=0,100
                  nintc(i,n,k,2)=tpic2*an(k)*
     &                 (betan(k)*(fr(i)/frtn(k))**
     &                 (-sn(k))+
     &                 (1.0-betan(k))*(fr(i)/frtn(k))**
     &                 (-sn(k)-1.0))*
     &                 fr(i)*fr(i)*dexp(-tau(n)*he0ffr(i))
               enddo
            enddo
         endif
      enddo
c     
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 3
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,4
         frsn(k,3)=dmax1(frtn(k),frthe1)
         if (frsn(k,3).lt.frtop1) then
            frmax=dmin1(frtop1,10.0*dmax1(frsn(k,3),frtop2))
            stepn(k,3)=(frmax-frsn(k,3))/dfloat(lint)
            do i=0,lint
               fr(i)=frsn(k,3)+stepn(k,3)*dfloat(i)
               he1ffr(i)=(betahe(1)*(fr(i)/frthe1)**(-she1)+
     &              (1.0-betahe(1))*(fr(i)/frthe1)**(-she1-1.0))
               do n=0,100
                  nintc(i,n,k,3)=tpic2*an(k)*
     &                 (betan(k)*(fr(i)/frtn(k))**
     &                 (-sn(k))+
     &                 (1.0-betan(k))*(fr(i)/frtn(k))**
     &                 (-sn(k)-1.0))*
     &                 fr(i)*fr(i)*dexp(-tau(n)*he1ffr(i))
               enddo
            enddo
         endif
      enddo
c
c-----------------------------------------------------------------------
c     tabulate the photo ionization integral cores for oxygen.
c     atomic data from osterbrock (1989), tab.2.7
c     ionisation potentials from lang (1980), tab.34
c-----------------------------------------------------------------------
c
      eto(0)=13.618
      frto(0)=ev2fr*eto(0)

      eto(1)=35.116
      frto(1)=ev2fr*eto(1)
      ao(1)=7.32d-18
      betao(1)=3.84
      so(1)=2.5

      eto(2)=54.934
      frto(2)=ev2fr*eto(2)
      ao(2)=3.65d-18
      betao(2)=2.01
      so(2)=3.0

      eto(3)=77.412
      frto(3)=ev2fr*eto(3)
      ao(3)=1.27d-18
      betao(3)=0.83
      so(3)=3.0

      eto(4)=113.896
      frto(4)=ev2fr*eto(4)
      ao(4)=0.78d-18
      betao(4)=2.6
      so(4)=3.0

      eto(5)=138.116
      frto(5)=ev2fr*eto(5)
      ao(5)=0.36d-18
      betao(5)=1.0
      so(5)=2.1
c
c-----------------------------------------------------------------------
c     o0
c
c     neutral oxygen has a rather complicated ionization cross section.
c     therefore it is treated seperately
c-----------------------------------------------------------------------
c
      frto0a=frto(0)
      frto0b=ev2fr*16.9
      frto0c=ev2fr*18.6

      ao0a=2.94d-18
      ao0b=3.85d-18
      ao0c=2.26d-18

      beto0a=2.66
      beto0b=4.38
      beto0c=4.31

      so0a=1.0
      so0b=1.5
      so0c=1.5
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 1
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     fill up the frequency array

      frso(0,1)=frto(0)
      if (frso(0,1).lt.frtop1) then
         frmax=dmin1(frthe0,frtop1)
         stepo(0,1)=(frmax-frso(0,1))/dfloat(lint)
         do i=0,lint
            fr(i)=frso(0,1)+dfloat(i)*stepo(0,1)
            h0ffr(i)=(betah0*(fr(i)/frth0)**(-sh0)+
     &           (1.0-betah0)*(fr(i)/frth0)**(-sh0-1.0))
         enddo
          
c     find the thresholds in frequqncy index
          
         it2=int((frto0b-frto0a)/stepo(0,1))
         it3=int((frto0c-frto0a)/stepo(0,1))
         
         
c     calculate the a
         
         do i=0,lint
            ao0fr(i)=ao0a*(beto0a*(fr(i)/frto0a)**(-so0a)+
     &           (1.0-beto0a)*(fr(i)/frto0a)**(-so0a-1.0))
         enddo
          
         if (frto0b.lt.frtop1) then
            do i=it2+1,lint
               ao0fr(i)=ao0fr(i)+
     &              ao0b*(beto0b*(fr(i)/frto0b)**(-so0b)+
     &              (1.0-beto0b)*(fr(i)/frto0b)**(-so0b-1.0))
            enddo
         endif
          
         if (frto0c.lt.frtop1) then
            do i=it3+1,lint
               ao0fr(i)=ao0fr(i)+
     &              ao0c*(beto0c*(fr(i)/frto0c)**(-so0c)+
     &              (1.0-beto0c)*(fr(i)/frto0c)**(-so0c-1.0))
            enddo
         endif
          
         do i=0,lint
            do n=0,100
               ointc(i,n,0,1)=tpic2*ao0fr(i)*
     &              fr(i)*fr(i)*dexp(-tau(n)*h0ffr(i))
            enddo
         enddo
      endif
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 2
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c

c     fill up the frequency array

      frso(0,2)=frthe0
      if (frso(0,2).lt.frtop1) then
         frmax=dmin1(frthe1,frtop1)
         stepo(0,2)=(frmax-frso(0,2))/dfloat(lint)
         do i=0,lint
            fr(i)=frso(0,2)+dfloat(i)*stepo(0,2)
            he0ffr(i)=(betahe(0)*(fr(i)/frthe0)**(-she0)+
     &           (1.0-betahe(0))*(fr(i)/frthe0)**(-she0-1.0))
         enddo

c     calculate the a

         do i=0,lint
            ao0fr(i)=ao0a*(beto0a*(fr(i)/frto0a)**(-so0a)+
     &           (1.0-beto0a)*(fr(i)/frto0a)**(-so0a-1.0))+
     &           ao0b*(beto0b*(fr(i)/frto0b)**(-so0b)+
     &           (1.0-beto0b)*(fr(i)/frto0b)**(-so0b-1.0))+
     &           ao0c*(beto0c*(fr(i)/frto0c)**(-so0c)+
     &           (1.0-beto0c)*(fr(i)/frto0c)**(-so0c-1.0))
         enddo
 
         do i=0,lint
            do n=0,100
               ointc(i,n,0,2)=tpic2*ao0fr(i)*
     &              fr(i)*fr(i)*dexp(-tau(n)*he0ffr(i))
            enddo
         enddo
      endif
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 3
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c

c     fill up the frequency array

      frso(0,3)=frthe1
      if (frso(0,3).lt.frtop1) then
         frmax=dmin1(frtop1,10.0*dmax1(frso(0,3),frtop2))
         stepo(0,3)=(frmax-frso(0,3))/dfloat(lint)
         do i=0,lint
            fr(i)=frso(0,3)+dfloat(i)*stepo(0,3)
            he1ffr(i)=(betahe(1)*(fr(i)/frthe1)**(-she1)+
     &           (1.0-betahe(1))*(fr(i)/frthe1)**(-she1-1.0))
         enddo

c     calculate the a

         do i=0,lint
            ao0fr(i)=ao0a*(beto0a*(fr(i)/frto0a)**(-so0a)+
     &           (1.0-beto0a)*(fr(i)/frto0a)**(-so0a-1.0))+
     &           ao0b*(beto0b*(fr(i)/frto0b)**(-so0b)+
     &           (1.0-beto0b)*(fr(i)/frto0b)**(-so0b-1.0))+
     &           ao0c*(beto0c*(fr(i)/frto0c)**(-so0c)+
     &           (1.0-beto0c)*(fr(i)/frto0c)**(-so0c-1.0))
         enddo
 
         do i=0,lint
            do n=0,100
               ointc(i,n,0,3)=tpic2*ao0fr(i)*
     &              fr(i)*fr(i)*dexp(-tau(n)*he1ffr(i))
            enddo
         enddo
      endif
c
c-----------------------------------------------------------------------
c     o1-5
c-----------------------------------------------------------------------
c
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 2
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=1,1
         frso(k,2)=dmax1(frto(k),frthe0)
         if (frso(k,2).lt.frtop1) then
            frmax=dmin1(frthe1,frtop1)
            stepo(k,2)=(frmax-frso(k,2))/dfloat(lint)
            do i=0,lint
               fr(i)=frso(k,2)+stepo(k,2)*dfloat(i)
               he0ffr(i)=(betahe(0)*(fr(i)/frthe0)**(-she0)+
     &              (1.0-betahe(0))*(fr(i)/frthe0)**(-she0-1.0))
               do n=0,100
                  ointc(i,n,k,2)=tpic2*ao(k)*
     &                 (betao(k)*(fr(i)/frto(k))**(-so(k))+
     &                 (1.0-betao(k))*(fr(i)/frto(k))**
     &                 (-so(k)-1.0))*
     &                 fr(i)*fr(i)*dexp(-tau(n)*he0ffr(i))
               enddo
            enddo
         endif
      enddo
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 3
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=1,5
          frso(k,3)=dmax1(frto(k),frthe1)
          if (frso(k,3).lt.frtop1) then
             frmax=dmin1(frtop1,10.0*dmax1(frso(k,3),frtop2))
             stepo(k,3)=(frmax-frso(k,3))/dfloat(lint)
             do i=0,lint
                fr(i)=frso(k,3)+stepo(k,3)*dfloat(i)
                he1ffr(i)=(betahe(1)*(fr(i)/frthe1)**(-she1)+
     &               (1.0-betahe(1))*(fr(i)/frthe1)**(-she1-1.0))
                do n=0,100
                   ointc(i,n,k,3)=tpic2*ao(k)*
     &                  (betao(k)*(fr(i)/frto(k))**(-so(k))+
     &                  (1.0-betao(k))*(fr(i)/frto(k))**
     &                  (-so(k)-1.0))*
     &                  fr(i)*fr(i)*dexp(-tau(n)*he1ffr(i))
                enddo
             enddo
          endif
       enddo
c
c-----------------------------------------------------------------------
c     tabulate the photo ionization integral cores for neon.
c     atomic data from osterbrock (1989), tab.2.7
c     ionisation potentials from lang (1980), tab.34
c-----------------------------------------------------------------------
c
      etne(0)=21.6
      frtne(0)=ev2fr*etne(0)
      ane(0)=5.35d-18
      betane(0)=3.77
      sne(0)=1.0

      etne(1)=41.1
      frtne(1)=ev2fr*etne(1)

      etne(2)=63.5
      frtne(2)=ev2fr*eto(2)

      etne(3)=97.0
      frtne(3)=ev2fr*etne(3)
      ane(3)=3.11d-18
      betane(3)=1.96
      sne(3)=3.0

      etne(4)=126.0
      frtne(4)=ev2fr*etne(4)
      ane(4)=1.40d-18
      betane(4)=1.47
      sne(4)=3.0

      etne(5)=158.0
      frtne(5)=ev2fr*etne(5)
      ane(5)=0.49d-18
      betane(5)=1.15
      sne(5)=3.0
c
c-----------------------------------------------------------------------
c     ne0,3,4
c-----------------------------------------------------------------------
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 1
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      k=0
      frsne(k,1)=dmax1(frtne(k),frth0)
      if (frsne(k,1).lt.frtop1) then
         frmax=dmin1(frthe0,frtop1)
         stepne(k,1)=(frmax-frsne(k,1))/dfloat(lint)
         do i=0,lint
            fr(i)=frsne(k,1)+stepne(k,1)*dfloat(i)
            h0ffr(i)=(betah0*(fr(i)/frth0)**(-sh0)+
     &           (1.0-betah0)*(fr(i)/frth0)**(-sh0-1.0))
            do n=0,100
               neintc(i,n,k,1)=tpic2*ane(k)*
     &              (betane(k)*(fr(i)/frtne(k))**(-sne(k))+
     &              (1.0-betane(k))*(fr(i)/frtne(k))**(-sne(k)-1.0))*
     &              fr(i)*fr(i)*dexp(-tau(n)*h0ffr(i))
            enddo
         enddo
      endif
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 2
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      k=0
      frsne(k,2)=dmax1(frtne(k),frthe0)
      if (frsne(k,2).lt.frtop1) then
         frmax=dmin1(frthe1,frtop1)
         stepne(k,2)=(frmax-frsne(k,2))/dfloat(lint)
         do i=0,lint
            fr(i)=frsne(k,2)+stepne(k,2)*dfloat(i)
            he0ffr(i)=(betahe(0)*(fr(i)/frthe0)**(-she0)+
     &           (1.0-betahe(0))*(fr(i)/frthe0)**(-she0-1.0))
            do n=0,100
               neintc(i,n,k,2)=tpic2*ane(k)*
     &              (betane(k)*(fr(i)/frtne(k))**
     &              (-sne(k))+
     &              (1.0-betane(k))*(fr(i)/frtne(k))**
     &              (-sne(k)-1.0))*
     &              fr(i)*fr(i)*dexp(-tau(n)*he0ffr(i))
            enddo
         enddo
      endif
c     
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 3
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      nk(0)=0
      nk(1)=3
      nk(2)=4
      nk(3)=5
      do lk=0,3
         k=nk(lk)
         frsne(k,3)=dmax1(frtne(k),frthe1)
         if (frsne(k,3).lt.frtop1) then
            frmax=dmin1(frtop1,10.0*dmax1(frsne(k,3),frtop2))
            stepne(k,3)=(frmax-frsne(k,3))/dfloat(lint)
            do i=0,lint
               fr(i)=frsne(k,3)+stepne(k,3)*dfloat(i)
               he1ffr(i)=(betahe(1)*(fr(i)/frthe1)**(-she1)+
     &              (1.0-betahe(1))*(fr(i)/frthe1)**(-she1-1.0))
               do n=0,100
                  neintc(i,n,k,3)=tpic2*ane(k)*
     &                 (betane(k)*(fr(i)/frtne(k))**
     &                 (-sne(k))+
     &                 (1.0-betane(k))*(fr(i)/frtne(k))**
     &                 (-sne(k)-1.0))*
     &                 fr(i)*fr(i)*dexp(-tau(n)*he1ffr(i))
               enddo
            enddo
         endif
      enddo
c
c-----------------------------------------------------------------------
c     ne1
c
c     ne1 has a rather complicated ionization cross section.
c     therefore it is treated seperately
c-----------------------------------------------------------------------
c
      frtne1a=frtne(1)
      frtne1b=ev2fr*16.9
      frtne1c=ev2fr*18.6

      ane1a=4.16d-18
      ane1b=2.71d-18
      ane1c=0.52d-18

      betne1a=2.72
      betne1b=2.15
      betne1c=2.13

      sne1a=1.5
      sne1b=1.5
      sne1c=1.5
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 2
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c

c     fill up the frequency array

      frsne(1,2)=frthe0
      if (frsne(1,2).lt.frtop1) then
         frmax=dmin1(frthe1,frtop1)
         stepne(1,2)=(frmax-frsne(1,2))/dfloat(lint)
         do i=0,lint
            fr(i)=frsne(1,2)+dfloat(i)*stepne(1,2)
            he0ffr(i)=(betahe(0)*(fr(i)/frthe0)**(-she0)+
     &           (1.0-betahe(0))*(fr(i)/frthe0)**(-she0-1.0))
         enddo

c     calculate the a

         do i=0,lint
            ane1fr(i)=ane1a*(betne1a*(fr(i)/frtne1a)**(-sne1a)+
     &           (1.0-betne1a)*(fr(i)/frtne1a)**(-sne1a-1.0))+
     &           ane1b*(betne1b*(fr(i)/frtne1b)**(-sne1b)+
     &           (1.0-betne1b)*(fr(i)/frtne1b)**(-sne1b-1.0))+
     &           ane1c*(betne1c*(fr(i)/frtne1c)**(-sne1c)+
     &           (1.0-betne1c)*(fr(i)/frtne1c)**(-sne1c-1.0))
         enddo
         
         do i=0,lint
            do n=0,100
               neintc(i,n,1,2)=tpic2*ane1fr(i)*
     &              fr(i)*fr(i)*dexp(-tau(n)*he0ffr(i))
            enddo
         enddo
      endif
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 3
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c

c     fill up the frequency array

      frsne(1,3)=frthe1
      if (frsne(1,3).lt.frtop1) then
          frmax=dmin1(frtop1,10.0*dmax1(frsne(1,3),frtop2))
          stepne(1,3)=(frmax-frsne(1,3))/dfloat(lint)
          do i=0,lint
              fr(i)=frsne(1,3)+dfloat(i)*stepne(1,3)
              he1ffr(i)=(betahe(1)*(fr(i)/frthe1)**(-she1)+
     &            (1.0-betahe(1))*(fr(i)/frthe1)**(-she1-1.0))
           enddo

c     calculate the a

          do i=0,lint
              ane1fr(i)=ane1a*(betne1a*(fr(i)/frtne1a)**(-sne1a)+
     &            (1.0-betne1a)*(fr(i)/frtne1a)**(-sne1a-1.0))+
     &            ane1b*(betne1b*(fr(i)/frtne1b)**(-sne1b)+
     &            (1.0-betne1b)*(fr(i)/frtne1b)**(-sne1b-1.0))+
     &            ane1c*(betne1c*(fr(i)/frtne1c)**(-sne1c)+
     &            (1.0-betne1c)*(fr(i)/frtne1c)**(-sne1c-1.0))
           enddo
 
          do i=0,lint
              do n=0,100
                  neintc(i,n,1,3)=tpic2*ane1fr(i)*
     &                fr(i)*fr(i)*dexp(-tau(n)*he1ffr(i))
               enddo
            enddo
      endif
c
c-----------------------------------------------------------------------
c     ne2
c
c     ne2 has a rather complicated ionization cross section.
c     therefore it is treated seperately
c-----------------------------------------------------------------------
c
c
      frtne2a=frtne(2)
      frtne2b=ev2fr*16.9
      frtne2c=ev2fr*18.6

      ane2a=1.80d-18
      ane2b=2.50d-18
      ane2c=1.48d-18

      betne2a=2.28
      betne2b=2.35
      betne2c=2.23

      sne2a=2.0
      sne2b=2.5
      sne2c=2.5
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 3
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c

c     fill up the frequency array

      frsne(2,3)=frthe1
      if (frsne(2,3).lt.frtop1) then
         frmax=dmin1(frtop1,10.0*dmax1(frsne(2,3),frtop2))
         stepne(2,3)=(frmax-frsne(2,3))/dfloat(lint)
         do i=0,lint
            fr(i)=frsne(2,3)+dfloat(i)*stepne(2,3)
            he1ffr(i)=(betahe(1)*(fr(i)/frthe1)**(-she1)+
     &           (1.0-betahe(1))*(fr(i)/frthe1)**(-she1-1.0))
         enddo
         
c     calculate the a
         
         do  i=0,lint
            ane2fr(i)=ane2a*(betne2a*(fr(i)/frtne2a)**(-sne2a)+
     &           (1.0-betne2a)*(fr(i)/frtne2a)**(-sne2a-1.0))+
     &           ane2b*(betne2b*(fr(i)/frtne2b)**(-sne2b)+
     &           (1.0-betne2b)*(fr(i)/frtne2b)**(-sne2b-1.0))+
     &           ane2c*(betne2c*(fr(i)/frtne2c)**(-sne2c)+
     &           (1.0-betne2c)*(fr(i)/frtne2c)**(-sne2c-1.0))
         enddo
         
         do i=0,lint
            do n=0,100
               neintc(i,n,2,3)=tpic2*ane2fr(i)*
     &              fr(i)*fr(i)*dexp(-tau(n)*he1ffr(i))
            enddo
         enddo
      endif
c
c-----------------------------------------------------------------------
c     tabulate the photo ionization integral cores for sulfur.
c     atomic data from Raymond (1976)
c     ionisation potentials from Raymond (1976)
c-----------------------------------------------------------------------
c
      ets(0)=10.4
      frts(0)=ev2fr*ets(0)
      as(0)=12.62d-18
      betas(0)=21.595
      ss(0)=3.0

      ets(1)=23.4
      frts(1)=ev2fr*ets(1)
      as(1)=8.20d-18
      betas(1)=1.695
      ss(1)=1.5

      ets(2)=35.0
      frts(2)=ev2fr*ets(2)
      as(2)=0.38d-18
      betas(2)=18.427
      ss(2)=2.0

      ets(3)=47.3
      frts(3)=ev2fr*ets(3)
      as(3)=0.29d-18
      betas(3)=6.837
      ss(3)=2.0

      ets(4)=73.0
      frts(4)=ev2fr*ets(4)
      as(4)=1.29d-18
      betas(4)=1.481
      ss(4)=3.02

      ets(5)=88.0
      frts(5)=ev2fr*ets(5)
      as(5)=0.195d-18
      betas(5)=0.4
      ss(5)=3.0
c
c-----------------------------------------------------------------------
c     s0-4
c-----------------------------------------------------------------------
c
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 1
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,1
         frss(k,1)=dmax1(frts(k),frth0)
         if (frss(k,1).lt.frtop1) then
            frmax=dmin1(frthe0,frtop1)
            steps(k,1)=(frmax-frss(k,1))/dfloat(lint)
            do i=0,lint
               fr(i)=frss(k,1)+steps(k,1)*dfloat(i)
               h0ffr(i)=(betah0*(fr(i)/frth0)**(-sh0)+
     &              (1.0-betah0)*(fr(i)/frth0)**(-sh0-1.0))
               do n=0,100
                  sintc(i,n,k,1)=tpic2*as(k)*
     &                 (betas(k)*(fr(i)/frts(k))**(-ss(k))+
     &                 (1.0-betas(k))*(fr(i)/frts(k))**(-ss(k)-1.0))*
     &                 fr(i)*fr(i)*dexp(-tau(n)*h0ffr(i))
               enddo
            enddo
         endif
      enddo
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 2
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,3
         frss(k,2)=dmax1(frts(k),frthe0)
         if (frss(k,2).lt.frtop1) then
            frmax=dmin1(frthe1,frtop1)
            steps(k,2)=(frmax-frss(k,2))/dfloat(lint)
            do i=0,lint
               fr(i)=frss(k,2)+steps(k,2)*dfloat(i)
               he0ffr(i)=(betahe(0)*(fr(i)/frthe0)**(-she0)+
     &              (1.0-betahe(0))*(fr(i)/frthe0)**(-she0-1.0))
               do n=0,100
                  sintc(i,n,k,2)=tpic2*as(k)*
     &                 (betas(k)*(fr(i)/frts(k))**
     &                 (-ss(k))+
     &                 (1.0-betas(k))*(fr(i)/frts(k))**
     &                 (-ss(k)-1.0))*
     &                 fr(i)*fr(i)*dexp(-tau(n)*he0ffr(i))
               enddo
            enddo
         endif
      enddo
c     
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     fequency range 3
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do k=0,5
         frss(k,3)=dmax1(frts(k),frthe1)
         if (frss(k,3).lt.frtop1) then
            frmax=dmin1(frtop1,10.0*dmax1(frss(k,3),frtop2))
            steps(k,3)=(frmax-frss(k,3))/dfloat(lint)
            do i=0,lint
               fr(i)=frss(k,3)+steps(k,3)*dfloat(i)
               he1ffr(i)=(betahe(1)*(fr(i)/frthe1)**(-she1)+
     &              (1.0-betahe(1))*(fr(i)/frthe1)**(-she1-1.0))
               do n=0,100
                  sintc(i,n,k,3)=tpic2*as(k)*
     &                 (betas(k)*(fr(i)/frts(k))**
     &                 (-ss(k))+
     &                 (1.0-betas(k))*(fr(i)/frts(k))**
     &                 (-ss(k)-1.0))*
     &                 fr(i)*fr(i)*dexp(-tau(n)*he1ffr(i))
               enddo
            enddo
         endif
      enddo

      return
      end

C=============================================================================
C
C
      SUBROUTINE Vector_Romberg (f, w, nc, nx, ny, itg)            
C
C
C=============================================================================
C                                                      
C      This subroutin performs an one-dimensional integration in      
C      the x-direction for every y and stores the results in itg.      
C      In principle it works the same as the function Romberg.            
C_______________________________________________________________________
C      Extra argument : itg      the results of the integrations            
C_______________________________________________________________________

C**      implicit none
c        implicit doubleprecision (a-h,o-z), integer (i-n)
        implicit real*8 (a-h,o-z)
      parameter (maxpow = 14)
      integer nc, nx, ny, px, x, y
c      DoublePrecision f (0:nc, 0:ny), w (0:nc, 0:ny),
c     &          itg (0:ny),romw (0:2**maxpow, -1:maxpow)
      real*8 f (0:nc, 0:ny), w (0:nc, 0:ny),
     &          itg (0:ny),romw (0:2**maxpow, -1:maxpow)

      common / rom / romw
      save
      
C                              nx = 2^px
      px = nint(log (dfloat(nx)) / log (2.0))

        do y = 0, ny
C                              One dimensional integration
          itg (y) = 0.0
          do x = 0, nx
              itg (y) = itg(y)+f(x, y)*w(x, y)*romw(x, px)
          enddo
      enddo
      return
        end
C========================================================================
C
C
      SUBROUTINE romberg_initialisation (nmax)
C
C
C========================================================================
C                                                      
C      In the initialisation procedure r and index are calculated      
C_______________________________________________________________________
C      The argument is :                                    
C      nmax      The maximum number of grid points in any direction      
C_______________________________________________________________________

C**      implicit none
c        implicit DoublePrecision (a-h,o-z), integer (i-n)

        implicit real*8 (a-h,o-z)
      parameter (maxpow = 14)
      integer pmax
c      DoublePrecision r (0:2**maxpow, -1:maxpow), a (0:maxpow), 
c     &                  b (0:maxpow), s (0:maxpow, 0:maxpow)
      real*8 r (0:2**maxpow, -1:maxpow), a (0:maxpow), 
     &                  b (0:maxpow), s (0:maxpow, 0:maxpow)

C      The main variables are :
C
C      r (i,j)      Romberg weight for point # i if the integration grid has
C            2^j+1 points. For a grid with one point, j = -1.
C      a, b      Coefficients in Richardson's extrapolation

      common / rom / r
      save
      
      pmax = nint(log (dfloat(nmax))/log(2.0d0))
C                              Set up extrapolation constants
      do k = 1, pmax
          b (k) = -1.0 / (4.0 ** k - 1.0)
          a (k) = - b (k) * 4.0 ** k
      enddo
      do i = 1, pmax
          s (i,0) = 0.0
      enddo
      do i = 0, pmax
          do j = 0, 2 ** i
            r (j, i) = 0.
          enddo
      enddo
C                        Calculate weight of Integral (#grid=2^k)
      do k = 0, pmax
C                              Select Integ (#2^k)
          s (k,0) = 1.0
C                              Apply Romberg procedure
          do j = 1, pmax
              do i = pmax, j, -1
                s (i,j) = a(j)*s(i,j-1) + b(j)*s(i-1,j-1)
              enddo
          enddo
C                  s (i,i) is the weight of Integ(#2^k) if #grid = 2^1
C              2 ** (i-k) is weight of Sum(#2^k) relative to Sum(#2^i)
          do i = k, pmax
            do j = 0, 2 ** k
                r(2**(i-k)*j,i)=s(i,i)*2**(i-k)+r(2**(i-k)*j,i)
            enddo
          enddo
          s (k,0) = 0.0
      enddo
      r (0, -1) = 1.0
C                            Weights at the edge have to be halved
      do i = 0, pmax
          r (0, i) = 0.5 * r (0, i)
          r (2**i, i) = 0.5 * r (2**i, i)
      enddo
       return
      end
C===============================================================================
C
C
c      DoublePrecision FUNCTION Romberg (f, w, nc, nx, ny)             
      double precision FUNCTION Romberg (f, w, nc, nx, ny)             
C
C
C===============================================================================
C                                                      
C      The following procedure performs a two-dimensional integral     
C      using first a trapezoidal rule integral and then applying the      
C      Romberg integration technique.                              
C      nx and ny have to be powers of two. For a one dimensional      
C      integral, ny = 0. First call Romberg_Initialisation before      
C      using the Romberg function                              
C
C       Author  Frank Robijn                                  
C      Version August 17, 1989                                    
C_______________________________________________________________________
C      The arguments are :                                    
C                                                      
C      f      Integrand values                              
C      w      Relative weights of grid points                        
C      nc      Leading dimension of f a w arrays                  
C      nx      Number of points - 1 in x direction                  
C      ny      Number of points - 1 in y direction                  
C_______________________________________________________________________

C**      implicit none
c        implicit DoublePrecision (a-h,o-z), integer (i-n)

        implicit real*8 (a-h,o-z)
      integer maxpow, nc
      parameter (maxpow = 14)
      integer nx, ny, px, py, x, y
c      DoublePrecision f(0:nc, 0:ny), w(0:nc, 0:ny)
c      DoublePrecisionrom w (0:2**maxpow, -1:maxpow), 
c     &                     integral
      real*8 f(0:nc, 0:ny), w(0:nc, 0:ny)
      real*8 romw(0:2**maxpow, -1:maxpow), 
     &                     integral

C      romw (x, px) is the Romberg weight for point x in a grid of 2^px

      common / rom / romw

C                              nx = 2^px, ny = 2^py
      px = nint(log (dfloat(nx)) / log (2.0))
      if (ny.gt.0) then
           py = nint(log (dfloat(ny)) / log (2.0))
      else
C                              One dimensional integral
           py = -1
      endif
      integral = 0.0
C                        Calculate integral
      do y = 0, ny
        do x = 0, nx
          integral = integral + f (x, y) * w (x, y) *
     &                  romw (x, px) * romw (y, py)
        enddo
      enddo
      romberg = integral
      return
      end
