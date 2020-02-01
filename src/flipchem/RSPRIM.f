C.................... RSPRIM.FOR ..................................
C.... This routine evaluates the ionization rates for photon impact
C.... It is based on a FLIP model routine that was modified in August 
C.... 2009 for the chemical equilibrium model by P. richards. 
      SUBROUTINE PRIMPR(IJ,Z,ZOX,ZN2,ZO2,HE,SZA,TN,JMAX,F107,F107A,
     > AP,YYYYDDD,UTSEC,GLATD,GLOND,LTHRS,N4S)
      IMPLICIT REAL(A-H,O-Z)
      INTEGER IVERT,YYYYDDD
      REAL Z,ZOX,ZN2,ZO2,HE,SZA,TN,UTSEC,GLATD,GLOND,LTHRS
     > ,CHI,ZZ,TNJ,TAU,FLUX,HEPLS,FBSBN
     > ,DISN,TAUGAM,FLUXG,ALTG,XNSIGF,DSPECT,GL,N4S
      REAL COLUM(3),OTHPR1(6),OTHPR2,OTHPR3(6)
      REAL COLUMN(3),XN(3),PROB(3,6,37),XSNPLS(37)
     > ,FNITE(37),CLNITE(3),AP(7)

      !.. common to hold the EUV and photoelectron production rates 
      COMMON/EUVPRD/EUVION(3,12),PEXCIT(3,12),PEPION(3,12),OTHPR2(6)
      COMMON/SIGS/ZFLUX(37),SIGABS(3,37),ZLAM(37),SIGION(3,37),
     > TPOT(3,10),NNI(3),LAMAX
      COMMON/SOL/UVFAC(59),EUV

      LMAX=0 
      F107SV=0.0 
      IPROBS=0
      !.. Fluxes for nighttime ion production in the 37 wavelength bins of
      !.. Torr et al GRL 1979. The fluxes are set to reproduce the production
      !.. rates in Strobel et al. PSS, p1027, 1980. Note that most bins are 
      !.. set to zero and that the Strobel production rates are scaled by 
      !.. FNFAC to stabilize the O+ solution below 200 km. Note also that
      !.. the wavelengths in FNITE go from largest (#3=HI) to smallest.
      DATA FNITE/9E5,0.0,9E5,2*0.0,9E6,13*0.0,3E5,8*0.0,3E5,8*0.0/
      DATA FNFAC/1.0/

      !.. UVFAC(58) is left over from FLIP routines for compatibility
      UVFAC(58)=-1.0 
      IF(ABS((F107-F107SV)/F107).GT.0.005) THEN
        !.. update UV flux factors
        CALL FACEUV(UVFAC,F107,F107A)
        CALL FACSR(UVFAC,F107,F107A)
        !.. call params to get solar flux data and cross sections
        CALL PARAMS(0,LMAX)
        F107SV=F107

        !..  find probability for formation of each state  ........
        IF(IPROBS.EQ.0) CALL PROBS(0,PROB,ZLAM,LMAX,NNI)
        IPROBS=1

      ENDIF

      !... initialization of production rates. 1.0E-15 stabilizes 
      !... e density evaluation at low altitudes in CMINOR
      DO 10 IS=1,3
      DO 10 IK=1,12
         EUVION(IS,IK)=1.0E-15
 10   CONTINUE

      DISN=0.0
      DO 687 I=1,6
      OTHPR2(I)=1.0E-15
 687  OTHPR1(I)=1.0E-15
C
C........ Nighttime He+ production is calculated and stored. Attenuated to
C........ avoid excess production at low altitudes
      OTHPR1(2)= 8E-11* EXP(-1.0E-11*ZN2) *HE
      DO 786 I=1,3
 786  COLUM(I)=1.0E+25
      TNJ=TN
      XN(1)=ZOX
      XN(2)=ZO2
      XN(3)=ZN2
      ZZ=Z*1.0E+5
      CHI=SZA
C
C*****  obtain reaction rates from subr rats to get their densities

C....... determine if sun is below the horizon ...
C---- Must now do calculation for night production - Feb 93
      ALTG=(6371.0+Z)*SIN(3.1416-CHI)-6371.0
C....      IF(CHI.GT.1.57.AND.ALTG.LT.85.) RETURN
      IF(Z.GT.1500) RETURN
C
C...... get column densities for scattered light at night  &&&&&&&&
      CALL SCOLUM(IJ,0.0E0,ZZ,TNJ,XN,YYYYDDD,UTSEC,GLATD,GLOND,LTHRS,
     >   F107A,F107,AP,CLNITE)
C
C...... evaluate the neutral column density  &&&&&&&&
      CALL SCOLUM(IJ,CHI,ZZ,TNJ,XN,YYYYDDD,UTSEC,GLATD,GLOND,LTHRS,
     >   F107A,F107,AP,COLUMN)
C........ Store the column densities for the 2-Stream program
      COLUM(1)=COLUMN(1)
      COLUM(2)=COLUMN(2)
      COLUM(3)=COLUMN(3)
C
C........ O2 dissociation by Schumann-Runge UV.
C........ OTHPR1(3)= dissociation rate. OTHPR1(5)= Energy
      CALL SCHUMN(IJ,Z,ZO2,COLUMN,OTHPR1(3),OTHPR1(5))
C
C---- Calculate hv + NO ion. freq. from Lyman-a (Brasseur & Solomon)
C---- OTHPR2(2) is photodissociation of NO in the SR bands. 
C---- A small night production from scattered light is included. FREQLY
C---- varies with solar activity using Richards et al. 1994 page 8981
C---- LY_a=2.5E11 (Lean), sigi(NO)=2.0E-18 (Brasseur & Solomon page 329)
      DATA O2LYXS,O2SRXS,FREQSR /1.0E-20,1.0E-21,5.0E-6/
       FREQLY=5.0E-7*(1+4.0E-3*(0.5*(F107+F107A)-80.0))
       OTHPR2(1)=FREQLY*(EXP(-O2LYXS*COLUMN(2))
     >    +0.001*EXP(-O2LYXS*CLNITE(2)))
       OTHPR2(2)=FREQSR*(EXP(-O2SRXS*COLUMN(2))
     >    +0.001*EXP(-O2SRXS*CLNITE(2)))
C
      !..  wavelength loop begins here  ----------
      !..  TAU, TAUN = optical depth for day, night 
      HEPLS=0.0
      DO 6 L=1,LMAX
        TAU=0.
        TAUN=0.0
        DO I=1,3
          TAUN=TAUN+SIGABS(I,L)*CLNITE(I)
          TAU=TAU+SIGABS(I,L)*COLUMN(I)
        ENDDO

        !.. evaluate nighttime flux and daytime flux
        FLUXN=FNFAC*(F107/75.0)*FNITE(L)*EXP(-TAUN)
        FLUX=ZFLUX(L)*EXP(-TAU) + FLUXN
        !..WRITE(9,'(I6,1P,22E10.2)') L,COLUMN(1),COLUMN(3),TAU,EXP(-TAU),
        !..>    FLUX, ZFLUX(L),FLUXN

        !.. he+ production. He+ X-S  = 0.25 N2  X-S. HEPRDN = nite He+
        IF(ZLAM(L).LT.500.) HEPLS=HEPLS+HE*0.25*SIGION(3,L)*FLUX

        !.. hv + N -> N+ + e. ion. freq. Richards et al. JGR 1994 page 8989
        DATA XSNPLS/6*0.0,.211,10.294,11.171,10.961,11.244,11.323,12.098
     >  ,13.265,12.423,11.951,11.212,11.798,11.758,11.778,11.772,11.503
     >  ,11.016,10.578,9.556,8.15,8.302,7.298,6.413,6.399,5.192,5.725
     >  ,4.787,3.778,2.3,.878,.286/

        OTHPR2(3)=OTHPR2(3)+XSNPLS(L)*1.0E-18*FLUX*N4S

        IF(ZLAM(L).GE.600.0) THEN
          !...... calculation of total euv absorption-ionization.....
          FBSBN=FLUX*(SIGABS(3,L)-SIGION(3,L))*XN(3)
          !.. Save energy absorbed in the photodissociative process
          OTHPR1(4)=OTHPR1(4)+1.24E+4*FBSBN/ZLAM(L)
          !.. production on atomic nitrogen by dissociation
          DISN=DISN+FBSBN
          !..      IF(J.EQ.1) WRITE(6,95) L,ZLAM(L),TAU,FLUX,FBSBN,DISN,HEPLS
          !95   FORMAT(I4,F9.1,1P,22E9.1)
          !.. take into account the large n2 absorption of lyman gamma(972.54)
          IF(NINT(ZLAM(L)).EQ.975) THEN
             TAUGAM=370E-18*COLUMN(3)
             IF(TAUGAM.GT.70.0)TAUGAM=70.0
             FLUXG=UVFAC(34) *0.82E+9 *EXP(-TAUGAM)
             DISN=DISN+FLUXG*370E-18*XN(3)
           ENDIF
         ENDIF

        !***** species loop begins here *****
        DO 304 I=1,3
          XNSIGF=XN(I)*SIGION(I,L)*FLUX
          K1=NNI(I)
 
          !.. dspect=# ions formed by w-l l by ionization of k state of species i
          DO 302 K=1,K1
            DSPECT=XNSIGF*PROB(I,K,L)
            !.. store ion production rates .....
            EUVION(I,K)=EUVION(I,K)+DSPECT

            !.. calculation of ion heating rate......
            EUVION(1,10)=EUVION(1,10)+DSPECT*TPOT(I,K)

 302      CONTINUE
 304    CONTINUE
 6    CONTINUE
C
      !..---   wavelength loop ends here   -----------
C
C.........Store UV disoc of N2 2 atoms produced for every dissociation
      OTHPR1(1)=2.0*DISN
C........ Transfer He+ production to storage
      OTHPR1(2)=OTHPR1(2)+HEPLS
C
 777  RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE SCOLUM(J,CHI,Z,TN,XN,YYYYDDD,UTSEC,GLATD,GLOND,LTHRS,
     >  F107A,F107,AP,COLUMN)
      IMPLICIT REAL(A-H,O-Z)
      INTEGER YYYYDDD
      REAL ZG,CHI,Z,TNJ,XN,COLUMN,ALTG,GE,GR,RE,RP,
     >  SH,XP,Y,ERFY2,CHAPFN,RG,HG,XG,EM,F,G,A,B,C,D,GL
      REAL UTSEC,GLATD,GLOND,LTHRS,F107,F107A,AP(7),CHIN
C++++ this routine evaluates the neutral column density
C++++ see smith & smith jgr 1972 p 3592, subr ambs is the msis
C++++ neutral atmosphere model for z<120 km, chi=solar zenith
C++++ angle, re & ge radius and grav con for earth
      DIMENSION XN(3),COLUMN(3),SN(5),M(3),DG(9),T(2)
      DATA A,B,C,D,F,G/1.0606963,0.55643831,1.0619896,1.7245609
     >  ,0.56498823,0.06651874/
      DATA SN/0.0,0.0,0.0,0.0,0.0/
        DATA    EM   ,   M(1) , M(2) , M(3) ,  RE   , GE
     1    / 1.662E-24 ,   16. ,  32. ,  28. ,6.357E8, 980/
      DATA T,ALTG,ERFY2/0.0,0.0,0.0D0,0.0D0/
      DATA DG/9*0.0/

      DO I=1,3
        SN(I)=0.0
        COLUMN(I)=1.E+25
      ENDDO

      !.. Avoids changing Tn at grazing incidence
      TNJ=TN

      !.. below lower boundary of FLIP
      IF(Z.LT.70*1.0E5) RETURN

      !.. is sza>90.0 degrees
      !IF(CHI.GE.1.5708) RETURN
      !.. is sza>90.0 degrees
      IF(CHI.LT.1.5708) GO  TO 2938

      !..Grazing incidence parameters not calculated in this version
      ALTG=(6371.0E5+Z)*SIN(3.1416-CHI)-6371.0E5
      IF(ALTG.GE.85*1.0E5) THEN
        ZG=ALTG*1.E-5
        !..  neutral density at grazing incidence
        CALL GTD7(YYYYDDD,UTSEC,ZG,GLATD,GLOND,LTHRS,F107A,F107,AP,48,
     >  DG,T)
        SN(1)=DG(2)
        SN(3)=DG(3)
        SN(2)=DG(4)
        TNJ=T(2)
      ELSE
        RETURN
      ENDIF
      !.. sn(1)=o , sn(2)=o2 , sn(3)=n2 , tnj=tn,  gr=gravity, rp=distance
      !.. to pt p, sh=scale height, rg=distance to pt g, hg=scale height at g

 2938 CONTINUE
      GR=GE*(RE/(RE+Z))**2
      RP=RE+Z
      DO 10 I=1,3
      SH=(1.38E-16*TNJ)/(EM*M(I)*GR)
      XP=RP/SH
      Y=SQRT(0.5*XP)*ABS(COS(CHI))
      IF(Y.GT.100.0) WRITE(6,100) I,Z/1.0E5,CHI*57.3,TNJ,EM,M(I),GR,RP
  100 FORMAT('WARNING, Y IN COLUMN(I) > 100',I4,1P,9E10.2)
      IF(Y.GT.8) ERFY2=F/(G+Y)
      IF(Y.LT.8) ERFY2=(A+B*Y)/(C+D*Y+Y*Y)
    4 IF(CHI.GT.1.5708)GO  TO 2
      CHAPFN=SQRT(0.5*3.1416*XP)*ERFY2

      COLUMN(I)=XN(I)*SH*CHAPFN
        GO TO 10
    2 RG=RP*SIN(3.1416-CHI)
      HG=1.38E-16*TNJ/
     1    (EM*M(I)*980.*(6371.E5/(6371.E5+ALTG))**2)
      XG=RG/HG
      COLUMN(I)=SQRT(0.5*3.1416*XG)*HG*(2.0*SN(I)-XN(I)*ERFY2)
10       CONTINUE
      RETURN
      END
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE PARAMS(ISW,LMAX)
C........ this program determines the cross sections, solar fluxes, and
C........ given in m. torr et al g.r.l 1979 p771, table 2 and 3. but
C........ the longer wavelengths come first
      COMMON/SIGS/ZFLUX(37),SIGABS(3,37),ZLAM(37),SIGION(3,37),
     > TPOT(3,10),NNI(3),LAMAX
      COMMON/SOL/UVFAC(59),EUV
      DIMENSION  X1(111),X2(111),X3(18),ZLX(37),ZFX(37)
C
C....... ionization potentials for o,o2 and n2 see kirby et al note the
C....... o2 2(pi)u and 2(sigma-)u , and n2 f2(sigma)u pots are guesses
C....... the sixth n2 potential is for dissociation
      DATA X3/13.6,16.9,18.6,28.5,40.,0.0,12.1,16.5,18.2,20.,
     > 25.,0.,15.6,16.7,18.8,25.,29.,37./
C........ wavelength data. average is taken for bands
      DATA ZLX/1025.,1031.91,1025.72,975.,977.02,925.,875.,825.,775.,
     > 789.36,770.41,765.15,725.,703.36,675.,625.,629.73,609.76,575.,
     > 584.33,554.31,525.,475.,465.22,425.,375.,368.07,325.,303.78,
     > 303.31,275.,284.15,256.3,225.,175.,125.,75./
C........ fluxes from table 3. these are for 74113. just replace this data
C........ for other years in the table. note!!!! flux doubled for lambda<250
C........ shortest wavelenghts have been tripled
      DATA ZFX/2.4665,2.1,3.5,1.4746,4.4,3.,3.537,1.625,.758,.702,
     > .26,.17,.141,.36,.23,.342,1.59,.53,.357,1.27,.72,.452,.285,
     > .29,.383,.314,.65,.965,6.9,.8,1.679,.21,.46,3.1,4.8,.45,1.2/
C........ absorption cross sections -- o first ,o2, then n2
      DATA X1/5*0.0,1.315,4.554,3.498,5.091,3.749,3.89,4,10.736,11.46
     > ,17.245,13.365,13.4,13.4,13.024,13.09,12.59,12.059,12.127,11.93
     > ,11.496,9.687,9.84,8.693,7.7,7.68,6.461,7.08,6.05,5.202,3.732
     > ,1.839,.73,  1.346,1.0,1.63,21.108,18.73,12.817,8.562,16.631
     > ,22.145,26.668,18.91,20.8,28.535,27.44,21.919,26.017,32.06
     > ,28.07,26.61,22.79,26.04,24.606,23.101,21.91,20.31,18.118
     > ,18.32,17.438,16.81,16.8,14.387,15.79,13.37,10.9,7.509,3.806
     > ,1.316,  3*0.0,50.988,2.24,9.68,20.249,16.992,33.578,16.487
     > ,14.18,120.49,24.662,26.54,31.755,23.339,23.37,22.79,22.787
     > ,22.4,24.13,24.501,23.471,23.16,21.675,16.395,16.91,13.857
     > ,11.7,11.67,10.493,10.9,10.21,8.392,4.958,2.261,0.72/
C....... ionization cross sections 
      DATA X2/5*0.0,1.315,4.554,3.498,5.091,3.749,3.89,4,10.736,11.46
     > ,17.245,13.365,13.4,13.4,13.024,13.09,12.59,12.059,12.127,11.93
     > ,11.496,9.687,9.84,8.693,7.7,7.68,6.461,7.08,6.05,5.202,3.732
     > ,1.839,.73,  .259,0.0,1.05,13.94,15.54,9.374,5.494,6.413,10.597
     > ,10.191,8.47,11.72,23.805,23.75,21.306,24.937,31.1,26.39
     > ,26.61,22.79,26.04,24.606,23.101,21.91,20.31,18.118
     > ,18.32,17.438,16.81,16.8,14.387,15.79,13.37,10.9,7.509,3.806
     > ,1.316,  8*0.0,14.274,8.86,8.5,65.8,15.06,25.48,29.235
     > ,23.339,23.37,22.79,22.787
     > ,22.4,24.13,24.501,23.471,23.16,21.675,16.395,16.91,13.857
     > ,11.7,11.67,10.493,10.9,10.21,8.392,4.958,2.261,0.72/
C
      NNI(1)=5
      NNI(2)=5
      NNI(3)=6
      LMAX=37
      IF(ISW.NE.0) WRITE(17,95)
 95   FORMAT(/5X,'EUV fluxes, Photoabsorption, and Photoionization ',
     >  'Cross sections',
     > /4X,'I',5X,'lam',5X,'flux',4X,'sigaOX',3X,'sigaO2'
     > ,3X,'sigaN2',3X,'sigiOX',3X,'sigiO2',3X,'sigiN2',3X,'UVfac')
C
      DO 20 L=1,LMAX
      ZLAM(L)=ZLX(L)
      FFAC=UVFAC(LMAX+1-L)
      IF(ZFX(L).LT.100) ZFLUX(L)=ZFX(L)*1.0E+9*FFAC
      !..- setting up ionization potentials
      IF(L.LE.6)THEN
         TPOT(1,L)=X3(L)
         TPOT(2,L)=X3(6+L)
         TPOT(3,L)=X3(12+L)
      ENDIF
      !..- setting up cross sections
      DO 10 IS=1,3
      IN=LMAX*(IS-1)+L
      SIGABS(IS,L)=X1(IN)*1.0E-18
      SIGION(IS,L)=X2(IN)*1.0E-18
      IF(SIGABS(IS,L).LT.SIGION(IS,L)) SIGABS(IS,L)=SIGION(IS,L)
 10   CONTINUE
C
      IF(ISW.EQ.0) GO TO 20
      WRITE(17,90) L,ZLAM(L),ZFLUX(L),(SIGABS(I,L),I=1,3)
     > ,(SIGION(I,L),I=1,3),FFAC
 20   CONTINUE
C
      IF(ISW.EQ.0) RETURN
      WRITE(17,94)
 94   FORMAT(/5X,' Ionization potentials for O, O2, N2'
     > ,/2X,'4S   2D   2P   4P   2P*  -   X2   a+A  b4   B2   dis  -'
     > ,'  X2   A2   B2   C2   F2   2s')
 60   WRITE(17,91) ((TPOT(I,J),J=1,6),I=1,3)
C
      RETURN
 90   FORMAT(1X,I4,F9.2,1P,22E9.2)
 91   FORMAT(22F5.1)
      END
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE PROBS(ISW,PROB,ZLAM,LMAX,NNI)
C.... program for finding branching ratios (probabilities for various ion
C.... and molecular states) of o,o2,n2
C.... ---refs--- m. torr et al grl 1979 page 771, kirby et al atomic data
C.... and nuclear tables 1979 23,page 63
      DIMENSION YO(37,5),PROB(3,6,37),ZLAM(37),NNI(3)
C...... coefficients of o ionization cross sections from torr et al
C..... table 2
      DATA YO/.19,.486,.952,1.311,1.539,1.77,1.628,1.92,1.925,2.259
     > ,2.559,2.523,3.073,3.34,3.394,3.421,3.65,3.92,3.62,3.61,3.88,4.25
     > ,5.128,4.89,6.739,4.0,3.89,3.749,5.091,3.498,4.554,1.315,5*0.0
     > ,.206,.529,1.171,1.762,2.138,2.62,2.325,2.842,2.849,3.446,3.936
     > ,3.883,4.896,5.37,5.459,5.427,5.67,6.02,5.91,6.17,6.29,6.159
     > ,11.453,6.57,3.997,12*0.0, .134,.345,.768,1.144,1.363,1.63
     > ,1.488,1.92,1.925,2.173,2.558,2.422,2.986,3.22,3.274,3.211,3.27
     > ,3.15,3.494,3.62,3.23,2.956,0.664,14*0.0,  .062,.163,.348,.508
     > ,.598,.71,.637,.691,.693,.815,.787,.859,.541,24*0.0, .049,.13
     > ,.278,.366,.412,.35,.383,.307,.308,28*0.0/
C
C....... production of o states from torr et al table 2 (yo array)
C....... need to reverse order of yo to correspond with lambda
      DO 10 L=1,LMAX
      LL=LMAX+1-L
      SUM=YO(LL,1)+YO(LL,2)+YO(LL,3)+YO(LL,4)+YO(LL,5)
      DO 20 I=1,5
      PROB(1,I,L)=0.0
 20   IF(SUM.NE.0.0) PROB(1,I,L)=YO(LL,I)/SUM
 10    CONTINUE
C
C....... call separate subroutines for o2 and n2 probabilities
      DO 30 L=1,LMAX
      CALL PROBO2(1,L,ZLAM(L),PROB,NNI(2))
      CALL PROBN2(1,L,ZLAM(L),PROB,NNI(3))
 30   CONTINUE
C
      IF(ISW.EQ.0) RETURN
      WRITE(17,95)
 95   FORMAT(/5X,' Photoionization branching ratios for O, O2, N2'
     > ,/3X,'Lam    4S   2D   2P   4P   2P*   -   X2   a+A  b4   B2 '
     > ,'  dis   -  X2   A2   B2   C2   F2   2s')
      DO 50 L=1,LMAX
 50   WRITE(17,90) ZLAM(L),((PROB(IS,J,L),J=1,6),IS=1,3)
 90   FORMAT(F8.2,22F5.2)
       RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         SUBROUTINE PROBN2(ISW,L,ZLAM,PROB,JPTS)
C...... the n2 probabilities are taken from kirby et al tables b and c
C...... the yield of n+ is determined first then the remaining portion
C...... of the cross section is distributed amongst the n2+ ion states
C...... (x,a,b). the dissociation yield is divided between the 3 higher
C...... energy states according to wight et al. j.phys. b, 1976
C...... the 2 other states of kirby et al are not included
C
      DIMENSION A(6),B(6),X(14),Y(14,6),PROB(3,6,37)
      DATA IPTS/14/
      DATA X/50.,210.,240.,280.,300.,332.,428.,500.,600.,660.,660.01,
     > 720.,747.,796./
      DATA Y/5*.32,.3,.46,.404,.308,.308,.308,.42,1.,1.,5*.55,
     > .52,.46,.506,2*.589,.692,.58,2*.0,5*.13,.12,.08,
     > .09,.103,.103,4*0.0,3*0.0,.05,.1,.15,.83,1.,6*.0,3*.0,
     > .3,.4,.79,.17,7*.0,3*1.,.65,.5,.06,8*.0/
     
      IPTS=14
     
C
C...... if zlam is too big set equal to x(max)
      YLAM=ZLAM
      !.. Prevent divide by zero
      IF(ZLAM.GT.X(14)) YLAM=X(14)-1
      IF(ZLAM.LT.X(1)) YLAM=X(1)+1

       YIELD=0.0
C...... determine yield of n+, and store in prob array
       CALL YLDISS(1,YLAM,YIELD)
C
       DO 10 I=1,IPTS
C kjh 6/22/92   NOTE:  I realize the following statement is strange
C   looking, but its purpose is to prevent the CRAY compiler from
C   vectorizing this loop.  (Which it does incorrectly).
      if(i.eq.25)write(6,*)' '
      IF(YLAM.GT.X(I).AND.YLAM.LE.X(I+1))  GO TO 20
 10   CONTINUE
 20   SUM=0.0
C...... fit straight line between points
      DO 30 J=1,JPTS
      A(J)=(Y(I+1,J)-Y(I,J))/(X(I+1)-X(I))
      B(J)=Y(I,J)-A(J)*X(I)
 30   CONTINUE
C...... determine probabilities of n2+ states
      DO 40 J=1,JPTS
      IF(J.LE.3) PROB(3,J,L)=(A(J)*YLAM+B(J))*(1-YIELD)
      IF(J.GT.3) PROB(3,J,L)=(A(J)*YLAM+B(J))*YIELD
      SUM=SUM+PROB(3,J,L)
 40   CONTINUE
C
      IF(SUM.EQ.0.0) RETURN
C....... normalise probabilities
      DO 50 J=1,JPTS
 50   PROB(3,J,L)=PROB(3,J,L)/SUM
      RETURN
      END
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         SUBROUTINE YLDISS(ISW,ZLAM,YIELD)
C..... determination of dissociative yield of n+, refer to kirby et al
C..... page 66 and table b
      DIMENSION X(9),Y(9)
C      DATA IPTS/9/
      DATA X/50.,210.,240.,302.,387.,477.,496.,509.,2000./
      DATA Y/.36,.36,.346,.202,.033,.041,.024,0.0,0.0/
      IPTS=9
C
C
       DO 10 I=1,IPTS
C kjh 6/22/92   NOTE:  I realize the following statement is strange
C   looking, but its purpose is to prevent the CRAY compiler from
C   vectorizing this loop.  (Which it does incorrectly).
      if(i.eq.25)write(6,*)' '
      IF(ZLAM.GE.X(I).AND.ZLAM.LT.X(I+1))  GO TO 20
 10   CONTINUE
 20   IF(ZLAM.GT.387.AND.ZLAM.LT.477) GO TO 40
C....... linear interpolation
      YIELD=(ZLAM-X(I))/(X(I+1)-X(I))*(Y(I+1)-Y(I))+Y(I)
      GO TO 30
C...... parabolic interpolation see formula page 66 kirby et al
 40   YIELD=.0329+8.13E-6*(ZLAM-442)**2
 30   RETURN
      END
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         SUBROUTINE PROBO2(ISW,L,ZLAM,PROB,JPTS)
C....... o2 branching ratios are taken from kirby et al table d
C....... columns 4 & 9 are combined. columns 5,6,7&8 are combined
      DIMENSION A(5),B(5),X(20),Y(20,5),PROB(3,6,37)
      DATA IPTS,X/20,304.,323.,454.,461.,504.,537.,556.,573.,584.,598.
     >  ,610.,637.,645.,662.,684.,704.,720.,737.,774.,1026./
      DATA Y/.365,.374,.432,.435,.384,.345,.356,.365,.306,.23,.235,.245,
     > .34,.27,.482,.675,.565,.565,1.,1.,.205,.21,.243,.245,.27,.29,
     > .23,.27,.33,.295,.385,.35,.305,.385,.518,.325,.435,.435,2*.0,
     > .125,.124,.12,.12,.126,.13,.225,.216,.21,.375,.305,.37,.33,.345,
     > 6*.0,.055,.167,.11,.105,.194,.234,.189,.149,.155,.103,.075,.036
     > ,.025,7*0.,.25,.125,.095,.95,.026,15*0./
C
C...... if zlam is too big set equal to x(max)
      YLAM=ZLAM
C...... if zlam is outside range of data values set equal to max or min
      IF(ZLAM.GT.X(20)) YLAM=X(20)
      IF(ZLAM.LE.X(1)) YLAM=X(1)+1.E-3
C
       DO 10 I=1,IPTS
C kjh 6/22/92   NOTE:  I realize the following statement is strange
C   looking, but its purpose is to prevent the CRAY compiler from
C   vectorizing this loop.  (Which it does incorrectly).
      if(i.eq.25)write(6,*)' '
      IF(YLAM.GT.X(I).AND.YLAM.LE.X(I+1))  GO TO 20
 10   CONTINUE
 20   SUM=0.0
C
      DO 30 J=1,JPTS
      A(J)=(Y(I+1,J)-Y(I,J))/(X(I+1)-X(I))
      B(J)=Y(I,J)-A(J)*X(I)
      SUM=SUM+A(J)*YLAM+B(J)
 30   CONTINUE
C
      DO 40 J=1,JPTS
      PROB(2,J,L)=(A(J)*YLAM+B(J))/SUM
 40   CONTINUE
C
      RETURN
      END
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE SCHUMN(J,Z,ZO2,COLUMN,SCHUPR,SCHUHT)
C......... production of o(1d) by schumann-runge bands
C......... The fluxes are from Torr et al. GRL 1980 p6063. Scaling is
C......... done using UVFAC which may be set according to F10.7 cm flux
C......... may be done in FACEUV
      IMPLICIT REAL(A-H,O-Z)
      REAL Z,ZO2,SCHUPR,SCHUHT,COLUMN,HSRX,FLD
      COMMON/SOL/UVFAC(59),EUV
      DIMENSION COLUMN(3),SRFLUX(8),SRXS(8),SRLAM(8)
      DATA SRFLUX/2.4,1.4,.63,.44,.33,.17,.12,.053/
      DATA SRXS/.5,1.5,3.4,6,10,13,15,12/
      DATA SRLAM/1725,1675,1625,1575,1525,1475,1425,1375/
C
C........ lmax=# of lambdas in sub. primpr: schuht=heating: schupr=o(1d) prod
      LMAX=37
C
      DO 505 LSR=1,8
C......... photoabsorption cross section
      SRXSCT=1.0E-18*SRXS(LSR)
      HSRX=SRXSCT*COLUMN(2)
      IF(HSRX.GT.70)HSRX=70
C........ attentuated solar flux
      FLD=UVFAC(LMAX+LSR)*1.E+11*SRFLUX(LSR)*EXP(-HSRX)
C....... modify for eclipse ...
      IF(J.LT.70) FLD=FLD*EUV
C............ neutral heating SCHUHT and photodissociation rate SCHUPR
      SCHUHT=SCHUHT+1.24E+4*(FLD*SRXSCT)*ZO2/SRLAM(LSR)
      SCHUPR=SCHUPR+FLD*SRXSCT
C....      IF(JTI.EQ.0) WRITE(6,90) LSR,SRXSCT,FLD,SCHUPR
 505  CONTINUE
C
      SCHUPR=ZO2*SCHUPR
 90   FORMAT(2X,I5,1P,9E9.1)
      JTI=1
      RETURN
      END
C:::::::::::::::::::::::::::::::: FACEUV :::::::::::::::::::::::
      SUBROUTINE FACEUV(UVFAC,F107,F107A)
C----- This routine uses the EUV scaling from Richards et al.[1994]
C----- The EUVAC flux model is based on the F74113 solar reference
C----- spectrum and Hinteregger's scaling factors. This subroutine
C----- just provides the scaling factors as a function of the proxy
C----- (F107+F107A)/2
      DIMENSION UVFAC(59),HFG200(37)
      DATA HFG200/2.202,1.855,2.605,3.334,1.333,17.522,4.176,4.0
     >  ,1.4,3.694,1.791,5.385,1.889,1.899,3.427,2.051,1.392,1.619
     >  ,1.439,2.941,1.399,2.416,1.512,1.365,1.570,1.462,2.537,1.393
     >  ,1.572,1.578,1.681,1.598,1.473,1.530,1.622,1.634,1.525/

      !..  Test to see if need to scale - see DATRD2 subroutine      
      IF(NINT(UVFAC(58)).EQ.-1.OR.NINT(UVFAC(58)).EQ.-3) THEN
         !........... EUV scaling
         F107AV=(F107+F107A)*0.5
         DO 50 I=1,37
         A=(HFG200(I)-1)/120.0
         B=1-A*80.0
         UVFAC(I)=A*F107AV+B
         IF(UVFAC(I).LT.0.8) UVFAC(I)=0.8
 50      CONTINUE
      ENDIF
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE FACSR(UVFAC,F107,F107A)
C........ The Schumann-Runge factors are scaled according to F10.7
C........ from Torr et al. GRL 1980 p6063
      DIMENSION UVFAC(59),SRFLUX(8),SRA(8),SRB(8)
      !............. Schumann-Runge scaling
      DATA SRFLUX/2.4,1.4,.63,.44,.33,.17,.12,.053/
      !...... first two SRA and SRB values out of order in Marsha's paper
      DATA SRA/25.5,20.7,13.2,11.6,11.3,7.86,7.68,4.56/
      DATA SRB/222.,129.,53.4,36.0,25.0,11.3,6.35,2.05/
C
C----  Test to see if need to scale - see DATRD2 subroutine      
      IF(NINT(UVFAC(58)).EQ.-2.OR.NINT(UVFAC(58)).EQ.-3) THEN
C
         DO 505 I=38,50
            LSR=I-37
            UVFAC(I)=1.0
            IF(LSR.LE.8)
     >      UVFAC(I)=(SRA(LSR)*1.0E7*F107+SRB(LSR)*1.0E9)/SRFLUX(LSR)
     >           /1.0E11
 505     CONTINUE
      ENDIF
      RETURN
      END
C::::::::::::::::::::::::: EPHEM ::::::::::::::::::::::::::::::::::::::
      SUBROUTINE EPHEM(IDAY,ETRAN)
      DIMENSION EPTRAN(26)
      DATA EPTRAN/4338,4374,4398,4404,4392,
     >  4374,4344,4320,4302,4296,4302
     > ,4320,4338,4356,4356,4344,4320
     > ,4290,4260,4236,4218,4224,4248,4284
     > ,4320,4350/
      ANUM=(MOD(IDAY,1000)+15.)/15.
      INUM=ANUM
      FRAC=ANUM-INUM
      ETRAN=(EPTRAN(INUM)*(1.0-FRAC)+EPTRAN(INUM+1)*FRAC)*10.
      RETURN
      END
C::::::::::::::::::::::::::: SOLDEC :::::::::::::::::::::::::
C...... This routine calculates the solar declination (DELTA radians)
C...... and ephemeris transit time (ETRAN in secs). ETRAN is the time
C...... the sun is overhead in Greenwich
C...... REFERENCE - page 484 "Explanatory Supplement to the Astronomical
C...... Almanac" Kenneth Seidelmann, University Science Books, 20 Edgehill 
C...... Road, Mill Valley, CA 94941, 1992
C...... IDAY is yyyyddd (eg 1988015 = feb 15), UT is the universal time
C...... in hours (0-24) 
       SUBROUTINE SOLDEC(IDAY,UT,DELTA,ETRAN)
      !.... Find day of year (0-366), and year (1996), then Julian day
       INTEGER JD,DAYNUM
       REAL T,YEAR,UT,L,G,LAMBDA,EPSIL,E,GHA,DELTA,SD,ETRAN,DTOR
       DOUBLE PRECISION DJD,DUT
       !.. degrees to radians
       DATA DTOR/57.29578/
      !..... Recover year and date and make sure UT is less than 24
       DAYNUM=MOD(IDAY,1000)
       YEAR=(IDAY-DAYNUM)/1000
       !.. Julian day
       JD=INT(365.25*(YEAR-1900)+DAYNUM)+2415020
       !.. extra precision needed         
       DJD=JD
       DUT=UT
       !.. # of centuries
       T=(DJD+DUT/24.0-2451545.0)/36525.0
      !..     WRITE(6,*) DAYNUM,YEAR,JD,UT
       !.. aberration
       L=AMOD(280.460+36000.770*T,360.0)
       !.. mean anomaly
       G=AMOD(357.528+35999.050*T,360.0)
      !...... LAMBDA= ecliptic longitude. DELTA=obliquity of the ecliptic.
       LAMBDA=AMOD(L+1.915*SIN(G/DTOR)+0.020*SIN(2.0*G/DTOR),360.0)
       EPSIL=23.4393-0.01300*T
      !...... Equation of time. Time difference between noon and overhead sun
       E=-1.915*SIN(G/DTOR)-0.020*SIN(2.0*G/DTOR)
     >   +2.466*SIN(2*LAMBDA/DTOR)-0.053*SIN(4*LAMBDA/DTOR)
       !.. Greenwich hour angle
       GHA=15*UT-180+E
       !.. solar declination
       DELTA=ASIN(SIN(EPSIL/DTOR)*SIN(LAMBDA/DTOR))
       !.. 
       SD=0.267/(1-0.017*COS(G/DTOR))
       !.. Ephemeris transit time in secs for FLIP
       ETRAN=(12-E/15)*3600
       RETURN
       END
C::::::::::::::::::::::: GETLTSZA ::::::::::::::::::::::::::::
      SUBROUTINE GETLTSZA(IDAY,SEC,GLATR,GLOND,SAT,SZA,DEC)
      IMPLICIT NONE
      !.. Day of year in form YYYYDDD
      INTEGER IDAY
      !.. UT in seconds
      REAL SEC
      !.. Geographic Latitude and Longitude
      REAL GLATR,GLOND
      !.. Solar Apparent time
      REAL SAT
      !.. Solar zenith angle & argument 
      REAL SZA,SZA_ARG
      !.. Solar declination
      REAL DEC
      !.. Hour angle
      REAL HH
      !.. ephemeris transit time
      REAL ETRAN

Cf2py intent(out) SAT
Cf2py intent(out) SZA
Cf2py intent(out) DEC

      !.. get ephemeris transit
      CALL EPHEM(IDAY,ETRAN)
      !..... Get solar declination and ephemeris transit time

      CALL SOLDEC(IDAY,SEC/3600.0,DEC,ETRAN)

      !...... Evaluate Solar apparent time
      SAT=(SEC-ETRAN+43200.0)/3600+GLOND/15.0
      IF(SAT.LT.0) SAT=SAT+24
      IF(SAT.GT.24) SAT=SAT-24

      !...... Evaluate Solar zenith angle
      !.. hour angle
      HH=(SAT-12.)*15.0*3.141592654/180.
      !.. Evaluate argument for ACOS and test because roundoff error 
      !.. may cause invalid ACOS argument      
      SZA_ARG=COS(GLATR)*COS(DEC)*COS(HH)+SIN(GLATR)*SIN(DEC)
      IF(SZA_ARG.GT.1.0) SZA_ARG=1.0
      !.. solar zenith angle
      SZA=ACOS(SZA_ARG)
      RETURN
      END
C::::::::::::::::::::::::::::::::: ACTUAL_DAY:::::::::::::::::::::::::
C...... Since IDAY does not change and UT continually increases, need to
C...... evaluate actual day (ID) and UT (SX) for MSIS, HWM, and IRI.
C...... Modified by P. Richards in November 2008 to allow for long runs
      !.... FLIP start day YYYYDDD
      !.... Current UT in secs (continuous)
      !.... actual day YYYYDDD (output) 
      !.... UT (0-24) in secs (output)
      SUBROUTINE ACTUAL_DAY(IDAY,
     >                       SEC,
     >                        ID,
     >                        SX)
      
      IMPLICIT NONE
      !.. Dimension for array
      INTEGER SDIM
      PARAMETER(SDIM=11)
      !.. Loop control, IWR= write control 
      INTEGER I,IWR
      !.. FLIP start YYYYDDD, ID= current YYYYDDD
      INTEGER IDAY,ID
      !.. Start year YYYY and start day DDD
      INTEGER YSTART,DSTART
      !.. # of years, days, and seconds in these years
      INTEGER NYRS,NDAYS,NSECS(SDIM)
      !.. # Number of days in run so far
      INTEGER IDAYS
      !.. UT in secs, actual UT in secs (output)
      REAL SEC,SX
      IWR=0

      !.. Check the run is not too long
      IF(SEC.GT.SDIM*3.16E+7) THEN
         WRITE(6,55) SDIM
         !CALL RUN_ERROR    !.. print ERROR warning in output files
         STOP
      ENDIF
 55   FORMAT(2X,'** ERROR: FLIP cannot do more than',I3,' year run **')

      !.. Determine the start year YYYY
      YSTART=IDAY/1000
      !.. Determine the start day DDD
      DSTART=MOD(IDAY,1000)

      !.. Determine the number of days and seconds simulated so far 
      !.. in order to determine the current year and the UT in that year
      NDAYS=0
      DO I=1,SDIM
        !.. leap year
        IF(MOD(YSTART+I-1,4).EQ.0) THEN
           NDAYS=NDAYS+366
        ELSE
           !.. normal year
           NDAYS=NDAYS+365
        ENDIF

        !.. Time elapsed since the start day in seconds
        NSECS(I)=(NDAYS-DSTART+1)*24*3600

        !.. Jump out of loop when the time exceeds the current FLIP run time
        IF(NSECS(I).GT.INT(SEC)) THEN
           NYRS=I-1
           GOTO 10
        ENDIF
      ENDDO

 10   CONTINUE

      !.. Determine year and UT in the first year and then later years
      IF(NYRS.LE.0) THEN
        IDAYS=INT(SEC)/3600/24
        ID=(YSTART)*1000+DSTART+IDAYS
      !.. after first year
      ELSE
        IDAYS=(INT(SEC)-NSECS(NYRS))/3600/24
        ID=(YSTART+NYRS)*1000+IDAYS+1
      ENDIF

      !.. UT 
      SX=AMOD(SEC,8.64E+4)

      !.. leap year
      IF(MOD(ID/1000,4).EQ.0) THEN
      IF(IDAYS.GT.366) THEN
          WRITE(6,*) '  ** ERROR > 366 days in leap year'
          !.. print ERROR warning in output files
          !CALL RUN_ERROR
        ENDIF
      ELSE
      IF(IDAYS.GT.365) THEN
          WRITE(6,*) '  ** ERROR > 365 days in normal year'
          !.. print ERROR warning in output files
          !CALL RUN_ERROR
        ENDIF
      ENDIF

      !..IWR+1  !.. IWR=0 switches off DEBUG write statements
      IWR=0
      IF(IWR.EQ.1) WRITE(6,191) 
 191  FORMAT(5X,'IDAY',8X,'ID    NYRS IDAYS     SX',9X,'SEC',8X,'Days')
      IF(IWR.GT.0) WRITE(6,'(2I11,2I5,9I11)') IDAY,ID,NYRS,IDAYS,
     >    INT(SX),INT(SEC),NDAYS,NSECS(NYRS)

      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C...... Empirical Te model from Brace and Theis GRL p275, 1978
      SUBROUTINE BRACE(H,NE,TN,TB)
      REAL H,NE,TN,TB,NEB
      NEB=NE
      IF(NEB.LT.1000) NEB=1000
      IF(NEB.GT.4E6)  NEB=4E6
      TB=1051+(17.07*H-2746)*
     >  EXP(-5.122E-4*H+6.094E-6*NEB-3.353E-8*H*NEB)
      IF(TB.LT.TN) TB=TN
      IF(TB.GT.4000) TB=4000
      RETURN
      END