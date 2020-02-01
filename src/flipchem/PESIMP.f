C.........................<PESIMP.FOR>......................3 APR 92
C...... This test driver demonstrates how to call the model and 
C...... how to calculate electron heating rate and 3371 excitation rate.
      SUBROUTINE SECIPRD(ALT,SZADEG,F107,F107A,AP,YYYYDDD,UTSEC,GLATD,
     >    GLOND,LTHRS,TE,TN,OXN,O2N,N2N,XNE,SECPRD,N2APRD)
      !.. loop control variables
      INTEGER I,K,IK
      !.. Array dimensions
      INTEGER IDIM
      !.. Maximum PE energy
      INTEGER EMAX
      PARAMETER (IDIM=501)
      !.. PE flux
      REAL PEFLUX(IDIM)
      !.. Solar EUV multiplication factors
      REAL UVFAC(59)
      !.. Total electron impact cross sections
      REAL SIGIT(3)
      !.. Secondary production rates
      REAL SECPRD(3,6)
      !.. ALT = altitude (km)  { 120 -> 500 }
      REAL ALT
      !.. solar zenith angle {0 -> 90 degrees}
      REAL SZADEG
      !.. F107 = Solar 10.7 cm flux
      REAL F107, F107A
      !.. electron, neutral temperatures (K)
      REAL TE,TN
      !.. XN, O, O2, N2, neutral densities  (cm-3)
      REAL XN(3),OXN,O2N,N2N
      !.. electron density  (cm-3)
      REAL XNE
      !.. N(2D) density for N(2D) + e -> 2.5 eV
      REAL XN2D
      !.. O+(2D) density for O+(2D) + e -> 3.3 eV
      REAL XOP2D
      !.. PE mult. factors for bins 1-9 Torr et al.
      REAL EUV(9)
      REAL SPRD(3,6)
      !.. Total exciation cross sections for O, N2, O2
      REAL SIGOX,SIGN2,SIGEE
      !.. Production of N2A
      REAL N2APRD
      REAL AP(7)
	INTEGER YYYYDDD
	REAL UTSEC,GLATD,GLOND,LTHRS

      DATA SPRD/.4,.56,.44, .4,.28,.44, .2,.06,.10, 0.,.05,.00, 0.,.05
     >             ,.00, 0.0,0.0,0.02/

      !.. Transfer neutral densities to the density array
      XN(1)=OXN
      XN(2)=O2N
      XN(3)=N2N
      N2APRD=0.0
      DO K=1,3
      DO IK=1,6
        SECPRD(K,IK)=1.0E-15
      ENDDO
      ENDDO

      !.. Cannot calculate PE if no densities
      IF(OXN.LT.1.0E5.OR.N2N.LT.1.0E5) RETURN
      IF(SZADEG.GT.105) RETURN

      !-- A crude scaling of solar EUV with F10.7 that works pretty well
      UVFAC(58)=-1.0
      IF(F107.GT.60.0) THEN
        CALL FACEUV(UVFAC,F107,F107A)
         DO 10 I=1,9
           !.. EUV is used for PEs
           EUV(I)=UVFAC(I)
 10      CONTINUE
      ENDIF

      !********************************************************************
      !.. Go and get the photoelectron fluxes
       XN2D=0        !.. N(2D) density for calculating N(2D) + e -> 2.5 eV
       XOP2D=0       !.. O+(2D) density for calculating O+(2D) + e -> 3.3 eV
       CALL FLXCAL(IDIM,ALT,SZADEG,F107,F107A,AP,YYYYDDD,UTSEC,GLATD,
     >    GLOND,LTHRS,TE,TN,EUV,XN,XNE,XN2D,XOP2D,PEFLUX,AFAC,EMAX)
      !***************************************************************

      !........ sample calculation of ion production rates. 
      DO I=1,EMAX
        E=I-0.5
        !.. total ion cross sections
        CALL TXSION(E,SIGIT)
        CALL SIGEXS(E,TE,XNE,SIGOX,SIGN2,SIGEE)

        IF(E.LT.250) N2APRD=N2APRD+0.22*PEFLUX(I)*SIGN2*XN(3)

        !.. Evaluate ionization branching ratios for O+
        CALL OXRAT(E,SPRD(1,1),SPRD(1,2),SPRD(1,3))

        !.. Calculate ion production rates
        DO K=1,3
        DO IK=1,6
          SECPRD(K,IK)=SECPRD(K,IK)+PEFLUX(I)*SIGIT(K)*XN(K)*SPRD(K,IK)
        ENDDO
        ENDDO

      EP=E+17
      PEFLX=PEFLUX(I)/12.57
C      WRITE(18,'(A,2F8.1,1P,22E10.2)') 'PESIMP',E,12398/EP,PEFLX,
C     >     PEFLUX(I),(SIGIT(K),K=1,3),T_XS_OX(EP),2.2*T_XS_OX(EP),
C     > T_XS_N2(EP)
      ENDDO

      RETURN
      END
C:::::::::::::::::::::::::: PHOTOELECTRON MODEL  ::::::::::::::::::::::::
C....... This subroutine evaluates the photoelectron flux using the concept
C.......  production frequencies developed by Phil Richards at Utah 
C....... State University March 1984. It supercedes the model described in
C....... JGR, p2155, 1983. Contact EAST::CSPARA::RICHARDS on SPAN network
C------- Some minor updates in April 1992 indicated by C----
C....... I would appreciate any feedback on bugs or clarity and if it 
C....... contributes substantially to a paper, I would appreciate the 
C....... appropriate acknowledgement.
C......       **************** WARNING ****************
C...... This program is constructed to produce reasonable agreement with
C...... the Atmosphere Explorer-E PES fluxes of John Doering (Lee et al.
C...... PSS 1980, page 947). It will NOT give good fluxes if the EUV 
C...... attenuation is greater than about a factor of 7 (AFAC < 0.14).
C...... The model accurately reproduces the measured fluxes very closely
C...... for the case in the test driver at 148 km SZA=53 when AFAC=0.19.
C...... You should compare the output against the Lee et al. 1980 fluxes
C...... periodically as a check. It is doubtful below 140km during the
C...... day and below 200km near sunset. Between 200km & 350km, it should
C...... be good for solar zenith angles < 90 degrees. Above 350 km there
C...... is considerable uncertainty due to neglect of transport but most
C...... models have similar uncertainties at high altitudes due to the 
C...... uncertainty in the conjugate photoelectron flux, and the pitch 
C...... angle distribution.
C
C------ ALT = altitude (km)  { 120 -> 500 }
C------ SZADEG = solar zenith angle  {0 -> 90 degrees ? }
C------ TE, TN = electron, neutral temperatures (K)
C------ EUV = multiplication factors for bins 1-9 Torr et al. solar EUV
C------ XN, XNE = O, O2, N2, and electron densities  (cm-3)
C------ XN2D, XOP2D = N(2D) and O+(2D) densities for electron quenching
C------ (cm-3). You may put these to ZERO if not available.
C------ PEFLUX = photoelectron flux to be returned (eV cm2 sec)-1
C------ AFAC = the solar EUV attenuation warning flag
      SUBROUTINE FLXCAL(IDIM,ALT,SZADEG,F107,F107A,AP,YYYYDDD,UTSEC,
     >  GLATD,GLOND,LTHRS,TE,TN,EUV,XN,XNE,XN2D,XOP2D,PEFLUX,AFAC,EMAX)

      INTEGER RDIM
      INTEGER EMAX
      PARAMETER (RDIM=501)
      REAL AP(7)
	INTEGER YYYYDDD
	REAL UTSEC,GLATD,GLOND,LTHRS
      DIMENSION RJOX(RDIM),RJN2(RDIM),XN(3),COLUM(3),EUV(9),PEFLUX(IDIM)
C
C....... photoelectron production frequencies by 1.0E9. Renormalized below
C------- Note that the EUV fluxes below 250A are doubled (see refs)
      DATA RJOX/10*19,15,18,14,10,13,9,13,9,7,11,6,26,6,31,6,5,22,4,4,5
     > ,3,5.4,3.4,3.4,5,2.9,2.5,3.2,2.3,1.9,1.8,1.8,1.8,1.5,2.6,1.5,2.5
     > ,2.8,2,2.6,2.1,3.2,1.3,2.5,1.5,1.8,1.3,.3,1,.4,4*.2,.3,.2,.3,.1
     > ,.2,.2,.1,.1,.1,.2,.2, 26*.1,100*.05,100*0.04,100*0.01,100*0.02/

      DATA RJN2/6*40,43,35,35,28,29,21,25,19,19,13,19,16,12,11,7,18,8,46
     > ,27,5*5,4.3,7.4,5.6,4.3,5.1,4.3,2.8,2.7,2.7,2.1,2.1,1.7,1.6,1.3
     > ,2.5,2,2.1,2.6,2.4,2,1.3,2.2,1.6,2,1,1.4,1.1,.5,4*.05,6*.05
     > ,33*.05,100*0.02,100*0.02,100*0.01,100*0.01/

      !.. convert solar zenith angle to radians
      SZA = SZADEG/57.29578

      !.. minimum photoelectron energy
      EMIN=1
      !.. maximum photoelectron energy. 300 is good
      EMAX=280
      !..  check upper and lower energy indices in range
      IF(EMAX.GT.RDIM) EMAX=RDIM
      IF(EMIN.LT.1) EMIN=1
C
C----- 2.5eV production from electron quenching of N2D
      PN2D=XN2D*XNE*6.0E-10*SQRT(TE/300.0)
C----- 3.3eV production from electron quenching of O+(2D)
      POP2D=XOP2D*XNE*6.6E-8*SQRT(300./TE)
C
C------ Initialize electron flux
      DO 122 IE=EMIN,EMAX
        PEFLUX(IE)=0.0
 122  CONTINUE
      CASEL=0.0

C...... evaluate the neutral column density  &&&&&&&&
      CALL SCOLUM(I,SZA,ALT*1.0E5,TN,XN,YYYYDDD,UTSEC,GLATD,GLOND,
     >   LTHRS,F107A,F107,AP,COLUM)
C
C.......... begin flux calculation loop............................
      DO 133 IE=1,EMAX
      I=EMAX+1-IE
      IF(I.LT.EMIN) GO TO 55
C
C....... evaluate energy of photon responsible for electron at energy EE
      EE=I-0.5
      EP=EE+17
      IF(EE.LT.22) EP=45
      IF(EE.GE.22.AND.EE.LT.28) EP=41
      IF(EE.GE.28.AND.EE.LT.38) EP=49
C
C..... evaluate total photoionization cross sections for photon energy EP
      !.. New OX cross section
      XSOXT=T_XS_OX(EP)
      !.. O2 XS is 2.2* O XS
      XSO2T=2.2*T_XS_OX(EP)
      !.. New N2 cross section
      XSN2T=T_XS_N2(EP)

C....... evaluate EUV attenuation factor AFAC
      TAU=COLUM(1)*XSOXT+COLUM(2)*XSO2T+COLUM(3)*XSN2T
      AFAC=EXP(-TAU)
C
C......... low energy cascade production from O(1D) and N2* impact
      CASOX=0.0
      IF(EE.LT.10) CASOX=PEFLUX(I+2)*SIGOX*XN(1)
      CASN2=0.0
      IF(EE.LT.6) CASN2=PEFLUX(I+1)*SIGN2*XN(3)
C
C......... cascade production from thermal electron degradation
      CASEL=0.0
      IF(I.LT.EMAX) CASEL=PEFLUX(I+1)*SIGEE*XNE
C
C....... Production from electron quenching of metastables
      EPN2D=0.0
      IF(NINT(EE).EQ.3) EPN2D=PN2D
      EPOP2D=0.0
      IF(NINT(EE).EQ.4) EPOP2D=POP2D
C
C........ evaluate cross sections (must be after cascade production)
      CALL SIGEXS(EE,TE,XNE,SIGOX,SIGN2,SIGEE)
C
C......... adjust production rate for different period of solar cycle
      CALL FACFLX(EE,EUV,FFAC)
C
C......... Production of pe's at energy EE, taking into account
C......... attenuation and EUV variation, and renormalize frequencies
C
      PRODOX=RJOX(I)*XN(1)*AFAC*FFAC*1.0E-9 
      PRODN2=RJN2(I)*XN(3)*AFAC*FFAC*1.0E-9 
C
C......... Sum all the production rates
      PROD=PRODOX+PRODN2+CASEL+CASOX+CASN2+EPN2D+EPOP2D
C
C....      WRITE(3,90) EE,PRODOX,PRODN2,CASEL,CASOX,CASN2,EPN2D,EPOP2D
 90   FORMAT(1X,F6.1,1P,11E8.1)
C
C........ total loss through collisions
      RLOSS=SIGOX*XN(1)+SIGN2*XN(3)+SIGEE*XNE
C
C........... evaluate photoelectron flux
      PEFLUX(I)=PROD/RLOSS
 133   CONTINUE
C
 55   RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE FACFLX(EE,EUV,FFAC)
C....... solar EUV factors. Correspond to the first 9 wavelengths
C....... TORR et al.[1979] GRL page 771 table 3. EUV(9) is for 304A
      DIMENSION EUV(9)
      FFAC=(7*EUV(9)+EUV(8)+0.2*EUV(6))/8.2
      IF(EE.GT.30.AND.EE.LE.38) FFAC=(2*EUV(7)+.5*EUV(5))/2.5
      IF(EE.GT.38.AND.EE.LE.45) FFAC=EUV(4)
      IF(EE.GT.45.AND.EE.LE.66) FFAC=EUV(3)
      IF(EE.GT.66.AND.EE.LE.108) FFAC=EUV(2)
      IF(EE.GT.108) FFAC=EUV(1)
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE SIGEXS(E,TE,XNE,SIGOX,SIGN2,SIGEE)
C..... Program for evaluating the total inelastic cross sections
C
C........ loss to thermal electrons ....
      ET=8.618E-5*TE
      SIGEE=(3.37E-12/E**0.94/XNE**0.03)*((E-ET)/
     >   (E-(0.53*ET)))**2.36
C
C...... cross section for o(1d)
      SIGO1D=0.0
      IF(E.GT.1.96) SIGO1D=4E-16*(1-1.96/E)**2/E
C...... total excitation cross section for O excluding O(1D)
      IF(E.LT.25) SIGO=(0.4*E-5)*1.4E-17
      IF(E.GE.25) SIGO=7.0E-17
      IF(SIGO.LT.0.0) SIGO=0.0
C
C...... total excitation cross section for N2......
      IF(E.LT.12) SIGN2=(15.5*E-104.8)*1.7E-18
      IF(E.LT.4.0) SIGN2=5.0E-9*(1-1.4/E)**9 * (1.4/E)**16
      IF(E.GT.11.5) SIGN2=1.4E-16
      IF(SIGN2.LT.0.0) SIGN2=0.0
C
C........ total ionization cross sections from Keiffer and Dunn ....
      SIGION=0.0
      AL=ALOG10(E)
      IF(AL.LT.2.7.AND.AL.GE.1.2) SIGION=-3.6E-16*(AL-1.2)*(AL-3)
      IF(AL.GT.2.7) SIGION=1.2E-14*EXP(-AL*1.6)
      IF(E.LT.50) SIGION=1.0E-16*(0.068*E-1.06)
      IF(SIGION.LE.0.0) SIGION=0.0
C
      SIGOX=SIGO1D+SIGO+0.5*SIGION
      SIGN2=SIGN2+SIGION
      RETURN
      END
C::::::::::::::::::::::: TXSION ::::::::::::::::::::::::::::::::::
C..... total ionization cross sections for O, O2, and N2
C..... ionization cross sections keiffer and dunn ........
C..... The N2+ and O2+ cross sections were modified in April 99 to agree
C..... with the Schram et al. cross sections at high energies
      SUBROUTINE TXSION(E,SIGIT)
      DIMENSION SIGIT(3)

      !... SIGTMP is used for N2+ and O2+ at the high energies
      SIGTMP=1.0E-13*EXP(-2.303*ALOG10(E))

      !... N2+ cross section
      SIGIT(3)=0.0
      IF(E.GT.15.0) SIGIT(3)=1.42E-14*(1-9.0/E)**7.1*E**(-0.7) 
      IF(SIGTMP.LT.SIGIT(3)) SIGIT(3)=SIGTMP
      !... This correction to convert units to cm**2. Keiffer and Dunn page 10
      SIGIT(3)=0.87972*SIGIT(3)

      !... O2+ cross section
      SIGIT(2)=0.0
      IF(E.GT.12.0) SIGIT(2)=1.08E-14*(1-7.0/E)**8.6*E**(-0.65)
      IF(SIGTMP.LT.SIGIT(2)) SIGIT(2)=SIGTMP
      !... This correction to convert units to cm**2. Keiffer and Dunn page 10
      SIGIT(2)=0.87972*SIGIT(2)

      !... O+ cross section from Brook et al. J. Phys. B. Vol 11 p 3115, 1978
      SIGIT(1)=0.0
      IF(E.GT.12.0) SIGIT(1)=7.33E-15*(1-2.0/E)**34.3*E**(-0.7)
      RETURN
      END
C:::::::::::::::::::::::::::::::::::: OXRAT ::::::::::::::::::::::::::::::::
      SUBROUTINE OXRAT(E,R4S,R2D,R2P)
C....... This subroutine returns the electron impact branching ratios
C....... for atomic oxygen from Burnett and Rountree Phys. Rev. A. 20
C....... 1979 page 1468
      R4S=1.0
      R2D=0.0
      R2P=0.0
      EV=E
      IF(E.GT.100.0) EV=100.0
      IF(EV.GT.17) R4S=-1.6E-3*EV+0.56
      IF(EV.GT.17) R2D=1.067E-3*EV+0.2933
      R2P=1-R4S-R2D
         IF(EV.LT.22) THEN
         R2P=0.0
         RTOT=R4S+R2D
         R4S=R4S/RTOT
         R2D=R2D/RTOT
         ENDIF
      RETURN
      END
C::::::::::::::::::::: T_XS_N2 :::::::::::::::::::::::::::
C.... This function calculates the N2 total photoionization
C.... cross section. P. Richards 2003-10-04
      REAL FUNCTION T_XS_N2(EP)
      IMPLICIT NONE
      REAL EP   !... photon energy
      REAL ESAVE
      DATA ESAVE/0.0/

      !.. Wavelength < 20 A, Auger ionization
      IF(EP.GE.600.0) THEN              
        T_XS_N2=0.5E-18
      !.. Wavelength < 31 A, Auger ionization
      ELSEIF(EP.GE.400.0) THEN              
        T_XS_N2=1.0E-18
      !.. Wavelength 31.62 to 23.70 A
      ELSEIF(EP.GE.392.0) THEN
        T_XS_N2=EXP(7.9864*ALOG(EP)-91.6604)
      !.. Wavelength 225 to 125 A
      ELSEIF(EP.GE.55.09) THEN
        T_XS_N2=EXP(-2.3711*ALOG(EP)-29.8142) 
      !.. Wavelength > 225 A
      ELSE
        T_XS_N2=EXP(-1.1077*ALOG(EP)-34.8787)  
      ENDIF

      !..IF(NINT(10*EP).NE.NINT(10*ESAVE)) WRITE(6,'(2F8.1,1P,2E10.2)') 
      !..> 12394.224/EP,EP, T_XS_N2/(3.39E-17*EXP(-0.0263*EP)), T_XS_N2
       ESAVE=EP

      !.. old parameterization
      !..T_XS_N2=3.39E-17*EXP(-0.0263*EP)

      RETURN
      END
C::::::::::::::::::::: T_XS_OX :::::::::::::::::::::::::::
C.... This function calculates the OX total photoionization
C.... cross section. P. Richards 2003-10-04
C.... Samson and Pareek Phys. Rev. A, 31, 1470, 1985

      REAL FUNCTION T_XS_OX(EP)
      IMPLICIT NONE
      !... photon energy
      REAL EP
      REAL ESAVE
      DATA ESAVE/0.0/

      !.. NEW parameterization
      IF(EP.GE.500.0) THEN                 
        !.. Wavelength shorter than 25 A, Auger ionization
        T_XS_OX=0.5E-18
      ELSEIF(EP.GE.165.26) THEN                 
        !.. Wavelength shorter than 75 A
        T_XS_OX=EXP(-2.5209*ALOG(EP)-28.8855)
      ELSEIF(EP.GE.55.09) THEN              
        !.. Wavelength between 78 and 256.26 A
        T_XS_OX=EXP(-1.7871*ALOG(EP)-32.6335)
      ELSE
        !.. Wavelength longer than 256.26 A
        T_XS_OX=EXP(-1.3077*ALOG(EP)-34.5556)   
      ENDIF

      !..IF(NINT(10*EP).NE.NINT(10*ESAVE)) WRITE(6,'(2F8.1,1P,2E10.2)') 
      !..> 12394.224/EP,EP, T_XS_OX/(27.2E-18*EXP(-3.09E-2*EP)), T_XS_OX
       ESAVE=EP

      !.. old parameterization
      !.. T_XS_OX=27.2E-18*EXP(-3.09E-2*EP)

      RETURN
      END
