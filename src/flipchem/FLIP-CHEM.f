C:::::::::::::::::::::::::::: CHEMION :::::::::::::::::::::::::::
C... This routine was written by Phil Richards in August 2009
C... It returns O+, O2+, NO+, N2+, N+, NO, and N(2D) densities.These
C... densities generally agree well with the FLIP model densities.
C... The electron density must be specified as an input parameter.
C... All the densities except O+ are calculated from chemical equilibrium.
C... O+ is calculated using a Newton iterative procedure so that the total 
C... ion density matches the input electron density.
C... N+ and NO densities can either be user specified or calculated by the model.
C... Because of the neglect of diffusion, this model will generally 
C... underestimate N+ at most altitudes. Similarly, NO will be very good except 
C... below about 130 km where it will be underestimated due to neglect of diffusion. 
C... There is an artificial floor on the NO density to prevent it from getting too 
C... low below 130 km.
C... If the Newton procedure fails to converge, all ions including O+ are 
C... calculated from chemical equilibrium and then normalized to reproduce 
C... the input electron density. This generally works well.
C... The Newton procedure usually works if the total calculated molecular ion 
C... densities do not approach the input electron density. Difficulties are most   
C... likely to happen below ~150 km and especially at night. A Newton solution is 
C... usually found when the O+ density is not too much smaller than the modeled 
C... molecular ion density.
C... The EUVAC model is used for solar EUV irradiances
      
    !.. Input: Turn file output on or off
    !.. Input: Day of year in the form YYYYDDD
    !.. Input: Altitude(km), latitude(D), longitude(D)
    !.. Input: Magnetic and solar activity indices
    !.. Input: Electron and ion temperatures
    !.. Input: O, O2, N2, and He densities (cm-3)
    !.. Input: User specified NO density (cm-3)
    !.. Input: N4S should be 0.5*MSIS N density (cm-3)
    !.. Input: electron density (cm-3)
    !.. Input: User specified N+ density (cm-3)
    !.. Input: UT time (secs), Local time (hrs), and solar zenith angle(D)
    !.. OUTPUT: O+, O2+, NO+ densities (cm-3)
    !.. OUTPUT: N2+ and N+ densities (cm-3)
    !.. OUTPUT: NO and N(2D) density (cm-3)
    !.. OUTPUT: Newton procedure fails if INEWT=0
      SUBROUTINE CHEMION(JPRINT,
     >                  YYYYDDD,
     >          ALT,GLATD,GLOND,
     >            AP,F107,F107A,
     >                 TE,TI,TN,
     >          OXN,O2N,N2N,HEN,
     >                  USER_NO,
     >                      N4S,
     >                       NE,
     >               USER_NPLUS,
     >         UTSEC,LTHRS,SZAD,
     >     OXPLUS,O2PLUS,NOPLUS,
     >             N2PLUS,NPLUS,
     >                  NNO,N2D,
     >                    INEWT)
      IMPLICIT NONE
      !.. loop control variables
      INTEGER I,J,K,ITERS
      !.. Day of year in the form YYYYDDD
      INTEGER YYYYDDD
      !.. Switch for different rates
      INTEGER IRATS
      !.. Signifies when the Newton procedure fails
      INTEGER INEWT
      !.. Turns on printing of production and loss
      INTEGER JPRINT
      !.. Variables for Newton procedure
      INTEGER ITS,JITER
      !.. Electron and ion temperatures
      REAL TE,TN,TI
      !.. Secondary ion production rates
      REAL SECPRD(3,6)
    !.. Geophysical parameters
      REAL UTSEC,F107,F107A,ALT,GLATD,GLOND,LTHRS,SZAD,AP(7)
    !.. Measured H+, He+, O+, N2+, NO+, O2+, N+, RPA ion density
      REAL HEPLUS,OXPLUS,N2PLUS,NOPLUS,O2PLUS,NPLUS,USER_NPLUS
      !.. O2,O,N2,NO,He,N4S, user specified NO
      REAL O2N,OXN,N2N,NNO,HEN,N4S,USER_NO
      !.. Ne, N(2P),N(2D),O+(2P),O+(2D) densities
      REAL NE,N2P,N2D,OP2D,OP2P
    !.. Total (photon & photoel) production rates O+(4S),O+(2P),O+(2D),O2+
      REAL TPROD1,PDISOP,TPROD2,TPROD3,TPROD5
    !.. Total Production rates from all sources for NO+, O2+, 
      REAL TPNOP,O2PPROD
    !.. Production rates hv(e*)+N2->N+, hv+N->N+, Lyman-a -> NO+ 
      REAL DISNP,PHOTN,PLYNOP
      !.. generic PE production
      REAL PSEC
      !.. Reaction rates array
      REAL RTS(99)
      !.. N2+ total production
      REAL SECPN2PLUS,EUVN2PLUS
      !.. used in Newton solver
      REAL H,DEX,FEX(2)
      !.. Sum of the major ions
      REAL SUMIONS
      !.. Production and loss of NO
      REAL PNO,LNO,PDNOSR
      !.. N2(A) density    
      REAL N2A
      !.. FLIP N2(v) factor. Not used here
      REAL VCON
      !.. Production and loss of N(2D)
      REAL DISN2D,UVDISN,PN2D,LN2D
      !.. altitude for chemistry calculation
      REAL ALTCHEM
      !.. PE production rate of N2(A)
      REAL N2APRD
      REAL PN4S,LN4S,DISN4S
      REAL OXPLUSAVE

Cf2py intent(out) OXPLUS
Cf2py intent(out) O2PLUS
Cf2py intent(out) NOPLUS
Cf2py intent(out) N2PLUS
Cf2py intent(out) NPLUS
Cf2py intent(out) NNO
Cf2py intent(out) N2D
Cf2py intent(out) INEWT

    !.. various ionization and excitation rates by EUV and PE
      REAL EUVION,PEXCIT,PEPION,OTHPR2
      COMMON/EUVPRD/EUVION(3,12),PEXCIT(3,12),PEPION(3,12),OTHPR2(6)

      !.. initialize parameters
      DATA VCON/1.0/K/0/
      DATA PNO,LNO,PDNOSR,PLYNOP,N2A/5*0.0/
      DATA DISN2D,UVDISN/0.0,0.0/
      DATA HEPLUS/0.0/

      !.. Initial altitude for O+ for imposing chemistry
      ALTCHEM=150
      !.. Counts the number of Newton iterations
      JITER=0
      !.. N(2P) density
      N2P=0.0

      !.. Get the reaction rates
      CALL RATS(0,TE,TI,TN,RTS)

      !.. PRIMPR calculates solar EUV production rates
      CALL PRIMPR(1,ALT,OXN,N2N,O2N,HEN,SZAD*0.01745,TN,1,F107,F107A,
     >   AP,YYYYDDD,UTSEC,GLATD,GLOND,LTHRS,N4S)

      !.. Calculate secondary Production from photoelectrons
      CALL SECIPRD(ALT,SZAD,F107,F107A,AP,YYYYDDD,UTSEC,GLATD,GLOND,
     > LTHRS,TE,TN,OXN,O2N,N2N,NE,SECPRD,N2APRD)

   
      !**********  Come back here if Newton fails
 5    CONTINUE 
        HEPLUS=0.0
        OXPLUS=0.0
        N2PLUS=0.0
        NOPLUS=0.0
        O2PLUS=0.0
        NPLUS=0.0
        N2P=0.0
        N2D=0.0
        OP2D=0.0
        OP2P=0.0
        N2A=0.0

      !.. Iterate through chemistry twice in case O+ is chemical equilibrium
      DO ITERS=1,2
        !.. k counts number of iterations for printing headers in routines
        K=K+1
        !.. O+(2P) Calculate and print densities, production, loss
        !.. Photoelectron production
        PSEC=SECPRD(1,3)
        !.. Add EUV and photoelectrons
        TPROD3=EUVION(1,3)+PSEC
        CALL COP2P(JPRINT,7,K,ALT,RTS,OXN,O2N,N2N,NE,OP2P,TPROD3,PSEC
     >    ,HEPLUS,N4S,NNO,TE)

        !.. O+(2D) Calculate and print densities, production, loss
        !.. Photoelectron production
        PSEC=SECPRD(1,2)
        !.. EUV
        TPROD2=EUVION(1,2)
        CALL COP2D(JPRINT,8,K,ALT,RTS,OXN,O2N,N2N,NE,OP2D,TPROD2,OP2P
     >    ,HEPLUS,N4S,NNO,PSEC)

        !.. O+(4S) Calculate and print densities, production, loss. 
        TPROD1=EUVION(1,1)
        PDISOP=EUVION(2,4)+EUVION(2,5)+SECPRD(2,4)+SECPRD(2,5)
        CALL COP4S(JPRINT,4,K,ALT,RTS,OXN,O2N,N2N,NE,OXPLUS,TPROD1,OP2D
     >    ,OP2P,SECPRD(1,1),PDISOP,N2PLUS,N2D,NNO,1.0,HEPLUS)

        CALL CN2A(JPRINT,14,K,ALT,RTS,OXN,O2N,N2N,NE,N2A,N2APRD,0.0,
     >     0.0,0.0)

        !.. N(2D) Calculate and print densities, production, loss. 
        DISN2D=0.2*SECPRD(3,1)
        CALL CN2D(JPRINT,16,K,ALT,RTS,OXN,O2N,N2N,NOPLUS,NE,PN2D,LN2D
     >    ,N2PLUS,DISN2D,UVDISN,NPLUS,N2P,N2D,OXPLUS,NNO,N2A)
        N2D=PN2D/LN2D

        !.. N2+ Calculate and print densities, production, loss. 
        CALL CN2PLS(JPRINT,9,K,ALT,RTS,OXN,O2N,N2N,NE,N2PLUS,EUVION(3,1)
     >   ,EUVION(3,2),EUVION(3,3),SECPRD(3,1),SECPRD(3,2),SECPRD(3,3)
     >   ,OP2D,OP2P,HEPLUS,NPLUS,NNO,N4S)

        !.. N+ Calculate and print densities, production, loss. 
      !.. Note that N(2D) is turned off in N+ solution 
        DISNP=EUVION(3,4)+EUVION(3,5)+EUVION(3,6)+SECPRD(3,6)
        !.. N+ prod
        PHOTN=OTHPR2(3)
        CALL CNPLS(JPRINT,10,K,ALT,RTS,OXN,O2N,N2N,NE,DISNP,NPLUS,
     >    OXPLUS,N2D,OP2P,HEPLUS,PHOTN,O2PLUS,N4S,OP2D,N2PLUS,NNO)
        !.. User specified N+
        IF(USER_NPLUS.GT.0) NPLUS=USER_NPLUS

        !.. NO Calculate and print densities, production, loss.
        CALL CNO(JPRINT,15,K,ALT,RTS,OXN,O2N,N2N,NE,PNO,LNO
     >    ,N2D,N4S,N2P,NNO,O2PLUS,OXPLUS,OTHPR2(2),OTHPR2(1),N2A,NPLUS)
        
        !.. NO chemical equilibrium density
        NNO=PNO/LNO
        !.. Set a floor on NO density, which is needed below ~150 km at night 
        IF(NNO.LT.1.0E8*EXP((100-ALT)/20)) NNO=1.0E8*EXP((100-ALT)/20)
        !.. substitute user specified value
        IF(USER_NO.GT.1.0) NNO=USER_NO

        !.. NO+ Calculate and print densities, production, loss. 
        CALL CNOP(JPRINT,11,K,ALT,RTS,OXN,O2N,N2N,NE,TPNOP,NOPLUS,OXPLUS
     >    ,N2PLUS,O2PLUS,N4S,NNO,NPLUS,N2P,PLYNOP,VCON,N2D,OP2D)

        !.. O2+ Calculate and print densities, production, loss. 
        !.. EUV + PE production
        TPROD5=EUVION(2,1)+EUVION(2,2)+EUVION(2,3)+SECPRD(2,1)+
     >       SECPRD(2,2)+SECPRD(2,3)
        CALL CO2P(JPRINT,12,K,ALT,RTS,OXN,O2N,N2N,NE,O2PPROD
     >    ,O2PLUS,TPROD5,OXPLUS,OP2D,N2PLUS,NPLUS,N4S,NNO,OP2P)
      ENDDO

      !.. This section for chemical equilibrium densities for all species 
      !.. including O+. It is used when the Newton procedure fails to get O+
      !.. Don't bother if molecular ions greater than Ne/2
      !.. It increments ALTCHEM to force this action. The ion densities are 
      !.. normalized to the input NE. N+ is excluded in case it is user specified
      INEWT=1
      SUMIONS=OXPLUS+NOPLUS+O2PLUS+NPLUS
      IF(ALT.LT.ALTCHEM.OR.NOPLUS+O2PLUS.GT.0.85*NE) THEN
         OXPLUS=OXPLUS*NE/SUMIONS
         NOPLUS=NOPLUS*NE/SUMIONS
         O2PLUS=O2PLUS*NE/SUMIONS
         NPLUS=NPLUS*NE/SUMIONS
         INEWT=0
         RETURN
      ENDIF

      !.. Newton solution for O+ density given the electron density (NE)
      !.. Go through twice to set up the derivative (F(x+h)-F(x))/h
      !.. First O+ guess for Newton. This is important for high altitudes
      !.. because Newton may converge on first try.
      OXPLUSAVE=OXPLUS
      IF(NE-NOPLUS-O2PLUS.GT.100) OXPLUS=NE-NOPLUS-O2PLUS
      !.. first guess at night
      !IF(SZAD.GT.90) OXPLUS=NE
 9    DO ITS=1,2
        !.. increment for dn/dt
        IF(ITS.EQ.1) H=OXPLUS*0.0001
        !.. increment N
        IF(ITS.EQ.2) OXPLUS=OXPLUS+H

        !.. N+ Calculate and print densities, production, loss. 
        CALL CNPLS(JPRINT,10,K,ALT,RTS,OXN,O2N,N2N,NE,DISNP,NPLUS,
     >    OXPLUS,N2D,OP2P,HEPLUS,PHOTN,O2PLUS,N4S,OP2D,N2PLUS,NNO)
        !.. User specified N+
        IF(USER_NPLUS.GT.0) NPLUS=USER_NPLUS

        !.. N2+ Calculate and print densities, production, loss. 
        CALL CN2PLS(JPRINT,9,K,ALT,RTS,OXN,O2N,N2N,NE,N2PLUS,EUVION(3,1)
     >    ,EUVION(3,2),EUVION(3,3),SECPRD(3,1),SECPRD(3,2),SECPRD(3,3)
     >    ,OP2D,OP2P,HEPLUS,NPLUS,NNO,N4S)

        !.. NO+ Calculate and print densities, production, loss. 
        CALL CNOP(JPRINT,11,K,ALT,RTS,OXN,O2N,N2N,NE,TPNOP,NOPLUS,OXPLUS
     >   ,N2PLUS,O2PLUS,N4S,NNO,NPLUS,N2P,PLYNOP,VCON,N2D,OP2D)

        !.. O2+ Calculate and print densities, production, loss. 
        !.. EUV + PE production
        TPROD5=EUVION(2,1)+EUVION(2,2)+EUVION(2,3)+SECPRD(2,1)+
     >     SECPRD(2,2)+SECPRD(2,3)
        CALL CO2P(JPRINT,12,K,ALT,RTS,OXN,O2N,N2N,NE,O2PPROD
     >    ,O2PLUS,TPROD5,OXPLUS,OP2D,N2PLUS,NPLUS,N4S,NNO,OP2P)

        !.. N(2D) Calculate and print densities, production, loss. 
        CALL CN2D(JPRINT,16,K,ALT,RTS,OXN,O2N,N2N,NOPLUS,NE,PN2D,LN2D
     >     ,N2PLUS,DISN2D,UVDISN,NPLUS,N2P,N2D,OXPLUS,NNO,N2A)
        N2D=PN2D/LN2D

        !.. calculation of F(x) from the total ion concentration
        FEX(ITS)=OXPLUS+NOPLUS+O2PLUS+N2PLUS+NPLUS-NE
      ENDDO
        
      !.. Test for convergence and add increment to O+ if not
      !.. for stopping the iterations
      JITER=JITER+1
      DEX=(FEX(2)-FEX(1))/H
      OXPLUS=OXPLUS-H-FEX(1)/DEX

      !.. If Newton fails, go back and calculate O+ chem equilibrium
      IF(JITER.GT.6.OR.OXPLUS.LT.0.0.OR.OXPLUS.GT.1.0E7) THEN
      !.. forces chemical equilibrium densities
        ALTCHEM=ALT+1
      !.. Go back to chemical eqilibrium
        GOTO 5
      ENDIF     

      !.. Test for convergence
      SUMIONS=OXPLUS+NOPLUS+O2PLUS+N2PLUS+NPLUS
      IF(ABS((NE-SUMIONS)/NE).GT.0.05) GOTO 9

      !.. Normalize ion densities to the input electron density
      OXPLUS=OXPLUS*NE/SUMIONS
      NOPLUS=NOPLUS*NE/SUMIONS
      O2PLUS=O2PLUS*NE/SUMIONS
      NPLUS=NPLUS*NE/SUMIONS

      !.. If O+ a very minor ion, Newton solution may not be good, force  
      !.. chemical equilibrium solution instead
      IF(OXPLUS/SUMIONS.LT.0.1) THEN
        !.. forces chemical equilibrium densities
        ALTCHEM=ALT+1
      !.. Go back to chemical eqilibrium
      GOTO 5
      ENDIF

      RETURN
      END
