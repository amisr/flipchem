C.................................... RATES.FOR ................. 
C.... This is the reaction rate subroutine for the FLIP model. It takes
C.... temperatures as input and puts the rates in array RTS. It includes
C.... reaction rates, efficiencies for various products, and Einstein
C.... coefficients. For a complete set of references see Fox and Sung JGR
C.... 2001, page 21,305. Rates different from Fox and Sung indicated by PGR
      SUBROUTINE RATS(J,TE,TI,TN,RTS)
      IMPLICIT REAL(A-H,L,N-Z)
      REAL TE,TI,TN,RTS,T13,TOT_NP_O2_RATE
      DIMENSION RTS(99)

      !.. zero out array
      DO 898 ITJS=1,99
 898  RTS(ITJS)=0.0

      !.. O + H+ -> O+ + H      Fox and Sung [2001]
      ZED=1+0.6*EXP(-228.0/TN)+0.2*EXP(-326.0/TN)
      RTS(1)=(8*6.4E-10/(9*ZED))
     >     *(EXP(-232.1/TI)+0.6*EXP(-228.0/TI)+0.2*EXP(-326.0/TI))

      !.. O+ + H -> O + H+    Anicich et al. [1993]
      RTS(2)=6.4E-10

      !.. O+ + N2 --> NO+ + N,   Hierl et al.[1997] 
      !.. The Hierl et al. [1997] lab rate is contaminated by N2(v) 
      !.. for T > 1300K. Therefore, the Hierl et al. rate is not really 
      !.. appropriate  in the ionosphere. The IDC model uses the Hierl et  
      !.. al. rate because it does not solve for N2(v). The FLIP model 
      !.. solves for N2(v) and uses the St. Maurice and Torr rate (JGR,1978,p969)
      !.. Hierl et al.[1997] 
      IF(TI.LE.1000) RTS(3)=1.2E-12*(300/TI)**0.45
      !.. Hierl et al.[1997] 
      IF(TI.GT.1000) RTS(3)=7.0E-13*(TI/1000)**2.12

      !.. O+ + O2 -> O2+ + O,   Lindinger et al. [1974] 
      !.. Hierl et al. lists different rates. Hierl et al. [1997] not 
      !.. used above 1600 because rates are contaminated by O2(v) for 
      !.. T > 1000K. We don't know the vibrational state in the 
      !.. thermosphere. This fit was done by PGR May 2009. It is similar 
      !.. to Fox and Sung but does not increase sharply above 1000K.
      IF(TI.LE.1600) RTS(4)=1.6E-11*(300/TI)**0.52
      IF(TI.GT.1600) RTS(4)=6.7E-12*(TI/1600)**0.6 

      !.. NO+ + e -> N + O    Walls and Dunn [1974)
      !.. Vejby-Christensen et al [1998] gives 4.0E-7*(300/TE)**0.5
      !.. Torr and Torr [1979] gives 4.3E-7*(300/TE)**0.83(-0.16,+.08)
      !.. Sheehan and St. Maurice gives 3.5E-7*(300/TE)**0.65
      RTS(5)=4.0E-7*(300/TE)**0.85

      !.. O2+ + e -> O + O   Mehr and Biondi (1969)
      IF(TE.LE.1200) RTS(6)=1.953E-7*(300/TE)**0.70
      IF(TE.GT.1200) RTS(6)=7.389E-8*(1200/TE)**0.56

      !..  O2 + N(4S)-> NO + O           Baulch et al.[1994]
      RTS(7)=1.5E-14*TN*EXP(-3270.0/TN)

      !..   N(2D) + e -> N(4S) + e     Berrington and Burke [1981]
      RTS(8)=3.86E-10*(TE/300.)**0.81

      !..   NO + N(4S) -> N2 + O      Lee et al. [1978]
      RTS(9)=3.4E-11 

      !.. N2+ + O -> NO+ + N   Scott et al.[1999]
      IF(TI.LE.1500) RTS(10)= 1.33E-10*(300/TI)**0.44
      IF(TI.GT.1500) RTS(10)= 6.55E-11*(1500/TI)**(-0.2)

      !.. N2+ + e -> N + N  Mehr and Biondi (1969)
      !.. Zipf (1980)
      RTS(11)=2.2E-7*(300/TE)**0.39

      !.. O+(2D) + e -> O+(4S) + e   McLaughlin and Bell (1998)
      !.. Henry [1969] gives 7.8E-8*(300/TE)**0.5
      RTS(12)=6.03E-8*(300/TE)**0.5

      !.. O+(2P) + e ->  O+(2D) + e   McLaughlin and Bell (1998)
      !.. RTS(13)+RTS(14) agrees with Walker et al (1975) and 
      !.. Chang et al (1993)
      RTS(13)=1.84E-7*(300/TE)**0.5

      !.. O+(2P) + e -> O+(4S) + e  McLaughlin and Bell (1998)
      RTS(14)=3.03E-8*(300/TE)**0.5

      !.. N(2D) + O ->  N(4S) + O  Fell et al.[1990]. Lin and Kaufman[1971]
      RTS(15)=6.9E-13

      !.. N(2D) + O2 -> NO + O  Herron[1999]. Shihira et al.[1994]
      RTS(16)=9.7E-12*EXP(-185.0/TN)

      !.. N2+ + O2 -> O2+ + N2   Scott et al.[1999]
      IF(TI.LT.1000) THEN
        RTS(17)=5.1E-11*(300/TI)**1.16
      ELSEIF(TI.LE.2000) THEN
        RTS(17)=1.26E-11*(TI/1000)**0.67
      ELSE
        RTS(17)=2.39E-11
      ENDIF

      !.. thermal electron excitation of O(1D); Rees et al 1967 pss, p1097 .....
      RTS(18)=1.1E-10*SQRT(TE)*EXP(-2.27E+4/TE)*(0.406
     >  +0.357E-4*TE-(0.333+0.183E-4*TE)*EXP(-1.37E4/TE)
     >  -(0.456+0.174E-4*TE)*EXP(-2.97E4/TE))

      !.. N2 + O+(2D) -> N2+ + O   
      !.. Johnson and Biondi
      !..RTS(19)=8.0E-10
      RTS(19)=1.50E-10*(300/Ti)**(-0.55)
      !.. Li et al by PGR
      
      !.. N2 + O+(2P) -> N2+ + 0    Fox 
      !.. Li et al from Fox wrong
      !.. RTS(20)=6.2E-10*EXP(-340/TI)
      !.. Li et al by PGR
      RTS(20)=2.0E-10*(300/Ti)**(-0.55)

      !.. O2+ + N(4S) -> NO+ + 0   Scott et al.[1999]
      RTS(21)=1.0E-10

      !.. N+ + O2 -> O+ + NO 
      !.. Torr and Torr gives 6.0E-10 for total N+ + O2 reaction rate
      !.. Dotan et al [1997] from Fox and Sung gives
      !IF(TI.LE.1000) TOT_NP_O2_RATE=2.02E-10*(300/TI)**-0.45
      !IF(TI.GT.1000) TOT_NP_O2_RATE=3.49E-10
      !.. does not seem to be correct. Probably vibrationally excited O2
      !.. Branching ratios for N+ + O2 from O'Keefe et al J. Chem. Phys. 1986
      !.. NO+ +O(3P) = .09, NO+ + O(1D) = .36, O2+ + N(4S) = 0.35, 
      !.. O2+ + N(2D) = 0.15, O+(4S) + NO = .05
      !.. Total N+ + O2 rate
      TOT_NP_O2_RATE=6.0E-10
      RTS(22)=0.05*TOT_NP_O2_RATE      

      !.. O2+ + NO -> NO+ + O2 Midey and Viggiano [1999]
      RTS(23)=4.5E-10 * 1.0000

      !.. O+ + NO -> O + NO+   Dotan and Viggiano [1999]
      IF(TI.LE.300) RTS(24)=7.0E-13*(300/TI)**0.66
      IF(TI.GT.300) RTS(24)=7.0E-13*(TI/300)**0.87

      !.. N+ + O2 -> O2+ + N(4S) 
      RTS(25)=0.35*TOT_NP_O2_RATE 

      !.. O+(2P) + O -> O+(4S) + O 
      !.. Fox appears to be wrong  
      !..RTS(26)=5.2E-10
      !.. (Chang et al., JGR 1993) c.f. 5.2E-11  (Rusch)    
      RTS(26)=4.0E-10

      !.. N2(A3sig) + O -> NO + N(2D) 
      !..see Campbell et al. 2006
      RTS(27)=2.0E-11
      !.. Torr and Torr value
      RTS(27)=0.000000

      !.. O+(2D) + O ->  O+(4S) + O  Torr and Torr [1980]
      RTS(28)=1.0E-11

      !.. O+ + N(2D) -> O + N+  Constantinides et al.[1979].Bates[1989]
      RTS(29)=1.3E-10

      !.. O2 + N+ -> O(3P) + NO+
      !.. Branching ratio from O'Keefe et al J. Chem. Phys. 1968
      RTS(30)=0.09*TOT_NP_O2_RATE 

      !.. O + N+ -> O+ + N   Constantinides et al.[1979].Bates[1989]
      RTS(31)=2.2E-12

      !.. Efficiency for   N2+ + e -> N(2D) + N(2D)
      RTS(32)=1.46

      !.. N2 + O(1D) -> O + NO
      RTS(33)=1.8E-11*EXP(107.0/TN)

      !.. O2 + O(1D) -> O + O2
      RTS(34)=3.2E-11*EXP(67/TN)

      !.. O2 + N(4S) -> O(1S) + NO. Kopp et al. 1977, JGR, p4715
      RTS(35)=2.5E-11


      !.. N2(A3sig) + O -> O(1S) + N2
      !..  see Campbell et al. 2006
      RTS(36)=2.5E-11*EXP(TN/298)**0.55
      !.. Torr et al.
      RTS(36)=2.0E-11

      !.. N(2P) + O -> products (N(2D,4S) and NO+) and O(3P,1D) 
      !.. from Piper et al 1993, J. Chem. Phys. vol 98 page 8560.
      RTS(37)=1.7E-11

      !.. N(2P) + O2 -> NO + O 
      RTS(38)=3.9E-12*EXP(-60/TN)

      !.. N(2P) quenching rates(O2+,NO) from Zipf et al jgr 1980 p687
      RTS(39)=2.2E-11
      RTS(40)=1.8E-10

      !.. N(2D) + NO -> N2 + O 
      RTS(41)=6.7E-11

      !.. efficiency N2+ + O -> N2 + O+(4S)   
      IF(TI.LE.1500) RTS(42)= 7.0E-12*(300/TI)**0.21
      IF(TI.GT.1500) RTS(42)= 4.83E-12*(1500/TI)**(-0.41)
      !.. converts to efficiency
      RTS(42)=RTS(42)/RTS(10)

      !.. O+(2D) + O2 -> O2+ + O   Fox
      RTS(43)=7.0E-10

      !.. He+ + N2 -> He + N2+
      RTS(44)=5.2E-10

      !.. He+ + N2 -> He + N+
      RTS(45)=7.8E-10

      !.. O(1S)+ e -> O(1D) + e  
      RTS(46)=8.5E-9

      !.. O(1S)+ e -> O(3P) + e  
      RTS(47)=1.56E-10*(TE/300)**0.94

      !.. O(1S) + O2 -> O2 + O 
      RTS(48)=4.4E-12*EXP(-815.0/TN)

      !.. NO+ + e -> N(4S) + O
      RTS(49)=0.15 * RTS(5)

      !.. NO+ + e -> N(2D) + O
      RTS(50)=0.85 * RTS(5)

      !.. O2+ + e -> O(1D) + O
      RTS(51)=1.11 * RTS(6)

      !.. O2+ + e -> O(1S) + O
      RTS(52)=0.05 * RTS(6)

      !.. Efficiency for   N2+ + e -> N(4S) + N(2D)
      RTS(53)=0.46

      !.. O(1D) -> O + 6300 + 6364 
      RTS(54)=0.00934

      !.. O(1S) -> O(1D) + 5577
      RTS(55)= 1.06

      !.. O(1S) -> O(3P) + hv (2972) RTS(56)= 4.5E-2 !.. old value
      RTS(56)= 0.10 * RTS(55)  !.. From Slanger, Spring AGU 2005

      !.. N(2P) -> N(2D) + hv
      RTS(57)=7.9E-2

      !.. N(2P) -> N(4S) + hv
      RTS(58)=5.0E-3

      !.. N+ + O2 -> NO+ + O(1S) Langford et al., PSS, 33,1225,1985
      RTS(59)=1.0E-3*TOT_NP_O2_RATE 

      !.. Efficiency for   N2(A3sig) + O -> O(1S) + N2
      RTS(60)=0.37

      !.. N(2D) -> N(4S) + hv
      RTS(61)=1.07E-5

      !.. hv(>600A) + N2 -> N(4S) + N   branching ratio
      RTS(62)=0.5

      !.. hv(>600A) + N2 -> N(2D) + N   branching ratio
      RTS(63)=0.4

      !.. hv(>600A) + N2 -> N(2P) + N   branching ratio
      RTS(64)=0.1

      !.. N+ + O2 -> O2+ + N(2D)
      !.. Branching ratio from O'Keefe et al J. Chem. Phys. 1968
      RTS(65)=0.15*TOT_NP_O2_RATE 

      !.. N+ + O2 -> NO+ + O(1D)
      !.. Branching ratio from O'Keefe et al J. Chem. Phys. 1968
      RTS(66)=0.36*TOT_NP_O2_RATE 

      !.. hv(Scum-Runge) + O2 -> O(1S) + O   branching ratio
      RTS(67)=0.001

      !.. Effic of O2(A3,DEL) + O -> O(1S)
      RTS(68)=0.1

      !.. O(1D) + O -> O + O   Abreu et al. PSS, p1143, 1986
      RTS(69)=6.47E-12*(TN/300)**0.14

      !.. hv + N2 -> N+(5S) -> 2143 A emission yield from the 2s sigma g state  
      !.. of N2. This was taken as 0.6 the value of Cleary and Barth JGR 1987, 
      !.. p13,635 because they did not double EUV below 250 A.
      RTS(70)=0.06

      !.. hv + N2 -> N+(1D) -> 6584 A emission (guess)
      RTS(71)=0.3

      !.. hv + N2 -> N+(1S) -> 5755 A emission (guess)
      RTS(72)=0.03

      !.. efficiency of production of N(2P) from e + N2+ reaction
      RTS(73)=0.08

      !.. Efficiency for production of O(1D) from N(2D) + O2 reaction
      !.. See Link and Swaminathan, PSS, 1992, page 699
      !??? check
      RTS(74)=0.1

      !.. He+ + O2 -> He + O2+
      RTS(75) = 9.2E-12

      !.. He+ + O2 -> He + O+(2D) + O(3P) 
      RTS(76) = 2.37E-10

      !.. O2+ + N(2D) -> NO+ + O
      RTS(77) = 1.8E-10

      !.. O2+ + N(2D) -> N+ + O2
      RTS(78) = 8.65E-10

      !.. N2+ + N(4S) -> N+ + N2
      RTS(79) = 1.0E-11

      !.. N2+ + NO -> NO+ + N2
      RTS(80) = 3.6E-10

      !.. N+ + NO -> N(4S) + NO+
      RTS(81) = 4.72E-10*(300/TI)**0.24

      !.. N+ + NO -> N2+ + O
      RTS(82) = 8.33E-11*(300/TI)**0.24

      !.. O+(2D) + NO -> NO+ + O
      RTS(83)=1.2E-9

      !.. O+(2D) + N -> N+ + O
      RTS(84)=1.5E-10

      !.. O+(2P) + O2 -> O+ + O2  Fox
      RTS(85)=1.3E-10

      !.. O+(2P) + O2 -> O2+ + O
      RTS(86)=1.3E-10

      !.. O+(2P) + N -> O+ + N(2D)
      RTS(87)=1.0E-11

      !.. O+(2P) + NO -> NO+ + O
      RTS(88)=1.2E-9

      !.. H+ + O2 -> O2+ + H
      RTS(89)=3.8E-9
 
      !.. O+(2D) + N2 -> NO+ + N  !.. Li et al. (1997).
      !.. From the ratio of the cross sections.
      !.. The branching ratio to O+(4S) + N2 not given by Li et al.
      RTS(90)=2.5E-11

      !.. He+ + O2 -> He + O+(4S) + O 
      RTS(91) = 2.39E-11

      !.. He+ + O2 -> He + O+(2P) + O 
      RTS(92) = 6.04E-10

      !.. He+ + O2 -> He + O+(4S) + O(1D)
      RTS(93) = 4.6E-11

      !.. He+ + NO -> He + N+ + O
      RTS(94) = 1.35E-9

      !.. He+ + NO -> He + O+ + N
      RTS(95) = 1.0E-10

      !.. N(2P) + e -> N(4S) + e
      RTS(96)=2.04E-10*(TE/300)**0.85
 
      !.. N(2P) + e -> N(2D) + e
      RTS(97)=9.5E-9

      !.. O(1D) + e -> O(3P) + e
      RTS(98)=2.87E-10*(TE/300)**0.91

      !.. N2+ + O -> O+ + N2  !.. McFarland et al.(1974)
      !.. From Fox and Sung 2001
      RTS(99)=0.07*1.0E-10*(300/Ti)**0.23

      !.. Replace some rates with those from a file
      !..CALL ALT_RATES(TE,TI,TN,RTS)

        RETURN
        END
C:::::::::::::::::::::::::::::::::::: RATCHK :::::::::::::::::::::::::::::
C.. this subroutine checks the reaction rates from subroutine rats
C.. against a standard set and prints a warning if any are different
C.. Modified in April 20008 to write warnings in multiple files
      !.. NUMU  number of files to write rate message
      !.. UNITNUM  Unit numbers to write message
      SUBROUTINE RATCHK(NUMU,UNITNUM)
      REAL RTS(99),RATIO,TSTAND
      !.. loop control
      INTEGER K
      !.. number of units to write message
      INTEGER NUMU
      !.. Unit numbers for writing warning
      INTEGER UNITNUM(22)
      !.. # of rates changed    
      INTEGER IRTS(99)
      !.. standard reaction rate values
      REAL R(99)
      DATA RTS/99*0.0/
      !.. rate values at the standard temperature
      DATA R/4.96E-10,6.40E-10,6.7E-13,8.56E-12,1.44E-07,8.41E-08
     > ,5.68E-13,1.02E-09,3.40E-11,7.83E-11,1.13E-07,3.30E-08,1.01E-07
     > ,1.66E-08,6.90E-13,8.06E-12,1.26E-11,2.08E-19,2.91E-10,3.88E-10
     > ,1.00E-10,0.0E-10,4.5E-10,1.99E-12,3.47E-10,4.0E-10,0.0,1.0E-11
     > ,1.30E-10,7.42E-11,2.20E-12,1.46E+00,2.00E-11,3.42E-11,2.50E-11
     > ,2.00E-11,1.70E-11,3.67E-12,2.20E-11,1.80E-10,6.70E-11,6.94E-02
     > ,7.00E-10,5.20E-10,7.80E-10,8.50E-09,4.83E-10,1.95E-12,2.16E-08
     > ,1.22E-07,9.34E-08,4.21E-09,4.60E-01,9.34E-03,1.06E+00,1.06E-01
     > ,7.90E-02,5.00E-03,6.00E-13,3.70E-01,1.07E-05,5.00E-01,4.00E-01
     > ,1.00E-01,1.49E-10,3.01E-10,1.00E-03,1.00E-01,7.66E-12,6.00E-02
     > ,3.00E-01,3.00E-02,8.00E-02,1.00E-01,9.20E-12,2.37E-10,1.80E-10
     > ,8.65E-10,1.00E-11,3.60E-10,3.54E-10,6.24E-11,1.20E-09,1.50E-10
     > ,1.30E-10,1.30E-10,1.00E-11,1.20E-09,3.80E-09,2.50E-11,2.39E-11
     > ,6.04E-10,4.60E-11,1.35E-09,1.00E-10,5.67E-10,9.50E-09,8.58E-10
     > ,5.31E-12/


      JRTS=0
      !.. call RATS with standard temperature TSTAND = 999
      TSTAND = 999.0
      CALL RATS(1,TSTAND,TSTAND,TSTAND,RTS)
      DO 20 I=1,99
      IF(R(I).GT.0) RATIO=RTS(I)/R(I)
      IF(R(I).LE.0.AND.RTS(I).EQ.0) RATIO=1.0
      IF(R(I).LE.0.AND.RTS(I).NE.0) RATIO=2.0
      IF(RATIO.GT.1.01.OR.RATIO.LT.0.99) THEN
        JRTS=JRTS+1
        IRTS(JRTS)=I
        ENDIF
 20   CONTINUE

      !.. Write warning in files
      IF(JRTS.GT.0) THEN
        DO K=1,NUMU
          WRITE(UNITNUM(K),88) (IRTS(I),I=1,JRTS)
        ENDDO
        DO K=1,NUMU
          WRITE(UNITNUM(K),89) (R(IRTS(I)),I=1,JRTS)
        ENDDO
        DO K=1,NUMU
          WRITE(UNITNUM(K),90) (RTS(IRTS(I)),I=1,JRTS)
        ENDDO
        DO K=1,NUMU
          !.. add a blank line
          WRITE(UNITNUM(K),*) ' '
        ENDDO
      ENDIF
 88   FORMAT(2X,'*** WARNING: Non-standard reaction rates'
     >  ,1X,'in FLIP, #s:',9I9)
 89   FORMAT(2X,'*** The FLIP standard reaction rates at'
     > ,1X,'Te=Ti=Tn=999 K:',1P,9E9.2)
 90   FORMAT(2X,'*** The  non-standard reaction rates at'
     > ,1X,'Te=Ti=Tn=999 K:',1P,9E9.2)
      RETURN
      END
C.................................... ALT_RATES ................. 
C.. This routine reads new rates off a file
C.. Written by P. Richards in May 2008
C.. Currently only capable of this functional form
C.. RTS(IRAT)=4.0E-7*(300/TE)**0.85 (COEFF*(TCON/TE)**POWER
C.. parameter TEMP = RCOEFF(I,4) converts between Te, Ti, Tn
C.. RATE  COEFF   TCON  TEMP POWER  others 
C..  3   1.2E-12  300.0  2.0  0.45  0 0 0 0 
      !.. electron temperature
      !.. ion temperature
      !.. neutral temperature
      !.. reaction rates
      SUBROUTINE ALT_RATES(TE,
     >                     TI,
     >                     TN,
     >                     RTS)
      IMPLICIT NONE
      !.. Loop control variable
      INTEGER I,NRATS
      !.. Rate number to process
      INTEGER IRAT
      !.. parameter to if file exists
      INTEGER FILEXIST
      !.. See input parameters
      REAL TE,TI,TN,RTS(99)
      !..
      REAL COEFFS(9), RCOEFF(99,9)
      FILEXIST=-99

      !.. Check to see if the rates file exists and store rate data for later calls
      IF(FILEXIST.EQ.-99) THEN
        CALL TFILE(3,FILEXIST)
        IF(FILEXIST.EQ.0) RETURN

        !.. Initialize rates numbers
        DO I=1,99
          RCOEFF(I,1)=0.0
        ENDDO
        NRATS=0
 10     READ(3,*,END=20,ERR=10) (COEFFS(I),I=1,9)

          NRATS=NRATS+1

          !.. WRITE(6,'(I6,1P,22E10.2)') NRATS,(COEFFS(I),I=1,9)
          !.. Rate number
          IRAT=NINT(COEFFS(1))

          DO I=1,9
            RCOEFF(NRATS,I)=COEFFS(I)
          ENDDO
        GO TO 10
      ENDIF

 20     CONTINUE


      !.. RATE    COEFF   TCON  TEMP POWER  others 
      !.. Currently only capable of this functional form
      !.. RTS(IRAT)=4.0E-7*(300/TE)**0.85
      !.. Now Reconstruct the rates
      DO I=1,NRATS
        IRAT=NINT(RCOEFF(I,1))  !.. Rate number
        IF(IRAT.GT.0.AND.RCOEFF(I,4).EQ.1) THEN
          RTS(IRAT)=RCOEFF(I,2)*(RCOEFF(I,3)/TE)**RCOEFF(I,5)
        ELSEIF(IRAT.GT.0.AND.RCOEFF(I,4).EQ.2) THEN
          RTS(IRAT)=RCOEFF(I,2)*(RCOEFF(I,3)/TI)**RCOEFF(I,5)
        ELSEIF(IRAT.GT.0.AND.RCOEFF(I,4).EQ.3) THEN
          RTS(NINT(RCOEFF(I,3)))=RCOEFF(I,2)*
     >      (RCOEFF(I,3)/TN)**RCOEFF(I,5)
        ELSE
          WRITE(6,*) ' ERROR: bad value in rates data file'
          STOP
        ENDIF
      ENDDO

      RETURN
      END
