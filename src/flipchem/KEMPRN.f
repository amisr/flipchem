C.................... KEMPRN.FOR ......................................   
C.. This file contains the chemistry routines for ions and neutrals
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CN2D(JPR,I,JPT,Z,RTS,ON,O2N,N2N,NOP,NE,P1,L1
     > ,N2PLS,DISN2D,UVDISN,NPLUS,N2P,N2D,OPLS,NNO,N2A)
C.......n(2d)
      IMPLICIT REAL(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=NOP*NE*RTS(50)
      PR(2)=N2PLS*NE*RTS(32)*RTS(11)
      PR(3)=N2PLS*ON*RTS(10)
      PR(4)=DISN2D
      PR(5)=RTS(63)*UVDISN
      PR(6)=RTS(65)*NPLUS*O2N
      PR(7)=N2P*RTS(57)
      PR(8)=RTS(27)*N2A*ON
      LR(1)=ON*RTS(15)
      LR(2)=O2N*RTS(16)
      LR(3)=NE*RTS(8)
      LR(4)=OPLS*RTS(29)
      LR(5)=RTS(61)
      LR(6)=RTS(41)*NNO
      P1=PR(1)+PR(2)+PR(3)+PR(4)+PR(5)+PR(6)+PR(7)+PR(8)
      L1=LR(1)+LR(2)+LR(3)+LR(4)+LR(5)+LR(6)
C....... EF is used to convert production rates to volume emission rates
      EF=1.0
C....... This line in for volume emission rates
C...      EF=RTS(61)*0.76/L1
      IF(JPT.EQ.1.AND.JPR.GT.0.AND.INT(EF+0.1).NE.1) WRITE(I,189)
      IF(JPT.EQ.1.AND.JPR.GT.0.AND.INT(EF+0.1).EQ.1) WRITE(I,191)
 189  FORMAT(/2X,'N(2D)',25X,'EMISSION',28X,':',20X,'Loss rate')
 191  FORMAT(/2X,'N(2D)',25X,'Production',36X,':',20X,'Loss rate')
      IF(JPT.EQ.1.AND.JPR.GT.0) WRITE(I,96)
 96   FORMAT(2X,'ALT   [N2D]   NO++e   N2++e   N2++O    e+N2   hv+N2'
     >  ,3X,'N++O2   N(2P)   N2A+O    +O     +O2      +e     +O+'
     >  ,5X,'RAD     +NO')
      IF(JPR.GT.0) WRITE(I,7) Z,N2D,(PR(K)*EF,K=1,8)
     > ,(LR(K)*N2D,K=1,6)
      RETURN
 7    FORMAT(F6.1,1P,22E8.1)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CNO(JPR,I,JPT,Z,RTS,ON,O2N,N2N,NE,P1,L1
     > ,N2D,N4S,N2P,NNO,O2P,OPLS,PDNOSR,PLYNOP,N2A,NPLUS)
C........no
      IMPLICIT REAL(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=RTS(16)*O2N*N2D
      PR(2)=RTS(7)*O2N*N4S
      PR(3)=RTS(38)*N2P*O2N
      PR(4)=RTS(27)*N2A*ON
      !.. Fox
      PR(5)=RTS(22)*NPLUS*O2N
      LR(1)=RTS(9)*N4S
      LR(2)=RTS(23)*O2P
      LR(3)=RTS(24)*OPLS
      LR(4)=RTS(41)*N2D
      LR(5)=PDNOSR
      LR(6)=PLYNOP
      P1=PR(1)+PR(2)+PR(3)+PR(4)+PR(5)
      L1=LR(1)+LR(2)+LR(3) + LR(4) + (LR(5)+LR(6))
      IF(JPT.EQ.1.AND.JPR.GT.0) WRITE(I,192)
 192   FORMAT(/2X,'NO',17X,'PRODUCTION',20X,':',10X,'LOSS RATES'/
     > ,4X,'ALT',3X,'[NO]',5X,'[NO]c',3X,'O2+N2D',
     > 3X,'O2+N4S   N2P+O2   N2A+O    N++O2    N4S+NO   O2P+NO   O++NO'
     > ,3X,'N2D+NO   hv<1910   Lyman-a')
      IF(JPR.GT.0) WRITE(I,7) Z,NNO,P1/L1,(PR(K),K=1,5)
     > ,(LR(K)*NNO,K=1,6)
      RETURN
 7    FORMAT(F6.1,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CN4S(JPR,I,JPT,Z,RTS,ON,O2N,N2N,NE,P1,L1,N4S,DISN4S
     >   ,N2D,N2P,OPLS,N2PLS,UVDISN,NOP,NPLUS,NNO,O2P,PDNOSR,VCON)
C........N(4S)
      IMPLICIT REAL(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=DISN4S
      PR(2)=RTS(15)*ON*N2D
      PR(3)=RTS(8)*NE*N2D
      PR(4)=VCON*RTS(3)*OPLS*N2N
      PR(5)=RTS(53)*RTS(11)*N2PLS*NE
      PR(6)=RTS(62)*UVDISN
      PR(7)=NOP*NE*RTS(49)
      PR(8)=N2D*RTS(61)
      PR(9)=N2P*RTS(58)
      PR(10)=RTS(25)*NPLUS*O2N
      PR(11)=PDNOSR*NNO
      !..Fox
      PR(12)=NPLUS*NNO*RTS(81)
      LR(1)=RTS(7)*O2N
      LR(2)=RTS(9)*NNO
      LR(3)=RTS(21)*O2P
      !..Fox
      LR(4)=RTS(79)*N2PLS
      P1=PR(1)+PR(2)+PR(3)+PR(4)+PR(5)+PR(6)+PR(7)+PR(8)+PR(9)+PR(10)
     >     +PR(11)+PR(12)
      L1=LR(1)+LR(2)+LR(3)+LR(4)
      IF(JPT.EQ.1.AND.JPR.GT.0) WRITE(I,193)
 193   FORMAT(/2X,'N(4S)',38X,'PRODUCTION',46X,':',7X,'LOSS RATES'/
     > ,3X,'ALT',2X,'[N4S]',2X,'hv->N+'
     > ,3X,'O+N2D',2X,'e+N2D',3X,'O++N2',3X,'N2++e',4X,'hv->2N'
     > ,2X,'NO++e',2X,'N(2D)',4X,'N(2P)   N+&X    hv+NO    +O2  '
     > ,2X,' +NO  ',2X,'+O2+ & N2+')
      !... for printing fit
      PR(10)=PR(10)+PR(12)
      !... for printing fit
      LR(3)=LR(3)+LR(4)
      IF(JPR.GT.0) WRITE(I,7) Z,N4S,(PR(K),K=1,11),(LR(K)*N4S,K=1,3)
      RETURN
 7    FORMAT(F6.1,1P,22E8.1)
      END
C::::::::::::::::::::::::::::::: CN2PLS :::::::::::::::::::::::::::::::
C..... Simplified chemistry of N2+.  PUN2P* = production of N2+ by euv 
C..... in the (X,A,B states). PEN2P* same for p.e.s (X,A,B states)
      SUBROUTINE CN2PLS(JPR,I,JPT,Z,RTS,ON,O2N,N2N,NE,N2PLS,PUN2PX,
     >  PUN2PA,PUN2PB,PEN2PX,PEN2PA,PEN2PB,OP2D,OP2P,HEPLUS,NPLUS,
     >  NNO,N4S)
      IMPLICIT REAL(A-H,L,N-Z)
      REAL PUN2PX,PUN2PA,PUN2PB,PEN2PX,PEN2PA,PEN2PB
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=PUN2PX
      PR(2)=PUN2PA
      PR(3)=PUN2PB
      PR(4)=PEN2PX
      PR(5)=PEN2PA
      PR(6)=PEN2PB
      PR(7)=RTS(19)*OP2D*N2N
      PR(8)=RTS(20)*OP2P*N2N
      PR(9)=RTS(44)*HEPLUS*N2N
      !..Fox
      PR(10)=RTS(82)*NPLUS*NNO
      LR(1)=RTS(10)*ON
      LR(2)=RTS(11)*NE
      LR(3)=RTS(17)*O2N
      LR(4)=RTS(99)*ON
      !..Fox
      LR(5)=RTS(79)*N4S
      !..Fox
      LR(6)=RTS(80)*NNO
      N2PLS=
     >  (PR(1)+PR(2)+PR(3)+PR(4)+PR(5)+PR(6)+PR(7)+PR(8)+PR(9)+PR(10))/
     >  (LR(1)+LR(2)+LR(3)+LR(4)+LR(5)+LR(6))
      IF(JPT.EQ.1.AND.JPR.GT.0) WRITE(I,95)
 95   FORMAT(/2X,'N2+',29X,'PRODUCTION',45X,':',12X,'LOSS RATES'/
     > ,3X,'ALT  [N2+]  EUV-X   EUV-A    EUV-B   PE-X'
     > ,5X,'PE-A    PE-B  O+2D+N2  O+2P+N2  He++N2  O+N2+'
     > ,2X,'e+N2+  O2+N2+  N2++O  Other')
      !.. for printing fit
      PR(9)=PR(9)+PR(10)
      !.. for printing fit
      LR(5)=LR(5)+LR(6)
      IF(JPR.GT.0) WRITE(I,7) Z,N2PLS,(PR(K),K=1,9)
     > ,(LR(K)*N2PLS,K=1,5)
      RETURN
 7    FORMAT(F6.1,1P,22E8.1)
      END
C:::::::::::::::::::::::::::::: CNOP ::::::::::::::::::::::::::::::::::
      SUBROUTINE CNOP(JPR,I,JPT,Z,RTS,ON,O2N,N2N,NE,P1,NOP,OPLS
     >  ,N2PLS,O2P,N4S,NNO,NPLUS,N2P,PLYNOP,VCON,N2D,OP2D)
C........no+
      IMPLICIT REAL(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=VCON*RTS(3)*N2N*OPLS
      PR(2)=N2PLS*ON*RTS(10)
      PR(3)=O2P*N4S*RTS(21)
      PR(4)=O2P*NNO*RTS(23)
      !.. N+ + O2 -> O2+ + N(2D,4S) or NO+ + O(1S)
      PR(5)=(RTS(30)+RTS(66)+RTS(59))*NPLUS*O2N
      PR(6)=RTS(37)*N2P*ON
      PR(7)=RTS(24)*OPLS*NNO
      PR(8)=PLYNOP*NNO
      !..Fox
      PR(9)=O2P*N2D*RTS(77)
      !..Fox
      PR(10)=N2PLS*NNO*RTS(80)
      !..Fox
      PR(11)=NPLUS*NNO*RTS(81)
      !..Fox
      PR(12)=RTS(83)*NNO*OP2D
      !.. -> NO+ + N, Li et al. [1997] 
      PR(13)=OP2D*RTS(90)*N2N
      LR(1)=NE*RTS(5)
      P1=PR(1)+PR(2)+PR(3)+PR(4)+PR(5)+PR(6)+PR(7)+PR(8)
     >    +PR(9)+PR(10)+PR(11)+PR(12)+PR(13)
      NOP=P1/LR(1)

      IF(JPT.EQ.1.AND.JPR.GT.0) WRITE(I,96)
 96   FORMAT(/2X,'NO+',31X,'PRODUCTION',48X,':',2X,'LOSS RATES'/
     > ,3X,'ALT',3X,'[NO+]',4X,'O++N2',3X,'N2++O',3X,'O2++N4S'
     > ,3X,'O2++NO',3X,'N++O2',4X,'N2P+O',3X,'O++NO   hv+NO'
     > ,5X,'O2++N2D   N2++NO   N++NO   OP2D+NO   OP2D+N2  NO++e')
      !PR(9)=PR(9)+PR(10)+PR(11)+PR(12)+PR(13)
      IF(JPR.GT.0) WRITE(I,7) Z,NOP,(PR(K),K=1,13),LR(1)*NOP
      RETURN
 7    FORMAT(F6.1,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CO2P(JPR,I,JPT,Z,RTS,ON,O2N,N2N,NE,P1
     > ,O2P,TPROD5,OPLS,OP2D,N2PLS,NPLUS,N4S,NNO,OP2P)
C........o2+
      IMPLICIT REAL(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
C......... TPROD5=euv @ p.e. production
      PR(1)=TPROD5
      PR(2)=RTS(4)*O2N*OPLS
      PR(3)=RTS(43)*OP2D*O2N
      PR(4)=RTS(17)*O2N*N2PLS
      PR(5)=RTS(25)*NPLUS*O2N
      !.. Fox
      PR(6)=RTS(86)*OP2P*O2N
      PR(7)=RTS(65)*NPLUS*O2N
      LR(1)=RTS(6)*NE
      LR(2)=RTS(21)*N4S
      LR(3)=RTS(23)*NNO
      P1=PR(1)+PR(2)+PR(3)+PR(4)+PR(5)+PR(6)+PR(7)
      O2P=P1/(LR(1)+LR(2)+LR(3))
      IF(JPT.EQ.1.AND.JPR.GT.0) WRITE(I,97)
  97  FORMAT(/2X,'O2+',22X,'PRODUCTION',24X,':',12X,'LOSS RATES'
     > /,3X,'ALT',3X,'[O2+]',3X,'hv+O2',3X,'O++O2',3X,'O+(2D)+O2'
     >  ,4X,'N2++O2   N++O2   O+(2P)+O2  O2++e   O2++N   O2++NO')
      IF(JPR.GT.0) WRITE(I,7) Z,O2P,(PR(K),K=1,7),(LR(K)*O2P,K=1,3)
      RETURN
 7    FORMAT(F6.1,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE COP4S(JPR,I,JPT,Z,RTS,ON,O2N,N2N,NE,OPLS,TPROD1,OP2D
     >  ,OP2P,PEPION,PDISOP,N2PLS,N2D,NNO,VCON,HEPLUS)
C...........o+(4s)
      IMPLICIT REAL(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
C.........pr(1)= euv production of o+(4s)
      PR(1)=TPROD1
      PR(2)=OP2D*NE*RTS(12)
      PR(3)=OP2P*ON*RTS(26)
      PR(4)=PEPION
      PR(5)=PDISOP
      PR(6)=RTS(99)*N2PLS*ON
      PR(7)=OP2P*NE*RTS(14)
      PR(8)=OP2P*0.047
      PR(9)=RTS(28)*ON*OP2D
      !.. Fox
      PR(10)=RTS(85)*OP2P*O2N
      !..Fox
      PR(11)=HEPLUS*O2N*(RTS(91)+RTS(93))
      !..Fox
      PR(12)=RTS(95)*NNO*HEPLUS
      !..Fox
      !PR(13)=RTS(22)*NPLUS*O2N
      LR(1)=N2N*VCON*RTS(3)
      LR(2)=O2N*RTS(4)
      LR(3)=NNO*RTS(24)
      !.... small loss?? ..Fox
      LR(4)=N2D*RTS(29)
      !.. total loss for printing
      !..LR(4)=(LR(1)+LR(2)+LR(3))
      PR(10)=PR(10)+PR(11)+PR(12)+PR(13)
      PRTOT=PR(1)+PR(2)+PR(3)+PR(4)+PR(5)+PR(6)+PR(7)+PR(8)+PR(9)+PR(10)
      LRTOT=LR(1)+LR(2)+LR(3)+LR(4)
      OPLS=PRTOT/LRTOT
      IF(JPT.EQ.1.AND.JPR.GT.0) WRITE(I,98)
 98   FORMAT(/2X,'O+',41X,'PRODUCTION',39X,':',10X,'LOSS RATES'/
     >  ,' ALT    [O+]   hv+O  O+(2D)+e O+(2P)+O   e*+O  O2-diss  '
     >  ,'N2++O  O+(2P)+e O+(2P) O+O+(2D)   Other  +N2     +O2    '
     >  ,'+NO   +N2D')
      IF(JPR.GT.0) WRITE(I,7) Z,OPLS,(PR(K),K=1,10)
     > ,(LR(K)*OPLS,K=1,4)
      RETURN
 7    FORMAT(F6.1,1P,22E8.1)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE COP2D(JPR,I,JPT,Z,RTS,ON,O2N,N2N,NE,OP2D,TPROD2,OP2P
     > ,HEPLUS,N4S,NNO,PSEC)
C.......o+(2d)
      IMPLICIT REAL (A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
      ! EUV  prod
      PR(1)=TPROD2
      PR(2)=OP2P*NE*RTS(13)
      PR(3)=OP2P*0.171
      !..Fox
      PR(4)=HEPLUS*O2N*RTS(76)
      PR(5)=PSEC
      LR(1)=RTS(19)*N2N
      !.. radiation at 3726 and 3729 A
      LR(2)=7.7E-5
      LR(3)=NE*RTS(12)
      LR(4)=ON*RTS(28)
      LR(5)=RTS(43)*O2N
      !..Fox
      LR(6)=RTS(83)*NNO
      !..Fox
      LR(7)=RTS(84)*N4S
      !.. -> NO+ + N, Li et al. [1997] 
      LR(8)=RTS(90)*N2N
      OP2D=(PR(1)+PR(2)+PR(3)+PR(4)+PR(5))/
     >   (LR(1)+LR(2)+LR(3)+LR(4)+LR(5)+LR(6)+LR(7)+LR(8))
      IF(JPT.EQ.1.AND.JPR.GT.0) WRITE(I,99)
 99   FORMAT(/2X,'O+(2D)',13X,'PRODUCTION',27X,':',18X,'LOSS RATES'/
     > ,3X,'ALT',3X,'[O+2D]',3X,'hv+O',4X,'e*+O',4X,'O+2P+e',3X,
     > 'O+2P>hv',2X,'He++O2     +N2    E3726_29    +e       +O      +O2'
     > '      +NO     +N  +N2>NO+')
      IF(JPR.GT.0) WRITE(I,7) Z,OP2D,PR(1),PR(5),PR(2),PR(3),PR(4)
     > ,(LR(K)*OP2D,K=1,8)
      RETURN
 7    FORMAT(F6.1,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE COP2P(JPR,I,JPT,Z,RTS,ON,O2N,N2N,NE,OP2P,TPROD3,PSEC
     > ,HEPLUS,N4S,NNO,TE)
C.......o+(2p)
      IMPLICIT REAL(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=0.0
      IF(TPROD3.GE.PSEC) PR(1)=TPROD3-PSEC
      PR(2)=PSEC
      !..Fox
      PR(3)=HEPLUS*O2N*RTS(92)
      LR(1)=RTS(26)*ON
      LR(2)=RTS(20)*N2N
      LR(3)=RTS(13)*NE
      LR(4)=0.218
      LR(5)=RTS(14)*NE
      !..Fox
      LR(6)=(RTS(85)+RTS(86))*O2N
      !..Fox
      LR(7)=RTS(87)*N4S
      !..Fox
      LR(8)=RTS(88)*NNO
      OP2P=(TPROD3+PR(3))
     >  /(LR(1)+LR(2)+LR(3)+LR(4)+LR(5)+LR(6)+LR(7)+LR(8))
      IF(JPT.EQ.1.AND.JPR.GT.0) WRITE(I,100)
 100   FORMAT(/2X,' O+(2P)',6X,'PRODUCTION',10X,':',12X,'LOSS RATES'/
     >,3X,'ALT   [O+2P]    hv+O     e*+O  He++O2      +O',7X,'+N2'
     > ,6x,'+e       RAD      +e      +O2      +N4S     +NO'
     > ,6x,'OX       N2        e      Te       E7320')
      IF(JPR.GT.0) WRITE(I,7) Z,OP2P,PR(1),PR(2),PR(3),
     >  (LR(K)*OP2P,K=1,8),ON,N2N,NE,TE,OP2P*0.218*0.781
      RETURN
 7    FORMAT(F6.1,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CNPLS(JPR,I,JPT,Z,RTS,ON,O2N,N2N,NE,DISNP
     > ,NPLUS,OPLS,N2D,OP2P,HEPLUS,PHOTN,O2P,N4S,OP2D,N2PLS,NNO)
C........n+
      IMPLICIT REAL(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=DISNP
      PR(2)=RTS(29)*OPLS*N2D
      PR(3)=0
      PR(4)=RTS(45)*HEPLUS*N2N
      PR(5)=PHOTN
      !..Fox
      PR(6)=O2P*N2D*RTS(78)
      !..Fox
      PR(7)=N2PLS*N4S*RTS(79)
      !..Fox
      PR(8)=OP2D*N4S*RTS(84)
      !..Fox
      PR(9)=RTS(94)*NNO*HEPLUS
      !..Fox
      LR(1)=RTS(30)*O2N
      !..Fox
      LR(2)=RTS(25)*O2N
      !..Fox
      LR(3)=RTS(22)*O2N
      !..Fox
      LR(4)=RTS(65)*O2N
      !..Fox
      LR(5)=RTS(66)*O2N
      !..Fox
      LR(6)=RTS(31)*ON

      CNPLUS=0.0
      IF(LR(1)+LR(2)+LR(3).GT.0.0)
     >CNPLUS=(PR(1)+PR(2)+PR(3)+PR(4)+PR(5)+PR(6)+PR(7)+PR(8)+PR(9))
     >   /(LR(1)+LR(2)+LR(3)+LR(4)+LR(5)+LR(6))
      NPLUS=CNPLUS

      IF(JPT.EQ.1.AND.JPR.GT.0) WRITE(I,101)
 101   FORMAT(/2X,'N+',20X,'PRODUCTION',71X,':',8X,'LOSS RATES'/
     > ,4X,'ALT   [N+]   [N+]c     hv+N2   O++N2D  O+2P+N2',3X
     > ,'He++N2',3X,' hv+N   O2++N2D  N2++N4S O+(2D)+N4S  He++NO'
     > ,3X,'N++O2    N++O2    N++O2    N++O2    N++O2    N++O')
      IF(JPR.GT.0) WRITE(I,7) Z,NPLUS,CNPLUS
     > ,(PR(K),K=1,9),(LR(K)*NPLUS,K=1,6)
      RETURN
 7    FORMAT(F6.1,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CN2A(JPR,I,JPT,Z,RTS,ON,O2N,N2N,NE
     > ,N2A,P3X1,P3X2,P3X3,P3X4)
C........n2(a3sigma) and total LBH
      IMPLICIT REAL(A-H,L,N-Z)
      REAL P3X1,P3X2,P3X3,P3X4
      DIMENSION RTS(99),LR(22),PR(22)
C....... pr(1,2,3)= electron impact excitation of n2(a,b,c) states
      PR(1)=P3X1
      PR(2)=P3X2
      PR(3)=P3X3
      LR(1)=RTS(36)*ON
      LR(2)=RTS(27)*ON
      LR(3)=0.57
      N2A=(PR(1)+PR(2)+PR(3))/(LR(1)+LR(2)+LR(3))
      IF(JPT.EQ.1.AND.JPR.GT.0) WRITE(I,102)
 102   FORMAT(/2X,'N2(A)',12X,'PRODUCTION',13X,':',5X,'LOSS RATES'
     > ,3X,':  Total LBH'
     > /,3X,'ALT',3X,'N2(A)',3X,'e*->N2A',3X,'e*->N2B',3X,'e*->N2C',2X
     >  ,'N2A>O1S',2X,'N2A>NO',2X,'RAD',5X,'LBH')
      IF(JPR.GT.0) WRITE(I,7) Z,N2A,(PR(K),K=1,3),(LR(K)*N2A,K=1,3)
     > ,P3X4
      RETURN
 7    FORMAT(F6.1,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CN2P(JPR,I,JPT,Z,RTS,ON,O2N,N2N,NE,P1,L1
     > ,N2P,P3X7,UVDISN,O2P,NNO,N2PLUS)
C....... N(2P). the rates are from Zipf et al JGR, 1980 p687
C.... 21-AUG-1992. Added N2+ recombination source
      IMPLICIT REAL(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=P3X7
      PR(2)=RTS(64)*UVDISN
      PR(3)=RTS(73)*RTS(11)*N2PLUS*NE
      LR(1)=RTS(37)*ON
      LR(2)=RTS(38)*O2N
      LR(3)=RTS(39)*O2P
      LR(4)=RTS(40)*NNO
      LR(5)=RTS(57)
      LR(6)=RTS(58)
      LR(7)=(RTS(96)+RTS(97))*NE    !..Fox
      P1=PR(1)+PR(2)+PR(3)
      L1=LR(1)+LR(2)+LR(3)+LR(4)+LR(5)+LR(6)+LR(7)
      IF(JPT.EQ.1.AND.JPR.GT.0) WRITE(I,103)
 103   FORMAT(/2X,'N(2P)',9X,'PRODUCTION',17X,':',20X,'LOSS RATES'/
     > ,3X,'ALT',3X,'[N2P]',3X,'e+N2',5X,'hv+N2',3X,'e+N2+',6X,'+O  '
     > ,3X,'+O2      +O2+      +NO       +2D     +4S      +e')
      IF(JPR.GT.0) WRITE(I,7) Z,N2P,(PR(K),K=1,3),(LR(K)*N2P,K=1,7)
      RETURN
 7    FORMAT(F6.1,1P,22E9.2)
      END
C:::::::::::::::::::::::::::::: CNOPV ::::::::::::::::::::::::::::::::::
C...... This routine calculates the vibrational distribution of no+
C...... Uses AFRL report by Winick et al. AFGL-TR-87-0334, Environmental
C...... Research papers, NO. 991, "An infrared spectral radiance code 
C...... for the auroral thermosphere (AARC)
C...... Written by P. Richards in February 2004
      SUBROUTINE CNOPV(JPR,I,JPT,Z,RTS,ON,O2N,N2N,NE,P1,NOP,OPLS
     >  ,N2PLS,O2P,N4S,NNO,NPLUS,N2P,PLYNOP,VCON,N2D,OP2D)
      IMPLICIT NONE
      INTEGER JPR,I,JPT,INV,IJ,IV,K
      !.. NO+(v) array dimensions
      PARAMETER (INV=20)
      !.. Altitude, rate coefficients
      !.. O, N2, O2, electron densities
      !.. total source output for finding [e] 
      !.. NO+, )+,N2+,O2+ densities
      !.. N(4S), NO, N+, N(2P) densities
      !.. Lyman-a source, N2(v) rate factor
      !.. N(2D), O+(2D) densities                
      !.. NO+(v) densities and total NO+ density
      !.. storage for NO+ sources and sinks
      !.. Einstein coeffs for delv=1,2
      !.. NO+(v) + e rate factors
      !.. Sources and sinks of NO+(v)
      !.. Temp total cascade source, sink
      !.. N2 queching rate :- coeff, source, sink
      REAL Z,RTS(99),
     >  ON,O2N,N2N,NE,
     >  P1,
     >  NOP,OPLS,N2PLS,O2P,
     >  N4S,NNO,NPLUS,N2P,
     >  PLYNOP,VCON,
     >  N2D,OP2D,
     >  NOPV(INV),NOPTOT,
     >  LR(22),PR(22),
     >  EINSCO1(INV),EINSCO2(INV),
     >  LRV(INV),
     >  PNOPV,LNOPV,
     >  PCASC,LRAD,
     >  K_N2_Q,P_N2_Q,L_N2_Q

      !.. Fractions of each source going to each vib. level. Assume
      !.. N+ + O2 fractions for each source. NEED UPDATE
      REAL PRV1(INV),PRV2(INV),PRV3(INV),PRV4(INV),PRV5(INV),PRV6(INV),
     >  PRV7(INV),PRV8(INV),PRV9(INV),PRV10(INV),PRV11(INV),PRV12(INV)
      DATA PRV1/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV2/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV3/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV4/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV5/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV6/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV7/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV8/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV9/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV10/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV12/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/

      DATA K_N2_Q/7.0E-12/  !.. Quenching rate coeff. by N2
      !.. Einstein coeff delv=1
      DATA EINSCO1/0.0,10.9,20.2,28.4,35.5,41.5,46.6,50.9,54.3,57.0,
     >     58.9,60.2,60.8,60.7,59.9,5*60.0/
      !.. Einstein coeff delv=1
      DATA EINSCO2/0.0,0.0,.697,1.93,3.61,5.74,8.24,11.1,14.2,17.7,
     >     21.3,25.1,29.0,33.2,37.4,5*40.0/
      !.. rate factors for NO+(v)+e -> N + O. Sheehan and St-Maurice 2004
      DATA LRV/1.0,19*0.3333/    
      
      !... Evaluate total production and loss rates
      !.. N2 + O+
      PR(1)=VCON*RTS(3)*N2N*OPLS
      !.. N2+ + O
      PR(2)=N2PLS*ON*RTS(10)
      !.. O2+ + N(4S)
      PR(3)=O2P*N4S*RTS(21)
      !.. O2+ + NO
      PR(4)=O2P*NNO*RTS(23)
      !.. N+ + O2 -> O2+ + N(2D,4S) or NO+ + O(1S)
      PR(5)=(RTS(30)+RTS(66)+RTS(59))*NPLUS*O2N
      !.. N2+ + O
      PR(6)=RTS(37)*N2P*ON
      !.. O+ + NO
      PR(7)=RTS(24)*OPLS*NNO
      !.. Lyman-a + NO
      PR(8)=PLYNOP*NNO
      !.. Fox: O2+ + N(2D)
      PR(9)=O2P*N2D*RTS(77)
      !.. Fox: N2+ + NO
      PR(10)=N2PLS*NNO*RTS(80)
      !.. Fox: N+ + NO
      PR(11)=NPLUS*NNO*RTS(81)
      !.. Fox: O+(2D) + NO
      PR(12)=RTS(83)*NNO*OP2D
      !.. -> NO+ + N, Li et al. [1997] 
      PR(13)=OP2D*RTS(90)*N2N
      !.. NO+ + e
      LR(1)=NE*RTS(5)

      !..Total source term used in main program to calculate [e]
      P1=PR(1)+PR(2)+PR(3)+PR(4)+PR(5)+PR(6)+PR(7)+PR(8)
     >    +PR(9)+PR(10)+PR(11)+PR(12)+PR(13)
      !.. NO+ density
      NOP=P1/LR(1)

      DO IJ=1,INV
        NOPV(IJ)=0.0
      ENDDO
      NOPTOT=0.0

      !.. loop down evaluating the vibrational populations. Must start 
      !.. less than INV-1 because need cascade from v+2
      DO IV=INV-4,1,-1
        !... chemical production for v=IV = total source * fraction 
        PNOPV=PR(1)*PRV1(IV)+PR(2)*PRV2(IV)+PR(3)*PRV3(IV)+
     >    PR(4)*PRV4(IV)+PR(5)*PRV5(IV)+PR(6)*PRV6(IV)+
     >    PR(7)*PRV7(IV)+PR(8)*PRV8(IV)+PR(9)*PRV9(IV)+
     >    PR(10)*PRV10(IV)+PR(11)*PRV11(IV)+PR(12)*PRV12(IV)

        !.. cascade production from v+1, v+2
        PCASC=NOPV(IV+1)*EINSCO1(IV+1)+NOPV(IV+2)*EINSCO2(IV+2)
        !.. total radiative loss
        LRAD=EINSCO1(IV)+EINSCO2(IV)

        !.. sink of quanta by N2 quenching
        L_N2_Q=K_N2_Q*N2N
        IF(IV.EQ.1) L_N2_Q=0.0

        !.. source of quanta by N2 quenching
        P_N2_Q=K_N2_Q*N2N*NOPV(IV+1)
        !.. recombination rate for level IV
        LNOPV=LR(1)*LRV(IV)

        !.. evaluate NO+(iV) population
        NOPV(IV)=(PNOPV+PCASC+P_N2_Q)/(LNOPV+LRAD+L_N2_Q)
        !.. total NO+ concentration
        NOPTOT=NOPTOT+NOPV(IV)

          !... diagnostic print. Set alt range to invoke
          IF(JPR.GT.0.AND.Z.GE.0.AND.Z.LT.10)
     >     WRITE(6,'(F10.1,I7,1P,22E10.2)') 
     >     Z,IV,PNOPV,P1,PCASC,P_N2_Q,LNOPV,LRAD,L_N2_Q,
     >     NOPV(IV),NOP,NOPTOT
      ENDDO

      !.. for printing only
      PR(9)=PR(9)+PR(10)+PR(11)+PR(12)
      IF(JPT.EQ.1.AND.JPR.GT.0) WRITE(I,96)
 96   FORMAT(/5X,'NO+',34X,'PRODUCTION',69X,':LOSS RATE'/
     > ,3X,'ALT NO+(v=0) NO+(v=1) NO+(v=2)  NO+(v=3) O++N2    N2++O',
     >  4X,'O2++N4S   O2++NO   N++O2    N2P+O   O++NO   hv+NO',
     >  4X,'Other_P   NO++e')
      IF(JPR.GT.0) WRITE(I,'(F6.1,1P,22E9.2)')  Z,
     >   (NOPV(K),K=1,4),(PR(K),K=1,9),LR(1)*NOP

      RETURN
      END
C:::::::::::::::::::::::::::::::: TFILE :::::::::::::::::::::::::
C------- This file checks to see if a file exists for the FORTRAN unit
C------- JUNIT and returns 0 for no and 1 for yes
       SUBROUTINE TFILE(JUNIT,IOPEN)
       CHARACTER ACHAR
       IOPEN=1
       READ(JUNIT,90,END=70,ERR=70) ACHAR 
 90    FORMAT(A)
       REWIND JUNIT
       RETURN

 70    CONTINUE
       !--- File did not exist
       IOPEN=0 
       RETURN
       END