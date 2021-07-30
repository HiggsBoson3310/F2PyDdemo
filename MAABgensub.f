      subroutine coulfg(ll,eps,rho,acc,f,fp,g,gp,k,ierr,actacc)
c
c  calculates coulomb functions f and g and their derivatives
c
c  input -
c        ll=angular momentum quantum number
c        eps=z-scaled energy in rydbergs
c        rho=z-scaled radial variable in atomic units
c        acc=accuracy required
c
c  output -
c        f=regular function
c        fp=derivative of f
c        g=irregular function
c        gp=derivative of g
c        k=number of terms needed in expansion
c        ierr=error code
c        actacc=accuracy actually achieved
c
c  convergence criterion -
c        value of wronskian converged to accuracy of 0.5*acc
c
c  error codes -
c        ierr=0, converged with actacc.lt.acc
c        ierr=1, converged with actacc.gt.acc
c        ierr=2, not converged with 101 terms in main summation
c
c  initialization
c
c
        implicit real*8 (a-h,o-z)
cd    delete previous card for single precision
c
ctj      data r2pi,ps0/.159154943d0,-.154431330d0/        !tj added d0 here. improves wronski-2/pi from 10^(-8) to 10^(-10)

        data r2pi,ps0/.159154943091895336d0,-.1544313298030657212d0/        !tj used r2pi value from fogo, and digamma function value from mathematica. improves wronski-2/pi from 10^(-8) to 10^(-16) to 10^(-13) depending on r


        ierr=0
        lp1=ll+1
        l2=2*ll
        l2p1=l2+1
        fl=ll
        flp1=lp1
        fl2p1=l2p1
        e2=0.5d0*eps
        r2=2.d0*rho
        acc2=2.d0*acc
c
c     initialize fa=factorial(2*ll+1)
c     and ps=psi(2*ll+2)+psi(1)
c
        fa=1.d0
        ps=ps0
c
c
c  calculate alpha(n) and beta(n) and initialize s and sp
c  continue calculation of fa and ps
c
c     s and sp for n=0
        x3=-l2
        x2=l2p1
        x1=-2.d0*r2**(-lp1)
        sp=x3*x1
        x1=r2*x1
        s=x1
c
c     initialize for coefficients in recursion formulae
        p1=fl*e2
        p2=p1
        q1=-e2
c
c     initialize alpha and beta
        alp1=1.d0
        alp2=1.d0+p2
        bet1=0.d0
        bet2=q1
c
        if(ll.eq.0)goto 20
c
c     s and sp for n=1
        x3=x3+2.d0
        x2=x2-1.d0
        x1=x1/x2
        sp=sp+x3*x1
        x1=r2*x1
        s=s+x1
c
c     loop for n=2 to 2*ll
        do 10 n=2,l2
c
c     continue calculation of fa and psi
        fn=n
        fa=fn*fa
        ps=ps+1.d0/fn
c
c     continue calculation of s and sp
        x3=x3+2.d0
        x2=x2-1.d0
        x1=x1/(x2*fn)
        sp=sp+x3*x1*alp2
        x1=r2*x1
        s=s+x1*alp2
c
c     compute coefficients in recursion formulae
        p1=p1-e2
        p2=p2+p1
        q1=q1-e2
c     now have p2=-n*(n-2*ll-1)*eps/4
c     and q1=-n*eps/2
c
c     new alpha and beta
        alp0=alp1
        alp1=alp2
        alp2=alp1+p2*alp0
        bet0=bet1
        bet1=bet2
 10     bet2=bet1+p2*bet0+q1*alp0
c
c     normalize s and sp, complete calculation of fa and ps
        s=s*fa
        sp=sp*fa
        fa=fl2p1*fa
        ps=ps+1.d0/fl2p1
c
c     complete calculation of alpha and beta
        p1=p1-e2
        p2=p2+p1
        q1=q1-e2
        alp0=alp1
        alp1=alp2
        bet0=bet1
        bet1=bet2
        bet2=bet1+p2*bet0+q1*alp0
c
 20     continue
c     now have alp1=alpha(2*ll+1)
c     and bet1=beta(2*ll+1), bet2=beta(2*ll+2)
c
c     value of a=a(eps,ll)
        a=alp1
        a4=4.d0*a
        cl=2.d0*a*dlog(dabs(r2))
cd    for single precision replace dlog by alog and dabs by abs
        clp=2.d0*a/rho
c
c  calculate a(n) and d(n), f and fp and
c  complete calculation of s and sp
c
c     calculate a0,a1,d0,d1
        a0=(2.d0**lp1)/fa
        a1=-a0/flp1
        ps=2.d0*ps*a
        d0=(bet1-ps)*a0
        d1=(bet2-ps-(2.d0+1.d0/flp1)*a)*a1
c
c     initialize f,fp, continue calculation of s,sp
c     - values for n=0
        fnplp1=flp1
        c1=rho**ll
        c1p=fnplp1*c1
        fp=c1p*a0
        sp=sp+c1p*d0
        c1=c1*rho
        f=c1*a0
        s=s+c1*d0
        w1=f*(clp*f+sp)-fp*s
c
c     - values for n=1
        fnplp1=fnplp1+1.d0
        c1p=fnplp1*c1
        fp=fp+c1p*a1
        sp=sp+c1p*d1
        c1=c1*rho
        f=f+c1*a1
        s=s+c1*d1
        w2=f*(clp*f+sp)-fp*s
        dw2=dabs(w2-w1)
cd    for single precision replace dabs by abs
c
c     initialize for coefficients in recursion formulae
        p1=-2.d0*flp1
        p2=p1
        q1=a4+2.d0*a*fl2p1
c
c     loop for n=2 to 100
        do 40 n=2,600
c
c     compute coefficients in recursion formulae
        p1=p1-2.d0
        p2=p2+p1
        q1=q1+a4
c     now have p2=-n*(n+2*ll+1)
c     and q1=2*a*(2*n+2*ll+1)
c
c     compute a2=a(n) and d2=d(n)
        a2=(2.d0*a1+eps*a0)/p2
        d2=(2.d0*d1+eps*d0+q1*a2)/p2
c
c     increment fp and sp
        fnplp1=fnplp1+1.d0
        fp=fp+a2*fnplp1*rho**(ll+2)
        sp=sp+d2*fnplp1*rho**(ll+2)
c
c     increment f and s
        f=f+a2*rho**(ll+3)
        s=s+d2*rho**(ll+3)
c
c     calculate wronskian
        w1=w2
        dw1=dw2
        w2=f*(clp*f+sp)-fp*s
        dw2=dabs(w2-w1)
cd    for single precision replace dabs by abs
c
c     convergence test
        k=n+1
        if(dw1.gt.acc2)goto 30
        if(dw2.gt.acc2)goto 30
        goto 50
c
c     new a0,a1,do,d1
30      a0=a1*rho
        a1=a2*rho
        d0=d1*rho
        d1=d2*rho
c
40      continue
c
c  not converged
c
        ierr=2
        actacc=dabs(0.25d0*w2-1.d0)
cd    for single precision replace dabs by abs
        goto 60
c
c  converged
c
50      actacc=dabs(0.25d0*w2-1.d0)
cd    for single precision replace dabs by abs
        if(actacc.gt.acc)ierr=1
c
c  complete calculation of g and gp
c
60      g=(s+cl*f)*r2pi
        gp=(sp+cl*fp+clp*f)*r2pi
c
        return
      end
c      
      subroutine seaton(l,eryd,r,zion,f,fp,g,gp)
	implicit real*8(a-h,o-z)
c -- assume for now that zion=1
	rho=zion*r
	epsr=eryd/(zion**2)
	acc=1.d-11
	call fogo(f,fp,g,gp,epsr,l,rho,ww,acc,actacc)
	s2=dsqrt(2.d0)
	f=f*s2
	fp=fp*s2
	g=g*s2
	gp=gp*s2
	return
	end
cccc
      SUBROUTINE FOGO(F0,FP0,G0,GP0,EPSR,LL,RO,WW,ACC,ACTACC)
*
*   INPUT EPSR (E RYDBERG) LL (L) ACC=1.E-11
*    OUTPUT ENERGY NORMALISED FUNCTIONS AND DERIVATIVES (RYDBERG)
*
*
        IMPLICIT REAL*8 (A-H,O-Z)
chg
	common/agqdt/aqdt,gqdt,fr0,fpr0,gr0,gpr0
chg
      DATA R2PI/.159154943 09189 5336D0/
       DATA PI/ 3.1415 92653 58979 32D0/
C     PRINT 1013
1013  FORMAT(//////,121(1H*) ,/)
*
       A=1.D0
! djh      IF(LL.EQ.0) GO TO 50
!       TESEN=1.D-10
       TESEN=1.D-10
       IF(LL.EQ.0) GO TO 50
       IF(DABS(EPSR).LT.TESEN)  GO TO 50
C
C
C CALCUL DE A= PRODUIT( 1. + EPS * P2 )
       DO 20 LP=1,LL
       B=1.D0 + LP*LP*EPSR
       A=A*B
C     PRINT *,' A B ',A,B
20     CONTINUE
50     CONTINUE
C     PRINT 1010 ,EPSR,LL,RO,A
1010  FORMAT(/////,' ENERGIE  ',D25.16,' RYD  L= ',I3,' RHO = ',
     =D25.16,/,' A = ',D25.16 ,//)
C
C
C
C CALCUL DE G ET DE B
C *******************
      IF(DABS(EPSR).LT.TESEN )  GO TO 8000
      IF(EPSR.GE.TESEN) GO TO 5000
C
C CAS ENERGIE NEGATIVE
      PNU=DSQRT(-1.D0/EPSR)
      PL1=DFLOAT(LL)
chg
	pl1=0.d0
chg
      IF(PNU.GT.PL1)  GO TO 61
      PRINT 1062,LL,EPSR
1062  FORMAT(//////,131(1H*),/,
     1'****** CAS ETUDIE IMPOSSIBLE  .  L = ',
     2I4,' ENERGIE =  ',D25.13,' RYDBERG ',/, 131(1H*) )
      STOP
61    CONTINUE
      X=PNU+LL+1.D0
      PSI=D1IGAM(X)
C     PRINT *,' X PSI ',X,PSI
      X=PNU-LL
      PSI=PSI +D1IGAM(X) -2.D0*DLOG(PNU)
      X2=DLOG(PNU)
      X1=D1IGAM(X)
C     PRINT *,' X PSI LOG ',X,X1,X2
      G=A*R2PI*PSI
C     PRINT 1020,G,PNU,PSI
1020  FORMAT('  G = ',D25.16,' KAPA ',D25.16,' PSI',D25.16)
      B=A
C     PRINT *, ' B = ',B
chg
	aqdt=b
	gqdt=g
chg
      GO TO 4000
5000    CONTINUE
C
C CAS ENERGIE POSITIVE
      PGA=DSQRT(1.D0/EPSR)
      SUM=0.D0
C CALCUL DE LA DEPENDANCE EN L
      IF(LL.EQ.0) GO TO 90
      DO 80 LP=1,LL
      S1=LP*EPSR
      SUM=SUM + S1/(1.D0 + LP*S1)
C     PRINT *,'  SUM S1  ',SUM,S1
80    CONTINUE
90    CONTINUE
      IF(EPSR . GT. 0.05D0 ) GO TO 92
C ENERGIE < 0.05 RYD . FORMULE 2
      X1=GAMI(PGA,2)
      GO TO 99
92    IF(EPSR . GT. 0.8D0 ) GO TO 96
C 0.05 RYD < ENERGIE < 0.8 RYD . FORMULE 1
      X1=GAMI(PGA,1)
      GO TO 99
C O.8 RYD < ENERGIE . FORMULE 3
96    X1=GAMI(PGA,3)
99    CONTINUE
      G=(SUM + X1)* A * 2.D0 * R2PI
C     PRINT 1040,G,PGA,SUM,X1
1040  FORMAT('  G = ',D25.16,' GAMA = ',D25.16,' SUM = ',
     +D25.16,' X = ',D25.16)
      B=A/(1.D0 - DEXP( -PGA/R2PI))
C     PRINT *, '  B = ',B
chg
	aqdt=b
	gqdt=g
chg
      GO TO 4000
8000   CONTINUE
C
C CAS ENERGIE=0
      G=0.D0
      B=1.D0
C     PRINT 1050,G
1050  FORMAT('  G = ',D25.16)
C     PRINT *,'  B = ',B
4000  CONTINUE
C
C CALCUL DE F(R0) ET DE G(R0) ANALYTIQUES EN ENERGIE
C **************************************************
C*********WRONSKIEN = 2 / PI
C         K TERMES DANS LE CALCUL
C         ACTACC PRECISION RELATIVE DU WRONSKIEN
C
      CALL COULFG(LL,EPSR,RO,ACC,FR0,FPR0,GR0,GPR0,KK
     1,IERR,ACTACC)
C     PRINT *
C     PRINT *
C     PRINT *,'  LL  EPSR   RO ',LL,EPSR,RO
C     PRINT *,' VALEURS DE F  FP  G  GP '
C     PRINT *,FR0,FPR0,GR0,GPR0
C     PRINT *,' K  IER WRONSKIEN ',KK,IERR,ACTACC
C
C CALCUL DE F0 ET G0 NORMEES EN RYDBERG
C *************************************
C
chg
	b=dabs(b)
chg
      HR0 =GR0 + G * FR0
      HPR0 = GPR0 + G * FPR0
      B = DSQRT ( B / 2.D0  )
      F0 = FR0 * B
      FP0 = FPR0 * B
      G0 = HR0 / ( 2.D0 * B )
       GP0 = HPR0 /( 2.D0 * B )
      WW=(F0*GP0-FP0*G0)*PI-1.D0
C       PRINT*,' FONCTION NORMEES EN RYDBERG'
C     PRINT *,F0,FP0,G0,GP0
C     PRINT*,' 2EME WRONSKIEN ',WW
C      PRINT *
      RETURN
      END
            subroutine seaton1(l,eryd,r,zion,f,fp,g,gp)
      implicit real*8(a-h,o-z)
      acc=1.d-12      ! tj changed from -10 to -16
      rl=l
      rr=r*zion
      eps=eryd/(zion**2)
      call coulfg(l,eps,rr,acc,f0,f0p,g0,g0p,k,ierr,actacc)

      if(.not.(eryd.lt.0))goto 23000
         ea=dabs(eryd)
         call ganda(a,gg,l,ea,zion,999)
         goto 23001
c     else
23000    continue
         gam=1.d0/dsqrt(eps)
         call gcmplx(a,gg,rl,gam)
23001 continue
      a5=dsqrt(dabs(a))
      f=a5*f0
      fp=a5*f0p
      g=(g0+gg*f0)/a5
      gp=(g0p+gg*f0p)/a5
c
c ** the next five lines changed on 1-22-88 by c.greene thanks to h. gao
c
	factor = dsqrt(zion)
	f=f/factor
	g=g/factor
	fp=fp*factor
	gp=gp*factor
c
      return
      end
c
c
      SUBROUTINE GANDA(A,G,L,E,ZION,NOPT)                               00000020
      IMPLICIT REAL*8(A-H,O-Z)                                          00000030
      DATA PI/3.1415926535897932D0/                                     00000040
      DPI = PI*2.D0                                                     00000050
C*** THIS PROGRAM RETURNS THE QDT TRANSFORMATION PARAMETERS A&G.        00000060
      IF(ZION.EQ.0) WRITE(6,90)                                         00000070
90    FORMAT(1X,'***** ZION = 0')                                       00000080
      IF(ZION.EQ.0) RETURN                                              00000090
      E = DABS(E)                                                       00000100
      XNUE = ZION/DSQRT(E)                                              00000110
C*** EVALUATE A(K,L) FIRST.                                             00000120
      A = 1.D0                                                          00000130
      IF(L.EQ.0) GO TO 109                                              00000140
      DO 100 I = 1,L                                                    00000150
      A = A* (1.D0 -I*I*E/(ZION**2) )                                   00000160
100   CONTINUE                                                          00000170
109   continue
C*** GIVE WARNINGS IN CASE A < OR = 0 .                                 00000180
C      IF(A.LE.0) WRITE(6,"(1X,'****** A < OR = 0', 2x, f10.3)") A                                         00000190
c5555  FORMAT(1X,'****** A < OR = 0', 2x, 'f10.3')                                    00000200
      IF(NOPT.EQ.1) RETURN                                              00000210
C*** CHECK WHETHER XNUE = INTEGER.                                      00000220
      N = XNUE + 1.D-02                                                 00000230
      Z = XNUE - N                                                      00000240
      IF(Z.EQ.0) G = 0.D0                                               00000250
      IF(Z.EQ.0) RETURN                                                 00000260
C*** G(K,L) IS NOW EVALUATED USING THE DIGAMMA PROGRAM.                 00000270
      G = A*(DIGAM(L+1+XNUE) + DIGAM(-L+XNUE)                           00000280
     1  - 2.D0 * DLOG(XNUE) )/DPI                                       00000290
      RETURN                                                            00000300
      END                                                               00000310
      SUBROUTINE GCMPLX(B,G,RL,GAM)                                     00000040
      IMPLICIT REAL*8(A-H,O-Y),COMPLEX*16(Z)                            00000050
CCCC  COMPLEX*16 ZQ1,ZQ2,ZQ3,ZSUM,ZG                                    00000060
      COMPLEX*16 CDLOG
      DIMENSION GG(2)                                                   00000070
      EQUIVALENCE (GG(1),ZG)                                            00000080
      DATA PI/3.1415926535897932D0/                                     00000090
C                                                                       00000100
      ZQ1 = DCMPLX(RL+1.D0,GAM)                                         00000110
      ZQ2 = DCMPLX(-RL,GAM)                                             00000120
      ZQ3 = DCMPLX(0.D0,GAM)                                            00000130
      ZSUM = ZDIGAM(ZQ1) + ZDIGAM(ZQ2) + CDLOG(ZQ3)*(-2.D0)             00000140
      PROD = 1.D0                                                       00000150
      L = RL                                                            00000160
      DO 100 I = 1,L                                                    00000170
      IF(L.EQ.0) GO TO 100                                              00000180
      PROD = PROD*( 1.D0 + I*I/(GAM*GAM) )                              00000190
100   CONTINUE                                                          00000200
      A1 = PROD                                                         00000210
      ZG = ZSUM*A1/(2.D0*PI)                                            00000220
      G = GG(1)                                                         00000230
      QQ = 1.D0 - DEXP(-2.D0*PI*GAM)                                    00000240
      B = A1 / QQ                                                       00000250
      RETURN                                                            00000260
      END

      FUNCTION GAMI (X,N)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION B(7)
      DIMENSION C(20)
      DIMENSION CCOR(4)
      DATA CCOR/ .20200 74006 59677 52D0,.03692 77526 92953 20D0,
     1 .0083492773817610 2D0, .00200839282608212D0 /
      DATA B /8.33333 33333 33333D-2,+8.33333 33333 33333D-3,
     1        3.96825 39682 53968D-3,+4.16666 66666 66667D-3,
     2        7.57575 7575757 576D-3,+2.10927 96092 79609D-2,
     3        8.33333 33333 33333D-2/
      DATA C /2.02056 90315 95943D-1, 3.69277 55143 36993D-2,
     1       8.34927 73819 22827D-3, 2.00839 28260 82214D-3,
     2       4.94188 60411 94646D-4, 1.22713 34757 84891D-4,
     3       3.05882 36307 02049D-5, 7.63719 76378 9976 D-6,
     4     1.90821 27165 5394 D-6, 4.76932 98678 781  D-7,
     5       1.19219 92596 531  D-7, 2.98035 0351465  D-8,
     6       7.45071 17898 4   D-9, 1.86265 97235 1 D-9,
     7       4.65662 90650    D-10, 1.16415 50173  D-10,
     8  2.910385044 D-11, 7.27595 984 D-12,
     9       1.81898 965    D-12,4.547473 8 D-13 /
      DATA EUL/.5772 15664 90153  29 D0/
c djh      DATA TES/1.D-10/
C CALCUL DE REEL(PSI(I*GAM)) - LOG(GAM)
C        N=1 SOMME 1/N(N2+GAM2)
C        N=2 DEVELOPPEMENT ASYMPTOTIQUE
C        N=3 CAS GAM INF A 2
      IF(N.GE.4) GO TO 8000
      A=DABS(X)
      AA=A*A
      IF(N.EQ.2) GO TO 50
      IF(N.EQ.3) GO TO 100
      S1=1.D0/(1.D0 + AA)
C     PRINT *, '  S1  ',S1
      DO 10 IK=2,100
      DS=AA + IK*IK
      DS=1.D0/(DS*IK)
      S1=S1+DS
10    CONTINUE
C     PRINT *,' X S1 DS TES ',X,S1,DS,TES
C     PRINT *,' ***** DEVELOPPEMENT NON CONVERGE'
      GO TO 30
CHG20    S1=S1+DS
30    CONTINUE
      RAP=DS/S1
      S2=DLOG(A)
C     PRINT *,' S1,LOG',S1,S2
      GAMI=S1*AA - EUL -S2
       RAP2=DABS(GAMI/(EUL+S2))
C     PRINT *,' CONVERG ANULATION ',RAP,RAP2
C     PRINT 1020,N,X,GAMI
1020  FORMAT(' OPTION N= ',I3,' X = ',D25.16,' GAMI = ',
     1D25.16,/)
      C1=0.D0
      X1=-X*X
      X2=-1.D0
      DO 5010 IC=1,4
      X2=X2*X1
      C2=X2*(C(IC)  - CCOR(IC))
      C1=C1 + C2
      GAMI= GAMI + C2
C     PRINT *,' C1 C2 GAMI ',C1,C2,GAMI
5010  CONTINUE
C     PRINT 1020,N,X,GAMI
      RETURN
50    AA=1.D0/AA
      A1=1.D0
      S1=0.D0
      DO 60 IK=1,7
      A1=A1*AA
      DS=B(IK)*A1
      S1= S1 + DS
      RAP=DS/S1
C     PRINT *,' S1 DS RAP ',S1,DS,RAP
60    CONTINUE
      GAMI=S1
C     PRINT 1020,N,X,GAMI
      RETURN
100   A1=-1.D0
      S1=0.D0
      DO 110 IK=1,20
      A1=-A1*AA
      DS=C(IK)*A1
      S1=S1 + DS
      RAP=DS/S1
C     PRINT *,' S1 DS RAP ', S1,DS,RAP
110   CONTINUE
      S2=DLOG(A)
      S3=1.D0/(1.D0 + AA)
      GAMI = 1.D0 -EUL -S2+S1 -S3
C     PRINT *,' S1 LOG S3 ',S1,S2,S3
      RAP=DABS(S1/GAMI)
      RAP=1.D0/RAP
C     PRINT *,'  ANULATION ',RAP
C     PRINT 1020,N,X,GAMI
      RETURN
8000   PRINT 8010,N
      GAMI=0.D0
8010      FORMAT('  FORMULE NON PROGRAMMEE N',
     +' ***********************',///)
      RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION D1IGAM(X)
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
      DIMENSION B(6)
      DATA PI/ 3.14159 26535 89793 D0/
      DATA B / +8.33333 33333 33333D-2, -8.33333 33333 33333D-3,
     1         +3.96825 39682 53968D-3, -4.16666 66666 66667D-3,
     2         +7.57575 75757 57576D-3, -2.10927 96092 79609D-2/
C CALCUL DE PSI(X)  X REEL
C PASSAGE X POSITIF PLUS GRAND QUE 15
C FORMULE ASYMPTOTIQUE A 7 TERMES
      A=DABS(X)
C     IF(DINT(X) + A)  ,4,
      XZ=DINT(X)+A
      IF(XZ.EQ.0.D0) GO TO 4
      V=A
      H=0.D0
      IF(A.GE.(15.D0)) GO TO 3
      N=14 - INT(A)
      H=1.D0/V
C     IF(N)  ,2,
      IF(N.EQ.0) GO TO 2
      DO 1 I=1,N
      V= V + 1.D0
1     H= H + 1.D0/V
2     V= V + 1.D0
3     R=1.D0/V**2
      D1IGAM=DLOG(V) -0.5D0/V -R*(B(1)+R*(B(2)+R*(B(3)+R*
     1 (B(4)+R*(B(5)+R*(B(6)+R*B(1))))))) - H
      IF(X . GE.(0.000D0) ) RETURN
      H=PI*A
      D1IGAM=D1IGAM + 1.D0/A + PI*DCOS(H)/DSIN(H)
      RETURN
C SORTIE ERREUR  :  X ENTIER NON POSITIF
4     PRINT 100,X
100   FORMAT(/////,131(1H*)//, '  *** DDIGAM
     1  ARGUMENT ENTIER NON NEGATIF =',D16.9,
     2  ' ***** ' )
      D1IGAM=0.D0
      RETURN
      END

      FUNCTION DIGAM(ARGG)                                              00000030
      IMPLICIT REAL*8(A-H,O-Z)                                          00000040
      DIGAM = 0.D0                                                      00000050
      ARG5 =-0.57721566490153286D0                                      00000060
      EN = 1.D0                                                         00000070
      ARG  = ARGG                                                       00000080
      ARG2 = ARGG                                                       00000090
1     IF(ARG2-40.D0) 2,3,3                                              00000100
2     DIGAM = DIGAM - 1.D0/ARG                                          00000110
      ARG = 1.D0+ ARG                                                   00000120
      ARG2 = 1.D0 + ARG2                                                00000130
      GO TO 1                                                           00000140
3     PSI = DLOG(ARG) - 1.D0/(2.D0*ARG) - 1.D0/(12.D0*ARG**2)           00000150
      PSI = PSI + 1.D0/(120.D0*ARG**4) - 1.D0/(252.D0*ARG**6)           00000160
      DIGAM = DIGAM + PSI                                               00000170
      RETURN                                                            00000180
      END                                                               00000190
c     COMPLEX FUNCTION ZDIGAM*16(ARG)                                   00000040
	function ZDIGAM(arg)
      IMPLICIT REAL*8(A-H,O-Y),COMPLEX*16(Z)                            00000050
      COMPLEX*16 ARG,ZDIGAM,cdlog                                       00000060
      COMPLEX*8 ARGG                                                    00000070
      REAL*8 INC                                                        00000080
      ZDIGAM = (0.D0,0.D0)                                              00000090
      ARG5 =-0.57721566490153286D0                                      00000100
      PI = 3.1415926535897932D0                                         00000110
      EN = 1.D0                                                         00000120
      ARGG = ARG                                                        00000130
      ARG2 = REAL(ARGG)                                                 00000140
      ARG3 = AIMAG(ARGG)                                                00000150
      IF(ARG3) 4,1,4                                                    00000160
1     IF(ARG2-40.D0) 2,3,3                                              00000170
2     ZDIGAM = ZDIGAM - 1.D0/ARG                                        00000180
      ARG = 1.D0+ ARG                                                   00000190
      ARG2 = 1.D0 + ARG2                                                00000200
      GO TO 1                                                           00000210
3     PSI = CDLOG(ARG)-1.D0/(2.D0*ARG)-1.D0/(12.D0*ARG**2)              00000220
      PSI=PSI +1.D0/(120.D0*ARG**4)-1.D0/(252.D0*ARG**6)                00000230
      ZDIGAM = ZDIGAM + PSI                                             00000240
      GO TO 12                                                          00000250
4     IF(ARG2) 5,7,6                                                    00000260
5     ZDIGAM = ZDIGAM - 1.D0/ARG                                        00000270
      ARG = ARG + 1.D0                                                  00000280
      ARG2 = ARG2 + 1.D0                                                00000290
      GO TO 4                                                           00000300
6     ARG = ARG - 1.D0                                                  00000310
      ARG2 = ARG2 - 1.D0                                                00000320
      ZDIGAM = ZDIGAM + 1.D0/ARG                                        00000330
      GO TO 4                                                           00000340
7     Y = CDABS(ARG)                                                    00000350
      ARG7 = PI*Y                                                       00000360
      ARG4 = 0.5D0/Y + (PI/2.D0)/DTANH(ARG7)                            00000370
      IF(Y-20.D0) 8,10,10                                               00000380
8     INC = Y*Y/(EN*(EN*EN+Y*Y))                                        00000390
      ARG5 = ARG5 + INC                                                 00000400
      IF(INC - 1.D-12) 11,11,9                                          00000410
9     EN = EN + 1.D0                                                    00000420
      GO TO 8                                                           00000430
10    ARG5 = 1.D0/(12.D0*Y**2) + 1.D0/(120.D0*Y**4)                     00000440
      ARG5 = ARG5 + 1.D0/(252.D0*Y**6)+ DLOG(Y)                         00000450
11    ZDIGAM = DCMPLX(ARG5,ARG4) + ZDIGAM                               00000460
C     XQ1 = REAL(ZDIGAM)                                                00000470
C     XQ2 = AIMAG(ZDIGAM)                                               00000480
12    continue                                                          00000490
	return
      END                                                               00000500      
      function rint(f,na,nb,nq,h)
c
c -- function rint performs a numerical integration over the function f,
c -----  assuming an equally spaced x-axis grid with step size h and a
c -----  total of nb mesh points.  Use na=1, nq=10, and nb>20 always.
c
	implicit none
      real*8 a,c,d,f,h, rint
      integer l,m,i,j,n,na,nb,nq

	dimension c(55),d(10),f(nb)
	c = (/1.d0,
     1  2.d0,1.d0,
     2  23.d0, 28.d0, 9.d0,
     3  25.d0, 20.d0, 31.d0, 8.d0,
     4  1413.d0, 1586.d0, 1104.d0, 1902.d0, 475.d0,
     5  1456.d0, 1333.d0, 1746.d0, 944.d0,1982.d0,459.d0,
     6  119585.d0, 130936.d0, 89437.d0, 177984.d0, 54851.d0,
     7  176648.d0, 36799.d0,
     8  122175.d0, 111080.d0, 156451.d0,46912.d0,220509.d0,
     9  29336.d0, 185153.d0, 35584.d0,
     a  7200319.d0, 7783754.d0, 5095890.d0,12489922.d0,-1020160.d0,
     b  16263486.d0, 261166.d0, 11532470.d0, 2082753.d0,
     c  7305728.d0, 6767167.d0, 9516362.d0,1053138.d0, 18554050.d0,
     d   -7084288.d0, 20306238.d0, -1471442.d0, 11965622.d0,
     e  2034625.d0/)
       d=(/  2.d0,    2.d0, 24.d0, 24.d0,1440.d0,1440.d0,120960.d0,
     &     120960.d0, 7257600.d0, 7257600.d0/)
	a=0.d0
	l=na
	m=nb
	i=nq*(nq+1)/2
	do 10 j=1,nq
	a=a+c(i)*(f(l)+f(m))
	l=l+1
	m=m-1
   10   i=i-1
	a=a/d(nq)
	do 20 n=l,m
   20   a=a+f(n)
	rint=a*h
	return
	end