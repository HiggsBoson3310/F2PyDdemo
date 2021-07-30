CACCUCOULN.  COULN, A PROGRAM FOR EVALUATING NEGATIVE ENERGY COULOMB     ACCU0000
C1   FUNCTIONS.  C.J. NOBLE, I.J. THOMPSON.                              ACCU0000
CREF. IN COMP. PHYS. COMMUN. 33 (1984) 413                               ACCU0000
C//NB1 JOB (4401,NB),CJNOBLE,MSGLEVEL=(2,0),MSGCLASS=T,                  ACCU0001
C// NOTIFY=NB                                                            ACCU0002
C// EXEC FVCLG,TIME.C=(0,10),REGION.C=800K,PARM.C='OPTIMIZE(0)',         ACCU0003
C//    TIME.L=(0,2),REGION.G=250K,TIME.G=(0,15)                          ACCU0004
C//C.SYSIN DD *                                                          ACCU0005
C                                                                       ACCU0006
C***********************************************************************ACCU0007
C                                                                       ACCU0008
C     COULN : DRIVER FOR NEGATIVE ENERGY COULOMB FUNCTION CODE COULN    ACCU0009
C             ( EXPONENTIALLY DECAYING WHITTAKER FUNCTION )             ACCU0010
C                                                                       ACCU0011
C             VERSION OPTIMIZED FOR IBM DOUBLE PRECISION AND FORTRAN    ACCU0012
C             H EXTENDED OR VS FORTRAN COMPILERS                        ACCU0013
C                                                                       ACCU0014
C     C.J. NOBLE AND I.J. THOMPSON, S.E.R.C. DARESBURY LABORATORY       ACCU0015
C     ( JANUARY, 1984 )                                                 ACCU0016
C                                                                       ACCU0017
C***********************************************************************ACCU0018
C                                                                       ACCU0019
      IMPLICIT REAL*8 (A-H,O-Z)                                         ACCU0020
      DIMENSION XG(200),XA(200),XC(200),XD(200)                         ACCU0021
C                                                                       ACCU0022
 1000 FORMAT(2I10,2D16.8)                                               ACCU0023
 1010 FORMAT(I10,2D16.8)                                                ACCU0024
 1020 FORMAT(1H1,//,' TEST RUN OUTPUT FOR EXPONENTIALLY DECAYING',      ACCU0025
     X       ' WHITTAKER FUNCTION ROUTINE COULN',//,' ACC =',D16.8,/)   ACCU0026
 1030 FORMAT(//,I10,' E =',D16.8,' K =',D16.8,' Z =',I4,' ETA =',D16.8, ACCU0027
     X       ' RHO =',D16.8,' RATIO =',D16.8,/)                         ACCU0028
 1040 FORMAT(' L =',I3,' FX =',D20.12,' FXP =',D20.12,' FNORM =',D20.12,ACCU0029
     X       ' EFX =',D8.3,' IERR =',I3,' NLAST =',I6)                  ACCU0030
C                                                                       ACCU0031
      IREAD=5                                                           ACCU0032
      IWRITE=6                                                          ACCU0033
C                                                                       ACCU0034
   10 READ(IREAD,1000) LMAX,IZ,R,ACC                                    ACCU0035
      IF(LMAX .LT. 0) GO TO 40                                          ACCU0036
      READ(IREAD,1010) NUME,EMIN,EINC                                   ACCU0037
C                                                                       ACCU0038
      WRITE(IWRITE,1020) ACC                                            ACCU0039
C                                                                       ACCU0040
      LMAX1=LMAX+1                                                      ACCU0041
      E=EMIN                                                            ACCU0042
      DO 30 I=1,NUME                                                    ACCU0043
C                                                                       ACCU0044
      FK=DSQRT(-E)                                                      ACCU0045
      Z=DFLOAT(IZ)                                                      ACCU0046
      ETA=-Z/FK                                                         ACCU0047
      RHO=Z*R                                                           ACCU0048
      RATIO=DABS(R/ETA)                                                 ACCU0049
C                                                                       ACCU0050
      WRITE(IWRITE,1030) I,E,FK,IZ,ETA,RHO,RATIO                        ACCU0051
C                                                                       ACCU0052
      DO 20 L1=1,LMAX1                                                  ACCU0053
      L=L1-1                                                            ACCU0054
C                                                                       ACCU0055
      CALL COULN(L,Z,E,R,FX,FXP,ACC,EFX,XG,XA,XC,XD,IERR,NLAST,FNORM)   ACCU0056
C                                                                       ACCU0057
      WRITE(IWRITE,1040) L,FX,FXP,FNORM,EFX,IERR,NLAST                  ACCU0058
C                                                                       ACCU0059
   20 CONTINUE                                                          ACCU0060
C                                                                       ACCU0061
      E=E+EINC                                                          ACCU0062
   30 CONTINUE                                                          ACCU0063
C                                                                       ACCU0064
      GO TO 10                                                          ACCU0065
C                                                                       ACCU0066
   40 STOP                                                              ACCU0067
      END                                                               ACCU0068
      SUBROUTINE COULN(LX,Z,E,R,FX,FXP,ACC,EFX,XG,XA,XC,XD,IERR,NLAST,  ACCU0069
     X                 FNORM)                                           ACCU0070
      IMPLICIT REAL*8 (A-H,O-Z)                                         ACCU0071
      REAL*8 NINE                                                       ACCU0072
C                                                                       ACCU0073
C***********************************************************************ACCU0074
C                                                                       ACCU0075
C     COULN COMPUTES THE EXPONENTIALLY DECAYING WHITTACKER FUNCTION     ACCU0076
C           AND ITS DERIVATIVE USING THE ASYMPTOTIC EXPANSION CONVERTED ACCU0077
C           TO A CONTINUED FRACTION                                     ACCU0078
C                                                                       ACCU0079
C                FX = W( KAPPA, MU, ZC )                                ACCU0080
C                                                                       ACCU0081
C                KAPPA = 1 / C = Z / DSQRT( -E )                        ACCU0082
C                MU    = LX + 1/2                                       ACCU0083
C                ZC    = 2 * DSQRT( -E ) * R                            ACCU0084
C                                                                       ACCU0085
C      REFERENCE : HANDBOOK OF MATHEMATICAL FUNCTIONS                   ACCU0086
C                  M. ABRAMOWITZ AND I.A. STEGUN                        ACCU0087
C                  ( DOVER, NEW YORK, 1965 )   EQUATION : 13.1.33       ACCU0088
C                                                                       ACCU0089
C***********************************************************************ACCU0090
C                                                                       ACCU0091
C     LX      ORBITAL ANGULAR MOMENTUM ( .GE. ZERO )                    ACCU0092
C     Z       CHARGE  ( NONZERO )                                       ACCU0093
C     E       ENERGY IN RYDBERGS ( .LT. ZERO )                          ACCU0094
C     R       RADIAL COORDINATE IN BOHR UNITS                           ACCU0095
C     ACC     RELATIVE ACCURACY REQUIRED                                ACCU0096
C                                                                       ACCU0097
C     FX      WHITTAKER FUNCTION                                        ACCU0098
C     FXP     DERIVATIVE OF WHITTAKER FUNCTION WITH RESPECT TO R        ACCU0099
C     FNORM   NORMALIZATION FACTOR MULTIPLYING BOTH FX AND FXP          ACCU0100
C     EFX     ESTIMATE OF RELATIVE ERROR IN VALUE OF FX RETURNED        ACCU0101
C                                                                       ACCU0102
C     IERR    = 0   FOR NORMAL RETURN WITH REQUESTED ACCURACY           ACCU0103
C             = 1   VALUE RETURNED WITH REDUCED ACCURACY                ACCU0104
C             = 2   INPUT ERROR, FX, FXP ZEROED                         ACCU0105
C             = 3   NORMAL RETURN, SERIES TERMINATES                    ACCU0106
C     NLAST   NUMBER OF TERMS USED IN THE ASYMPTOTIC EXPANSION          ACCU0107
C                                                                       ACCU0108
C     XG,XA,XC,XD    WORK ARRAYS, LENGTH OF EACH = KKKMAX * NTERM       ACCU0109
C                    COEFFICIENTS ARE COMPUTED IN UP TO KKKMAX BLOCKS   ACCU0110
C                    EACH BLOCK CONSISTS OF NTERM COEFFICIENTS          ACCU0111
C                                                                       ACCU0112
C     MACHINE DEPENDENT CONSTANTS :                                     ACCU0113
C                                                                       ACCU0114
C     TMAX     MAGNITUDE OF LARGEST COEFFICIENT TO BE USED              ACCU0115
C     SCAL     MAGNITUDE OF SMALLEST NUMBER TO BE USED                  ACCU0116
C              DETERMINES THE SCALING OF THE SERIES COEFFICIENTS        ACCU0117
C     IXMAX    ARGUMENT RANGE PARAMETER OF THE INTRINSIC EXPONENTIAL    ACCU0118
C              FUNCTION                                                 ACCU0119
C                                                                       ACCU0120
C     ROUTINES CALLED :       DFRACT, HSUM                              ACCU0121
C     INTRINSIC FUNCTIONS :   IDINT, MAX0, DABS, DEXP, DFLOAT, DLOG,    ACCU0122
C                             DMIN1, DSQRT                              ACCU0123
C                                                                       ACCU0124
C***********************************************************************ACCU0125
C                                                                       ACCU0126
      DIMENSION XG(2),XA(2),XC(2),XD(2)                                 ACCU0127
C                                                                       ACCU0128
      DATA ZERO/0.0D0/,ONE/1.0D0/,TWO/2.0D0/,NINE/9.0D0/,HALF/0.5D0/    ACCU0129
      DATA SCAL/1.0D-50/,TMAX/1.0D+60/,EPS/1.0D-10/,AK5/1.0D+05/        ACCU0130
      DATA IXMAX/150/,KKKMAX/20/,NTERM/10/                              ACCU0131
C                                                                       ACCU0132
      IERR=0                                                            ACCU0133
      ITERM=0                                                           ACCU0134
      FNORM=ONE                                                         ACCU0135
      XMAX=DFLOAT(IXMAX)                                                ACCU0136
C                                                                       ACCU0137
      IF(E .GE. ZERO .OR. R .LE. ZERO .OR. Z .EQ. ZERO .OR. LX .LT. 0)  ACCU0138
     X   GO TO 140                                                      ACCU0139
C                                                                       ACCU0140
      AK=DSQRT(-E)                                                      ACCU0141
      ZI=ONE/Z                                                          ACCU0142
      C=AK*ZI                                                           ACCU0143
      X=Z*R                                                             ACCU0144
      XI=ONE/X                                                          ACCU0145
      XI2=XI*XI                                                         ACCU0146
      CX=AK*R                                                           ACCU0147
      SCALE=SCAL*DEXP(DMIN1(CX,-DLOG(SCAL)))                            ACCU0148
      RSCALE=ONE/SCALE                                                  ACCU0149
      C1=ONE/C                                                          ACCU0150
      NC=IDINT(C1+HALF)                                                 ACCU0151
C                                                                       ACCU0152
      NFIRST=1                                                          ACCU0153
      NLAST=NTERM                                                       ACCU0154
C                                                                       ACCU0155
C     CHECK FOR FINITE SERIES SPECIAL CASES                             ACCU0156
C                                                                       ACCU0157
      IF(DABS(DFLOAT(NC)-C1) .LT. EPS .AND. NC .GE. -LX) ITERM=1        ACCU0158
C                                                                       ACCU0159
      L=LX                                                              ACCU0160
      LM=0                                                              ACCU0161
      IF(ITERM .EQ. 1) LM=NC+1                                          ACCU0162
      IF(X .GE. ZERO) LM=MAX0(LM,IDINT(DSQRT(X)+HALF))                  ACCU0163
      IF(LX .GT. LM) L=LM                                               ACCU0164
C                                                                       ACCU0165
      L1=L+1                                                            ACCU0166
      FL=DFLOAT(L)                                                      ACCU0167
      FL1=DFLOAT(L1)                                                    ACCU0168
C                                                                       ACCU0169
      IF(ITERM .EQ. 0) GO TO 20                                         ACCU0170
      IF(NC .GE. L1) GO TO 10                                           ACCU0171
      NLAST=L+NC                                                        ACCU0172
      GO TO 20                                                          ACCU0173
   10 NLAST=NC-L1                                                       ACCU0174
C                                                                       ACCU0175
C     COMPUTE EXPANSION COEFFICIENTS FOR ASYMPTOTIC SERIES L = LM       ACCU0176
C     ( IBID. EQUATION : 13.5.2 )                                       ACCU0177
C                                                                       ACCU0178
   20 Q=SCALE                                                           ACCU0179
      XG0=Q                                                             ACCU0180
      G1=FL-C1                                                          ACCU0181
      G2=-FL1-C1                                                        ACCU0182
      G3=-HALF*C1                                                       ACCU0183
C                                                                       ACCU0184
      DO 90 KKK=1,KKKMAX                                                ACCU0185
C                                                                       ACCU0186
      IF(NLAST .LT. NFIRST) GO TO 50                                    ACCU0187
      DO 30 N=NFIRST,NLAST                                              ACCU0188
      AN=DFLOAT(N)                                                      ACCU0189
      Q=(G1+AN)*(G2+AN)*G3*Q/AN                                         ACCU0190
      IF(DABS(Q) .GT. TMAX) GO TO 100                                   ACCU0191
      XG(N)=Q                                                           ACCU0192
   30 CONTINUE                                                          ACCU0193
      IF(ITERM .EQ. 1) GO TO 50                                         ACCU0194
C                                                                       ACCU0195
C     OBTAIN COEFFICIENTS OF THE CONTINUED FRACTION USING THE NESBET    ACCU0196
C     ALGORITHM                                                         ACCU0197
C                                                                       ACCU0198
      CALL DFRACT(NLAST,XG0,XG,XA,XC,XD)                                ACCU0199
C                                                                       ACCU0200
      KMAX=(NLAST+1)/2                                                  ACCU0201
      KMAX1=KMAX-1                                                      ACCU0202
      NC=KMAX+KMAX-1                                                    ACCU0203
C                                                                       ACCU0204
C     EVALUATE THE CONTINUED FRACTION USING THE RATIONAL FORM           ACCU0205
C                                                                       ACCU0206
      NCP=NC+1                                                          ACCU0207
      NQ=NLAST+1                                                        ACCU0208
      SN1=XC(NC)                                                        ACCU0209
      SN2=XC(NCP)                                                       ACCU0210
      AI=DFLOAT(KMAX)                                                   ACCU0211
      SNP1=AI*SN1                                                       ACCU0212
      SNP2=AI*SN2                                                       ACCU0213
C                                                                       ACCU0214
      SD1=XD(NC)                                                        ACCU0215
      SD2=XD(NCP)                                                       ACCU0216
      SDP1=AI*SD1                                                       ACCU0217
      SDP2=AI*SD2                                                       ACCU0218
C                                                                       ACCU0219
      DO 40 K=1,KMAX1                                                   ACCU0220
      NC=NC-2                                                           ACCU0221
      NCP=NC+1                                                          ACCU0222
      AI=DFLOAT(KMAX-K)                                                 ACCU0223
C                                                                       ACCU0224
      SN1=XC(NC)+XI*SN1                                                 ACCU0225
      SN2=XC(NCP)+XI*SN2                                                ACCU0226
      SNP1=AI*XC(NC)+XI*SNP1                                            ACCU0227
      SNP2=AI*XC(NCP)+XI*SNP2                                           ACCU0228
C                                                                       ACCU0229
      SD1=XD(NC)+XI*SD1                                                 ACCU0230
      SD2=XD(NCP)+XI*SD2                                                ACCU0231
      SDP1=AI*XD(NC)+XI*SDP1                                            ACCU0232
      SDP2=AI*XD(NCP)+XI*SDP2                                           ACCU0233
C                                                                       ACCU0234
      IF(DABS(SN1)+DABS(SN2)+DABS(SD1)+DABS(SD2) .GT.  NINE*TMAX*X)     ACCU0235
     X   GO TO 100                                                      ACCU0236
      IF(DABS(SNP1)+DABS(SNP2)+DABS(SDP1)+DABS(SDP2) .GT. AK5*TMAX*X*X) ACCU0237
     X   GO TO 100                                                      ACCU0238
   40 CONTINUE                                                          ACCU0239
C                                                                       ACCU0240
      SN1=XG0+XI*SN1                                                    ACCU0241
      SN2=XG0+XI*SN2                                                    ACCU0242
      SNP1=-XI2*SNP1                                                    ACCU0243
      SNP2=-XI2*SNP2                                                    ACCU0244
C                                                                       ACCU0245
      SD1=ONE+XI*SD1                                                    ACCU0246
      SD2=ONE+XI*SD2                                                    ACCU0247
      SDP1=-XI2*SDP1                                                    ACCU0248
      SDP2=-XI2*SDP2                                                    ACCU0249
C                                                                       ACCU0250
C     OBTAIN FUNCTIONS AND ERROR ESTIMATES                              ACCU0251
C                                                                       ACCU0252
      F=SN1/SD1                                                         ACCU0253
      F1=SN2/SD2                                                        ACCU0254
      EF=DABS(F-F1)                                                     ACCU0255
C                                                                       ACCU0256
      FP=(SNP1-F*SDP1)/SD1                                              ACCU0257
      FP1=(SNP2-F1*SDP2)/SD2                                            ACCU0258
      EFP=DABS(FP-FP1)                                                  ACCU0259
      GO TO 60                                                          ACCU0260
C                                                                       ACCU0261
C     EVALUATE FINITE SERIES CASES USING HORNERS ALGORITHM              ACCU0262
C                                                                       ACCU0263
   50 CALL HSUM(NLAST,1,X,XG0,XG,F,FP)                                  ACCU0264
      EFX=ZERO                                                          ACCU0265
      EFXP=ZERO                                                         ACCU0266
      IERR=3                                                            ACCU0267
C                                                                       ACCU0268
   60 XARG=-CX+C1*DLOG(TWO*CX)                                          ACCU0269
      IF(XARG .LT. XMAX) GO TO 70                                       ACCU0270
      XARG=XARG-XMAX                                                    ACCU0271
      FNORM=DEXP(-XMAX)                                                 ACCU0272
      GO TO 80                                                          ACCU0273
C                                                                       ACCU0274
   70 IF(XARG .GT. -XMAX) GO TO 80                                      ACCU0275
      XARG=XARG+XMAX                                                    ACCU0276
      FNORM=DEXP(XMAX)                                                  ACCU0277
C                                                                       ACCU0278
   80 FAC=DEXP(XARG)                                                    ACCU0279
      FX=FAC*F                                                          ACCU0280
      FXP=Z*((ONE/CX-C)*F+FP)*FAC                                       ACCU0281
      IF(ITERM .EQ. 1) GO TO 110                                        ACCU0282
C                                                                       ACCU0283
      EFX=EF/DMAX1(DABS(F),SCAL)                                        ACCU0284
      EFXP=EFP/DMAX1(DABS(FP),SCAL)                                     ACCU0285
      IF(DABS(EFX) .LT. ACC .AND. DABS(EFXP) .LT. ACC) GO TO 110        ACCU0286
C                                                                       ACCU0287
      NFIRST=NLAST+1                                                    ACCU0288
      NLAST=NLAST+NTERM                                                 ACCU0289
C                                                                       ACCU0290
   90 CONTINUE                                                          ACCU0291
C                                                                       ACCU0292
  100 IERR=1                                                            ACCU0293
      NLAST=NLAST-NTERM                                                 ACCU0294
      IF(KKK .EQ. 1) GO TO 140                                          ACCU0295
C                                                                       ACCU0296
  110 FX=RSCALE*FX                                                      ACCU0297
      FXP=RSCALE*FXP                                                    ACCU0298
C                                                                       ACCU0299
C     UPWARDS RECURRENCE ON L IF LX GREATER THAN LM                     ACCU0300
C                                                                       ACCU0301
      IF(LX .LE. LM) GO TO 130                                          ACCU0302
C                                                                       ACCU0303
      VP=FX                                                             ACCU0304
      VQ=FXP                                                            ACCU0305
      DO 120 L=L1,LX                                                    ACCU0306
      AL=DFLOAT(L)                                                      ACCU0307
      AL2=DFLOAT(L*L)                                                   ACCU0308
      FX=((ONE-AL2*XI)*VP+AL*VQ*ZI)/(ONE-C*AL)                          ACCU0309
      FXP=Z*((ONE-AL2*XI)*FX-(ONE+C*AL)*VP)/AL                          ACCU0310
      VP=FX                                                             ACCU0311
      VQ=FXP                                                            ACCU0312
  120 CONTINUE                                                          ACCU0313
C                                                                       ACCU0314
  130 RETURN                                                            ACCU0315
C                                                                       ACCU0316
  140 IERR=2                                                            ACCU0317
      FX=ZERO                                                           ACCU0318
      FXP=ZERO                                                          ACCU0319
      RETURN                                                            ACCU0320
C                                                                       ACCU0321
      END                                                               ACCU0322
      SUBROUTINE DFRACT(N,G0,G,A,C,D)                                   ACCU0323
      IMPLICIT REAL*8 (A-H,O-Z)                                         ACCU0324
C                                                                       ACCU0325
C***********************************************************************ACCU0326
C                                                                       ACCU0327
C  DFRACT CONVERTS POLYNOMIAL G TO THE CORRESPONDING CONTINUED          ACCU0328
C         FRACTION, IN STIELTJES FORM WITH COEFFICIENTS A               ACCU0329
C                                                                       ACCU0330
C   G(Z) = G0+G1*Z+G2*Z**2+G3*Z**3+...+GN*Z**N                          ACCU0331
C                                                                       ACCU0332
C   GN(Z)= CN(Z)/DN(Z)                                                  ACCU0333
C        = G0/1+ A1*Z/1+A2*Z/1+A3*Z/1+.../1+AN*Z                        ACCU0334
C                                                                       ACCU0335
C  DATA:                                                                ACCU0336
C   N     ORDER OF G, EVEN OR ODD       INPUT                           ACCU0337
C   G0    CONSTANT TERM                 INPUT                           ACCU0338
C   G     VECTOR G(K), K=1,N            INPUT                           ACCU0339
C   A     VECTOR A(K), K=1,N            OUTPUT                          ACCU0340
C   C     OUTPUT: NUMERATOR   POLYNOMIALS, CONVERGENTS N AND N-1        ACCU0341
C   D     OUTPUT: DENOMINATOR POLYNOMIALS, CONVERGENTS N AND N-1        ACCU0342
C         CALLER PROVIDES SPACE FOR A,C,D                               ACCU0343
C   C AND D ARRAYS CONTAIN NT=2*((N+1)/2) ELEMENTS EACH                 ACCU0344
C   STORAGE CONVENTIONS FOR C AND D:                                    ACCU0345
C     C0=G0 AND D0=1.0                                                  ACCU0346
C     COEFFICIENTS IN SEQUENCE, WITH M=N-1 AND K=(N+1)/2 :              ACCU0347
C       CN1,CM1,CN2,CM2,...CNI,CMI,...  I=1,K                           ACCU0348
C       DN1,DM1,DN2,DM2,...DNI,DMI,...  I=1,K                           ACCU0349
C     NOTE THAT DNK=0.0 IF N IS ODD                                     ACCU0350
C                                                                       ACCU0351
C   ALGORITHM: RK NESBET, 82.10.27                                      ACCU0352
C                                                                       ACCU0353
C***********************************************************************ACCU0354
C                                                                       ACCU0355
      DIMENSION G(1),A(1),C(1),D(1)                                     ACCU0356
      DATA ZERO/0.0D0/,ONE/1.0D0/                                       ACCU0357
C                                                                       ACCU0358
      NT=2*((N+1)/2)                                                    ACCU0359
      DO 10  I=1,NT                                                     ACCU0360
      C(I)=ZERO                                                         ACCU0361
   10 D(I)=ZERO                                                         ACCU0362
      CT=ZERO                                                           ACCU0363
      DT=ONE                                                            ACCU0364
      DN= G0                                                            ACCU0365
      DO 40 K=1,N                                                       ACCU0366
      DD=DN                                                             ACCU0367
      DN=G(K)                                                           ACCU0368
      DO 20 I=1,K,2                                                     ACCU0369
      IF((I+1)/2.GE.K) GO TO 20                                         ACCU0370
      DN=DN+G(K-(I+1)/2)*D(I)                                           ACCU0371
   20 CONTINUE                                                          ACCU0372
      A(K)=-DN/DD                                                       ACCU0373
C                                                                       ACCU0374
C     PLANT AK=0.0 AND RETURN IF SEQUENCE TRUNCATES                     ACCU0375
C                                                                       ACCU0376
      IF(A(K) .EQ. ZERO) RETURN                                         ACCU0377
      DO 30 I=1,K,2                                                     ACCU0378
      CI=C(I)                                                           ACCU0379
      C(I)=CI+A(K)*CT                                                   ACCU0380
      CT=C(I+1)                                                         ACCU0381
      C(I+1)=CI                                                         ACCU0382
      DI=D(I)                                                           ACCU0383
      D(I)=DI+A(K)*DT                                                   ACCU0384
      DT=D(I+1)                                                         ACCU0385
      D(I+1)=DI                                                         ACCU0386
   30 CONTINUE                                                          ACCU0387
      CT=G0                                                             ACCU0388
      DT=ONE                                                            ACCU0389
   40 CONTINUE                                                          ACCU0390
C                                                                       ACCU0391
      RETURN                                                            ACCU0392
      END                                                               ACCU0393
      SUBROUTINE HSUM(N,M,X,C0,C,S,SP)                                  ACCU0394
      IMPLICIT REAL*8 (A-H,O-Z)                                         ACCU0395
C                                                                       ACCU0396
C***********************************************************************ACCU0397
C                                                                       ACCU0398
C     HSUM EVALUATES THE SUMS                                           ACCU0399
C                                                                       ACCU0400
C           S  = C0 + C(1)/X + C(2)/X**2 + ... + C(N)/X**N              ACCU0401
C           SP =  - ( C(1)/X + 2*C(2)/X**2 + ... + N*C(N)/X**N ) / X    ACCU0402
C                                                                       ACCU0403
C***********************************************************************ACCU0404
C                                                                       ACCU0405
      DIMENSION C(2)                                                    ACCU0406
      DATA ZERO/0.0D0/,ONE/1.0D0/                                       ACCU0407
C                                                                       ACCU0408
      IF(N .EQ. 0) GO TO 30                                             ACCU0409
C                                                                       ACCU0410
      N1=N-1                                                            ACCU0411
      NC=N                                                              ACCU0412
C                                                                       ACCU0413
      XI=ONE/X                                                          ACCU0414
C                                                                       ACCU0415
      NQ=N1*M+1                                                         ACCU0416
      S=C(NQ)                                                           ACCU0417
      AI=DFLOAT(N)                                                      ACCU0418
      SP=AI*S                                                           ACCU0419
      IF(N .EQ. 1) GO TO 20                                             ACCU0420
C                                                                       ACCU0421
      DO 10 K=1,N1                                                      ACCU0422
      NC=NC-1                                                           ACCU0423
      NQ=NQ-M                                                           ACCU0424
      AI=DFLOAT(NC)                                                     ACCU0425
      S=C(NQ)+XI*S                                                      ACCU0426
      SP=AI*C(NQ)+XI*SP                                                 ACCU0427
   10 CONTINUE                                                          ACCU0428
C                                                                       ACCU0429
   20 S=C0+XI*S                                                         ACCU0430
      SP=-XI*XI*SP                                                      ACCU0431
      RETURN                                                            ACCU0432
C                                                                       ACCU0433
   30 S=C0                                                              ACCU0434
      SP=ZERO                                                           ACCU0435
      RETURN                                                            ACCU0436
C                                                                       ACCU0437
      END                                                               ACCU0438
C//G.SYSIN DD *                                                          ACCU0439
C        20         2       10.00D+00       1.000D-12                    ACCU0440
C         3     -10.000D+00     2.99000D+00                              ACCU0441
C        10         2        1.00D+00       1.000D-12                    ACCU0442
C         2      -0.210D+00     0.10000D+00                              ACCU0443
C        30        22        2.50D+00       1.000D-12                    ACCU0444
C         1 -0.25397958D+01     0.05000D+00                              ACCU0445
C        30         2        1.50D+00       1.000D-12                    ACCU0446
C         1 -0.48000000D+01     0.05000D+00                              ACCU0447
C       -10                                                              ACCU0448
C//                                                                      ACCU0449
C                                                                        ACCU****
