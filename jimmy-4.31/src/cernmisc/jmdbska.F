*
* $Id: jmdbska.F 22 2006-11-19 13:14:22Z jmb $
*
* $Log$
* Revision 1.1  2005/07/29 16:15:46  jmb
* Include the various CERNLIB functions jimmy needs
*
* Revision 1.1.1.1  1996/04/01 15:02:08  mclareni
* Mathlib gen
*
*
      SUBROUTINE JMDBSKA(X,IA,JA,NL,B)
      IMPLICIT DOUBLE PRECISION (A-H,J,O-Z)
      SAVE
      LOGICAL LEX
      INTEGER JA

      CHARACTER NAME*(*),ENAM*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'BSKA/JMDBSKA', ENAM = 'EBSKA/DEBKA')
 
      PARAMETER (Z1 = 1, Z2 = 2, Z3 = 3, Z4 = 4)
      PARAMETER (Z12 = Z1/Z2, Z13 = Z1/Z3, Z14 = Z1/Z4, Z23 = Z2/Z3)
      PARAMETER (Z34 = Z3/Z4)
 
      DIMENSION B(0:*)
 
      PARAMETER (PI = 3.14159 26535 89793D0, PIV = PI/4)
 
      LEX=.FALSE.
      GO TO 9
 
      ENTRY JMDEBKA(X,IA,JA,NL,B)
      LEX=.TRUE.
 
    9 MODE=10*IA+JA
      N=NL-1
      U=2/X
      IF(LEX) THEN
       IF(X .LE. 0) THEN
        N=0
        WRITE(ERRTXT,101) X
        CALL JMMLPT(ENAM,'C341.1',ERRTXT)
       ELSEIF(NL .LT. 0 .OR. NL .GT. 100) THEN
        N=0
        WRITE(ERRTXT,103) NL
        CALL JMMLPT(ENAM,'C341.3',ERRTXT)
       ELSEIF(IA .EQ. 0) THEN
        A=0
        B(0)=JMDEBSK0(X)
        B(1)=JMDEBSK1(X)
       ELSEIF(MODE .EQ. 12) THEN
        A=Z12
        B(0)=SQRT(PIV*U)
        B(1)=B(0)*(1+A*U)
       ELSEIF(MODE .EQ. 13) THEN
        A=Z13
        B(0)=JMDEBKR3(X,1)
        B(1)=JMDEBKR3(X,2)+A*U*B(0)
       ELSEIF(MODE .EQ. 14) THEN
        A=Z14
        B(0)=JMDEBKR4(X,1)
        B(1)=JMDEBKR4(X,3)+A*U*B(0)
       ELSEIF(MODE .EQ. 23) THEN
        A=Z23
        B(0)=JMDEBKR3(X,2)
        B(1)=JMDEBKR3(X,1)+A*U*B(0)
       ELSEIF(MODE .EQ. 34) THEN
        A=Z34
        B(0)=JMDEBKR4(X,3)
        B(1)=JMDEBKR4(X,1)+A*U*B(0)
       ELSE
        N=0
        WRITE(ERRTXT,102) IA,JA
        CALL JMMLPT(ENAM,'C341.2',ERRTXT)
       ENDIF
      ELSE
       IF(X .LE. 0) THEN
        N=0
        WRITE(ERRTXT,101) X
        CALL JMMLPT(NAME,'C341.1',ERRTXT)
       ELSEIF(NL .LT. 0 .OR. NL .GT. 100) THEN
        N=0
        WRITE(ERRTXT,103) NL
        CALL JMMLPT(NAME,'C341.3',ERRTXT)
       ELSEIF(IA .EQ. 0) THEN
        A=0
        B(0)=JMDBESK0(X)
        B(1)=JMDBESK1(X)
       ELSEIF(MODE .EQ. 12) THEN
        A=Z12
        B(0)=EXP(-X)*SQRT(PIV*U)
        B(1)=B(0)*(1+A*U)
       ELSEIF(MODE .EQ. 13) THEN
        A=Z13
        B(0)=JMDBSKR3(X,1)
        B(1)=JMDBSKR3(X,2)+A*U*B(0)
       ELSEIF(MODE .EQ. 14) THEN
        A=Z14
        B(0)=JMDBSKR4(X,1)
        B(1)=JMDBSKR4(X,3)+A*U*B(0)
       ELSEIF(MODE .EQ. 23) THEN
        A=Z23
        B(0)=JMDBSKR3(X,2)
        B(1)=JMDBSKR3(X,1)+A*U*B(0)
       ELSEIF(MODE .EQ. 34) THEN
        A=Z34
        B(0)=JMDBSKR4(X,3)
        B(1)=JMDBSKR4(X,1)+A*U*B(0)
       ELSE
        N=0
        WRITE(ERRTXT,102) IA,JA
        CALL JMMLPT(NAME,'C341.2',ERRTXT)
       ENDIF
      ENDIF
      DO 1 IJ = 1,N
      A=A+1
    1 B(IJ+1)=B(IJ-1)+A*U*B(IJ)
      RETURN
  101 FORMAT('NON-POSITIVE ARGUMENT  X = ',E15.6)
  102 FORMAT('PAIR (IA,JA) = (',I5,I5,')  ILLEGAL')
  103 FORMAT('ILLEGAL  NL =',I5)
      END
