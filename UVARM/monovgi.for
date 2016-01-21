C     Vincente Pericoli
C     UC Davis
C     14 Jan 2016
C     
C     This subroutine calculates stress triaxiality and monotonic VGI
C
C-----------------------------------------------------------------------
C     Submit via command line using:
C     abaqusfortran job=INPUT_FILE cpus=2 mp_mode=mpi user=monovgi.for
C     
C     need to setup abaqusfortran batch file to execute 
C     ifortvars_intel64.bat before launching abaqus
C-----------------------------------------------------------------------

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)
C
C     Required by ABAQUS:    
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
C     My Variables:
      DIMENSION XNEWVAL(2)
      PARAMETER ZTOL=0.00001
C      
C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.
C
C-----------------------------------------------------------------------
C     Define Variables
C     UVAR(1) : Stress Triaxiality
C     UVAR(2) : PEEQ
C     UVAR(3) : Monotonic VGI
C     XNEWVAL(1) : Current Stress Triaxiality
C     XNEWVAL(2) : Current PEEQ
C-----------------------------------------------------------------------
C
C     Obtain stress values, and calculate Triaxiality
C     note that SINV is all stress invariant components 
C     i.e. (MISES, TRESC, PRESS, INV3)
C     JRCD is an implicit integer.
      CALL GETVRM('SINV',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,MATLAYO,
     & LACCFLA)
C      
      IF (ARRAY(1).LE.ZTOL) THEN
        XNEWVAL(1) = 0.0
      ELSE
        XNEWVAL(1) = -ARRAY(3)/ARRAY(1)
      ENDIF
C
C     Obtain the PEEQ
C     This is basically a hack so that we can easily save/access the
C     PEEQ from the previous time increment.
      CALL GETVRM('PE',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,MATLAYO,
     & LACCFLA)
      XNEWVAL(2) = ARRAY(7)
C
C     Calculate the Monotonic VGI
      UVAR(3) = UVAR(3) + 0.5*(XNEWVAL(2) - UVAR(2)) * 
     &                    EXP(ABS(1.5*XNEWVAL(1)))
C
C     Update user-variables
      UVAR(1) = XNEWVAL(1)
      UVAR(2) = XNEWVAL(2)
      RETURN
      END