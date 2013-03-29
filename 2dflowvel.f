************************************************************************
*   2-D UNSTEADY NONLINEAR TIDAL & WIND-DRIVEN COASTAL CIRCULATION     *
************************************************************************
*       This program solves the non-linear, depth averaged conserv-    *
*       ation equations using a numerical scheme due to Koutitas       *
*       (1988) as translated into FORTRAN by R. Slingerland.           *
*                                                                      *
************************************************************************
*                The set up constants are:                             *
*   DT,DX  --- Time (s) and space (m) descretisation steps             *
*   CS     --- Wind friction coefficient (dimensionless)               *
*   CF     --- Chezy Bed friction coefficient (units?)                 *
*   WX,WY  --- Wind velocity components in x and y (m/s)               *
*   f      --- Coriolis parameter (1/s)                                *
*   IM     --- Max number of grid points along x direction             *
*   JM     --- Max number of grid points along y direction             *
*   NM     --- Max number of time steps desired                        *
*   KB     --- Number of coastal and open boundary nodes               *
*   DAT    --- Dependent variables are saved every dat timesteps       *
*   AMPL   --- Amplitude of the incident long waves (m)                *
*   PER    --- Period of the incident long waves (s)                   *
*   IS(J)  --- Starting node number of the computation field in the jth*
*              row                                                     *
*   IE(J)  --- Ending node number of the computation field in the jth  *
*              row                                                     *
*   IB(K),                                                             *
*   JB(K),                                                             *
*   NB(K)  --- IB and JB are the coordinates of the kth boundary       *
*              node; nb is the cell type                               *
*   H(i,j) --- Still water depths (m)                                  *
*   Z      --- free surface elevation with respect to the SWL (m)      *
*   U,V    --- Vertically averaged velocities (m/s)                    *
************************************************************************
************************************************************************
      Implicit None
      Character*80 formt
      Real*8 u(50,50),un(50,50),v(50,50),vn(50,50),h(50,50),snx,sny,
     #z(25,25),uu,vv,hm,tsx,tsy,ts,tbx,tby,cf,cs,wx,wy,dx,dt,ad,
     #f,g,t,dat,time,ubar(50,50),vbar(50,50),ampl,per,zr1,zr2,zi(50,50),
     #c,pi,L
      Integer*4 i,j,k,im,jm,nm,is(50),ie(50),ib(50),jb(50),nb(50),kb,
     #n,kk,p,q
C-----------------------------------------------------------------------
C    Open the input and output files
C-----------------------------------------------------------------------
      OPEN(10, FILE='input.dat')
      OPEN(11, FILE='H.DATA')
      OPEN(12, FILE='U.DATA')
      OPEN(13, FILE='V.DATA')
      OPEN(14, FILE='Z.DATA')
C-----------------------------------------------------------------------
C  Read the Input Variables
C-----------------------------------------------------------------------
      READ(10,*)DT,DX,CS,CF,WX,WY,F,IM,JM,NM,KB,DAT,AMPL,PER
      DO J = 1, JM - 1
         READ (10, *) IS(J), IE(J)
      END DO
      DO K = 1, KB
         READ (10, *) IB(K), JB(K), NB(K)
      END DO
      DO J = 1, JM - 1
         DO I = 1, IM
            H(I,J) = 415.0 - 17.0*J
C-----------------------------------------------------------------------
C  If bathymetry is to be read in, uncomment the next line.
C-----------------------------------------------------------------------
C          read(10,*)(h(i,j),i=is(j)-1,ie(j)+1)
         END DO
      END DO
C-----------------------------------------------------------------------
C  Define Wind Stresses and Zero Out Velocity Variables
C-----------------------------------------------------------------------
      TSX = CS*WX*SQRT(WX**2+WY**2)
      TSY = CS*WY*SQRT(WX**2+WY**2)
      SNX = SIGN(1.,TSX)
      SNY = SIGN(1.,TSY)
      TS = SQRT(TSX**2+TSY**2)
      G = 9.81
      PI = 3.1459
      T = 0
      KK = 100
      DO I = 1, IM
         DO J = 1, JM
            U(I,J) = 0.0
            UN(I,J) = 0.0
            V(I,J) = 0.0
            VN(I,J) = 0.0
            Z(I,J) = 0.0
            ZI(I,J) = 0.0
         END DO
      END DO
C-----------------------------------------------------------------------
C  Start the Main Time Loop
C-----------------------------------------------------------------------
      DO 33 N = 1, NM
         T = T + DT
C-----------------------------------------------------------------------
C  Solve the Continuity Equation for Updated WSE (z) at n+1/2
C-----------------------------------------------------------------------
         DO J = 2, JM - 2
            DO I = IS(J), IE(J)
               Z(I,J) = Z(I,J) - DT/2.0/DX*(U(I+1,J)*(H(I,J)+H(I+1,J))-U
     #            (I,J)*(H(I,J)+H(I-1,J))+V(I,J+1)*(H(I,J+1)+H(I,J))-V(I
     #            ,J)*(H(I,J)+H(I,J-1)))
            END DO
         END DO
C-----------------------------------------------------------------------
C  Solve x-directed Momentum Equation for u at n+1
C-----------------------------------------------------------------------
         DO J = 2, JM - 2
            DO I = IS(J) + 1, IE(J)
               VV = (V(I,J)+V(I-1,J)+V(I,J+1)+V(I-1,J+1))/4.0
               HM = (H(I,J)+H(I-1,J))/2.0
               TBX = CF*U(I,J)*SQRT(VV**2+U(I,J)**2)
               AD = ((U(I,J)+U(I+1,J))**2-(U(I,J)+U(I-1,J))**2)/8/DX + 
     #            VV*(U(I,J+1)-U(I,J-1))/2/DX
               UN(I,J) = U(I,J) - DT*(AD+G*(Z(I,J)-Z(I-1,J))/DX-F*VV-(
     #            TSX-TBX)/HM)
            END DO
         END DO
C-----------------------------------------------------------------------
C  Solve y-directed Momentum Equation for v at n+1
C-----------------------------------------------------------------------
         DO J = 3, JM - 2
            DO I = IS(J), IE(J)
               UU = (U(I,J)+U(I+1,J)+U(I,J-1)+U(I+1,J-1))/4.0
               HM = (H(I,J)+H(I,J-1))/2.0
               TBY = CF*V(I,J)*SQRT(UU**2+V(I,J)**2)
               AD = ((V(I,J+1)+V(I,J))**2-(V(I,J)+V(I,J-1))**2)/8/DX + 
     #            UU*(V(I+1,J)-V(I-1,J))/2/DX
               VN(I,J) = V(I,J) - DT*(AD+G*(Z(I,J)-Z(I,J-1))/DX+F*UU-(
     #            TSY-TBY)/HM)
            END DO
         END DO
C-----------------------------------------------------------------------
C  Treat Boundary Nodes as in Koutitas' wind circ model
C-----------------------------------------------------------------------
         DO 26 K = 1, KB
            I = IB(K)
            J = JB(K)
            GO TO (13,14,15,16,17,18,19,20,21,22,23) NB(K) - 1
   13       CONTINUE
            UN(I,J) = 0.0
            GO TO 25
   14       CONTINUE
            VN(I,J) = 0.0
            GO TO 25
   15       CONTINUE
            UN(I,J) = 0.0
            VN(I,J) = 0.0
            GO TO 25
   16       CONTINUE
            UN(I,J) = -Z(I,J)*SQRT(G/H(I,J))
            VN(I-1,J) = VN(I,J)
            GO TO 24
   17       CONTINUE
            UN(I+1,J) = Z(I,J)*SQRT(G/H(I,J))
            VN(I+1,J) = VN(I,J)
            GO TO 24
   18       CONTINUE
            VN(I,J+1) = Z(I,J)*SQRT(G/H(I,J))
            UN(I,J+1) = UN(I,J)
            GO TO 24
   19       CONTINUE
            VN(I,J) = -Z(I,J)*SQRT(G/H(I,J))
            UN(I,J-1) = UN(I,J)
            GO TO 24
   20       CONTINUE
            UN(I,J) = -Z(I,J)*SQRT(G/H(I,J))
            VN(I,J) = 0.0
            GO TO 24
   21       CONTINUE
            UN(I+1,J) = Z(I,J)*SQRT(G/H(I,J))
            VN(I,J) = 0.0
            GO TO 24
   22       CONTINUE
            UN(I,J) = 0.0
            VN(I,J+1) = Z(I,J)*SQRT(G/H(I,J))
            GO TO 24
   23       CONTINUE
            UN(I,J) = 0.0
            VN(I,J) = -Z(I,J)*SQRT(G/H(I,J))
   24       CONTINUE
C            ZR1 = Z(I,J) - ZI(I,J)
C           ZR2 = Z(I,J+1) - ZI(I,J+1)
C           ZI(I,J) = AMPL*DSIN(2.0*PI*(N-1)*DT/PER)
C           C = SQRT(G*H(I,J))
C           L = C*PER
C           ZI(I,J+1) = AMPL*DSIN(2.0*PI*(N-1)*DT/PER-DX/L)
C           ZR1 = ZR1 + DT/DX*C*(ZI(I,J+1)-ZI(I,J))
C           Z(I,J) = ZI(I,J) + ZR1
C           Flather BC (courtesy of Paul Meijer (Utrecht))
            un(i,j)= 0.0
            vn(i,j)= -sqrt(g/h(i,j))*( z(i,j) -
     &        ampl*dsin(2.0*pi*(n-1)*dt/per))
   25       CONTINUE
   26    CONTINUE
C-----------------------------------------------------------------------
C  Update u and v with their New Values
C-----------------------------------------------------------------------
         DO J = 1, JM
            DO I = 1, IM
               U(I,J) = UN(I,J)
               V(I,J) = VN(I,J)
            END DO
         END DO
C-----------------------------------------------------------------------
C  Translate u and v to the Mesh Centers and Store Along with z
C-----------------------------------------------------------------------
         DO J = 1, JM - 2
            DO I = 2, IM - 2
               UBAR(I,J) = (U(I,J)+U(I+1,J))/2.0
               VBAR(I,J) = (V(I,J)+V(I,J+1))/2.0
            END DO
         END DO
         TIME = N/DAT - INT(N/DAT)
         IF (TIME .EQ. 0.0) THEN
            KK = KK + 1
C-----------------------------------------------------------------------
C  Write the Results to a Regular File
C-----------------------------------------------------------------------
            P = 1
            Q = 7
            FORMT = '(25x,i6)'
            DO WHILE(P .LE. JM)
               Q = MIN0(JM,Q)
               WRITE (12, FORMT) N
               WRITE (13, FORMT) N
               DO I = 1, IM
                  WRITE (11, 100) (H(I,J), J = P, Q)
                  WRITE (12, 100) (UBAR(I,J), J = P, Q)
                  WRITE (13, 100) (VBAR(I,J), J = P, Q)
                  WRITE (14, 100) (Z(I,J), J = P, Q)
  100 format(4x,7(1x,1pe9.2))
               END DO
               P = Q + 1
               Q = Q + 7
            END DO
         ENDIF
   33 CONTINUE
      STOP 
      END
