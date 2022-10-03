
	PROGRAM SHOCK_WAVE


      REAL (KIND=8) PR_1,PR_2,RO_1,RO_2,T_1,T_2,A_1,MACH,TETA,BETA,GAMMA
      REAL (KIND=8) PI,R,V_1,ST,SB,CT,CB,DTDT,FT,ET,TET,MT,HT,Vn,Un


      GAMMA = 1.4

      PI = 3.1416

      R = 8.3145D+03/29.

      MACH  = 5.0
      BETA  = 15.0
      TETA = PI*BETA/360. + PI/4.
      BETA = PI*BETA/180.
      
      PR_1 = 1.01325D+05
       T_1 = 288.0

      RO_1 = PR_1/(R*T_1) 
      A_1  = DSQRT (GAMMA*PR_1/RO_1)
      V_1  = MACH*A_1

      SB = DSIN (BETA)
      CB = DCOS (BETA)

      MT = (1./MACH)**2

      HT = 0.5*(GAMMA+1.)

      ET = 1.

      DO WHILE (ET > 1.0D-02)

      ST = DSIN (TETA) 
      CT = DCOS (TETA)

      DTDT =   2.*ST*CT - HT*CT*SB/ DCOS (TETA-BETA) - DSIN (TETA-BETA)*HT*ST*SB/ DCOS (TETA-BETA)**2

      FT = ST**2-HT*ST*SB/ DCOS (TETA-BETA) - MT

      TET = TETA - FT/DTDT

      ET = DABS ((TET-TETA)/TETA)

      TETA = TET

      END DO

      Vn = V_1*DSIN(TETA)
      Un = 1.-(DSIN(TETA-BETA)*DCOS(TETA))/(DSIN(TETA)*DCOS(TETA-BETA))

      RO_2 = RO_1/(1.-Un)
      PR_2 = PR_1+RO_1*Un*Vn**2
       T_2 = PR_2/(RO_2*R)

      TETA = 180.*TETA/PI
      BETA = 180.*BETA/PI
      
      write(1,"(1x,G11.2)")MACH
      write(1,"(1x,G11.2)")BETA
      write(1,"(1x,G11.2)")TETA
      write(1,"(1x,G11.2)")RO_2
      write(1,"(1x,G11.2)")PR_2
      write(1,"(1x,G11.2)")T_2

      
     
	END PROGRAM SHOCK_WAVE