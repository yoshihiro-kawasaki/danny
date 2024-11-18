 program main
 EXTERNAL FEX, JEX
 DOUBLE PRECISION ATOL, RTOL, RWORK, T, TOUT, Y
 DIMENSION Y(12), RWORK(500), IWORK(30)
     DATA LRW/500/, LIW/30/
     NEQ = 12
     DO 10 I = 1,NEQ
 10    Y(I) = 0.0D0
     Y(1) = 1.0D0
     T = 0.0D0
     TOUT = 0.1D0
     ITOL = 1
     RTOL = 1.0D-4
     ATOL = 1.0D-6
     ITASK = 1
     ISTATE = 1
     IOPT = 0
     MF = 121
     DO 40 IOUT = 1,5
       CALL DLSODES (FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL,
    1     ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
       WRITE(6,30)T,IWORK(11),RWORK(11),(Y(I),I=1,NEQ)
 30    FORMAT(//' At t =',D11.3,4X,
    1    ' No. steps =',I5,4X,' Last step =',D11.3/
    2    '  Y array =  ',4D14.5/13X,4D14.5/13X,4D14.5)
       IF (ISTATE .LT. 0) GO TO 80
       TOUT = TOUT*10.0D0
 40    CONTINUE
     LENRW = IWORK(17)
     LENIW = IWORK(18)
     NST = IWORK(11)
     NFE = IWORK(12)
     NJE = IWORK(13)
     NLU = IWORK(21)
     NNZ = IWORK(19)
     NNZLU = IWORK(25) + IWORK(26) + NEQ
     WRITE (6,70) LENRW,LENIW,NST,NFE,NJE,NLU,NNZ,NNZLU
 70  FORMAT(//' Required RWORK size =',I4,'   IWORK size =',I4/
    1   ' No. steps =',I4,'   No. f-s =',I4,'   No. J-s =',I4,
    2   '   No. LU-s =',I4/' No. of nonzeros in J =',I5,
    3   '   No. of nonzeros in LU =',I5)
     STOP
 80  WRITE(6,90)ISTATE
 90  FORMAT(///' Error halt.. ISTATE =',I3)
     STOP

end program main

    SUBROUTINE FEX (NEQ, T, Y, YDOT)
        DOUBLE PRECISION T, Y, YDOT
        DOUBLE PRECISION RK1, RK2, RK3, RK4, RK5, RK6, RK7, RK8, RK9,
    1   RK10, RK11, RK12, RK13, RK14, RK15, RK16, RK17
        DIMENSION Y(12), YDOT(12)
        DATA RK1/0.1D0/, RK2/10.0D0/, RK3/50.0D0/, RK4/2.5D0/, RK5/0.1D0/,
    1   RK6/10.0D0/, RK7/50.0D0/, RK8/2.5D0/, RK9/50.0D0/, RK10/5.0D0/,
    2   RK11/50.0D0/, RK12/50.0D0/, RK13/50.0D0/, RK14/30.0D0/,
    3   RK15/100.0D0/, RK16/2.5D0/, RK17/100.0D0/, RK18/2.5D0/,
    4   RK19/50.0D0/, RK20/50.0D0/
        YDOT(1)  = -RK1*Y(1)
        YDOT(2)  = RK1*Y(1) + RK11*RK14*Y(4) + RK19*RK14*Y(5)
    1           - RK3*Y(2)*Y(3) - RK15*Y(2)*Y(12) - RK2*Y(2)
        YDOT(3)  = RK2*Y(2) - RK5*Y(3) - RK3*Y(2)*Y(3) - RK7*Y(10)*Y(3)
    1           + RK11*RK14*Y(4) + RK12*RK14*Y(6)
        YDOT(4)  = RK3*Y(2)*Y(3) - RK11*RK14*Y(4) - RK4*Y(4)
        YDOT(5)  = RK15*Y(2)*Y(12) - RK19*RK14*Y(5) - RK16*Y(5)
        YDOT(6)  = RK7*Y(10)*Y(3) - RK12*RK14*Y(6) - RK8*Y(6)
        YDOT(7)  = RK17*Y(10)*Y(12) - RK20*RK14*Y(7) - RK18*Y(7)
        YDOT(8)  = RK9*Y(10) - RK13*RK14*Y(8) - RK10*Y(8)
        YDOT(9)  = RK4*Y(4) + RK16*Y(5) + RK8*Y(6) + RK18*Y(7)
        YDOT(10) = RK5*Y(3) + RK12*RK14*Y(6) + RK20*RK14*Y(7)
    1           + RK13*RK14*Y(8) - RK7*Y(10)*Y(3) - RK17*Y(10)*Y(12)
    2           - RK6*Y(10) - RK9*Y(10)
        YDOT(11) = RK10*Y(8)
        YDOT(12) = RK6*Y(10) + RK19*RK14*Y(5) + RK20*RK14*Y(7)
    1           - RK15*Y(2)*Y(12) - RK17*Y(10)*Y(12)
    RETURN
    END

     SUBROUTINE JEX (NEQ, T, Y, J, IA, JA, PDJ)
     DOUBLE PRECISION T, Y, PDJ
     DOUBLE PRECISION RK1, RK2, RK3, RK4, RK5, RK6, RK7, RK8, RK9,
    1   RK10, RK11, RK12, RK13, RK14, RK15, RK16, RK17
     DIMENSION Y(12), IA(*), JA(*), PDJ(12)
     DATA RK1/0.1D0/, RK2/10.0D0/, RK3/50.0D0/, RK4/2.5D0/, RK5/0.1D0/,
    1   RK6/10.0D0/, RK7/50.0D0/, RK8/2.5D0/, RK9/50.0D0/, RK10/5.0D0/,
    2   RK11/50.0D0/, RK12/50.0D0/, RK13/50.0D0/, RK14/30.0D0/,
    3   RK15/100.0D0/, RK16/2.5D0/, RK17/100.0D0/, RK18/2.5D0/,
    4   RK19/50.0D0/, RK20/50.0D0/
     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), J
 1   PDJ(1) = -RK1
     PDJ(2) = RK1
     RETURN
 2   PDJ(2) = -RK3*Y(3) - RK15*Y(12) - RK2
     PDJ(3) = RK2 - RK3*Y(3)
     PDJ(4) = RK3*Y(3)
     PDJ(5) = RK15*Y(12)
     PDJ(12) = -RK15*Y(12)
     RETURN
 3   PDJ(2) = -RK3*Y(2)
     PDJ(3) = -RK5 - RK3*Y(2) - RK7*Y(10)
     PDJ(4) = RK3*Y(2)
     PDJ(6) = RK7*Y(10)
     PDJ(10) = RK5 - RK7*Y(10)
     RETURN
 4   PDJ(2) = RK11*RK14
     PDJ(3) = RK11*RK14
     PDJ(4) = -RK11*RK14 - RK4
     PDJ(9) = RK4
     RETURN
 5   PDJ(2) = RK19*RK14
     PDJ(5) = -RK19*RK14 - RK16
     PDJ(9) = RK16
     PDJ(12) = RK19*RK14
     RETURN
 6   PDJ(3) = RK12*RK14
     PDJ(6) = -RK12*RK14 - RK8
     PDJ(9) = RK8
     PDJ(10) = RK12*RK14
     RETURN
 7   PDJ(7) = -RK20*RK14 - RK18
     PDJ(9) = RK18
     PDJ(10) = RK20*RK14
     PDJ(12) = RK20*RK14
     RETURN
 8   PDJ(8) = -RK13*RK14 - RK10
     PDJ(10) = RK13*RK14
     PDJ(11) = RK10
 9   RETURN
 10  PDJ(3) = -RK7*Y(3)
     PDJ(6) = RK7*Y(3)
     PDJ(7) = RK17*Y(12)
     PDJ(8) = RK9
     PDJ(10) = -RK7*Y(3) - RK17*Y(12) - RK6 - RK9
     PDJ(12) = RK6 - RK17*Y(12)
 11  RETURN
 12  PDJ(2) = -RK15*Y(2)
     PDJ(5) = RK15*Y(2)
     PDJ(7) = RK17*Y(10)
     PDJ(10) = -RK17*Y(10)
     PDJ(12) = -RK15*Y(2) - RK17*Y(10)
     RETURN
     END