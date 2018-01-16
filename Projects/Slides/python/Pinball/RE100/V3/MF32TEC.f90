       	program MF32TEC

        implicit none

        REAL, DIMENSION(:)  , ALLOCATABLE :: X,Y,Vx,Vy,P,VORT
        REAL :: U1,U2,U3,U4,V1,V2,V3,V4,X1,X2,X3,X4,Y1,Y2,Y3,Y4
        REAL :: DX,DY,AREA,CIRC,Vort_max,Vort_min
        INTEGER, DIMENSION(:,:)  , ALLOCATABLE :: IX
        INTEGER, DIMENSION(:,:)  , ALLOCATABLE :: ICON
        INTEGER :: NDE,NEL,i,j,K1,K2,K3,K4

!———————————————
       open(unit=4,file='../Code_Input/GRID6',err=1000,status='OLD')
!———————————————
       open(unit=9,file='../Code_Output/Flow.dat',err=1000,status='OLD')

       open(76,file='TEC.plt',form='FORMATTED',status='UNKNOWN')


       read(4,*)
       read(4,*) NDE

       do i=1,NDE
        read(4,*)
       enddo
       read(4,*)
       read(4,*) NEL
       Rewind(4)

       ALLOCATE (X(NDE),Y(NDE),Vx(NDE),Vy(NDE),P(NDE),VORT(NDE))
       ALLOCATE (IX(6,NDE))
       ALLOCATE (ICON(4*NEL,3))
       



       write(76,*) 'TITLE = "TECPLOT: FE-Data"'
       write(76,*) 'VARIABLES = "X", "Y", "Vx", "Vy", "P", "VORT"' 
       write(76,*) 'ZONE NODES=',NDE ,' ELEMENTS=',4*NEL,  &
                   ' DATAPACKING=POINT, ZONETYPE=FETRIANGLE'


       read(4,*)
       read(4,*) NDE

       do i=1,NDE
        read(4,*) X(i),Y(i)
        read(9,*) Vx(i), Vy(i), P(i)
       enddo

       read(4,*)
       read(4,*) NEL

       do i=1,NEL
        read(4,*) (IX(j,i),j=1,6)
       enddo


!————————————————————
!—— Calculate Vorticity
!————————————————————



      do i=1,NEL
       ICON(4*I-3,1:3) = (/abs(IX(1,I)) , abs(IX(6,I)) , abs(IX(4,I)) /)
       ICON(4*I-2,1:3) = (/abs(IX(4,I)) , abs(IX(6,I)) , abs(IX(5,I)) /)
       ICON(4*I-1,1:3) = (/abs(IX(5,I)) , abs(IX(6,I)) , abs(IX(3,I)) /)
       ICON(4*I-0,1:3) = (/abs(IX(4,I)) , abs(IX(5,I)) , abs(IX(2,I)) /)
      enddo


 
       DO i = 1, NDE
         VORT(i)  = 0.0
       ENDDO

      

       DX = 0.0
       DY = 0.0
       DO i = 1, NEL*4

         K1= ICON(i,1)
         K2= ICON(i,2)
         K3= ICON(i,3)
         K4= ICON(i,1)

         X1 = X(K1)
         X2 = X(K2)
         X3 = X(K3)
         X4 = X(K4)
         DX = MAX(DX, MAX(X1,X2,X3,X4)-MIN(X1,X2,X3,X4))
         Y1 = Y(K1)
         Y2 = Y(K2)
         Y3 = Y(K3)
         Y4 = Y(K4)
         DY = MAX(DY, MAX(Y1,Y2,Y3,Y4)-MIN(Y1,Y2,Y3,Y4))
         U1 = Vx(K1)
         V1 = Vy(K1)
         U2 = Vx(K2)
         V2 = Vy(K2)
         U3 = Vx(K3)
         V3 = Vy(K3)
         U4 = Vx(K4)
         V4 = Vy(K4)

         AREA = ( (X3-X1)*(Y4-Y2) + (X2-X4)*(Y3-Y1) )
         CIRC = ( (U1+U2)*(X2-X1) + (V1+V2)*(Y2-Y1)   &
                + (U2+U3)*(X3-X2) + (V2+V3)*(Y3-Y2)   &
                + (U3+U4)*(X4-X3) + (V3+V4)*(Y4-Y3)   &
                + (U4+U1)*(X1-X4) + (V4+V1)*(Y1-Y4) ) / AREA

         VORT(K1) = VORT(K1) + CIRC
         VORT(K2) = VORT(K2) + CIRC
         VORT(K3) = VORT(K3) + CIRC
         IF(K4.NE.K1) VORT(K4) = VORT(K4) + CIRC

        if (Vort_max.lt.VORT(k1)) Vort_max=VORT(k1)
        if (Vort_max.lt.VORT(k2)) Vort_max=VORT(k2)
        if (Vort_max.lt.VORT(k3)) Vort_max=VORT(k3)
        if (Vort_min.gt.VORT(k1)) Vort_min=VORT(k1)
        if (Vort_min.gt.VORT(k2)) Vort_min=VORT(k2)
        if (Vort_min.gt.VORT(k3)) Vort_min=VORT(k3)
      ENDDO

       do i=1,NDE
        write(76,'(6E12.4)') X(i),Y(i),Vx(i),Vy(i),P(i),VORT(i)
       enddo


       do i=1,NEL
      WRITE(76,*) abs(IX(1,I)),abs(IX(6,I)),abs(IX(4,I))
      WRITE(76,*) abs(IX(4,I)),abs(IX(6,I)),abs(IX(5,I))
      WRITE(76,*) abs(IX(5,I)),abs(IX(6,I)),abs(IX(3,I))
      WRITE(76,*) abs(IX(4,I)),abs(IX(5,I)),abs(IX(2,I))
       enddo
       write(76,*) 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXX'

       stop

 1000  write(*,*)'Error in read'

        stop
        end


