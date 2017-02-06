SUBROUTINE matrixInversion(N,M,AA,BB,DET)                                
  PARAMETER(EPSMACH=1.D-20)                                    
  REAL*8 AA(N,N),BB(N,M)                                                   
  REAL*8,POINTER ::  PC(:),PL(:),CS(:)                                                  
  REAL*8 PV,PAV,DET,TT    
  
  !Initializations :                       
  ALLOCATE(PC(1:N),stat=ialloc)
  ALLOCATE(PL(1:N),stat=ialloc)
  ALLOCATE(CS(1:N),stat=ialloc)                                         
  DET=1.D0

  PC(:) = 0.
  PL(:) = 0.
  CS(:) = 0.

  !main loop                                                                   
kloop:  DO K=1,N                                                                  
   !Searching greatest pivot :                                               
     PV=AA(K,K)                                                              
     IK=K                                                                    
     JK=K                                                                    
     PAV=ABS(PV)                                                            
     DO j=K,N                                                                
        DO i=K,N                                                              
           IF (ABS(AA(I,J)) > PAV) THEN                                      
              PV=AA(I,J)                                                        
              PAV=ABS(PV)                                                      
              IK=I                                                              
              JK=J                                                              
           ENDIF
        ENDDO
     ENDDO
     
     !Search terminated, the pivot is in location I=IK, J=JK.
     !Memorizing pivot location: :                                        
     
     PC(K)=JK                                                                
     PL(K)=IK                                                                
     
     !Determinant DET is actualised
     !If DET=0, ERROR MESSAGE and STOP
     !Machine dependent EPSMACH equals here 1.DE-20                                        
     
     IF (IK.NE.K) DET=-DET                                                   
     IF (JK.NE.K) DET=-DET                                                   
     DET=DET*PV                                                              
     IF (ABS(DET) < EPSMACH) THEN                                          
        !Error message and Stop                                             
        PRINT 10                                                              
        STOP                                                                  
     ENDIF
     
     !POSITIONING PIVOT IN K,K:                                              
     IF(IK /= K) THEN                                                        
        DO I=1,N                                                              
           !EXCHANGE LINES IK and K:                                       
           TT=AA(IK,I) 
           AA(IK,I)=AA(K,I)
           AA(K,I)=TT                                                     
        ENDDO
     ENDIF
     
     !Pivot is at correct line                                                

     IF(JK /= K) THEN                                                        
        DO I=1,N                                                              
           !Exchange columns JK and K of matrix AA 
           TT=AA(I,JK)                                                         
           AA(I,JK)=AA(I,K)                                                    
           AA(I,K)=TT                                                          
        ENDDO
     ENDIF
     !Pivot is at correct column and located in K,K 
     
     !Store column K in vector CS                             
     !then set column K to zero                                             

     DO I=1,N                                                                
        CS(I)=AA(I,K)                                                         
        AA(I,K)=0.D0                                                          
     ENDDO
                                                                                    
     CS(K)=0.                                                                
     AA(K,K)=1.                                                              
     
     !Modify line K :                                            
     IF(ABS(PV) < EPSMACH) THEN                                            
        WRITE(*,*) '  PIVOT TOO SMALL - STOP'                               
        STOP                                                                  
     ENDIF

     DO I=1,N                                                                
        AA(K,I)=AA(K,I)/PV                                                    
     ENDDO

     !Modify other lines of matrix AA:                                        
     DO J=1,N                                                                
        IF (J == K) CONTINUE                                                  
        DO I=1,N                                                              
           !Modify line J of matrix AA :                                            
           AA(J,I)=AA(J,I)-CS(J)*AA(K,I)                                       
        ENDDO
     END DO

    !Line K is ready.                                                
    
      END DO kloop                                                                     

      !End of K loop                                                              
      
      !The matrix AA is inverted - Rearrange AA                         
      
      !Exchange lines                                                            
      DO I=N,1,-1                                                               
         IK=PC(I)                                                                
         IF (IK.EQ.I) CONTINUE                                                   
         !EXCHANGE LINES I AND PC(I) OF AA:                                         
         DO J=1,N                                                                
            TT=AA(I,J)                                                            
            AA(I,J)=AA(IK,J)                                                      
            AA(IK,J)=TT                                                           
         ENDDO
         
         !NO MORE EXCHANGE NEEDED                                                      
          !GO TO NEXT LINE                                                  
      ENDDO
      
      !EXCHANGE COLUMNS                                                          
      
      DO J=N,1,-1                                                               
         JK=PL(J)                                                                
         IF (JK.EQ.J) CONTINUE                                                   
         !EXCHANGE COLUMNS J AND PL(J) OF AA :                                       
         DO I=1,N                                                                
            TT=AA(I,J)                                                            
            AA(I,J)=AA(I,JK)                                                      
            AA(I,JK)=TT                                                           
         ENDDO
          !NO MORE EXCHANGE NEEDED                                                      
          !GO TO NEXT COLUMN   
      ENDDO
      !REARRANGEMENT TERMINATED.                                                        

det = 1.0D0

      RETURN                                                                    
10    FORMAT(///'  DETERMINANT EQUALS ZERO, NO SOLUTION!')                    
      
END SUBROUTINE matrixInversion
