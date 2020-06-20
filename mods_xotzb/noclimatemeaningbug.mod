*DECLARE DUMPCTL1
*B GKR1F404.70
            IF (STP2im(A_IM).gt.0)THEN
*I GLW2F402.12
            ELSE
              IF (END_DUMPim(A_IM).NE."              " .AND.
     &           .NOT. LKEEPATM ) THEN
!               Filename to be deleted is not blank and is not to 
!               be kept until the ocean dump for the current step
!               is written
                WRITE (8,610) END_DUMPim(A_IM) ! Delete request
                CLOSE(8)
                OPEN(8,FILE=FILENAME)
              ENDIF
            ENDIF
*B GKR1F404.142
            IF (STP2im(O_IM).gt.0)THEN
*I GKR1F404.143
              IF (END_DUMPim(O_IM).NE."              ") THEN
!               Filename to be deleted is not blank
                WRITE(8,610) END_DUMPim(O_IM)
                CLOSE(8)
                OPEN(8,FILE=FILENAME)
              ENDIF

              IF (LDELATM .AND. LASTATMim(A_IM).NE."              ")THEN
!               There is an atmos dump to delete and the filename
!               to be deleted is not blank
                  WRITE(8,610) LASTATMim(A_IM)     
                  CLOSE(8)             
              ENDIF 
              OPEN(8,FILE=FILENAME)
            ENDIF
            ELSE
