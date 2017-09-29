! -------------------------------------------------------------------------------------------------
! This function determine the lines of files
! -------------------------------------------------------------------------------------------------
        FUNCTION count_lines(filename) result(nlines)

        IMPLICIT NONE
        CHARACTER(LEN=*)    :: filename
        INTEGER             :: nlines 
        INTEGER             :: io

        OPEN(10,FILE=filename, IOSTAT=io, STATUS='old')
        IF (io/=0) STOP 'Cannot open file! '

        nlines = 0
        DO
           READ(10,*,IOSTAT=io)
           IF (io/=0) exit
           nlines = nlines + 1
        END DO
        CLOSE(10)

        END FUNCTION count_lines
