PROGRAM U1_lattice
    IMPLICIT NONE

    CONTAINS


!--------JACKKNIFE SUBROUTINE---------------------------------------------------------
! This subroutine takes an array "data" of lengh N_data that contains all the measures
! from which we want to extract the mean value and its sigma, using the Jackknife method.
! These value are then stored inside the variable "average" and "error".

    SUBROUTINE jackknife(N_data, data, average, error)

        IMPLICIT NONE
        INTEGER:: i, N_data
        REAL(KIND = 8):: average, error, sig_reduced, sig, sig0, data_avg, sum_data
        REAL(KIND = 8), DIMENSION(N_data):: data, y

        sig0 = 0.0_8
        sum_data = 0.0_8
        DO i=1, N_data
            sum_data = sum_data + data(i)
        END DO
        data_avg = sum_data/REAL(N_data, 8)
        sig_reduced = 0.0_8
        DO i = 1, N_data
            y(i) = sum_data - data(i)
            y(i) = y(i)/REAL(N_data - 1, 8)
            sig_reduced = sig_reduced + (y(i) - data_avg)*(y(i) - data_avg)
        END DO
        sig=DSQRT((REAL(N_data-1,8)/REAL(N_data,8))*sig_reduced)
        IF (sig0 .lt. sig) THEN
            sig0 = sig
        END IF
        error = sig0
        average = data_avg

    END SUBROUTINE jackknife

!--------RAND1----------------------------------------------------------------------------------
! This is the Rand1() function taken from "Numerical Recipies in Fortran 90, Vol 2"

    FUNCTION rand1(idum)

        IMPLICIT NONE
        INTEGER, PARAMETER :: K4B = selected_int_kind(9)
        INTEGER(K4B), INTENT(INOUT) :: idum
        REAL :: rand1
        INTEGER(K4B), PARAMETER :: IA=16807, IM=2147483647, IQ=127773, IR=2836
        REAL, SAVE :: am
        INTEGER(K4B), SAVE :: ix = -1, iy = -1, k
        IF (idum <= 0 .or. iy < 0) THEN
            am = nearest(1.0,-1.0)/IM
            iy = ior(ieor(888889999,abs(idum)),1)
            ix = ieor(777755555,abs(idum))
            idum = abs(idum) + 1
        END IF
        ix = ieor(ix,ishft(ix,13)) 
        ix = ieor(ix,ishft(ix,-17))
        ix = ieor(ix,ishft(ix,5))
        k = iy/IQ 
        iy = IA * (iy-k*IQ) - IR * k
        IF (iy < 0) THEN 
            iy = iy + IM
        END IF
        rand1 = am * ior(iand(IM,ieor(ix,iy)),1)

    END FUNCTION rand1

END PROGRAM U1_lattice