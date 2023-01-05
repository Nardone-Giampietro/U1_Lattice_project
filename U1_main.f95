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
        IF (sig0.lt.sig) THEN
            sig0=sig
        END IF
        error=sig0
        average=data_avg

    END SUBROUTINE jackknife

END PROGRAM U1_lattice