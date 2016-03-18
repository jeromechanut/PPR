
    !
    ! POLYFN.f90: polynomial basis vectors, reconstruction
    ! coefficients and utilities.
    !
    ! Darren Engwirda 
    ! 17-Mar-2016
    ! engwirda [at] mit [dot] edu
    !

    module polyfn

    implicit none

    integer, parameter :: max_poly_ndof =   +7

    public :: poly_basis
    public :: poly_value

    contains

    !----------------------------------------------------------------!
    ! POLY_BASIS: evaluate the basis vector BVEC (SVAL).             !
    !----------------------------------------------------------------!
    !   NDOF - no. degrees-of-freedom.                               !
    !   SVAL - local coordinate, -1.d0 <= SVAL <= +1.d0.             !
    !   BVEC - NDOF-by-1 basis vector.                               !
    !----------------------------------------------------------------!
    subroutine poly_basis(nord,ndof,sval,bvec)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: nord,ndof
        real*8 , intent( in) :: sval
        real*8 , intent(out) :: bvec(:)

    !------------------------------------ push basis order !
    !   -1 : integral                                      !
    !   +0 : function                                      !
    !   +1 : 1st derivative                                !
    !   +2 : 2nd derivative                                !
    !------------------------------------------------------!

        select case (nord)
        case (-1)
    !------------------------------------ -1th-order basis !
            select case (ndof)
            case (+1)
                bvec(1) = sval**1 / 1.d0
            case (+2)
                bvec(1) = sval**1 / 1.d0
                bvec(2) = sval**2 / 2.d0
            case (+3)
                bvec(1) = sval**1 / 1.d0
                bvec(2) = sval**2 / 2.d0
                bvec(3) = sval**3 / 3.d0
            case (+4)
                bvec(1) = sval**1 / 1.d0
                bvec(2) = sval**2 / 2.d0
                bvec(3) = sval**3 / 3.d0
                bvec(4) = sval**4 / 4.d0
            case (+5)
                bvec(1) = sval**1 / 1.d0
                bvec(2) = sval**2 / 2.d0
                bvec(3) = sval**3 / 3.d0
                bvec(4) = sval**4 / 4.d0
                bvec(5) = sval**5 / 5.d0
            case (+6)
                bvec(1) = sval**1 / 1.d0
                bvec(2) = sval**2 / 2.d0
                bvec(3) = sval**3 / 3.d0
                bvec(4) = sval**4 / 4.d0
                bvec(5) = sval**5 / 5.d0
                bvec(6) = sval**6 / 6.d0
            case (+7)
                bvec(1) = sval**1 / 1.d0
                bvec(2) = sval**2 / 2.d0
                bvec(3) = sval**3 / 3.d0
                bvec(4) = sval**4 / 4.d0
                bvec(5) = sval**5 / 5.d0
                bvec(6) = sval**6 / 6.d0
                bvec(7) = sval**7 / 7.d0
            case default
                bvec =    0.d0
            end select

        case (+0)
    !------------------------------------ +0th-order basis !
            select case (ndof)
            case (+1)
                bvec(1) = 1.d0
            case (+2)
                bvec(1) = 1.d0
                bvec(2) = sval**1
            case (+3)
                bvec(1) = 1.d0
                bvec(2) = sval**1
                bvec(3) = sval**2
            case (+4)
                bvec(1) = 1.d0
                bvec(2) = sval**1
                bvec(3) = sval**2
                bvec(4) = sval**3
            case (+5)
                bvec(1) = 1.d0
                bvec(2) = sval**1
                bvec(3) = sval**2
                bvec(4) = sval**3
                bvec(5) = sval**4
            case (+6)
                bvec(1) = 1.d0
                bvec(2) = sval**1
                bvec(3) = sval**2
                bvec(4) = sval**3
                bvec(5) = sval**4
                bvec(6) = sval**5
            case (+7)
                bvec(1) = 1.d0
                bvec(2) = sval**1
                bvec(3) = sval**2
                bvec(4) = sval**3
                bvec(5) = sval**4
                bvec(6) = sval**5
                bvec(7) = sval**6
            case default
                bvec =    0.d0
            end select

        case (+1)
    !------------------------------------ +1st-order basis !
            select case (ndof)
            case (+1)
                bvec(1) =           0.d0
            case (+2)
                bvec(1) =           0.d0
                bvec(2) =           1.d0
            case (+3)
                bvec(1) =           0.d0
                bvec(2) =           1.d0
                bvec(3) = sval**1 * 2.d0
            case (+4)
                bvec(1) =           0.d0
                bvec(2) =           1.d0
                bvec(3) = sval**1 * 2.d0
                bvec(4) = sval**2 * 3.d0
            case (+5)
                bvec(1) =           0.d0
                bvec(2) =           1.d0
                bvec(3) = sval**1 * 2.d0
                bvec(4) = sval**2 * 3.d0
                bvec(5) = sval**3 * 4.d0
            case (+6)
                bvec(1) =           0.d0
                bvec(2) =           1.d0
                bvec(3) = sval**1 * 2.d0
                bvec(4) = sval**2 * 3.d0
                bvec(5) = sval**3 * 4.d0
                bvec(6) = sval**4 * 5.d0
            case (+7)
                bvec(1) =           0.d0
                bvec(2) =           1.d0
                bvec(3) = sval**1 * 2.d0
                bvec(4) = sval**2 * 3.d0
                bvec(5) = sval**3 * 4.d0
                bvec(6) = sval**4 * 5.d0
                bvec(7) = sval**5 * 6.d0
            case default
                bvec =              0.d0
            end select

        case (+2)
    !------------------------------------ +2nd-order basis !
            select case (ndof)
            case (+1)
                bvec(1) =           0.d0
            case (+2)
                bvec(1) =           0.d0
                bvec(2) =           0.d0
            case (+3)
                bvec(1) =           0.d0
                bvec(2) =           0.d0
                bvec(3) =           2.d0
            case (+4)
                bvec(1) =           0.d0
                bvec(2) =           0.d0
                bvec(3) =           2.d0
                bvec(4) = sval**1 * 6.d0
            case (+5)
                bvec(1) =           0.d0
                bvec(2) =           0.d0
                bvec(3) =           2.d0
                bvec(4) = sval**1 * 6.d0
                bvec(5) = sval**2 *12.d0
            case (+6)
                bvec(1) =           0.d0
                bvec(2) =           0.d0
                bvec(3) =           2.d0
                bvec(4) = sval**1 * 6.d0
                bvec(5) = sval**2 *12.d0
                bvec(6) = sval**3 *20.d0
            case (+7)
                bvec(1) =           0.d0
                bvec(2) =           0.d0
                bvec(3) =           2.d0
                bvec(4) = sval**1 * 6.d0
                bvec(5) = sval**2 *12.d0
                bvec(6) = sval**3 *20.d0
                bvec(7) = sval**4 *30.d0
            case default
                bvec =              0.d0
            end select

        case default
            bvec = 0.d0
        end select

    end subroutine poly_basis

    !----------------------------------------------------------------!
    ! POLY_VALUE: evaluate FUNC(SVAL) = DOT(FHAT,BVEC(SVAL)).        !
    !----------------------------------------------------------------!
    !   NDOF - no. degrees-of-freedom.                               !
    !   FHAT - NDOF-by-1 array of coefficients.                      !
    !   SVAL - local coordinate, -1.d0 <= SVAL <= +1.d0.             !
    !   BVEC - NDOF-by-1 basis vector.                               !
    !----------------------------------------------------------------!
    function poly_value(nord,ndof,fhat,sval, &
        &               bvec) result  (rval)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: nord,ndof
        real*8 , intent( in) :: fhat(:)
        real*8 , intent( in) :: sval
        real*8 , intent( in), optional :: bvec(:)

    !------------------------------------------- variables !
        real*8  :: rval
        integer :: idof
        real*8  :: lvec(max_poly_ndof)

        rval = 0.d0

        if (present(bvec)) then

    !--------------------------- use existing basis vector !

            do  idof = +1,ndof
            
                rval = & 
        &       rval + fhat(idof) * bvec(idof)
            
            end do

        else

    !--------------------------- assemble new basis vector !

            call poly_basis(nord,ndof, &
            &               sval,lvec)

            do  idof = +1,ndof
            
                rval = &
        &       rval + fhat(idof) * lvec(idof)
            
            end do

        end if

    end function poly_value

    end module polyfn



