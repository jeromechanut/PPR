
    !
    ! MATRIX.f90: matrix factorisation, adapted from LAPACK.
    !
    ! Darren Engwirda 
    ! 17-Mar-2016
    ! engwirda [at] mit [dot] edu
    !

    module matrix

    implicit none

    public :: trisolve
    public :: lufactor
    public :: lusolver

    contains

    !----------------------------------------------------------------!
    ! ROWSWAPS: apply a set of row-permutations.                     !
    !----------------------------------------------------------------!
    subroutine rowswaps(ncol,aa,adim,k1,k2,ipiv,incx)

        implicit none

        integer, intent(in)    :: ncol
        integer, intent(in)    :: adim,incx
        integer, intent(in)    :: k1,k2
        integer, intent(in)    :: ipiv(*)
        real*8 , intent(inout) :: aa(adim,*)

        integer :: i0,i1,i2,ix,ip,ii,kk,inc
        real*8  :: ta

    !----------------------------------------- quick exits !
        if (ncol.eq.0) return
        if (incx.eq.0) return

    !----------------------------------------- set indexes !
        if (incx.gt.0) then
            i0 = k1
            i1 = k1
            i2 = k2
            inc = 1
        else
            i0 = 1 + (1 - k2) * incx
            i1 = k2
            i2 = k1
            inc =-1
        end if

    !---- swap AA(II,:) and AA(IPIV(II),:) for IPIV(K1:K2) !
        ix = i0
        do  ii = i1, i2, inc
            ip = ipiv(ix)
            if (ip.ne.ii) then
            do  kk = 1, ncol
                ta        = aa(ii,kk)
                aa(ii,kk) = aa(ip,kk)
                aa(ip,kk) = ta
            end do
            end if
            ix = ix + incx
        end do
        
        return

    end subroutine rowswaps

    !----------------------------------------------------------------!
    ! TRISOLVE: solve a set of (triangular) linear systems.          !
    !----------------------------------------------------------------!
    subroutine trisolve(uplo,tran,diag,nrow,nrhs,aa,adim,xx,xdim)

        implicit none

    !------------------------------------------- arguments !
        character , intent(in) :: uplo,tran,diag
        integer, intent(in)    :: nrow,nrhs
        integer, intent(in)    :: adim,xdim
        real*8 , intent(in)    :: aa(adim,*)
        real*8 , intent(inout) :: xx(xdim,*)

    !------------------------------------------- variables !
        integer :: ii,jj,kk
        logical :: unit_diag

        if (nrow.le.0) return

        if (diag.eq.'N') then
            unit_diag =.FALSE.
        else
            unit_diag = .TRUE.
        end if

        if (tran.eq.'N') then
        
        if (uplo.eq.'U') then
    !------------------------------------- solve U * X = X !

            do  jj = nrow, 1, -1
    !-------------------------------------------- diagonal !
                if (.not.unit_diag) then
                do  kk = 1, nrhs, +1
                    xx(jj,kk) = xx(jj,kk) / aa(jj,jj)
                end do
                end if
    !-------------------------------------------- backsubs !
                do  kk = 1, nrhs, +1
                do  ii = jj-1, 1, -1
                    xx(ii,kk) = xx(ii,kk) - xx(jj,kk) &
            &                             * aa(ii,jj)
                end do
                end do

            end do

        else
    !------------------------------------- solve L * X = X !

            do  jj = 1, nrow
    !-------------------------------------------- diagonal !
                if (.not.unit_diag) then
                do  kk = 1, nrhs, +1
                    xx(jj,kk) = xx(jj,kk) / aa(jj,jj)
                end do
                end if
    !-------------------------------------------- backsubs !
                do  kk =    1, nrhs, +1
                do  ii = jj+1, nrow, +1
                    xx(ii,kk) = xx(ii,kk) - xx(jj,kk) &
            &                             * aa(ii,jj)
                end do
                end do
            end do

        end if

        else

        if (uplo.eq.'U') then
    !----------------------------------- solve U^T * X = X !

            do  jj = 1, nrow
    !-------------------------------------------- backsubs !
                do  kk = 1, nrhs, +1
                do  ii = 1, jj+0, +1
                    xx(jj,kk) = xx(jj,kk) - xx(ii,kk) &
            &                             * aa(ii,jj)
                end do
                end do
    !-------------------------------------------- diagonal !
                if (.not.unit_diag) then
                do  kk = 1, nrhs, +1
                    xx(jj,kk) = xx(jj,kk) / aa(jj,jj)
                end do
                end if
            end do

        else
    !----------------------------------- solve L^T * X = X !

            do  jj = nrow, 1, -1
    !-------------------------------------------- backsubs !
                do  kk =    1, nrhs, +1
                do  ii = nrow, jj+1, -1
                    xx(jj,kk) = xx(jj,kk) - xx(ii,kk) &
            &                             * aa(ii,jj)
                end do
                end do
    !-------------------------------------------- diagonal !
                if (.not.unit_diag) then
                do  kk = 1, nrhs, +1
                    xx(jj,kk) = xx(jj,kk) / aa(jj,jj)
                end do
                end if
            end do

        end if
        
        end if
        
        return

    end subroutine trisolve

    !----------------------------------------------------------------!
    ! LUFACTOR: form the LU factorisation of a matrix.               !
    !----------------------------------------------------------------!
    !   AA   - NROW-by-NCOL array, storing the coefficient matrix on !
    !          entry and the LU factorisation on exit.               !
    !          AA is allocated with leading dimension ADIM, such th- !
    !          at AA = AA(1:ADIM,:).                                 !
    !   IPIV - NROW-by-1 array of row permutations due to pivoting.  !
    !   INFO - Error flag. ZERO if factorisation successful, negati- !
    !          ve if inputs sized incorrectly.                       !
    !----------------------------------------------------------------!
    subroutine lufactor(nrow,ncol,aa,adim,ipiv,info)

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)    :: nrow,ncol
        integer, intent(in)    :: adim
        integer, intent(inout) :: ipiv(*)
        real*8 , intent(inout) :: aa(adim,*)
        integer, intent(out)   :: info

    !------------------------------------------- variables !
        integer :: ii,jj,kk,jp
        real*8  :: ta,da

    !----------------------------------------- error tests !
        info = 0
        if (nrow.lt.0) then
            info = -1
        end if
        if (ncol.lt.0) then
            info = -2
        end if
        if (adim.lt.max(1, nrow)) then
            info = -4
        end if

    !----------------------------------------- quick exits !
        if (info.lt.0) return
        if (nrow.eq.0) return
        if (ncol.eq.0) return

    !------------- right-looking column-wise factorisation !

        do  jj = 1, min(nrow, ncol)

    !-------------------------------------- find pivot row !
            jp = jj
            da = abs(aa(jj,jj))
            ta = abs(aa(jj,jj))
            do  ii = jj+1, nrow
                if (ta.lt.abs(aa(ii,jj))) then
                    jp  = ii
                    ta  = abs(aa(ii,jj))
                end if
            end do
            if (ta.lt.(1000.d+0*da)) then
    !-------------------------------------- bias onto diag !
                jp  = jj
            end if

            ipiv(jj)= jp

            if (abs(aa(jp,jj)).gt.0.0d+0) then

    !-------------------------------------- swap pivot row !
                if (jp.ne.jj) then
                do  kk =  1 , ncol
                    ta        = aa(jp,kk)
                    aa(jp,kk) = aa(jj,kk)
                    aa(jj,kk) = ta
                end do
                end if

    !-------------------------------------- scal pivot col !
                do  ii = jj+1, nrow
                    aa(ii,jj) = aa(ii,jj) / aa(jj,jj)
                end do

            else if (info.eq.0) then

    !-------------------------- store first singular pivot !
                info = jj

            end if

            if (jj.lt.min(nrow,ncol)) then
            
    !-------------------------- update trailing sub-matrix !
                do  kk = jj+1, ncol
                do  ii = jj+1, nrow
                    aa(ii,kk) = aa(ii,kk) - aa(ii,jj) &
                &                         * aa(jj,kk)
                end do
                end do
            
            end if

        end do !------------------ jj = 1, min(nrow, ncol) !

        return

    end subroutine lufactor

    !----------------------------------------------------------------!
    ! LUSOLVER: solve a linear system via LU factorisation.          !
    !----------------------------------------------------------------!
    !   TRAN - Transpose flag {'N' or 'T'}.                          !
    !   AA   - NROW-by-NCOL matrix, storing the in-place LU facoris- !
    !          ation returned bu LUFACTOR.                           !
    !          AA is allocated with leading dimension ADIM, such th- !
    !          at AA = AA(1:ADIM,:).                                 !
    !   IPIV - NROW-by-1 array of row permutations due to pivoting.  !
    !   XX   - NROW/NCOL-by-NRHS matrix, storing the set of RHS vec- !
    !          tors on entry - the set of solution vectors on exit.  !
    !          XX is allocated with leading dimension XDIM, such th- !
    !          at XX = XX(1:XDIM,:).                                 !
    !----------------------------------------------------------------!
    subroutine lusolver(tran,nrow,ncol,nrhs,aa,adim,ipiv,xx,xdim)

        implicit none

    !------------------------------------------- arguments !
        character , intent(in) :: tran
        integer, intent(in)    :: nrow,ncol,nrhs
        integer, intent(in)    :: adim,xdim
        integer, intent(in)    :: ipiv(*)
        real*8 , intent(in)    :: aa(adim,*)
        real*8 , intent(inout) :: xx(xdim,*)

        if (tran.eq.'N') then
    !--------------- solve OP(A) * X = X, with OP(A) = A^1 !

    !------------------------------ apply row-permutations !
            call rowswaps(nrhs, xx, xdim, &
        &           +1, nrow, ipiv, +1)

    !----------------------------------- solve L^1 * X = X !
            call trisolve('L', 'N', 'D' , &
        &           nrow, nrhs, aa, adim, xx, xdim)

    !----------------------------------- solve U^1 * X = X !
            call trisolve('U', 'N', 'N' , &
        &           nrow, nrhs, aa, adim, xx, xdim)

        else
    !--------------- solve OP(A) * X = X, with OP(A) = A^t !

    !----------------------------------- solve U^t * X = X !
            call trisolve('U', 'N', 'N' , &
        &           nrow, nrhs, aa, adim, xx, xdim)

    !----------------------------------- solve L^t * X = X !
            call trisolve('L', 'N', 'N' , &
        &           nrow, nrhs, aa, adim, xx, xdim)

    !------------------------------ apply row-permutations !
            call rowswaps(nrhs, xx, xdim, &
        &           +1, ncol, ipiv, -1)

        end if
        
        return

    end subroutine lusolver

    end module matrix



