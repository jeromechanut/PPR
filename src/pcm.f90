
    !
    ! PCM.f90: piecewise constant reconstruction method.
    !
    ! Darren Engwirda 
    ! 17-Mar-2016
    ! engwirda [at] mit [dot] edu
    !

    !----------------------------------------------------------------!
    ! PCM: piecewise constant method.                                !
    !----------------------------------------------------------------!
    !   XPOS - NPOS-by-1 array of grid coordinates. Grid spacing can !
    !          be non-uniform.                                       !
    !   FDAT - NDOF-by-NVAR-by-NPOS-1 array of cell-wise moments for !
    !          NVAR discrete variables.                              !
    !   FHAT - NDOF-by-NVAR-by-NPOS-1 array of piece-wise polynomial !
    !          coefficients. See POLY_VALUE etc for additional info. !
    !----------------------------------------------------------------!
    
    subroutine pcm(npos,nvar,ndof,fdat,fhat)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent(out) :: fhat(:,:,:)
        real*8 , intent( in) :: fdat(:,:,:)

    !------------------------------------------- variables !
        integer:: ipos,ivar,idof

        do  ipos = +1, npos - 1
        do  ivar = +1, nvar + 0
        do  idof = +1, ndof + 0

            fhat(idof,ivar,ipos) = fdat(idof,ivar,ipos)

        end do
        end do
        end do
        
        return

    end subroutine pcm
    
    

