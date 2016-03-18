
    !
    ! P1E.f90: degree-1 interpolation for edge values/slopes.
    !
    ! Darren Engwirda 
    ! 17-Mar-2016
    ! engwirda [at] mit [dot] edu
    !

	!----------------------------------------------------------------!
    ! P1E: approximate edge values with degree-1 polynomials.        !
    !----------------------------------------------------------------!
    !   DELX - NPOS-1-by-1 array of grid coordinates. Grid spacing   !
    !          can be non-uniform.                                   !
    !   FDAT - NDOF-by-NVAR-by-NPOS-1 array of cell-wise moments for !
    !          NVAR discrete variables.                              !
    !   EDGE - 2-by-NVAR-by-NPOS-1 array of interpolated edge values !
    !          for each grid-cell. EDGE(1,:,IPOS) and EDGE(2,:,IPOS) !
    !          store the lower/upper edge values for grid-cell IPOS. !
    !   DFDX - 2-by-NVAR-by-NPOS-1 array of interpolated edge slopes !
    !          for each grid-cell. DFDX(1,:,IPOS) and DFDX(2,:,IPOS) !
    !          store the lower/upper edge slopes for grid-cell IPOS. !
    !----------------------------------------------------------------!
    subroutine p1e(npos,nvar,ndof,delx,fdat, &
        &          edge,dfdx)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: edge(:,:,:)
        real*8 , intent(out) :: dfdx(:,:,:)

    !------------------------------------------- variables !
        integer :: ipos,ivar, &
        &          head,tail
        real*8  :: fnxx,dfxx,dd10
        real*8  :: delh(-1:+0)

        head = +2; tail = npos-1

        if (npos.lt.2) return
        if (npos.eq.2) then
    !----- default to reduced order if insufficient points !
        do  ivar = 1,nvar
            edge(1,1,ivar) = fdat(1,1,ivar)
            edge(2,1,ivar) = fdat(1,1,ivar)
            edge(1,2,ivar) = fdat(1,1,ivar)
            edge(2,2,ivar) = fdat(1,1,ivar)

            dfdx(1,1,ivar) = 0.d0
            dfdx(2,1,ivar) = 0.d0
            dfdx(1,2,ivar) = 0.d0
            dfdx(2,2,ivar) = 0.d0
        end do
        end if
        if (npos.le.2) return
   
    ! Reconstruct edge-centred 2nd-order polynomials. Com- !
    ! pute values/slopes at edges directly. Full-order ex- !
    ! trapolation at endpoints.                            !

        do  ipos = head , tail

            if (size(delx).gt.+1) then  ! variable spacing

            delh(-1) = delx(ipos -1)
            delh(+0) = delx(ipos +0)

            dd10 = delh(-1)+delh(+0)

    !------------- interpolate values/slopes at ii-th edge !

            do  ivar = +1, nvar

                fnxx = &
        &     + delh(+0) * fdat(1,ivar,ipos-1) &
        &     + delh(-1) * fdat(1,ivar,ipos+0)

                dfxx = &
        &     - 4.0d+000 * fdat(1,ivar,ipos-1) &
        &     + 4.0d+000 * fdat(1,ivar,ipos+0)

                fnxx = fnxx / dd10
                dfxx = dfxx / dd10

                edge(2,ivar,ipos-1) = fnxx
                edge(1,ivar,ipos+0) = fnxx

                dfdx(2,ivar,ipos-1) = dfxx
                dfdx(1,ivar,ipos+0) = dfxx

            end do

            else

    !------------- interpolate values/slopes at ii-th edge !

            do  ivar = +1, nvar

                fnxx = &
        &     + delx(+1) * fdat(1,ivar,ipos-1) &
        &     + delx(+1) * fdat(1,ivar,ipos+0)

                dfxx = &
        &     - 4.0d+000 * fdat(1,ivar,ipos-1) &
        &     + 4.0d+000 * fdat(1,ivar,ipos+0)

                fnxx = fnxx / delx(1)
                dfxx = dfxx / delx(1)

                edge(2,ivar,ipos-1) = fnxx
                edge(1,ivar,ipos+0) = fnxx

                dfdx(2,ivar,ipos-1) = dfxx
                dfdx(1,ivar,ipos+0) = dfxx

            end do
            
            end if

        end do

    !- impose low-order value/slope B.C.'s about endpoints !

        do  ivar = +1, nvar

            edge(1,ivar,head-1) = fdat(1,ivar,head-1)
            dfdx(1,ivar,head-1) = 0.d0

            edge(2,ivar,tail+0) = fdat(1,ivar,tail+0)
            dfdx(2,ivar,tail+0) = 0.d0

        end do
        
        return

    end subroutine p1e



