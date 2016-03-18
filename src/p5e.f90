
    !
    ! P5E.f90: degree-5 interpolation for edge values/slopes.
    !
    ! Darren Engwirda 
    ! 17-Mar-2016
    ! engwirda [at] mit [dot] edu
    !

    !----------------------------------------------------------------!
    ! P5E: approximate edge values with degree-5 polynomials.        !
    !----------------------------------------------------------------!
    !   DELX - NPOS-1-by-1 array of grid coordinates. Grid spacing   !
    !          can be non-uniform.                                   !
    !   FDAT - NDOF-by-NVAR-by-NPOS-1 array of cell-wise moments for !
    !          NVAR discrete variables.                              !
    !   BCLO - NVAR-by-1 array of boundary conditions for NVAR disc- !
    !          rete variables at the lower endpoint.                 !
    !   BCUP - NVAR-by-1 array of boundary conditions for NVAR disc- !
    !          rete variables at the upper endpoint.                 !
    !   EDGE - 2-by-NVAR-by-NPOS-1 array of interpolated edge values !
    !          for each grid-cell. EDGE(1,:,IPOS) and EDGE(2,:,IPOS) !
    !          store the lower/upper edge values for grid-cell IPOS. !
    !   DFDX - 2-by-NVAR-by-NPOS-1 array of interpolated edge slopes !
    !          for each grid-cell. DFDX(1,:,IPOS) and DFDX(2,:,IPOS) !
    !          store the lower/upper edge slopes for grid-cell IPOS. !
    !----------------------------------------------------------------!
    subroutine p5e(npos,nvar,ndof,delx,fdat, &
        &          bclo,bcup,edge,dfdx,opts)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        type (rcon_ends), intent(in) :: bclo(:)
        type (rcon_ends), intent(in) :: bcup(:)
        real*8 , intent(out) :: edge(:,:,:)
        real*8 , intent(out) :: dfdx(:,:,:)
        class(rcon_opts), intent(in) :: opts

    !------------------------------------------- variables !
        integer :: ipos,ivar,idof,info
        integer :: head,tail
        real*8  :: fnxx,dfxx,xhat
        real*8  :: xmap(-3:+3)
        real*8  :: fhat(+6, nvar)
        real*8  :: ivec(+6,-3:+3)
        integer :: ipiv(+6)
        real*8  :: cmat(+6,+6)
        
        character, parameter :: NTRA = 'N'
        integer  , parameter :: NSIZ = +6
   
        head = +4; tail = npos-3

        if (npos.le.6) then
    !----- default to reduced order if insufficient points !
            call p3e (npos,nvar,ndof, &
        &             delx,fdat,bclo, &
        &             bcup,edge,dfdx, &
        &             opts)
        end if
        if (npos.le.6) return

    !------ impose value/slope B.C.'s about lower endpoint !

        call lbc(npos,nvar,ndof,delx, &
        &        fdat,bclo,edge,dfdx)

    !------ impose value/slope B.C.'s about upper endpoint !

        call ubc(npos,nvar,ndof,delx, &
        &        fdat,bcup,edge,dfdx)

    ! Reconstruct edge-centred 6th-order polynomials. Com- !
    ! pute values/slopes at edges directly. Full-order ex- !
    ! trapolation at endpoints.                            !

        do  ipos = head , tail

            if (size(delx).eq.+1) then
            
    !--------------- reconstruction: constant grid-spacing !

            do  ivar = +1, nvar

                fnxx = &
        &+( 1.d0 / 60.d0) * fdat(+1,ivar,ipos-3) &
        &-( 8.d0 / 60.d0) * fdat(+1,ivar,ipos-2) &
        &+(37.d0 / 60.d0) * fdat(+1,ivar,ipos-1) &
        &+(37.d0 / 60.d0) * fdat(+1,ivar,ipos+0) &
        &-( 8.d0 / 60.d0) * fdat(+1,ivar,ipos+1) &
        &+( 1.d0 / 60.d0) * fdat(+1,ivar,ipos+2)

            edge(2,ivar,ipos-1) = fnxx
            edge(1,ivar,ipos+0) = fnxx

                dfxx = &
        &-( 1.d0 / 90.d0) * fdat(+1,ivar,ipos-3) &
        &+( 5.d0 / 36.d0) * fdat(+1,ivar,ipos-2) &
        &-(49.d0 / 36.d0) * fdat(+1,ivar,ipos-1) &
        &+(49.d0 / 36.d0) * fdat(+1,ivar,ipos+0) &
        &-( 5.d0 / 36.d0) * fdat(+1,ivar,ipos+1) &
        &+( 1.d0 / 90.d0) * fdat(+1,ivar,ipos+2)

            dfxx = dfxx / delx(+1)

            dfdx(2,ivar,ipos-1) = dfxx
            dfdx(1,ivar,ipos+0) = dfxx

            end do
            
    !--------------- reconstruction: variable grid-spacing !

            else

            xhat = +.5d0 * (delx(ipos-1)&
        &                +  delx(ipos+0))

            xmap(-3) = -(delx(ipos-3)&
        &              + delx(ipos-2)&
        &              + delx(ipos-1)) / xhat
            xmap(-2) = -(delx(ipos-2)&
        &              + delx(ipos-1)) / xhat
            xmap(-1) = - delx(ipos-1)  / xhat
            xmap(+0) = + 0.d0
            xmap(+1) = + delx(ipos+0)  / xhat
            xmap(+2) = +(delx(ipos+0)&
        &              + delx(ipos+1)) / xhat
            xmap(+3) = +(delx(ipos+0)&
        &              + delx(ipos+1)&
        &              + delx(ipos+2)) / xhat

    !--------------------------- calc. integral basis vec. !

            call poly_basis(-1,+6,xmap(  -3), &
        &                         ivec(:,-3))
            call poly_basis(-1,+6,xmap(  -2), &
        &                         ivec(:,-2))
            call poly_basis(-1,+6,xmap(  -1), &
        &                         ivec(:,-1))
            call poly_basis(-1,+6,xmap(  +0), &
        &                         ivec(:,+0))
            call poly_basis(-1,+6,xmap(  +1), &
        &                         ivec(:,+1))
            call poly_basis(-1,+6,xmap(  +2), &
        &                         ivec(:,+2))
            call poly_basis(-1,+6,xmap(  +3), &
        &                         ivec(:,+3))

    !--------------------------- linear system: lhs matrix !

            do  idof = +1 , +6

            cmat(1,idof) = ivec(idof,-2) &
        &                - ivec(idof,-3)
            cmat(2,idof) = ivec(idof,-1) &
        &                - ivec(idof,-2)
            cmat(3,idof) = ivec(idof,+0) &
        &                - ivec(idof,-1)
            cmat(4,idof) = ivec(idof,+1) &
        &                - ivec(idof,+0)
            cmat(5,idof) = ivec(idof,+2) &
        &                - ivec(idof,+1)
            cmat(6,idof) = ivec(idof,+3) &
        &                - ivec(idof,+2)

            end do

    !--------------------------- linear system: rhs vector !

            do  ivar = +1, nvar

            fhat(+1,ivar) = delx(ipos-3) &
        &       * fdat(+1,ivar,ipos-3) / xhat
            fhat(+2,ivar) = delx(ipos-2) &
        &       * fdat(+1,ivar,ipos-2) / xhat
            fhat(+3,ivar) = delx(ipos-1) &
        &       * fdat(+1,ivar,ipos-1) / xhat
            fhat(+4,ivar) = delx(ipos+0) &
        &       * fdat(+1,ivar,ipos+0) / xhat
            fhat(+5,ivar) = delx(ipos+1) &
        &       * fdat(+1,ivar,ipos+1) / xhat
            fhat(+6,ivar) = delx(ipos+2) &
        &       * fdat(+1,ivar,ipos+2) / xhat
        
            end do

    !------------------------- factor/solve linear systems !

            call lufactor(NSIZ,NSIZ,cmat,NSIZ, &
        &                 ipiv,info)

            call lusolver(NTRA,NSIZ,NSIZ,nvar, &
        &                 cmat,NSIZ,ipiv,fhat, &
        &                 NSIZ)

    !------------------------- eval. f(x), dfdx(x) at edge !

            do  ivar = +1, nvar

            edge(2,ivar,ipos-1) = fhat(1,ivar)
            dfdx(2,ivar,ipos-1) = fhat(2,ivar) &
        &                       / xhat
        
            edge(1,ivar,ipos+0) = fhat(1,ivar)       
            dfdx(1,ivar,ipos+0) = fhat(2,ivar) &
        &                       / xhat

            end do

            end if

        end do

        return

    end subroutine p5e
    
    
