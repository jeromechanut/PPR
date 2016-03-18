
    !
    ! PBC.f90: boundary conditions for edge interpolation.
    !
    ! Darren Engwirda 
    ! 17-Mar-2016
    ! engwirda [at] mit [dot] edu
    !

    !----------------------------------------------------------------!
    ! LBC_IMPL: setup a single boundary condition type at the lower  !
    ! endpoint. Edges/slopes at first three edges are set.           !
    !----------------------------------------------------------------!
    subroutine lbc_impl(npos,nvar,ndof,delx,fdat, &
        &               bcon,bopt,edge,dfdx)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        integer, intent( in) :: bopt
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: edge(:,:,:)
        real*8 , intent(out) :: dfdx(:,:,:)
        type(rcon_ends), intent(in) :: bcon(:)

    !------------------------------------------- variables !
        integer :: ivar,idof,isel, &
        &          head,tail,info
        real*8  :: xhat
        real*8  :: delh(-1:+1)
        real*8  :: xmap(-1:+2)
        real*8  :: bvec(+3,-1:+2)
        real*8  :: gvec(+3,-1:+2)
        real*8  :: cmat(+3,+3)
        real*8  :: fhat(+3, nvar)
        integer :: ipiv(+3)
        real*8  :: eval(-1:+2)
        real*8  :: gval(-1:+2)

        character, parameter :: NTRA = 'N'
        integer  , parameter :: NSIZ = +3

        head = +2; tail = npos - 2

        if (size(delx).gt.+1) then

    !------------------ mean grid spacing about ii-th cell !

        xhat = delx(head) * 0.5d+0

    !------------------ grid spacing for all stencil cells !

        delh(-1) = delx(head-1)
        delh(+0) = delx(head+0)
        delh(+1) = delx(head+1)

        else
        
    !------------------ mean grid spacing about ii-th cell !

        xhat = delx(  +1) * 0.5d+0

    !------------------ grid spacing for all stencil cells !

        delh(-1) = delx(    +1)
        delh(+0) = delx(    +1)
        delh(+1) = delx(    +1)
        
        end if

    !---------- local coordinate mapping for stencil edges !

        xmap(-1) =-(delh(-1)+0.5d0*delh(0)) / xhat
        xmap(+0) = -1.d0
        xmap(+1) = +1.d0
        xmap(+2) = (delh(+1)+0.5d0*delh(0)) / xhat

    !------------ linear system: lhs reconstruction matrix !

        select case (bopt)
        case ( bcon_loose)

            call poly_basis(-1,+3,xmap(-1),bvec(:,-1))
            call poly_basis(-1,+3,xmap(+0),bvec(:,+0))
            call poly_basis(-1,+3,xmap(+1),bvec(:,+1))
            call poly_basis(-1,+3,xmap(+2),bvec(:,+2))

            do  idof = +1 , +3

                cmat(1,idof) = bvec(idof,+0) &
        &                    - bvec(idof,-1)
                cmat(2,idof) = bvec(idof,+1) &
        &                    - bvec(idof,+0)
                cmat(3,idof) = bvec(idof,+2) &
        &                    - bvec(idof,+1)

            end do

        case ( bcon_value)

            call poly_basis(+0,+3,xmap(-1),gvec(:,-1))
            call poly_basis(-1,+3,xmap(-1),bvec(:,-1))
            call poly_basis(-1,+3,xmap(+0),bvec(:,+0))
            call poly_basis(-1,+3,xmap(+1),bvec(:,+1))

            do  idof = +1 , +3

                cmat(1,idof) = bvec(idof,+0) &
        &                    - bvec(idof,-1)
                cmat(2,idof) = bvec(idof,+1) &
        &                    - bvec(idof,+0)

                cmat(3,idof) = gvec(idof,-1)

            end do

        case ( bcon_slope)

            call poly_basis(+1,+3,xmap(-1),gvec(:,-1))
            call poly_basis(-1,+3,xmap(-1),bvec(:,-1))
            call poly_basis(-1,+3,xmap(+0),bvec(:,+0))
            call poly_basis(-1,+3,xmap(+1),bvec(:,+1))

            do  idof = +1 , +3

                cmat(1,idof) = bvec(idof,+0) &
        &                    - bvec(idof,-1)
                cmat(2,idof) = bvec(idof,+1) &
        &                    - bvec(idof,+0)

                cmat(3,idof) = gvec(idof,-1)

            end do

        end select

    !------------ linear system: rhs reconstruction vector !

        isel  = +0

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bopt)  then

            isel = isel + 1

            select case (bopt)
            case ( bcon_loose)

            fhat(1,isel) = &
        &   delh(-1) * fdat(1,ivar,head-1) / xhat
            fhat(2,isel) = &
        &   delh(+0) * fdat(1,ivar,head+0) / xhat
            fhat(3,isel) = &
        &   delh(+1) * fdat(1,ivar,head+1) / xhat

            case ( bcon_value)

            fhat(1,isel) = &
        &   delh(-1) * fdat(1,ivar,head-1) / xhat
            fhat(2,isel) = &
        &   delh(+0) * fdat(1,ivar,head+0) / xhat

            fhat(3,isel) = bcon(ivar)%value

            case ( bcon_slope)

            fhat(1,isel) = &
        &   delh(-1) * fdat(1,ivar,head-1) / xhat
            fhat(2,isel) = &
        &   delh(+0) * fdat(1,ivar,head+0) / xhat

            fhat(3,isel) = &
        &       bcon(ivar)%slope  * xhat

            end select

        end if

        end do

    !------------------------- factor/solve linear systems !

        call lufactor(NSIZ,NSIZ,cmat,NSIZ, &
        &             ipiv,info)

        call lusolver(NTRA,NSIZ,NSIZ,nvar, &
        &             cmat,NSIZ,ipiv,fhat, &
        &             NSIZ)

    !------------- extrapolate values/slopes at lower edge !

        isel  = +0

        call poly_basis(+0,+3,xmap(-1),bvec(:,-1))
        call poly_basis(+0,+3,xmap(+0),bvec(:,+0))
        call poly_basis(+0,+3,xmap(+1),bvec(:,+1))
        
        call poly_basis(+1,+3,xmap(-1),gvec(:,-1))
        call poly_basis(+1,+3,xmap(+0),gvec(:,+0))
        call poly_basis(+1,+3,xmap(+1),gvec(:,+1))

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bopt)  then

            isel = isel  + 1

            eval(-1) = poly_value(+0,+3 , &
        &                   fhat(:,isel), &
        &                   xmap(   -1) , &
        &                   bvec(:, -1))
            eval(+0) = poly_value(+0,+3 , &
        &                   fhat(:,isel), &
        &                   xmap(   +0) , &
        &                   bvec(:, +0))
            eval(+1) = poly_value(+0,+3 , &
        &                   fhat(:,isel), &
        &                   xmap(   +1) , &
        &                   bvec(:, +1))

            gval(-1) = poly_value(+1,+3 , &
        &                   fhat(:,isel), &
        &                   xmap(   -1) , &
        &                   gvec(:, -1))
            gval(+0) = poly_value(+1,+3 , &
        &                   fhat(:,isel), &
        &                   xmap(   +0) , &
        &                   gvec(:, +0))
            gval(+1) = poly_value(+1,+3 , &
        &                   fhat(:,isel), &
        &                   xmap(   +1) , &
        &                   gvec(:, +1))

            edge(+1,ivar,head-1) = eval(-1)
            edge(+2,ivar,head-1) = eval(+0)
            edge(+1,ivar,head+0) = eval(+0)
            edge(+2,ivar,head+0) = eval(+0)
            edge(+1,ivar,head+1) = eval(+0)

            dfdx(+1,ivar,head-1) = gval(-1) &
        &                        / xhat
            dfdx(+2,ivar,head-1) = gval(+0) &
        &                        / xhat
            dfdx(+1,ivar,head+0) = gval(+0) &
        &                        / xhat
            dfdx(+2,ivar,head+0) = gval(+1) &
        &                        / xhat
            dfdx(+1,ivar,head+1) = gval(+1) &
        &                        / xhat

        end if

        end do
        
        return

    end subroutine lbc_impl

    !----------------------------------------------------------------!
    ! LBC: set boundary conditions at lower endpoint - edges/slopes  !
    ! at first two edges are set.                                    !
    !----------------------------------------------------------------!
    !   DELX - NPOS-1-by-1 array of grid coordinates. Grid spacing   !
    !          can be non-uniform.                                   !
    !   FDAT - NDOF-by-NVAR-by-NPOS-1 array of cell-wise moments for !
    !          NVAR discrete variables.                              !
    !   BCON - NVAR-by-1 array of boundary conditions for NVAR disc- !
    !          rete variables, specified as an array of RCON_ENDS.   !
    !   EDGE - 2-by-NVAR-by-NPOS-1 array of interpolated edge values !
    !          for each grid-cell. EDGE(1,:,IPOS) and EDGE(2,:,IPOS) !
    !          store the lower/upper edge values for grid-cell IPOS. !
    !   DFDX - 2-by-NVAR-by-NPOS-1 array of interpolated edge slopes !
    !          for each grid-cell. DFDX(1,:,IPOS) and DFDX(2,:,IPOS) !
    !          store the lower/upper edge slopes for grid-cell IPOS. !
    !----------------------------------------------------------------!
    subroutine lbc(npos,nvar,ndof,delx,fdat, &
        &          bcon,edge,dfdx)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        type(rcon_ends), intent(in) :: bcon(:)
        real*8 , intent(out) :: edge(:,:,:)
        real*8 , intent(out) :: dfdx(:,:,:)

    !------------------------------------------- variables !
        integer :: ivar
        integer :: num_loose,num_value, &
        &          num_slope

        num_loose= +0
        num_value= +0
        num_slope= +0

        do  ivar = +1, nvar

            select case (bcon(ivar)%bcopt)
    !-------------------------------------- count BC types !
            case (bcon_loose)
                num_loose = num_loose + 1
            case (bcon_value)
                num_value = num_value + 1
            case (bcon_slope)
                num_slope = num_slope + 1

            end select

        end do

        if (num_loose.gt.0) then
    !---------------------------- setup "unset" conditions !
            call lbc_impl(npos,nvar,ndof, &
        &                 delx,fdat,bcon, &
        &                 bcon_loose ,  &
        &                 edge,dfdx)

        end if

        if (num_value.gt.0) then
    !---------------------------- setup "value" conditions !
            call lbc_impl(npos,nvar,ndof, &
        &                 delx,fdat,bcon, &
        &                 bcon_value ,  &
        &                 edge,dfdx)

        end if

        if (num_slope.gt.0) then
    !---------------------------- setup "slope" conditions !
            call lbc_impl(npos,nvar,ndof, &
        &                 delx,fdat,bcon, &
        &                 bcon_slope ,  &
        &                 edge,dfdx)

        end if
        
        return

    end subroutine lbc

    !----------------------------------------------------------------!
    ! UBC_IMPL: setup a single boundary condition type at the upper  !
    ! endpoint. Edges/slopes at final three edges are set.           !
    !----------------------------------------------------------------!
    subroutine ubc_impl(npos,nvar,ndof,delx,fdat, &
        &               bcon,bopt,edge,dfdx)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        integer, intent( in) :: bopt
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        type(rcon_ends), intent(in) :: bcon(:)
        real*8 , intent(out) :: edge(:,:,:)
        real*8 , intent(out) :: dfdx(:,:,:)

    !------------------------------------------- variables !
        integer :: ivar,idof,isel, &
        &          head,tail,info
        real*8  :: xhat
        real*8  :: delh(-1:+1)
        real*8  :: xmap(-1:+2)
        real*8  :: bvec(+3,-1:+2)
        real*8  :: gvec(+3,-1:+2)
        real*8  :: cmat(+3,+3)
        real*8  :: fhat(+3, nvar)
        integer :: ipiv(+3)
        real*8  :: eval(-1:+2)
        real*8  :: gval(-1:+2)

        character, parameter :: NTRA = 'N'
        integer  , parameter :: NSIZ = +3

        head = +2; tail = npos - 2

        if (size(delx).gt.+1) then

    !------------------ mean grid spacing about ii-th cell !

        xhat = delx(tail) * 0.5d+0

    !------------------ grid spacing for all stencil cells !

        delh(-1) = delx(tail-1)
        delh(+0) = delx(tail+0)
        delh(+1) = delx(tail+1)

        else

    !------------------ mean grid spacing about ii-th cell !

        xhat = delx(  +1) * 0.5d+0

    !------------------ grid spacing for all stencil cells !

        delh(-1) = delx(    +1)
        delh(+0) = delx(    +1)
        delh(+1) = delx(    +1)

        
        end if

    !---------- local coordinate mapping for stencil edges !

        xmap(-1) =-(delh(-1)+0.5d0*delh(0)) / xhat
        xmap(+0) = -1.d0
        xmap(+1) = +1.d0
        xmap(+2) = (delh(+1)+0.5d0*delh(0)) / xhat

    !------------ linear system: lhs reconstruction matrix !

        select case (bopt)
        case ( bcon_loose)

            call poly_basis(-1,+3,xmap(-1),bvec(:,-1))
            call poly_basis(-1,+3,xmap(+0),bvec(:,+0))
            call poly_basis(-1,+3,xmap(+1),bvec(:,+1))
            call poly_basis(-1,+3,xmap(+2),bvec(:,+2))

            do  idof = +1 , +3

                cmat(1,idof) = bvec(idof,+0) &
        &                    - bvec(idof,-1)
                cmat(2,idof) = bvec(idof,+1) &
        &                    - bvec(idof,+0)
                cmat(3,idof) = bvec(idof,+2) &
        &                    - bvec(idof,+1)

            end do

        case ( bcon_value)

            call poly_basis(+0,+3,xmap(-1),gvec(:,+2))
            call poly_basis(-1,+3,xmap(+0),bvec(:,+0))
            call poly_basis(-1,+3,xmap(+1),bvec(:,+1))
            call poly_basis(-1,+3,xmap(+2),bvec(:,+2))

            do  idof = +1 , +3

                cmat(1,idof) = gvec(idof,+0)

                cmat(2,idof) = bvec(idof,+1) &
        &                    - bvec(idof,+0)
                cmat(3,idof) = bvec(idof,+2) &
        &                    - bvec(idof,+1)

            end do

        case ( bcon_slope)

            call poly_basis(+1,+3,xmap(+2),gvec(:,+2))
            call poly_basis(+0,+3,xmap(+0),bvec(:,+0))
            call poly_basis(+0,+3,xmap(+1),bvec(:,+1))
            call poly_basis(+0,+3,xmap(+2),bvec(:,+2))

            do  idof = +1 , +3

                cmat(1,idof) = gvec(idof,+2)

                cmat(2,idof) = bvec(idof,+1) &
        &                    - bvec(idof,+0)
                cmat(3,idof) = bvec(idof,+2) &
        &                    - bvec(idof,+1)

            end do

        end select

    !------------ linear system: rhs reconstruction vector !

        isel  = +0

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bopt)  then

            isel = isel + 1

            select case (bopt)
            case ( bcon_loose)

            fhat(1,isel) = &
        &   delh(-1) * fdat(1,ivar,tail-1) / xhat
            fhat(2,isel) = &
        &   delh(+0) * fdat(1,ivar,tail+0) / xhat
            fhat(3,isel) = &
        &   delh(+1) * fdat(1,ivar,tail+1) / xhat

            case ( bcon_value)

            fhat(1,isel) = bcon(ivar)%value

            fhat(2,isel) = &
        &   delh(+0) * fdat(1,ivar,tail+0) / xhat
            fhat(3,isel) = &
        &   delh(+1) * fdat(1,ivar,tail+1) / xhat

            case ( bcon_slope)

            fhat(1,isel) = &
        &       bcon(ivar)%slope  * xhat

            fhat(2,isel) = &
        &   delh(+0) * fdat(1,ivar,tail+0) / xhat
            fhat(3,isel) = &
        &   delh(+1) * fdat(1,ivar,tail+1) / xhat

            end select

        end if

        end do

    !------------------------- factor/solve linear systems !

        call lufactor(NSIZ,NSIZ,cmat,NSIZ, &
        &             ipiv,info)

        call lusolver(NTRA,NSIZ,NSIZ,nvar, &
        &             cmat,NSIZ,ipiv,fhat, &
        &             NSIZ)

    !------------- extrapolate values/slopes at upper edge !

        isel  = +0

        call poly_basis(+0,+3,xmap(+0),bvec(:,+0))
        call poly_basis(+0,+3,xmap(+1),bvec(:,+1))
        call poly_basis(+0,+3,xmap(+2),bvec(:,+2))
        
        call poly_basis(+1,+3,xmap(+0),gvec(:,+0))
        call poly_basis(+1,+3,xmap(+1),gvec(:,+1))
        call poly_basis(+1,+3,xmap(+2),gvec(:,+2))
        
        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bopt)  then

            isel = isel  + 1

            eval(+0) = poly_value(+0,+3 , &
        &                   fhat(:,isel), &
        &                   xmap(   +0) , &
        &                   bvec(:, +0))
            eval(+1) = poly_value(+0,+3 , &
        &                   fhat(:,isel), &
        &                   xmap(   +1) , &
        &                   bvec(:, +1))
            eval(+2) = poly_value(+0,+3 , &
        &                   fhat(:,isel), &
        &                   xmap(   +2) , &
        &                   bvec(:, +2))

            gval(+0) = poly_value(+1,+3 , &
        &                   fhat(:,isel), &
        &                   xmap(   +0) , &
        &                   gvec(:, +0))
            gval(+1) = poly_value(+1,+3 , &
        &                   fhat(:,isel), &
        &                   xmap(   +1) , &
        &                   gvec(:, +1))
            gval(+2) = poly_value(+1,+3 , &
        &                   fhat(:,isel), &
        &                   xmap(   +2) , &
        &                   gvec(:, +2))

            edge(+2,ivar,tail-1) = eval(+0)
            edge(+1,ivar,tail+0) = eval(+0)
            edge(+2,ivar,tail+0) = eval(+1)
            edge(+1,ivar,tail+1) = eval(+1)
            edge(+2,ivar,tail+1) = eval(+2)

            dfdx(+2,ivar,tail-1) = gval(+0) & 
        &                        / xhat
            dfdx(+1,ivar,tail+0) = gval(+0) &
        &                        / xhat
            dfdx(+2,ivar,tail+0) = gval(+1) & 
        &                        / xhat
            dfdx(+1,ivar,tail+1) = gval(+1) &
        &                        / xhat
            dfdx(+2,ivar,tail+1) = gval(+2) &
        &                        / xhat

        end if

        end do
        
        return

    end subroutine ubc_impl

    !----------------------------------------------------------------!
    ! UBC: set boundary conditions at upper endpoint - edges/slopes  !
    ! at final two edges are set.                                    !
    !----------------------------------------------------------------!
    !   DELX - NPOS-1-by-1 array of grid coordinates. Grid spacing   !
    !          can be non-uniform.                                   !
    !   FDAT - NDOF-by-NVAR-by-NPOS-1 array of cell-wise moments for !
    !          NVAR discrete variables.                              !
    !   BCON - NVAR-by-1 array of boundary conditions for NVAR disc- !
    !          rete variables, specified as an array of RCON_ENDS.   !
    !   EDGE - 2-by-NVAR-by-NPOS-1 array of interpolated edge values !
    !          for each grid-cell. EDGE(1,:,IPOS) and EDGE(2,:,IPOS) !
    !          store the lower/upper edge values for grid-cell IPOS. !
    !   DFDX - 2-by-NVAR-by-NPOS-1 array of interpolated edge slopes !
    !          for each grid-cell. DFDX(1,:,IPOS) and DFDX(2,:,IPOS) !
    !          store the lower/upper edge slopes for grid-cell IPOS. !
    !----------------------------------------------------------------!
    subroutine ubc(npos,nvar,ndof,delx,fdat, &
        &          bcon,edge,dfdx)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        type(rcon_ends), intent(in) :: bcon(:)
        real*8 , intent(out) :: edge(:,:,:)
        real*8 , intent(out) :: dfdx(:,:,:)

    !------------------------------------------- variables !
        integer :: ivar
        integer :: num_loose,num_value, &
        &          num_slope

        num_loose= +0
        num_value= +0
        num_slope= +0

        do  ivar = +1, nvar

            select case (bcon(ivar)%bcopt)
    !-------------------------------------- count BC types !
            case (bcon_loose)
                num_loose = num_loose + 1
            case (bcon_value)
                num_value = num_value + 1
            case (bcon_slope)
                num_slope = num_slope + 1

            end select

        end do

        if (num_loose.gt.0) then
    !---------------------------- setup "unset" conditions !
            call ubc_impl(npos,nvar,ndof, &
        &                 delx,fdat,bcon, &
        &                 bcon_loose ,  &
        &                 edge,dfdx)

        end if

        if (num_value.gt.0) then
    !---------------------------- setup "value" conditions !
            call ubc_impl(npos,nvar,ndof, &
        &                 delx,fdat,bcon, &
        &                 bcon_value ,  &
        &                 edge,dfdx)

        end if

        if (num_slope.gt.0) then
    !---------------------------- setup "slope" conditions !
            call ubc_impl(npos,nvar,ndof, &
        &                 delx,fdat,bcon, &
        &                 bcon_slope ,  &
        &                 edge,dfdx)

        end if

        return

    end subroutine ubc
    
    
    
