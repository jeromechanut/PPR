
    !
    ! PPM.f90: piecewise parabolic reconstruction method.
    !
    ! Darren Engwirda 
    ! 17-Mar-2016
    ! engwirda [at] mit [dot] edu
    !

    !----------------------------------------------------------------!
    ! PPMFN: piecewise parabolic interpolation.                      !
    !----------------------------------------------------------------!
    subroutine ppmfn(ff00,ffll,ffrr,fell,ferr, &
        &            dfds,fhat,mono,ilim)

        implicit none

    !------------------------------------------- arguments !
        real*8 , intent(inout) :: ff00
        real*8 , intent(inout) :: ffll,ffrr
        real*8 , intent(inout) :: fell,ferr
        real*8 , intent(in)    :: dfds(-1:+1)
        real*8 , intent(out)   :: fhat(+1:+3)
        integer, intent(in)    :: ilim
        integer, intent(out)   :: mono
          
    !------------------------------------------- variables !
        real*8  :: turn
        
        mono    = 0
        
        if (ilim.eq.null_limit) then
        
    !-------------------------------- "null" slope-limiter !
            fhat( 1 ) = &
        & + (3.0d+0 / 2.0d+0) * ff00 &
        & - (1.0d+0 / 4.0d+0) *(ferr+fell)
            fhat( 2 ) = &
        & + (1.0d+0 / 2.0d+0) *(ferr-fell)
            fhat( 3 ) = &
        & - (3.0d+0 / 2.0d+0) * ff00 &
        & + (3.0d+0 / 4.0d+0) *(ferr+fell)
    
        else
        
    !-------------------------------- "mono" slope-limiter !
        
            if((ffrr - ff00) * & 
        &      (ff00 - ffll) .le. 0.d+0) then

    !----------------------------------- "flatten" extrema !
    
                mono = +1

                fhat(1) = ff00
                fhat(2) = 0.d0
                fhat(3) = 0.d0
                 
                return
                  
            end if

    !----------------------------------- limit edge values !
    
            if((ffll - fell) * &
        &      (fell - ff00) .le. 0.d+0) then

                mono = +1
           
                fell = ff00 - dfds(0)

            end if

            if((ffrr - ferr) * &
        &      (ferr - ff00) .le. 0.d+0) then

                mono = +1

                ferr = ff00 + dfds(0)
             
            end if
       
    !----------------------------------- limit cell values !
    
            fhat( 1 ) = &
        & + (3.0d+0 / 2.0d+0) * ff00 &
        & - (1.0d+0 / 4.0d+0) *(ferr+fell)
            fhat( 2 ) = &
        & + (1.0d+0 / 2.0d+0) *(ferr-fell)
            fhat( 3 ) = &
        & - (3.0d+0 / 2.0d+0) * ff00 &
        & + (3.0d+0 / 4.0d+0) *(ferr+fell)
       
            if (abs(fhat(3)) .gt. & 
        &       abs(fhat(2))*.5d+0) then

            turn = -0.5d+0 * fhat(2) &
        &                  / fhat(3)

            if ((turn .ge. -1.d+0)&
        &  .and.(turn .le. +0.d+0)) then

                mono = +2

    !--------------------------- push TURN onto lower edge !
    
                ferr = +3.0d+0 * ff00 &
        &              -2.0d+0 * fell

            end if

            if ((turn .gt. +0.d+0)&
        &  .and.(turn .le. +1.d+0)) then

                mono = +2

    !--------------------------- push TURN onto upper edge !
    
                fell = +3.0d+0 * ff00 &
        &              -2.0d+0 * ferr

            end if
          
            end if
       
            if (mono .ge. +2) then

    !------------------ re-assemble coefficients on demand !

            fhat( 1 ) = &
        & + (3.0d+0 / 2.0d+0) * ff00 &
        & - (1.0d+0 / 4.0d+0) *(ferr+fell)
            fhat( 2 ) = &
        & + (1.0d+0 / 2.0d+0) *(ferr-fell)
            fhat( 3 ) = &
        & - (3.0d+0 / 2.0d+0) * ff00 &
        & + (3.0d+0 / 4.0d+0) *(ferr+fell)
       
            end if
        
        end if
        
        return
    
    end subroutine ppmfn

    !----------------------------------------------------------------!
    ! PPM: piecewise parabolic method.                               !
    !----------------------------------------------------------------!
    !   DELX - NPOS-1-by-1 array of grid coordinates. Grid spacing   !
    !          can be non-uniform.                                   !
    !   FDAT - NDOF-by-NVAR-by-NPOS-1 array of cell-wise moments for !
    !          NVAR discrete variables.                              !
    !   FHAT - NDOF-by-NVAR-by-NPOS-1 array of piece-wise polynomial !
    !          coefficients. See EVAL_FUNC, etc for additional info. !
    !   EDGE - 2-by-NVAR-by-NPOS-1 array of interpolated edge values !
    !          for each grid-cell. EDGE(1,:,IPOS) and EDGE(2,:,IPOS) !
    !          store the lower/upper edge values for grid-cell IPOS. !
    !   OSCL - 2-by-NVAR-by-NPOS-1 array of oscillation indicators.  !
    !   ILIM - slope-limiter selection: NULL_LIMIT, MONO_LIMIT,      !
    !                                   WENO_LIMIT                   !
    !   HALO - stencil width.                                        !
    !----------------------------------------------------------------!
    subroutine ppm(npos,nvar,ndof,delx,fdat, &
        &          fhat,edge,oscl,ilim,halo)

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,nvar,ndof
        integer, intent(in)  :: ilim,halo
        real*8 , intent(out) :: fhat(:,:,:)
        real*8 , intent(out) :: oscl(:,:,:)
        real*8 , intent(in)  :: delx(:)
        real*8 , intent(in)  :: fdat(:,:,:)
        real*8 , intent(in)  :: edge(:,:,:)

    !------------------------------------------- variables !
        integer :: ipos,ivar
        integer :: mono
        integer :: iill,iirr
        integer :: head,tail
        real*8  :: ff00
        real*8  :: ffll,ffrr
        real*8  :: fell,ferr
        real*8  :: dfds(-1:+1,nvar)
        real*8  :: uhat(3)
        real*8  :: lhat(3)
        real*8  :: wval(2)
        
        head = +1; tail = npos - 1

        if (npos.eq.2) then
    !----- default to reduced order if insufficient points !
        do  ivar = +1, nvar
            fhat(1,ivar,+1) = fdat(1,ivar,+1)
            fhat(2,ivar,+1) = 0.d0
            fhat(3,ivar,+1) = 0.d0
        end do
        end if
        if (npos.le.2) return

        do  ipos = +1, npos-1

            if (ilim.eq.weno_limit) then
    !------------------- oscillator indicator on each cell !
            call oscli(npos,nvar,ndof, &
        &              delx,fdat,ipos, &
        &              oscl)
            end if

        end do

    !------------------- reconstruct function on each cell !

        uhat = +0.d+0
        lhat = +0.d+0

        do  ipos = +1, npos-1

    !----------------------- calc. piecewise linear slopes !
            call pls(npos,nvar,ndof,delx,&
        &            fdat, &
        &            dfds,ipos,mono_limit)

            iill = max(head,ipos-1)
            iirr = min(tail,ipos+1)

            do  ivar = +1, nvar

    !----------------------------- cell mean + edge values !
                ff00 = fdat(1,ivar,ipos)
                ffll = fdat(1,ivar,iill)
                ffrr = fdat(1,ivar,iirr)
    
                fell = edge(1,ivar,ipos)
                ferr = edge(2,ivar,ipos)

    !----------------------------- calc. cell-wise profile !
                select case(ilim)
                case (null_limit)
    
    !----------------------------- calc. unlimited profile !              
                call ppmfn(ff00,ffll,ffrr ,&
        &                  fell,ferr, &
        &                  dfds(   :,ivar),&
        &                  uhat     , &
        &                  mono,null_limit)
                
    !----------------------------- pref. unlimited profile !
                wval(1) = +1.d+0
                wval(2) = +0.d+0
                
                case (mono_limit)
                
    !----------------------------- calc. monotonic profile !              
                call ppmfn(ff00,ffll,ffrr ,&
        &                  fell,ferr, &
        &                  dfds(   :,ivar),&
        &                  lhat     , &
        &                  mono,mono_limit)
                
    !----------------------------- pref. monotonic profile !
                wval(1) = +0.d+0
                wval(2) = +1.d+0
                
                case (weno_limit)
  
    !----------------------------- calc. unlimited profile !              
                call ppmfn(ff00,ffll,ffrr ,&
        &                  fell,ferr, &
        &                  dfds(   :,ivar),&
        &                  uhat     , &
        &                  mono,null_limit)
  
    !----------------------------- calc. monotonic profile !              
                call ppmfn(ff00,ffll,ffrr ,&
        &                  fell,ferr, &
        &                  dfds(   :,ivar),&
        &                  lhat     , &
        &                  mono,mono_limit)
                
                if (mono.gt.+0) then
  
    !----------------------------- calc. WENO-type weights !      
                call wenoi(npos,delx,oscl ,&
        &                  ipos,ivar,halo ,&
        &                  wval)
                
                else
                
    !----------------------------- pref. unlimited profile !
                wval(1) = +1.d+0
                wval(2) = +0.d+0
                
                end if
                
                end select
 
    !----------------------------- blend "null" and "mono" !               
                fhat(1,ivar,ipos) = &
        &           wval(1) * uhat(1) + &
        &           wval(2) * lhat(1)
                fhat(2,ivar,ipos) = &
        &           wval(1) * uhat(2) + &
        &           wval(2) * lhat(2)
                fhat(3,ivar,ipos) = &
        &           wval(1) * uhat(3) + &
        &           wval(2) * lhat(3)

            end do

        end do
        
        return

    end subroutine ppm



