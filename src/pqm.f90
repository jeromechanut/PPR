
    !
    ! PQM.f90: piecewise quartic reconstruction method.
    !
    ! Darren Engwirda 
    ! 17-Mar-2016
    ! engwirda [at] mit [dot] edu
    !

    !----------------------------------------------------------------!
    ! POLYROOT: find roots of a degree-2 polynomial.                 !
    !----------------------------------------------------------------!
    !   aa * xx**2 + bb * xx + cc = 0                                !
    !----------------------------------------------------------------!
    function polyroot(aa,bb,cc,xx) result(haveroot)

        implicit none

    !------------------------------------------- arguments !
        real*8, intent( in) :: aa,bb,cc
        real*8, intent(out) :: xx(1:2)

    !------------------------------------------- variables !
        logical :: haveroot
        real*8 :: sq,ia,a0,b0,c0,x0

        real*8, parameter :: rt = +1.d-14

        a0 = abs(aa)
        b0 = abs(bb)
        c0 = abs(cc)

        sq = bb * bb - 4.0d+0 * aa * cc

        if (sq .ge. 0.0d+0) then

            sq = sqrt (sq)

            xx(1) =  - bb + sq
            xx(2) =  - bb - sq

            x0 = max(abs(xx(1)), &
        &            abs(xx(2)))

            if (a0 .gt. (rt*x0)) then

    !-------------------------------------- degree-2 roots !
    
            haveroot =  .TRUE.

            ia = 0.5d+0   / aa

            xx(1) = xx(1) * ia
            xx(2) = xx(2) * ia
            
            else
            
            if (b0 .gt. (rt*c0)) then

    !-------------------------------------- degree-1 roots !

            haveroot =  .TRUE.
            
            xx(1) =  - cc / bb
            xx(2) =  - cc / bb
            
            else
            
            haveroot = .FALSE.
            
            end if
            
            end if

        else

            haveroot = .FALSE.

        end if

        return

    end function polyroot

    !----------------------------------------------------------------!
    ! PQMFN: piecewise quartic interpolation.                        !
    !----------------------------------------------------------------!
    subroutine pqmfn(ff00,ffll,ffrr,fell,ferr, &
        &            dell,derr,dfds,fhat,mono, &
        &            ilim)

        implicit none

    !------------------------------------------- arguments !
        real*8 , intent(inout) :: ff00
        real*8 , intent(inout) :: ffll,ffrr
        real*8 , intent(inout) :: fell,ferr
        real*8 , intent(inout) :: dell,derr
        real*8 , intent(in)    :: dfds(-1:+1)
        real*8 , intent(out)   :: fhat(+1:+5)
        integer, intent(in)    :: ilim
        integer, intent(out)   :: mono
          
    !------------------------------------------- variables !
        integer :: bind
        real*8  :: aval,bval,cval
        real*8  :: iflx(+1:+2)
        real*8  :: dflx(+1:+2)
        
        mono    = 0
        
        if (ilim.eq.null_limit) then
        
    !-------------------------------- "null" slope-limiter !
        
            fhat(1) = &
        & + (30.d+0 / 16.d+0) * ff00 &
        & - ( 7.d+0 / 16.d+0) *(ferr+fell) &
        & + ( 1.d+0 / 16.d+0) *(derr-dell)
            fhat(2) = &
        & + ( 3.d+0 /  4.d+0) *(ferr-fell) &
        & - ( 1.d+0 /  4.d+0) *(derr+dell)
            fhat(3) = &
        & - (30.d+0 /  8.d+0) * ff00 &
        & + (15.d+0 /  8.d+0) *(ferr+fell) &
        & - ( 3.d+0 /  8.d+0) *(derr-dell)
            fhat(4) = &
        & - ( 1.d+0 /  4.d+0) *(ferr-fell  &
        &                      -derr-dell)
            fhat(5) = &
        & + (30.d+0 / 16.d+0) * ff00 &
        & - (15.d+0 / 16.d+0) *(ferr+fell) &
        & + ( 5.d+0 / 16.d+0) *(derr-dell)

        else
        
    !-------------------------------- "mono" slope-limiter !
        
            if((ffrr - ff00) * & 
        &      (ff00 - ffll) .le. 0.d+0) then

    !----------------------------------- "flatten" extrema !
    
                mono = +1

                fhat(1) = ff00
                fhat(2) = 0.d0
                fhat(3) = 0.d0
                fhat(4) = 0.d0
                fhat(5) = 0.d0
                 
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

    !----------------------------------- limit edge slopes !
        
            if((dell * dfds(-1)) .lt. 0.d+0) then

                mono = +1

                dell = dfds(-1)
     
            end if

            if((derr * dfds(+1)) .lt. 0.d+0) then

                mono = +1

                derr = dfds(+1)

            end if
        
    !----------------------------------- limit cell values !
        
            fhat(1) = &
        & + (30.d+0 / 16.d+0) * ff00 &
        & - ( 7.d+0 / 16.d+0) *(ferr+fell) &
        & + ( 1.d+0 / 16.d+0) *(derr-dell)
            fhat(2) = &
        & + ( 3.d+0 /  4.d+0) *(ferr-fell) &
        & - ( 1.d+0 /  4.d+0) *(derr+dell)
            fhat(3) = &
        & - (30.d+0 /  8.d+0) * ff00 &
        & + (15.d+0 /  8.d+0) *(ferr+fell) &
        & - ( 3.d+0 /  8.d+0) *(derr-dell)
            fhat(4) = &
        & - ( 1.d+0 /  4.d+0) *(ferr-fell  &
        &                      -derr-dell)
            fhat(5) = &
        & + (30.d+0 / 16.d+0) * ff00 &
        & - (15.d+0 / 16.d+0) *(ferr+fell) &
        & + ( 5.d+0 / 16.d+0) *(derr-dell)
 
    !------------------ calc. inflexion via 2nd-derivative !
            aval = 12.d+0 * fhat(5)
            bval =  6.d+0 * fhat(4)
            cval =  2.d+0 * fhat(3)

            if ( polyroot(aval,bval,cval,iflx) ) then

                bind = +0

                if ( ( iflx(1) .gt. -1.d+0 ) &
        &      .and. ( iflx(1) .lt. +1.d+0 ) ) then

    !------------------ check for non-monotonic inflection !
                dflx(1) =       fhat(2) &
        &     + iflx(1)       * fhat(3) * 2.d+0 &
        &     +(iflx(1) ** 2) * fhat(4) * 3.d+0 &
        &     +(iflx(1) ** 3) * fhat(5) * 4.d+0

                if (dflx(1)*dfds(+0) .lt. 0.d+0) then

                    if (abs(dell) &
        &          .lt. abs(derr) ) then

                        bind = -1

                    else

                        bind = +1

                    end if

                end if

                end if
                
                if ( ( iflx(2) .gt. -1.d+0 ) &
        &      .and. ( iflx(2) .lt. +1.d+0 ) ) then

    !------------------ check for non-monotonic inflection !
                dflx(2) =       fhat(2) &
        &     + iflx(2)       * fhat(3) * 2.d+0 &
        &     +(iflx(2) ** 2) * fhat(4) * 3.d+0 &
        &     +(iflx(2) ** 3) * fhat(5) * 4.d+0

                if (dflx(2)*dfds(+0) .lt. 0.d+0) then

                    if (abs(dell) &
        &          .lt. abs(derr) ) then

                        bind = -1

                    else

                        bind = +1

                    end if

                end if

                end if
                
    !------------------ pop non-monotone inflexion to edge !
              
                if (bind .eq. -1) then

    !------------------ pop inflection points onto -1 edge !
    
                    mono = +2

                    derr = &
        &       - ( 5.d+0 / 1.d+0) * ff00 &
        &       + ( 3.d+0 / 1.d+0) * ferr &
        &       + ( 2.d+0 / 1.d+0) * fell
                    dell = &
        &       + ( 5.d+0 / 3.d+0) * ff00 &
        &       - ( 1.d+0 / 3.d+0) * ferr &
        &       - ( 4.d+0 / 3.d+0) * fell

                    if (dell*dfds(-1) .lt. 0.d+0) then

                        dell = 0.d+0

                        ferr = &
        &           + ( 5.d+0 / 1.d+0) * ff00 &
        &           - ( 4.d+0 / 1.d+0) * fell
                        derr = &
        &           + (10.d+0 / 1.d+0) * ff00 &
        &           - (10.d+0 / 1.d+0) * fell

                    end if

                    if (derr*dfds(+1) .lt. 0.d+0) then

                        derr = 0.d+0

                        fell = &
        &           + ( 5.d+0 / 2.d+0) * ff00 &
        &           - ( 3.d+0 / 2.d+0) * ferr
                        dell = &
        &           - ( 5.d+0 / 3.d+0) * ff00 &
        &           + ( 5.d+0 / 3.d+0) * ferr

                    end if

                end if

                if (bind .eq. +1) then

    !------------------ pop inflection points onto -1 edge !
                    mono = +2

                    derr = &
        &       - ( 5.d+0 / 3.d+0) * ff00 &
        &       + ( 4.d+0 / 3.d+0) * ferr &
        &       + ( 1.d+0 / 3.d+0) * fell
                    dell = &
        &       + ( 5.d+0 / 1.d+0) * ff00 &
        &       - ( 2.d+0 / 1.d+0) * ferr &
        &       - ( 3.d+0 / 1.d+0) * fell

                    if (dell*dfds(-1) .lt. 0.d+0) then

                        dell = 0.d+0

                        ferr = &
        &           + ( 5.d+0 / 2.d+0) * ff00 &
        &           - ( 3.d+0 / 2.d+0) * fell
                        derr = &
        &           + ( 5.d+0 / 3.d+0) * ff00 &
        &           - ( 5.d+0 / 3.d+0) * fell

                    end if

                    if (derr*dfds(+1) .lt. 0.d+0) then

                        derr = 0.d+0

                        fell = &
        &           + ( 5.d+0 / 1.d+0) * ff00 &
        &           - ( 4.d+0 / 1.d+0) * ferr
                        dell = &
        &           - (10.d+0 / 1.d+0) * ff00 &
        &           + (10.d+0 / 1.d+0) * ferr

                    end if

                end if
                
            end if ! if polyroot(aval,bval,cval,iflx)

            if (mono .ge. +2) then

    !------------------ re-assemble coefficients on demand !

            fhat(1) = &
        & + (30.d+0 / 16.d+0) * ff00 &
        & - ( 7.d+0 / 16.d+0) *(ferr+fell) &
        & + ( 1.d+0 / 16.d+0) *(derr-dell)
            fhat(2) = &
        & + ( 3.d+0 /  4.d+0) *(ferr-fell) &
        & - ( 1.d+0 /  4.d+0) *(derr+dell)
            fhat(3) = &
        & - (30.d+0 /  8.d+0) * ff00 &
        & + (15.d+0 /  8.d+0) *(ferr+fell) &
        & - ( 3.d+0 /  8.d+0) *(derr-dell)
            fhat(4) = &
        & - ( 1.d+0 /  4.d+0) *(ferr-fell  &
        &                      -derr-dell)
            fhat(5) = &
        & + (30.d+0 / 16.d+0) * ff00 &
        & - (15.d+0 / 16.d+0) *(ferr+fell) &
        & + ( 5.d+0 / 16.d+0) *(derr-dell)
 
            end if
 
        end if

        return

    end subroutine pqmfn

    !----------------------------------------------------------------!
    ! PQM: piecewise quartic method.                                 !
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
    !   DFDX - 2-by-NVAR-by-NPOS-1 array of interpolated edge slopes !
    !          for each grid-cell. DFDX(1,:,IPOS) and DFDX(2,:,IPOS) !
    !          store the lower/upper edge values for grid-cell IPOS. !
    !   OSCL - 2-by-NVAR-by-NPOS-1 array of oscillation indicators.  !
    !   ILIM - slope-limiter selection: NULL_LIMIT, MONO_LIMIT,      !
    !                                   WENO_LIMIT                   !
    !   HALO - stencil width.                                        !
    !----------------------------------------------------------------!
    subroutine pqm(npos,nvar,ndof,delx,fdat, &
        &          fhat,edge,dfdx,oscl,ilim, &
        &          halo)

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,nvar,ndof
        integer, intent(in)  :: ilim,halo
        real*8 , intent(out) :: fhat(:,:,:)
        real*8 , intent(out) :: oscl(:,:,:)
        real*8 , intent(in)  :: delx(:)
        real*8 , intent(in)  :: fdat(:,:,:)
        real*8 , intent(in)  :: edge(:,:,:)
        real*8 , intent(in)  :: dfdx(:,:,:)

    !------------------------------------------- variables !
        integer :: ipos,ivar,mono
        integer :: iill,iirr
        integer :: head,tail
        real*8  :: xhat
        real*8  :: ff00,ffll,ffrr
        real*8  :: fell,ferr
        real*8  :: dell,derr
        real*8  :: dfds(-1:+1,nvar)
        real*8  :: uhat(5)
        real*8  :: lhat(5)
        real*8  :: wval(2)
        logical :: uniform
        
        head = +1; tail = npos - 1

        if (npos.eq.2) then
    !----- default to reduced order if insufficient points !
        do  ivar = +1, nvar
            fhat(1,ivar,+1) = fdat(1,ivar,+1)
            fhat(2,ivar,+1) = 0.d0
            fhat(3,ivar,+1) = 0.d0
            fhat(4,ivar,+1) = 0.d0
            fhat(5,ivar,+1) = 0.d0
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

        uniform = (size(delx) .eq. +1)

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

            if (uniform) then
            xhat = delx(    +1) * 0.5d+0
            else
            xhat = delx(ipos+0) * 0.5d+0
            end if

            do  ivar = +1, nvar

    !----------------------------- cell mean + edge values !
                ff00 = fdat(1,ivar,ipos)
                ffll = fdat(1,ivar,iill)
                ffrr = fdat(1,ivar,iirr)
               
                dell = dfdx(1,ivar,ipos)
                derr = dfdx(2,ivar,ipos)
                dell = dell * xhat
                derr = derr * xhat

                fell = edge(1,ivar,ipos)
                ferr = edge(2,ivar,ipos)

    !----------------------------- calc. cell-wise profile !
                select case(ilim)
                case (null_limit)
    
    !----------------------------- calc. unlimited profile !              
                call pqmfn(ff00,ffll,ffrr ,&
        &                  fell,ferr,dell ,&
        &                  derr, &
        &                  dfds(   :,ivar),&
        &                  uhat, &
        &                  mono,null_limit)
                
    !----------------------------- pref. unlimited profile !
                wval(1) = +1.d+0
                wval(2) = +0.d+0
                
                case (mono_limit)
                
    !----------------------------- calc. monotonic profile !              
                call pqmfn(ff00,ffll,ffrr ,&
        &                  fell,ferr,dell ,&
        &                  derr, &
        &                  dfds(   :,ivar),&
        &                  lhat, &
        &                  mono,mono_limit)
                
    !----------------------------- pref. monotonic profile !
                wval(1) = +0.d+0
                wval(2) = +1.d+0
                
                case (weno_limit)
  
    !----------------------------- calc. unlimited profile !              
                call pqmfn(ff00,ffll,ffrr ,&
        &                  fell,ferr,dell ,&
        &                  derr, &
        &                  dfds(   :,ivar),&
        &                  uhat, &
        &                  mono,null_limit)
  
    !----------------------------- calc. monotonic profile !              
                call pqmfn(ff00,ffll,ffrr ,&
        &                  fell,ferr,dell ,&
        &                  derr, &
        &                  dfds(   :,ivar),&
        &                  lhat, &
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
                fhat(4,ivar,ipos) = &
        &           wval(1) * uhat(4) + &
        &           wval(2) * lhat(4)
                fhat(5,ivar,ipos) = &
        &           wval(1) * uhat(5) + &
        &           wval(2) * lhat(5)

            end do

        end do
        
        return

    end subroutine pqm


