
    !
    ! PLM.f90: piecewise linear reconstruction method.
    !
    ! Darren Engwirda 
    ! 17-Mar-2016
    ! engwirda [at] mit [dot] edu
    !

    !----------------------------------------------------------------!
    ! PLSV: calc. piecewise linear slopes.                           !
    !----------------------------------------------------------------!
    !   ==> variable grid-spacing variant.                           !
    !----------------------------------------------------------------!
    subroutine plsv(npos,nvar,ndof,delx,fdat, &
        &           dfds,ipos,ilim)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        integer, intent( in) :: ilim,ipos
        real*8 , intent( in) :: delx(:)
        real*8 , intent(out) :: dfds(-1:+1,nvar)
        real*8 , intent( in) :: fdat(:,:,:)

    !------------------------------------------- variables !
        integer :: ivar,head,tail
        real*8  :: hhll,hh00,hhrr
        real*8  :: ffll,ff00,ffrr
        real*8  :: fell,ferr,scal

    !------ construct centred slopes, extrapolate at edges !

        head = 1 ; tail = npos - 1

        if ((ipos.gt.head).and.(ipos.lt.tail)) then
        
    !------------------------------ centred approximations !
    
            hhll = delx(ipos-1)
            hh00 = delx(ipos+0)
            hhrr = delx(ipos+1)
    
            do  ivar = +1, nvar
    
    !------------------------------ compute each component !
    
                ffll = fdat(1,ivar,ipos-1)
                ff00 = fdat(1,ivar,ipos+0)
                ffrr = fdat(1,ivar,ipos+1)
    
                dfds(-1,ivar) = ff00-ffll
                dfds(+1,ivar) = ffrr-ff00
     
                if (dfds(-1,ivar) * &
        &           dfds(+1,ivar) .gt. 0.0d+0) then
    
    !---------------------------- calc. ll//rr edge values !
    
                fell = (hh00*ffll + hhll*ff00) & 
        &            / (hhll+hh00)        
                ferr = (hhrr*ff00 + hh00*ffrr) &
        &            / (hh00+hhrr)
    
    !---------------------------- calc. centred derivative !
                dfds(+0,ivar) = &
        &               0.5d+0 * (ferr - fell)
     
    !---------------------------- monotonic slope-limiting !
                scal = min(abs(dfds(-1,ivar)), &
        &                  abs(dfds(+1,ivar))) &
        &            / max(abs(dfds(+0,ivar)), &
                            epsilon(ff00))
                scal = min(scal, 1.0d+0)
    
                dfds(+0,ivar) = scal * dfds(+0,ivar)
    
                else
     
     !--------------------------- flatten if local extrema ! 
              
                dfds(+0,ivar) = +0.0d+0
                
                end if
                
    !---------------------------- scale onto local co-ord. !
                
                dfds(-1,ivar) = dfds(-1,ivar) &
        &               / (hhll + hh00) * hh00
                dfds(+1,ivar) = dfds(+1,ivar) &
        &               / (hh00 + hhrr) * hh00
                
            end do

        else

    !-------------------------------------- lower-endpoint !

        if (ipos.eq.head) then
        
            hhll = +0.0d+0
            hh00 = delx(ipos+0)
            hhrr = delx(ipos+1)
    
            do  ivar = +1, nvar
    
    !------------------------------ compute each component !
    
                ff00 = fdat(1,ivar,ipos+0)
                ffrr = fdat(1,ivar,ipos+1)
    
                dfds(+1,ivar) = &
        &           (ffrr - ff00) * hh00 &
        &         / (hh00 + hhrr)
    
                dfds(-1,ivar) = + 0.0d+0
                dfds(+0,ivar) = + 0.0d+0
                    
            end do
                
        end if

    !-------------------------------------- upper-endpoint !

        if (ipos.eq.tail) then
  
            hhll = delx(ipos-1)
            hh00 = delx(ipos+0)
            hhrr = +0.0d+0
            
            do  ivar = +1, nvar
    
    !------------------------------ compute each component !
    
                ffll = fdat(1,ivar,ipos-1)
                ff00 = fdat(1,ivar,ipos+0)
                
                dfds(-1,ivar) = & 
        &           (ff00 - ffll) * hh00 &
        &         / (hhll + hh00)
        
                dfds(+0,ivar) = + 0.0d+0
                dfds(+1,ivar) = + 0.0d+0
                
            end do

        end if
        
        end if

        return

    end subroutine plsv
    
    !----------------------------------------------------------------!
    ! PLSC: calc. piecewise linear slopes.                           !
    !----------------------------------------------------------------!
    !   ==> constant grid-spacing variant.                           !
    !----------------------------------------------------------------!
    subroutine plsc(npos,nvar,ndof,delx,fdat, &
        &           dfds,ipos,ilim)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        integer, intent( in) :: ilim,ipos
        real*8 , intent( in) :: delx(+1)
        real*8 , intent(out) :: dfds(-1:+1,nvar)
        real*8 , intent( in) :: fdat(:,:,:)

    !------------------------------------------- variables !
        integer :: ivar,head,tail
        real*8  :: ffll,ff00,ffrr
        real*8  :: fell,ferr,scal

    !------ construct centred slopes, extrapolate at edges !

        head = 1 ; tail = npos - 1

        if ((ipos.gt.head).and.(ipos.lt.tail)) then
        
    !------------------------------ centred approximations !
    
            do  ivar = +1, nvar
    
    !------------------------------ compute each component !
    
                ffll = fdat(1,ivar,ipos-1)
                ff00 = fdat(1,ivar,ipos+0)
                ffrr = fdat(1,ivar,ipos+1)
    
                dfds(-1,ivar) = ff00-ffll
                dfds(+1,ivar) = ffrr-ff00
     
                if (dfds(-1,ivar) * &
        &           dfds(+1,ivar) .gt. 0.0d+0) then
    
    !---------------------------- calc. ll//rr edge values !
    
                fell = (ffll + ff00) * 0.5d+0        
                ferr = (ff00 + ffrr) * 0.5d+0
    
    !---------------------------- calc. centred derivative !
                dfds(+0,ivar) = &
        &               0.5d+0 * (ferr - fell)
     
    !---------------------------- monotonic slope-limiting !
                scal = min(abs(dfds(-1,ivar)), &
        &                  abs(dfds(+1,ivar))) &
        &            / max(abs(dfds(+0,ivar)), &
                            epsilon(ff00))
                scal = min(scal, 1.0d+0)
    
                dfds(+0,ivar) = scal * dfds(+0,ivar)
    
                else
     
     !--------------------------- flatten if local extrema ! 
              
                dfds(+0,ivar) = +0.0d+0
                
                end if
                
    !---------------------------- scale onto local co-ord. !
                
                dfds(-1,ivar) = &
        &           0.5d+0 * dfds(-1,ivar) 
                dfds(+1,ivar) = &
        &           0.5d+0 * dfds(+1,ivar)
                
            end do

        else

    !-------------------------------------- lower-endpoint !

        if (ipos.eq.head) then
        
            do  ivar = +1, nvar
    
    !------------------------------ compute each component !
    
                ff00 = fdat(1,ivar,ipos+0)
                ffrr = fdat(1,ivar,ipos+1)
    
                dfds(+1,ivar) = &
        &         (ffrr - ff00) * 0.5d+0
    
                dfds(-1,ivar) = + 0.0d+0
                dfds(+0,ivar) = + 0.0d+0
                    
            end do
                
        end if

    !-------------------------------------- upper-endpoint !

        if (ipos.eq.tail) then
  
            do  ivar = +1, nvar
    
    !------------------------------ compute each component !
    
                ffll = fdat(1,ivar,ipos-1)
                ff00 = fdat(1,ivar,ipos+0)
                
                dfds(-1,ivar) = & 
        &         (ff00 - ffll) * 0.5d+0
        
                dfds(+0,ivar) = + 0.0d+0
                dfds(+1,ivar) = + 0.0d+0
                
            end do

        end if
        
        end if

        return

    end subroutine plsc

    !----------------------------------------------------------------!
    ! PLS: calc. piecewise linear slopes.                           !
    !----------------------------------------------------------------!
    subroutine pls(npos,nvar,ndof,delx,fdat, &
        &          dfds,ipos,ilim)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        integer, intent( in) :: ilim,ipos
        real*8 , intent( in) :: delx(:)
        real*8 , intent(out) :: dfds(-1:+1,nvar)
        real*8 , intent( in) :: fdat(:,:,:)
  
        if (size(delx).gt.+1) then
        
    !------------------------------- variable grid-spacing !
        
            call plsv(npos,nvar,ndof,delx, &
        &             fdat,dfds, &
        &             ipos,ilim)
        
        else
        
    !------------------------------- constant grid-spacing !
        
            call plsc(npos,nvar,ndof,delx, &
        &             fdat,dfds, &
        &             ipos,ilim)
        
        end if
  
        return
        
    end subroutine pls

    !----------------------------------------------------------------!
    ! PLM: piecewise linear method.                                  !
    !----------------------------------------------------------------!
    !   DELX - NPOS-1-by-1 array of grid coordinates. Grid spacing   !
    !          can be non-uniform.                                   !
    !   FDAT - NDOF-by-NVAR-by-NPOS-1 array of cell-wise moments for !
    !          NVAR discrete variables.                              !
    !   FHAT - NDOF-by-NVAR-by-NPOS-1 array of piece-wise polynomial !
    !          coefficients. See EVAL_FUNC, etc for additional info. !
    !   ILIM - slope-limiter selection: NULL_LIMIT, MONO_LIMIT,      !
    !                                   WENO_LIMIT                   !
    !----------------------------------------------------------------!    
    subroutine plm(npos,nvar,ndof,delx,fdat, &
        &          fhat,ilim)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        integer, intent( in) :: ilim
        real*8 , intent( in) :: delx(:)
        real*8 , intent(out) :: fhat(:,:,:)
        real*8 , intent( in) :: fdat(:,:,:)

    !------------------------------------------- variables !
        integer :: ipos,ivar
        real*8  :: dfds(-1:+1,nvar)

        if (npos.eq.2) then
    !----------------------- reduce order if small stencil !
        do  ivar = +1, nvar
            fhat(1,ivar,+1) = fdat(1,ivar,1)
            fhat(2,ivar,+1) = 0.d0
        end do
        end if
        if (npos.le.2) return

        do  ipos = +1, npos-1
        
    !----------------------- calc. piecewise linear slopes !
            
            call pls(npos,nvar,ndof, &
            &        delx,fdat,      &
            &        dfds,ipos,ilim)

    !----------------------- store piecewise linear coeff. !
            
            do  ivar = +1, nvar

            fhat(1,ivar,ipos) = fdat(1,ivar,ipos)
            fhat(2,ivar,ipos) = dfds(     0,ivar)

            end do
            
        end do
        
        return

    end subroutine plm
    
    
    
