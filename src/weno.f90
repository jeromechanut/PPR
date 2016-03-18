
    !
    ! WENO.f90: oscillation indicators and WENO-type weights.
    !
    ! Darren Engwirda 
    ! 17-Mar-2016
    ! engwirda [at] mit [dot] edu
    !

    !----------------------------------------------------------------!
    ! WENOV: weighted non-oscillatory indicator weights.             !
    !----------------------------------------------------------------!
    !   ==> variable grid-spacing variant.                           !
    !----------------------------------------------------------------!
    subroutine wenov(npos,delx,oscl,ipos,ivar, &
        &            halo,omin,omax)

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,halo
        integer, intent(in)  :: ipos,ivar
        real*8 , intent(in)  :: delx(:)
        real*8 , intent(in)  :: oscl(:,:,:)
        real*8 , intent(out) :: omin,omax

    !------------------------------------------- variables !
        integer :: hpos
        integer :: head,tail
        integer :: imin,imax
        real*8  :: deli,delh
        real*8  :: hh00
        real*8  :: dfx1,dfx2
        real*8  :: oval
        
    !------------------- calc. lower//upper stencil bounds !    

        head = 1 ; tail = npos - 1

        imin = max(ipos-halo,head)
        imax = min(ipos+halo,tail)
      
    !------------------ find min/max indicators on stencil !
 
        dfx1 = oscl(1,ivar,ipos)
        dfx2 = oscl(2,ivar,ipos)
      
        hh00 = delx(ipos+0)
      
        oval = (hh00**1*dfx1)**2 &
        &    + (hh00**2*dfx2)**2
      
        omin = oval
        omax = oval
      
    !---------------------------------------- "lower" part !
      
        delh = 0.d0

        do  hpos = ipos-1, imin, -1
        
    !------------------ calc. derivatives centred on IPOS. !     
   
            deli = delx(hpos+0) &
        &        + delx(hpos+1)
        
            delh = delh + deli*.5d0
            
            dfx1 = oscl(1,ivar,hpos)
            dfx2 = oscl(2,ivar,hpos)

            dfx1 = dfx1 + dfx2*delh
    
    !------------------ indicator: NORM(H^N * D^N/DX^N(F)) !

            oval = (hh00**1*dfx1)**2 &
        &        + (hh00**2*dfx2)**2

            omin = min(omin, oval)
            omax = max(omax, oval)
        
        end do
      
    !---------------------------------------- "upper" part !     
      
        delh = 0.d0
        
        do  hpos = ipos+1, imax, +1
        
    !------------------ calc. derivatives centred on IPOS. !     
   
            deli = delx(hpos+0) &
        &        + delx(hpos-1)
        
            delh = delh - deli*.5d0
            
            dfx1 = oscl(1,ivar,hpos)
            dfx2 = oscl(2,ivar,hpos)

            dfx1 = dfx1 + dfx2*delh
    
    !------------------ indicator: NORM(H^N * D^N/DX^N(F)) !

            oval = (hh00**1*dfx1)**2 &
        &        + (hh00**2*dfx2)**2

            omin = min(omin, oval)
            omax = max(omax, oval)
        
        end do

        return

    end subroutine wenov
    
    !----------------------------------------------------------------!
    ! WENOC: weighted non-oscillatory indicator weights.             !
    !----------------------------------------------------------------!
    !   ==> constant grid-spacing variant.                           !
    !----------------------------------------------------------------!
    subroutine wenoc(npos,delx,oscl,ipos,ivar, &
        &            halo,omin,omax)

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,halo
        integer, intent(in)  :: ipos,ivar
        real*8 , intent(in)  :: delx(1)
        real*8 , intent(in)  :: oscl(:,:,:)
        real*8 , intent(out) :: omin,omax

    !------------------------------------------- variables !
        integer :: hpos
        integer :: head,tail
        integer :: imin,imax
        real*8  :: delh
        real*8  :: dfx1,dfx2
        real*8  :: oval
    
    !------------------- calc. lower//upper stencil bounds !    

        head = 1 ; tail = npos - 1

        imin = max(ipos-halo,head)
        imax = min(ipos+halo,tail)
      
    !------------------ find min/max indicators on stencil !
 
        dfx1 = oscl(1,ivar,ipos)
        dfx2 = oscl(2,ivar,ipos)
      
        oval = (2.d0**1*dfx1)**2 &
        &    + (2.d0**2*dfx2)**2
      
        omin = oval
        omax = oval
      
    !---------------------------------------- "lower" part !
      
        delh = 0.d0

        do  hpos = ipos-1, imin, -1
        
    !------------------ calc. derivatives centred on IPOS. !     
   
            delh = delh + 2.d0
            
            dfx1 = oscl(1,ivar,hpos)
            dfx2 = oscl(2,ivar,hpos)

            dfx1 = dfx1 + dfx2*delh
    
    !------------------ indicator: NORM(H^N * D^N/DX^N(F)) !

            oval = (2.d0**1*dfx1)**2 &
        &        + (2.d0**2*dfx2)**2

            omin = min(omin, oval)
            omax = max(omax, oval)
        
        end do
      
    !---------------------------------------- "upper" part !     
      
        delh = 0.d0
        
        do  hpos = ipos+1, imax, +1
        
    !------------------ calc. derivatives centred on IPOS. !     
   
            delh = delh - 2.d0
            
            dfx1 = oscl(1,ivar,hpos)
            dfx2 = oscl(2,ivar,hpos)

            dfx1 = dfx1 + dfx2*delh
    
    !------------------ indicator: NORM(H^N * D^N/DX^N(F)) !

            oval = (2.d0**1*dfx1)**2 &
        &        + (2.d0**2*dfx2)**2

            omin = min(omin, oval)
            omax = max(omax, oval)
        
        end do

        return

    end subroutine wenoc
    
    !----------------------------------------------------------------!
    ! WENOI: weighted non-oscillatory indicator weights.             !
    !----------------------------------------------------------------!
    !   DELX - NPOS-1-by-1 array of "old" grid coordinates. Grid sp- !
    !          acing can be non-uniform.                             !
    !   IPOS - gridcell index for which the indicator is calculated. !
    !   IVAR -                                                       !
    !   OSCL - 2-by-NVAR-by-NPOS-1 array of oscillation coefficients !
    !          OSCL(1,:,:) is the 1st derivative at the cell centre. !
    !          OSCL(2,:,:) is the 2nd derivative at the cell center. !
    !   HALO -                                                       !
    !----------------------------------------------------------------!
    subroutine wenoi(npos,delx,oscl,ipos,ivar, &
        &            halo,wval)

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,halo
        integer, intent(in)  :: ipos,ivar
        real*8 , intent(in)  :: delx(:)
        real*8 , intent(in)  :: oscl(:,:,:)
        real*8 , intent(out) :: wval(2)
    
    !------------------------------------------- variables !
        real*8 :: omin,omax,wsum
        
        real*8 , parameter :: zero = +1.d-16
    
        if (size(delx).gt.+1) then
        
    !------------------- use variable grid spacing variant !
        
        call wenov(npos,delx,oscl, &
        &          ipos,ivar,halo, &
        &          omin,omax)
        
        else
        
    !------------------- use constant grid spacing variant !
        
        call wenoc(npos,delx,oscl, &
        &          ipos,ivar,halo, &
        &          omin,omax)
        
        end if

    !------------------ compute WENO-style profile weights !

        omax = omax + zero
        omin = omin + zero

        wval(1) = +1.0d+6 / omax ** 3
        wval(2) = +1.0d+0 / omin ** 3

        wsum = wval(1) + wval(2)
        wval(1) = wval(1) / wsum
        wval(2) = wval(2) / wsum

        return

    end subroutine wenoi
    
    !----------------------------------------------------------------!
    ! OSCLV: cell-wise oscillation indicators.                       !
    !----------------------------------------------------------------!
    !   ==> variable grid-spacing variant.                           !
    !----------------------------------------------------------------!
    subroutine osclv(npos,nvar,ndof,delx,fdat, &
        &            ipos,oscl)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        integer, intent( in) :: ipos
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: oscl(:,:,:)

    !------------------------------------------- variables !
        integer :: head,tail
        integer :: ivar
        real*8  :: hhll,hhcc,hhrr
        real*8  :: hhmm,hhrc,hhlc
        real*8  :: cmat(2,3)

        head = +1; tail = npos-1

        if ((ipos.gt.head).and.(ipos.lt.tail)) then

    !--------------------------------------- centred point !

            hhll = delx(ipos-1)
            hhcc = delx(ipos+0)
            hhrr = delx(ipos+1)

            hhrc = hhrr + hhcc
            hhlc = hhll + hhcc
            hhmm = hhll + hhcc + hhrr

            cmat(1,1) = -(hhcc+2.d0*hhrr)/(hhlc*hhmm)
            cmat(1,2) = -(hhll-hhrr)* &
            &          (3.d0*hhcc+2.d0*(hhll+hhrr))/&
            &            (hhlc*hhrc*hhmm)
            cmat(1,3) = +(hhcc+2.d0*hhll)/(hhrc*hhmm)

            cmat(2,1) = +3.d0/(hhlc*hhmm)
            cmat(2,2) = -3.d0*(2.d0*hhcc+hhll+hhrr)/&
            &            (hhlc*hhrc*hhmm)
            cmat(2,3) = +3.d0/(hhrc*hhmm)

            do  ivar = 1, nvar

                oscl(1,ivar,ipos) = +1.d0 * ( &
            & + cmat(1,1)*fdat(1,ivar,ipos-1) &
            & + cmat(1,2)*fdat(1,ivar,ipos+0) &
            & + cmat(1,3)*fdat(1,ivar,ipos+1) )

                oscl(2,ivar,ipos) = +2.d0 * ( &
            & + cmat(2,1)*fdat(1,ivar,ipos-1) &
            & + cmat(2,2)*fdat(1,ivar,ipos+0) &
            & + cmat(2,3)*fdat(1,ivar,ipos+1) )

            end do

        else

        if (ipos.eq.head) then

    !-------------------------------------- lower endpoint !

            hhll = delx(ipos+0)
            hhcc = delx(ipos+1)
            hhrr = delx(ipos+2)

            hhrc = hhrr + hhcc
            hhlc = hhll + hhcc
            hhmm = hhll + hhcc + hhrr

            cmat(1,1) = -(hhcc+2.d0*hhrr)/(hhlc*hhmm)
            cmat(1,2) = -(hhll-hhrr)* &
            &          (3.d0*hhcc+2.d0*(hhll+hhrr))/&
            &            (hhlc*hhrc*hhmm)
            cmat(1,3) = +(hhcc+2.d0*hhll)/(hhrc*hhmm)

            cmat(2,1) = +3.d0/(hhlc*hhmm)
            cmat(2,2) = -3.d0*(2.d0*hhcc+hhll+hhrr)/&
            &            (hhlc*hhrc*hhmm)
            cmat(2,3) = +3.d0/(hhrc*hhmm)

            do  ivar = 1, nvar

                oscl(1,ivar,ipos) = +1.d0 * ( &
            & + cmat(1,1)*fdat(1,ivar,ipos+0) &
            & + cmat(1,2)*fdat(1,ivar,ipos+1) &
            & + cmat(1,3)*fdat(1,ivar,ipos+2) )

                oscl(2,ivar,ipos) = +2.d0 * ( &
            & + cmat(2,1)*fdat(1,ivar,ipos+0) &
            & + cmat(2,2)*fdat(1,ivar,ipos+1) &
            & + cmat(2,3)*fdat(1,ivar,ipos+2) )

                oscl(1,ivar,ipos) = &
            &   oscl(1,ivar,ipos) - &
            &   0.5d0 * hhlc * oscl(2,ivar,ipos)

            end do

        end if

        if (ipos.eq.tail) then

    !-------------------------------------- upper endpoint !

            hhll = delx(ipos-2)
            hhcc = delx(ipos-1)
            hhrr = delx(ipos-0)

            hhrc = hhrr + hhcc
            hhlc = hhll + hhcc
            hhmm = hhll + hhcc + hhrr

            cmat(1,1) = -(hhcc+2.d0*hhrr)/(hhlc*hhmm)
            cmat(1,2) = -(hhll-hhrr)* &
            &          (3.d0*hhcc+2.d0*(hhll+hhrr))/&
            &            (hhlc*hhrc*hhmm)
            cmat(1,3) = +(hhcc+2.d0*hhll)/(hhrc*hhmm)

            cmat(2,1) = +3.d0/(hhlc*hhmm)
            cmat(2,2) = -3.d0*(2.d0*hhcc+hhll+hhrr)/&
            &            (hhlc*hhrc*hhmm)
            cmat(2,3) = +3.d0/(hhrc*hhmm)

            do  ivar = 1, nvar

                oscl(1,ivar,ipos) = +1.d0 * ( &
            & + cmat(1,1)*fdat(1,ivar,ipos-2) &
            & + cmat(1,2)*fdat(1,ivar,ipos-1) &
            & + cmat(1,3)*fdat(1,ivar,ipos+0) )

                oscl(2,ivar,ipos) = +2.d0 * ( &
            & + cmat(2,1)*fdat(1,ivar,ipos-2) &
            & + cmat(2,2)*fdat(1,ivar,ipos-1) &
            & + cmat(2,3)*fdat(1,ivar,ipos+0) )

                oscl(1,ivar,ipos) = &
            &   oscl(1,ivar,ipos) + &
            &   0.5d0 * hhrc * oscl(2,ivar,ipos)

            end do

        end if

        end if
        
        return

    end subroutine osclv
    
    !----------------------------------------------------------------!
    ! OSCLV: cell-wise oscillation indicators.                       !
    !----------------------------------------------------------------!
    !   ==> constant grid-spacing variant.                           !
    !----------------------------------------------------------------!
    subroutine osclc(npos,nvar,ndof,delx,fdat, &
        &            ipos,oscl)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        integer, intent( in) :: ipos
        real*8 , intent( in) :: delx(1)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: oscl(:,:,:)

    !------------------------------------------- variables !
        integer :: head,tail,ivar
        
        head = +1; tail = npos -1

        if ((ipos.gt.head).and.(ipos.lt.tail)) then

    !--------------------------------------- centred point !

            do  ivar = 1, nvar

                oscl(1,ivar,ipos) = &
        &     + .25d+0 * fdat(1,ivar,ipos+1) &
        &     - .25d+0 * fdat(1,ivar,ipos-1)
            
                oscl(2,ivar,ipos) = &
        &     + .25d+0 * fdat(1,ivar,ipos+1) &
        &     - .50d+0 * fdat(1,ivar,ipos+0) &
        &     + .25d+0 * fdat(1,ivar,ipos-1)

            end do

        else

        if (ipos.eq.head) then

    !-------------------------------------- lower endpoint !

            do  ivar = 1, nvar

                oscl(1,ivar,ipos) = &
        &     + .25d+0 * fdat(1,ivar,ipos+2) &
        &     - .25d+0 * fdat(1,ivar,ipos+0)
            
                oscl(2,ivar,ipos) = &
        &     + .25d+0 * fdat(1,ivar,ipos+2) &
        &     - .50d+0 * fdat(1,ivar,ipos+1) &
        &     + .25d+0 * fdat(1,ivar,ipos+0)

                oscl(1,ivar,ipos) = & !! check this
                oscl(1,ivar,ipos) - &
                oscl(2,ivar,ipos) * 2.d+0

            end do

        end if

        if (ipos.eq.tail) then

    !-------------------------------------- upper endpoint !

            do  ivar = 1, nvar

                oscl(1,ivar,ipos) = &
        &     + .25d+0 * fdat(1,ivar,ipos+0) &
        &     - .25d+0 * fdat(1,ivar,ipos-2)
            
                oscl(2,ivar,ipos) = &
        &     + .25d+0 * fdat(1,ivar,ipos+0) &
        &     - .50d+0 * fdat(1,ivar,ipos-1) &
        &     + .25d+0 * fdat(1,ivar,ipos-2)

                oscl(1,ivar,ipos) = & !! check this
                oscl(1,ivar,ipos) + &
                oscl(2,ivar,ipos) * 2.d+0

            end do

        end if

        end if
        
        return

    end subroutine osclc
    
    !----------------------------------------------------------------!
    ! OSCLI: cell-wise oscillation indicators.                       !
    !----------------------------------------------------------------!
    !   DELX - NPOS-1-by-1 array of "old" grid coordinates. Grid sp- !
    !          acing can be non-uniform.                             !
    !   FDAT - NDOF-by-NVAR-by-NPOS-1 array of cell-wise moments for !
    !          NVAR discrete variables.                              !
    !   IPOS - gridcell index for which the indicator is calculated. !
    !   OSCL - 2-by-NVAR-by-NPOS-1 array of oscillation coefficients !
    !          OSCL(1,:,:) is the 1st derivative at the cell centre. !
    !          OSCL(2,:,:) is the 2nd derivative at the cell center. !
    !----------------------------------------------------------------!
    subroutine oscli(npos,nvar,ndof,delx,fdat, &
        &            ipos,oscl)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        integer, intent( in) :: ipos
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: oscl(:,:,:)

    !------------------------------------------- variables !
        integer :: ivar
        
        if (npos.lt.4) then
    !------------------------------- at least 3 grid-cells !
        do  ivar = +1, nvar
            oscl(1,ivar,ipos) = +0.d0
            oscl(2,ivar,ipos) = +0.d0
        end do
        end if
        if (npos.lt.4) return
        if (nvar.lt.1) return
        if (ndof.lt.1) return

        if (size(delx).gt.+1) then

    !------------------------------- variable grid-spacing !

            call osclv(npos,nvar,ndof,delx, &
        &              fdat,ipos,oscl)
        
        else

    !------------------------------- constant grid-spacing !
        
            call osclc(npos,nvar,ndof,delx, &
        &              fdat,ipos,oscl)
        
        end if

        return
        
    end subroutine oscli



