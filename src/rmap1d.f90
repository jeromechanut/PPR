
    !-------------------------------------------------------
    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !-------------------------------------------------------
    !

    !    
    ! RMAP1D.f90: high-order integral re-mapping operators.
    !
    ! Darren Engwirda 
    ! 07-Sep-2016
    ! ​de2363 [at] columbia [dot] edu
    !
    !

    subroutine rmap1d(npos,nnew,nvar,ndof,xpos, &
        &             xnew,fdat,fnew,bclo,bcup, &
        &             work,opts,tCPU)

    !
    ! NPOS  no. edges in old grid.
    ! NNEW  no. edges in new grid.
    ! NVAR  no. discrete variables to remap.
    ! NDOF  no. degrees-of-freedom per cell.
    ! XPOS  old grid edge positions. XPOS is a length NPOS
    !       array .
    ! XNEW  new grid edge positions. XNEW is a length NNEW
    !       array .
    ! FDAT  grid-cell moments on old grid. FNEW has SIZE = 
    !       NDOF-by-NVAR-by-NNEW-1 .
    ! FNEW  grid-cell moments on new grid. FNEW has SIZE = 
    !       NDOF-by-NVAR-by-NNEW-1 .
    ! BCLO  boundary condition at lower endpoint .
    ! BCHI  boundary condition at upper endpoint . 
    ! WORK  method work-space. See RCON-WORK for details .
    ! OPTS  method parameters. See RCON-OPTS for details .
    ! TCPU  method tcpu-timer.
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nnew
        integer, intent( in) :: nvar,ndof
        class(rmap_work), intent(inout):: work
        class(rmap_opts), intent(inout):: opts
        real*8 , intent( in) :: xpos(:)
        real*8 , intent( in) :: xnew(:)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: fnew(:,:,:)
        type (rcon_ends), intent(in) :: bclo(:)
        type (rcon_ends), intent(in) :: bcup(:)
        type (rmap_tics), &
        &   intent(inout) , optional :: tCPU

        real*8 , parameter :: RTOL = +1.d-14

    !------------------------------------------- variables !
        integer :: ipos
        real*8  :: diff,spac,same,xtol
        real*8  :: delx(1:1)
        logical :: uniform
        
#       ifdef __PPR_TIMER__
        integer(kind=8) :: ttic,ttoc,rate
#       endif

        if (ndof.lt.1) return
        if (npos.lt.2) return
        if (nnew.lt.2) return
        if (nvar.lt.1) return

    !------------- calc. grid-spacing and check uniformity !

        same = (xpos(npos)& 
             -  xpos(  +1)) / (npos-1)

        uniform = .true.
             
        xtol = same * RTOL

        do  ipos = +1 , npos-1, +1
        
            spac = xpos(ipos+1) &
        &        - xpos(ipos+0)
  
            diff = abs(spac - same)
        
            if (diff.gt.xtol) then
            
                uniform = .false.
            
            end if
  
            work% &
        &   cell_spac(ipos) = spac
        
        end do

       !uniform = .false.
    
    !------------- reconstruct FHAT over all cells in XPOS !

        if (.not. uniform) then

    !------------------------------------ variable spacing !
        call rcon1d(npos,nvar,ndof, &
        &           work%cell_spac, &
        &           fdat,bclo,bcup, &
        &           work%cell_func, &
        &           work,opts,tCPU)

        else
        
    !------------------------------------ constant spacing !
        delx(1) = work%cell_spac(1)
        
        call rcon1d(npos,nvar,ndof, &
        &           delx,    &
        &           fdat,bclo,bcup, &
        &           work%cell_func, &
        &           work,opts,tCPU)
        
        end if

    !------------- remap FDAT from XPOS to XNEW using FHAT !

        __TIC__

        select case(opts%cell_meth)
        case(pcm_method)
    !------------------------------------ 1st-order method !
        call p0m   (npos,nnew,nvar, &
        &           ndof,xpos,xnew, &
        &           work%cell_func, &
        &           fnew,xtol)
    
        case(plm_method)
    !------------------------------------ 2nd-order method !
        call p1m   (npos,nnew,nvar, &
        &           ndof,xpos,xnew, &
        &           work%cell_func, &
        &           fnew,xtol)

        case(ppm_method)
    !------------------------------------ 3rd-order method !
        call p2m   (npos,nnew,nvar, &
        &           ndof,xpos,xnew, &
        &           work%cell_func, &
        &           fnew,xtol)

        case(pqm_method)
    !------------------------------------ 5th-order method !
        call p4m   (npos,nnew,nvar, &
        &           ndof,xpos,xnew, &
        &           work%cell_func, &
        &           fnew,xtol)
        
        end select
        
        __TOC__(tCPU,rmap_time)
        
        return
    
    end  subroutine
    
    !------------ P0M : one-dimensional degree-0 remapping !
    
    pure subroutine p0m(npos,nnew,nvar,ndof, &
        &               xpos,xnew,fhat,fnew, &
        &               XTOL)

    !
    ! NPOS  no. edges in old grid.
    ! NNEW  no. edges in new grid.
    ! NVAR  no. discrete variables to remap.
    ! NDOF  no. degrees-of-freedom per cell.
    ! XPOS  old grid edge positions. XPOS is a length NPOS
    !       array .
    ! XNEW  new grid edge positions. XNEW is a length NNEW
    !       array .
    ! FHAT  reconstruction over old grid. FHAT has SIZE =
    !       MDOF-by-NVAR-by-NPOS-1 .
    ! FNEW  reconstruction over new grid. FNEW has SIZE = 
    !       NDOF-by-NVAR-by-NNEW-1 .
    ! XTOL  min. grid-cell thickness . 
    !

        implicit none    
    
    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nnew
        integer, intent( in) :: nvar,ndof
        real*8 , intent( in) :: xpos(:)
        real*8 , intent( in) :: xnew(:)
        real*8 , intent( in) :: fhat(:,:,:)
        real*8 , intent(out) :: fnew(:,:,:)
        real*8 , intent( in) :: XTOL
        
    !------------------------------------------- variables !
        integer :: kpos,jpos,ipos
        integer :: ivar,idof,pos0,pos1
        real*8  :: xmid,xhat,khat
        real*8  :: xxlo,xxhi,sslo,sshi
        real*8  :: intf,ivec(1:1)
        
    !------------- remap FDAT from XPOS to XNEW using FHAT !

        pos0 = +1 ; pos1 = +1

        do  kpos = +1, nnew-1

    !------ first cell in XPOS overlapping with XNEW(KPOS) !

            pos1 = max(pos1,1)

            do  pos0 = pos1, npos-1, +1
                if (xpos(pos0+1)&
            &       .gt. xnew(kpos+0)) then
                    exit
                end if
            end do

    !------ final cell in XPOS overlapping with XNEW(KPOS) !
    
            do  pos1 = pos0, npos-1, +1
                if (xpos(pos1+0)&
            &       .ge. xnew(kpos+1)) then
                    exit
                end if    
            end do
            
            pos1 = pos1 - 1

    !------------- integrate FHAT across overlapping cells !

            khat = xnew(kpos+1) &
            &    - xnew(kpos+0)
            khat = max (khat , XTOL)

            do  idof = +1,ndof
            do  ivar = +1,nvar
                
                fnew(idof,ivar,kpos) = 0.d0
            
            end do
            end do

            do  ipos = pos0 , pos1 , +1

    !------------------------------- integration endpoints !
    
                xxlo = max (xpos(ipos+0) , &
            &               xnew(kpos+0))
                xxhi = min (xpos(ipos+1) , &
            &               xnew(kpos+1))

    !------------------------------- local endpoint coords !
    
                xmid = xpos(ipos+1) * .5d0 &
            &        + xpos(ipos+0) * .5d0    
                xhat = xpos(ipos+1) * .5d0 &
            &        - xpos(ipos+0) * .5d0
     
                sslo = &
            &  (xxlo-xmid) / max(xhat,XTOL)
                sshi = &
            &  (xxhi-xmid) / max(xhat,XTOL)

    !------------------------------- integral basis vector !
    
                ivec(1) = sshi ** 1 / 1.d0 &
            &           - sslo ** 1 / 1.d0
        
    !--------- integrate FHAT across the overlap XXLO:XXHI !

                xhat =  xhat / khat
    
                do  ivar = 1 , nvar , +1
    
                intf =  dot_product (  &
            &   ivec,fhat(1:1,ivar,ipos-0))

                intf =  intf * xhat
        
    !--------- accumulate integral contributions from IPOS !
    
                fnew(  +1,ivar,kpos) = &
            &   fnew(  +1,ivar,kpos) + intf

                end do

            end do

        end do
        
        return
    
    end  subroutine

    !------------ P1M : one-dimensional degree-1 remapping !
    
    pure subroutine p1m(npos,nnew,nvar,ndof, &
        &               xpos,xnew,fhat,fnew, &
        &               XTOL)

    !
    ! NPOS  no. edges in old grid.
    ! NNEW  no. edges in new grid.
    ! NVAR  no. discrete variables to remap.
    ! NDOF  no. degrees-of-freedom per cell.
    ! XPOS  old grid edge positions. XPOS is a length NPOS
    !       array .
    ! XNEW  new grid edge positions. XNEW is a length NNEW
    !       array .
    ! FHAT  reconstruction over old grid. FHAT has SIZE =
    !       MDOF-by-NVAR-by-NPOS-1 .
    ! FNEW  reconstruction over new grid. FNEW has SIZE = 
    !       NDOF-by-NVAR-by-NNEW-1 .
    ! XTOL  min. grid-cell thickness . 
    !

        implicit none    
    
    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nnew
        integer, intent( in) :: nvar,ndof
        real*8 , intent( in) :: xpos(:)
        real*8 , intent( in) :: xnew(:)
        real*8 , intent( in) :: fhat(:,:,:)
        real*8 , intent(out) :: fnew(:,:,:)
        real*8 , intent( in) :: XTOL
        
    !------------------------------------------- variables !
        integer :: kpos,jpos,ipos
        integer :: ivar,idof,pos0,pos1
        real*8  :: xmid,xhat,khat
        real*8  :: xxlo,xxhi,sslo,sshi
        real*8  :: intf,ivec(1:2)
        
    !------------- remap FDAT from XPOS to XNEW using FHAT !

        pos0 = +1 ; pos1 = +1

        do  kpos = +1, nnew-1

    !------ first cell in XPOS overlapping with XNEW(KPOS) !

            pos1 = max(pos1,1)

            do  pos0 = pos1, npos-1, +1
                if (xpos(pos0+1)&
            &       .gt. xnew(kpos+0)) then
                    exit
                end if
            end do

    !------ final cell in XPOS overlapping with XNEW(KPOS) !
    
            do  pos1 = pos0, npos-1, +1
                if (xpos(pos1+0)&
            &       .ge. xnew(kpos+1)) then
                    exit
                end if    
            end do
            
            pos1 = pos1 - 1

    !------------- integrate FHAT across overlapping cells !

            khat = xnew(kpos+1) &
            &    - xnew(kpos+0)
            khat = max (khat , XTOL)

            do  idof = +1,ndof
            do  ivar = +1,nvar
                
                fnew(idof,ivar,kpos) = 0.d0
            
            end do
            end do

            do  ipos = pos0 , pos1 , +1

    !------------------------------- integration endpoints !
    
                xxlo = max (xpos(ipos+0) , &
            &               xnew(kpos+0))
                xxhi = min (xpos(ipos+1) , &
            &               xnew(kpos+1))

    !------------------------------- local endpoint coords !
    
                xmid = xpos(ipos+1) * .5d0 &
            &        + xpos(ipos+0) * .5d0    
                xhat = xpos(ipos+1) * .5d0 &
            &        - xpos(ipos+0) * .5d0
     
                sslo = &
            &  (xxlo-xmid) / max(xhat,XTOL)
                sshi = &
            &  (xxhi-xmid) / max(xhat,XTOL)

    !------------------------------- integral basis vector !
    
                ivec(1) = sshi ** 1 / 1.d0 &
            &           - sslo ** 1 / 1.d0
        
                ivec(2) = sshi ** 2 / 2.d0 &
            &           - sslo ** 2 / 2.d0
                
    !--------- integrate FHAT across the overlap XXLO:XXHI !

                xhat =  xhat / khat
    
                do  ivar = 1 , nvar , +1
    
                intf =  dot_product (  &
            &   ivec,fhat(1:2,ivar,ipos-0))

                intf =  intf * xhat
        
    !--------- accumulate integral contributions from IPOS !
    
                fnew(  +1,ivar,kpos) = &
            &   fnew(  +1,ivar,kpos) + intf

                end do

            end do

        end do
        
        return
    
    end  subroutine

    !------------ P2M : one-dimensional degree-2 remapping !
    
    pure subroutine p2m(npos,nnew,nvar,ndof, &
        &               xpos,xnew,fhat,fnew, &
        &               XTOL)

    !
    ! NPOS  no. edges in old grid.
    ! NNEW  no. edges in new grid.
    ! NVAR  no. discrete variables to remap.
    ! NDOF  no. degrees-of-freedom per cell.
    ! XPOS  old grid edge positions. XPOS is a length NPOS
    !       array .
    ! XNEW  new grid edge positions. XNEW is a length NNEW
    !       array .
    ! FHAT  reconstruction over old grid. FHAT has SIZE =
    !       MDOF-by-NVAR-by-NPOS-1 .
    ! FNEW  reconstruction over new grid. FNEW has SIZE = 
    !       NDOF-by-NVAR-by-NNEW-1 .
    ! XTOL  min. grid-cell thickness . 
    !

        implicit none    
    
    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nnew
        integer, intent( in) :: nvar,ndof
        real*8 , intent( in) :: xpos(:)
        real*8 , intent( in) :: xnew(:)
        real*8 , intent( in) :: fhat(:,:,:)
        real*8 , intent(out) :: fnew(:,:,:)
        real*8 , intent( in) :: XTOL
        
    !------------------------------------------- variables !
        integer :: kpos,jpos,ipos
        integer :: ivar,idof,pos0,pos1
        real*8  :: xmid,xhat,khat
        real*8  :: xxlo,xxhi,sslo,sshi
        real*8  :: intf,ivec(1:3)
        
    !------------- remap FDAT from XPOS to XNEW using FHAT !

        pos0 = +1 ; pos1 = +1

        do  kpos = +1, nnew-1

    !------ first cell in XPOS overlapping with XNEW(KPOS) !

            pos1 = max(pos1,1)

            do  pos0 = pos1, npos-1, +1
                if (xpos(pos0+1)&
            &       .gt. xnew(kpos+0)) then
                    exit
                end if
            end do

    !------ final cell in XPOS overlapping with XNEW(KPOS) !
    
            do  pos1 = pos0, npos-1, +1
                if (xpos(pos1+0)&
            &       .ge. xnew(kpos+1)) then
                    exit
                end if    
            end do
            
            pos1 = pos1 - 1

    !------------- integrate FHAT across overlapping cells !

            khat = xnew(kpos+1) &
            &    - xnew(kpos+0)
            khat = max (khat , XTOL)

            do  idof = +1,ndof
            do  ivar = +1,nvar
                
                fnew(idof,ivar,kpos) = 0.d0
            
            end do
            end do

            do  ipos = pos0 , pos1 , +1

    !------------------------------- integration endpoints !
    
                xxlo = max (xpos(ipos+0) , &
            &               xnew(kpos+0))
                xxhi = min (xpos(ipos+1) , &
            &               xnew(kpos+1))

    !------------------------------- local endpoint coords !
    
                xmid = xpos(ipos+1) * .5d0 &
            &        + xpos(ipos+0) * .5d0    
                xhat = xpos(ipos+1) * .5d0 &
            &        - xpos(ipos+0) * .5d0
     
                sslo = &
            &  (xxlo-xmid) / max(xhat,XTOL)
                sshi = &
            &  (xxhi-xmid) / max(xhat,XTOL)

    !------------------------------- integral basis vector !
    
                ivec(1) = sshi ** 1 / 1.d0 &
            &           - sslo ** 1 / 1.d0
        
                ivec(2) = sshi ** 2 / 2.d0 &
            &           - sslo ** 2 / 2.d0
            
                ivec(3) = sshi ** 3 / 3.d0 &
            &           - sslo ** 3 / 3.d0
                
    !--------- integrate FHAT across the overlap XXLO:XXHI !
    
                xhat =  xhat / khat

                do  ivar = 1 , nvar , +1
    
                intf =  dot_product (  &
            &   ivec,fhat(1:3,ivar,ipos-0))

                intf =  intf * xhat
        
    !--------- accumulate integral contributions from IPOS !
    
                fnew(  +1,ivar,kpos) = &
            &   fnew(  +1,ivar,kpos) + intf

                end do

            end do

        end do
        
        return
    
    end  subroutine

    !------------ P4M : one-dimensional degree-4 remapping !
    
    pure subroutine p4m(npos,nnew,nvar,ndof, &
        &               xpos,xnew,fhat,fnew, &
        &               XTOL)

    !
    ! NPOS  no. edges in old grid.
    ! NNEW  no. edges in new grid.
    ! NVAR  no. discrete variables to remap.
    ! NDOF  no. degrees-of-freedom per cell.
    ! XPOS  old grid edge positions. XPOS is a length NPOS
    !       array .
    ! XNEW  new grid edge positions. XNEW is a length NNEW
    !       array .
    ! FHAT  reconstruction over old grid. FHAT has SIZE =
    !       MDOF-by-NVAR-by-NPOS-1 .
    ! FNEW  reconstruction over new grid. FNEW has SIZE = 
    !       NDOF-by-NVAR-by-NNEW-1 .
    ! XTOL  min. grid-cell thickness . 
    !

        implicit none    
    
    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nnew
        integer, intent( in) :: nvar,ndof
        real*8 , intent( in) :: xpos(:)
        real*8 , intent( in) :: xnew(:)
        real*8 , intent( in) :: fhat(:,:,:)
        real*8 , intent(out) :: fnew(:,:,:)
        real*8 , intent( in) :: XTOL
        
    !------------------------------------------- variables !
        integer :: kpos,jpos,ipos
        integer :: ivar,idof,pos0,pos1
        real*8  :: xmid,xhat,khat
        real*8  :: xxlo,xxhi,sslo,sshi
        real*8  :: intf,ivec(1:5)
        
    !------------- remap FDAT from XPOS to XNEW using FHAT !

        pos0 = +1 ; pos1 = +1

        do  kpos = +1, nnew-1

    !------ first cell in XPOS overlapping with XNEW(KPOS) !

            pos1 = max(pos1,1)

            do  pos0 = pos1, npos-1, +1
                if (xpos(pos0+1)&
            &       .gt. xnew(kpos+0)) then
                    exit
                end if
            end do

    !------ final cell in XPOS overlapping with XNEW(KPOS) !
    
            do  pos1 = pos0, npos-1, +1
                if (xpos(pos1+0)&
            &       .ge. xnew(kpos+1)) then
                    exit
                end if    
            end do
            
            pos1 = pos1 - 1

    !------------- integrate FHAT across overlapping cells !

            khat = xnew(kpos+1) &
            &    - xnew(kpos+0)
            khat = max (khat , XTOL)

            do  idof = +1,ndof
            do  ivar = +1,nvar
                
                fnew(idof,ivar,kpos) = 0.d0
            
            end do
            end do

            do  ipos = pos0 , pos1 , +1

    !------------------------------- integration endpoints !
    
                xxlo = max (xpos(ipos+0) , &
            &               xnew(kpos+0))
                xxhi = min (xpos(ipos+1) , &
            &               xnew(kpos+1))

    !------------------------------- local endpoint coords !
    
                xmid = xpos(ipos+1) * .5d0 &
            &        + xpos(ipos+0) * .5d0    
                xhat = xpos(ipos+1) * .5d0 &
            &        - xpos(ipos+0) * .5d0
     
                sslo = &
            &  (xxlo-xmid) / max(xhat,XTOL)
                sshi = &
            &  (xxhi-xmid) / max(xhat,XTOL)

    !------------------------------- integral basis vector !
    
                ivec(1) = sshi ** 1 / 1.d0 &
            &           - sslo ** 1 / 1.d0
        
                ivec(2) = sshi ** 2 / 2.d0 &
            &           - sslo ** 2 / 2.d0
            
                ivec(3) = sshi ** 3 / 3.d0 &
            &           - sslo ** 3 / 3.d0
            
                ivec(4) = sshi ** 4 / 4.d0 &
            &           - sslo ** 4 / 4.d0
            
                ivec(5) = sshi ** 5 / 5.d0 &
            &           - sslo ** 5 / 5.d0
                
    !--------- integrate FHAT across the overlap XXLO:XXHI !
    
                xhat =  xhat / khat

                do  ivar = 1 , nvar , +1
    
                intf =  dot_product (  &
            &   ivec,fhat(1:5,ivar,ipos-0))

                intf =  intf * xhat
        
    !--------- accumulate integral contributions from IPOS !
    
                fnew(  +1,ivar,kpos) = &
            &   fnew(  +1,ivar,kpos) + intf

                end do

            end do

        end do
        
        return
    
    end  subroutine

    
    
