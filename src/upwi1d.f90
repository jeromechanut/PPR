
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
    ! UPWI1D.f90: upwind-biased flux-reconstruction scheme.
    !
    ! Darren Engwirda 
    ! 08-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !

    subroutine upwi1d(npos,nvar,ndof,spac,tDEL, &
        &             mask,uvel,qbar,qedg,bclo, &
        &             bchi,work,opts)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! SPAC  grid-cell spacing array. LENGTH(SPAC) == +1 if 
    !       spacing is uniform .
    ! TDEL  time-step .
    ! MASK  logical grid-cell masking array.
    ! UVEL  edge-centred velocity vectors. UVEL has SIZE = 
    !       NPOS-by-1 .
    ! QBAR  cell-centred integral moments. QBAR has SIZE =
    !       NDOF-by-NVAR-by-NPOS-1 .
    ! QEDG  edge-centred upwind-type eval. QEDG has SIZE = 
    !       NVAR-by-NPOS .
    ! BCLO  boundary condition at lower endpoint .
    ! BCHI  boundary condition at upper endpoint . 
    ! WORK  method work-space. See RCON-WORK for details .
    ! OPTS  method parameters. See RCON-OPTS for details .
    !
    
        implicit none
    
    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,nvar,ndof
        class(rmap_work), intent(inout):: work
        class(rmap_opts), intent(inout):: opts
        real*8 , intent(in)  :: spac(:)
        real*8 , intent(in)  :: tDEL
        logical, intent(in)  :: mask(:)
        real*8 , intent(in)  :: qbar(:,:,:)
        real*8 , intent(in)  :: uvel(:)
        real*8 , intent(out) :: qedg(:,:)
        class(rcon_ends), intent(in) :: bclo(:)
        class(rcon_ends), intent(in) :: bchi(:)
        
    !------------------------------------------- variables !      
        integer :: head,tail,nprt
    
        head = +0 ; tail = +0 ; qedg = 0.d+0
    
        do while (.true.)

    !--------------------------------- 1. find active part !
           
            do  head = tail+1, npos-1
            if (mask(head) .eqv..true.) exit
            end do
            
            do  tail = head+1, npos-1
            if (mask(tail).neqv..true.) exit
            end do
            tail = tail - 1

            if (head.ge.npos) exit
           
    !--------------------------------- 2. rcon active part !
            
            nprt = tail - head + 1
            
            if (size(spac).ne.+1) then
            
            call rcon1d(nprt+1,nvar,ndof , &
            &    spac(    head:tail), &
            &    qbar(:,:,head:tail), &
            &    bclo,bchi,work%cell_func, &
            &    work,opts )
            
            else
            
            call rcon1d(nprt+1,nvar,ndof , &
            &    spac,qbar(:,:,head:tail), &
            &    bclo,bchi,work%cell_func, &
            &    work,opts )
            
            end if
            
    !--------------------------------- 3. int. active part !

            select case(opts%cell_meth)
                case(pcm_method)       !! 1st-order scheme
   
                if (size(spac).ne.+1) then
                
                call p0u(nprt+1,nvar , &
                &    spac(  head:tail+0) , &
                &    tDEL, &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                else
                
                call p0u(nprt+1,nvar , &
                &    spac,tDEL , &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                end if
                
                case(plm_method)       !! 2nd-order scheme
                
                if (size(spac).ne.+1) then
                
                call p1u(nprt+1,nvar , &
                &    spac(  head:tail+0) , &
                &    tDEL, &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                else
                
                call p1u(nprt+1,nvar , &
                &    spac,tDEL , &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                end if
                
                case(ppm_method)       !! 3rd-order scheme
                
                if (size(spac).ne.+1) then
                
                call p2u(nprt+1,nvar , &
                &    spac(  head:tail+0) , &
                &    tDEL, &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                else
                
                call p2u(nprt+1,nvar , &
                &    spac,tDEL , &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                end if
                
                case(pqm_method)       !! 5th-order scheme
  
                if (size(spac).ne.+1) then
                
                call p4u(nprt+1,nvar , &
                &    spac(  head:tail+0) , &
                &    tDEL, &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                else
                
                call p4u(nprt+1,nvar , &
                &    spac,tDEL , &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                end if
            
            end select

        end do    
  
        return
        
    end  subroutine
    
    ! P0U: upwind-type flux reconstruction, degree-0 impl. !

    pure subroutine p0u(npos,nvar,SPAC,tDEL, &
        &               uvel,QHAT,qedg)
    
    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! SPAC  grid spacing vector. SIZE(SPAC)==+1 if uniform .
    ! TDEL  time-step .
    ! UVEL  edge-centred velocity vectors. UVEL has SIZE = 
    !       NPOS-by-1 .
    ! QHAT  cell-centred polynomial recon. QHAT has SIZE =  
    !       NDOF-by-NVAR-by-NPOS-1 .
    ! QEDG  edge-centred upwind-type eval. QEDG has SIZE = 
    !       NVAR-by-NPOS .
    !

        implicit none
        
    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,nvar
        real*8 , intent(in)  :: SPAC(:)
        real*8 , intent(in)  :: tDEL
        real*8 , intent(in)  :: uvel(:)
        real*8 , intent(in)  :: QHAT(:,:,:)
        real*8 , intent(out) :: qedg(:,:)  
  
    !------------------------------------------- variables !      
        integer :: ipos,ivar
        real*8  :: uCFL,SS11,SS22,QQ11,QQ22
        real*8  :: Qdel,Sdel,Qval,Smag
        real*8  :: BV11(1:1),BV22(1:1)
        
        real*8 , parameter :: ZERO = 1.d-14

    !----------- single-cell, lagrangian-type upwind rcon. !
        
        do  ipos = +2 , npos - 1
        
            if (size(SPAC).ne.+1) then          
            uCFL = uvel(ipos) &
        &        * tDEL / SPAC(ipos-1)
            else
            uCFL = uvel(ipos) &
        &        * tDEL / SPAC(    +1)
            end if
        
            if (uCFL .gt. +0.0d0) then
            
    !----------- integrate profile over upwind cell IPOS-1 !
  
            SS11 = +1.d0 - 2.d0 * uCFL
            SS22 = +1.d0

            BV11(1) = SS11 ** 1 / 1.d0
            
            BV22(1) = SS22 ** 1 / 1.d0
            
            do  ivar = 1, nvar
            
                QQ11 = dot_product( &
        &           BV11, &
        &       QHAT(1:1,ivar,ipos-1))
        
                QQ22 = dot_product( &
        &           BV22, &
        &       QHAT(1:1,ivar,ipos-1))
                
                Sdel = SS22 - SS11
                Smag = abs(Sdel)
                Qdel = QQ22 - QQ11
                
                if (Smag.gt.ZERO) then
                    Qval = Qdel & 
        &                / Sdel
                else
                    Qval = QQ22 &
        &                / SS22
                end if
                
                qedg(ivar,ipos) = Qval
            
            end do
                      
            else &
        &   if (uCFL .lt. -0.0d0) then
            
    !----------- integrate profile over upwind cell IPOS+0 !
            
            SS11 = -1.d0 - 2.d0 * uCFL
            SS22 = -1.d0
        
            BV11(1) = SS11 ** 1 / 1.d0
            
            BV22(1) = SS22 ** 1 / 1.d0
            
            do  ivar = 1, nvar
                
                QQ11 = dot_product( &
        &           BV11, &
        &       QHAT(1:1,ivar,ipos-0))
        
                QQ22 = dot_product( &
        &           BV22, &
        &       QHAT(1:1,ivar,ipos-0))
                
                Sdel = SS22 - SS11
                Smag = abs(Sdel)
                Qdel = QQ22 - QQ11
                
                if (Smag.gt.ZERO) then
                    Qval = Qdel & 
        &                / Sdel
                else
                    Qval = QQ22 &
        &                / SS22
                end if
                
                qedg(ivar,ipos) = Qval
                
            end do
            
            end if
        
        end do
               
        return
    
    end  subroutine

    ! P1U: upwind-type flux reconstruction, degree-1 impl. !
    
    pure subroutine p1u(npos,nvar,SPAC,tDEL, &
        &               uvel,QHAT,qedg)
    
    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! SPAC  grid spacing vector. SIZE(SPAC)==+1 if uniform .
    ! TDEL  time-step .
    ! UVEL  edge-centred velocity vectors. UVEL has SIZE = 
    !       NPOS-by-1 .
    ! QHAT  cell-centred polynomial recon. QHAT has SIZE =  
    !       NDOF-by-NVAR-by-NPOS-1 .
    ! QEDG  edge-centred upwind-type eval. QEDG has SIZE = 
    !       NVAR-by-NPOS .
    !

        implicit none
        
    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,nvar
        real*8 , intent(in)  :: SPAC(:)
        real*8 , intent(in)  :: tDEL
        real*8 , intent(in)  :: uvel(:)
        real*8 , intent(in)  :: QHAT(:,:,:)
        real*8 , intent(out) :: qedg(:,:)  
  
    !------------------------------------------- variables !      
        integer :: ipos,ivar
        real*8  :: uCFL,SS11,SS22,QQ11,QQ22
        real*8  :: Qdel,Sdel,Qval,Smag
        real*8  :: BV11(1:2),BV22(1:2)
        
        real*8 , parameter :: ZERO = 1.d-14

    !----------- single-cell, lagrangian-type upwind rcon. !
        
        do  ipos = +2 , npos - 1
        
            if (size(SPAC).ne.+1) then          
            uCFL = uvel(ipos) &
        &        * tDEL / SPAC(ipos-1)
            else
            uCFL = uvel(ipos) &
        &        * tDEL / SPAC(    +1)
            end if
        
            if (uCFL .gt. +0.0d0) then
            
    !----------- integrate profile over upwind cell IPOS-1 !
  
            SS11 = +1.d0 - 2.d0 * uCFL
            SS22 = +1.d0

            BV11(1) = SS11 ** 1 / 1.d0
            BV11(2) = SS11 ** 2 / 2.d0
            
            BV22(1) = SS22 ** 1 / 1.d0
            BV22(2) = SS22 ** 2 / 2.d0
            
            do  ivar = 1, nvar
            
                QQ11 = dot_product( &
        &           BV11, &
        &       QHAT(1:2,ivar,ipos-1))
        
                QQ22 = dot_product( &
        &           BV22, &
        &       QHAT(1:2,ivar,ipos-1))
                
                Sdel = SS22 - SS11
                Smag = abs(Sdel)
                Qdel = QQ22 - QQ11
                
                if (Smag.gt.ZERO) then
                    Qval = Qdel & 
        &                / Sdel
                else
                    Qval = QQ22 &
        &                / SS22
                end if
                
                qedg(ivar,ipos) = Qval
            
            end do
                      
            else &
        &   if (uCFL .lt. -0.0d0) then
            
    !----------- integrate profile over upwind cell IPOS+0 !
            
            SS11 = -1.d0 - 2.d0 * uCFL
            SS22 = -1.d0
        
            BV11(1) = SS11 ** 1 / 1.d0
            BV11(2) = SS11 ** 2 / 2.d0
            
            BV22(1) = SS22 ** 1 / 1.d0
            BV22(2) = SS22 ** 2 / 2.d0
            
            do  ivar = 1, nvar
                
                QQ11 = dot_product( &
        &           BV11, &
        &       QHAT(1:2,ivar,ipos-0))
        
                QQ22 = dot_product( &
        &           BV22, &
        &       QHAT(1:2,ivar,ipos-0))
                
                Sdel = SS22 - SS11
                Smag = abs(Sdel)
                Qdel = QQ22 - QQ11
                
                if (Smag.gt.ZERO) then
                    Qval = Qdel & 
        &                / Sdel
                else
                    Qval = QQ22 &
        &                / SS22
                end if
                
                qedg(ivar,ipos) = Qval
                
            end do
            
            end if
        
        end do
               
        return
    
    end  subroutine

    ! P2U: upwind-type flux reconstruction, degree-2 impl. !

    pure subroutine p2u(npos,nvar,SPAC,tDEL, &
        &               uvel,QHAT,qedg)
    
    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! SPAC  grid spacing vector. SIZE(SPAC)==+1 if uniform .
    ! TDEL  time-step .
    ! UVEL  edge-centred velocity vectors. UVEL has SIZE = 
    !       NPOS-by-1 .
    ! QHAT  cell-centred polynomial recon. QHAT has SIZE =  
    !       NDOF-by-NVAR-by-NPOS-1 .
    ! QEDG  edge-centred upwind-type eval. QEDG has SIZE = 
    !       NVAR-by-NPOS .
    !

        implicit none
        
    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,nvar
        real*8 , intent(in)  :: SPAC(:)
        real*8 , intent(in)  :: tDEL
        real*8 , intent(in)  :: uvel(:)
        real*8 , intent(in)  :: QHAT(:,:,:)
        real*8 , intent(out) :: qedg(:,:)  
  
    !------------------------------------------- variables !      
        integer :: ipos,ivar
        real*8  :: uCFL,SS11,SS22,QQ11,QQ22
        real*8  :: Qdel,Sdel,Qval,Smag
        real*8  :: BV11(1:3),BV22(1:3)
        
        real*8 , parameter :: ZERO = 1.d-14

    !----------- single-cell, lagrangian-type upwind rcon. !
        
        do  ipos = +2 , npos - 1
        
            if (size(SPAC).ne.+1) then          
            uCFL = uvel(ipos) &
        &        * tDEL / SPAC(ipos-1)
            else
            uCFL = uvel(ipos) &
        &        * tDEL / SPAC(    +1)
            end if
        
            if (uCFL .gt. +0.0d0) then
            
    !----------- integrate profile over upwind cell IPOS-1 !
  
            SS11 = +1.d0 - 2.d0 * uCFL
            SS22 = +1.d0

            BV11(1) = SS11 ** 1 / 1.d0
            BV11(2) = SS11 ** 2 / 2.d0
            BV11(3) = SS11 ** 3 / 3.d0
            
            BV22(1) = SS22 ** 1 / 1.d0
            BV22(2) = SS22 ** 2 / 2.d0
            BV22(3) = SS22 ** 3 / 3.d0
            
            do  ivar = 1, nvar
            
                QQ11 = dot_product( &
        &           BV11, &
        &       QHAT(1:3,ivar,ipos-1))
        
                QQ22 = dot_product( &
        &           BV22, &
        &       QHAT(1:3,ivar,ipos-1))
                
                Sdel = SS22 - SS11
                Smag = abs(Sdel)
                Qdel = QQ22 - QQ11
                
                if (Smag.gt.ZERO) then
                    Qval = Qdel & 
        &                / Sdel
                else
                    Qval = QQ22 &
        &                / SS22
                end if
                
                qedg(ivar,ipos) = Qval
            
            end do
                      
            else &
        &   if (uCFL .lt. -0.0d0) then
            
    !----------- integrate profile over upwind cell IPOS+0 !
            
            SS11 = -1.d0 - 2.d0 * uCFL
            SS22 = -1.d0
        
            BV11(1) = SS11 ** 1 / 1.d0
            BV11(2) = SS11 ** 2 / 2.d0
            BV11(3) = SS11 ** 3 / 3.d0
            
            BV22(1) = SS22 ** 1 / 1.d0
            BV22(2) = SS22 ** 2 / 2.d0
            BV22(3) = SS22 ** 3 / 3.d0
            
            do  ivar = 1, nvar
                
                QQ11 = dot_product( &
        &           BV11, &
        &       QHAT(1:3,ivar,ipos-0))
        
                QQ22 = dot_product( &
        &           BV22, &
        &       QHAT(1:3,ivar,ipos-0))
                
                Sdel = SS22 - SS11
                Smag = abs(Sdel)
                Qdel = QQ22 - QQ11
                
                if (Smag.gt.ZERO) then
                    Qval = Qdel & 
        &                / Sdel
                else
                    Qval = QQ22 &
        &                / SS22
                end if
                
                qedg(ivar,ipos) = Qval
                
            end do
            
            end if
        
        end do
               
        return
    
    end  subroutine

    ! P4U: upwind-type flux reconstruction, degree-4 impl. !

    pure subroutine p4u(npos,nvar,SPAC,tDEL, &
        &               uvel,QHAT,qedg)
    
    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! SPAC  grid spacing vector. SIZE(SPAC)==+1 if uniform .
    ! TDEL  time-step .
    ! UVEL  edge-centred velocity vectors. UVEL has SIZE = 
    !       NPOS-by-1 .
    ! QHAT  cell-centred polynomial recon. QHAT has SIZE =  
    !       NDOF-by-NVAR-by-NPOS-1 .
    ! QEDG  edge-centred upwind-type eval. QEDG has SIZE = 
    !       NVAR-by-NPOS .
    !

        implicit none
        
    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,nvar
        real*8 , intent(in)  :: SPAC(:)
        real*8 , intent(in)  :: tDEL
        real*8 , intent(in)  :: uvel(:)
        real*8 , intent(in)  :: QHAT(:,:,:)
        real*8 , intent(out) :: qedg(:,:)  
  
    !------------------------------------------- variables !      
        integer :: ipos,ivar
        real*8  :: uCFL,SS11,SS22,QQ11,QQ22
        real*8  :: Qdel,Sdel,Qval,Smag
        real*8  :: BV11(1:5),BV22(1:5)
        
        real*8 , parameter :: ZERO = 1.d-14

    !----------- single-cell, lagrangian-type upwind rcon. !
        
        do  ipos = +2 , npos - 1
        
            if (size(SPAC).ne.+1) then          
            uCFL = uvel(ipos) &
        &        * tDEL / SPAC(ipos-1)
            else
            uCFL = uvel(ipos) &
        &        * tDEL / SPAC(    +1)
            end if
        
            if (uCFL .gt. +0.0d0) then
            
    !----------- integrate profile over upwind cell IPOS-1 !
  
            SS11 = +1.d0 - 2.d0 * uCFL
            SS22 = +1.d0

            BV11(1) = SS11 ** 1 / 1.d0
            BV11(2) = SS11 ** 2 / 2.d0
            BV11(3) = SS11 ** 3 / 3.d0
            BV11(4) = SS11 ** 4 / 4.d0
            BV11(5) = SS11 ** 5 / 5.d0
            
            BV22(1) = SS22 ** 1 / 1.d0
            BV22(2) = SS22 ** 2 / 2.d0
            BV22(3) = SS22 ** 3 / 3.d0
            BV22(4) = SS22 ** 4 / 4.d0
            BV22(5) = SS22 ** 5 / 5.d0
            
            do  ivar = 1, nvar
            
                QQ11 = dot_product( &
        &           BV11, &
        &       QHAT(1:5,ivar,ipos-1))
        
                QQ22 = dot_product( &
        &           BV22, &
        &       QHAT(1:5,ivar,ipos-1))
                
                Sdel = SS22 - SS11
                Smag = abs(Sdel)
                Qdel = QQ22 - QQ11
                
                if (Smag.gt.ZERO) then
                    Qval = Qdel & 
        &                / Sdel
                else
                    Qval = QQ22 &
        &                / SS22
                end if
                
                qedg(ivar,ipos) = Qval
            
            end do
                      
            else &
        &   if (uCFL .lt. -0.0d0) then
            
    !----------- integrate profile over upwind cell IPOS+0 !
            
            SS11 = -1.d0 - 2.d0 * uCFL
            SS22 = -1.d0
        
            BV11(1) = SS11 ** 1 / 1.d0
            BV11(2) = SS11 ** 2 / 2.d0
            BV11(3) = SS11 ** 3 / 3.d0
            BV11(4) = SS11 ** 4 / 4.d0
            BV11(5) = SS11 ** 5 / 5.d0
            
            BV22(1) = SS22 ** 1 / 1.d0
            BV22(2) = SS22 ** 2 / 2.d0
            BV22(3) = SS22 ** 3 / 3.d0
            BV22(4) = SS22 ** 4 / 4.d0
            BV22(5) = SS22 ** 5 / 5.d0
            
            do  ivar = 1, nvar
                
                QQ11 = dot_product( &
        &           BV11, &
        &       QHAT(1:5,ivar,ipos-0))
        
                QQ22 = dot_product( &
        &           BV22, &
        &       QHAT(1:5,ivar,ipos-0))
                
                Sdel = SS22 - SS11
                Smag = abs(Sdel)
                Qdel = QQ22 - QQ11
                
                if (Smag.gt.ZERO) then
                    Qval = Qdel & 
        &                / Sdel
                else
                    Qval = QQ22 &
        &                / SS22
                end if
                
                qedg(ivar,ipos) = Qval
                
            end do
            
            end if
        
        end do
               
        return
    
    end  subroutine



