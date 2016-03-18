
    !
    ! ADVCFN.f90: finite-volume type scalar advection.
    !
    ! Darren Engwirda 
    ! 17-Mar-2016
    ! engwirda [at] mit [dot] edu
    !
    
    module advcfn

    use polyfn
    use rmapfn

    implicit none
    
    public :: advc1d
    
    contains

    !----------------------------------------------------------------!
    ! ADVC1D: 1-dimensional scalar advection.                        !
    !----------------------------------------------------------------!
    subroutine advc1d(npos,nvar,ndof,delx,delt, &
        &             uvel,uflx,qdat,qflx, &
        &             bclo,bchi,work,opts  )

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,nvar,ndof
        class(rmap_work), intent(inout):: work
        class(rmap_opts), intent(inout):: opts
        real*8 , intent(in)  :: delx(:)
        real*8 , intent(in)  :: delt
        real*8 , intent(in)  :: uvel(:)
        real*8 , intent(in)  :: uflx(:)
        real*8 , intent(in)  :: qdat(:,:,:)
        real*8 , intent(out) :: qflx(:,:)
        type (rcon_ends), intent(in) :: bclo(:)
        type (rcon_ends), intent(in) :: bchi(:)

    !------------------------------------------- variables !
        integer :: ipos,ivar,mdof
        integer :: head,tail
        real*8  :: uCFL
        real*8  :: SS11,SS22
        real*8  :: QQ11,QQ22
        real*8  :: BV11(max_poly_ndof)
        real*8  :: BV22(max_poly_ndof)
        real*8  :: Qdel,Sdel,Qval

        if (npos.lt.2) return
        if (nvar.lt.1) return

        head = +2 ; tail = npos - 1

    !------------- reconstruct QHAT over all cells in XPOS !

        call rcon1d(npos,nvar,ndof, &
        &           delx,qdat, &
        &           bclo,bchi, &
        &           work%cell_func, &
        &           work,opts)

        mdof = ndof1d(opts%cell_meth)
    
    !------------- compute transport for each edge in XPOS !
    
        do  ivar = 1, nvar
                
            qflx(ivar,head-1) = 0.d0
            qflx(ivar,tail+1) = 0.d0
        
        end do
    
        do  ipos = head, tail
        
            if (uvel(ipos) .gt. 0.d0) then
            
    !----------- integrate profile over upwind cell IPOS-1 !
            
                uCFL = uvel(ipos) &
        &            * delt / delx(ipos-1)
            
                SS11 = +1.d0 - 2.d0 * uCFL
                SS22 = +1.d0

                call poly_basis(-1,mdof,SS11,BV11)
                call poly_basis(-1,mdof,SS22,BV22)
                
                do  ivar = 1, nvar
                
                QQ11 = poly_value(-1,mdof, &
        &           work%cell_func(:,ivar,ipos-1), &
        &               SS11,BV11)
            
                QQ22 = poly_value(-1,mdof, &
        &           work%cell_func(:,ivar,ipos-1), &
        &               SS22,BV22)
                
                Qdel = QQ22 - QQ11
                Sdel = SS22 - SS11
                
                if (abs(Qdel).gt.epsilon(Sdel)) then
                    Qval = Qdel & 
        &                / Sdel
                else
                    Qval = QQ22
                end if
                
                qflx(ivar,ipos) = +uflx(ipos) * Qval
                
                end do
                      
            else
            
    !----------- integrate profile over upwind cell IPOS+0 !
            
                uCFL = uvel(ipos) &
        &            * delt / delx(ipos+0)
                     
                SS11 = -1.d0 - 2.d0 * uCFL
                SS22 = -1.d0
            
                call poly_basis(-1,mdof,SS11,BV11)
                call poly_basis(-1,mdof,SS22,BV22)
                
                do  ivar = 1, nvar
                
                QQ11 = poly_value(-1,mdof, &
        &           work%cell_func(:,ivar,ipos+0), &
        &               SS11,BV11)
            
                QQ22 = poly_value(-1,mdof, &
        &           work%cell_func(:,ivar,ipos+0), &
        &               SS22,BV22)
                
                Qdel = QQ22 - QQ11
                Sdel = SS22 - SS11
                
                if (abs(Qdel).gt.epsilon(Sdel)) then
                    Qval = Qdel & 
        &                / Sdel
                else
                    Qval = QQ22
                end if
                
                qflx(ivar,ipos) = +uflx(ipos) * Qval
                
                end do
            
            end if
        
        end do
  
        return
    
    end subroutine advc1d
    
    end module advcfn 



