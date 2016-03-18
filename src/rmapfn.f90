
    !
    ! RMAPFN.f90: piecewise polynomial reconstructions and 
    ! conservative integral re-mapping.
    !
    ! Darren Engwirda 
    ! 17-Mar-2016
    ! engwirda [at] mit [dot] edu
    !

    module rmapfn

    use matrix
    use polyfn

    implicit none

    !------------------------------------ method selection !
    
    integer, parameter :: p1e_method = +100
    integer, parameter :: p3e_method = +101
    integer, parameter :: p5e_method = +102

    integer, parameter :: pcm_method = +200
    integer, parameter :: plm_method = +201
    integer, parameter :: ppm_method = +202
    integer, parameter :: pqm_method = +203

    integer, parameter :: null_limit = +300
    integer, parameter :: mono_limit = +301
    integer, parameter :: weno_limit = +302

    integer, parameter :: bcon_loose = +400
    integer, parameter :: bcon_value = +401
    integer, parameter :: bcon_slope = +402
   
    public :: ndof1d
    public :: rcon1d
    public :: rmap1d

    type rcon_opts
    !------------------------------- parameters for RCON1D !
        integer :: edge_meth
        integer :: cell_meth
        integer :: cell_lims
    end type rcon_opts

    type rcon_ends
    !------------------------------- end-conditions struct !
        integer :: bcopt
        real*8  :: value
        real*8  :: slope
    end type rcon_ends

    type rcon_work
    !------------------------------- work-space for RCON1D !
        real*8, allocatable  :: edge_func(:,:,:)
        real*8, allocatable  :: edge_dfdx(:,:,:)
        real*8, allocatable  :: cell_oscl(:,:,:)
    contains
        procedure :: init => init_rcon_work
        procedure :: free => free_rcon_work
    end type rcon_work

    type, extends(rcon_opts) :: rmap_opts
    !------------------------------- parameters for RMAP1D !
    end type rmap_opts

    type, extends(rcon_work) :: rmap_work
    !------------------------------- work-space for RMAP1D !
        real*8, allocatable  :: cell_spac(:)
        real*8, allocatable  :: cell_func(:,:,:)
    contains
        procedure :: init => init_rmap_work
        procedure :: free => free_rmap_work
    end type rmap_work

    contains

    !----------------------------------------------------------------!
    ! INIT_RCON: init work-space for RCON1D.                         !
    !----------------------------------------------------------------!
    !   NPOS - no. edges in grid.                                    !
    !   NVAR - no. variables to reconstruct.                         !
    !   OPTS - method options, see RMAP_OPTS for additional details. !
    !----------------------------------------------------------------!
    subroutine init_rcon_work(this,npos,nvar,opts)

        implicit none

    !------------------------------------------- arguments !
        class(rcon_work) , intent(inout) :: this
        integer, intent(in):: npos
        integer, intent(in):: nvar
        class(rcon_opts) , optional      :: opts

    !------------------------------------------- variables !
        integer :: okay

        allocate(this%edge_func(2,nvar,npos), &
        &        this%edge_dfdx(2,nvar,npos), &
        &        this%cell_oscl(2,nvar,npos), &
        &        stat=okay)
        
    end subroutine init_rcon_work

    !----------------------------------------------------------------!
    ! FREE_RCON: free work-space for RCON1D.                         !
    !----------------------------------------------------------------!
    subroutine free_rcon_work(this)

        implicit none

    !------------------------------------------- arguments !
        class(rcon_work), intent(inout) :: this

    !------------------------------------------- variables !
        integer :: okay

        deallocate(this%edge_func, &
        &          this%edge_dfdx, &
        &          this%cell_oscl, &
        &          stat=okay)
        
    end subroutine free_rcon_work

    !----------------------------------------------------------------!
    ! INIT_RMAP: init work-space for RMAP1D.                         !
    !----------------------------------------------------------------!
    !   NPOS - no. edges in grid.                                    !
    !   NVAR - no. variables to reconstruct.                         !
    !   OPTS - method options, see RMAP_OPTS for additional details. !
    !----------------------------------------------------------------!
    subroutine init_rmap_work(this,npos,nvar,opts)

        implicit none

    !------------------------------------------- arguments !
        class(rmap_work) , intent(inout) :: this
        integer, intent(in) :: npos
        integer, intent(in) :: nvar
        class(rcon_opts) , optional      :: opts

    !------------------------------------------- variables !
        integer :: okay,ndof

        ndof = ndof1d(opts%cell_meth)

        allocate(this%edge_func(2,nvar,npos), &
        &        this%edge_dfdx(2,nvar,npos), &
        &        this%cell_oscl(2,nvar,npos), &
        &        this%cell_spac(       npos), &
        &        this%cell_func(ndof,nvar,npos), &
        &        stat=okay)
        
    end subroutine init_rmap_work

    !----------------------------------------------------------------!
    ! FREE_RMAP: free work-space for RMAP1D.                         !
    !----------------------------------------------------------------!
    subroutine free_rmap_work(this)

        implicit none

    !------------------------------------------- arguments !
        class(rmap_work), intent(inout) :: this

    !------------------------------------------- variables !
        integer :: okay

        deallocate(this%edge_func, &
        &          this%edge_dfdx, &
        &          this%cell_oscl, &
        &          this%cell_func, &
        &          this%cell_spac, &
        &          stat=okay)

    end subroutine free_rmap_work

    !----------------------------- include implementations !

        include 'pbc.f90'

        include 'p1e.f90'
        include 'p3e.f90'
        include 'p5e.f90'

        include 'weno.f90'

        include 'pcm.f90'
        include 'plm.f90'
        include 'ppm.f90'
        include 'pqm.f90'


    !----------------------------------------------------------------!
    ! NDOF1D : no. degrees-of-freedom per polynomial.                !
    !----------------------------------------------------------------!
    pure function ndof1d(meth) result(rdof)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: meth

    !------------------------------------------- variables !
        integer  :: rdof

        select case(meth)
    !-------------------------------- edge reconstructions !
        case (p1e_method)
            rdof = +2
        case (p3e_method)
            rdof = +4
        case (p5e_method)
            rdof = +6
    !-------------------------------- cell reconstructions !
        case (pcm_method)
            rdof = +1
        case (plm_method)
            rdof = +2
        case (ppm_method)
            rdof = +3
        case (pqm_method)
            rdof = +5
        
        case default 
            rdof = +0

        end select
        
    end function ndof1d

    !----------------------------------------------------------------!
    ! RCON1D : one-dimensional integral reconstruction.              !
    !----------------------------------------------------------------!
    !   DELX - NPOS-1-by-1 array of grid coordinates. Grid spacing   !
    !          can be non-uniform.                                   !
    !   FDAT - NDOF-by-NVAR-by-NPOS-1 array of cell-wise moments for !
    !          NVAR discrete variables.                              !
    !   BCLO - NVAR-by-1 array of boundary conditions for NVAR disc- !
    !          rete variables at the lower endpoint.                 !
    !   BCUP - NVAR-by-1 array of boundary conditions for NVAR disc- !
    !          rete variables at the upper endpoint.                 !
    !   FHAT - MDOF-by-NVAR-by-NPOS-1 array of polynomial coefficie- !
    !          nts for each grid-cell. See POLY_VALUE, etc for addi- !
    !          tional details.                                       !
    !   WORK - work-space. See RCON_WORK for additional details.     !
    !   OPTS - parameters. See RCON_OPTS for additional details.     !
    !----------------------------------------------------------------!
    subroutine rcon1d(npos,nvar,ndof,delx,fdat, &
        &             bclo,bcup,fhat,work,opts)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        class(rcon_work), intent(inout):: work
        class(rcon_opts), intent(in)   :: opts
        real*8 , intent( in) :: delx(:)
        real*8 , intent(out) :: fhat(:,:,:)
        real*8 , intent( in), contiguous :: fdat(:,:,:)
        type (rcon_ends), intent(in) :: bclo(:)
        type (rcon_ends), intent(in) :: bcup(:)

    !------------------------------------------- variables !
        integer :: halo

        if (ndof.lt.1) return
        if (npos.lt.2) return
        if (nvar.lt.1) return
     
    !-------------------------- compute edge values/slopes !

        if ( (opts%cell_meth.eq.ppm_method) &
    &  .or.  (opts%cell_meth.eq.pqm_method) ) then

        select case (opts%edge_meth)
            case(p1e_method)
    !------------------------------------ 2nd-order method !
            halo = +1
            call p1e(npos,nvar,ndof, &
        &            delx,fdat,      &
        &            work%edge_func, &
        &            work%edge_dfdx)

            case(p3e_method)
    !------------------------------------ 4th-order method !           
            halo = +2
            call p3e(npos,nvar,ndof, &
        &            delx,fdat,      &
        &            bclo,bcup,      &
        &            work%edge_func, &
        &            work%edge_dfdx, &
        &            opts)

            case(p5e_method)
    !------------------------------------ 6th-order method !           
            halo = +3
            call p5e(npos,nvar,ndof, &
        &            delx,fdat,      &
        &            bclo,bcup,      &
        &            work%edge_func, &
        &            work%edge_dfdx, &
        &            opts)

        end select

        end if

    !-------------------------- compute grid-cell profiles !

        select case (opts%cell_meth)
            case(pcm_method)
    !------------------------------------ 1st-order method !
            call pcm(npos,nvar,ndof, &
        &            fdat,fhat)

            case(plm_method)
    !------------------------------------ 2nd-order method !
            call plm(npos,nvar,ndof, &
        &            delx,fdat,fhat, &
        &            opts%cell_lims)

            case(ppm_method)
    !------------------------------------ 3rd-order method !
            call ppm(npos,nvar,ndof, &
        &            delx,fdat,fhat, &
        &            work%edge_func, &
        &            work%cell_oscl, &
        &            opts%cell_lims, &
        &            halo)

            case(pqm_method)
    !------------------------------------ 5th-order method !
            call pqm(npos,nvar,ndof, &
        &            delx,fdat,fhat, &
        &            work%edge_func, &
        &            work%edge_dfdx, &
        &            work%cell_oscl, &
        &            opts%cell_lims, &
        &            halo)

        end select

    end subroutine rcon1d

    !----------------------------------------------------------------!
    ! RMAP1D : one-dimensional integral remapping.                   !
    !----------------------------------------------------------------!
    !   XPOS - NPOS-by-1 array of "old" grid coordinates. Grid spac- !
    !          ing can be non-uniform.                               !
    !   XNEW - NNEW-by-1 array of "new" grid coordinates. Grid spac- !
    !          ing can be non-uniform.                               !
    !   FDAT - NDOF-by-NVAR-by-NPOS-1 array of old cell-wise moments !
    !          for NVAR discrete variables.                          !
    !   FNEW - NDOF-by-NVAR-by-NPOS-1 array of new cell-wise moments !
    !          for NVAR discrete variables.                          !
    !   BCLO - NVAR-by-1 array of boundary conditions for NVAR disc- !
    !          rete variables at the lower endpoint.                 !
    !   BCUP - NVAR-by-1 array of boundary conditions for NVAR disc- !
    !          rete variables at the upper endpoint.                 !
    !   WORK - work-space. See RMAP_WORK for additional details.     !
    !   OPTS - parameters. See RMAP_OPTS for additional details.     !
    !----------------------------------------------------------------!
    subroutine rmap1d(npos,nnew,nvar,ndof,xpos, &
        &             xnew,fdat,fnew,bclo,bcup, &
        &             work,opts)

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

    !------------------------------------------- variables !
        integer :: kpos,ipos,ivar
        integer :: mdof
        integer :: pos0,pos1
        real*8  :: diff,spac
        real*8  :: same,xtol
        real*8  :: delx(1:1)
        real*8  :: xmid,xhat,khat
        real*8  :: xxlo,xxhi
        real*8  :: sslo,sshi
        real*8  :: fnlo,fnhi
        real*8  :: bvlo(max_poly_ndof)
        real*8  :: bvhi(max_poly_ndof)
        logical :: uniform

        if (ndof.lt.1) return
        if (npos.lt.2) return
        if (nnew.lt.2) return
        if (nvar.lt.1) return

    !------------- calc. grid-spacing and check uniformity !

        same = (xpos(npos)& 
             -  xpos(  +1)) / (npos-1)

        uniform = .TRUE.
             
        xtol = same * 1.0d-12

        do  ipos = +1, npos-1
        
            spac = xpos(ipos+1) &
        &        - xpos(ipos+0)
  
            diff = abs(spac - same)
        
            if (diff.gt.xtol) then
            
                uniform = .FALSE.
            
            end if
  
            work% &
        &   cell_spac(ipos) = spac
        
        end do
    
    !------------- reconstruct FHAT over all cells in XPOS !

        if (.not. uniform) then

        call rcon1d(npos,nvar,ndof, &
        &           work%cell_spac, &
        &           fdat,bclo,bcup, &
        &           work%cell_func, &
        &           work,opts)

        else
        
        delx(1) = work%cell_spac(1)
        
        call rcon1d(npos,nvar,ndof, &
        &           delx,    &
        &           fdat,bclo,bcup, &
        &           work%cell_func, &
        &           work,opts)
        
        end if

        mdof = ndof1d(opts%cell_meth)

    !------------- remap FDAT from XPOS to XNEW using FHAT !

        pos0 = +1 ; pos1 = +1

        do  kpos = +1, nnew-1

            pos1 = max(pos1,1)

            do  pos0 = pos1, npos-1
            
    !------ first cell in XPOS overlapping with XNEW(KPOS) !
            
                if (xpos(pos0+1)&
        &      .gt. xnew(kpos+0)) then
                    exit
                end if
            
            end do

            do  pos1 = pos0, npos-1
            
    !------ final cell in XPOS overlapping with XNEW(KPOS) !
    
                if (xpos(pos1+0)&
        &      .ge. xnew(kpos+1)) then
                    exit
                end if
                
            end do
            
            pos1 = pos1-1

    !------------- integrate FHAT across overlapping cells !

            khat = xnew(kpos+1)-xnew(kpos+0)

            do  ivar = +1,nvar
                fnew(  +1,ivar,kpos) = 0.d0
            end do

            do  ipos = pos0 , pos1

                xmid =(xpos(ipos+1)&
            &        + xpos(ipos+0)) * .5d0
                xhat =(xpos(ipos+1)&
            &        - xpos(ipos+0)) * .5d0

    !------------------------------- integration endpoints !
    
                xxlo = max(xpos(ipos+0),  &
            &              xnew(kpos+0))
                xxhi = min(xpos(ipos+1),  &
            &              xnew(kpos+1))

    !------------------------------- local endpoint coords !
    
                sslo = (xxlo - xmid) / xhat
                sshi = (xxhi - xmid) / xhat

    !------------------------------- integral basis vector !
    
                call poly_basis(-1,mdof,sslo,bvlo)
                call poly_basis(-1,mdof,sshi,bvhi)

                do  ivar = +1, nvar

    !--------- integrate FHAT across the overlap XXLO:XXHI !
    
                fnlo = poly_value(-1,mdof, &
            &       work%cell_func(:,ivar,ipos), &
            &           sslo,bvlo)
            
                fnhi = poly_value(-1,mdof, &
            &       work%cell_func(:,ivar,ipos), &
            &           sshi,bvhi)

    !--------- accumulate integral contributions from IPOS !
    
                fnew(+1,ivar,kpos) = &
            &       fnew(+1,ivar,kpos) + &
            &           (fnhi-fnlo) * xhat / khat

                end do

            end do

        end do
        
        return

    end subroutine rmap1d

    end module rmapfn



