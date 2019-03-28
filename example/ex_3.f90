
!   gfortran -cpp -O3 -flto ex_3.f90 -o ex_3
!   ./ex_3

!   "Convergence" testing: remap a smooth (Gaussian) profile 
!   between perturbed grids and compute error in the L2-nrm. 
!   Re-run with N = 50, 100, 200, ... etc (and N-tmp at 80%) 
!   to assess the order of spatial accuracy. With NUll/WENO-
!   LIMIT, all methods should achieve their best-case O(h^p) 
!   scaling.
!

#   include "../src/ppr_1d.f90"

    program ex

        use ppr_1d

        implicit none

        integer, parameter :: npos = 50 ! no. edge (old grid) 
        integer, parameter :: ntmp = 40 ! no. edge (new grid)
        integer, parameter :: nvar = 1  ! no. variables to remap
        integer, parameter :: ndof = 1  ! no. FV DoF per cell
        integer :: ipos

    !------------------------------ position of cell edges !
        real*8  :: xpos(npos), xtmp(ntmp)
        real*8  :: xdel, xmid
        
    !-------------------------------- finite-volume arrays !

    !   Arrays represent a "block" of finite-volume tracers
    !   to remap. The 1st dim. is the no. of DoF per cell,
    !   NDOF=1 is a standard finite-volume scheme where the
    !   data is specified as cell means. NDOF>1 is reserved
    !   for future use with DG-style schemes. NVAR is the
    !   number of tracers to remap. Processing tracers in a
    !   batch is typically more efficient than one-by-one. 
    !   The last dim. is the no. cells (layers) in the grid.

        real*8  :: fdat(ndof,nvar,npos-1)
        real*8  :: ftmp(ndof,nvar,ntmp-1)
        real*8  :: fnew(ndof,nvar,npos-1)

    !------------------------------ method data-structures !
        type(rmap_work) :: work
        type(rmap_opts) :: opts
        type(rcon_ends) :: bc_l(nvar)
        type(rcon_ends) :: bc_r(nvar)

    !------------------------------ define a simple domain !

        xpos(   1) = 0.0d+00
        xpos(npos) = 1.0d+00
        xtmp(   1) = 0.0d+00
        xtmp(ntmp) = 1.0d+00

        xdel = (xpos(npos)-xpos(1))/(npos-1)

        do ipos = +2, npos-1

            xpos(ipos) = (ipos-1) * xdel         

        end do

        xdel = (xtmp(ntmp)-xtmp(1))/(ntmp-1)

        do ipos = +2, ntmp-1

            xtmp(ipos) = (ipos-1) * xdel         

        end do

    !------------------------------ setup some simple data !

        do ipos = +1, npos-1

            xmid = xpos(ipos+0) * 0.5d+00 &
    &            + xpos(ipos+1) * 0.5d+00

            fdat(1,1,ipos) = &
    &   .8d+0 * exp( -75.0d+0 * (xmid - 0.275d+0) ** 2 ) &
    & + .9d+0 * exp(-100.0d+0 * (xmid - 0.500d+0) ** 2 ) &
    & + 1.d+0 * exp(-125.0d+0 * (xmid - 0.725d+0) ** 2 )

        end do

    !------------------------------ specify method options !

        opts%edge_meth = p3e_method     ! 3rd-order edge interp.
        opts%cell_meth = ppm_method     ! PPM method in cells
        opts%cell_lims = null_limit     ! monotone limiter
        
    !------------------------------ set BC.'s at endpoints !

        bc_l%bcopt = bcon_loose         ! "loose" = extrapolate
        bc_r%bcopt = bcon_loose

    !------------------------------ init. method workspace !

        call work%init(npos,nvar,opts)

    !------------------------------ re-map from new-to-tmp !

        fnew = fdat

        do ipos = +1, +1000

        call rmap1d(npos,ntmp,nvar,ndof, &
    &               xpos,xtmp,fnew,ftmp, &
    &               bc_l,bc_r,work,opts)

    !------------------------------ re-map from tmp-to-new !

        call rmap1d(ntmp,npos,nvar,ndof, &
    &               xtmp,xpos,ftmp,fnew, &
    &               bc_l,bc_r,work,opts)

        end do

    !------------------------------ clear method workspace !

        call work%free()

    !------------------------------ dump results to stdout !

        print*,"Cell data: [INIT.] [FINAL]"

        do ipos = +1, npos-1

            print *, fdat(1,1,ipos) &
    &              , fnew(1,1,ipos)

        end do

    !------------------------------ calc. errors in L2 nrm !

        xdel = (xpos(npos)-xpos(1))/(npos-1)

        print*, "Error (L2): " , &
    &   sqrt(sum(xdel * (fdat(1,1,:) - fnew(1,1,:)) ** 2))

    end program



