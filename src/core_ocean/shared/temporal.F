!***********************************************************************
!
!  routine ocn_vel_vmix_tend_implicit_vegetation_drag
!
!> \brief   Computes tendencies for implicit momentum vertical mixing
!>          with spatially-variable bottom drag and vegetation drag
!> \template  Mark Petersen, Phillip J. Wolfram, Zhendong Cao
!> \date    September 2011, September 2019, April 2020
!> \author Zhendong Cao
!> \date April 2020
!> \details
!>  This routine computes the tendencies for implicit vertical mixing for momentum
!>  using computed coefficients from spatially-variable bottom drag and vegetation drag.
!>  Except for bottom drag coefficient, routine should be identifcal to
!>  ocn_vel_vmix_tend_implicit above.  Cd uses Mannings' n values for the
!>  Cd=g*n^2*h^{-1/3}. Vegetation drag ....
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vmix_tend_implicit_vegetation_drag(meshPool, & !{{{
                                         manningN_veg, diam_veg, dens_veg, hght_veg, &
                                         dt, kineticEnergyCell, vertViscTopOfEdge, layerThickness, &
                                         layerThicknessEdge, normalVelocity, ssh, bottomDepth, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         kineticEnergyCell        !< Input: kinetic energy at cell

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertViscTopOfEdge !< Input: vertical mixing coefficients

      real (kind=RKIND), intent(in) :: &
         dt            !< Input: time step

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness !< Input: thickness at cell center

       real (kind=RKIND), dimension(:), intent(in) :: &
         manningN_veg !< Input: Vegetation Manning coefficient at cell centeres

       real (kind=RKIND), dimension (:), intent(in) :: &
         diam_veg !< Input: Vegetation stem diameter (m)

       real (kind=RKIND), dimension (:), intent(in) :: &
         dens_veg !< Input: Vegetation stem density (m^-2)

       real (kind=RKIND), dimension (:), intent(in) :: &
         hght_veg !< Input: Vegetation stem height (m)

       real (kind=RKIND), dimension(:), pointer :: ssh

       real (kind=RKIND), dimension(:), pointer :: bottomDepth

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         normalVelocity             !< Input: velocity

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         layerThicknessEdge        !< Input: thickness at edge

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, k, cell1, cell2, N, nEdges
      integer, pointer :: nVertLevels
      real (kind=RKIND) :: implicitCd
      integer, dimension(:), pointer :: nEdgesArray

      integer, dimension(:), pointer :: maxLevelEdgeTop

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(:), allocatable :: A, B, C, velTemp

      real (kind=RKIND), dimension(:), allocatable :: Cd_veg, R_veg, dab

      real (kind=RKIND) :: cff

      err = 0

      if(.not.velVmixOn) return

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      nEdges = nEdgesArray( 1 )

      allocate(A(nVertLevels),B(nVertLevels),C(nVertLevels),velTemp(nVertLevels))

      allocate(Cd_veg(nVertLevels), R_veg(nVertLevels),dab(nVertLevels+1))

      A(1)=0.0_RKIND

      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
        N = maxLevelEdgeTop(iEdge)
        if (N .gt. 0) then

         ! Compute A(k), B(k), C(k)
         ! layerThicknessEdge is computed in compute_solve_diag, and is not available yet,
         ! so recompute layerThicknessEdge here.
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = 1, N
            layerThicknessEdge(k,iEdge) = 0.5_RKIND * (layerThickness(k,cell1) + layerThickness(k,cell2))
         end do

         ! average cell-based implicit bottom drag to edges and convert Mannings n to Cd
         implicitCd = gravity*(0.5_RKIND*(manningN_veg(cell1) + manningN_veg(cell2)))**2.0 * &
           0.5_RKIND * (ssh(cell1) + ssh(cell2) + bottomDepth(cell1) + bottomDepth(cell2))**(-1.0_RKIND/3.0_RKIND)
         ! compute the relative height of the vegetation at vertical layer, R_veg
         ! code adopted from COAWST
         dab = 0.0_RKIND
         do k = N,1,-1
           dab(k) = dab(k+1) + layerThicknessEdge(k,iEdge)
           cff = MIN(1,(dab(k)-(hght_veg(cell1)+hght_veg(cell2)/2.0)/layerThickness(k,iEdge))
           R_veg(k) = MIN(1-cff, 1)
         enddo

         ! compute the vegetation drag
         do k = 1, N
           Cd_veg(k) = implicitCd*(diam_veg(cell1)+diam_veg(cell2))/2.0_RKIND * &
                                  (dens_veg(cell1)+dens_veg(cell2))/2.0_RKIND * &
                                  layerThicknessEdge(k,iEdge)*R_veg(k)
         enddo


         ! A is lower diagonal term
         do k = 2, N
            A(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k,iEdge) &
               / (layerThicknessEdge(k-1,iEdge) + layerThicknessEdge(k,iEdge)) &
               / layerThicknessEdge(k,iEdge)
         enddo

         ! C is upper diagonal term
         do k = 1, N-1
            C(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k+1,iEdge) &
               / (layerThicknessEdge(k,iEdge) + layerThicknessEdge(k+1,iEdge)) &
               / layerThicknessEdge(k,iEdge)
         enddo

         ! B is diagonal term
!         B(1) = 1.0_RKIND - C(1)
!         do k = 2, N-1
!            B(k) = 1.0_RKIND - A(k) - C(k)
!         enddo

         ! Apply bottom drag boundary condition on the viscous term
         ! second line uses sqrt(2.0*kineticEnergyEdge(k,iEdge))
         ! use implicitCd from spatially variable bottom drag
!         B(N) = 1.0_RKIND - A(N) + dt*implicitCd &
!              * sqrt(kineticEnergyCell(k,cell1) + kineticEnergyCell(k,cell2)) / layerThicknessEdge(k,iEdge)
         ! Considering that A(1)=0 and C(k)=0, we define B in a more uniform way
         do k = 1, N
           B(k) = 1.0_RKIND - A(k) - C(k) + &
             dt*Cd_veg(k) * sqrt(kineticEnergyCell(k,cell1) + kineticEnergyCell(k,cell2)) / layerThicknessEdge(k,iEdge)
         enddo

         call tridiagonal_solve(A(2:N),B,C(1:N-1),normalVelocity(:,iEdge),velTemp,N)

         normalVelocity(1:N,iEdge) = velTemp(1:N)
         normalVelocity(N+1:nVertLevels,iEdge) = 0.0_RKIND

        end if
      end do
      !$omp end do

      deallocate(A,B,C,velTemp)
      deallocate(Cd_veg,R_veg, dab)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vmix_tend_implicit_vegetation_drag!}}}


