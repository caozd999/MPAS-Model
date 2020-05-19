!***********************************************************************
!
!  routine ocn_vel_vmix_tend_implicit_vegetation_drag
!
!> \brief   Computes tendencies for implicit momentum vertical mixing
!>          induced by aquatic vegetation
!> \author  Zhendong Cao
!> \date    May 2020
!> \details
!>  This routine computes the tendencies for implicit vertical mixing for momentum
!>  using vegetation drag.
!>  Except for bottom drag coefficient, routine should be identifcal to
!>  ocn_vel_vmix_tend_implicit above.
!
subroutine ocn_vel_vmix_tend_implicit_vegetation_drag(meshPool, bottomDrag, dt, kineticEnergyCell, & !{{{
                                         vertViscTopOfEdge, layerThickness, &
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
         bottomDrag !< Input: bottomDrag at cell centeres

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

      real(kind=RKIND), pointer :: config_idealized_profile_Cd, &
                                   config_idealized_profile_Cd_marsh, &
                                   config_vegetation_density, &
                                   config_vegetation_diameter, &
                                   config_vegetation_height

      real(kind=RKIND) :: Cd_veg, implicitCd, cff,cff1,cff2, mask_veg

      real (kind=RKIND), dimension(:), allocatable :: R_veg, dab

      err = 0

      if(.not.velVmixOn) return

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

! read in vegetation properties
      call mpas_pool_get_config(ocnConfigs, 'config_idealized_profile_Cd',config_idealized_profile_Cd)
      call mpas_pool_get_config(ocnConfigs, 'config_idealized_profile_Cd_marsh',config_idealized_profile_Cd_marsh)
      call mpas_pool_get_config(ocnConfigs, 'config_vegetation_diameter',config_vegetation_diameter)
      call mpas_pool_get_config(ocnConfigs, 'config_vegetation_height',config_vegetation_height)
      call mpas_pool_get_config(ocnConfigs, 'config_vegetation_density',config_vegetation_density)

      allocate(R_veg(nVertLevels), dab(nVertLevels+1))
      allocate(A(nVertLevels),B(nVertLevels),C(nVertLevels),velTemp(nVertLevels))

      nEdges = nEdgesArray( 1 )

      A(1)=0.0_RKIND

      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
        N = maxLevelEdgeTop(iEdge)
        if (N .gt. 0) then

         ! layerThicknessEdge is computed in compute_solve_diag, and is not available yet,
         ! so recompute layerThicknessEdge here.
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = 1, N
            layerThicknessEdge(k,iEdge) = 0.5_RKIND * (layerThickness(k,cell1) + layerThickness(k,cell2))
         end do

          ! compute the relative height of the vegetation at each vertical layer, R_veg
          dab = 0.0_RKIND
          do k = N,1,-1
            dab(k) = dab(k+1) + layerThicknessEdge(k, iEdge)
            cff = MIN(1.0_RKIND, (dab(k)-config_vegetation_height)/layerThicknessEdge(k,iEdge))
            R_veg(k) = MIN(1-cff, 1.0_RKIND)
          enddo

          ! Define bottom drag coefficient 
          ! cff1 for regular bottom drag; cff2 for vegetation bottom drag
!          cff1 = config_idealized_profile_Cd
          cff1 = gravity * config_idealized_profile_Cd**2.0_RKIND * &
            (0.5_RKIND * (ssh(cell1) + ssh(cell2) + bottomDepth(cell1) + bottomDepth(cell2)))**(-1.0_RKIND/3.0_RKIND)
          cff2 = config_idealized_profile_Cd_marsh * config_vegetation_density * config_vegetation_diameter * &
                         layerThicknessEdge(N,iEdge) * R_veg(N) ! should we assume R_veg(N)=1?

          ! define mask_veg on the cell edge and compute bottom drag coefficient
          ! mask_veg=1 :: the cell on (at least) one side of the edge is vegetated
          ! mask_veg=0 :: cells on both sides of the edge are non-vegetated
          mask_veg = 1.0_RKIND
          if(bottomDrag(cell1).eq.config_idealized_profile_Cd .AND. &
             bottomDrag(cell2).eq.config_idealized_profile_Cd) then
            mask_veg = 0.0_RKIND
            implicitCd = cff1
          elseif (bottomDrag(cell1).eq.config_idealized_profile_Cd_marsh .AND. &
              bottomDrag(cell2).eq.config_idealized_profile_Cd_marsh) then
            implicitCd = cff2
          else
            implicitCd = 0.5_RKIND *(cff1+cff2)
          endif

          ! compute A,B,C
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
         B(1) = 1.0_RKIND - C(1)
         do k = 2, N-1
            B(k) = 1.0_RKIND - A(k) - C(k)
         enddo

         ! Update B and C by vegetation drag
         do k=1, N-1
            Cd_veg = config_idealized_profile_Cd_marsh * config_vegetation_density * &
                        config_vegetation_diameter*R_veg(k)*mask_veg
            B(k) = B(k) - dt*Cd_veg * &
                   sqrt(kineticEnergyCell(k,cell1) + kineticEnergyCell(k,cell2)) &
                   / (layerThicknessEdge(k,iEdge)+layerThicknessEdge(k+1,iEdge))

            C(k) = C(k) + dt*Cd_veg * &
                   sqrt(kineticEnergyCell(k,cell1) + kineticEnergyCell(k,cell2)) &
                   / (layerThicknessEdge(k,iEdge)+layerThicknessEdge(k+1,iEdge))
         enddo

         ! Apply bottom drag boundary condition on the viscous term
         ! second line uses sqrt(2.0*kineticEnergyEdge(k,iEdge))
         ! use implicitCd from spatially variable bottom drag
         B(N) = 1.0_RKIND - A(N) + dt*implicitCd &
              * sqrt(kineticEnergyCell(k,cell1) + kineticEnergyCell(k,cell2)) / layerThicknessEdge(k,iEdge)

         call tridiagonal_solve(A(2:N),B,C(1:N-1),normalVelocity(:,iEdge),velTemp,N)

         normalVelocity(1:N,iEdge) = velTemp(1:N)
         normalVelocity(N+1:nVertLevels,iEdge) = 0.0_RKIND

        end if
      end do
      !$omp end do

      deallocate(A,B,C,velTemp)
      deallocate(dab,R_veg)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vmix_tend_implicit_vegetation_drag!}}}

