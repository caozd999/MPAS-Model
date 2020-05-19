! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_forcing_explicit_vegetation_drag
!
!> \brief MPAS ocean explicit vegetation drag
!> \author Zhendong Cao
!> \date   April 2020
!> \details
!>  This module contains the routine for computing
!>  tendencies from explicit vegetation drag.
!
!-----------------------------------------------------------------------

module ocn_vel_forcing_explicit_vegetation_drag

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_constants, only : gravity
   use ocn_constants
   use ocn_forcing

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_vel_forcing_explicit_vegetation_drag_tend, &
             ocn_vel_forcing_explicit_vegetation_drag_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------
   logical :: explicitVegetationDragOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_forcing_explicit_vegetation_drag_tend
!
!> \brief   Computes tendency term from explicit vegetation drag
!> \author  Zhendong Cao
!> \date    23 April 2020
!> \details
!>  This routine computes the explicit vegetation drag tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_explicit_vegetation_drag_tend(meshPool, ssh, bottomDrag, & !{{{
                            normalVelocity, layerThickness, kineticEnergyCell, layerThicknessEdge, tend, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity,       &!< Input: velocity
         kineticEnergyCell,    &!< Input: kinetic energy at cell
         layerThickness         !< Input: layer thickness
      real (kind=RKIND), dimension(:), intent(in) ::  &
         ssh,                  &!< Input: sea surface height
         bottomDrag             !< Input: bottom drag coefficient
      type (mpas_pool_type), intent(in) :: &
         meshPool               !< Input: mesh information
      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         layerThicknessEdge,   &!< Input/Output: thickness at edge
         tend                   !< Input/Output: velocity tendency
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

      integer :: N, iEdge, k, cell1, cell2, nEdges
      integer, dimension(:), pointer :: nEdgesArray
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, pointer :: nVertLevels
      real(kind=RKIND), dimension (:), pointer :: bottomDepth

      real(kind=RKIND), pointer :: config_idealized_profile_Cd, &
                                   config_idealized_profile_Cd_marsh, &
                                   config_vegetation_density, &
                                   config_vegetation_diameter, &
                                   config_vegetation_height
      real(kind=RKIND) :: Cd_veg, explicitCd, cff,cff1,cff2, mask_veg
      real (kind=RKIND), dimension(:), allocatable :: R_veg, dab


      err = 0
      if ( .not. explicitVegetationDragOn ) return

      call mpas_timer_start('vel explicit vegetation drag')

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)

! read in vegetation properties
      call mpas_pool_get_config(ocnConfigs, 'config_idealized_profile_Cd',config_idealized_profile_Cd)
      call mpas_pool_get_config(ocnConfigs, 'config_idealized_profile_Cd_marsh',config_idealized_profile_Cd_marsh)
      call mpas_pool_get_config(ocnConfigs, 'config_vegetation_diameter',config_vegetation_diameter)
      call mpas_pool_get_config(ocnConfigs, 'config_vegetation_height',config_vegetation_height)
      call mpas_pool_get_config(ocnConfigs, 'config_vegetation_density',config_vegetation_density)

      allocate(R_veg(nVertLevels), dab(nVertLevels+1))

      nEdges = nEdgesArray( 1 )

      !$omp do schedule(runtime) private(k, cell1, cell2)

      do iEdge = 1, nEdges
        N =  maxLevelEdgeTop(iEdge)
        if (N .gt. 0) then
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)

          do k = 1, N
            layerThicknessEdge(k, iEdge) = 0.5_RKIND * (layerThickness(k,cell1) + layerThickness(k,cell2))
          enddo

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
            explicitCd = cff1
          elseif (bottomDrag(cell1).eq.config_idealized_profile_Cd_marsh .AND. &
              bottomDrag(cell2).eq.config_idealized_profile_Cd_marsh) then
            explicitCd = cff2
          else
            explicitCd = 0.5_RKIND *(cff1+cff2)
          endif

          ! compute the bottom drag
          ! Explicit bottom drag term:
          ! du/dt = ... - c |u| u / h
          ! appied to bottom layer only.
          tend(N,iEdge) = tend(N,iEdge) - 0.5_RKIND*explicitCd* &
                          sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) * normalVelocity(N,iEdge)  &
                          / layerThicknessEdge(N,iEdge)

          ! compute the vegetation drag effects on tendency at different vertical layer
          ! Cd_veg = Cd_marsh * density * diameter * layerThickness * R_veg
          ! where Cd_marsh = config_idealized_profile_Cd_marsh
          if (N .ge. 2) then
             do k = 1, N-1
               !Cd_veg = config_idealized_profile_Cd_marsh*config_vegetation_density * &
               !         config_vegetation_diameter*layerThicknessEdge(k,iEdge)*R_veg(k)*mask_veg
               !tend(k,iEdge) = tend(k,iEdge) - Cd_veg* &
               !               sqrt(kineticEnergyCell(k,cell1) + kineticEnergyCell(k,cell2)) * normalVelocity(k,iEdge)  &
               !                / layerThicknessEdge(k,iEdge)
               ! simplify the above commented equations to the bottom ones to get rid of layerThicknessEdge
               config_idealized_profile_Cd_marsh = 1.0_RKIND
               Cd_veg = config_idealized_profile_Cd_marsh * config_vegetation_density * &
                        config_vegetation_diameter*R_veg(k)*mask_veg
               tend(k,iEdge) = tend(k,iEdge) -0.5_RKIND* Cd_veg* &
                           sqrt(kineticEnergyCell(k,cell1) + kineticEnergyCell(k,cell2)) * normalVelocity(k,iEdge)
             enddo
           endif

        endif
      enddo
      !$omp end do

      deallocate(R_veg, dab)

      call mpas_timer_stop('vel explicit vegetation drag')

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_explicit_vegetation_drag_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_forcing_explicit_vegetation_drag_init
!
!> \brief   Initializes ocean explicit vegetation drag forcing
!> \author  Zhendong Cao
!> \date    April 2020
!> \details
!>  This routine initializes quantities related to explicit vegetation drag
!>  in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_explicit_vegetation_drag_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      logical, pointer :: config_use_explicit_vegetation_drag
      logical, pointer :: config_idealized_profile

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_idealized_profile', config_idealized_profile)
      call mpas_pool_get_config(ocnConfigs, 'config_use_explicit_vegetation_drag', config_use_explicit_vegetation_drag)

      if (config_idealized_profile .and. config_use_explicit_vegetation_drag ) then
        explicitVegetationDragOn = .TRUE.
      endif
   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_explicit_vegetation_drag_init!}}}

!***********************************************************************

end module ocn_vel_forcing_explicit_vegetation_drag

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
