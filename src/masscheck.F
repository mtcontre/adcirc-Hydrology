!***********************************************************************
!                 M O D U L E   M A S S   C H E C K
!***********************************************************************
!     Computes and stores element-wise mass conservation errors
!     sb 06102021
!***********************************************************************
      module masscheck
      use sizes, only : sz
      use global, only : setMessageSource, unsetMessageSource,
     &  DEBUG, ERROR, allMessage, screenMessage

!      real(sz), allocatable, target :: errElem(:)   ! Errors for each element to be used to compute nodal attributions

      contains

!!----------------------------------------------------------------------
!!                  S U B R O U T I N E
!!             I N I T   M A S S   C H E C K
!!----------------------------------------------------------------------
!!     Initialize variables in the module
!!----------------------------------------------------------------------
!      subroutine initMassCheck()
!      use sizes, only : mne
!      use global, only : nodecode, noutge
!      implicit none
!   
!      call setMessageSource("initMassCheck")
!#if defined(MASSCHECK_TRACE) || defined(ALL_TRACE)
!      call allMessage(DEBUG,"Enter.")
!#endif
!
!      allocate(errElem(mne))
!
!#if defined(MASSCHECK_TRACE) || defined(ALL_TRACE)
!      call allMessage(DEBUG,"Return.")
!#endif
!      call unsetMessageSource()
!!----------------------------------------------------------------------
!!      end subroutine
!!----------------------------------------------------------------------

!----------------------------------------------------------------------
!                  S U B R O U T I N E
!     C O M P U T E   M A S S   C O N S E R V   E R R O R
!----------------------------------------------------------------------
!     Computes element-wise nodal mass conservation errors.
!     + errMassSum: Accumulated errors observed for each for.63 output period
!     + errMassMax: Signed maximum errors observed for each for.63 output period
!     + errMassSumFull: Accumulated errors observed for an entire simulation lifetime
!     + errMassMaxFull: Signed maximum errors observed for an entire simulation lifetime
!     Note: Errors are distribuited to nodes once computed element-wise.
!     sb 06102021
!----------------------------------------------------------------------
      subroutine computeMassConservError(timeloc)
      use global, only : 
     &  dtdp,eta1,eta2,uu1,uu2,vv1,vv2,qx1,qx2,qy1,qy2,
     &  ifnlfa,nodecode,noutge
      use mesh, only : np,x,y,dp,ne,nm,areas ! areas = 2 x Area
!      use gwce, only : a00,b00,c00
      use write_output, only :
     &  errMassSum,errMassMax,errMassMax_time,
     &  errMassSumFull,errMassMaxFull,errMassMaxFull_time,
     &  errMassSumDescript,errMassMaxDescript

      implicit none

      real(8), intent(in) :: timeloc ! time in s since cold start
      integer :: i ! node index
      integer :: j ! element index
      integer :: e ! edge index of element
      integer :: nm1,nm2,nm3 ! element node global index
      integer :: nmi         ! another element node global index
      integer :: l1,l2 ! edge node local index
      integer :: i1,i2 ! edge node global index
      real(sz) :: x1,y1,x2,y2 ! coordinates of edge nodes
      real(sz) :: H1,H2 ! total water depth at edge nodes
      real(sz) :: dx,dy ! edge tangential vector
      real(sz) :: L ! edge length of element
      real(sz) :: nx,ny ! edge normal vector
      real(sz) :: V1,V2 ! velocities normal inward to edge
      real(sz) :: P ! mass change computed from elevation change
      real(sz) :: Q ! mass flow computed from flux through edges
      real(sz) :: QK0,QK1,QK2 ! intermidiate variables for K = -2, -1 and 0
      real(sz) :: Err ! mass conservation error
      real(sz) :: Ptot ! P summed over all elements
      real(sz) :: Qtot ! Q summed over all elements
      real(sz) :: Etot ! Err summed over all elements

      call setMessageSource("computeMassConservError")
#if defined(MASSCHECK_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      !
      ! Leave subroutine if the latter proccess is not needed.
      if (noutge.eq.0) then
#if defined(MASSCHECK_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         return
      endif
      !
      ! Reset arrays for each fort.63 output period
      if (errMassSumDescript % spoolCounter.eq.0) then
         errMassSum(:) = 0.d0
      endif
      if (errMassMaxDescript % spoolCounter.eq.0) then
         errMassMax(:) = 0.d0
         errMassMax_time(:) = 0.d0
      endif
      !
      ! Initialize
      Ptot = 0.d0
      Qtot = 0.d0
      Etot = 0.d0
      !
      ! Compute element-wise error
      do j=1,ne
         nm1 = nm(j,1)
         nm2 = nm(j,2)
         nm3 = nm(j,3)

         ! skip if completely dry
         if( (nodecode(nm1) + nodecode(nm2) + 
     &        nodecode(nm3)) .eq. 0) then
            cycle
         endif

         P = 0.d0
         Q = 0.d0
         do e=1,3
            ! prep node/edge properties
            l1 = mod(e  ,3) + 1
            l2 = mod(e+1,3) + 1
            i1 = nm(j,l1)
            i2 = nm(j,l2)
            x1 = x(i1)
            x2 = x(i2)
            y1 = y(i1)
            y2 = y(i2)
            dx = x2 - x1
            dy = y2 - y1
            L = sqrt(dx*dx + dy*dy)
            nx = -dy / L
            ny = dx / L
            ! partial calculation of mass change
            P = P + (eta2(i1) - eta1(i1))
            ! partial calculation of mass flow
            V1 = UU1(i1)*nx + VV1(i1)*ny
            V2 = UU1(i2)*nx + VV1(i2)*ny
            H1 = dp(i1) + ifnlfa*eta1(i1)
            H2 = dp(i2) + ifnlfa*eta1(i2)
            QK1 = 2.0*V1*H1 + V1*H2 + V2*H1 + 2.0*V2*H2

            V1 = UU2(i1)*nx + VV2(i1)*ny
            V2 = UU2(i2)*nx + VV2(i2)*ny
            H1 = dp(i1) + ifnlfa*eta2(i1)
            H2 = dp(i2) + ifnlfa*eta2(i2)
            QK2 = 2.0*V1*H1 + V1*H2 + V2*H1 + 2.0*V2*H2

            Q = Q + L*(QK1 + QK2)

         enddo

         P = P*areas(j)*0.5d0/3.d0 ! complete mass change calculation
         Q = Q*dtdp/12.d0          ! complete mass flow calculation
         Err = P - Q               ! compute mass conservation error

         Ptot = Ptot + P
         Qtot = Qtot + Q
         Etot = Etot + Err

         do i=1,3
            nmi = nm(j,i)

            ! errMassSum
            errMassSum(nmi) = errMassSum(nmi) + Err/3.d0

            ! errMassMax
            if(abs(errMassMax(nmi)).lt.abs(Err)) then
               errMassMax(nmi) = Err
               errMassMax_time(nmi) = timeLoc
            endif

            ! errMassSumFull
            errMassSumFull(nmi) = errMassSumFull(nmi) + Err/3.d0

            ! errMassMaxFull
            if(abs(errMassMaxFull(nmi)).lt.abs(Err)) then
               errMassMaxFull(nmi) = Err
               errMassMaxFull_time(nmi) = timeLoc
            endif
         enddo
      enddo

!      write(*,101), Ptot, Qtot, Etot, Etot/Ptot
!101   FORMAT(1X,'tot:',E11.3,E11.3,E11.3,F10.2)
#if defined(MASSCHECK_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
!----------------------------------------------------------------------
      end subroutine
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      end module masscheck
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
          
