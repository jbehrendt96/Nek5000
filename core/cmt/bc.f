C> @file bc.f Boundary condition routines

C> \ingroup bcond
C> @{
C> Determining rind state for Dirichlet boundary conditions
      subroutine InviscidBC(flux)
!-------------------------------------------------------------------------------
! JH091514 A fading copy of RFLU_ModAUSM.F90 from RocFlu
!-------------------------------------------------------------------------------

!#ifdef SPEC
!      USE ModSpecies, ONLY: t_spec_type
!#endif
      include 'SIZE'
      include 'INPUT' ! if3d
      include 'CMTDATA' ! do we need this without outflsub?
!     include 'TSTEP' ! for ifield?
      include 'GEOM'

! ==============================================================================
! Arguments
! ==============================================================================
      real flux(lx1*lz1,2*ldim,nelt,toteq)

! ==============================================================================
! Locals
! ==============================================================================

      integer e,f,fdim,i,k,nxz,nface,ifield,eq

      common /nekcb/ cb
      character*3 cb
      parameter(lstate=5) ! KG is expecting 4.
                                      ! hopefully, this won't cause problems
      COMMON /SCRNS/ wminus(lstate,lxz),wplus(lstate,lxz),
     >               jaminus(3,lxz),japlus(3,lxz),
     >               uminus(toteq,lxz),uplus(toteq,lxz),
     >               flx(toteq,lxz)
      real wminus,wplus,jaminus,japlus,uminus,uplus,flx
      external kennedygruber,llf_euler

      fdim=ldim-1
      nface = 2*ldim
      nxz   = lx1*lz1
      ifield= 1 ! You need to figure out the best way of dealing with
                ! this variable
!     if (outflsub)then
!        call maxMachnumber
!     endif

      do e=1,nelt
      do f=1,nface

         cb=cbc(f,e,ifield)
         if (cb.ne.'E  '.and.cb.ne.'P  ') then ! cbc bndy
! BC routines to fill face points
! facind + userbc for wminus, uminus
! dirichlet routines for wplus, uplus
            if (cb.eq.'v  ' .or. cb .eq. 'V  ') then
              call inflow(f,e,wminus,wplus,uminus,uplus,lstate)
            elseif (cb.eq.'O  ') then
              call outflow(f,e,wminus,wplus,uminus,uplus,lstate)
            elseif (cb .eq. 'W  ' .or. cb .eq.'I  '.or.cb .eq.'SYM')then
              call wallbc_inviscid(f,e,wminus,wplus,uminus,uplus,lstate)
            endif 

! convert surface normals into metric terms for two-point fluxes (just
! face Jacobian times normal vector)

            do i=1,nxz
               jaminus(1,i)=unx(i,1,f,e)
               jaminus(2,i)=uny(i,1,f,e)
            enddo
            if(if3d) then
               do i=1,nxz
                  jaminus(3,i)=unz(i,1,f,e)
               enddo
            else
               do i=1,nxz
                  jaminus(3,i)=0.0
               enddo
            endif
! boundaries are watertight lol
            call copy(japlus,jaminus,3*nxz)

! two-point flux
            call sequential_flux(flx,wminus,wplus,uminus,uplus,jaminus,
     >                           japlus,kennedygruber,lstate,nxz)
            do eq=1,toteq
            do i=1,nxz
            flux(i,f,e,eq)=flux(i,f,e,eq)+flx(eq,i)*jface(i,1,f,e)
            enddo
            enddo
! stabilization flux
            call sequential_flux(flx,wminus,wplus,uminus,uplus,jaminus,
     >                           japlus,llf_euler,lstate,nxz)
            do eq=1,toteq
            do i=1,nxz
            flux(i,f,e,eq)=flux(i,f,e,eq)+flx(eq,i)*jface(i,1,f,e)
            enddo
            enddo

         endif
      enddo
      enddo

C> @}
      return
      end

C> \ingroup bcond
C> @{
C> Mask to make sure Fsharp doesn't clobber boundary faces, where gs_op is null
C> This routine intents to take a real array for all face points, bmask, and
C> only zero out faces on boundaries. It is thus not limited to an array only
C> of indicators.
      subroutine bcmask_cmt(bmsk)
      include 'SIZE'
      include 'INPUT' ! do we need this?
      include 'CMTDATA' ! do we need this without outflsub?
!     include 'TSTEP' ! for ifield?
      include 'DG'

! ==============================================================================
! Arguments
! ==============================================================================
      integer nstate,nflux
      real bmsk(lx1*lz1,2*ldim,nelt)

! ==============================================================================
! Locals
! ==============================================================================

      integer e,f,fdim,i,k,nxz,nface,ifield

      common /nekcb/ cb
      character*3 cb

      fdim=ldim-1
      nface = 2*ldim
      nxz   = lx1*lz1
      ifield= 1 ! You need to figure out the best way of dealing with
                ! this variable

      do e=1,nelt
      do f=1,nface

         cb=cbc(f,e,ifield)
         if (cb.ne.'E  '.and.cb.ne.'P  ') then ! cbc bndy
            call rzero(bmsk(1,f,e),nxz)
         endif
      enddo
      enddo

C> @}
      return
      end

!-----------------------------------------------------------------------

C> \ingroup bcond
C> @{
C> Determining IGU contribution to boundary flux. 0 for artificial
C> viscosity, and strictly interior for physical viscosity.
      subroutine bcflux(flux,agradu,qminus)
! Need proper indexing and nekasgn & cmtasgn calls
      include 'SIZE'
      include 'INPUT'
      include 'DG'
!     include 'NEKUSE'
      include 'TSTEP' ! wait how do we know what ifield is?
      integer e,eq,f
      real flux  (lx1*lz1,2*ldim,nelt,toteq)
      real agradu(lx1*lz1,2*ldim,nelt,toteq)
!     real qminus(lx1*lz1,2*ldim,nelt,nqq) ! include CMTDATA?
      real qminus(*) ! 'scuse me. comin' through
      common /nekcb/ cb
      character*3 cb

      nfaces=2*ldim
      nxz=lx1*lz1
      ifield=1

      do e=1,nelt
         do f=1,nfaces
            if (cbc(f, e, ifield).ne.'E  '.and.
     >          cbc(f, e, ifield).ne.'P  ') then ! cbc bndy
               cb=cbc(f,e,ifield)
               if (cb .eq. 'I  ') then ! NEUMANN CONDITIONS GO HERE
!-------------------------------------------------------------
! JH112216 HARDCODING ADIABATIC WALL. DO SMARTER SOON
                  call rzero(flux(1,f,e,1),nxz)
                  do eq=2,ldim+1
                     call copy(flux(1,f,e,eq),agradu(1,f,e,eq),nxz)
                  enddo
! METHOD "B", ADIABATIC NO-SLIP
                  call rzero(flux(1,f,e,toteq),nxz)
! METHOD "A", ADIABATIC NO-SLIP augments with viscous work. triage below
! because, predictably, NOW I need to computate AgradU on surfaces and I don't
! have general code for that.
                  call a5adiabatic_wall(flux(1,1,1,toteq),f,e,agradu,
     >                                  qminus)
! JH112216 HARDCODING ADIABATIC WALL. DO SMARTER SOON
!-------------------------------------------------------------
!                 cbu=cb
!                 do eq=1,toteq
!                    call userflux(flux(1,f,e,eq)) ! replace this with userbc
!                 enddo
               else  ! if (cb .eq. 'SYM') then ! NEED PHYSICAL VISC TEST
! JH031617 But this code block basically guarantees that artificial viscosity
!          does not contribute to viscous fluxes at boundaries.
                  do eq=1,toteq
                     call rzero(flux(1,f,e,eq),nxz)
                  enddo
               endif
            endif
         enddo
      enddo

C> @}
      return
      end

      subroutine a5adiabatic_wall(eflx,f,e,dU,wstate)
      include 'SIZE'
      include 'INPUT'
      include 'GEOM' ! for UNX under ADIABATIC WALL METHOD "A"
      include 'CMTDATA'
      real eflx  (lx1*lz1,2*ldim,nelt) ! better be zero on entry
      real dU    (lx1*lz1,2*ldim,nelt,toteq)
      real wstate(lx1*lz1,2*ldim,nelt,nqq)
      common /scrns/ flxscr(lx1*lz1)
      real flxscr
      integer e,f

      nxz=lx1*lz1

      call rzero(eflx(1,f,e),nxz)
      call rzero(hface,nxz)

      call a51dUadia(flxscr,f,e,dU,wstate)
      call add2col2(eflx(1,f,e),flxsxcr,unx(1,1,f,e),nxz)
      call a52dUadia(flxscr,f,e,dU,wstate)
      call add2col2(eflx(1,f,e),flxsxcr,uny(1,1,f,e),nxz)
      if (if3d) then
         call a53dUadia(flxscr,f,e,dU,wstate)
         call add2col2(eflx(1,f,e),flxsxcr,unz(1,1,f,e),nxz)
      endif
      return
      end

      subroutine a51dUadia(flux,f,ie,dU,wstate)
! same as A51 for volume flux, but
! 1. uses surface storage of quantities in wstate <-qminus (intent(in))
! 2. SETS K=0. ADIABATIC WALLS HAVE VISCOUS HEATING, BUT DON'T CONDUCT
      include 'SIZE'
      include 'CMTDATA'
      real wstate(lx1*lz1,2*ldim,nelt,nqq)
      real dU    (lx1*lz1,2*ldim,nelt,toteq,3)
      real flux  (lx1*ly1*lz1)
      real K,E,kmcvmu,lambdamu
      integer f
      npt=lx1*lz1

      do i=1,npt
         dU1x=dU(i,f,ie,1,1)
         dU2x=dU(i,f,ie,2,1)
         dU3x=dU(i,f,ie,3,1)
         dU4x=dU(i,f,ie,4,1)
         dU5x=dU(i,f,ie,5,1)
         dU1y=dU(i,f,ie,1,2)
         dU2y=dU(i,f,ie,2,2)
         dU3y=dU(i,f,ie,3,2)
         dU4y=dU(i,f,ie,4,2)
         dU5y=dU(i,f,ie,5,2)
         dU1z=dU(i,f,ie,1,3)
         dU2z=dU(i,f,ie,2,3)
         dU3z=dU(i,f,ie,3,3)
         dU4z=dU(i,f,ie,4,3)
         dU5z=dU(i,f,ie,5,3)
         rho   =wstate(i,f,ie,jrho)
         cv    =wstate(i,f,ie,jcvf)/rho
         lambda=wstate(i,f,ie,jlamf)
         mu    =wstate(i,f,ie,jmuf)
         K     =0.0 ! ADIABATIC HARDCODING
         u1    =wstate(i,f,ie,jux)
         u2    =wstate(i,f,ie,juy)
         u3    =wstate(i,f,ie,juz)
         E     =wstate(i,f,ie,ju5)/rho
         lambdamu=lambda+mu
         kmcvmu=K-cv*mu
         flux(i)=
     >(K*dU5x+cv*lambda*u1*dU4z-kmcvmu*u3*dU4x+cv*lambda*u1*dU3y
     1   -kmcvmu*u2*dU3x+cv*mu*u3*dU2z+cv*mu*u2*dU2y+(cv*lambda-
     2   K+2*cv*mu)*u1*dU2x-cv*lambdamu*u1*u3*dU1z-cv*lambdamu
     3   *u1*u2*dU1y+(K*u3**2-cv*mu*u3**2+K*u2**2-cv*mu*u2**2-cv*la
     4   mbda*u1**2+K*u1**2-2*cv*mu*u1**2-E*K)*dU1x)/(cv*rho)
      enddo
      return
      end

      subroutine a52dUadia(flux,f,ie,dU,wstate)
! same as A52 for volume flux, but
! 1. uses surface storage of quantities in wstate <-qminus (intent(in))
! 2. SETS K=0. ADIABATIC WALLS HAVE VISCOUS HEATING, BUT DON'T CONDUCT
      include 'SIZE'
      include 'CMTDATA'
      real wstate(lx1*lz1,2*ldim,nelt,nqq)
      real dU    (lx1*lz1,2*ldim,nelt,toteq,3)
      real flux  (lx1*ly1*lz1)
      real K,E,kmcvmu,lambdamu
      integer f
      npt=lx1*lz1

      do i=1,npt
         dU1x=dU(i,f,ie,1,1)
         dU2x=dU(i,f,ie,2,1)
         dU3x=dU(i,f,ie,3,1)
         dU4x=dU(i,f,ie,4,1)
         dU5x=dU(i,f,ie,5,1)
         dU1y=dU(i,f,ie,1,2)
         dU2y=dU(i,f,ie,2,2)
         dU3y=dU(i,f,ie,3,2)
         dU4y=dU(i,f,ie,4,2)
         dU5y=dU(i,f,ie,5,2)
         dU1z=dU(i,f,ie,1,3)
         dU2z=dU(i,f,ie,2,3)
         dU3z=dU(i,f,ie,3,3)
         dU4z=dU(i,f,ie,4,3)
         dU5z=dU(i,f,ie,5,3)
         rho   =wstate(i,f,ie,jrho)
         cv    =wstate(i,f,ie,jcvf)/rho
         lambda=wstate(i,f,ie,jlamf)
         mu    =wstate(i,f,ie,jmuf)
         K     =0.0 ! ADIABATIC HARDCODING
         u1    =wstate(i,f,ie,jux)
         u2    =wstate(i,f,ie,juy)
         u3    =wstate(i,f,ie,juz)
         E     =wstate(i,f,ie,ju5)/rho
         lambdamu=lambda+mu
         kmcvmu=K-cv*mu
         flux(i)=
     >(K*dU5y+cv*lambda*u2*dU4z-kmcvmu*u3*dU4y+cv*mu*u3*dU3z+(cv
     1   *lambda-K+2*cv*mu)*u2*dU3y+cv*mu*u1*dU3x-kmcvmu*u1*dU2y+
     2   cv*lambda*u2*dU2x-cv*lambdamu*u2*u3*dU1z+(K*u3**2-cv*mu
     3   *u3**2-cv*lambda*u2**2+K*u2**2-2*cv*mu*u2**2+K*u1**2-cv*mu*
     4   u1**2-E*K)*dU1y-cv*lambdamu*u1*u2*dU1x)/(cv*rho)
      enddo
      return
      end

      subroutine a53dUadia(flux,f,ie,dU,wstate)
! same as A53 for volume flux, but
! 1. uses surface storage of quantities in wstate <-qminus (intent(in))
! 2. SETS K=0. ADIABATIC WALLS HAVE VISCOUS HEATING, BUT DON'T CONDUCT
      include 'SIZE'
      include 'CMTDATA'
      real wstate(lx1*lz1,2*ldim,nelt,nqq)
      real dU    (lx1*lz1,2*ldim,nelt,toteq,3)
      real flux  (lx1*ly1*lz1)
      real K,E,kmcvmu,lambdamu
      integer f
      npt=lx1*lz1

      do i=1,npt
         dU1x=dU(i,f,ie,1,1)
         dU2x=dU(i,f,ie,2,1)
         dU3x=dU(i,f,ie,3,1)
         dU4x=dU(i,f,ie,4,1)
         dU5x=dU(i,f,ie,5,1)
         dU1y=dU(i,f,ie,1,2)
         dU2y=dU(i,f,ie,2,2)
         dU3y=dU(i,f,ie,3,2)
         dU4y=dU(i,f,ie,4,2)
         dU5y=dU(i,f,ie,5,2)
         dU1z=dU(i,f,ie,1,3)
         dU2z=dU(i,f,ie,2,3)
         dU3z=dU(i,f,ie,3,3)
         dU4z=dU(i,f,ie,4,3)
         dU5z=dU(i,f,ie,5,3)
         rho   =wstate(i,f,ie,jrho)
         cv    =wstate(i,f,ie,jcvf)/rho
         lambda=wstate(i,f,ie,jlamf)
         mu    =wstate(i,f,ie,jmuf)
         K     =0.0 ! ADIABATIC HARDCODING
         u1    =wstate(i,f,ie,jux)
         u2    =wstate(i,f,ie,juy)
         u3    =wstate(i,f,ie,juz)
         E     =wstate(i,f,ie,ju5)/rho
         lambdamu=lambda+mu
         kmcvmu=K-cv*mu
         flux(i)=
     >(K*(dU5z-E*dU1z)+cv*u3*(lambda*dU4z+2*mu*dU4z+lambda*dU3y+lambda
     1   *dU2x)-K*u3*dU4z+cv*mu*u2*(dU4y+dU3z)+cv*mu*u1*(dU4x+dU2z)-
     2   K*u2*dU3z-K*u1*dU2z-cv*(lambda+2*mu)*u3**2*dU1z+K*u3**2*dU1z+
     3   K*u2**2*dU1z-cv*mu*u2**2*dU1z+K*u1**2*dU1z-cv*mu*u1**2*dU1z-c
     4   _v*(lambda+mu)*u2*u3*dU1y-cv*(lambda+mu)*u1*u3*dU1x)/(cv*rho)
      enddo
      return
      end
