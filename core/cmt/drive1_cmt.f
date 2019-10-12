C> @file drive1_cmt.f high-level driver for CMT-nek
C> \defgroup convhvol Volume integral for inviscid fluxes
C> \defgroup bcond Surface integrals due to boundary conditions
C> \defgroup diffhvol Volume integral for viscous fluxes
C> \defgroup vfjac Jacobians for viscous fluxes
C> \defgroup isurf Inviscid surface terms
C> \defgroup vsurf Viscous surface terms
C> \defgroup faceops utility functions for manipulating face data
C> Branch from subroutine nek_advance in core/drive1.f
C> Advance CMT-nek one time step within nek5000 time loop
      subroutine cmt_nek_advance
c     Solve the Euler equations

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
      include 'SOLN'
      include 'GEOM'
      include 'CTIMER'
      include 'CMTDATA'
      include 'CMTTIMERS'
      
      integer e,eq
      character*32 dumchars

      ftime_dum = dnekclock()
      nxyz1=lx1*ly1*lz1
      n = nxyz1*lelt*toteq
      nfldpart = ldim*npart

      if(istep.eq.1) then
         call cmt_ics
         if (ifrestart) then
            time_cmt=time
         else
            time_cmt=0.0 !time !0.0 ! until we can get settime to behave
         endif
         call cmt_flow_ics
         call init_cmt_timers
         call cmtchk ! need more ifdefs
         call compute_mesh_h(meshh,xm1,ym1,zm1)
         call compute_grid_h(gridh,xm1,ym1,zm1)
         call compute_primitive_vars(1) ! get good mu
         call limiter
         call entropy_viscosity      ! for high diffno
         call compute_transport_props! at t=0

      endif

      nstage = 3
      do stage=1,nstage
         if (stage.eq.1) call copy(res3(1,1,1,1,1),U(1,1,1,1,1),n)

         rhst_dum = dnekclock()
         call compute_rhs_and_dt
         rhst = rhst + dnekclock() - rhst_dum
c particle equations of motion are solved (also includes forcing)
c In future this subroutine may compute the back effect of particles
c on the fluid and suitably modify the residue computed by 
c compute_rhs_dt for the 5 conserved variables
         call usr_particles_solver

! JH111815 soon....
! JH082316 someday...maybe?
!        do eq=1,toteq
!           call fbinvert(res1(1,1,1,1,eq))
!        enddo

         do e=1,nelt
            do eq=1,toteq
            do i=1,nxyz1
! JH071218 res1 is premultiplied by B^{-1}
               u(i,1,1,eq,e) = tcoef(1,stage)*res3(i,1,1,eq,e)+
     >                         tcoef(2,stage)*u(i,1,1,eq,e)-
     >                         tcoef(3,stage)*res1(i,1,1,e,eq)
            enddo
            enddo
         enddo ! nelt
!-----------------------------------------------------------------------
! JH080918 Now with solution limiters of Zhang & Shu (2010)
!                                    and   Lv & Ihme (2015) 
!          Also, FINALLY rewritten to consider solution at the
!          END OF RK STAGES AND END OF TIME STEP AS THE SOLUTION OF INTEREST
! JH081018 OK I can't do that for some reason. CHECK SOLN COMMONS BETWEEN
!          cmt_nek_advance and istep=istep+1
!-----------------------------------------------------------------------
!        call compute_primitive_vars(0)
!        call limiter
!        call compute_primitive_vars(1)

      enddo ! RK stage loop

      ftime = ftime + dnekclock() - ftime_dum

!     call print_cmt_timers ! NOT NOW

 101  format(4(2x,e18.9))
      return
      end

c-----------------------------------------------------------------------

C> Compute right-hand-side of the semidiscrete conservation law
C> Store it in res1
      subroutine compute_rhs_and_dt
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
      include 'CMTDATA'
      include 'CTIMER'

! hdsize needs to be big enough for 15 full fields
      common /CMTSURFLX/ fatface(heresize),graduf(hdsize)
      real graduf

      integer e,eq
      real wkj(lx1+lxd)
      character*32  dumchars

      nxyz=lx1*ly1*lz1

      call compute_mesh_h(meshh,xm1,ym1,zm1)
      call compute_grid_h(gridh,xm1,ym1,zm1)

      call cmt_metrics(istep)

!     filter the conservative variables before start of each
!     time step
!     if(IFFLTR)  call filter_cmtvar(IFCNTFILT)
!        primitive vars = rho, u, v, w, p, T, phi_g

      call compute_primitive_vars(0)
!! JH090518 Shock detector is not ready for prime time. Lean on EVM for
!!          sane default 
!!     if (stage.eq.1)
!!    >call shock_detector(t(1,1,1,1,5),vtrans(1,1,1,1,jrho),scrent)
      call limiter
      call compute_primitive_vars(1)

!-----------------------------------------------------------------------
! JH072914 We can really only proceed with dt once we have current
!          primitive variables. Only then can we compute CFL and/or dt.
!-----------------------------------------------------------------------
      if(stage.eq.1) then
         call setdtcmt
         call set_tstep_coef
!-----------------------------------------------------------------------
! JH081018 a whole bunch of this stuff should really be done AFTER the
!          RK loop at the END of the time step, but I lose custody
!          of commons in SOLN between cmt_nek_advance and the rest of
!          the time loop.
         call copy(t(1,1,1,1,2),vtrans(1,1,1,1,jrho),nxyz*nelt)
!! JH070119 Tait mixture model extension. Need T(:,2) for mass fraction
!!          of one of the two species. put mixture density (for
!!          post-processing only) into T(:,4)           
!c JB080119 more species change mix density, T(:,5) for 3 species
!         call copy(t(1,1,1,1,4),vtrans(1,1,1,1,jrho),nxyz*nelt)
!c        call copy(t(1,1,1,1,5),vtrans(1,1,1,1,jrho),nxyz*nelt)
         call cmtchk

         if (mod(istep,iostep).eq.0.or.istep.eq.1)then ! migrate to iostep2
            call out_fld_nek ! solution checkpoint for restart
! T2 S1 rho
! T3 S2 wave visc
! T4 S3 epsebdg
            call outpost2(vx,vy,vz,pr,t,ldimt,'CMT')
            call mass_balance(if3d)
! dump out particle information. 
#ifdef LPM
            call lpm_usr_particles_io(istep)
#endif
         end if
      endif

      call entropy_viscosity ! accessed through uservp. computes
                             ! entropy residual and max wave speed
      call compute_transport_props ! everything inside rk stage

      ntot = lx1*ly1*lz1*lelt*toteq
      call rzero(res1,ntot)
      call rzero(fatface,heresize)

!     !Total_eqs = 5 (we will set this up so that it can be a user 
!     !defined value. 5 will be its default value)
!     !eq = 1 -------- Mass balance
!     !eq = 2 -------- x  momentum 
!     !eq = 3 -------- y  momentum 
!     !eq = 4 -------- z  momentum 
!     !eq = 5 -------- Energy Equation 

C> Restrict via \f$\mathbf{E}\f$ to get primitive and conserved variables
C> on interior faces \f$\mathbf{U}^-\f$ and neighbor faces
C> \f$\mathbf{U}^+\f$; store in CMTSURFLX

      call cmt_usrsurf

C> res1+=\f$\oint \mathbf{H}^{c\ast}\cdot\mathbf{n}dA\f$ on face points
      nstate=nqq
      nfq=lx1*lz1*2*ldim*nelt
      iwm =1
      iwp =iwm+nstate*nfq
! KEPEC hopefully streamlined
!     iflx=iwp+nfq ! ok this needs to be segregated and CMTSURFLX redeclared.
! KEPEC done naively via duplicated
      iflx=iwp+nstate*nfq ! ok this needs to be segregated and CMTSURFLX redeclared.
                   ! W+ depends on flux function and may not always be 1
      do eq=1,toteq
         ieq=(eq-1)*ndg_face+iflx
         call surface_integral_full(res1(1,1,1,1,eq),fatface(ieq))
      enddo
      dumchars='after_inviscid'
!     call dumpresidue(dumchars,999)

!     call gtu_wrapper(fatface) ! for penalty methods. not yet

! JH091319 overwrite wminus with -[[U]] for viscous terms BR1
      call fillujumpu

      do e=1,nelt
!-----------------------------------------------------------------------
! JH082216 Since the dawn of CMT-nek we have called this particular loop
!***********************************************************************
!*         "THE" ELEMENT LOOP                                          *
!***********************************************************************
!          since it does several operations, mostly for volume integrals,
!          for all equations, one element at a time. If we get memory
!          under control and GPUs really need to act on gigabytes all
!          at once, then this and its dependents can still have their
!          loop order flipped and things like totalh declared for
!          15 full fields or more.
! JH060418 totalh is now 15 elements. interchanged with equation loop
!-----------------------------------------------------------------------
! Get user defined forcing from userf defined in usr file
         call cmtusrf(e)
         call compute_gradients_contra(e) ! gradU
         i_cvars=1
         do eq=1,toteq
            call br1auxflux(e,gradu(1,1,eq),fatface(i_cvars)) ! SEE HEAT.USR
            i_cvars=i_cvars+nfq
         enddo
         call convective_cmt(e)        ! convh & totalh -> res1
         do eq=1,toteq
            call    viscous_cmt(e,eq) ! diffh -> half_iku_cmt -> res1
                                             !       |
                                             !       -> diffh2graduf
! Compute the forcing term in each of the 5 eqs
            if (1.eq.2) then
               call compute_forcing(e,eq)
            endif
         enddo
      enddo
      dumchars='after_elm'
!     call dumpresidue(dumchars,999)

! COMPARE TO HEAT.USR. IGU should be the same
C> res1+=\f$\int_{\Gamma} \{\{\mathbf{A}\nabla \mathbf{U}\}\} \cdot \left[v\right] dA\f$
!     call igu_cmt(flux(iwp),graduf,flux(iwm))
      call br1primary(fatface(iwm),graduf)
      do eq=1,toteq
         ieq=(eq-1)*ndg_face+iwm
!Finally add viscous surface flux functions of derivatives to res1.
         call surface_integral_full(res1(1,1,1,1,eq),fatface(ieq))
      enddo

! one last
      do eq=1,toteq
         call col2(res1(1,1,1,1,eq),jacmi,nelt*lx1*ly1*lz1)
      enddo

      dumchars='end_of_rhs'
!     call dumpresidue(dumchars,999)

      return
      end
!-----------------------------------------------------------------------
C> Compute coefficients for Runge-Kutta stages \cite{TVDRK}
      subroutine set_tstep_coef

      real tcoef(3,3),dt_cmt,time_cmt
      COMMON /TIMESTEPCOEF/ tcoef,dt_cmt,time_cmt

      tcoef(1,1) = 0.0
      tcoef(2,1) = 1.0 
      tcoef(3,1) = dt_cmt
      tcoef(1,2) = 3.0/4.0
      tcoef(2,2) = 1.0/4.0 
      tcoef(3,2) = dt_cmt/4.0 
      tcoef(1,3) = 1.0/3.0
      tcoef(2,3) = 2.0/3.0 
      tcoef(3,3) = dt_cmt*2.0/3.0 

      return
      end
!-----------------------------------------------------------------------

      subroutine cmt_flow_ics
      include 'SIZE'
      include 'CMTDATA'
      include 'SOLN'

      integer e
      nxyz1 = lx1*ly1*lz1
      n     = nxyz1*lelt*toteq
      if (ifrestart)then
         do e=1,nelt
            call copy(U(1,1,1,2,e),vx(1,1,1,e),nxyz1) 
            call copy(U(1,1,1,3,e),vy(1,1,1,e),nxyz1) 
            call copy(U(1,1,1,4,e),vz(1,1,1,e),nxyz1) 
            call copy(U(1,1,1,5,e),t(1,1,1,e,1),nxyz1) 
            call copy(U(1,1,1,1,e),pr(1,1,1,e),nxyz1) 
         enddo
         call copy(tlag(1,1,1,1,1,2),t(1,1,1,1,2),nxyz1*nelt) ! s_{n-1}
         call copy(tlag(1,1,1,1,2,1),t(1,1,1,1,3),nxyz1*nelt) ! s_n
      endif
      call rzero(res1,n)
!     call copy(res2,t(1,1,1,1,5),n) ! art visc hardcoding. old entropy resid
      call rzero(res2,n) ! Actually,...
      return
      end
!-----------------------------------------------------------------------

      subroutine print_cmt_timers
      include 'SIZE'
      include 'CMTTIMERS'
      include 'TSTEP'
      include 'PARALLEL'

c we need our own IO features. Until then we use the default nek routines
      if ((mod(istep,flio_freq).eq.0.and.istep.gt.0)
     $                               .or.istep.eq.nstep)then
         dmtime1 = ftime/istep
         dtime_ = glsum(dmtime1,1)
         if(nio.eq.0) write(6,*) 'fluid rhs compute time(Avg)  '
     $               ,dtime_/np
      endif
      return 
      end
!-----------------------------------------------------------------------

      subroutine init_cmt_timers
      include 'CMTTIMERS'

      rhst    = 0.00
      ftime   = 0.00

      return
      end
!-----------------------------------------------------------------------
