C> @file surface_fluxes.f Routines for surface terms on RHS.

C> \ingroup isurf
C> @{
C> Restrict and copy face data and compute inviscid numerical flux 
C> \f$\oint \mathbf{H}^{c\ast}\cdot\mathbf{n}dA\f$ on face points
      subroutine fluxes_full_field(fstab,parameter_vector,fsharp)
!-----------------------------------------------------------------------
! JH060314 First, compute face fluxes now that we have the primitive variables
! JH091514 renamed from "surface_fluxes_inviscid" since it handles all fluxes
!          that we compute from variables stored for the whole field (as
!          opposed to one element at a time).
! JH070918 redone for two-point fluxes 
!-----------------------------------------------------------------------
      include 'SIZE'
      include 'DG'
      include 'SOLN'
      include 'CMTDATA'
      include 'INPUT'

      integer lfq,heresize,hdsize
      parameter (lfq=lx1*lz1*2*ldim*lelt,
     >                   heresize=nqq*3*lfq,! guarantees transpose of Q+ fits
     >                   hdsize=toteq*3*lfq) ! might not need ldim
! JH070214 OK getting different answers whether or not the variables are
!          declared locally or in common blocks. switching to a different
!          method of memory management that is more transparent to me.
      common /CMTSURFLX/ fatface(heresize),notyet(hdsize)
      real fatface,notyet
      external fstab,parameter_vector,fsharp
      integer eq
      character*32 cname
      nfq=lx1*lz1*2*ldim*nelt
      nstate = nqq
! where different things live
      iwm =1
      iwp =iwm+nstate*nfq
      iflx=iwp+nstate*nfq

      call rzero(fatface(iflx),nfq*toteq)

! fill parameter vector z for two-point flux applied at faces
! start with primitive variables at faces

! just multipy by {{phi}}
! until we can verify correct multiphase two-point fluxes
!     call faceu(1,fatface(iwm+nfq*(irho-1)))
      call fillq(iux, vx,    fatface(iwm),fatface(iwp))
      call fillq(iuy, vy,    fatface(iwm),fatface(iwp))
      call fillq(iuz, vz,    fatface(iwm),fatface(iwp))
      call fillq(iph, phig,  fatface(iwm),fatface(iwp))

! stabilization first.
      call fillq(isnd,csound,fatface(iwm),fatface(iwp))
      call fstab(fatface(iwm),fatface(iwp),fatface(iflx),nstate)

! need total energy, not internal
      i_cvars=(iu5-1)*nfq+1
!     call faceu(toteq,fatface(i_cvars)) ! should persist from llf_euler
      call invcol2(fatface(i_cvars),fatface(iwm+nfq*(iph-1)),nfq)


      call fillq(irho,vtrans,fatface(iwm),fatface(iwp))
      call fillq(ipr, pr,    fatface(iwm),fatface(iwp))

! q- -> z-. Kennedy-Gruber, Pirozzoli, and most energy-
!           conserving fluxes have z=q, so I just divide total energy by
!           U1 here since Kennedy-Gruber needs E
      call parameter_vector(fatface(iwm),nfq,nstate)

! z- -> z^, which is {{z}} for Kennedy-Gruber, Pirozzoli, and some parts
!           of other energy-conserving fluxes.
      call dg_face_avg(fatface(iwm),nfq,nstate,dg_hndl)

! z^ -> F#. Some parameter-vector stuff can go here too as long as it's all
!           local to a given element.
      call fsharp(fatface(iwm),fatface(iflx),nstate,toteq)

!     i_cvars=iwm!(iu1-1)*nfq+1
!     do eq=1,toteq
!        call faceu(eq,fatface(i_cvars))
!        i_cvars=i_cvars+nfq
!     enddo

C> @}
! Now do all fluxes for all boundaries, both F# and stabilized
      call InviscidBC(fatface(iflx))

      return
      end

!-----------------------------------------------------------------------

      subroutine faceu(ivar,yourface)
! get faces of conserved variables stored contiguously
      include 'SIZE'
      include 'CMTDATA'
      include 'DG'
      integer e
      real yourface(lx1,lz1,2*ldim,nelt)

      do e=1,nelt
         call full2face_cmt(1,lx1,ly1,lz1,iface_flux(1,e),
     >                      yourface(1,1,1,e),u(1,1,1,ivar,e))
      enddo

      return
      end

!-----------------------------------------------------------------------

      subroutine fillq(ivar,field,wminus,yourface)
      include 'SIZE'
      include 'DG'

      integer ivar! intent(in)
      real field(lx1,ly1,lz1,nelt)! intent(in)
!     real, intent(out)wminus(7,lx1*lz1*2*ldim*nelt) ! gs_op no worky
      real wminus(lx1*lz1*2*ldim*nelt,*)! intent(out)
      real yourface(lx1,lz1,2*ldim,*)
      integer e,f

      nxz  =lx1*lz1
      nface=2*ldim

      call full2face_cmt(nelt,lx1,ly1,lz1,iface_flux,yourface,field)

      do i=1,ndg_face
         wminus(i,ivar)=yourface(i,1,1,1)
      enddo

      return
      end

!-----------------------------------------------------------------------

      subroutine dg_face_avg(mine,nf,nstate,handle)

! JH110818 A lot of entropy-stable fluxes have a product of averages.
!          mine: starting parameter vector of nstate quantities,
!                quantity outermost, on piles of faces with nf points

      integer handle,nf,nstate ! intent(in)
      real mine(*)

      ntot=nf*nstate
      call cmult(mine,0.5,ntot)
!-----------------------------------------------------------------------
! operation flag is second-to-last arg, an integer
!                                                1 ==> +
      call fgslib_gs_op_fields(handle,mine,nf,nstate,1,1,0)
      return
      end

!-----------------------------------------------------------------------

      subroutine face_state_commo(mine,yours,nf,nstate,handle)

! JH060414 if we ever want to be more intelligent about who gets what,
!          who gives what and who does what, this is the place where all
!          that is done. At the very least, gs_op may need the transpose
!          flag set to 1. Who knows. Everybody duplicates everything for
!          now.
! JH070714 figure out gs_op_fields, many, vec, whatever (and the
!          corresponding setup) to get this done for the transposed
!          ordering of state variables. I want variable innermost, not
!          grid point.

      integer handle,nf,nstate ! intent(in)
      real yours(*),mine(*)

      ntot=nf*nstate
      call copy(yours,mine,ntot)
!-----------------------------------------------------------------------
! operation flag is second-to-last arg, an integer
!                                                1 ==> +
      call fgslib_gs_op_fields(handle,yours,nf,nstate,1,1,0)
      call sub2 (yours,mine,ntot)
      return
      end

!-----------------------------------------------------------------------

      subroutine face_flux_commo(flux1,flux2,nf,neq,handle)
! JH060514 asymmetric transposed gs_op, gs_unique magic may be needed if
!          we ever decide to avoid redundancy. For now, this routine
!          doesn't need to do anything.
      integer ntot,handle
      real flux1(*),flux2(*)
! JH061814 It doesn't need to do anything, but a sanity check would be
!          wise.
      return
      end

!-------------------------------------------------------------------------------

      subroutine sequential_flux(flux,wminus,wplus,uminus,uplus,
     >                           jaminus,japlus,
     >                           fluxfunction,nstate,npt)
! Calls two-point flux functions one point at a time for npt points for which both
! points are given. Mostly intended to allow quantity-innermost volume flux
! functions to be used where needed for boundary points too, after *bc routines
! provide Dirichlet ``rind'' states in wplus and uplus.
      include 'SIZE'
      real flux(toteq,npt),wminus(nstate,npt),wplus(nstate,npt),
     >   jaminus(3,npt),japlus(3,npt),uminus(toteq,npt),uplus(toteq,npt)
      external fluxfunction

      do i=1,npt
         call fluxfunction(flux(1,i),uminus(1,i),uplus(1,i),wminus(1,i),
     >                     wplus(1,i),jaminus(1,i),japlus(1,i))
      enddo

      return
      end

!-------------------------------------------------------------------------------
! OBSOLETE ROUTINES
!-------------------------------------------------------------------------------

      subroutine InviscidFlux(wminus,wplus,flux,nstate,nflux)
!-------------------------------------------------------------------------------
! JH091514 A fading copy of RFLU_ModAUSM.F90 from RocFlu
!-------------------------------------------------------------------------------

!#ifdef SPEC
!      USE ModSpecies, ONLY: t_spec_type
!#endif
      include 'SIZE'
      include 'INPUT' ! do we need this?
      include 'GEOM' ! for unx
      include 'CMTDATA' ! do we need this without outflsub?
      include 'TSTEP' ! for ifield?
      include 'DG'

! ==============================================================================
! Arguments
! ==============================================================================
      integer nstate,nflux
      real wminus(lx1*lz1,2*ldim,nelt,nstate),
     >     wplus(lx1*lz1,2*ldim,nelt,nstate),
     >     flux(lx1*lz1,2*ldim,nelt,nflux)

! ==============================================================================
! Locals
! ==============================================================================

      integer e,f,fdim,i,k,nxz,nface
      parameter (lfd=lxd*lzd)
! JH111815 legacy rocflu names.
!
! nx,ny,nz : outward facing unit normal components
! fs       : face speed. zero until we have moving grid
! jaco_c   : fdim-D GLL grid Jacobian
! nm       : jaco_c, fine grid
!
! State on the interior (-, "left") side of the face
! rl       : density
! ul,vl,wl : velocity
! tl       : temperature
! al       : sound speed
! pl       : pressure, then phi
! el      : rho*cp
! State on the exterior (+, "right") side of the face
! rr       : density
! ur,vr,wr : velocity
! tr       : temperature
! ar       : sound speed
! pr       : pressure
! er      : rho*cp

      COMMON /SCRNS/ nx(lfd), ny(lfd), nz(lfd), rl(lfd), ul(lfd),
     >               vl(lfd), wl(lfd), pl(lfd), tl(lfd), al(lfd),
     >               el(lfd),rr(lfd), ur(lfd), vr(lfd), wr(lfd),
     >               pr(lfd),tr(lfd), ar(lfd),er(lfd),phl(lfd),fs(lfd),
     >               jaco_f(lfd),flx(lfd,toteq),jaco_c(lx1*lz1)
      real nx, ny, nz, rl, ul, vl, wl, pl, tl, al, el, rr, ur, vr, wr,
     >                pr,tr, ar,er,phl,fs,jaco_f,flx,jaco_c

!     REAL vf(3)
      real nTol
      character*132 deathmessage
      common /nekcb/ cb
      character*3 cb

      nTol = 1.0E-14

      fdim=ldim-1
      nface = 2*ldim
      nxz   = lx1*lz1
      nxzd  = lxd*lzd
      ifield= 1 ! You need to figure out the best way of dealing with
                ! this variable

!     if (outflsub)then
!        call maxMachnumber
!     endif
      do e=1,nelt
      do f=1,nface
! diagnostic
           call facind (kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,f)

! JH021717 Finally corrected BC wrongheadedness. Riemann solver handles
!          all fluxes with even the slightest Dirichlet interpretation,
!          and BC fill wplus sanely for the Riemann solver to provide
!          a good flux for weak BC.
! JH111715 now with dealiased surface integrals. I am too lazy to write
!          something better

         if (lxd.gt.lx1) then
            call map_faced(nx,unx(1,1,f,e),lx1,lxd,fdim,0)
            call map_faced(ny,uny(1,1,f,e),lx1,lxd,fdim,0)
            call map_faced(nz,unz(1,1,f,e),lx1,lxd,fdim,0)

            call map_faced(rl,wminus(1,f,e,irho),lx1,lxd,fdim,0)
            call map_faced(ul,wminus(1,f,e,iux),lx1,lxd,fdim,0)
            call map_faced(vl,wminus(1,f,e,iuy),lx1,lxd,fdim,0)
            call map_faced(wl,wminus(1,f,e,iuz),lx1,lxd,fdim,0)
            call map_faced(pl,wminus(1,f,e,ipr),lx1,lxd,fdim,0)
            call map_faced(tl,wminus(1,f,e,ithm),lx1,lxd,fdim,0)
            call map_faced(al,wminus(1,f,e,isnd),lx1,lxd,fdim,0)
            call map_faced(el,wminus(1,f,e,icpf),lx1,lxd,fdim,0)

            call map_faced(rr,wplus(1,f,e,irho),lx1,lxd,fdim,0)
            call map_faced(ur,wplus(1,f,e,iux),lx1,lxd,fdim,0)
            call map_faced(vr,wplus(1,f,e,iuy),lx1,lxd,fdim,0)
            call map_faced(wr,wplus(1,f,e,iuz),lx1,lxd,fdim,0)
            call map_faced(pr,wplus(1,f,e,ipr),lx1,lxd,fdim,0)
            call map_faced(tr,wplus(1,f,e,ithm),lx1,lxd,fdim,0)
            call map_faced(ar,wplus(1,f,e,isnd),lx1,lxd,fdim,0)
            call map_faced(er,wplus(1,f,e,icpf),lx1,lxd,fdim,0)

            call map_faced(phl,wminus(1,f,e,iph),lx1,lxd,fdim,0)

            call invcol3(jaco_c,area(1,1,f,e),w2m1,nxz)
            call map_faced(jaco_f,jaco_c,lx1,lxd,fdim,0) 
            call col2(jaco_f,wghtf,nxzd)
         else

            call copy(nx,unx(1,1,f,e),nxz)
            call copy(ny,uny(1,1,f,e),nxz)
            call copy(nz,unz(1,1,f,e),nxz)

            call copy(rl,wminus(1,f,e,irho),nxz)
            call copy(ul,wminus(1,f,e,iux),nxz)
            call copy(vl,wminus(1,f,e,iuy),nxz)
            call copy(wl,wminus(1,f,e,iuz),nxz)
            call copy(pl,wminus(1,f,e,ipr),nxz)
            call copy(tl,wminus(1,f,e,ithm),nxz)
            call copy(al,wminus(1,f,e,isnd),nxz)
            call copy(el,wminus(1,f,e,icpf),nxz)

            call copy(rr,wplus(1,f,e,irho),nxz)
            call copy(ur,wplus(1,f,e,iux),nxz)
            call copy(vr,wplus(1,f,e,iuy),nxz)
            call copy(wr,wplus(1,f,e,iuz),nxz)
            call copy(pr,wplus(1,f,e,ipr),nxz)
            call copy(tr,wplus(1,f,e,ithm),nxz)
            call copy(ar,wplus(1,f,e,isnd),nxz)
            call copy(er,wplus(1,f,e,icpf),nxz)

            call copy(phl,wminus(1,f,e,iph),nxz)

            call copy(jaco_f,jface(1,1,f,e),nxz) 
         endif
         call rzero(fs,nxzd) ! moving grid stuff later

         call AUSM_FluxFunction(nxzd,nx,ny,nz,jaco_f,fs,rl,ul,vl,wl,pl,
     >                          al,tl,rr,ur,vr,wr,pr,ar,tr,flx,el,er)

         do j=1,toteq
            call col2(flx(1,j),phl,nxzd)
         enddo

         if (lxd.gt.lx1) then
            do j=1,toteq
               call map_faced(flux(1,f,e,j),flx(1,j),lx1,lxd,fdim,1)
            enddo
         else
            do j=1,toteq
               call copy(flux(1,f,e,j),flx(1,j),nxz)
            enddo
         endif

      enddo
      enddo

      end

!-----------------------------------------------------------------------

      subroutine surface_integral_full(vol,flux)
! Integrate surface fluxes for an entire field. Add contribution of flux
! to volume according to add_face2full_cmt
      include 'SIZE'
      include 'GEOM'
      include 'DG'
      include 'CMTDATA'
      real vol(lx1*ly1*lz1*nelt),flux(*)
      integer e,f

! weak form until we get the time loop rewritten
!     onem=-1.0
!     ntot=lx1*lz1*2*ldim*nelt
!     call cmult(flux,onem,ntot)
! weak form until we get the time loop rewritten
      call add_face2full_cmt(nelt,lx1,ly1,lz1,iface_flux,vol,flux)

      return
      end

!-------------------------------------------------------------------------------

      subroutine diffh2graduf(e,eq,graduf)
! peels off diffusiveH into contiguous face storage via restriction operator R
! for now, stores {{gradU}} for igu
      include  'SIZE'
      include  'DG' ! iface
      include  'CMTDATA'
      include  'GEOM'
      integer e,eq
      real graduf(lx1*lz1*2*ldim,nelt,toteq)
      common /scrns/ hface(lx1*lz1,2*ldim)
     >              ,normal(lx1*ly1,2*ldim)
      real hface, normal

      integer f

      nf    = lx1*lz1*2*ldim*nelt
      nfaces=2*ldim
      nxz   =lx1*lz1
      nxzf  =nxz*nfaces
      nxyz  = lx1*ly1*lz1

      call rzero(graduf(1,e,eq),nxzf) !   . dot nhat -> overwrites beginning of flxscr
      do j =1,ldim
         if (j .eq. 1) call copy(normal,unx(1,1,1,e),nxzf)
         if (j .eq. 2) call copy(normal,uny(1,1,1,e),nxzf)
         if (j .eq. 3) call copy(normal,unz(1,1,1,e),nxzf)
         call full2face_cmt(1,lx1,ly1,lz1,iface_flux,hface,diffh(1,j)) 
         call addcol3(graduf(1,e,eq),hface,normal,nxzf)
      enddo
      call col2(graduf(1,e,eq),area(1,1,1,e),nxzf)

      return
      end

!-----------------------------------------------------------------------

      subroutine igu_cmt(flxscr,gdudxk,wminus)
! Hij^{d*}
      include 'SIZE'
      include 'CMTDATA'
      include 'DG'

      real flxscr(lx1*lz1*2*ldim*nelt,toteq)
      real gdudxk(lx1*lz1*2*ldim,nelt,toteq)
      real wminus(lx1*lz1,2*ldim,nelt,nqq)
      real const
      integer e,eq,f

      nxz = lx1*lz1
      nfaces=2*ldim
      nxzf=nxz*nfaces
      nfq =lx1*lz1*nfaces*nelt
      ntot=nfq*toteq

      call copy (flxscr,gdudxk,ntot) ! save AgradU.n
      const = 0.5
      call cmult(gdudxk,const,ntot)
!-----------------------------------------------------------------------
! supa huge gs_op to get {{AgradU}}
! operation flag is second-to-last arg, an integer
!                                                   1 ==> +
      call fgslib_gs_op_fields(dg_hndl,gdudxk,nfq,toteq,1,1,0)
!-----------------------------------------------------------------------
      call sub2  (flxscr,gdudxk,ntot) ! overwrite flxscr with
                                      !           -
                                      ! (AgradU.n)  - {{AgradU.n}}
! [v] changes character on boundaries, so we actually need
! 1. (AgradU.n)- on Dirichlet boundaries
      call igu_dirichlet(flxscr,gdudxk)
! 2. (Fbc.n)- on Neumann boundaries
      call bcflux(flxscr,gdudxk,wminus)
      call chsign(flxscr,ntot) ! needs to change with sign changes

      return
      end

!-----------------------------------------------------------------------

      subroutine igu_dirichlet(flux,agradu)
! Acts on ALL boundary faces because I'm lazy. SET NEUMANN BC AFTER THIS
! CALL. BCFLUX IS PICKIER ABOUT THE BOUNDARY FACES IT ACTS ON.
      include 'SIZE'
      include 'TOTAL'
      integer e,eq,f
      real flux(lx1*lz1,2*ldim,nelt,toteq)
      real agradu(lx1*lz1,2*ldim,nelt,toteq)
      character*3 cb2

      nxz=lx1*lz1
      nfaces=2*ldim

      ifield=1
      do e=1,nelt
         do f=1,nfaces
            cb2=cbc(f, e, ifield)
            if (cb2.ne.'E  '.and.cb2.ne.'P  ') then ! cbc bndy.
! all Dirichlet conditions result in IGU being
! strictly one-sided, so we undo 0.5*QQT
! UNDER THE ASSUMPTIONS THAT
! 1. agradu's actual argument is really gdudxk AND
! 2. IT HAS ALREADY BEEN MULTIPLIED BY 0.5
! 3. gs_op has not changed it at all.
! overwriting flux with it and and multiplying it 2.0 should do the trick
               do eq=1,toteq
!                  call copy(flux(1,f,e,eq),agradu(1,f,e,eq),nxz)
!! in fact, that copy should not be necessary at all. TEST WITHOUT IT
!                  call cmult(flux(1,f,e,eq),2.0,nxz)
! JH112216 This may be better because agradu (without the factor of 1/2) is
!          needed for some Neumann conditions (like adiabatic walls)
                   call cmult(agradu(1,f,e,eq),2.0,nxz)
                   call copy(flux(1,f,e,eq),agradu(1,f,e,eq),nxz)
               enddo
            endif
         enddo
      enddo

      return
      end

!-----------------------------------------------------------------------

      subroutine strong_sfc_flux(flux,vflx,e,eq)
! JH052818 dedicated routine for stripping faces from an array vflx of
!          physical fluxes on GLL faces and storing them in the flux
!          pile of faces for a call to surface_integral_full
      include 'SIZE'
      include 'INPUT'
      include 'GEOM' ! for unx
      include 'DG'
      
      parameter (lf=lx1*lz1*2*ldim)
      COMMON /SCRNS/ yourface(lf),normal(lf)
      real yourface,normal
      real flux(lx1*lz1,2*ldim,nelt,toteq) ! out
      real vflx(lx1*ly1*lz1,ldim)          ! in
      integer e,eq
      integer f

      nxz =lx1*lz1
      nxzf=nxz*2*ldim
      nfaces=2*ldim

      call rzero(flux(1,1,e,eq),nxzf)

      do j=1,ldim
         if (j .eq. 1) call copy(normal,unx(1,1,1,e),nxzf)
         if (j .eq. 2) call copy(normal,uny(1,1,1,e),nxzf)
         if (j .eq. 3) call copy(normal,unz(1,1,1,e),nxzf)
         call full2face_cmt(1,lx1,ly1,lz1,iface_flux(1,e),yourface,
     >                      vflx(1,j))
         call col2(yourface,normal,nxzf)
!        call add2(flux(1,1,e,eq),yourface,nxzf) ! needed for +sign in RK loop
         call sub2(flux(1,1,e,eq),yourface,nxzf)
      enddo

!     call col2(flux(1,1,e,eq),area(1,1,1,e),nxzf)
      do f=1,nfaces
         call col2(flux(1,f,e,eq),w2m1,nxz)
      enddo

      return
      end

!-----------------------------------------------------------------------
! JUNKYARD
!-----------------------------------------------------------------------

      subroutine fluxes_full_field_old
!-----------------------------------------------------------------------
! JH060314 First, compute face fluxes now that we have the primitive variables
! JH091514 renamed from "surface_fluxes_inviscid" since it handles all fluxes
!          that we compute from variables stored for the whole field (as
!          opposed to one element at a time).
! JH070918 redone for two-point fluxes 
!-----------------------------------------------------------------------
      include 'SIZE'
      include 'DG'
      include 'SOLN'
      include 'CMTDATA'
      include 'INPUT'

      integer lfq,heresize,hdsize
      parameter (lfq=lx1*lz1*2*ldim*lelt,
     >                   heresize=nqq*3*lfq,! guarantees transpose of Q+ fits
     >                   hdsize=toteq*3*lfq) ! might not need ldim
! JH070214 OK getting different answers whether or not the variables are
!          declared locally or in common blocks. switching to a different
!          method of memory management that is more transparent to me.
      common /CMTSURFLX/ fatface(heresize),notyet(hdsize)
      real fatface,notyet
      integer eq
      character*32 cname

      nstate = nqq
      nfq=lx1*lz1*2*ldim*nelt
! where different things live
      iwm =1
      iwp =iwm+nstate*nfq
      iflx=iwp+nstate*nfq
 
      call fillq(irho,vtrans,fatface(iwm),fatface(iwp))
      call fillq(iux, vx,    fatface(iwm),fatface(iwp))
      call fillq(iuy, vy,    fatface(iwm),fatface(iwp))
      call fillq(iuz, vz,    fatface(iwm),fatface(iwp))
      call fillq(ipr, pr,    fatface(iwm),fatface(iwp))
      call fillq(ithm,t,     fatface(iwm),fatface(iwp))
      call fillq(isnd,csound,fatface(iwm),fatface(iwp))
      call fillq(iph, phig,  fatface(iwm),fatface(iwp))
      call fillq(icvf,vtrans(1,1,1,1,icv),fatface(iwm),fatface(iwp))
      call fillq(icpf,vtrans(1,1,1,1,icp),fatface(iwm),fatface(iwp))
      call fillq(imuf, vdiff(1,1,1,1,imu), fatface(iwm),fatface(iwp))
      call fillq(ikndf,vdiff(1,1,1,1,iknd),fatface(iwm),fatface(iwp))
      call fillq(ilamf,vdiff(1,1,1,1,ilam),fatface(iwm),fatface(iwp))

      i_cvars=(iu1-1)*nfq+1
      do eq=1,toteq
         call faceu(eq,fatface(i_cvars))
! JH080317 at least get the product rule right until we figure out how
!          we want the governing equations to look
         call invcol2(fatface(i_cvars),fatface(iwm+nfq*(iph-1)),nfq)
         i_cvars=i_cvars+nfq
      enddo
      call face_state_commo(fatface(iwm),fatface(iwp),nfq,nstate
     >                     ,dg_hndl)
      call InviscidBC(fatface(iwm),fatface(iwp),nstate)
!     call InviscidFlux(fatface(iwm),fatface(iwp),fatface(iflx)
      call InviscidFluxRot(fatface(iwm),fatface(iwp),fatface(iflx)
     >                 ,nstate,toteq)
      return
      end

!-----------------------------------------------------------------------

      subroutine InviscidFluxRot(wminus,wplus,flux,nstate,nflux)
      include 'SIZE'
      include 'INPUT' ! do we need this?
      include 'GEOM' ! for unx
      include 'CMTDATA' ! do we need this without outflsub?
      include 'TSTEP' ! for ifield?
      include 'DG'

! ==============================================================================
! Arguments
! ==============================================================================
      integer nstate,nflux
      real wminus(lx1*lz1,2*ldim,nelt,nstate),
     >     wplus(lx1*lz1,2*ldim,nelt,nstate),
     >     flux(lx1*lz1,2*ldim,nelt,nflux)

! ==============================================================================
! Locals
! ==============================================================================

      integer e,f,fdim,i,k,nxz,nface
      parameter (lfd=lxd*lzd)
! JH111815 legacy rocflu names.
!
! nx,ny,nz : outward facing unit normal components
! fs       : face speed. zero until we have moving grid
! jaco_c   : fdim-D GLL grid Jacobian
! nm       : jaco_c, fine grid
!
! State on the interior (-, "left") side of the face
! rl       : density
! ul,vl,wl : velocity
! tl       : temperature
! al       : sound speed
! pl       : pressure, then phi
! el      : rho*cp
! State on the exterior (+, "right") side of the face
! rr       : density
! ur,vr,wr : velocity
! tr       : temperature
! ar       : sound speed
! pr       : pressure
! er      : rho*cp

      COMMON /SCRNS/ nx(lfd), ny(lfd), nz(lfd), rl(lfd), ul(lfd),
     >               vl(lfd), wl(lfd), pl(lfd), tl(lfd), al(lfd),
     >               el(lfd),rr(lfd), ur(lfd), vr(lfd), wr(lfd),
     >               pr(lfd),tr(lfd), ar(lfd),er(lfd),phl(lfd),fs(lfd),
     >               jaco_f(lfd),flx(lfd,toteq),jaco_c(lx1*lz1)
      real nx, ny, nz, rl, ul, vl, wl, pl, tl, al, el, rr, ur, vr, wr,
     >                pr,tr, ar,er,phl,fs,jaco_f,flx,jaco_c

!     REAL vf(3)
      real nTol
      character*132 deathmessage
      common /nekcb/ cb
      character*3 cb

      nTol = 1.0E-14

      fdim=ldim-1
      nface = 2*ldim
      nxz   = lx1*lz1
      nxzd  = lxd*lzd
      ifield= 1 ! You need to figure out the best way of dealing with
                ! this variable

      do e=1,nelt
      do f=1,nface

! JH111218 Test 2-point fluxes in rotated form.

         call copy(nx,unx(1,1,f,e),nxz)
         call copy(ny,uny(1,1,f,e),nxz)
         call copy(nz,unz(1,1,f,e),nxz)

         call copy(rl,wminus(1,f,e,irho),nxz)

         if(if3d) then
            call vdot3(ul,wminus(1,f,e,iux),wminus(1,f,e,iuy),
     >                    wminus(1,f,e,iuz),nx,ny,nz,nxz)
            call vdot3(vl,wminus(1,f,e,iux),wminus(1,f,e,iuy),
     >wminus(1,f,e,iuz),t1x(1,1,f,e),t1y(1,1,f,e),t1z(1,1,f,e),nxz)
            call vdot3(wl,wminus(1,f,e,iux),wminus(1,f,e,iuy),
     >wminus(1,f,e,iuz),t2x(1,1,f,e),t2y(1,1,f,e),t2z(1,1,f,e),nxz)
            call vdot3(ur,wplus(1,f,e,iux),wplus(1,f,e,iuy),
     >                    wplus(1,f,e,iuz),nx,ny,nz,nxz)
            call vdot3(vr,wplus(1,f,e,iux),wplus(1,f,e,iuy),
     >wplus(1,f,e,iuz),t1x(1,1,f,e),t1y(1,1,f,e),t1z(1,1,f,e),nxz)
            call vdot3(wr,wplus(1,f,e,iux),wplus(1,f,e,iuy),
     >wplus(1,f,e,iuz),t2x(1,1,f,e),t2y(1,1,f,e),t2z(1,1,f,e),nxz)
         else
            call vdot2(ul,wminus(1,f,e,iux),wminus(1,f,e,iuy),
     >                    nx,ny,nxz)
            call vdot2(vl,wminus(1,f,e,iux),wminus(1,f,e,iuy),
     >                 t1x(1,1,f,e),t1y(1,1,f,e),nxz)
            call rzero(wl,nxz)
            call vdot2(ur,wplus(1,f,e,iux),wplus(1,f,e,iuy),
     >                    nx,ny,nxz)
            call vdot2(vr,wplus(1,f,e,iux),wplus(1,f,e,iuy),
     >                 t1x(1,1,f,e),t1y(1,1,f,e),nxz)
            call rzero(wr,nxz)
         endif

         call copy(pl,wminus(1,f,e,ipr),nxz)
         call copy(tl,wminus(1,f,e,ithm),nxz)
         call copy(al,wminus(1,f,e,isnd),nxz)
         call copy(el,wminus(1,f,e,icpf),nxz)

         call copy(rr,wplus(1,f,e,irho),nxz)
         call copy(pr,wplus(1,f,e,ipr),nxz)
         call copy(tr,wplus(1,f,e,ithm),nxz)
         call copy(ar,wplus(1,f,e,isnd),nxz)
         call copy(er,wplus(1,f,e,icpf),nxz)

         call copy(phl,wminus(1,f,e,iph),nxz)

         call copy(jaco_f,jface(1,1,f,e),nxz) 
         do i=1,nxz
            el(i)=el(i)+0.5*(ul(i)**2+vl(i)**2+wl(i)**2)
            er(i)=er(i)+0.5*(ur(i)**2+vr(i)**2+wr(i)**2)
         enddo

         call KGrotFluxFunction(nxzd,jaco_f,rl,ul,vl,wl,pl,
     >                          al,tl,rr,ur,vr,wr,pr,ar,tr,flx,el,er)

         do j=1,toteq
            call col2(flx(1,j),phl,nxzd)
         enddo
! rotate back
         call copy(flux(1,f,e,1),flx(1,1),nxz)
         if (if3d) then
            call vdot3(flux(1,f,e,2),flx(1,2),flx(1,3),flx(1,4),
     >                                 nx,t1x(1,1,f,e),t2x(1,1,f,e),nxz)
            call vdot3(flux(1,f,e,3),flx(1,2),flx(1,3),flx(1,4),
     >                                 ny,t1y(1,1,f,e),t2y(1,1,f,e),nxz)
            call vdot3(flux(1,f,e,4),flx(1,2),flx(1,3),flx(1,4),
     >                                 nz,t1z(1,1,f,e),t2z(1,1,f,e),nxz)
         else
            call vdot2(flux(1,f,e,2),flx(1,2),flx(1,3),
     >                                 nx,t1x(1,1,f,e),nxz)
            call vdot2(flux(1,f,e,3),flx(1,2),flx(1,3),
     >                                 ny,t1y(1,1,f,e),nxz)
         endif
         call copy(flux(1,f,e,5),flx(1,5),nxz)

      enddo
      enddo

      end
