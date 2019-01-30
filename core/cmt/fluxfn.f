C> @file fluxfn.f Riemann solvers, other rocflu miscellany and two-point fluxes
      SUBROUTINE KGrotFluxFunction(ntot,nm,rl,ul,vl,wl,pl,
     >                         al,tl,rr,ur,vr,wr,pr,ar,tr,flx,el,er)
! JH111218 sanity check of order of operations in two-point form of surface fluxes
!          contravariant (like in the papers) should match rotated (like in fluxo)
!          This is the rotated form of the flux. I don't want to rely on it because
!          the contravariant form is formally correct and far cheaper.
! ==============================================================================
! Arguments
! ==============================================================================
      integer ntot
      REAL al(ntot),ar(ntot),nm(ntot),
     >     pl(ntot),pr(ntot),rl(ntot),rr(ntot),ul(ntot),
     >     ur(ntot),vl(ntot),vr(ntot),wl(ntot),wr(ntot),el(ntot),
     >     er(ntot),tl(ntot),tr(ntot)! INTENT(IN) ::
      REAL flx(ntot,5)
! Locals
      real rav,uav(3),pav,eav,jav(3)

      do i=1,ntot
         rav=0.5*(rl(i)+rr(i))
         pav=0.5*(pl(i)+pr(i))
         uav(1)=0.5*(ul(i)+ur(i))
         uav(2)=0.5*(vl(i)+vr(i))
         uav(3)=0.5*(wl(i)+wr(i))
         eav=0.5*(el(i)+er(i))
         flx(i,1)=rav*uav(1)
         flx(i,2)=flx(i,1)*uav(1)+pav
         flx(i,3)=flx(i,1)*uav(2)
         flx(i,4)=flx(i,1)*uav(3)
         flx(i,5)=flx(i,1)*eav+pav*uav(1)
      enddo
      do j=1,5
         call col2(flx(1,j),nm,ntot)
      enddo
      return
      end

! Purpose: Compute convective fluxes using AUSM+ scheme.
!
! Description: None.
!
! Input: 
!   nx          x-component of face normal
!   ny          y-component of face normal
!   nz          z-component of face normal
!   nm          Magnitude of face normal
!   fs          Face speed
!   rl          Density of left state
!   ul          x-component of velocity of left state
!   vl          y-component of velocity of left state
!   wl          z-component of velocity of left state   
!   Hl		Total enthalpy of left state
!   al		Speed of sound of left state
!   pl          Pressure of left state
!   rr          Density of right state
!   ur          x-component of velocity of right state
!   vr          y-component of velocity of right state
!   wr          z-component of velocity of right state  
!   pr          Pressure of right state
!   Hr		Total enthalpy of right state
!   ar		Speed of sound of right state
!
! Output: 
!   flx         Fluxes
!   vf          Face velocities ! NOT USED IN CMT-NEK YET
!
! Notes: 
!   1. Liou M.-S., Progress towards an improved CFD method: AUSM+, AIAA Paper
!      95-1701, 1995
!   2. Do not use computation of face speed of sound which leads to exact 
!      capturing of isolated normal shock waves because of robustness problems
!      for unsteady flows and because that formulation is not applicable to 
!      anything but calorically and thermally perfect gases.
!
! ******************************************************************************

C> \ingroup isurf
C> @{
C> Computes inviscid numerical surface flux from AUSM+ Riemann solver
      SUBROUTINE AUSM_FluxFunction(ntot,nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,
     >                         al,tl,rr,ur,vr,wr,pr,ar,tr,flx,el,er)

!     IMPLICIT NONE ! HAHAHHAHHAHA
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************
      real     MixtJWL_Enthalpy
      external MixtJWL_Enthalpy

! ==============================================================================
! Arguments
! ==============================================================================
      integer ntot
      REAL al(ntot),ar(ntot),fs(ntot),nm(ntot),nx(ntot),ny(ntot),
     >     nz(ntot),pl(ntot),pr(ntot),rl(ntot),rr(ntot),ul(ntot),
     >     ur(ntot),vl(ntot),vr(ntot),wl(ntot),wr(ntot),el(ntot),
     >     er(ntot),tl(ntot),tr(ntot)! INTENT(IN) ::
      REAL flx(ntot,5)!,vf(3) ! INTENT(OUT) ::

! ==============================================================================
! Locals
! ==============================================================================

      REAL af,mf,mfa,mfm,mfp,ml,mla,mlp,mr,mra,mrm,pf,ql,qr,vml,vmr,
     >        wtl,wtr,Hl,Hr

! ******************************************************************************
! Start, compute face state
! ******************************************************************************

      do i=1,ntot
!        Change the Enthalpy 
         Hl = MixtJWL_Enthalpy(rl(i),pl(i),ul(i),vl(i),wl(i),el(i))
         Hr = MixtJWL_Enthalpy(rr(i),pr(i),ur(i),vr(i),wr(i),er(i))

         ql = ul(i)*nx(i) + vl(i)*ny(i) + wl(i)*nz(i) - fs(i)
         qr = ur(i)*nx(i) + vr(i)*ny(i) + wr(i)*nz(i) - fs(i)

         af = 0.5*(al(i)+ar(i)) ! NOTE not using original formulation, see note
         ml  = ql/af
         mla = ABS(ml)

         mr  = qr/af
         mra = ABS(mr)    

         IF ( mla .le. 1.0 ) THEN 
            mlp = 0.25*(ml+1.0)*(ml+1.0) + 0.125*(ml*ml-1.0)*(ml*ml-1.0)
            wtl = 0.25*(ml+1.0)*(ml+1.0)*(2.0-ml) +
     >            0.1875*ml*(ml*ml-1.0)*(ml*ml-1.0)
         ELSE
            mlp = 0.5*(ml+mla)
            wtl = 0.5*(1.0+ml/mla)
         END IF ! mla

         IF ( mra .le. 1.0 ) THEN 
            mrm = -0.25*(mr-1.0)*(mr-1.0)-0.125*(mr*mr-1.0)*(mr*mr-1.0)
            wtr = 0.25*(mr-1.0)*(mr-1.0)*(2.0+mr) -
     >            0.1875*mr*(mr*mr-1.0)*(mr*mr-1.0)
         ELSE
            mrm = 0.5*(mr-mra)
            wtr = 0.5*(1.0-mr/mra)
         END IF ! mla

         mf  = mlp + mrm
         mfa = ABS(mf)
         mfp = 0.5*(mf+mfa)
         mfm = 0.5*(mf-mfa)

         pf = wtl*pl(i) + wtr*pr(i)

! ******************************************************************************
! Compute fluxes
! ******************************************************************************

!        vf(1) = mfp*ul + mfm*ur ! I'm sure we'll need this someday
!        vf(2) = mfp*vl + mfm*vr
!        vf(3) = mfp*wl + mfm*wr

         flx(i,1)=(af*(mfp*rl(i)      +mfm*rr(i)   )        )*nm(i)
         flx(i,2)=(af*(mfp*rl(i)*ul(i)+mfm*rr(i)*ur(i))+pf*nx(i))*
     >            nm(i)
         flx(i,3)=(af*(mfp*rl(i)*vl(i)+mfm*rr(i)*vr(i))+pf*ny(i))*
     >            nm(i)
         flx(i,4)=(af*(mfp*rl(i)*wl(i)+mfm*rr(i)*wr(i))+pf*nz(i))*
     >            nm(i)
         flx(i,5)=(af*(mfp*rl(i)*Hl   +mfm*rr(i)*Hr) + pf*fs(i))*
     >            nm(i)
      enddo
C> @}
      return
      END

!-----------------------------------------------------------------------
! NOT LONG FOR THIS WORLD

      SUBROUTINE CentralInviscid_FluxFunction(ntot,nx,ny,nz,fs,ul,pl,
     >                                     ur,pr,flx)
! JH081915 More general, more obvious
! JH111815 HEY GENIUS THIS MAY BE SECOND ORDER AND THUS KILLING YOUR
!          CONVERGENCE. REPLACE WITH AUSM AND SHITCAN IT
! JH112015 This isn't why walls aren't converging. There's something
!          inherently second-order about your wall pressure. Think!
      real nx(ntot),ny(ntot),nz(ntot),fs(ntot),ul(ntot,5),pl(ntot),
     >     ur(ntot,5),pr(ntot) ! intent(in)
      real flx(ntot,5)! intent(out),dimension(5) ::

      do i=1,ntot
         rl =ul(i,1)
         rul=ul(i,2)
         rvl=ul(i,3)
         rwl=ul(i,4)
         rel=ul(i,5)

         rr =ur(i,1)
         rur=ur(i,2)
         rvr=ur(i,3)
         rwr=ur(i,4)
         rer=ur(i,5)

         ql = (rul*nx(i) + rvl*ny(i) + rwl*nz(i))/rl - fs(i)
         qr = (rur*nx(i) + rvr*ny(i) + rwr*nz(i))/rr - fs(i)

         flx(i,1) = 0.5*(ql* rl+ qr*rr               )
         flx(i,2) = 0.5*(ql* rul+pl(i)*nx(i) + qr* rur     +pr(i)*nx(i))
         flx(i,3) = 0.5*(ql* rvl+pl(i)*ny(i) + qr* rvr     +pr(i)*ny(i))
         flx(i,4) = 0.5*(ql* rwl+pl(i)*nz(i) + qr* rwr     +pr(i)*nz(i))
         flx(i,5) = 0.5*(ql*(rel+pl(i))+pl(i)*fs(i)+qr*(rer+pr(i))+
     >               pr(i)*fs(i))
      enddo

      return
      end

!-----------------------------------------------------------------------
      subroutine llf_euler(flx,ul,ur,wl,wr,nrm,dum)
! local Lax-Friedrichs done for the Euler equations (i.e., wave speed is
! hardcoded for gas dynamics) on an arbitrary pair of points sharing
! a normal vector nrm. the dummy argument dum is to match the signature
! of other two-point fluxes
      include 'SIZE' ! for ldim
      parameter (isnd=5) ! not sure yet
      real flx(5)
      real ul(5),ur(5),wl(5),wr(5),nrm(3),dum(3)
      real lambda
      rl=ul(1)
      rr=ur(1)
      al=wl(isnd)
      ar=wr(isnd)
      unrml=0.0
      unrmr=0.0
      do i=1,ldim
         unrml=unrml+nrm(i)*wl(i)
         unrmr=unrmr+nrm(i)*wr(i)
      enddo
      unrml=abs(unrml)+al
      unrmr=abs(unrmr)+ar
      lambda=0.5*max(unrml,unrmr)
      do i=1,toteq ! bake the minus sign into the formula here
!        flx(i)=-lambda*(ur(i)-ul(i))
         flx(i)=lambda*(ul(i)-ur(i))
      enddo
      return
      end

!-----------------------------------------------------------------------

      subroutine llf_euler_vec(wminus,uplus,flux,nstate) ! fstab
! local Lax-Friedrichs done for the Euler equations (i.e., wave speed is
! hardcoded for gas dynamics) on interior faces only (so gs_op is used
! to get neighbor values)
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'CMTDATA'
      include 'DG'

      real wminus(lx1*lz1*2*ldim*nelt,nstate)!*)
      real uplus(lx1*lz1*2*ldim*nelt,1) ! jumps are computed one at a time for now, so I only need
                                      ! U+ for a single conserved variable
      real flux(lx1*lz1*2*ldim*nelt,toteq) ! intent(inout); incremented

      parameter (lfq=lx1*lz1*2*ldim*lelt)
      common /SCRNS/ nx(lx1*lz1,2*ldim,lelt),ny(lx1*lz1,2*ldim,lelt),
     >               nz(lx1*lz1,2*ldim,lelt),jscr(lfq)
      real nx,ny,nz,jscr

      integer e,f,eq

      nfaces=2*ldim
      nxz=lx1*lz1
      nf=nxz*nfaces*nelt

      do e=1,nelt
         do f=1,nfaces
            call copy(nx(1,f,e),unx(1,1,f,e),nxz)
            call copy(ny(1,f,e),uny(1,1,f,e),nxz)
            if (if3d) call copy(nz(1,f,e),unz(1,1,f,e),nxz)
         enddo
      enddo

! overwrite wminus(:,jsnd) with max wave speed local to each GLL point.
! use uplus as scratch for awhile
      if (if3d) then
         call vdot3(uplus,wminus(1,jux),wminus(1,juy),wminus(1,juz),
     >              nx,ny,nz,nf)
      else
         call vdot2(uplus,wminus(1,jux),wminus(1,juy),nx,ny,nf)
      endif
      do i=1,nf
         wminus(i,jsnd)=0.5*(abs(uplus(i,1))+wminus(i,jsnd))
      enddo
! get max wave speed
      call fgslib_gs_op(dg_hndl,wminus(1,jsnd),1,4,0)

! do BC twice? not sure.
      call copy(jscr,jface,nf)
      call bcmask_cmt(jscr)

! now get flux from jump in the conserved variable one equation at a time
      do eq=1,toteq
         call faceu(eq,wminus(1,ju5)) ! peel U once or twice (thrice)?
         call face_state_commo(wminus(1,ju5),uplus,nf,1,dg_hndl)
         call sub2(uplus,wminus(1,ju5),nf)
! might math.f this someday
         do i=1,nf
            flux(i,eq)=flux(i,eq)-wminus(i,jsnd)*jscr(i)*uplus(i,1)
         enddo
      enddo

      return
      end

!-----------------------------------------------------------------------
! JH060618 Hooray for two-point fluxes! Needed for entropy-conserving DGSEM
!          in CMT-nek. These are volume functions here. calls to GS for
!          surface fluxes need dedicated subroutines.

      subroutine kennedygruber(flx,ul,ur,wl,wr,jal,jar)
      include 'SIZE' ! for ldim
      include 'CMTDATA' ! for i* indices
      real flx(5)
      real ul(5),ur(5),wl(4),wr(4),jal(3),jar(3)
      real rav,uav(3),pav,eav,jav(3)
      real rl,rr,pl,pr,qav,rq
      rl=ul(1)
      rr=ur(1)
      pl=wl(ipr)
      pr=wr(ipr)
      call rzero(jav,3)
      call rzero(uav,3)
      do j=1,ldim
         jav(j)=0.5*(jal(j)+jar(j))
         uav(j)=0.5*( wl(j)+ wr(j))
      enddo
      rav=0.5*(rl +rr )
      pav=0.5*(pl +pr )
!     eav=0.5*(ul(5)/rl+ur(5)/rr)
      eav=0.5*(ul(5)/ul(1)+ur(5)/ur(1)) !{{E}}
      qav=0.0
      do j=1,ldim
         qav=qav+uav(j)*jav(j) ! {{u}}.{{Ja}}
      enddo
      rq=rav*qav
      flx(1)=rq
      flx(2)=rq*uav(1)+jav(1)*pav
      flx(3)=rq*uav(2)+jav(2)*pav
      flx(4)=rq*uav(3)+jav(3)*pav
      flx(5)=rq*eav+pav*qav ! ({{rho}}{{E}}+{{p}}){{u}}
      return
      end

      subroutine kepec_ch(flx,ul,ur,wl,wr,jal,jar)
! kinetic energy preserving entropy-conservative flux a la Chandrashekar
! straight from his 2013 paper (well, GWK16 equation 3.19-20),
! no playing around with pressure
      include 'SIZE' ! for ldim
      include 'CMTDATA' ! for i* indices, and /CMTGASREF/
      external logmean
      real logmean
      real flx(5)
      real ul(5),ur(5),wl(4),wr(4),jal(3),jar(3)
      real rav,uav(3),pav,hhat,jav(3)
      real rl,rr,pl,pr,qav,rq
      rl=ul(1) !lol need phi
      rr=ur(1)
      pl=wl(ipr)
      pr=wr(ipr)
      bl=0.5*rl/pl
      br=0.5*rr/pr
      call rzero(jav,3)
      call rzero(uav,3)
      vsqav=0.0
      do j=1,ldim
         jav(j)=0.5*(jal(j)+jar(j))
         uav(j)=0.5*( wl(j)+ wr(j))
         vsqav=vsqav+0.5*(wl(j)**2+wr(j)**2)
      enddo
      rav=0.5*(rl +rr ) !temporary
      bav=0.5*(bl+br) !temporary
      phat=0.5*rav/bav
      rav=logmean(rl,rr)
      bav=logmean(bl,br)
      qav=0.0
      do j=1,ldim
         qav=qav+uav(j)*jav(j) ! {{u}}.{{Ja}}
      enddo
      rq=rav*qav
      flx(1)=rq
      flx(2)=rq*uav(1)+jav(1)*phat
      flx(3)=rq*uav(2)+jav(2)*phat
      flx(4)=rq*uav(3)+jav(3)*phat
      flx(5)=flx(1)*(0.5/(bav*(gmaref-1.0))-0.5*vsqav)+uav(1)*flx(2)+
     >       uav(2)*flx(3)+uav(3)*flx(4)
      return
      end

      subroutine kennedygruber_vec(z,flux,nstate,nflux) ! fsharp
! JH111218 Kennedy-Gruber fluxes, but acting on vectors of face nodes
!  instead of two arbitrary points. Gratuitously assuming watertight geometry.
! JH122518 Ask Paul if multiplication of elements of z is faster via gs_op
!          given that z+=z-={{stuff}}
      include 'SIZE'
      include 'INPUT' ! for if3d
      include 'GEOM' ! for normal vectors at faces
      include 'CMTDATA' ! for jface

! ==============================================================================
! Arguments
! ==============================================================================
      integer nstate,nflux
      real z(lx1*lz1*2*ldim*nelt,nstate),
     >     flux(lx1*lz1*2*ldim*nelt,nflux)
      
      parameter (lfq=lx1*lz1*2*ldim*lelt)
      common /SCRNS/ scrf(lfq),scrg(lfq),scrh(lfq),fdot(lfq),jscr(lfq),
     >                 nx(lx1*lz1,2*ldim,lelt),ny(lx1*lz1,2*ldim,lelt),
     >                 nz(lx1*lz1,2*ldim,lelt)
      real scrf,scrg,scrh,fdot,jscr,nx,ny,nz
      integer e,f

!-----------------------------------------------------------------------
! I need to get this stuff out of this routine and one level up
!-----------------------------------------------------------------------
      nfaces=2*ldim
      nxz=lx1*lz1
      nf=nxz*nfaces*nelt

! I don't know what to do with volume fraction phi, and this is my first guess
      call col3(jscr,jface,z(1,jph),nf) ! Jscr=JA*{{\phi_g}}

! zero out jscr at boundary faces; gs_op is degenerate there.
      call bcmask_cmt(jscr)

      do e=1,nelt
         do f=1,nfaces
            call copy(nx(1,f,e),unx(1,1,f,e),nxz)
            call copy(ny(1,f,e),uny(1,1,f,e),nxz)
            if (if3d) call copy(nz(1,f,e),unz(1,1,f,e),nxz)
         enddo
      enddo
!-----------------------------------------------------------------------
! I need to get that stuff out of this routine and one level up
!-----------------------------------------------------------------------

! mass. scrF={{rho}}{{u}}, scrG={{rho}}{{v}}, scrH={{rho}}{{w}}
      call col3(scrf,z(1,jrhof),z(1,jux),nf)
      call col3(scrg,z(1,jrhof),z(1,juy),nf)
      if (if3d) then
         call col3(scrh,z(1,jrhof),z(1,juz),nf)
         call vdot3(fdot,scrf,scrg,scrh,nx,ny,nz,nf)
      else
         call vdot2(fdot,scrf,scrg,nx,ny,nf)
      endif
      call add2col2(flux(1,1),fdot,jscr,nf)

! x-momentum
      call col3(fdot,scrf,z(1,jux),nf) ! F={{rho}}{{u}}{{u}}
      call add2(fdot,z(1,jpr),nf) ! F+={{p}}
      call col2(fdot,nx,nf) ! F contribution to f~
      call addcol4(fdot,scrf,z(1,juy),ny,nf) ! G={{rho}}{{v}}{{u}} .ny -> f~
      if (if3d) call addcol4(fdot,scrf,z(1,juz),nz,nf) ! H={{rho}}{{w}}{{u}} .nz -> f~
      call add2col2(flux(1,2),fdot,jscr,nf)

! y-momentum
      call col3(fdot,scrg,z(1,juy),nf) ! G={{rho}}{{v}}{{v}}
      call add2(fdot,z(1,jpr),nf) ! G+={{p}}
      call col2(fdot,ny,nf)
      call addcol4(fdot,scrg,z(1,jux),nx,nf)

      if (if3d) then
         call addcol4(fdot,scrg,z(1,juz),nz,nf)
         call add2col2(flux(1,3),fdot,jscr,nf)
! z-momentum
         call col3(fdot,scrh,z(1,juz),nf)
         call add2(fdot,z(1,jpr),nf)
         call col2(fdot,nz,nf)
         call addcol4(fdot,scrh,z(1,jux),nx,nf)
         call addcol4(fdot,scrh,z(1,juy),ny,nf)
         call add2col2(flux(1,4),fdot,jscr,nf)
      else ! 2D only. couldn't resist deleting one if(if3d)
         call add2col2(flux(1,3),fdot,jscr,nf)
      endif

! energy ({{rho}}{{E}}+{{p}}){{u}}.n
      call col2(scrf,z(1,ju5),nf)
      call col2(scrg,z(1,ju5),nf)
      call add2col2(scrf,z(1,jux),z(1,jpr),nf)
      call add2col2(scrg,z(1,juy),z(1,jpr),nf)
      if (if3d) then
         call col2(scrh,z(1,ju5),nf)
         call add2(scrh,z(1,juz),z(1,jpr),nf)
         call vdot3(fdot,scrf,scrg,scrh,nx,ny,nz,nf)
      else
         call vdot2(fdot,scrf,scrg,nx,ny,nf)
      endif
      call add2col2(flux(1,5),fdot,jscr,nf)

      return
      end

!-----------------------------------------------------------------------
! Parameter vectors. need to work out some way of specifying indices
! instead of hardcoding them
      subroutine trivial
      return
      end

      subroutine rhoe_to_e(fatface,nf,ns)
! Routine to convert primitive variables on face to parameter vector for
! entropy-stable numerical fluxes consistent with volume fluxes.
! converts U5=phi*rho*E on faces to E=e+1/2*ui*ui for Kennedy-Gruber fluxes
! consider putting indices for quantities (ju5, jrho, etc.) in the argument
! list instead of including CMTDATA
      include 'SIZE'
      include 'CMTDATA'
      real fatface(nf,ns)
      call invcol2(fatface(1,ju5),fatface(1,jrhof),nf)
      call invcol2(fatface(1,ju5),fatface(1,jph),nf)
      return
      end

!-----------------------------------------------------------------------
! JUNKYARD
!-----------------------------------------------------------------------

      subroutine KEPEC_duplicated(wminus,wplus,flux)
! Chandrashekar KEPEC flux done the expensive, naive way.
      include 'SIZE'
      include 'INPUT' ! do we need this?
      include 'GEOM' ! for unx
      include 'CMTDATA' ! do we need this without outflsub?
      include 'TSTEP' ! for ifield?
      include 'DG'
      external logmean

! ==============================================================================
! Arguments
! ==============================================================================
      integer nstate,nflux
      real wminus(lx1*lz1,2*ldim,nelt,nqq),
     >     wplus(lx1*lz1,2*ldim,nelt,nqq),
     >     flux(lx1*lz1,2*ldim,nelt,toteq)

! ==============================================================================
! Locals
! ==============================================================================

      integer e,f,fdim,i,k,nxz,nface
      parameter (lfd=lx1*lz1)
! nx,ny,nz : outward facing unit normal components
! jaco_c   : fdim-D GLL grid Jacobian
!
      COMMON /SCRNS/ jaco_c(lx1*lz1)
      real jaco_c
      real ul(5),ur(5),wl(4),wr(4),jal(3),jar(3),flx(5)

      nface = 2*ldim
      nxz   = lx1*lz1

      do e=1,nelt
      do f=1,nface

         call rzero(jal,3)
         do i=1,nxz
            jal(1)=unx(i,1,f,e)
            jal(2)=uny(i,1,f,e)
            if (if3d) jal(3)=unz(i,1,f,e)
            ul(1)=wminus(i,f,e,ju1)
            ul(2)=wminus(i,f,e,ju2)
            ul(3)=wminus(i,f,e,ju3)
            ul(4)=wplus(i,f,e,ju4)
            ur(1)=wplus(i,f,e,ju1)
            ur(2)=wplus(i,f,e,ju2)
            ur(3)=wplus(i,f,e,ju3)
            ur(4)=wplus(i,f,e,ju4)
            wl(1)=wminus(i,f,e,jux)
            wl(2)=wminus(i,f,e,jux)
            wl(3)=wminus(i,f,e,jux)
            wl(ipr)=wminus(i,f,e,jpr)
            wr(1)=wplus(i,f,e,jux)
            wr(2)=wplus(i,f,e,jux)
            wr(3)=wplus(i,f,e,jux)
            wr(ipr)=wplus(i,f,e,jpr)
            call copy(jar,jal,3)
            call kepec_ch(flx,ul,ur,wl,wr,jal,jar)
! max wave speed
            al=
