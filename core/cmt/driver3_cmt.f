! JH100919
! Now with some code introduced by Jacob's work on Tait-equation water
! mixture modeling. Commented out until shocks are working in essplit
C> @file driver3_cmt.f routines for primitive variables, usr-file interfaces
C> and properties

C> Compute primitive variables (velocity, thermodynamic state) from 
C> conserved unknowns U
      subroutine compute_primitive_vars(ilim)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'GEOM'
      include 'CMTDATA'
      include 'SOLN'
      include 'DEALIAS' ! until we are comfortable with setup_convect

      parameter (lxyz=lx1*ly1*lz1)
      common /ctmp1/ energy(lx1,ly1,lz1),scr(lx1,ly1,lz1)
      integer e, eq
      common /posflags/ ifailr,ifaile,ifailt,ilimflag
      integer ifailr,ifaile,ifailt,ilimflag
      integer ilim

      nxyz= lx1*ly1*lz1
      ntot=nxyz*nelt
      ifailr=-1
      ifaile=-1
      ifailt=-1
      ilimflag=ilim

      do e=1,nelt
! JH020918 long-overdue sanity checks
         dmin=vlmin(u(1,1,1,irg,e),nxyz)
         if (dmin .lt. 0.0 .and. ilim .ne. 0) then
            ifailr=lglel(e)
            write(6,*) nid,'***NEGATIVE DENSITY***',dmin,lglel(e)
         endif
         call invcol3(vx(1,1,1,e),u(1,1,1,irpu,e),u(1,1,1,irg,e),nxyz)
         call invcol3(vy(1,1,1,e),u(1,1,1,irpv,e),u(1,1,1,irg,e),nxyz)
!        if (if3d)
         call invcol3(vz(1,1,1,e),u(1,1,1,irpw,e),u(1,1,1,irg,e),nxyz)
! first kinetic energy
         if (if3d) then
            call vdot3(scr,
     >             u(1,1,1,irpu,e),u(1,1,1,irpv,e),u(1,1,1,irpw,e),
     >             u(1,1,1,irpu,e),u(1,1,1,irpv,e),u(1,1,1,irpw,e),nxyz)
         else
            call vdot2(scr,u(1,1,1,irpu,e),u(1,1,1,irpv,e),
     >                     u(1,1,1,irpu,e),u(1,1,1,irpv,e),nxyz)
         endif
         call invcol2(scr,u(1,1,1,irg,e),nxyz)
         call cmult(scr,0.5,nxyz)
! then to internal energy
         call sub3(energy,u(1,1,1,iret,e),scr,nxyz)
! now mass-specific
         call invcol2(energy,u(1,1,1,irg,e),nxyz)
! don't forget to get density where it belongs
         call invcol3(vtrans(1,1,1,e,jrho),u(1,1,1,irg,e),phig(1,1,1,e),
     >                nxyz)
! JH020718 long-overdue sanity checks
         emin=vlmin(energy,nxyz)
         if (emin .lt. 0.0 .and. ilim .ne. 0) then
            ifaile=lglel(e)
            write(6,*) stage,nid, ' HAS NEGATIVE ENERGY ',emin,lglel(e)
         endif
!! JH070219 Tait mixture model mass fractions. just one for now
!c JB080119 go throug hmultiple species
!         call invcol3(t(1,1,1,e,2),u(1,1,1,imfrac,e),
!     >                u(1,1,1,irg,e),
!     >                nxyz)
!c        do iscal = 1,NPSCAL
!c        call invcol3(t(1,1,1,e,1+iscal),u(1,1,1,imfrac+iscal-1,e),
!c    >                u(1,1,1,irg,e),
!c    >                nxyz)
!c        enddo
         call tdstate(e,energy) ! compute state, fill ifailt
      enddo

! Avoid during EBDG testing
! JH070219 Tait mixture model: man up and test T(:,2) for positivity
!          someday.
      call poscheck(ifailr,'density    ')
      call poscheck(ifaile,'energy     ')
      call poscheck(ifailt,'temperature')

      return
      end

!-----------------------------------------------------------------------

C> Compute thermodynamic state for element e from internal energy.
C> usr file.
      subroutine tdstate(e,energy)!,energy)
c compute the gas properties. We will have option to add real gas models
c We have perfect gas law. Cvg is stored full field
      include 'SIZE'
      include 'CMTDATA'
      include 'SOLN'
      include 'PARALLEL'
      include 'NEKUSE'
      integer   e,eg
      real energy(lx1,ly1,lz1)

      common /posflags/ ifailr,ifaile,ifailt,ilimflag
      integer ifailr,ifaile,ifailt,ilimflag

      eg = lglel(e)
      do k=1,lz1
      do j=1,ly1
      do i=1,lx1
         call nekasgn(i,j,k,e)
         call cmtasgn(i,j,k,e)
         e_internal=energy(i,j,k) !cmtasgn should do this, but can't
         call cmt_userEOS(i,j,k,eg)
! JH020718 long-overdue sanity checks
         if (temp .lt. 0.0 .and. ilimflag .ne. 0) then
            ifailt=eg
            write(6,'(i6,a26,e12.4,3i2,i8,3e15.6)') ! might want to be less verbose
     >      nid,' HAS NEGATIVE TEMPERATURE ', x,i,j,k,eg,temp,rho,pres
         endif
         vtrans(i,j,k,e,jen)= e_internal
         vtrans(i,j,k,e,jcv)= cv*rho
         vtrans(i,j,k,e,jcp)= cp*rho
         t(i,j,k,e,1)       = temp
         pr(i,j,k,e)        = pres
         csound(i,j,k,e)    = asnd
      enddo
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine cmtasgn (ix,iy,iz,e)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA'
      include 'NEKUSE'

      integer e,eqnum
!     do eqnum=1,toteq
!        varsic(eqnum)=u(ix,iy,iz,eqnum,e)  
!     enddo
      phi  = phig  (ix,iy,iz,e)
      rho  = vtrans(ix,iy,iz,e,jrho)
      pres = pr    (ix,iy,iz,e)
      if (rho.ne.0) then
         cv   = vtrans(ix,iy,iz,e,jcv)/rho
         cp   = vtrans(ix,iy,iz,e,jcp)/rho
         e_internal = vtrans(ix,iy,iz,e,jen)
      endif
      asnd = csound(ix,iy,iz,e)
      mu     = vdiff(ix,iy,iz,e,jmu)
      udiff  = vdiff(ix,iy,iz,e,jknd)
! MAKE SURE WE''RE NOT USING UTRANS FOR ANYTHING IN pre-v16 code!!
      lambda = vdiff(ix,iy,iz,e,jlam)

      return
      end

!-----------------------------------------------------------------------

      subroutine cmt_ics
! overlaps with setics. -DCMT will require IFDG as well
      include 'SIZE'
      include 'SOLN'
      include 'GEOM'
      include 'PARALLEL'
      include 'CMTDATA'
      include 'NEKUSE'
      nxyz2=lx2*ly2*lz2       ! Initialize all fields:
      ntot2=nxyz2*nelv
      nxyz1=lx1*ly1*lz1
      ntott=nelt*nxyz1
      ntotv=nelv*nxyz1
      ltott=lelt*nxyz1
      ntotcv=lelt*nxyz1*toteq
      call rone(phig,ltott)
      call rzero(csound,ltott)
      call rzero(vtrans,ltott*ldimt1)
      call rzero(vdiff ,ltott*ldimt1)
      call rzero(u,ntotcv)
! JH100919 where does particle stuff live these days?
!     call usr_particles_init
      call cmtuic
      if(ifrestart) call my_full_restart !  Check restart files. soon...

C print min values
      xxmax = glmin(xm1,ntott)
      yymax = glmin(ym1,ntott)
      zzmax = glmin(zm1,ntott)

      vxmax = glmin(vx,ntotv)
      vymax = glmin(vy,ntotv)
      vzmax = glmin(vz,ntotv)
      prmax = glmin(pr,ntot2)

      ntot = nxyz1*nelt
      ttmax = glmin(t ,ntott)

      if (nio.eq.0) then
         write(6,19) xxmax,yymax,zzmax
   19    format('Cxyz min  ',5g25.18)
      endif
      if (nio.eq.0) then
         write(6,20) vxmax,vymax,vzmax,prmax,ttmax
   20    format('Cuvwpt min',5g25.18)
      endif

c print max values
      xxmax = glmax(xm1,ntott)
      yymax = glmax(ym1,ntott)
      zzmax = glmax(zm1,ntott)

      vxmax = glmax(vx,ntotv)
      vymax = glmax(vy,ntotv)
      vzmax = glmax(vz,ntotv)
      prmax = glmax(pr,ntot2)

      ntot = nxyz1*nelt
      ttmax = glmax(t ,ntott)

      if (nio.eq.0) then
         write(6,16) xxmax,yymax,zzmax
   16    format('Cxyz max  ',5g25.18)
      endif

      if (nio.eq.0) then
         write(6,17) vxmax,vymax,vzmax,prmax,ttmax
   17    format('Cuvwpt max',5g25.18)
      endif

c     ! save velocity on fine mesh for dealiasing
!     call setup_convect(2) ! check what this does again. might be a good
!                           ! idea, or it might be counterproductive
      if(nio.eq.0) then
        write(6,*) 'done :: set initial conditions, CMT-nek'
        write(6,*) ' '
      endif
      return
      end

!-----------------------------------------------------------------------

      subroutine cmtuic
! overlaps with setics. -DCMT will require IFDG as well
! need to make sure setics has no effect.
! JH070219 cmtuic now sets U and U alone. EVERYTHING else should come
!          from compute_primitive_variables
      include 'SIZE'
      include 'SOLN'
      include 'PARALLEL'
      include 'CMTDATA'
      include 'NEKUSE'
      integer e,eg
      do e=1,nelt
         eg = lglel(e)
         do k=1,lz1
         do j=1,ly1
         do i=1,lx1           
            call nekasgn (i,j,k,e)
            call cmtasgn (i,j,k,e)
            call useric  (i,j,k,eg)
            phig(i,j,k,e)  = phi ! only sane way to run CMT-nek without
                                 ! particles is to have useric set phi=1
            u(i,j,k,irg,e) = phi*rho
            u(i,j,k,irpu,e)= phi*rho*ux
            u(i,j,k,irpv,e)= phi*rho*uy
            u(i,j,k,irpw,e)= phi*rho*uz
            u(i,j,k,iret,e)=phi*rho*(e_internal+0.5*(ux**2+uy**2+uz**2))
!            u(i,j,k,imfrac,e)=phi*rho*ps(1)
!c JB080119 multiple species
!               t(i,j,k,e,2) = ps(l)
!c           do l = 2,NPSCAL
!c               t(i,j,k,e,l) = ps(l-1)
!c           enddo
         enddo
         enddo
         enddo
      enddo
      return
      end

!-----------------------------------------------------------------------

      subroutine poscheck(ifail,what)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA'
      include 'PARALLEL'
      include 'INPUT'
!JH020918 handles reporting, I/O and exit from failed positivity checks
!         in compute_primitive_variables
      character*11 what

      ifail0=iglmax(ifail,1)
      if(ifail0 .ne. -1) then
         if (nio .eq. 0)
     >   write(6,*) 'dumping solution after negative ',what,'@ eg=',
     >             ifail0
         ifxyo=.true.
!        call out_fld_nek
         call outpost2(vx,vy,vz,pr,t,ldimt,'EBL')
         call exitt
      endif

      return
      end
