C> @file outflow_bc.f Dirichlet states for outflow boundary conditions
C> wrapper for other BC routines. Just one for now. More to come.
      subroutine outflow(f,e,wminus,wplus,uminus,uplus,nvar)
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'CMTBCDATA'
      integer nvar,f,e
      real wminus(nvar,lx1*lz1),wplus(nvar,lx1*lz1),
     >     uminus(toteq,lx1*lz1),uplus(toteq,lx1*lz1)

      call outflow_df(f,e,wminus,wplus,uminus,uplus,nvar)
! diagnostic
      do i=1,lxz
         write(57,'(a6,2i4,10e15.7)')'f,e,w-',f,e,(wminus(j,i),j=1,nvar)
         write(58,'(a6,2i4,10e15.7)')'f,e,w+',f,e,(wplus(j,i),j=1,nvar)
         write(59,'(a6,2i4,5e15.7)')'f,e,u-',f,e,(uminus(j,i),j=1,toteq)
         write(60,'(a6,2i4,5e15.7)')'f,e,u+',f,e,(uplus(j,i),j=1,toteq)
      enddo
! diagnostic

      return
      end

!--------------------------------------------------------------------

C> \ingroup isurf
C> @{
      subroutine outflow_df(f,e,wm,wp,um,up,nvar)
C> more conventional Dolejsi & Feistauer (2015) Section 8.3.2.2
C> ``physical'' boundary conditions. Also encountered in
C> Hartmann & Houston (2006). A poor default.
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'CMTDATA'

      integer f,e,nvar ! intent(in)
      real wm(nvar,lx1*lz1),wp(nvar,lx1*lz1),
     >     um(toteq,lx1*lz1),up(toteq,lx1*lz1)
      real mach
      integer eq

      nxz=lx1*lz1
      ieg=lglel(e)

      call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,f)    
      l=0
      do iz=k0,k1
      do iy=j0,j1
      do ix=i0,i1
         call nekasgn(ix,iy,iz,e)
         call cmtasgn(ix,iy,iz,e)
         call userbc (ix,iy,iz,f,ieg) ! get pinfty and e_internal=\rho e
         l=l+1
         do eq=1,toteq
            um(eq,l)=u(ix,iy,iz,eq,e)
         enddo
         wm(iux,l)=vx(ix,iy,iz,e)
         wm(iuy,l)=vy(ix,iy,iz,e)
         wm(iuz,l)=vz(ix,iy,iz,e)
         wm(ipr,l)=pr(ix,iy,iz,e)
         wm(ithm,l)=t(ix,iy,iz,e,1)
         wm(irho,l)=vtrans(ix,iy,iz,e,jrho)
         wm(isnd,l)=csound(ix,iy,iz,e)
         wm(iph,l)=phig(ix,iy,iz,e)

         wp(iux,l)= wm(iux,l)
         wp(iuy,l)= wm(iuy,l)
         wp(iuz,l)= wm(iuz,l)
         wp(iph,l)  = wm(iph,l)
         wp(irho,l)= wm(irho,l)
         call copy(up(1,l),um(1,l),4)

         snx  = unx(l,1,f,e)
         sny  = uny(l,1,f,e)
         snz  = unz(l,1,f,e)
         mach = abs(wm(iux,l)*snx+wm(iuy,l)*sny+wm(iuz,l)*snz)/
     >          wm(isnd,l)

         if (mach.lt.1.0) then ! userbc should have set this to pinfty

            wp(ipr,l)  = pres ! userbc should have set this to pinfty
            wp(isnd,l) = asnd ! userbc should have set this to a(pinfty,rho-)
            wp(ithm,l) = temp   ! userbc should have set this to T(pinfty,rho-)
!           up(5,l)=wm(jrho,l)*e_internal ! userbc plz set e_internal(temp)
            up(5,l)=e_internal ! here AND ONLY HERE is e_internal density-weighted
     >          +0.5*wm(irho,l)*(wm(iux,l)**2+wm(iuy,l)**2+wm(iux,l)**2)
            up(5,l)=up(5,l)*wm(iph,l)

         else ! supersonic outflow

            wp(ipr,l)  = wm(ipr,l)
            wp(isnd,l) = wm(isnd,l)
            wp(ithm,l) = wm(ithm,l)
            up(5,l)=um(5,l)

         endif

         call copy(up(1,l),um(1,l),4)

      enddo
      enddo
      enddo

      return
      end

!--------------------------------------------------------------------

      subroutine outflow_rflu(nvar,f,e,facew,wbc)
      include 'SIZE'
      include 'NEKUSE'
      include 'CMTDATA'
      include 'CMTBCDATA'
      include 'INPUT'
      include 'GEOM'
      include 'PARALLEL'
      include 'DG'

      integer i,bcOpt
      integer  f,e,fdim
      real facew(lx1*lz1,2*ldim,nelt,nvar)
      real wbc(lx1*lz1,2*ldim,nelt,nvar)
      real sxn,syn,szn,rhou,rhov,rhow,pl,rhob,rhoub,rhovb,rhowb,rhoeb

      nxz=lx1*lz1
      nxzd=lxd*lzd
      fdim=ldim-1
      ieg=lglel(e)

      call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,f)    
      l=0
      do iz=k0,k1
      do iy=j0,j1
      do ix=i0,i1
         call nekasgn(ix,iy,iz,e)     ! gives us phi- and rho-
         call cmtasgn(ix,iy,iz,e)
         call userbc (ix,iy,iz,f,ieg) ! just for molarmass, and
                                      ! pres
         l=l+1
         sxn = unx(l,1,f,e)
         syn = uny(l,1,f,e)
         szn = unz(l,1,f,e)
         rhou= facew(l,f,e,iu2)/phi
         rhov= facew(l,f,e,iu3)/phi
         rhow= facew(l,f,e,iu4)/phi
         rhoe= facew(l,f,e,iu5)/phi
         pl= facew(l,f,e,ipr) ! P- here
         wbc(l,f,e,jcpf)=facew(l,f,e,jcpf)
         wbc(l,f,e,jcvf)=facew(l,f,e,jcvf)
         cp=facew(l,f,e,jcpf)/rho
         cv=facew(l,f,e,jcvf)/rho
c        fs = 0.0
         if(outflsub)then
            pres= pinfty
            idbc=1
         else
            pres= facew(l,f,e,ipr)
            idbc=0
         endif
         call BcondOutflowPerf(idbc,pres,sxn,syn,szn,cp,molarmass,
     >                         rho,rhou,rhov,rhow,rhoe,pl,
     >                         rhob,rhoub,rhovb,rhowb,rhoeb )
         wbc(l,f,e,jrho)=rhob
         wbc(l,f,e,jux)=rhoub/rhob
         wbc(l,f,e,juy)=rhovb/rhob
         wbc(l,f,e,juz)=rhowb/rhob
! dammit fix this. tdstate to the rescue?
         wbc(l,f,e,jthm)=(rhoeb-0.5*(rhoub**2+rhovb**2+rhowb**2)/rhob)/
     >                   cv
! dammit fix that
         wbc(l,f,e,ju1)=rhob*phi
         wbc(l,f,e,ju2)=rhoub*phi
         wbc(l,f,e,ju3)=rhovb*phi
         wbc(l,f,e,ju4)=rhowb*phi
         wbc(l,f,e,ju5)=rhoeb*phi
         wbc(l,f,e,jph)=phi
         wbc(l,f,e,jpr)=pres
! dammit fix this. tdstate to the rescue?
         wbc(l,f,e,jsnd)=sqrt(cp/cv*pres/rho)
! dammit fix that
      enddo
      enddo
      enddo

      return
      end

!******************************************************************************
!
! Purpose: set outflow boundary condition for one cell.
!
! Description: the subsonic boundary condition is based on non-reflecting,
!              characteristics method of Whitfield and Janus: Three-Dimensional
!              Unsteady Euler Equations Solution Using Flux Vector Splitting.
!              AIAA Paper 84-1552, 1984. The supersonic boundary condition
!              consists of simple extrapolation.
!
! Input: bcOpt    = boundary treatment: subsonic, supersonic, or mixed
!        pout     = given static outlet pressure
!        sx/y/zn  = components of ortho-normalized face vector (outward facing)
!        cpgas    = specific heat at constant pressure (boundary cell)
!        mol      = molecular mass at boundary cell
!        rho      = density at boundary cell
!        rhou/v/w = density * velocity components at boundary cell
!        rhoe     = density * total energy at boundary cell
!        press    = static pressure at boundary cell
!
! Output: rhob      = density at boundary
!         rhou/v/wb = density * velocity components at boundary
!         rhoeb     = density * total energy at boundary
!
! Notes: this condition is valid only for thermally and calorically
!        perfect gas (supersonic outflow valid for all gases).
!
!******************************************************************************
!
! $Id: BcondOutflowPerf.F90,v 1.1.1.1 2014/05/05 21:47:47 tmish Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

      SUBROUTINE BcondOutflowPerf( bcOpt,pout,sxn,syn,szn,cpgas,mol,
     >                       rho,rhou,rhov,rhow,rhoe,press,
     >                       rhob,rhoub,rhovb,rhowb,rhoeb )

      IMPLICIT NONE
      integer bcopt_subsonic,bcopt_mixed
      parameter (BCOPT_SUBSONIC=1)
      parameter (BCOPT_MIXED=0)
      real MixtPerf_C_DGP,MixtPerf_Eo_DGPUVW,MixtPerf_G_CpR,MixtPerf_R_M
     >   , MixtPerf_P_DEoGVm2
      external MixtPerf_C_DGP,MixtPerf_Eo_DGPUVW,MixtPerf_G_CpR,
     >     MixtPerf_R_M,MixtPerf_P_DEoGVm2

! ... parameters
      INTEGER bcOpt

      REAL pout
      REAL rho, rhou, rhov, rhow, rhoe, press
      REAL sxn, syn, szn, cpgas, mol
      REAL rhob, rhoub, rhovb, rhowb, rhoeb

! ... local variables
      REAL csound, rgas, gamma, gam1, u, v, w, mach, rrhoc, deltp,
     >            ub, vb, wb, vnd

!******************************************************************************
! gas properties; velocity components; Mach number

      rgas  = MixtPerf_R_M( mol )
      gamma = MixtPerf_G_CpR( cpgas,rgas )
      gam1  = gamma - 1.0

      u      = rhou/rho
      v      = rhov/rho
      w      = rhow/rho
      csound = MixtPerf_C_DGP( rho,gamma,press )
      mach   = SQRT(u*u+v*v+w*w)/csound

! subsonic flow ---------------------------------------------------------------

      IF (mach .lt. 1.0 .AND.
     >(bcOpt .eq. BCOPT_SUBSONIC .OR. bcOpt .eq. BCOPT_MIXED)) THEN
         rrhoc = 1.0/(rho*csound)
         deltp = press - pout
         rhob  = rho - deltp/(csound*csound)
         ub    = u + sxn*deltp*rrhoc
         vb    = v + syn*deltp*rrhoc
         wb    = w + szn*deltp*rrhoc

! - special treatment to prevent "deltp" from changing the sign
!   of velocity components. This may happen for very small u, v, w.

         vnd = ub*sxn + vb*syn + wb*szn
         IF ( vnd .lt. 0.0 ) THEN ! inflow at outflow boundary
            ub = SIGN(1.0,u)*MAX(ABS(ub),ABS(u))
            vb = SIGN(1.0,v)*MAX(ABS(vb),ABS(v))
            wb = SIGN(1.0,w)*MAX(ABS(wb),ABS(w))
         END IF ! vnd

         rhoub = rhob*ub
         rhovb = rhob*vb
         rhowb = rhob*wb
         rhoeb = rhob*MixtPerf_Eo_DGPUVW( rhob,gamma,pout,ub,vb,wb )

! supersonic flow -------------------------------------------------------------

      ELSE
         rhob  = rho
         rhoub = rhou
         rhovb = rhov
         rhowb = rhow
         rhoeb = rhoe
      END IF ! mach

      end

!******************************************************************************
!
! RCS Revision history:
!
! $Log: BcondOutflowPerf.F90,v $
! Revision 1.1.1.1  2014/05/05 21:47:47  tmish
! Initial checkin for rocflu macro.
!
! Revision 1.3  2008/12/06 08:43:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:24  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/03/26 20:21:09  haselbac
! Fix mistake in declarations
!
! Revision 1.1  2004/12/01 16:48:04  haselbac
! Initial revision after changing case
!
! Revision 1.5  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.2  2002/06/22 00:49:50  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.1  2002/06/10 21:19:34  haselbac
! Initial revision
!
!******************************************************************************
