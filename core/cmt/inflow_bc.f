C> @file Dirichlet states for inflow boundary conditions
C> wrapper for other BC routines. Just one for now. More to come.
      subroutine inflow(f,e,wminus,wplus,uminus,uplus,nvar)
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'CMTBCDATA'
      integer nvar,f,e
      real wminus(nvar,lx1*lz1),wplus(nvar,lx1*lz1),
     >     uminus(toteq,lx1*lz1),uplus(toteq,lx1*lz1)

      call inflow_df(f,e,wminus,wplus,uminus,uplus,nvar)
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

      subroutine inflow_rflu(nvar,f,e,wminus,wbc)
! JH112718 Poorly thought-out routine to get rind state from RocFlu's
!          implementation of Holmes' inflow conditions. Backburnered
!          for awhile
      include 'SIZE'
      include 'INPUT'
      include 'NEKUSE'
      include 'CMTDATA'
      include 'GEOM'
      include 'PARALLEL'
      include 'DG'
      include 'PERFECTGAS'

      integer f,e,fdim ! intent(in)
      integer i,bcOptType
      real wminus(lx1*lz1,2*ldim,nelt,nvar) ! intent(in)
      real wbc(lx1*lz1,nvar) ! intent(in)
      real snx,sny,snz,rhou,rhov,rhow,pl,rhob,rhoub,rhovb,rhowb
     >     ,rhoeb, mach

      nxz=lx1*lz1
      nxzd=lxd*lzd
      fdim=ldim-1
      ieg=lglel(e)

      call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,f)    
      l=0
      do iz=k0,k1
      do iy=j0,j1
      do ix=i0,i1
         call nekasgn(ix,iy,iz,e)
         call cmtasgn(ix,iy,iz,e)
         call userbc (ix,iy,iz,f,ieg) ! get molarmass asnd phi t0in p0in cp cv
c                                     !     ux,uy,uz
         l=l+1

         bcOptType=0
         snx  = unx(l,1,f,e)
         sny  = uny(l,1,f,e)

!        rho  = facew(l,f,e,ju1)/facew(l,f,e,jph)
!        rhou = facew(l,f,e,ju2)/facew(l,f,e,jph)
!        rhov = facew(l,f,e,ju3)/facew(l,f,e,jph)
!        rhow = facew(l,f,e,ju4)/facew(l,f,e,jph)
!        rhoe = facew(l,f,e,ju5)/facew(l,f,e,jph)
         rho  = u(ix,iy,iz,1,e)/phig(ix,iy,iz,e)
         rhou = u(ix,iy,iz,2,e)/phig(ix,iy,iz,e)
         rhov = u(ix,iy,iz,3,e)/phig(ix,iy,iz,e)
         rhow = u(ix,iy,iz,4,e)/phig(ix,iy,iz,e)
         rhoe = u(ix,iy,iz,5,e)/phig(ix,iy,iz,e)

         if (if3d) then
            mach = sqrt(ux**2+uy**2+uz**2)/asnd
            snz  = unz(l,1,f,e)
         else
            mach = sqrt(ux**2+uy**2)/asnd
            snz=0.0
         endif
         if (mach.lt.1.0) bcOptType=1

         call BcondInflowPerf(bcOptType,0,p0in,t0in,ux,uy,uz
     >                       ,mach,snx,sny,snz,cp
     >                       ,molarmass,rho,rhou,rhov,rhow,rhob,rhoub
     >                       ,rhovb,rhowb,rhoeb,pres,asnd,temp)
         
         wbc(l,jrho) = rhob
         wbc(l,jux)  = ux
         wbc(l,juy)  = uy
         wbc(l,juz)  = uz
         wbc(l,jsnd) = asnd ! overwritten by Bcond
         wbc(l,jpr)  = pres ! overwritten by Bcond
         wbc(l,jthm) = temp ! overwritten by Bcond
         wbc(l,jcpf) = rho*cp
         wbc(l,jcvf) = rho*cv
         wbc(l,ju1)  = rhob*phi
         wbc(l,ju2)  = rhoub*phi
         wbc(l,ju3)  = rhovb*phi
         wbc(l,ju4)  = rhowb*phi
         wbc(l,ju5)  = rhoeb*phi
      enddo
      enddo
      enddo

      return
      end

!--------------------------------------------------------------------

C> \ingroup isurf
C> @{
      subroutine inflow_df(f,e,wm,wp,um,up,nvar)
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
         call userbc (ix,iy,iz,f,ieg) ! get rho, molarmass, ux, uy and uz
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

         wp(iux,l)= ux   ! Dirichlet, userbc
         wp(iuy,l)= uy   ! Dirichlet, userbc
         wp(iuz,l)= uz   ! Dirichlet, userbc
         wp(iph,l)  = phi  ! Dirichlet, userbc.
         wp(irho,l)=rho  ! Dirichlet, userbc
         up(1,l)= rho*phi
         up(2,l)= rho*ux*phi
         up(3,l)= rho*uy*phi
         up(4,l)= rho*uz*phi

         snx  = unx(l,1,f,e)
         sny  = uny(l,1,f,e)
         snz  = unz(l,1,f,e)
         mach = abs(ux*snx+uy*sny+uz*snz)/asnd ! better be uz=0.0 in 2D

         if (mach.lt.1.0) then

            wp(ipr,l)  = wm(ipr,l)
            wp(isnd,l) = asnd ! userbc should have set this to a(p-,rho+)
            wp(ithm,l) = temp ! userbc should have set this to T(p-,rho+)

         else ! supersonic inflow

            wp(ipr,l)  = pres ! Inflow state, userbc
            wp(isnd,l) = asnd ! Inflow state, userbc
            wp(ithm,l) = temp ! Inflow state, userbc

         endif

! userbc should have set e_internal to Dirichlet state (supersonic) or p-/(gm-1)
! here and only here is e_internal density-weighted.
         up(5,l) = e_internal+0.5*rho*(ux**2+uy**2+uz**2)
         up(5,l) = up(5,l)*phi

      enddo
      enddo
      enddo

      return
      end

!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
!******************************************************************************
!
! Purpose: set inflow boundary condition for one cell.
!
! Description: The subsonic boundary condition is based on the extrapolation of 
!   the Riemann invariant from the interior field (Holmes, D.G.: Inviscid 2D 
!   Solutions on Unstructured, Adaptive Grids. VKI Lecture Series 1989-06, 1989). 
!   The supersonic inflow boundary condition computes the conservative variables
!   from given velocity components, density and pressure.
!
! Input: bcOptType  = boundary treatment: subsonic, supersonic, or mixed
!        bcOptFixed = whether _computed_ inflow angle should be fixed or not
!        ptot       = given total pressure
!        ttot       = given total temperature
!        sx/y/zn    = components of normalized face vector (outward facing)
!        cpgas      = specific heat at constant pressure (boundary cell)
!        mm         = molecular mass at boundary cell
!        rl         = given density
!        ru/v/wl    = given velocity components
!
! Output: rr      = density at boundary | velocity compts inout from userbc
!         ru/v/wr = density * velocity components at boundary
!         rer     = density * total internal energy at boundary
!         pr      = pressure at boundary
!
! Notes: 
!   1. This condition is valid only for thermally and calorically perfect  
!      gas.
!   2. Important to avoid division by MakeNonZero(sl) when computing eta 
!      because that computation can become undefined for reservoir inflow
!      conditions, i.e., for a vanishing velocity vector. 
!
!******************************************************************************
!
! $Id: BcondInflowPerf.F90,v 1.1.1.1 2014/05/05 21:47:47 tmish Exp $
!
! Copyright: (c) 2002-2004 by the University of Illinois
!
!******************************************************************************

      SUBROUTINE BcondInflowPerf(bcOptType,bcOptFixed,ptot,ttot,ur,vr,wr
     >                          ,mach,sxn,syn,szn,cpgas,mm,rl,rul
     >                          ,rvl,rwl,rr,rur,rvr,rwr,rer,pr,ar,tr)
      IMPLICIT NONE

      integer BCOPT_SUBSONIC, BCOPT_MIXED, BCOPT_FIXED_NO
! faked flags
      parameter (BCOPT_SUBSONIC = 1)
      parameter (BCOPT_MIXED = 2)
      parameter (BCOPT_FIXED_NO = 2)
      real  MixtPerf_C_Co2GUVW, MixtPerf_C_DGP, MixtPerf_C_GRT,
     >  MixtPerf_Co2_CGUVW, MixtPerf_C2_GRT,MixtPerf_D_PRT,
     >  MixtPerf_Eo_DGPUVW, MixtPerf_Eo_DGPVm, MixtPerf_G_CpR,
     >  MixtPerf_P_GMaPo, MixtPerf_P_GPoTTo, MixtPerf_Po_GPTTo,
     >  MixtPerf_Po_CGPUVW, MixtPerf_R_M, MixtPerf_T_CGR,
     >  MixtPerf_T_GMaTo,  MixtPerf_Vm_C2Co2G
      external  MixtPerf_C_Co2GUVW, MixtPerf_C_DGP, MixtPerf_C_GRT,
     >  MixtPerf_Co2_CGUVW, MixtPerf_C2_GRT,MixtPerf_D_PRT,
     >  MixtPerf_Eo_DGPUVW, MixtPerf_Eo_DGPVm, MixtPerf_G_CpR,
     >  MixtPerf_P_GMaPo, MixtPerf_P_GPoTTo, MixtPerf_Po_GPTTo,
     >  MixtPerf_Po_CGPUVW, MixtPerf_R_M, MixtPerf_T_CGR,
     >  MixtPerf_T_GMaTo,  MixtPerf_Vm_C2Co2G


! ... parameters
      INTEGER bcOptFixed,bcOptType!, INTENT(IN) ::
!          .      .     .                  .   .   .   from userbc
      REAL cpgas, mach, mm, sxn, syn, szn, ur, vr, wr,!, INTENT(IN) ::
     >                        ptot, rl, rul, rvl, rwl, ttot

      REAL rer, rr, rur, rvr, rwr, pr!, INTENT(OUT) ::

      REAL al, ar, a02, cp, disc, eta, g, gm1, igm1, ql, rgas, Rm,
     >            sl, sr, tr, ul,  vl,  wl

!******************************************************************************
! gas properties

      rgas = MixtPerf_R_M(mm)
      g    = MixtPerf_G_CpR(cpgas,rgas)

! subsonic or mixed -----------------------------------------------------------   

      IF ( bcOptType.eq.BCOPT_SUBSONIC.OR.bcOptType.eq.BCOPT_MIXED) THEN
         gm1  = g - 1.0
         igm1 = 1.0/gm1

         ul = rul/rl
         vl = rvl/rl
         wl = rwl/rl

         a02 = MixtPerf_C2_GRT(g,rgas,ttot)

         al = MixtPerf_C_Co2GUVW(a02,g,ul,vl,wl) ! make al consistent with a02
         ql = ul*sxn + vl*syn + wl*szn

! - subsonic

         IF ( ABS(ql) .lt. al ) THEN
            sl = SQRT(ul*ul + vl*vl + wl*wl)
         
            IF ( bcOptFixed .eq. BCOPT_FIXED_NO ) THEN 
               IF ( sl .gt. 1.0E-6 ) THEN ! Avoid ill-defined angle computation
                  eta = ql/sl        
               ELSE 
                  eta = -1.0
               END IF ! sl
            ELSE 
               eta = -1.0
            END IF ! bcOptFixed

            Rm   = al - 0.5*gm1*ql
            disc = 0.5*gm1*eta**2*
     >       (a02/(Rm*Rm)*(1.0 + 0.5*gm1*eta**2) - 1.0)     
        
            IF ( disc .lt. 0.0 ) THEN ! discriminant cannot be negative
               ar = SQRT(a02)
               tr = ttot
               pr = ptot
               sr = 0.0
            ELSE
               ar = Rm/(1.0 + 0.5*gm1*eta*eta)*(1.0+SQRT(disc))
               tr = MixtPerf_T_CGR(ar,g,rgas)
               pr = MixtPerf_P_GPoTTo(g,ptot,tr,ttot)
               sr = MixtPerf_Vm_C2Co2G(ar*ar,a02,g)
            END IF ! disc    
                     
            rr = MixtPerf_D_PRT( pr,rgas,tr )

            rer = rr*MixtPerf_Eo_DGPVm(rr,g,pr,sr)
            rur = rr*ur
            rvr = rr*vr
            rwr = rr*wr

! - supersonic

         ELSE
            IF ( bcOptType .eq. BCOPT_MIXED ) THEN
               pr = mixtPerf_P_GMaPo(g,mach,ptot)
               tr = mixtPerf_T_GMaTo(g,mach,ttot)
               rr = mixtPerf_D_PRT(pr,rgas,tr)
               ar = mixtPerf_C_GRT(g,rgas,tr)

               rer = rr*MixtPerf_Eo_DGPUVW(rr,g,pr,ur,vr,wr)
               rur = rr*ur
               rvr = rr*vr
               rwr = rr*wr
            END IF ! bcOptType
         END IF ! ql < al

! supersonic ------------------------------------------------------------------

      ELSE ! bcOptType .eq. BCOPT_SUPERSONIC
         pr = mixtPerf_P_GMaPo(g,mach,ptot)
         tr = mixtPerf_T_GMaTo(g,mach,ttot)
         rr = mixtPerf_D_PRT(pr,rgas,tr)
         ar = mixtPerf_C_GRT(g,rgas,tr)

         rer = rr*MixtPerf_Eo_DGPUVW(rr,g,pr,ur,vr,wr)
         rur = rr*ur
         rvr = rr*vr
         rwr = rr*wr
      END IF ! bcOptType

      return
      end

!******************************************************************************
!
! RCS Revision history:
!
! $Log: BcondInflowPerf.F90,v $
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
! Revision 1.1  2004/12/01 16:47:56  haselbac
! Initial revision after changing case
!
! Revision 1.7  2004/01/29 22:52:36  haselbac
! Added bcOptFixed, fixed bug, clean-up
!
! Revision 1.6  2003/12/04 03:22:56  haselbac
! Fixed bug in formulation, added partial fix for eta
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
