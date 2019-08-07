C> @file intpdiff.f interpolation and differentiation routines not already provided
C> by nek5000
!--------------------------------------------------------------------
! JH061914 propose name change to intpdiff since that is what is in
!          here.
! JH081816 went to local_grad. redimensioned gradu. wish I could use
!          gradm11, but stride
! JH080719 compute_gradients superceded by compute_gradients_contra
!          until nek5000 community teaches me the best way of handling
!          lxd. for now, contravariant metrics are stored in rx,
!          not in /giso1/. Hence compute_gradients_contra
!--------------------------------------------------------------------

      subroutine compute_gradients(e)
      include 'SIZE'
      include 'INPUT'
      include 'DXYZ'
      include 'GEOM'
      include 'SOLN'
      include 'CMTDATA'
      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ud(ldd),tu(ldd)

      integer eq,e

!     !  Compute d/dx, d/dy and d/dz of all the cons vars

      nxy1  = lx1*ly1
      nyz1  = ly1*lz1
      nxyz1 = lx1*ly1*lz1
      m0 = lx1-1

      do eq=1,toteq
         call invcol3(ud,u(1,1,1,eq,e),phig(1,1,1,e),nxyz1)

         if (if3d) then
            call local_grad3(ur,us,ut,ud,m0,1,dxm1,dxtm1)
            do i=1,nxyz1
               gradu(i,eq,1) = jacmi(i,e)*(rxm1(i,1,1,e)*ur(i)+
     >                                     sxm1(i,1,1,e)*us(i)+
     >                                     txm1(i,1,1,e)*ut(i))
            enddo
            do i=1,nxyz1
               gradu(i,eq,2) = jacmi(i,e)*(rym1(i,1,1,e)*ur(i)+
     >                                     sym1(i,1,1,e)*us(i)+
     >                                     tym1(i,1,1,e)*ut(i))
            enddo
            do i=1,nxyz1
               gradu(i,eq,3) = jacmi(i,e)*(rzm1(i,1,1,e)*ur(i)+
     >                                     szm1(i,1,1,e)*us(i)+
     >                                     tzm1(i,1,1,e)*ut(i))
            enddo

         else

            call local_grad2(ur,us   ,ud,m0,1,dxm1,dxtm1)
            do i=1,nxyz1
               gradu(i,eq,1) = jacmi(i,e)*(rxm1(i,1,1,e)*ur(i)+
     >                                     sxm1(i,1,1,e)*us(i))
            enddo
            do i=1,nxyz1
               gradu(i,eq,2) = jacmi(i,e)*(rym1(i,1,1,e)*ur(i)+
     >                                     sym1(i,1,1,e)*us(i))
            enddo

         endif

      enddo ! equation loop

      return
      end

!--------------------------------------------------------------------
! JH080719 compute_gradients_contra computes gradU using
!          contravariant metrics stored in rx
!          FIGURE OUT IF GRADU*=JACMI or NOT!!!
!--------------------------------------------------------------------
      subroutine compute_gradients_contra(e)
      include 'SIZE'
      include 'INPUT'
      include 'DXYZ'
      include 'GEOM'
      include 'SOLN'
      include 'CMTDATA'
      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ud(ldd),tu(ldd)

      integer eq,e

!     !  Compute d/dx, d/dy and d/dz of all the cons vars

      nxyz1 = lx1*ly1*lz1
      m0 = lx1-1

      do eq=1,toteq
         call invcol3(ud,u(1,1,1,eq,e),phig(1,1,1,e),nxyz1)

         if (if3d) then
            call local_grad3(ur,us,ut,ud,m0,1,dxm1,dxtm1)
            do j=1,ldim !xyz
               j0=j+0 ! r_xyz
               j3=j+3 ! s_xyz
               j6=j+6 ! t_xyz
               do i=1,nxyz1
                  gradu(i,eq,j) = rx(i,j0,e)*ur(i)+
     >                            rx(i,j3,e)*us(i)+
     >                            rx(i,j6,e)*ut(i)
               enddo
            enddo
         else
            call local_grad2(ur,us,ud,m0,1,dxm1,dxtm1)
            do j=1,ldim !xy
               j0=j+0 ! r_xy
               j2=j+2 ! s_xy
               do i=1,nxyz1
                  gradu(i,eq,j) = rx(i,j0,e)*ur(i)+
     >                            rx(i,j2,e)*us(i)
               enddo
            enddo
         endif

      return
      end

!-----------------------------------------------------------------------

      subroutine set_dealias_face

!-----------------------------------------------------------------------
! JH111815 needed for face Jacobian and surface integrals
!-----------------------------------------------------------------------

      include 'SIZE'
      include 'INPUT' ! for if3d
      include 'GEOM'  ! for ifgeom
      include 'TSTEP' ! for istep
      include 'WZ'    ! for wxm1
      include 'DG'    ! for facewz

      integer ilstep
      save    ilstep
      data    ilstep /-1/

      if (.not.ifgeom.and.ilstep.gt.1) return  ! already computed
      if (ifgeom.and.ilstep.eq.istep)  return  ! already computed
      ilstep = istep

      call zwgl(zptf,wgtf,lxd)

      if (if3d) then
         k=0
         do j=1,lyd
         do i=1,lxd
            k=k+1
            wghtf(k)=wgtf(i)*wgtf(j)
         enddo
         enddo
      else
         call copy(wghtf,wgtf,lxd)
      endif

      return
      end
!-----------------------------------------------------------------------

      subroutine cmt_metrics(istp)
! compute freestream-preserving metrics $Ja^i$ for transforming fluxes
! F to F~ in a contravariant frame according to
! Kopriva (2006)
! Follows methodology in FLUXO (github.com/project-fluxo/fluxx.git)

! I REALLY want to duplicate get_dgll_ptr, make my own version of
! /dgradl/, set lxd=1, and replace rx with ja(lx1*ly1*lz1,ldim*ldim,lelv) 
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'DG'
      include 'WZ'
      include 'MASS'
      include 'CMTDATA'

      parameter (ldg=lx1**3,lwkd=4*lx1*lx1)
      common /dgradl/d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     > ,wkd(lwkd)
      real jgl,jgt
      common /ctmp1/ ur(ldg),us(ldg),ut(ldg),ud(ldg)
      parameter (ngeo=2,ngeoref=3*(ngeo-1)+1)
      parameter (ld=2*lxd,ldw=2*(ld**ldim)) ! sigh
      common /igrad/ pd    (0:ld*ld) ! reap 
     $             , pdg   (0:ld*ld) ! where you
     $             , pjgl  (0:ld*ld) ! do not sow
      integer pd , pdg , pjgl
      common /ctmp0/ w(ld**ldim,2)
      common /scrns/xyz(ngeo,ngeo,ngeo,3),xrref(ngeoref**ldim,3,3),
     >              xr(lx1*ly1*lz1,3,3),jref(ngeoref**ldim)
     >             ,screlm(lx1*ly1*lz1),jaface(lx1*lz1,2*ldim,ldim,ldim)
      real jref,jaface
      integer e,f
      integer ndir(6),nsgn(6)
      data ndir /2,1,2,1,3,3/
      data nsgn /-1,1,1,-1,-1,1/

      integer ilstep
      save    ilstep
      data    ilstep /-1/

      if (.not.ifgeom.and.ilstep.gt.1) return  ! already computed
      if (ifgeom.and.ilstep.eq.istp)  return  ! already computed
      ilstep = istp
      write(6,*) 'welcome to cmt_metrics'

      n=ngeoref**ldim
      nxyz=lx1*ly1*lz1

      call get_int_gll2gll(jgeo2ref,ngeo,ngeoref)
      call get_int_gll2gll(jgeo2n,ngeo,lx1)

      if (lx1.lt.ngeoref) then

         ij = lx1 + ld*(ngeoref-1)
         jref2n = pjgl(ij)
         if (jref2n.eq.0) then
            nstore   = pjgl(0)
            pjgl(ij) = nstore+1
            nstore   = nstore + lx1*ngeoref
            pjgl(0)  = nstore
            jref2n   = pjgl(ij)
            nwrkd = lx1 + ngeoref
            call lim_chk(nstore,ldg ,'jgl  ','ldg  ','getgll2gll')
            call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','getgll2gll')

            call proj_legmodal(jgl(jref2n),ngeoref,lx1)
            call transpose(jgt(jref2n),ngeoref,jgl(jref2n),lx1)

         endif

      else!, lx1>ngeoref
         call get_int_gll2gll(jref2n,ngeoref,lx1)
      endif

      call get_dgll_ptr(igeo,ngeo,ngeo)
      call get_dgll_ptr(ilx1,lx1,lx1)

! JH060418 strong-form derivative matrix for 2-point fluxes + surface
      call copy (dstrong,d(ilx1),lx1**2)
      call cmult(dstrong,2.0,    lx1**2)
      dstrong(1,1)     = 2.0*d(ilx1)         +1.0/wxm1(1)
      dstrong(lx1,lx1) = 2.0*d(ilx1+lx1**2-1)-1.0/wxm1(lx1)
      call transpose(dstrongt,lx1,dstrong,lx1)

      do e=1,nelt
         m=ngeo-1
! check for curved faces here someday and set ngeo on the fly,
! tucking all the get* calls above into the element loop
         if (m .eq. 1) then 
            call xyztriv(xyz(1,1,1,1),xyz(1,1,1,2),xyz(1,1,1,3),
     >               ngeo,ngeo,ngeo,e)
         elseif (m .eq. 2) then
            call exitti('need an xyzquad call in metrics plz$',m)
         endif

         if (if3d) then
!        l=0
            do i=1,ldim ! x_i
               call local_grad3(ur,us,ut,xyz(1,1,1,i),m,1,d(igeo),
     >                                                   dt(igeo))
! interpolate for Jacobian
               call specmpn(xrref(1,i,1),ngeoref,ur,ngeo,
     >                   jgl(jgeo2ref),jgt(jgeo2ref),if3d,w,ldw)
               call specmpn(xrref(1,i,2),ngeoref,us,ngeo,
     >                   jgl(jgeo2ref),jgt(jgeo2ref),if3d,w,ldw)
               call specmpn(xrref(1,i,3),ngeoref,ut,ngeo,
     >                   jgl(jgeo2ref),jgt(jgeo2ref),if3d,w,ldw)
! interpolate for metrics
               if (lx1 .gt. ngeo) then
                  call specmpn(xr(1,i,1),lx1,ur,ngeo,
     >                   jgl(jgeo2n),jgt(jgeo2n),if3d,w,ldw)
                  call specmpn(xr(1,i,2),lx1,us,ngeo,
     >                   jgl(jgeo2n),jgt(jgeo2n),if3d,w,ldw)
                  call specmpn(xr(1,i,3),lx1,ut,ngeo,
     >                   jgl(jgeo2n),jgt(jgeo2n),if3d,w,ldw)
               else
                  call exitti('copy rstxyzm1 into rx and copy to xr$',
     >                        ngeo)
               endif
            enddo ! x_i

! compute Jacobian
            call rzero(jref,n) 
            call addcol4(jref,xrref(1,1,1),xrref(1,2,2),xrref(1,3,3),n)
            call addcol4(jref,xrref(1,1,3),xrref(1,2,1),xrref(1,3,2),n)
            call addcol4(jref,xrref(1,1,2),xrref(1,2,3),xrref(1,3,1),n)
            call subcol4(jref,xrref(1,1,1),xrref(1,2,3),xrref(1,3,2),n)
            call subcol4(jref,xrref(1,1,2),xrref(1,2,1),xrref(1,3,3),n)
            call subcol4(jref,xrref(1,1,3),xrref(1,2,2),xrref(1,3,1),n)
! project back
            call specmpn(jacm1(1,1,1,e),lx1,jref,ngeoref,jgl(jref2n),
     >                jgt(jref2n),if3d,w,ldw)
            call chkjac(jacm1(1,1,1,e),nxyz,e,xm1(i,1,1,e),ym1(i,1,1,e),
     >               zm1(i,1,1,e),ldim,ierr)
            if (ierr .ne. 0) then
               call exitti('failed jacobian check in element $',e)
            endif

!call intp_rstd(rx(1,1,e),rxm1(1,1,1,e),lx1,lxd,if3d,0) 
!call intp_rstd(rx(1,2,e),rym1(1,1,1,e),lx1,lxd,if3d,0) 
!call intp_rstd(rx(1,3,e),rzm1(1,1,1,e),lx1,lxd,if3d,0) 
!call intp_rstd(rx(1,4,e),sxm1(1,1,1,e),lx1,lxd,if3d,0) 
!call intp_rstd(rx(1,5,e),sym1(1,1,1,e),lx1,lxd,if3d,0) 
!call intp_rstd(rx(1,6,e),szm1(1,1,1,e),lx1,lxd,if3d,0) 
!call intp_rstd(rx(1,7,e),txm1(1,1,1,e),lx1,lxd,if3d,0) 
!call intp_rstd(rx(1,8,e),tym1(1,1,1,e),lx1,lxd,if3d,0) 
!call intp_rstd(rx(1,9,e),tzm1(1,1,1,e),lx1,lxd,if3d,0)
! COMPUTE R and store temporarily in ja=rx
            j0=1
            do j=1,3 ! R_j in   x          y            z directions
               call vcross(rx(1,j0,e),rx(1,j0+1,e),rx(1,j0+2,e),
     >                     xr(1,1,j), xr(1,2,j),   xr(1,3,j),
     >                    xm1(1,1,1,e),ym1(1,1,1,e),zm1(1,1,1,e),nxyz)
               j0=j0+3
            enddo
            call cmult(rx(1,1,e),0.5,lxd*lyd*lzd*ldim*ldim)

! xr becomes scratch for Ja. I don't feel like writing my own
! nested tensor product, so I've elected to do 50% more
! derivatives needlessly
            j0=1
            m2=lx1-1
! j=1, R_j; rst loop hardcoded until I write my own tensor prod.
            do i=1,3 !      xyz
               j1=j0+(i-1)
               call local_grad3(ur,us,ut,rx(1,j1,e),m2,1,d(ilx1),
     >                                                  dt(ilx1))
! j0=1,2,3 ut=dR1/dxi_3 goes into -ja2
               call chsign(ut,nxyz)
               call copy(xr(1,i,2),ut,nxyz)
! j0=1,2,3 us=dR1/dxi_2 goes into +ja3
               call copy(xr(1,i,3),us,nxyz)
            enddo

            j0=j0+3
! j=2
            do i=1,3 !      xyz
               j1=j0+(i-1)
               call local_grad3(ur,us,ut,rx(1,j1,e),m2,1,d(ilx1),
     >                                                  dt(ilx1))
! j0=4,5,6 ut=dR2/dxi_3 goes into +ja1
               call copy(xr(1,i,1),ut,nxyz)
! j0=4,5,6 ur=dR2/dxi_1 goes into -ja3
               call sub2(xr(1,i,3),ur,nxyz)
            enddo

            j0=j0+3
! j=3
            do i=1,3 !      xyz
               j1=j0+(i-1)
               call local_grad3(ur,us,ut,rx(1,j1,e),m2,1,d(ilx1),
     >                                                  dt(ilx1))
! j0=7,8,9 us=dR3/dxi_2 goes into -ja1
               call sub2(xr(1,i,1),us,nxyz)
! j0=7,8,9 ur=dR3/dxi_1 goes into +ja2
               call add2(xr(1,i,2),ur,nxyz)
            enddo

! overwrite ja=rx
! REAAAALLY want to get rid of lxd here
            j0=0
            do j=1,3
               do i=1,3
                  j0=j0+1
                  call copy(rx(1,j0,e),xr(1,i,j),nxyz)
               enddo
            enddo

         else ! 2D 

            do i=1,ldim
               call local_grad2(ur,us,xyz(1,1,1,i),m,1,d(igeo),dt(igeo))
! interpolate for Jacobian
               call specmpn(xrref(1,i,1),ngeoref,ur,ngeo,
     >                   jgl(jgeo2ref),jgt(jgeo2ref),if3d,w,ldw)
               call specmpn(xrref(1,i,2),ngeoref,us,ngeo,
     >                   jgl(jgeo2ref),jgt(jgeo2ref),if3d,w,ldw)
            enddo
            call rzero(jref,n) 
            call addcol3(jref,xrref(1,1,1),xrref(1,2,2),n)
            call subcol3(jref,xrref(1,1,2),xrref(1,2,1),n)
! project back
            call copy(screlm,jacm1(1,1,1,e),nxyz)
            call specmpn(jacm1(1,1,1,e),lx1,jref,ngeoref,jgl(jref2n),
     >                jgt(jref2n),if3d,w,ldw)
            call chkjac(jacm1(1,1,1,e),nxyz,e,xm1(1,1,1,e),ym1(1,1,1,e),
     >               zm1(1,1,1,e),ldim,ierr)
            if (ierr .ne. 0) then
               call exitti('failed jacobian check in element $',e)
            endif
! cross product form is fine in 2D; no need to recompute metrics, but
! contravariant Ja do NOT have quadrature weights in them.
            call copy(rx(1,1,e),rxm1(1,1,1,e),nxyz) 
            call copy(rx(1,2,e),rym1(1,1,1,e),nxyz) 
            call copy(rx(1,3,e),sxm1(1,1,1,e),nxyz) 
            call copy(rx(1,4,e),sym1(1,1,1,e),nxyz) 

         endif ! if3d

         j0=0
         do j=1,ldim    ! r_j
            do i=1,ldim ! x_i
               j0=j0+1
               call full2face_cmt(1,lx1,ly1,lz1,iface_flux(1,e),
     >                            jaface(1,1,i,j),rx(1,j0,e))
            enddo
         enddo

         do f=1,2*ldim
            if (nsgn(f).lt.0) then
               wght=1.0/wxm1(1)
            else
               wght=1.0/wxm1(lx1)
            endif
            l=0
            do iz=1,lz1
            do ix=1,lx1
               anew=0.0
               l=l+1
               do i=1,ldim
                  anew=anew+jaface(l,f,i,ndir(f))**2
               enddo
               anew=sqrt(anew)
               jface(ix,iz,f,e)=anew*wght
               unx(ix,iz,f,e)=jaface(l,f,1,ndir(f))/anew*nsgn(f)
               uny(ix,iz,f,e)=jaface(l,f,2,ndir(f))/anew*nsgn(f)
               if (if3d) then
                  unz(ix,iz,f,e)=jaface(l,f,3,ndir(f))/anew*nsgn(f)
               endif
! TODO: check again on deformed meshes.

! we're taking over this town
               area(ix,iz,f,e)=anew*wxm1(ix)*wzm1(iz)

            enddo
            enddo
         enddo

! might as well recompute the mass matrix
         call col3(bm1(1,1,1,e),jacm1(1,1,1,e),w3m1,nxyz)
! and we need this for strong-form discontinuous surface fluxes computed
! from contravariant volume fluxes
! JH061318 Actually, we don't need this for strong form. Dstrong for life.
         do j=1,lz1
         do i=1,lx1
            w2m1(i,j)=wxm1(i)*wzm1(j)
         enddo
         enddo

      enddo ! e=1,nelt

      call invers2(jacmi,jacm1,nxyz*nelt)

      return
      end

!-----------------------------------------------------------------------

      subroutine xyztriv(xl,yl,zl,nxl,nyl,nzl,e)
c     reorder corner vertices lexicographically. better living
c     through paranoia. yes, all those J better be identities.

      include 'SIZE'
      include 'INPUT'

      real xl(nxl,nyl,nzl),yl(nxl,nyl,nzl),zl(nxl,nyl,nzl)
      integer e

c   Preprocessor Corner notation:      Symmetric Corner notation:
c
c           4+-----+3    ^ s                    3+-----+4    ^ s
c           /     /|     |                      /     /|     |
c          /     / |     |                     /     / |     |
c        8+-----+7 +2    +----> r            7+-----+8 +2    +----> r
c         |     | /     /                     |     | /     /
c         |     |/     /                      |     |/     /
c        5+-----+6    t                      5+-----+6    t

      integer indx(8)
      save    indx
      data    indx / 1,2,4,3,5,6,8,7 /

      parameter (ldw=4*lx1*ly1*lz1)
      common /ctmp0/ xcb(2,2,2),ycb(2,2,2),zcb(2,2,2),w(ldw)

      common /cxyzl/ zgml(2,3),jx (2*2),jy (2*2),jz (2*2)
     $                          ,jxt(2*2),jyt(2*2),jzt(2*2)
     $                          ,zlin(2)
      real jx,jy,jz,jxt,jyt,jzt

      zlin(1) = -1
      zlin(2) =  1
      do j=1,3
         call copy(zgml(1,j),zlin,nxl)
      enddo

      k = 1
      do i=1,nxl
         call fd_weights_full(zgml(i,1),zlin,1,0,jxt(k))
         call fd_weights_full(zgml(i,2),zlin,1,0,jyt(k))
         call fd_weights_full(zgml(i,3),zlin,1,0,jzt(k))
         k=k+2
      enddo
      call transpose(jx,nxl,jxt,2)

      ldim2 = 2**ldim
      do ix=1,ldim2          ! Convert prex notation to lexicographical
         i=indx(ix)
         xcb(ix,1,1)=xc(i,e)
         ycb(ix,1,1)=yc(i,e)
         zcb(ix,1,1)=zc(i,e)
      enddo

c     Map R-S-T space into physical X-Y-Z space.

      ! NOTE:  Assumes nxl=nyl=nzl !

      call tensr3(xl,nxl,xcb,2,jx,jyt,jzt,w)
      call tensr3(yl,nxl,ycb,2,jx,jyt,jzt,w)
      call tensr3(zl,nxl,zcb,2,jx,jyt,jzt,w)

      return
      end

!-----------------------------------------------------------------------

      subroutine get_int_gll2gll (ip,mx,md) ! GLL-->GLL pointer

c     Get pointer to jgl() for interpolation pair (mx,md)

      include 'SIZE'

      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgradl/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt
c
      parameter (ld=2*lxd)
      common /igrad/ pd    (0:ld*ld) ! reap 
     $             , pdg   (0:ld*ld) ! where you
     $             , pjgl  (0:ld*ld) ! do not sow
      integer pd , pdg , pjgl
c
      ij = md + ld*(mx-1)
      ip = pjgl(ij)
c
      if (ip.eq.0) then
c
         nstore   = pjgl(0)
         pjgl(ij) = nstore+1
         nstore   = nstore + md*mx
         pjgl(0)  = nstore
         ip       = pjgl(ij)
c
         nwrkd = mx + md
         call lim_chk(nstore,ldg ,'jgl  ','ldg  ','getgll2gll')
         call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','getgll2gll')
c
         call gen_int_gll2gll(jgl(ip),jgt(ip),md,mx,wkd)
      endif
c
      return
      end

!-----------------------------------------------------------------------

      subroutine gen_int_gll2gll(jgl,jgt,mp,np,w)

c     Generate interpolation from np GLL points to mp GLL points

c        jgl  = interpolation matrix, mapping from velocity nodes to pressure
c        jgt  = transpose of interpolation matrix
c        w    = work array of size (np+mp)

c        np   = number of points on GLL grid
c        mp   = number of points on GL  grid


      real jgl(mp,np),jgt(np*mp),w(*)

      iz = 1
      id = iz + np

      call zwgll (w(iz),jgt,np)
      call zwgll (w(id),jgt,mp)

      n  = np-1
      do i=1,mp
         call fd_weights_full(w(id+i-1),w(iz),n,0,jgt)
         do j=1,np
            jgl(i,j) = jgt(j)                  !  Interpolation matrix
         enddo
      enddo

      call transpose(jgt,np,jgl,mp)

      return
      end

!-----------------------------------------------------------------------

      subroutine proj_legmodal(promat,nin,nout)

      parameter (lm=84)
      parameter (lm2=lm*lm)
      real      rmult(lm),vin(lm2),vout(lm2)
      real      z(lm),w(lm)
      integer   indr(lm),indc(lm),ipiv(lm)

      real promat(nout,nin)

      if (nin.le.nout) then
         call exitti('proj_legmodal is not for N_in < N_out$',
     >              nin-nout)
      endif

      call zwgll(z,w,nin)
      call vandermonde_legendre(vin,z,nin)
      call zwgll(z,w,nout)
      call vandermonde_legendre(vout,z,nout)
      call gaujordf  (vin,nin,nin,indr,indc,ipiv,ierr,rmult)
      call discard_rows(vin,vin,nout,nin)
      call mxm  (vout ,nout,vin,nout,promat,nin)
      return
      end

!-----------------------------------------------------------------------

      subroutine discard_rows(trunc,matrix,nsmall,nlarge)
      real trunc(nsmall,nlarge),matrix(nlarge,nlarge)
      if (nsmall .ge. nlarge) then
         call exitti('discard_rows needs more rows than $',nsmall)
      endif
      do j=1,nlarge
         do i=1,nsmall
            trunc(i,j)=matrix(i,j)
         enddo
      enddo
      return
      end

!-----------------------------------------------------------------------

      subroutine vandermonde_legendre(v,z,nx)
      parameter (lm=84)
      parameter (lm2=lm*lm)
      real      vtmp(lm2),Lj(lm)
      real v(*),z(nx)

      if (nx.gt.lm) then
         call exitti('OK vandermonde_legendre can''t handle n=$',
     >               n)
      endif

      kj = 0
      n  = nx-1
      do j=1,nx
         call legendre_poly(Lj,z(j),n) ! is PNLEG better?
         do k=1,nx
            kj = kj+1
            coef=2.0/(2.0*(k-1)+1)
            vtmp(kj) = Lj(k)/sqrt(coef)
         enddo
      enddo
      call transpose (v,nx,vtmp,nx)

      return
      end

!-----------------------------------------------------------------------

      subroutine gradm11_t(grad,uxyz,csgn,e) ! grad is incremented, not overwritten
c     source: . gradm1, from navier5.f 
c             . gradm1, from navier5.f 
c     Compute divergence^T of ux,uy,uz -- mesh 1 to mesh 1 (vel. to vel.)
!     single element, but references jacmi and the metrics

      include 'SIZE'
      include 'DXYZ'
      include 'GEOM'
      include 'INPUT'

      parameter (lxyz=lx1*ly1*lz1)
      real grad(lxyz),uxyz(lxyz,ldim)

      common /ctmp1/ ur(lxyz),us(lxyz),ut(lxyz),ud(lxyz),tmp(lxyz)
      real ur,us,ut,tmp

      integer e

      nxyz = lx1*ly1*lz1
      call rzero(ud,nxyz)

      N = lx1-1
      if (if3d) then

         do i=1,lxyz
            ur(i) = jacmi(i,e)*(uxyz(i,1)*rxm1(i,1,1,e)
     >                        + uxyz(i,2)*rym1(i,1,1,e)
     >                        + uxyz(i,3)*rzm1(i,1,1,e) )
            us(i) = jacmi(i,e)*(uxyz(i,1)*sxm1(i,1,1,e)
     >                        + uxyz(i,2)*sym1(i,1,1,e)
     >                        + uxyz(i,3)*szm1(i,1,1,e) )
            ut(i) = jacmi(i,e)*(uxyz(i,1)*txm1(i,1,1,e)
     >                        + uxyz(i,2)*tym1(i,1,1,e)
     >                        + uxyz(i,3)*tzm1(i,1,1,e) )
         enddo
         call local_grad3_t(ud,ur,us,ut,N,1,dxm1,dxtm1,tmp)
      else
         do i=1,lxyz
            ur(i) =jacmi(i,e)*(uxyz(i,1)*rxm1(i,1,1,e)
     >                       + uxyz(i,2)*rym1(i,1,1,e) )
            us(i) =jacmi(i,e)*(uxyz(i,1)*sxm1(i,1,1,e)
     >                       + uxyz(i,2)*sym1(i,1,1,e) )
         enddo
         call local_grad2_t(ud,ur,us,N,1,dxm1,dxtm1,tmp)
      endif
      call cmult(ud,csgn,nxyz)
      call add2(grad,ud,nxyz)

      return
      end

!-----------------------------------------------------------------------

      subroutine gradm1_t(u,ux,uy,uz)
! JH082516 torn bleeding from Lu's dgf3.f. someday you are going to have
!          to vectorize cmt-nek properly.
c     source: . gradm1, from navier5.f 
c             . gradm1, from navier5.f 
c
c     Compute divergence of ux,uy,uz -- mesh 1 to mesh 1 (vel. to vel.)
c
      include 'SIZE'
      include 'DXYZ'
      include 'GEOM'
      include 'INPUT'
c
      parameter (lxyz=lx1*ly1*lz1)
      real ux(lxyz,*),uy(lxyz,*),uz(lxyz,*),u(lxyz,*)

      common /ctmp1/ ur(lxyz),us(lxyz),ut(lxyz),tmp(lxyz)
      real ur,us,ut,tmp

      integer e

      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelt

      N = lx1-1
      do e=1,nelt
         if (if3d) then

            do i=1,lxyz
               ur(i) = jacmi(i,e)*(ux(i,e)*rxm1(i,1,1,e)
     $                           + uy(i,e)*rym1(i,1,1,e)
     $                           + uz(i,e)*rzm1(i,1,1,e) )
               us(i) = jacmi(i,e)*(ux(i,e)*sxm1(i,1,1,e)
     $                           + uy(i,e)*sym1(i,1,1,e)
     $                           + uz(i,e)*szm1(i,1,1,e) )
               ut(i) = jacmi(i,e)*(ux(i,e)*txm1(i,1,1,e)
     $                           + uy(i,e)*tym1(i,1,1,e)
     $                           + uz(i,e)*tzm1(i,1,1,e) )
            enddo
            call local_grad3_t(u,ur,us,ut,N,e,dxm1,dxtm1,tmp)
         else
            do i=1,lxyz
               ur(i) =jacmi(i,e)*(ux(i,e)*rxm1(i,1,1,e)
     $                          + uy(i,e)*rym1(i,1,1,e) )
               us(i) =jacmi(i,e)*(ux(i,e)*sxm1(i,1,1,e)
     $                          + uy(i,e)*sym1(i,1,1,e) )
            enddo
            call local_grad2_t(u,ur,us,N,e,dxm1,dxtm1,tmp)
         endif
      enddo
c
      return
      end
