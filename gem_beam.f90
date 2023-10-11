
subroutine ppushb(n)

  use gem_com
  use gem_equil
  implicit none
  real :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
  real :: wx0,wx1,wy0,wy1,wz0,wz1,dum,vxdum,dum1,bstar
  INTEGER :: m,i,j,k,l,n,ns
  INTEGER :: np_old,np_new
  real :: rhog,vfac,kap,vpar,pidum,kaptp,kapnp,xnp
  real :: b,th,r,enerb,cost,sint,qr,laps,sz,ter,dtp
  real :: xt,xs,yt,xdot,ydot,zdot,pzdot,edot,pzd0,vp0
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,grdgtp
  real :: grp,gxdgyp,rhox(4),rhoy(4),psp,pzp,vncp,vparspp,psip2p,bdcrvbp,curvbzp,dipdrp
  integer :: mynopi
  real :: fdum,gdum,fisrcp,dnisrcp,avwixepsp,fovg,avwixezp,dnisrczp
  !MVPT
  real :: dapsdtp, delbxhp, delbyhp, dahdzp, aparhp,zdot0,zdot1

  mynopi = 0
  nopb = 0

#ifdef OPENMP
!$omp parallel do default(shared) &
!$omp private(i,j,k,l,r,th,wx0,wx1,wy0,wy1,wz0,wz1,dbdrp,dbdtp,grcgtp,bfldp,radiusp,dydrp,qhatp,grp,gxdgyp,curvbzp,bdcrvbp) &
!$omp private(grdgtp,fp,jfnp,psipp,psp,ter,kaptp,kapnp,xnp,vncp,vparspp,psip2p,dipdrp,b,pzp,rhog,rhox,rhoy) &
!$omp private(phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp,xs,xt,yt,vfac,vp0,vpar,bstar,enerb,kap,dum1,vxdum,xdot,ydot,zdot) &
!$omp private(pzd0,pzdot,edot,dum,laps,qr) &
!$omp private(fdum,gdum,fisrcp,dnisrcp,avwixepsp,fovg,dtp,avwixezp,dnisrczp) &
!$omp private(dapsdtp,delbxhp,delbyhp,dahdzp,aparhp,zdot0,zdot1) &    !MVPT
!$omp reduction(+: mynopi)
#endif

#ifdef OPENACC
!$acc parallel                                                                                                                                                                         
!$acc loop gang vector private(rhox,rhoy)                                                                                                                                              
#endif  

  do m=1,mmb
     r=x2b(m)-0.5*lx+lr0
     k = int(z2b(m)/delz)
     wz0 = ((k+1)*delz-z2b(m))/delz
     wz1 = 1-wz0
     th = wz0*thfnz(k)+wz1*thfnz(k+1)

     i = int((r-rin)/dr)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     k = int((th+pi)/dth)
     wz0 = (-pi+(k+1)*dth-th)/dth
     wz1 = 1.-wz0
     dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
          +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
     dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
          +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1) 
     grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
          +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 
     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
     radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
          +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
     dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
          +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1) 
     qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
          +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1) 
     grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
          +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1) 
     gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
          +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1) 

     curvbzp = wx0*wz0*curvbz(i,k)+wx0*wz1*curvbz(i,k+1) &
          +wx1*wz0*curvbz(i+1,k)+wx1*wz1*curvbz(i+1,k+1) 
     bdcrvbp = wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
          +wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1) 
     grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
          +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1) 

     fp = wx0*f(i)+wx1*f(i+1)        
     jfnp = wz0*jfn(k)+wz1*jfn(k+1)
     psipp = wx0*psip(i)+wx1*psip(i+1)        
     psp = wx0*psi(i)+wx1*psi(i+1)        
     vncp = 0.
     psip2p = wx0*psip2(i)+wx1*psip2(i+1)        
     dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)        
     b=1.-tor+tor*bfldp
     pzp = mbeam*u2b(m)/b*fp/br0-qbeam*psp/br0

     rhog=sqrt(2.*b*mub(m)*mbeam)/(qbeam*b)*iflr

     rhox(1) = rhog*(1-tor)+rhog*grp*tor
     rhoy(1) = rhog*gxdgyp/grp*tor
     rhox(2) = -rhox(1)
     rhoy(2) = -rhoy(1)
     rhox(3) = 0
     rhoy(3) = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor
     rhox(4) = 0
     rhoy(4) = -rhoy(3)
     !    calculate avg. e-field...
     !    do 1,2,4 point average, where lr is the no. of points...

     phip=0.
     exp1=0.
     eyp=0.
     ezp=0.
     delbxp=0.
     delbyp=0.
     dpdzp = 0.
     dadzp = 0.
     aparp = 0.

     !MVPT
     delbxhp = 0.
     delbyhp = 0.
     dahdzp = 0.
     aparhp = 0.

#ifdef OPENACC
!$acc loop seq                                                                                                                                                                         
#endif     
     !  4 pt. avg. done explicitly for vectorization...
     do l=1,lr(1)
        !
        xs=x2b(m)+rhox(l) !rwx(1,l)*rhog
        yt=y2b(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
        !
        !   particle can go out of bounds during gyroavg...
        xt=mod(xs+800.*lx,lx)
        yt=mod(yt+800.*ly,ly)
        xt = min(xt,lx-1.0e-8)
        yt = min(yt,ly-1.0e-8)

        include "ppushbngp.h"
     enddo
     exp1 = exp1/4.
     eyp = eyp/4.
     ezp = ezp/4.
     delbxp = delbxp/4.
     delbyp = delbyp/4.
     dpdzp = dpdzp/4.
     dadzp = dadzp/4.
     aparp = aparp/4.

     !MVPT
     delbxhp = delbxhp/4.
     delbyhp = delbyhp/4.
     dahdzp = dahdzp/4.
     aparhp = aparhp/4.
     
     !
     vfac = 0.5*(mbeam*u2b(m)**2 + 2.*mub(m)*b)
     vp0 = 1./b**2*lr0/q0*qhatp*fp/radiusp*grcgtp
     vp0 = vp0*vncp*vexbsw

     vpar = u2b(m)-qbeam/mbeam*aparp*nonlinb
     bstar = b*(1+mbeam*vpar/(qbeam*b)*bdcrvbp)
     enerb=(mub(m)+mbeam*vpar*vpar/b)/qbeam*b/bstar*tor

     dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
     vxdum = (eyp/b+vpar/b*delbxp)*dum1
     xdot = vxdum*nonlinb -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
     ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlinb &
          +iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
          (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0   &
          +enerb/(bfldp**2)*psipp*lr0/q0/radiusp**2*(dbdrp*grp**2+dbdtp*grdgtp) &
          -mbeam*vpar**2/(qbeam*bstar*b)*(psip2p*grp**2/radiusp+curvbzp)*lr0/(radiusp*q0) &
          -dipdrp/radiusp*mbeam*vpar**2/(qbeam*bstar*b)*grcgtp*lr0/q0*qhatp 

     zdot0 = vpar*b/bstar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp
     zdot1 = q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp &
          -1./b**2*q0*br0*fp/radiusp*grcgtp*vncp*vexbsw/jfnp &
          -dipdrp/radiusp*mbeam*vpar**2/(qbeam*bstar*b)*q0*br0*grcgtp/jfnp
     zdot = zdot0+zdot1
     
     pzd0 = tor*(-mub(m)/mbeam/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar &
          +mub(m)*vpar/(qbeam*bstar*b)*dipdrp/radiusp*dbdtp*grcgtp
     pzdot = pzd0 + (qbeam/mbeam*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
          +qbeam/mbeam*(-xdot*delbyp+ydot*delbxp+zdot*dadzp))*nonlinb    &
          +qbeam/mbeam*(delbyp*eyp+delbxp*exp1)* 1./bstar*lr0/q0*qhatp*fp/radiusp*grcgtp *nonlinb  &
          +vpar/bstar * ( exp1 *(-1)/b**2*fp/radiusp*dbdtp*grcgtp                                     &
                         +eyp *(-dydrp*fp/(radiusp*b**2)*dbdtp*grcgtp+fp/(radiusp*b**2)*dbdrp*lr0/q0*qhatp*grcgtp) &
                         +eyp /b*(psip2p*grp**2/radiusp+curvbzp)*(-1)*lr0/(radiusp*q0) &
                         -eyp *dipdrp/(radiusp*b)*lr0/q0*qhatp*grcgtp) *nonlinb        &
          -mub(m)/(mbeam*bstar)*(-dbdtp*(-delbyp+delbxp*dydrp)+dbdrp*delbxp*lr0/q0*qhatp)*fp/radiusp*grcgtp               
     
     x3b(m) = x2b(m) + 0.5*dt*xdot
     y3b(m) = y2b(m) + 0.5*dt*ydot
     z3b(m) = z2b(m) + 0.5*dt*zdot
     u3b(m) = u2b(m) + 0.5*dt*pzdot


     if(itube==1)go to 333
     if(abs(pzp-pzib(m))>pzcrit(1).or.abs(vfac-ekib(m))>0.5*ekib(m))then
#ifdef OPENACC
        !$acc atomic                                                                                                                                                                   
#endif
        mynopi = mynopi+1
        x3b(m) = xiib(m)
        z3b(m) = z0ib(m)
        r = x3b(m)-lx/2+lr0
        k = int(z3b(m)/delz)
        wz0 = ((k+1)*delz-z3b(m))/delz
        wz1 = 1-wz0
        th = wz0*thfnz(k)+wz1*thfnz(k+1)

        i = int((r-rin)/dr)
        wx0 = (rin+(i+1)*dr-r)/dr
        wx1 = 1.-wx0
        k = int((th+pi)/dth)
        wz0 = (-pi+(k+1)*dth-th)/dth
        wz1 = 1.-wz0
        b = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
             +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
        u3b(m) = u0ib(m)     
        u2b(m) = u3b(m)
        x2b(m) = x3b(m)
        z2b(m) = z3b(m)
     end if

333  continue
     laps=anint((z3b(m)/lz)-.5)*(1-peritr)
     r=x3b(m)-0.5*lx+lr0
     i = int((r-rin)/dr)
     i = min(i,nr-1)
     i = max(i,0)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     qr = wx0*sf(i)+wx1*sf(i+1)
     y3b(m)=mod(y3b(m)-laps*2*pi*qr*lr0/q0*sign(1.0,q0)+8000.*ly,ly)
     if(x3b(m)>lx.and.iperidf==0)then
        x3b(m) = lx-1.e-8
        z3b(m)=lz-z3b(m)
        x2b(m) = x3b(m)
        z2b(m) = z3b(m)
        wb(m) = 0.0
     end if
     if(x3b(m)<0..and.iperidf==0)then
        x3b(m) = 1.e-8
        z3b(m)=lz-z3b(m)
        x2b(m) = x3b(m)
        z2b(m) = z3b(m)
!        wb(m) = 0.0        
     end if
     z3b(m)=mod(z3b(m)+8.*lz,lz)
     x3b(m)=mod(x3b(m)+800.*lx,lx)         
     x3b(m) = min(x3b(m),lx-1.0e-8)
     y3b(m) = min(y3b(m),ly-1.0e-8)
     z3b(m) = min(z3b(m),lz-1.0e-8)

  enddo
#ifdef OPENMP
  !$omp end parallel do
#endif
#ifdef OPENACC
!$acc end parallel                                                                                                                                                                     
#endif  
  call MPI_ALLREDUCE(mynopi,nopb,1,MPI_integer, &
       MPI_SUM, MPI_COMM_WORLD,ierr)

  np_old=mmb
  initpmove_start_tm=initpmove_start_tm+MPI_WTIME()
  call test_init_pmove(z3b(:),np_old,lz,ierr)
  initpmove_end_tm=initpmove_end_tm+MPI_WTIME()

  pmove_start_tm=pmove_start_tm+MPI_WTIME()    
  call test_pmove(x2b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(x3b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(y2b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(y3b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(z2b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(z3b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(u2b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(u3b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(mub(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(wb(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call test_pmove(xiib(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(z0ib(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(pzib(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(ekib(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(u0ib(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  pmove_end_tm=pmove_end_tm+MPI_WTIME()  
  call end_pmove(ierr)
  mmb=np_new

  !      return
end subroutine ppushb

!-----------------------------------------------------------------------

subroutine cpushb(n)

  use gem_com
  use gem_equil
  implicit none
  INTEGER :: n
  real :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
  real :: wx0,wx1,wy0,wy1,wz0,wz1,w3old
  INTEGER :: m,i,j,k,ns,l,mynowi
  INTEGER :: np_old,np_new
  real :: vfac,vpar,vxdum,dum,xdot,ydot,zdot,pzdot,edot,pzd0,vp0
  real :: rhog,xt,yt,zt,kap,xs,pidum,dum1,kaptp,kapnp,xnp
  real :: b,th,r,enerb,cost,sint,qr,laps,sz,ter,bstar,dtp
  real :: myke,mypfl_es(1:nsubd),mypfl_em(1:nsubd),myavewi
  real :: myefl_es(1:nsubd),myefl_em(1:nsubd),mynos
  real :: ketemp,pfltemp
  real :: efltemp,nostemp
  real :: sbuf(10),rbuf(10)
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,grdgtp
  real :: grp,gxdgyp,rhox(4),rhoy(4),psp,pzp,vncp,vparspp,psip2p,bdcrvbp,curvbzp,dipdrp
  real :: fdum,gdum,fisrcp,dnisrcp,avwixepsp,fovg,avwixezp,dnisrczp
  !MVPT
  real :: dapsdtp, delbxhp, delbyhp, dahdzp, aparhp,zdot0,zdot1

  
  sbuf(1:10) = 0.
  rbuf(1:10) = 0.
  myavewi = 0.
  myke=0.    
  mypfl_es=0.    
  mypfl_em=0.    
  myefl_es=0. 
  myefl_em=0. 
  mynos=0.   
  ketemp=0.
  pfltemp=0.
  efltemp=0.
  nostemp=0.
  mynowi = 0

#ifdef OPENMP
!$omp parallel do default(shared) &
!$omp private(i,j,k,l,r,th,wx0,wx1,wy0,wy1,wz0,wz1,dbdrp,dbdtp,grcgtp,bfldp,radiusp,dydrp,qhatp,grp,gxdgyp,curvbzp,bdcrvbp,grdgtp) &
!$omp private(fp,jfnp,psipp,psp,ter,kaptp,kapnp,xnp,b,psip2p,dipdrp,pzp,vncp,vparspp,rhog,rhox,rhoy) &
!$omp private(phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp,xs,xt,yt,vfac,vp0,vpar,bstar,enerb,xdot,ydot,zdot) &
!$omp private(pzd0,pzdot,edot,dum,dum1,vxdum,kap,w3old,laps,qr) &
!$omp private(fdum,gdum,fisrcp,dnisrcp,avwixepsp,fovg,dtp,avwixezp,dnisrczp) &
!$omp private(dapsdtp,delbxhp,delbyhp,dahdzp,aparhp,zdot0,zdot1) &    !MVPT
!$omp reduction(+: myke,mypfl_es,mypfl_em,myefl_es,myefl_em,mynos)
#endif

#ifdef OPENACC
!$acc parallel                                                                                                                                                                         
!$acc loop gang vector private(rhox,rhoy)                                                                                                                                              
#endif

  do m=1,mmb
     r=x3b(m)-0.5*lx+lr0

     k = int(z3b(m)/delz)
     wz0 = ((k+1)*delz-z3b(m))/delz
     wz1 = 1-wz0
     th = wz0*thfnz(k)+wz1*thfnz(k+1)

     i = int((r-rin)/dr)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     k = int((th+pi)/dth)
     wz0 = (-pi+(k+1)*dth-th)/dth
     wz1 = 1.-wz0
     dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
          +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
     dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
          +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1) 
     grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
          +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 
     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
     radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
          +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
     dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
          +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1) 
     qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
          +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1) 
     grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
          +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1) 
     gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
          +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1) 

     curvbzp = wx0*wz0*curvbz(i,k)+wx0*wz1*curvbz(i,k+1) &
          +wx1*wz0*curvbz(i+1,k)+wx1*wz1*curvbz(i+1,k+1) 
     bdcrvbp = wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
          +wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1) 
     grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
          +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1) 

     fp = wx0*f(i)+wx1*f(i+1)        
     jfnp = wz0*jfn(k)+wz1*jfn(k+1)
     psipp = wx0*psip(i)+wx1*psip(i+1)        
     psp = wx0*psi(i)+wx1*psi(i+1)        
     b=1.-tor+tor*bfldp
     psip2p = wx0*psip2(i)+wx1*psip2(i+1)        
     dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)        
     pzp = mbeam*u3b(m)/b*fp/br0-qbeam*psp/br0
     vncp = 0.

     rhog=sqrt(2.*b*mub(m)*mbeam)/(qbeam*b)*iflr

     rhox(1) = rhog*(1-tor)+rhog*grp*tor
     rhoy(1) = rhog*gxdgyp/grp*tor
     rhox(2) = -rhox(1)
     rhoy(2) = -rhoy(1)
     rhox(3) = 0
     rhoy(3) = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor
     rhox(4) = 0
     rhoy(4) = -rhoy(3)
     !    calculate avg. e-field...
     !    do 1,2,4 point average, where lr is the no. of points...

     phip=0.
     exp1=0.
     eyp=0.
     ezp=0.
     delbxp = 0.
     delbyp = 0.
     dpdzp = 0.
     dadzp = 0.
     aparp = 0.

     !MVPT
     delbxhp = 0.
     delbyhp = 0.
     dahdzp = 0.
     aparhp = 0.

#ifdef OPENACC
!$acc loop seq                                                                                                                                                                         
#endif

     !  4 pt. avg. written out explicitly for vectorization...
     do l=1,lr(1)
        xs=x3b(m)+rhox(l) !rwx(1,l)*rhog
        yt=y3b(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
        !   BOUNDARY
        xt=mod(xs+800.*lx,lx)
        yt=mod(yt+800.*ly,ly)
        xt = min(xt,lx-1.0e-8)
        yt = min(yt,ly-1.0e-8)

        include "cpushbngp.h"
     enddo

     exp1=exp1/4.
     eyp=eyp/4.
     ezp=ezp/4.
     delbxp=delbxp/4.
     delbyp=delbyp/4.
     dpdzp = dpdzp/4.
     dadzp = dadzp/4.
     aparp = aparp/4.

     !MVPT
     delbxhp = delbxhp/4.
     delbyhp = delbyhp/4.
     dahdzp = dahdzp/4.
     aparhp = aparhp/4.


     vfac = 0.5*(mbeam*u3b(m)**2 + 2.*mub(m)*b)
     vp0 = 1./b**2*lr0/q0*qhatp*fp/radiusp*grcgtp
     vp0 = vp0*vncp*vexbsw

     vpar = u3b(m)-qbeam/mbeam*aparp *nonlinb
     bstar = b*(1+mbeam*vpar/(qbeam*b)*bdcrvbp)
     enerb=(mub(m)+mbeam*vpar*vpar/b)/qbeam*b/bstar*tor

     dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
     vxdum = (eyp/b+vpar/b*delbxp)*dum1
     xdot = vxdum*nonlinb -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
     ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlinb     &
          +iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
          (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0 &
          +enerb/(bfldp**2)*psipp*lr0/q0/radiusp**2*(dbdrp*grp**2+dbdtp*grdgtp) &
          -mbeam*vpar**2/(qbeam*bstar*b)*(psip2p*grp**2/radiusp+curvbzp)*lr0/(radiusp*q0) &
          -dipdrp/radiusp*mbeam*vpar**2/(qbeam*bstar*b)*grcgtp*lr0/q0*qhatp  

     !MVPT
     zdot0 = vpar*b/bstar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp
     zdot1 = q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp &
          -1./b**2*q0*br0*fp/radiusp*grcgtp*vncp*vexbsw/jfnp &
          -dipdrp/radiusp*mbeam*vpar**2/(qbeam*bstar*b)*q0*br0*grcgtp/jfnp
     zdot = zdot0+zdot1

     pzd0 = tor*(-mub(m)/mbeam/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar &
          +mub(m)*vpar/(qbeam*bstar*b)*dipdrp/radiusp*dbdtp*grcgtp
     pzdot = pzd0 + (qbeam/mbeam*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
          +qbeam/mbeam*(-xdot*delbyp+ydot*delbxp+zdot*dadzp))*nonlinb    &
          +qbeam/mbeam*(delbyp*eyp+delbxp*exp1)* 1./bstar*lr0/q0*qhatp*fp/radiusp*grcgtp *nonlinb  &
          +vpar/bstar * ( exp1 *(-1)/b**2*fp/radiusp*dbdtp*grcgtp                                     &
                         +eyp *(-dydrp*fp/(radiusp*b**2)*dbdtp*grcgtp+fp/(radiusp*b**2)*dbdrp*lr0/q0*qhatp*grcgtp) &
                         +eyp /b*(psip2p*grp**2/radiusp+curvbzp)*(-1)*lr0/(radiusp*q0) &
                         -eyp *dipdrp/(radiusp*b)*lr0/q0*qhatp*grcgtp) *nonlinb        &
          -mub(m)/(mbeam*bstar)*(-dbdtp*(-delbyp+delbxp*dydrp)+dbdrp*delbxp*lr0/q0*qhatp)*fp/radiusp*grcgtp
     
     x3b(m) = x2b(m) + dt*xdot
     y3b(m) = y2b(m) + dt*ydot
     z3b(m) = z2b(m) + dt*zdot
     u3b(m) = u2b(m) + dt*pzdot

     laps=anint((z3b(m)/lz)-.5)*(1-peritr)
     r=x3b(m)-0.5*lx+lr0
     i = int((r-rin)/dr)
     i = min(i,nr-1)
     i = max(i,0)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     qr = wx0*sf(i)+wx1*sf(i+1)
     y3b(m)=mod(y3b(m)-laps*2*pi*qr*lr0/q0*sign(1.0,q0)+8000.*ly,ly)
     if(x3b(m)>lx.and.iperidf==0)then
        x3b(m) = lx-1.e-8
        z3b(m)=lz-z3b(m)
        x2b(m) = x3b(m)
        z2b(m) = z3b(m)
        wb(m) = 0.
     end if
     if(x3b(m)<0..and.iperidf==0)then
        x3b(m) = 1.e-8
        z3b(m)=lz-z3b(m)
        x2b(m) = x3b(m)
        z2b(m) = z3b(m)
!        wb(m) = 0.0
     end if
     z3b(m)=mod(z3b(m)+8.*lz,lz)
     x3b(m)=mod(x3b(m)+800.*lx,lx)         
     x3b(m) = min(x3b(m),lx-1.0e-8)
     y3b(m) = min(y3b(m),ly-1.0e-8)
     z3b(m) = min(z3b(m),lz-1.0e-8)

     !     particle diagnostics done here because info is available...
     k = int(x3b(m)/(lx/nsubd))
     k = min(k,nsubd-1)
     k = k+1
#ifdef OPENACC
     !$acc atomic
     mypfl_es(k)=mypfl_es(k) + (eyp/b)*dum1/grp *wb(m)
     !$acc atomic
     mypfl_em(k)=mypfl_em(k) + (vpar*delbxp/b)*dum1/grp *wb(m)
     !$acc atomic
     myefl_es(k)=myefl_es(k) + vfac*(eyp/b)*dum1/grp *wb(m)
     !$acc atomic
     myefl_em(k)=myefl_em(k) + vfac*(vpar*delbxp/b)*dum1/grp *wb(m)
     !$acc atomic
     myke=myke + vfac *wb(m)
     !$acc atomic
     mynos=mynos + wb(m)
#endif
     
#ifdef OPENMP     
     mypfl_es(k)=mypfl_es(k) + (eyp/b)*dum1/grp *wb(m)
     mypfl_em(k)=mypfl_em(k) + (vpar*delbxp/b)*dum1/grp *wb(m)
     myefl_es(k)=myefl_es(k) + vfac*(eyp/b)*dum1/grp *wb(m)
     myefl_em(k)=myefl_em(k) + vfac*(vpar*delbxp/b)*dum1/grp *wb(m)
     myke=myke + vfac *wb(m)
     mynos=mynos + wb(m)
#endif
     
     !     xn+1 becomes xn...
     u2b(m)=u3b(m)
     x2b(m)=x3b(m)
     y2b(m)=y3b(m)
     z2b(m)=z3b(m)

     !     100     continue
  enddo
#ifdef OPENMP
!$omp end parallel do                                                                                                                                                                  
#endif

#ifdef OPENACC
!$acc end parallel                                                                                                                                                                     
#endif

  sbuf(1)=myke
  sbuf(2)=myefl_es(nsubd/2)
  sbuf(3)=mypfl_es(nsubd/2)
  sbuf(4)=mynos
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,           &
       MPI_COMM_WORLD,ierr)

  ketemp=rbuf(1)
  efltemp=rbuf(2)
  pfltemp=rbuf(3)
  nostemp=rbuf(4)
  avewb(n)=nostemp/( tmmb )
  kebeam(n)=ketemp/( tmmb )

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  sbuf(1:nsubd) = myefl_es(1:nsubd)
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,  &
       MPI_COMM_WORLD,ierr)
  do k = 1,nsubd
     eflb_es(k,n)=rbuf(k)/( float(tmmb) )*totvol/vol(k)*cn0b
  end do

  sbuf(1:nsubd) = myefl_em(1:nsubd)
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,  &
       MPI_COMM_WORLD,ierr)
  do k = 1,nsubd
     eflb_em(k,n)=rbuf(k)/( float(tmmb) )*totvol/vol(k)*cn0b
  end do

  sbuf(1:nsubd) = mypfl_es(1:nsubd)
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,  &
       MPI_COMM_WORLD,ierr)
  do k = 1,nsubd
     pflb_es(k,n)=rbuf(k)/( float(tmmb) )*totvol/vol(k)*cn0b
  end do

  sbuf(1:nsubd) = mypfl_em(1:nsubd)
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,  &
       MPI_COMM_WORLD,ierr)
  do k = 1,nsubd
     pflb_em(k,n)=rbuf(k)/( float(tmmb) )*totvol/vol(k)*cn0b
  end do

  !      pfl(1,n)=pfltemp/( float(tmm(1)) )
  !      efl(1,n)=mims(ns)/tets(1)*efltemp/( float(tmm(1)) )

  np_old=mmb
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  initpmove_start_tm=initpmove_start_tm+MPI_WTIME()  
  call test_init_pmove(z3b(:),np_old,lz,ierr)
  initpmove_end_tm=initpmove_end_tm+MPI_WTIME()
  
  pmove_start_tm=pmove_start_tm+MPI_WTIME()    
  call test_pmove(x2b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(x3b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(y2b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(y3b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(z2b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(z3b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(u2b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(u3b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(mub(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(wb(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(xiib(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(z0ib(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(pzib(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(ekib(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(u0ib(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  pmove_end_tm=pmove_end_tm+MPI_WTIME()  
  call end_pmove(ierr)
  mmb=np_new
  !     write(*,*)MyId,mm(ns)

  !      return
end subroutine cpushb

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine grid1b(ip,n)

  !    source quantities are are calculated: n_i
  !    right now only ion quantitities are calculated...

  use gem_com
  use gem_equil
  use omp_lib
  implicit none
  real :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
  real :: enerb,vxdum,dum,xdot,ydot
  INTEGER :: m,n,i,j,k,l,ns,ip
  real :: wx0,wx1,wy0,wy1,wz0,wz1,vte
  real :: sz,wght,r,th,cost,sint,b,qr,dv
  real :: xt,yt,rhog,pidum,vpar,xs,dely,vfac
  real :: lbfs(0:imx,0:jmx)
  real :: rbfs(0:imx,0:jmx)
  real :: lbfr(0:imx,0:jmx)
  real :: rbfr(0:imx,0:jmx)
  real :: myden(0:imx,0:jmx,0:1),myjpar(0:imx,0:jmx,0:1)
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
  real :: grp,gxdgyp,rhox(4),rhoy(4)

  if(idg.eq.1)write(*,*)'enter ion grid1',mmb

     denb(:,:,:)=0.
     jparb(:,:,:)=0.
     myden = 0.
     myjpar = 0.

#ifdef OPENMP     
!$omp parallel do &
!$omp default(none) & 
!$omp private(dv,i,j,k,l,r,wx0,wx1,wz0,wz1,th,dbdrp,dbdtp,grcgtp,bfldp,radiusp,dydrp,qhatp,grp,gxdgyp,fp,jfnp,psipp,b) &
!$omp private(rhog,rhox,rhoy,vfac,wght,vpar,xs,xt,yt,aparp) &
!$omp shared(mmb,rin,u3b,x3b,y3b,z3b,mub,wb,thfnz,dbdr,dbdth,grcgt,bfld,radius,dydr,qhat,gr,gxdgy,f,jfn,psip) &
!$omp shared(mbeam,qbeam,tor,iflr,emass,amie,qel,lr0,q0,lx,ly,lz,dx,dy,dz,delz,dr,dth,lr,pi,gclr,tclr,kcnt,apar,nonlinb) &
!$omp reduction(+: myden,myjpar)
#endif

#ifdef OPENACC
     !$acc data copy(myden,myjpar)
     !$acc parallel
     !$acc loop gang vector private(rhox,rhoy)
#endif
     do m=1,mmb
        dv=float(lr(1))*(dx*dy*dz)
        r=x3b(m)-0.5*lx+lr0

        k = int(z3b(m)/delz)
        wz0 = ((k+1)*delz-z3b(m))/delz
        wz1 = 1-wz0
        th = wz0*thfnz(k)+wz1*thfnz(k+1)

        i = int((r-rin)/dr)
        wx0 = (rin+(i+1)*dr-r)/dr
        wx1 = 1.-wx0
        k = int((th+pi)/dth)
        wz0 = (-pi+(k+1)*dth-th)/dth
        wz1 = 1.-wz0
        dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
             +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
        dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
             +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1) 
        grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
             +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 
        bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
             +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
        radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
             +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
        dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
             +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1) 
        qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
             +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1) 
        grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
             +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1) 
        gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
             +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1) 
        fp = wx0*f(i)+wx1*f(i+1)        
        jfnp = wz0*jfn(k)+wz1*jfn(k+1)
        psipp = wx0*psip(i)+wx1*psip(i+1)        
        !         b=1.-lr0/br0*cost
        b=1.-tor+tor*bfldp

        rhog=sqrt(2.*b*mub(m)*mbeam)/(qbeam*b)*iflr

        rhox(1) = rhog*(1-tor)+rhog*grp*tor
        rhoy(1) = rhog*gxdgyp/grp*tor
        rhox(2) = -rhox(1)
        rhoy(2) = -rhoy(1)
        rhox(3) = 0
        rhoy(3) = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor
        rhox(4) = 0
        rhoy(4) = -rhoy(3)

        vfac=0.5*(mbeam*u3b(m)**2 + 2.*mub(m)*b )
        wght=wb(m)/dv

        aparp = 0.

#ifdef OPENACC
        !$acc loop seq
#endif        
        !  4 pt. avg. done explicitly for vectorization...
        do l=1,lr(1)
        !
           xs=x3b(m)+rhox(l) !rwx(1,l)*rhog
           yt=y3b(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
           !
           !   particle can go out of bounds during gyroavg...
           xt=mod(xs+800.*lx,lx)
           yt=mod(yt+800.*ly,ly)
           xt = min(xt,lx-1.0e-8)
           yt = min(yt,ly-1.0e-8)
           i=int(xt/dx+0.5)
           j=int(yt/dy+0.5)
           k=int(z3b(m)/dz+0.5)-gclr*kcnt
           
           aparp = aparp+apar(i,j,k)           
        enddo
        aparp = aparp/4.

        vpar = u3b(m) -qbeam/mbeam*aparp*nonlinb

#ifdef OPENACC
        !$acc loop seq
#endif        
        
        !    now do 1,2,4 point average, where lr is the no. of points...
        do l=1,lr(1)
           xs=x3b(m)+rhox(l) !rwx(1,l)*rhog
           yt=y3b(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
           xt=max(1e-8,xt) !mod(xs+800.*lx,lx)
           yt=mod(yt+800.*ly,ly)
           xt = min(xt,lx-1.e-8)
           yt = min(yt,ly-1.e-8)

           i=int(xt/dx+0.5)
           j=int(yt/dy+0.5)
           k=int(z3b(m)/dz+0.5)-gclr*kcnt

#ifdef OPENACC
           !$acc atomic update
           myden(i,j,k) = myden(i,j,k) + wght
           !$acc atomic update           
           myjpar(i,j,k) = myjpar(i,j,k)+wght*vpar
#endif           
#ifdef OPENMP
           myden(i,j,k) = myden(i,j,k) + wght
           myjpar(i,j,k) = myjpar(i,j,k)+wght*vpar
#endif           
        enddo
     enddo
#ifdef OPENMP
     !$omp end parallel do
#endif

#ifdef OPENACC
     !$acc end parallel
     !$acc end data
#endif

     if(idg.eq.1)write(*,*)myid,'pass ion grid1'
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     !   enforce periodicity
     call enforce(myden(:,:,:))
     call enforce(myjpar)

     do i=0,im
        do j=0,jm
           do k=0,mykm
              denb(i,j,k)=qbeam*myden(i,j,k)/n0b/jac(i,k)*cn0b
              jparb(i,j,k) = qbeam*myjpar(i,j,k)/n0b/jac(i,k)*cn0b
           enddo
        enddo
     enddo
     denb0 = denb
!     call ngt0(denb)
!     call ngt0(jparb)
     if(ibmst<=n .and. n<ibmend)then
        denbeq = denbeq+denb/float(ibmend-ibmst)
        jparbeq = jparbeq+jparb/float(ibmend-ibmst)
     end if
     
     do i=0,im
        do j=0,jm
           do k=0,mykm
              rho(i,j,k)=rho(i,j,k)+(denb(i,j,k)-denbeq(i,j,k)*ibmcnl)*isbeam
              jion(i,j,k) = jion(i,j,k)+(jparb(i,j,k)-jparbeq(i,j,k)*ibmcnl)*isbeam
           enddo
        enddo
     enddo


  !      return
end subroutine grid1b


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine loadb

  use gem_com
  use gem_equil
  implicit none
  INTEGER :: i,k,m,idum,ns,m1
  INTEGER :: np_old,np_new
  integer*8 :: j,maxtry
  real :: vpar,vperp2,r,qr,th,b,cost,ter,vbeam,ptch
  real :: avgv,myavgv,avgw,myavgw
  real :: dumx,dumy,dumz,jacp,rdum,drdum,jv,kappan,wn,s0
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,psp
  real :: grp,gxdgyp,zoldp
  real :: wx0,wx1,wz0,wz1
  real :: denprof(0:imx)
  
  maxtry=1000**3
  !radial density profile, flat in [rina,dumx], linear in [dumx,dumy]
  vbeam = 20.0
  
  dumx = 0.4
  dumy = 0.7
  rdum = 0.05
  drdum = 0.25
  s0=0.45
  kappan=10.0
  wn=0.2
  denprof = 1.0
  do i = 0,imx
     r = rina+lxa/imx*i
!     if(r>=dumx .and. r<=dumy)denprof(i)=1.0-(r-dumx)/(dumy-dumx)
!     if(r>dumy)denprof(i)=0.
!     denprof(i)=exp(-(r-rdum)**2/drdum**2)
     denprof(i) = exp(-kappan*wn*a/rmaj0*tanh((r-s0)/(wn)))
  end do
  myavgv=0.
  avgv=0.
  avgw = 0.
  myavgw = 0.

  !      ns = 1
  m = 0
  do j = 1,maxtry

     !     load a slab of ions...

     dumx=lx*ran2(iseed) !revers(MyId*cnt+j,2) !ran2(iseed)
     dumy=ly*ran2(iseed) !revers(MyId*cnt+j,3) !ran2(iseed)
     dumz=lz*ran2(iseed) !revers(MyId*cnt+j,5) !ran2(iseed)
     dumz = min(dumz,lz-1.e-8)
     r = lr0+dumx-0.5*lx
     th = (dumz-lz/2)/(q0*br0)
     i = int((r-rin)/dr)
     k = int((pi+th)/dth)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     wz0 = (-pi+(k+1)*dth-th)/dth
     wz1 = 1.-wz0
     jacp = wx0*wz0*jacob(i,k)+wx0*wz1*jacob(i,k+1) &
          +wx1*wz0*jacob(i+1,k)+wx1*wz1*jacob(i+1,k+1) 
     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
     b=1.-tor+tor*bfldp
     jv = 1.0 !b if(mu,vpar), 1 if(v,ptch)
     if(ran2(iseed)<0.8*jv*(jacp/jacmax)*denprof(int((r-rin)/dx))/denprof(1) )then
        m = m+1
        if(m>mmb)goto 170
        x2b(m)=min(dumx,lx-1e-8)
        y2b(m)=min(dumy,ly-1e-8)

        k = int((th+pi)/dth)
        wz0 = (-pi+(k+1)*dth-th)/dth
        wz1 = 1-wz0
        z2b(m) = wz0*zfnth(k)+wz1*zfnth(k+1)
        z2b(m)=min(z2b(m),lz-1e-8)

        r=x2b(m)-0.5*lx+lr0
        i = int((r-rin)/dr)
        wx0 = (rin+(i+1)*dr-r)/dr
        wx1 = 1.-wx0
        k = int((th+pi)/dth)
        wz0 = (-pi+(k+1)*dth-th)/dth
        wz1 = 1.-wz0
        bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
             +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
        psp = wx0*psi(i)+wx1*psi(i+1)
        fp = wx0*f(i)+wx1*f(i+1)        
        b=1.-tor+tor*bfldp

        ptch = 2*(ran2(iseed)-0.5)
        u2b(m) = vbeam*ptch
        vperp2 = vbeam**2*(1-ptch**2)
        mub(m) = 0.5*mbeam*vperp2/b
        ekib(m) = mub(m)*b+0.5*mbeam*u2b(m)**2
        pzib(m) = mbeam*u2b(m)/b*fp/br0-qbeam*psp/br0
        z0ib(m) = z2b(m)
        xiib(m) = x2b(m)
        u0ib(m) = u2b(m)
        wb(m) = 1.0
     end if
  enddo
170 continue

  if(idg.eq.1)write(*,*)'all reduce'
  do m=1,mmb
     u2b(m)=u2b(m) !-avgv
     x3b(m)=x2b(m)
     y3b(m)=y2b(m)
     z3b(m)=z2b(m)
     u3b(m)=u2b(m)
  enddo
  do i = 0,last
     if(myid==i)then
        open(9, file='plot',status='unknown',position='append')
        write(9,*)'loadbb before init_pmove myid,j,mmb',myid,j,mmb
        close(9)
     end if
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end do

  np_old=mmb
  call init_pmove(z2b(:),np_old,lz,ierr)
  if(idg.eq.1)write(*,*)'pass init_pmove'
  !     
  call pmove(x2b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(x3b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(y2b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(y3b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z2b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z3b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(u2b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(u3b(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(mub(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(wb(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(xiib(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z0ib(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(pzib(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(ekib(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(u0ib(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !     
  call end_pmove(ierr)
  mmb=np_new
  
#ifdef OPENACC
!$acc update device(mub,x2b,y2b,z2b,u2b,wb,x3b,y3b,z3b,u3b,xiib,z0ib,pzib,ekib,u0ib)
#endif
  do i = 0,last
     if(myid==i)then
        open(9, file='plot',status='unknown',position='append')
        write(9,*)'myid,j,mmb',myid,j,mmb
        close(9)
     end if
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end do
  !      return
end subroutine loadb

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine jieb(ip,n)

  !    source quantities are are calculated: n_i
  !    right now only ion quantitities are calculated...

  use gem_com
  use gem_equil
  implicit none
  real :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
  real :: enerb,vxdum,dum,xdot,ydot,zdot,xdt0,ydt0,pidum,dum1,dum2,fdum,gdum,fovg
  INTEGER :: m,n,i,j,k,l,ns,ip,nonfi,nonfe
  real :: wx0,wx1,wy0,wy1,wz0,wz1,vte
  real :: sz,wght,wght0,wght1,wght2,r,th,cost,sint,b,qr,dv,kap,ter
  real :: kapnp,kaptp,xnp,psip2p,bdcrvbp,curvbzp,dipdrp,bstar
  real :: xt,yt,rhog,vpar,xs,dely,vfac,vp0
  real :: lbfs(0:imx,0:jmx)
  real :: rbfs(0:imx,0:jmx)
  real :: lbfr(0:imx,0:jmx)
  real :: rbfr(0:imx,0:jmx)
  real :: myjpar(0:imx,0:jmx,0:1),myjpex(0:imx,0:jmx,0:1)
  real :: myjpey(0:imx,0:jmx,0:1),myupar(0:imx,0:jmx,0:1)
  real :: myupex(0:imx,0:jmx,0:1),myupey(0:imx,0:jmx,0:1)
  real :: myupazd(0:imx,0:jmx,0:1)
  real :: mydnidt(0:imx,0:jmx,0:1),mydnedt(0:imx,0:jmx,0:1)
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,grdgtp
  real :: grp,gxdgyp,rhox(4),rhoy(4),vncp,vparspp
  real :: sn(4),cn(4),bdgrhop,e1gxp,e2gxp,be1gxp,be2gxp,e1gyp,e2gyp,be1gyp,be2gyp,bdgvrx,bdgvry
  
  sn(1) = 1;   cn(1) = 0
  sn(2) = -1;  cn(2) = 0
  sn(3) = 0;   cn(3) = 1
  sn(4) = 0;   cn(4) = -1

  nonfi = 1 
  nonfe = 1 

     myjpar = 0.
     myjpex = 0.
     myjpey = 0.
#ifdef OPENMP
!$omp parallel do &
!$omp private(i,j,k,l,dv,r,th,wx0,wx1,wz0,wz1,dbdrp,dbdtp,grcgtp,bfldp,radiusp,dydrp,qhatp,grp,gxdgyp,curvbzp,bdcrvbp,grdgtp,fp,jfnp) &
!$omp private(psipp,ter,kaptp,kapnp,xnp,vncp,vparspp,psip2p,dipdrp,b,rhog,vfac,vp0,vpar,kap,wght,wght0,wght1,rhox,rhoy) &
!$omp private(exp1,eyp,delbxp,delbyp,aparp,xs,xt,yt,bstar,enerb,dum1,vxdum,xdot,ydot,zdot,fdum,gdum,fovg) &
!$omp private(bdgrhop,e1gxp,e2gxp,be1gxp,be2gxp,e1gyp,e2gyp,be1gyp,be2gyp,bdgvrx,bdgvry) &
!$omp reduction(+: myjpar,myjpex,myjpey)
#endif

#ifdef OPENACC
     !$acc data copy(myjpar,myjpex,myjpey)
     !$acc parallel
     !$acc loop gang vector private(rhox,rhoy)
#endif
     
     do m=1,mmb
        dv=float(lr(1))*(dx*dy*dz)

        r=x3b(m)-0.5*lx+lr0

        k = int(z3b(m)/delz)
        wz0 = ((k+1)*delz-z3b(m))/delz
        wz1 = 1-wz0
        th = wz0*thfnz(k)+wz1*thfnz(k+1)

        i = int((r-rin)/dr)
        wx0 = (rin+(i+1)*dr-r)/dr
        wx1 = 1.-wx0
        k = int((th+pi)/dth)
        wz0 = (-pi+(k+1)*dth-th)/dth
        wz1 = 1.-wz0
        dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
             +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
        dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
             +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1) 
        grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
             +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 
        bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
             +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
        radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
             +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
        dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
             +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1) 
        qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
             +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1) 
        grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
             +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1) 
        gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
             +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1) 

        curvbzp = wx0*wz0*curvbz(i,k)+wx0*wz1*curvbz(i,k+1) &
             +wx1*wz0*curvbz(i+1,k)+wx1*wz1*curvbz(i+1,k+1) 
        bdcrvbp = wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
             +wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1) 
        grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
             +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1) 

        e1gxp = wx0*wz0*e1gx(i,k)+wx0*wz1*e1gx(i,k+1) &
             +wx1*wz0*e1gx(i+1,k)+wx1*wz1*e1gx(i+1,k+1) 
        e1gyp = wx0*wz0*e1gy(i,k)+wx0*wz1*e1gy(i,k+1) &
             +wx1*wz0*e1gy(i+1,k)+wx1*wz1*e1gy(i+1,k+1) 
        e2gxp = wx0*wz0*e2gx(i,k)+wx0*wz1*e2gx(i,k+1) &
             +wx1*wz0*e2gx(i+1,k)+wx1*wz1*e2gx(i+1,k+1) 
        e2gyp = wx0*wz0*e2gy(i,k)+wx0*wz1*e2gy(i,k+1) &
             +wx1*wz0*e2gy(i+1,k)+wx1*wz1*e2gy(i+1,k+1) 
        be1gxp = wx0*wz0*bdge1gx(i,k)+wx0*wz1*bdge1gx(i,k+1) &
             +wx1*wz0*bdge1gx(i+1,k)+wx1*wz1*bdge1gx(i+1,k+1) 
        be1gyp = wx0*wz0*bdge1gy(i,k)+wx0*wz1*bdge1gy(i,k+1) &
             +wx1*wz0*bdge1gy(i+1,k)+wx1*wz1*bdge1gy(i+1,k+1) 
        be2gxp = wx0*wz0*bdge2gx(i,k)+wx0*wz1*bdge2gx(i,k+1) &
             +wx1*wz0*bdge2gx(i+1,k)+wx1*wz1*bdge2gx(i+1,k+1) 
        be2gyp = wx0*wz0*bdge2gy(i,k)+wx0*wz1*bdge2gy(i,k+1) &
             +wx1*wz0*bdge2gy(i+1,k)+wx1*wz1*bdge2gy(i+1,k+1) 

        
        fp = wx0*f(i)+wx1*f(i+1)        
        jfnp = wz0*jfn(k)+wz1*jfn(k+1)
        psipp = wx0*psip(i)+wx1*psip(i+1)        

        vncp = 0.
        psip2p = wx0*psip2(i)+wx1*psip2(i+1)        
        dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)        

        !         b=1.-lr0/br0*cost
        b=1.-tor+tor*bfldp

        rhog=sqrt(2.*b*mub(m)*mbeam)/(qbeam*b)*iflr
        vfac = 0.5*(mbeam*u3b(m)**2 + 2.*mub(m)*b)

        vp0 = 1./b**2*lr0/q0*qhatp*fp/radiusp*grcgtp
        vp0 = vp0*vncp*vexbsw
        
        wght=wb(m)/dv

        rhox(1) = rhog*(1-tor)+rhog*grp*tor
        rhoy(1) = rhog*gxdgyp/grp*tor
        rhox(2) = -rhox(1)
        rhoy(2) = -rhoy(1)
        rhox(3) = 0
        rhoy(3) = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor
        rhox(4) = 0
        rhoy(4) = -rhoy(3)

        
        exp1=0.
        eyp=0.
        delbxp = 0.
        delbyp = 0.
        aparp = 0.
#ifdef OPENACC
        !$acc loop seq
#endif        
        do l=1,lr(1)
           xs=x3b(m)+rhox(l) !rwx(1,l)*rhog
           yt=y3b(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
           xt=mod(xs+800.*lx,lx)
           yt=mod(yt+800.*ly,ly)
           xt = min(xt,lx-1.0e-8)
           yt = min(yt,ly-1.0e-8)

           i=int(xt/dx+0.5)
           j=int(yt/dy+0.5)
           k=int(z3b(m)/dz+0.5)-gclr*kcnt
           exp1=exp1 + ex(i,j,k)
           eyp=eyp + ey(i,j,k)
           delbxp = delbxp+delbx(i,j,k)
           delbyp = delbyp+delby(i,j,k)
           aparp = aparp+apar(i,j,k)
        enddo

        exp1=exp1/4.
        eyp=eyp/4.
        delbxp=delbxp/4.
        delbyp=delbyp/4.
        aparp = aparp/4.

        vpar = u3b(m)-qbeam/mbeam*aparp*nonlinb
        bstar = b*(1+mbeam*vpar/(qbeam*b)*bdcrvbp)
        enerb=(mub(m)+mbeam*vpar*vpar/b)/qbeam*b/bstar*tor
        dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
        vxdum = (eyp/b+vpar/b*delbxp)*dum1

        xdot = vxdum*nonlinb  &
             -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
        ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlinb  &
             +iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
             (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0 &
             +enerb/(bfldp**2)*psipp*lr0/q0/radiusp**2*(dbdrp*grp**2+dbdtp*grdgtp) &
             -mbeam*vpar**2/(qbeam*bstar*b)*(psip2p*grp**2/radiusp+curvbzp)*lr0/(radiusp*q0) &
             -dipdrp/radiusp*mbeam*vpar**2/(qbeam*bstar*b)*grcgtp*lr0/q0*qhatp  
        zdot =  vpar*b/bstar*(1.-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
             +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp &
             -1./b**2*q0*br0*fp/radiusp*grcgtp*vncp*vexbsw/jfnp &
             -dipdrp/radiusp*mbeam*vpar**2/(qbeam*bstar*b)*q0*br0*grcgtp/jfnp

        bdgrhop = -rhog*psipp/(2*radiusp*b*b)*dbdtp*grcgtp

#ifdef OPENACC
        !$acc loop seq
#endif        

        !    now do 1,2,4 point average, where lr is the no. of points...
        do l=1,lr(1)
           xs=x3b(m)+rhox(l) !rwx(1,l)*rhog
           yt=y3b(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
           xt=mod(xs+800.*lx,lx)
           yt=mod(yt+800.*ly,ly)
           xt = min(xt,lx-1.0e-8)
           yt = min(yt,ly-1.0e-8)

           i=int(xt/dx+0.5)
           j=int(yt/dy+0.5)
           k=int(z3b(m)/dz+0.5)-gclr*kcnt

           bdgvrx = bdgrhop*(e2gxp*cn(l)+e1gxp*sn(l))   &
                +rhog*(be2gxp*cn(l)+be1gxp*sn(l))
           bdgvry = bdgrhop*(e2gyp*cn(l)+e1gyp*sn(l))   &
                +rhog*(be2gyp*cn(l)+be1gyp*sn(l))

#ifdef OPENACC
           !$acc atomic update
           myjpar(i,j,k) = myjpar(i,j,k)+wght*zdot
           !$acc atomic update           
           myjpex(i,j,k) = myjpex(i,j,k)+(wght*xdot+wght*vpar*bdgvrx)
           !$acc atomic update           
           myjpey(i,j,k) = myjpey(i,j,k)+(wght*ydot+wght*vpar*bdgvry)
#endif           
#ifdef OPENMP 
           myjpar(i,j,k) = myjpar(i,j,k)+wght*zdot
           myjpex(i,j,k) = myjpex(i,j,k)+(wght*xdot+wght*vpar*bdgvrx)
           myjpey(i,j,k) = myjpey(i,j,k)+(wght*ydot+wght*vpar*bdgvry)
#endif           
        enddo
     enddo
#ifdef OPENACC
     !$acc end parallel
     !$acc end data
#endif

#ifdef OPENMP
     !$omp end parallel do
#endif

     !   enforce periodicity
     call enforce(myjpar)
     call enforce(myjpex)
     call enforce(myjpey)

     do i=0,im
        do j=0,jm
           do k=0,mykm
              jparb(i,j,k) = qbeam*myjpar(i,j,k)/n0b/jac(i,k)*cn0b
              jpexb(i,j,k) = qbeam*myjpex(i,j,k)/n0b/jac(i,k)*cn0b
              jpeyb(i,j,k) = qbeam*myjpey(i,j,k)/n0b/jac(i,k)*cn0b
           enddo
        enddo
     enddo
!     call ngt0(jparb)
!     call ngt0(jpexb)
!     call ngt0(jpeyb)
     if(ibmst<=n .and. n<ibmend)then
        jparbeqp = jparbeqp+jparb/float(ibmend-ibmst)
        jpexbeq = jpexbeq+jpexb/float(ibmend-ibmst)
        jpeybeq = jpeybeq+jpeyb/float(ibmend-ibmst)        
     end if
     
     do i = 0,im
        do j = 0,jm
           do k = 0,mykm
              jion(i,j,k) = jion(i,j,k)+(jparb(i,j,k)-jparbeqp(i,j,k)*ibmcnl)*isbeam
              jionx(i,j,k) = jionx(i,j,k)+(jpexb(i,j,k)-jpexbeq(i,j,k)*ibmcnl)*isbeam
              jiony(i,j,k) = jiony(i,j,k)+(jpeyb(i,j,k)-jpeybeq(i,j,k)*ibmcnl)*isbeam
              drhoidt(i,j,k) = drhoidt(i,j,k)
           end do
        end do
     end do

999 continue
  !      return
end subroutine jieb
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine jpar0b(ip,n,it,itp)

  !    source quantities are are calculated: n_i
  !    right now only ion quantitities are calculated...

  use gem_com
  use gem_equil
  implicit none
  real :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
  real :: enerb,vxdum,dum,xdot,ydot,avede0
  INTEGER :: m,n,i,j,k,l,ns,ip,it,itp
  real :: wx0,wx1,wy0,wy1,wz0,wz1,vte
  real :: sz,wght,wght0,wght1,r,th,cost,sint,b,qr,dv,xnp
  real :: xt,yt,zt,rhog,pidum,vpar,ppar,xs,dely,vfac,ter
  real :: lbfs(0:imx,0:jmx)
  real :: rbfs(0:imx,0:jmx)
  real :: lbfr(0:imx,0:jmx)
  real :: rbfr(0:imx,0:jmx)
  real :: myupa(0:imx,0:jmx,0:1),myupa0(0:imx,0:jmx,0:1),myden0(0:imx,0:jmx,0:1)
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
  real :: grp,gxdgyp
  real :: zdot,fdum,gdum,fovg

  myupa = 0.
  myupa0 = 0.
  myden0 = 0.
  !      upa0(:,:,:) = apar(ip,:,:,:)
  !      return
  if(it.eq.1)then
     upa0(:,:,:) = 0.
     upa00(:,:,:) = 0.
     den0apa(:,:,:) = 0.
     return
  end if

!$omp parallel do &
!$omp private(i,j,k,r,th,vpar,wx0,wx1,wz0,wz1,ter,xnp,bfldp,b,vfac,grcgtp,radiusp,dbdrp,jfnp,psipp,fp,enerb,zdot) &
!$omp private(dv,wght0,wght1,xt,yt,zt,aparp,ppar,fdum,gdum,fovg) &
!$omp reduction(+: myupa0,myupa,myden0)
  do m=1,mmb
     xt=x3b(m)
     yt=y3b(m)
     zt=z3b(m)
     i=int(xt/dx)
     j=int(yt/dy)
     k=0 !int(z3e(m)/dz)-gclr*kcnt
     aparp = w000(m)*apar(i,j,k)  &
          + w100(m)*apar(i+1,j,k) &
          + w010(m)*apar(i,j+1,k) &
          + w110(m)*apar(i+1,j+1,k) &
          + w001(m)*apar(i,j,k+1) &
          + w101(m)*apar(i+1,j,k+1) &
          + w011(m)*apar(i,j+1,k+1) &
          + w111(m)*apar(i+1,j+1,k+1)
     ppar = u3e(m)
     vpar = u3b(m)-qbeam/mbeam*aparp

     r=x3b(m)-0.5*lx+lr0

     k = int(z3b(m)/delz)
     wz0 = ((k+1)*delz-z3b(m))/delz
     wz1 = 1-wz0
     th = wz0*thfnz(k)+wz1*thfnz(k+1)

     i = int((r-rin)/dr)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     k = int((th+pi)/dth)
     wz0 = (-pi+(k+1)*dth-th)/dth
     wz1 = 1.-wz0

     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
     b=1.-tor+tor*bfldp
     
     vfac = 0.5*(mbeam*u3b(m)**2 + 2.*mub(m)*b)

     if(itp==1)then
        grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
             +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 
        radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
             +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
        bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
             +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
        dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
             +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
        jfnp = wz0*jfn(k)+wz1*jfn(k+1)
        psipp = wx0*psip(i)+wx1*psip(i+1)        
        fp = wx0*f(i)+wx1*f(i+1)
        b=1.-tor+tor*bfldp
        enerb=(mub(m)+mbeam*vpar*vpar/b)/qbeam*tor
        zdot =  vpar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
             +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp
     end if

     dv=(dx*dy*dz)
     wght0 = 1/dv
     wght1 = 1/dv

     xt=x3b(m)
     yt=y3b(m)
     zt=z3b(m)
     i=int(xt/dx)
     j=int(yt/dy)
     k=0 !int(z3e(m)/dz)-gclr*kcnt

     if(itp==0)then
        wght1 = wght1*aparp*ppar*ppar
        myupa0(i,j,k)      =myupa0(i,j,k)+wght1*w000(m)
        myupa0(i+1,j,k)    =myupa0(i+1,j,k)+wght1*w100(m)
        myupa0(i,j+1,k)    =myupa0(i,j+1,k)+wght1*w010(m)
        myupa0(i+1,j+1,k)  =myupa0(i+1,j+1,k)+wght1*w110(m)
        myupa0(i,j,k+1)    =myupa0(i,j,k+1)+wght1*w001(m)
        myupa0(i+1,j,k+1)  =myupa0(i+1,j,k+1)+wght1*w101(m)
        myupa0(i,j+1,k+1)  =myupa0(i,j+1,k+1)+wght1*w011(m)
        myupa0(i+1,j+1,k+1)=myupa0(i+1,j+1,k+1)+wght1*w111(m)
     end if
     if(itp==1)then
        wght0 = wght0*aparp*ppar*zdot
        myupa(i,j,k)      =myupa(i,j,k)+wght0*w000(m)
        myupa(i+1,j,k)    =myupa(i+1,j,k)+wght0*w100(m)
        myupa(i,j+1,k)    =myupa(i,j+1,k)+wght0*w010(m)
        myupa(i+1,j+1,k)  =myupa(i+1,j+1,k)+wght0*w110(m)
        myupa(i,j,k+1)    =myupa(i,j,k+1)+wght0*w001(m)
        myupa(i+1,j,k+1)  =myupa(i+1,j,k+1)+wght0*w101(m)
        myupa(i,j+1,k+1)  =myupa(i,j+1,k+1)+wght0*w011(m)
        myupa(i+1,j+1,k+1)=myupa(i+1,j+1,k+1)+wght0*w111(m)

        wght0 = 1./dv*aparp*vpar
        myden0(i,j,k)      =myden0(i,j,k)+wght0*w000(m)
        myden0(i+1,j,k)    =myden0(i+1,j,k)+wght0*w100(m)
        myden0(i,j+1,k)    =myden0(i,j+1,k)+wght0*w010(m)
        myden0(i+1,j+1,k)  =myden0(i+1,j+1,k)+wght0*w110(m)
        myden0(i,j,k+1)    =myden0(i,j,k+1)+wght0*w001(m)
        myden0(i+1,j,k+1)  =myden0(i+1,j,k+1)+wght0*w101(m)
        myden0(i,j+1,k+1)  =myden0(i,j+1,k+1)+wght0*w011(m)
        myden0(i+1,j+1,k+1)=myden0(i+1,j+1,k+1)+wght0*w111(m)
     end if
  enddo

  !   enforce periodicity
  if(itp==1)call enforce(myupa(:,:,:))
  if(itp==1)call enforce(myden0(:,:,:))
  if(itp==0)call enforce(myupa0(:,:,:))
  !      call filter(myupa(:,:,:))

  if(itp==1)then
     do  i=0,im
        do  j=0,jm
           do  k=0,mykm
              upa0(i,j,k)= myupa(i,j,k)/n0b/jac(i,k)*cn0b
              den0apa(i,j,k)= myden0(i,j,k)/n0b/jac(i,k)*cn0b
           end do
        end do
     end do
  end if
  if(itp==0)then
     do  i=0,im
        do  j=0,jm
           do  k=0,mykm
              upa00(i,j,k)= myupa0(i,j,k)/n0b/jac(i,k)*cn0b
           end do
        end do
     end do
  end if

  upa0(:,:,:) = upa0(:,:,:)+ &
       (isg*phi(:,:,:)+dene(:,:,:))*apar(:,:,:)*nonline*0.

  upa00(:,:,:) = upa00(:,:,:)+ &
       (isg*phi(:,:,:)+dene(:,:,:))*apar(:,:,:)*nonline*0.
999 continue

  !      return
end subroutine jpar0b
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine collb(ip,n)

  use gem_com
  use gem_equil

  implicit none

  integer :: i,ip,k,m,n,ncol,icol
  real :: edum,vdum,dum,dum1,ptch,vti,r,qr,th,cost,b
  real :: h_x,h_coll,x,eps,dtcol,uold,hee,nue,ter
  real :: wx0,wx1,wz0,wz1

  ncol = 1
  if(ip.eq.1)dtcol = dt/ncol*2
  if(ip.eq.0)dtcol = dt/ncol
  if(rneui==0.0)return
  do k = 1,mmb
     r=x3b(k)-0.5*lx+lr0

     m = int(z3b(k)/delz)
     wz0 = ((m+1)*delz-z3b(k))/delz
     wz1 = 1-wz0
     th = wz0*thfnz(m)+wz1*thfnz(m+1)
     i = int((r-rin)/dr)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     m = int((th+pi)/dth)
     wz0 = (-pi+(m+1)*dth-th)/dth
     wz1 = 1.-wz0
     b = wx0*wz0*bfld(i,m)+wx0*wz1*bfld(i,m+1) &
          +wx1*wz0*bfld(i+1,m)+wx1*wz1*bfld(i+1,m+1) 
     uold = u3b(k)
     edum = b*mub(k)+0.5*mbeam*u3b(k)*u3b(k)
     vdum = sqrt(2.*edum/mbeam)
     ptch = u3b(k)/vdum
     nue=rneui
     dum1 = 1. !/(eps+0.1)**1.5  
     !         dum = dtcol*rneu*dum1
     dum = dtcol*nue*dum1
     !         if(x<0.3)dum=0.0
     do icol = 1,ncol
        ptch =ptch-dum*ptch+sqrt(dum*(1.-ptch*ptch)) &
             *sign(1.0,ran2(iseed)-0.5)
        ptch = min(ptch,0.999)
        ptch = max(ptch,-0.999)
     end do
     u3b(k) = vdum*ptch
     mub(k) = 0.5*mbeam*vdum*vdum*(1.-ptch*ptch)/b
     u2b(k) = u3b(k)
  end do
  !      return
end subroutine collb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ngt0(u)
  use gem_com
  use gem_equil
  use gem_fft_wrapper
  implicit none
  INTEGER :: i,j,k,k1,l,m,n
  real :: u(0:imx,0:jmx,0:1)
  real :: dum,dum1,r1,r2,r

  r1 = 0.25*a
  r2 = 0.3*a
  m = floor((r2-rin)/dx)
  do k = 0,1
     do i = 0,nxpp
        dum = 0.
        do j = 0,jm-1
           dum = dum+u(i,j,k)
        end do
        dum = dum/jmx
        do j = 0,jmx
           u(i,j,k) = u(i,j,k)-dum
        end do
     end do
  end do

  do k = 0,1
     do j = 0,jmx
        dum = u(m,j,k)
        do i = 0,imx
           r = rin+i*dx
           if(r<r1)u(i,j,k) = 0.
           if(r1<=r .and. r<=r2)u(i,j,k) = (r-r1)/(r2-r1)*dum
        end do
     end do
  end do
  
  !      return
end subroutine ngt0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
