module qfactor_common
	integer :: nx ,ny, nz, nxm1, nym1, nzm1,r0max(0:2) , qx, qy, qz, factor, &
	nbridges, q1, q2, mmaxsteps, maxsteps, CS_ND, q1tq2
	real:: xmax, ymax, zmax, xmin, ymin, zmin, pmin(0:2),pmax(0:2), NaN, &
	xreg(0:1),yreg(0:1),zreg(0:1), &
	point0(0:2), point1(0:2), point2(0:2), ev1(0:2), ev2(0:2), ev3(0:2), step, mstep, cut_coordinate, delta
	real(8),parameter:: pi=3.141592653589793D0
	integer, allocatable:: rsboundary(:,:), reboundary(:,:), tFlag(:, :)
	real, allocatable:: Bfield(:, :, :1, :), CurlB(:,:, :, :), &
	rsF(:, :, :), reF(:, :, :), bnr(:, :), twistF(:, :), q(:, :), qtmp(:, :), length(:, :)
	logical::twistFlag, vqflag, q0flag, csflag
end module

!CS_ND: normal direction of the cross section

!#########################################################################
! trilinear interpolation
subroutine interpolateB(vp, bp)
use qfactor_common
implicit none
real::w(0:1,0:2), weigh(0:1,0:1,0:1), vp(0:2), bp(0:2)
integer:: round(0:1,0:2), i, j, k

round(0,:)=floor(vp)
w(1,:)=vp-round(0,:)

do i=0,2
	if (vp(i) .le. 0.0) then
		round(0,i)=0
		w(1,i)=0.0
	endif

	if (vp(i) .ge. pmax(i)) then
		round(0,i)=r0max(i)
		w(1,i)=1.0
	endif
enddo

round(1,:)=round(0,:)+1
w(0,:)=1.0-w(1,:)
 
forall(i=0:1,j=0:1,k=0:1) weigh(i,j,k)=w(i,0)*w(j,1)*w(k,2)

forall(i=0:2) bp(i)=sum(Bfield(round(:,0) ,round(:,1) , round(:,2), i)*weigh)
end
!#########################################################################
 !d vp/dt=Bp/Bt  
subroutine RK4(dt, vp0, vp1)
implicit none
real ::  dt, vp0(0:2), vp1(0:2), bp(0:2), &
k1(0:2), k2(0:2), k3(0:2), k4(0:2)
call interpolateB(vp0, bp)
k1=bp/sqrt(sum(DBLE(bp)*bp))
call interpolateB(vp0+0.5*dt*k1, bp)
k2=bp/sqrt(sum(DBLE(bp)*bp))
call interpolateB(vp0+0.5*dt*k2, bp)
k3=bp/sqrt(sum(DBLE(bp)*bp))
call interpolateB(vp0+dt*k3, bp)
k4=bp/sqrt(sum(DBLE(bp)*bp))
vp1=vp0+dt/6.0*(k1+2*k2+2*k3+k4)
end
!#########################################################################
subroutine RK4_boundary(dt, vp0, vp1, b_dim)
implicit none
real ::  dt, vp0(0:2), vp1(0:2), bp(0:2), &
k1(0:2), k2(0:2), k3(0:2), k4(0:2)
integer:: b_dim

call interpolateB(vp0, bp)
k1=bp/bp(b_dim)
k1(b_dim)=1.0
call interpolateB(vp0+0.5*dt*k1, bp)
k2=bp/bp(b_dim)
k1(b_dim)=1.0
call interpolateB(vp0+0.5*dt*k2, bp)
k3=bp/bp(b_dim)
k3(b_dim)=1.0
call interpolateB(vp0+dt*k3, bp)
k4=bp/bp(b_dim)
k4(b_dim)=1.0
vp1=vp0+dt/6.0*(k1+2*k2+2*k3+k4)
end

!#########################################################################
subroutine interpolateCurlB(vp, Bp, CurlBp)
use qfactor_common
implicit none
real::w(0:1,0:2), weigh(0:1,0:1,0:1), vp(0:2), bp(0:2), CurlBp(0:2)
integer:: round(0:1,0:2), i, j, k

round(0,:)=floor(vp)
w(1,:)=vp-round(0,:)

do i=0,2
	if (vp(i) .le. 0.0) then
		round(0,i)=0
		w(1,i)=0.0
	endif

	if (vp(i) .ge. pmax(i)) then
		round(0,i)=r0max(i)
		w(1,i)=1.0
	endif
enddo

round(1,:)=round(0,:)+1
w(0,:)=1.0-w(1,:)
 
forall(i=0:1,j=0:1,k=0:1) weigh(i,j,k)=w(i,0)*w(j,1)*w(k,2)
forall(i=0:2) Bp(i)=sum(Bfield(round(:,0) ,round(:,1), round(:,2), i)*weigh)
forall(i=0:2) CurlBp(i)=sum(CurlB(round(:,0) ,round(:,1), round(:,2), i)*weigh)
end

!#########################################################################
subroutine curl_B(z)
use qfactor_common
implicit none
integer:: z
real:: dx, dy ,dz
real, allocatable:: dBxdy(:, :), dBxdz(:, :), dBydx(:, :), dBydz(:, :),&
 dBzdx(:, :), dBzdy(:, :), b2dx(:, :), b2dy(:, :), b2dz(:, :)

allocate(dBxdy(0:nxm1, 0:nym1))
allocate(dBxdz(0:nxm1, 0:nym1))
allocate(dBydx(0:nxm1, 0:nym1))
allocate(dBydz(0:nxm1, 0:nym1))
allocate(dBzdx(0:nxm1, 0:nym1))
allocate(dBzdy(0:nxm1, 0:nym1))

allocate(b2dx(0:nxm1, 0:nym1))
allocate(b2dy(0:nxm1, 0:nym1))
allocate(b2dz(0:nxm1, 0:nym1))

b2dx=Bfield(:,:,z,0)
b2dy=Bfield(:,:,z,1)
b2dz=Bfield(:,:,z,2)

dx=2.0
dy=2.0
dz=2.0

dBxdy(:,1:nym1-1) =(b2dx(:,2:nym1)-b2dx(:,0:nym1-2))/dy
dBxdy(:,0) = (-3.0*b2dx(:, 0) + 4.0*b2dx(:, 1) - b2dx(:, 2))/dy
dBxdy(:,nym1) = (3.0*b2dx(:, nym1) - 4.0*b2dx(:, nym1-1) + b2dx(:, nym1-2))/dy

dBzdy(:,1:nym1-1) =(b2dz(:,2:nym1)-b2dz(:,0:nym1-2))/dy
dBzdy(:,0) = (-3.0*b2dz(:, 0) + 4.0*b2dz(:, 1) - b2dz(:, 2))/dy
dBzdy(:,nym1) = (3.0*b2dz(:, nym1) - 4.0*b2dz(:, nym1-1) + b2dz(:, nym1-2))/dy

dBydx(1:nxm1-1,:) =(b2dy(2:nxm1,:)-b2dy(0:nxm1-2,:))/dx
dBydx(0,:) = (-3.0*b2dy(0,:) + 4.0*b2dy(1,:) - b2dy(2,:))/dx
dBydx(nxm1,:) = (3.0*b2dy(nxm1,:)  - 4.0*b2dy(nxm1-1,:)  + b2dy(nxm1-2,:) )/dx

dBzdx(1:nxm1-1,:) =(b2dz(2:nxm1,:)-b2dz(0:nxm1-2,:))/dx
dBzdx(0,:) = (-3.0*b2dz(0,:) + 4.0*b2dz(1,:) - b2dz(2,:))/dx
dBzdx(nxm1,:) = (3.0*b2dz(nxm1,:)  - 4.0*b2dz(nxm1-1,:)  + b2dz(nxm1-2,:) )/dx

	curlB(:,:,z,2) = dbydx - dbxdy

if(z==0) then
	dBydz(:, :) = (-3.0*Bfield(:,:,0,1) + 4.0*Bfield(:,:,1,1) - Bfield(:,:,2,1))/dz
	dBxdz(:, :) = (-3.0*Bfield(:,:,0,0) + 4.0*Bfield(:,:,1,0) - Bfield(:,:,2,0))/dz
	curlB(:,:,z,0) = dbzdy - dbydz
	curlB(:,:,z,1) = dbxdz - dbzdx
	deallocate(dBxdy, dBxdz, dBydx, dBydz, dBzdx, dBzdy, b2dx, b2dy, b2dz)
	return
endif

if(z==nzm1) then
	dBydz(:, :) = (3.0*Bfield(:,:,nzm1,1) - 4.0*Bfield(:,:,nzm1-1,1) + Bfield(:,:,nzm1-2,1))/dz
	dBxdz(:, :) = (3.0*Bfield(:,:,nzm1,0) - 4.0*Bfield(:,:,nzm1-1,0) + Bfield(:,:,nzm1-2,0))/dz
	curlB(:,:,z,0) = dbzdy - dbydz 
	curlB(:,:,z,1) = dbxdz - dbzdx
	deallocate(dBxdy, dBxdz, dBydx, dBydz, dBzdx, dBzdy, b2dx, b2dy, b2dz)
	return
endif

	dBydz(:, :) = (Bfield(:,:,z+1,1) -Bfield(:,:,z-1,1))/dz
	dBxdz(:, :) = (Bfield(:,:,z+1,0) -Bfield(:,:,z-1,0))/dz
  curlB(:,:,z,0) = dbzdy - dbydz
  curlB(:,:,z,1) = dbxdz - dbzdx   
  deallocate(dBxdy, dBxdz, dBydx, dBydz, dBzdx, dBzdy, b2dx, b2dy, b2dz)
END


!#########################################################################
subroutine trace_bline(vp0, rs, re, rbs, rbe, twist, line_length)

!+
! NAME :
!   trace_bline
!   
!PURPOSE:
!     Uses RK4 solver to integrate a line of force for a 3D magnetic field

!CATEGORY:

!INPUTS:
!     vp0 = the position to start the line at

!OUTPUTS:
!     rs   - pixel coordinates of the start point of field lines
!     re  - pixel coordinates of the end point of field lines
!     rbs,rbe- categorize field lines based on where they thread the boundary of the field cube
!     twist     - twist number of the field line

!SIDE EFFECTS:

!MODIFICATION HISTORY:
!       
!       January 13, 2014 R. Liu, adapted from Wieglemann's NLFFF package, trace_bline
!       
!       June 18, 2014 R. Liu, use pointer instead of structure for better efficency
!       
!       June 24, 2014 R. Liu, use common blocks for much better efficiency
!
!      March 20, 2015 Jun Chen, rboundary
!
!      June 15, 2015 Jun Chen, Fortran Edition, modify re, rs
!
!   This software is provided without any warranty. Permission to use,
!   copy, modify. Distributing modified or unmodified copies is granted,
!   provided this disclaimer and information are included unchanged.       
!-
use qfactor_common
implicit none
real ::  dt, dt1, vp0(0:2), vp1(0:2), vp(0:2), rs(0:2), re(0:2),&
 bp(0:2), bp1(0:2), rx ,ry, rz, CurlBp(0:2),twist, line_length
real, allocatable:: loop(:, :),alpha(:)
integer:: n, is, ie, rbs, rbe, t, m_ind, remainder, lenOver4, ind0
real(8)::alpha0,alpha1,alpha2,alpha3,alpha4,x1,x2,x3,x4,a0,a1,a2,a3,a4,loop0(0:2)
!!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 
allocate(loop(0:2, mmaxsteps:maxsteps) )
mstep=-step
ie=0
is=0
loop(:,0)=vp0
!!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 
!Integrate in +B direction
vp= vp0

do while( minval(vp-pmin)>=0 .and. maxval(vp-pmax)<=0 .and. ie < maxsteps)
	call RK4(step, vp, vp1) 
	vp=vp1
	ie=ie+1	
	loop(:,ie)=vp1
end do
re=vp1

!!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  
!Integrate in -B direction

vp= vp0
do while( minval(vp-pmin)>=0 .and. maxval(vp-pmax)<=0 .and. is > mmaxsteps)
	call RK4(mstep, vp, vp1)	
	vp=vp1	
	is=is-1	
	loop(:,is)=vp1
end do
rs=vp1

!!###############################################################
!rboundary
rbs=0
rx=rs(0)
ry=rs(1)
rz=rs(2)
  
if ((rx .ge. xmin) .and. (rx .le. xmax) .and. &
      (ry .ge. ymin) .and. (ry .le. ymax) .and. &
      (rz .lt. zmin)) then 
      rbs=1
endif

if ((rx .ge. xmin) .and. (rx .le. xmax) .and. &
      (ry .ge. ymin) .and. (ry .le. ymax) .and. &
      (rz .gt. zmax)) then 
      rbs=2
endif

if ((rx .ge. xmin) .and. (rx .le. xmax) .and. &
      (ry .lt. ymin) .and. &
      (rz .ge. zmin) .and.(rz .le. zmax)) then 
      rbs=3
endif

if ((rx .ge. xmin) .and. (rx .le. xmax) .and. &
      (ry .gt. ymax)  .and. &
      (rz .ge. zmin) .and.(rz .le. zmax)) then 
      rbs=4
endif

if ((rx .lt. xmin) .and. &
      (ry .ge. ymin) .and. (ry .le. ymax) .and. &
      (rz .ge. zmin) .and.(rz .le. zmax)) then 
      rbs=5
endif

if ((rx .gt. xmax) .and. &
      (ry .ge. ymin) .and. (ry .le. ymax) .and. &
      (rz .ge. zmin) .and.(rz .le. zmax)) then 
      rbs=6
endif

!!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,      
  rbe=0
  rx=re(0)
  ry=re(1)
  rz=re(2)
  
if ((rx .ge. xmin) .and. (rx .le. xmax) .and. &
      (ry .ge. ymin) .and. (ry .le. ymax) .and. &
      (rz .lt. zmin)) then 
      rbe=1
endif

if ((rx .ge. xmin) .and. (rx .le. xmax) .and. &
      (ry .ge. ymin) .and. (ry .le. ymax) .and. &
      (rz .gt. zmax)) then 
      rbe=2
endif

if ((rx .ge. xmin) .and. (rx .le. xmax) .and. &
      (ry .lt. ymin) .and. &
      (rz .ge. zmin) .and.(rz .le. zmax)) then 
      rbe=3
endif

if ((rx .ge. xmin) .and. (rx .le. xmax) .and. &
      (ry .gt. ymax)  .and. &
      (rz .ge. zmin) .and.(rz .le. zmax)) then 
      rbe=4
endif

if ((rx .lt. xmin) .and. &
      (ry .ge. ymin) .and. (ry .le. ymax) .and. &
      (rz .ge. zmin) .and.(rz .le. zmax)) then 
      rbe=5
endif

if ((rx .gt. xmax) .and. &
      (ry .ge. ymin) .and. (ry .le. ymax) .and. &
      (rz .ge. zmin) .and.(rz .le. zmax)) then 
      rbe=6
endif
!###############################################################
!modify re, rs
dt=1.0/16
ie=ie-1
vp=loop(:,ie)
do while( minval(vp-pmin)>=0 .and. maxval(vp-pmax)<=0 .and. ie < maxsteps)
	call RK4(dt, vp, vp1)
	vp=vp1	
	ie=ie+1	
	loop(:,ie)=vp1
end do

vp=loop(:,ie-1)
vp1=loop(:,ie)
if(rbe .eq. 0)then
	re=vp1
else
	m_ind=(6-rbe)/2
	
	if( mod(rbe,2) .eq. 1) then
	 	dt=pmin(m_ind)-vp(m_ind)
	 	dt1=vp1(m_ind)-pmin(m_ind)
	else
	 	dt=pmax(m_ind)-vp(m_ind)
	 	dt1=vp1(m_ind)-pmax(m_ind)
	endif
	
    call interpolateB(vp, bp)  
    call interpolateB(vp1, bp1)
    if((abs(bp(m_ind))*20 .le. sqrt(sum(DBLE(bp)*bp))) &
    	.or.(abs(bp1(m_ind))*20 .le. sqrt(sum(DBLE(bp1)*bp1))))then	   
			if((dt+dt1) .eq. 0.0) then
				re=vp
			else
				re=(vp*dt1+vp1*dt)/(dt+dt1)
			endif
		else
    	call RK4_boundary(dt, vp, re, m_ind)
    endif
endif
!!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 
dt=-1.0/16
is=is+1
vp=loop(:,is)
do while( minval(vp-pmin)>=0 .and. maxval(vp-pmax)<=0 .and. is >mmaxsteps)
	call RK4(dt, vp, vp1)    !d (x,y,z)/dt=(Bx,By,Bz)/Bs
	vp=vp1	
	is=is-1	
	loop(:,is)=vp1
end do

vp=loop(:,is+1)
vp1=loop(:,is)
if(rbs .eq. 0)then
	rs=vp1
else
	m_ind=(6-rbs)/2
	
	if( mod(rbs,2) .eq. 1) then
	 	dt=pmin(m_ind)-vp(m_ind)
	 	dt1=vp1(m_ind)-pmin(m_ind)
	else
	 	dt=pmax(m_ind)-vp(m_ind)
	 	dt1=vp1(m_ind)-pmax(m_ind)
	endif
	
    call interpolateB(vp, bp)
    call interpolateB(vp1, bp1)  
	if((abs(bp(m_ind))*20 .le. sqrt(sum(DBLE(bp)*bp))) & 
	.or.(abs(bp1(m_ind))* 20 .le. sqrt(sum(DBLE(bp1)*bp1))))then	  
	  if((dt+dt1) .eq. 0.0) then
	  	rs=vp
	  else
			rs=(vp*dt1+vp1*dt)/(dt+dt1)
		endif
	else
    	call RK4_boundary(dt, vp, rs, m_ind)
    endif
endif
!!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 
n=ie-is+1
	if ( ( sum( ( loop(:, ie-1)-re )**2.0 ) .lt. 0.001 ) .and.  (ie .ne. is+1 ) ) then 
		ie=ie-1
	else	
		loop(:, ie)=re
	endif
	
	if ( ( sum( ( loop(:, is+1)-rs )**2.0 ) .lt. 0.001 ) .and.  (is .ne. ie-1 ) ) then 
		is=is+1
	else	
		loop(:, is)=rs
	endif
	n=ie-is+1

!#######################################################################
!compute length of field line
line_length=0.0
do t=is, ie-1
	line_length=line_length+sqrt(sum( (loop(:,t)-loop(:,t+1))**2.0 ))
enddo

!#######################################################################

!compute twist
twist=0.0
if (twistFlag) then
    if ((n == 2).or. (n >= maxsteps)) then
				if (n == 2) then  
					twist=0.0
				else 
					twist=NAN   		
				endif
		else
			twist=0.0		
			
!		  do t=is, ie-1
!		  		vp=loop(:,t)
!		  		call interpolateCurlB(vp,bp,curlbp)
!	  		  twist=twist+sqrt(sum( (vp-loop(:,t+1))**2.0 ))*sum(curlbp*bp)/(sqrt(sum(DBLE(bp)*bp))**2)
!		  enddo 	   
 
!========================================================
!Newton-Cotes integral
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 	
remainder=mod(ie-is,4)
lenOver4=(ie-is)/4

allocate(alpha(mmaxsteps:maxsteps) )
do t=is, ie
	vp=loop(:,t)
	call interpolateCurlB(vp, bp, curlbp)
	alpha(t)=sum(curlbp*bp)/(sum(DBLE(bp)*bp))
enddo
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 	
do t=0,lenOver4-1
ind0=is+4*t
loop0=loop(:,ind0)
x1=sqrt(sum( (loop(:,ind0+1)-loop0)**2))
x2=sqrt(sum( (loop(:,ind0+2)-loop0)**2))
x3=sqrt(sum( (loop(:,ind0+3)-loop0)**2))
x4=sqrt(sum( (loop(:,ind0+4)-loop0)**2))

a0=(x4*(30*x1*x2*x3-10*(x2*x3+x1*(x2+x3))*x4+5*(x1+x2+x3)*x4**2-3*x4**3))/(60*x1*x2*x3)
a1=(x4**3*(-10*x2*x3+5*(x2+x3)*x4-3*x4**2))/(60*x1*(x1-x2)*(x1-x3)*(x1-x4))
a2=(x4**3*(-10*x1*x3+5*(x1+x3)*x4-3*x4**2))/(60*x2*(-x1+x2)*(x2-x3)*(x2-x4))
a3=(x4**3*(-10*x1*x2+5*(x1+x2)*x4-3*x4**2))/(60*x3*(-x1+x3)*(-x2+x3)*(x3-x4))
a4=x4-a3-a2-a1-a0
!(x4*(-30*x1*x2*x3+20*(x2*x3+x1*(x2+x3))*x4-15*(x1+x2+x3)*x4**2+12*x4**3))/(60*(-x1+x4)*(-x2+x4)*(-x3+x4))
alpha0=alpha(ind0)
alpha1=alpha(ind0+1)
alpha2=alpha(ind0+2)
alpha3=alpha(ind0+3)
alpha4=alpha(ind0+4)
twist=twist+a0*alpha0+a1*alpha1+a2*alpha2+a3*alpha3+a4*alpha4
enddo

ind0=is+4*lenOver4
loop0=loop(:,ind0)
select case(remainder)
	case(0) 
		twist=twist
	case(1) 
		x1=sqrt(sum( (loop(:,ind0+1)-loop0)**2))
		a0=x1/2
		a1=a0 !x1/2
		alpha0=alpha(ind0)
		alpha1=alpha(ind0+1)
		twist=twist+a0*alpha0+a1*alpha1
	case(2) 
		x1=sqrt(sum( (loop(:,ind0+1)-loop0)**2))
		x2=sqrt(sum( (loop(:,ind0+2)-loop0)**2))

		a0=((3*x1-x2)*x2)/(6*x1)
		a1=-(x2**3/(6*x1**2-6*x1*x2))
		a2=x2-a1-a0 !((3*x1-2*x2)*x2)/(6*(x1-x2))

		alpha0=alpha(ind0)
		alpha1=alpha(ind0+1)
		alpha2=alpha(ind0+2)
		twist=twist+a0*alpha0+a1*alpha1+a2*alpha2
	case(3) 
		x1=sqrt(sum( (loop(:,ind0+1)-loop0)**2))
		x2=sqrt(sum( (loop(:,ind0+2)-loop0)**2))
		x3=sqrt(sum( (loop(:,ind0+3)-loop0)**2))

		a0=(x3*(6*x1*x2-2*(x1+x2)*x3+x3**2))/(12*x1*x2)
		a1=((2*x2-x3)*x3**3)/(12*x1*(x1-x2)*(x1-x3))
		a2=(x3**3*(-2*x1+x3))/(12*(x1-x2)*x2*(x2-x3))
		a3=x3-a2-a1-a0 !(x3*(-6*x1*x2+4*(x1+x2)*x3-3*x3**2))/(12*(x1-x3)*(-x2+x3))

		alpha0=alpha(ind0)
		alpha1=alpha(ind0+1)
		alpha2=alpha(ind0+2)
		alpha3=alpha(ind0+3)
		twist=twist+a0*alpha0+a1*alpha1+a2*alpha2+a3*alpha3
end select		
deallocate(alpha)
!========================================================
 			twist=twist/(4.0*pi)	 
 		endif 
endif  

deallocate(loop)

END


!#########################################################################
subroutine initialize()
use qfactor_common
implicit none
 character(len=30):: file_b3d, file_head
 integer:: k,rec1, times(8), twistFlag_int, csFlag_int
 character(len=30):: time_str(3)
 real(8):: length01,length02
 real::vp(0:2),bp(0:2)
!!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 
file_b3d='b3d.bin'
file_head='head.txt'
open(unit=8, file=file_head, status='unknown')
read(8, *) nx, ny ,nz ,nbridges, factor, xreg, yreg, zreg, twistFlag_int, step, &
  qx, qy, qz, q1, q2, csFlag_int, point0, point1, point2
 close(8, status='KEEP')
 
twistFlag= twistFlag_int .eq. 1
 csFlag= csFlag_int .eq. 1 
 vqflag=(qx .ne. 1).and.(qy .ne. 1).and. (qz .ne. 1) .and. ( .not. csflag)
q0flag=( zreg(0) .eq. 0) .and. ( zreg(1) .eq. 0) .and. ( .not. csflag)

if (csflag) then
	length01=sqrt(sum(DBLE(point1-point0)*(point1-point0)))
	length02=sqrt(sum(DBLE(point2-point0)*(point2-point0)))
	 !ev*: Unit basis vector
	ev1=(point1-point0)/length01
	ev2=(point2-point0)/length02
	ev3(0)=dble(ev1(1))*ev2(2)-dble(ev1(2))*ev2(1)
	ev3(1)=dble(ev1(2))*ev2(0)-dble(ev1(0))*ev2(2)
	ev3(2)=dble(ev1(0))*ev2(1)-dble(ev1(1))*ev2(0)
else
	if(vqflag) CS_ND=2	
	
	if (qx .eq. 1) then 
		CS_ND=0; cut_coordinate=xreg(0)
	endif
	if (qy .eq. 1) then 
		CS_ND=1; cut_coordinate=yreg(0)
	endif
	if (qz .eq. 1) then 
		CS_ND=2; cut_coordinate=zreg(0)
	endif
endif

NaN=-log(-1.0)
q1tq2=q1*q2
delta=1.0/factor
nxm1=nx-1
nym1=ny-1
nzm1=nz-1
xmin=0
ymin=0
zmin=0
pmin=[xmin,ymin,zmin]
xmax=nxm1
ymax=nym1
zmax=nzm1
pmax=[xmax,ymax,zmax]
r0max=[nxm1,nym1,nzm1]-1

maxsteps=nint(4*(nx+ny+nz)/step)
mmaxsteps=-maxsteps

CALL OMP_set_num_threads(nbridges)

allocate(Bfield(0:nxm1, 0:nym1, 0:nzm1, 0:2))
allocate(q(0:q1-1, 0:q2-1))
allocate(qtmp(0:q1-1, 0:q2-1))
allocate(reF(0:2, -1:q1, -1:q2))
allocate(reboundary(-1:q1, -1:q2))
allocate(bnr(0:q1-1, 0:q2-1))
allocate(length(0:q1-1, 0:q2-1))
if(twistFlag)	allocate(twistF(0:q1-1, 0:q2-1))
!!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 
! read Bx, By, Bz
!max recl could be 2^29-1 ~ 5.36e8 ~ correspond to 2GiB

if(nx*ny*nz< 2**29-1) then
	open(unit=8, file=file_b3d, access="direct",recl=nx*ny*nz)
	read(8,rec=1) Bfield(:,:,:,0)
	read(8,rec=2) Bfield(:,:,:,1)
	read(8,rec=3) Bfield(:,:,:,2)
else	
	open(unit=8, file=file_b3d, access="direct",recl=nx*ny)
	do rec1=0, nz-1 
		read(8,rec=1+rec1) Bfield(:, :,rec1,0)
		read(8,rec=1+rec1+nz) Bfield(:, :,rec1,1)
		read(8,rec=1+rec1+2*nz) Bfield(:, :,rec1,2)
	enddo
endif

 close(8, status='KEEP')
 
 !call system('rm '//file_b3d)
!!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  
if(twistFlag)then
	allocate(curlB(0:nxm1, 0:nym1, 0:nzm1, 0:2))
	
	!$OMP PARALLEL DO  PRIVATE(k)
	DO k = 0, nzm1
		call curl_B(k) 
	END DO	
	!$OMP END PARALLEL DO
endif
!!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 
print*, '  _____________________________________'
print*, '        schedule         time'
call date_and_time( time_str(1) , time_str(2),  time_str(3), times )		 
print 300, 0.0, times(5), times(6), times(7)
300 format( '         ', F6.2, '%        ' I2.2, ':', I2.2, ':', I2.2)  
end
