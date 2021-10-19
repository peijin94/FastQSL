module qfactor_common
	integer:: qx, qy, qz, factor, &
	nbridges, q1, q2, q1m1, q2m1, q1tq2, CS_ND
	real:: xreg(0:1),yreg(0:1),zreg(0:1), NaN, cut_coordinate, delta, &
	point0(0:2), point1(0:2), point2(0:2), ev1(0:2), ev2(0:2), ev3(0:2)	
	integer, allocatable:: rsboundary(:,:), reboundary(:,:), tFlag(:, :)
	real, allocatable:: rsF(:, :, :), reF(:, :, :), bnr(:, :), &
	q_perp(:,:), q(:, :), qtmp(:, :), length(:, :), twistF(:, :)
	logical::twistFlag, vflag, q0flag, csflag, scottFlag
end module qfactor_common


module trace_common
	integer:: nx ,ny, nz, nxm1, nym1, nzm1, r0max(0:2), mmaxsteps, maxsteps
	real:: xmax, ymax, zmax, xmin, ymin, zmin, pmin(0:2),pmax(0:2), step, tol
	logical:: RK4flag
	real(8),parameter:: pi=3.141592653589793D0
end module trace_common


module field_common
	real, allocatable:: Bfield(:, :, :, :), CurlB(:,:, :, :), gradUVBfield(:, :, :, :, :)
end module field_common


module rkf45_common
	real:: c2,c3,c4,c5,c6, &
	a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65, &
	b1,b3,b4,b5,b6,&
	ce1,ce3,ce4,ce5,ce6;
end module rkf45_common


! trilinear interpolation
subroutine interpolateB(vp, bp)
use trace_common
use field_common
implicit none
real::w(0:1,0:2), weigh(0:1,0:1,0:1), vp(0:2), bp(0:2)
integer:: round(0:1,0:2), i, j, k
!----------------------------------------------------------------------------
round(0,:)=floor(vp)
w(1,:)=vp-round(0,:)

do i=0,2
	if (vp(i) .lt. 0.0) then
		round(0,i)=0	
		w(1,i)=0.0		
	else if (vp(i) .ge. pmax(i)) then
		round(0,i)=r0max(i)
		w(1,i)=1.0
	endif
enddo

round(1,:)=round(0,:)+1
w(0,:)=1.0-w(1,:)

forall(i=0:1,j=0:1,k=0:1) weigh(i,j,k)=w(i,0)*w(j,1)*w(k,2)
forall(i=0:2) Bp(i)=sum(Bfield(i, round(:,0) ,round(:,1), round(:,2))*weigh)
end subroutine interpolateB


subroutine interpolateUVB(vp, uvbp)
!Unit vector B
implicit none
real::vp(0:2), uvbp(0:2), bp(0:2)
call interpolateB(vp, bp)
uvbp= bp/norm2(bp) 
end subroutine interpolateUVB


subroutine interpolateAlpha(vp, uvbp, alpha, alphaFlag)
use trace_common
use field_common
implicit none
real::w(0:1,0:2), weigh(0:1,0:1,0:1), vp(0:2), bp(0:2), uvbp(0:2), CurlBp(0:2), alpha
integer:: round(0:1,0:2), i, j, k
logical:: alphaFlag
!----------------------------------------------------------------------------
round(0,:)=floor(vp)
w(1,:)=vp-round(0,:)

do i=0,2
	if (vp(i) .lt. 0.0) then
		round(0,i)=0
		w(1,i)=0.0
	else if (vp(i) .ge. pmax(i)) then
		round(0,i)=r0max(i)
		w(1,i)=1.0
	endif
enddo

round(1,:)=round(0,:)+1



w(0,:)=1.0-w(1,:)



forall(i=0:1,j=0:1,k=0:1) weigh(i,j,k)=w(i,0)*w(j,1)*w(k,2)
forall(i=0:2) Bp(i)=sum(Bfield(i, round(:,0) ,round(:,1), round(:,2))*weigh)
uvbp=bp/norm2(bp)

if (alphaFlag) then
	forall(i=0:2) CurlBp(i)=sum(CurlB(i, round(:,0) ,round(:,1), round(:,2))*weigh)
	alpha=sum(DBLE(curlbp)*bp)/(sum(DBLE(bp)*bp))
endif
end subroutine interpolateAlpha


!d vp/dt=Bp/Bt
subroutine RK4(dt, vp0, vp1, alpha, alphaFlag)
implicit none
real ::  dt, alpha, vp0(0:2), vp1(0:2), &
k1(0:2), k2(0:2), k3(0:2), k4(0:2), bp(0:2)
logical:: alphaFlag
!----------------------------------------------------------------------------
call interpolateAlpha(vp0, k1, alpha, alphaFlag)
call interpolateUVB(vp0+dt*k1/3., k2)
call interpolateUVB(vp0+dt*(-k1/3.+k2), k3)
call interpolateUVB(vp0+dt*(k1-k2+k3), k4)
vp1=vp0+dt/8.0*(k1+3.*(k2+k3)+k4)
end subroutine RK4


subroutine RK4_Boundary(dt, vp0, vp1, b_dim)
implicit none
real ::  dt, vp0(0:2), vp1(0:2), bp(0:2), &
k1(0:2), k2(0:2), k3(0:2), k4(0:2)
integer:: b_dim
!----------------------------------------------------------------------------
call interpolateB(vp0, bp)
k1=bp/bp(b_dim)
k1(b_dim)=1.0
call interpolateB(vp0+dt*k1/3., bp)
k2=bp/bp(b_dim)
k2(b_dim)=1.0
call interpolateB(vp0+dt*(-k1/3.+k2), bp)
k3=bp/bp(b_dim)
k3(b_dim)=1.0
call interpolateB(vp0+dt*(k1-k2+k3), bp)
k4=bp/bp(b_dim)
k4(b_dim)=1.0
vp1=vp0+dt/8.0*(k1+3.*(k2+k3)+k4)
end subroutine RK4_Boundary


subroutine RKF45(dt, vp0, vp1, alpha, alphaFlag, tol_this)
use trace_common
use rkf45_common
implicit none
real:: k1(0:2),k2(0:2),k3(0:2),k4(0:2),k5(0:2),k6(0:2),vp0(0:2),vp1(0:2)
real:: dt, dt0, dt1, alpha, step_error, min_step_error, tol_this
real:: dvp(0:2), scale_dt
logical:: continue_flag, alphaFlag
integer:: rb, rb_index
!----------------------------------------------------------------------------

continue_flag=.true.

call interpolateAlpha(vp0, k1, alpha, alphaFlag)
do while ( continue_flag ) 

	call interpolateUVB(vp0+dt*a21*k1, k2)
	call interpolateUVB(vp0+dt*(a31*k1+ a32*k2),k3)   
	call interpolateUVB(vp0+dt*(a41*k1+ a42*k2+ a43*k3), k4)   
	call interpolateUVB(vp0+dt*(a51*k1+ a52*k2+ a53*k3+ a54*k4),k5)   
	call interpolateUVB(vp0+dt*(a61*k1+ a62*k2+ a63*k3+ a64*k4+ a65*k5),k6)

	vp1 = vp0 + dt*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6)
	dvp=dt*(ce1*k1+ce3*k3+ce4*k4+ce5*k5+ce6*k6)
	step_error = norm2(dvp)
	
	continue_flag =.false.
	if (abs(dt) .gt. step) then
		
		if (.not.(minval(vp1-pmin)>=0 .and. maxval(vp1-pmax)<=0 ))  then	
			call vp_rboundary(vp1, rb, rb_index)		
			if (rb .ne. 7) then 
				if( mod(rb, 2) .eq. 1) then
				 	dt0=pmin(rb_index)-vp0(rb_index)
				 	dt1=vp1(rb_index)-pmin(rb_index)
				else
				 	dt0 =pmax(rb_index)-vp0(rb_index)
				 	dt1=vp1(rb_index)-pmax(rb_index)
				endif
				dt=dt*dt0/(dt0+dt1)*0.98
				
				if (abs(dt)  .lt. step) then	
					dt=sign(step,dt)			
					continue_flag=.false.
				else
					continue_flag=.true.
				endif
			else
				dt=dt*0.618
				if (abs(dt)  .lt. step) dt=sign(step,dt)
				continue_flag=.true.
			endif
			
			cycle

		endif
		
		if  (step_error .gt. tol_this) then
			dt=dt*0.618
			if (abs(dt)  .lt. step) dt=sign(step,dt)
			continue_flag=.true.
			cycle
		endif
	endif
	
enddo

!11.0932=(0.618)^-5
min_step_error=(tol_this/11.0932)/((100./abs(dt))**5.)
if (step_error .lt. min_step_error) then
	dt=sign(100.,dt)
else
	scale_dt=((tol_this/step_error)**0.2)*0.618
	dt=dt*scale_dt
endif

if (abs(dt)  .lt. step) dt=sign(step,dt)

end subroutine RKF45


subroutine modify_foot(vp, vp1, sign_dt, rb)
use trace_common
implicit none
real :: dt, dt1, alpha
real :: vp(0:2), vp1(0:2), vp1_orig(0:2), vp_tmp0(0:2), vp_tmp(0:2), uvbp(0:2), uvbp1(0:2)
integer:: sign_dt, rb, rb_index, maxsteps2, it
!----------------------------------------------------------------------------
call vp_rboundary(vp1, rb, rb_index)
if (rb .eq. 0) return

vp1_orig=vp1
dt=1.0/16* sign_dt
vp_tmp=vp
it=0
maxsteps2=norm2(vp1-vp)/abs(dt)*2
do while( minval(vp_tmp-pmin)>=0 .and. maxval(vp_tmp-pmax)<=0)
	call RK4(dt, vp_tmp, vp1, alpha, .false.)
	it=it+1
	if (it .ge. maxsteps2) then
		vp_tmp0=vp
		vp1=vp1_orig
		exit
	endif	
	vp_tmp0=vp_tmp
	vp_tmp =vp1
end do
call vp_rboundary(vp1, rb, rb_index)
if (rb .eq. 7) return
	
if( mod(rb, 2) .eq. 1) then
	dt =pmin(rb_index)-vp_tmp0(rb_index)
	dt1= vp1(rb_index)-   pmin(rb_index)
else
	dt =pmax(rb_index)-vp_tmp0(rb_index)
 	dt1= vp1(rb_index)-   pmax(rb_index)
endif
	
call interpolateUVB(vp_tmp0, uvbp)  
call interpolateUVB(vp1, uvbp1)
	
if ((abs(uvbp (rb_index)) .le. 0.05) &
.or.(abs(uvbp1(rb_index)) .le. 0.05))then
	if((dt+dt1) .ne. 0.0) &
	vp1=(vp_tmp0*dt1+vp1*dt)/(dt+dt1)
else
	call RK4_Boundary(dt, vp_tmp0, vp1, rb_index)
endif

end subroutine modify_foot


subroutine trace_bline(vp0, rs, re, rbs, rbe, line_length, twist, twistFlag)
!+
! NAME :
!   trace_bline
!   
!PURPOSE:
!     Uses RK4 solver to integrate a line of force for a 3D magnetic field

!INPUTS:
!     vp0        - the position to start the line at
!     twistFlag  - compute the twist or not, additional 15.4% time is required

!OUTPUTS:
!     rs   - pixel coordinates of the start point of field lines
!     re  - pixel coordinates of the end point of field lines
!     rbs,rbe- categorize field lines based on where they thread the boundary of the field cube
!     twist     - twist number of the field line
!     line_length  -length of the field line

!MODIFICATION HISTORY:
!       
!     January 13, 2014 R. Liu, adapted from Wieglemann's NLFFF package, trace_bline
!       
!     June 18, 2014 R. Liu, use pointer instead of structure for better efficency
!       
!     June 24, 2014 R. Liu, use common blocks for much better efficiency
! 
!     March 20, 2015 Jun Chen, rboundary
! 
!     June 15, 2015 Jun Chen, Fortran Edition, modify footpoints
! 
!   This software is provided without any warranty. Permission to use,
!   copy, modify. Distributing modified or unmodified copies is granted,
!   provided this disclaimer and information are included unchanged.       
!-
use trace_common
implicit none
real ::  dt, vp0(0:2), vp1(0:2), vp(0:2), vp_tmp(0:2), rs(0:2), re(0:2), uvbp(0:2), uvbp_z,&
twist, line_length, alpha, alpha0, dL, dL0, dtwist
integer:: it, rb, rbs, rbe, sign_dt
real:: tol_this, incline
logical::twistFlag, z0Flag
!----------------------------------------------------------------------------
twist=0.
line_length=0.
!----------------------------------------------------------------------------
z0flag= vp0(2) .eq. 0.0
if (z0flag) then 
	call interpolateUVB(vp0,uvbp)
	uvbp_z=uvbp(2)
	incline=abs(uvbp_z)
	if (incline .le. 0.05)then
		tol_this=tol/(8.e3)
	else 
		tol_this=tol*(incline**3.)
	endif
else
	tol_this=tol
endif

do sign_dt=-1,1,2
	if (z0flag) then
		if ( uvbp_z*sign_dt .le. 0) then 
			if (sign_dt .eq. -1) then
				rs =vp0
				rbs=1
			else
				re =vp0
				rbe=1
			endif		
			cycle
		endif
	endif
	
	vp= vp0
	it=0	
	dt=step*sign_dt
	dL=0.
	do while( minval(vp-pmin)>=0 .and. maxval(vp-pmax)<=0 .and. it < maxsteps)
		line_length=line_length+dL
				
		if (RK4flag) then  
		 	call RK4  (dt, vp, vp1, alpha, twistflag) 
		 	! alpha @ vp
		else			
			call RKF45(dt, vp, vp1, alpha, twistflag, tol_this)
		endif
		
		dL0=dL
		dL=norm2(vp1-vp)
				
		if (twistflag) then 
			if (it .ne. 0) then
				dtwist=(alpha0+alpha)/2.*dL0
				twist=twist+dtwist
			endif
			
			alpha0=alpha
		endif
		
		it=it+1	
		vp_tmp=vp
		vp=vp1
	end do
	call modify_foot(vp_tmp, vp1, sign_dt, rb)
	dL=norm2(vp1-vp_tmp)
	
	line_length=line_length+dL
	if (twistflag) then 
		call interpolateAlpha(vp1, uvbp, alpha, twistflag)
		dtwist=(alpha0+alpha)/2.*dL
		twist=twist+dtwist
	endif
	
	if (sign_dt .eq. -1) then
		rs =vp1
		rbs=rb
	else
		re =vp1
		rbe=rb
	endif
enddo

if (twistflag) twist=twist/(4.0*pi)

END subroutine trace_bline


subroutine curlB_calculate(k)
use trace_common
use field_common
implicit none
integer:: i,j,k
real:: dBxdy, dBxdz, dBydx, dBydz, dBzdx, dBzdy
!----------------------------------------------------------------------------
do j=0,nym1
do i=0,nxm1
!----------------------------------------------------------------------------
	if (i .eq. 0) then
		dBydx =-3.0*Bfield(1,0,j,k)   +4.0*Bfield(1,1,j,k)     -Bfield(1,2,j,k)
		dBzdx =-3.0*Bfield(2,0,j,k)   +4.0*Bfield(2,1,j,k)     -Bfield(2,2,j,k)
	else if (i .eq. nxm1) then
		dBydx = 3.0*Bfield(1,nxm1,j,k)-4.0*Bfield(1,nxm1-1,j,k)+Bfield(1,nxm1-2,j,k)
		dBzdx = 3.0*Bfield(2,nxm1,j,k)-4.0*Bfield(2,nxm1-1,j,k)+Bfield(2,nxm1-2,j,k)
	else
		dBydx = Bfield(1,i+1,j,k)-Bfield(1,i-1,j,k)
		dBzdx = Bfield(2,i+1,j,k)-Bfield(2,i-1,j,k)	
	endif
!----------------------------------------------------------------------------
	if (j .eq. 0) then
		dBxdy =-3.0*Bfield(0,i,0,k)   +4.0*Bfield(0,i,1,k)     -Bfield(0,i,2,k)
		dBzdy =-3.0*Bfield(2,i,0,k)   +4.0*Bfield(2,i,1,k)     -Bfield(2,i,2,k)
	else if (j .eq. nym1) then
		dBxdy = 3.0*Bfield(0,i,nym1,k)-4.0*Bfield(0,i,nym1-1,k)+Bfield(0,i,nym1-2,k)
		dBzdy = 3.0*Bfield(2,i,nym1,k)-4.0*Bfield(2,i,nym1-1,k)+Bfield(2,i,nym1-2,k)
	else
		dBxdy = Bfield(0,i,j+1,k)-Bfield(0,i,j-1,k)
		dBzdy = Bfield(2,i,j+1,k)-Bfield(2,i,j-1,k)	
	endif
!----------------------------------------------------------------------------
	if (k .eq. 0) then
		dBydz =-3.0*Bfield(1,i,j,0)   +4.0*Bfield(1,i,j,1)     -Bfield(1,i,j,2)
		dBxdz =-3.0*Bfield(0,i,j,0)   +4.0*Bfield(0,i,j,1)     -Bfield(0,i,j,2)
	else if (k .eq. nzm1) then
		dBydz = 3.0*Bfield(1,i,j,nzm1)-4.0*Bfield(1,i,j,nzm1-1)+Bfield(1,i,j,nzm1-2)
		dBxdz = 3.0*Bfield(0,i,j,nzm1)-4.0*Bfield(0,i,j,nzm1-1)+Bfield(0,i,j,nzm1-2)
	else
		dBydz = Bfield(1,i,j,k+1)-Bfield(1,i,j,k-1)
		dBxdz = Bfield(0,i,j,k+1)-Bfield(0,i,j,k-1)
	endif
!----------------------------------------------------------------------------
	curlB(0,i,j,k) = dbzdy - dbydz
	curlB(1,i,j,k) = dbxdz - dbzdx
	curlB(2,i,j,k) = dbydx - dbxdy
	curlB(:,i,j,k)=curlB(:,i,j,k)/2.
enddo
enddo

END subroutine curlB_calculate


subroutine vp_rboundary(vp, rb, rb_index)
use trace_common
implicit none
real:: vp(0:2),rx,ry,rz
integer::rb, rb_index
!----------------------------------------------------------------------------
if (minval(vp-pmin)>=0 .and. maxval(vp-pmax)<=0 ) then
	rb=0
	return
endif

rx=vp(0)
ry=vp(1)
rz=vp(2)
  
if ((rx .ge. xmin) .and. (rx .le. xmax) .and. &
    (ry .ge. ymin) .and. (ry .le. ymax) .and. &
    (rz .lt. zmin)) then 
      rb=1
      rb_index=2
      return
endif

if ((rx .ge. xmin) .and. (rx .le. xmax) .and. &
    (ry .ge. ymin) .and. (ry .le. ymax) .and. &
    (rz .gt. zmax)) then 
      rb=2
      rb_index=2
      return
endif

if ((rx .ge. xmin) .and. (rx .le. xmax) .and. &
    (ry .lt. ymin) .and. &
    (rz .ge. zmin) .and.(rz .le. zmax)) then 
      rb=3
      rb_index=1
      return
endif

if ((rx .ge. xmin) .and. (rx .le. xmax) .and. &
    (ry .gt. ymax)  .and. &
    (rz .ge. zmin) .and.(rz .le. zmax)) then 
      rb=4
      rb_index=1
      return
endif

if ((rx .lt. xmin) .and. &
    (ry .ge. ymin) .and. (ry .le. ymax) .and. &
    (rz .ge. zmin) .and.(rz .le. zmax)) then 
      rb=5
      rb_index=0
      return
endif

if ((rx .gt. xmax) .and. &
    (ry .ge. ymin) .and. (ry .le. ymax) .and. &
    (rz .ge. zmin) .and.(rz .le. zmax)) then 
      rb=6
      rb_index=0
      return
endif

rb=7

end subroutine vp_rboundary


subroutine initialize()
use qfactor_common
use trace_common
use field_common
use rkf45_common
implicit none
integer::k, s, twistFlag_int, csFlag_int, RK4flag_int, scottFlag_int
real:: length01, length02
real, allocatable:: Bfield_tmp(:, :, :, :)
!----------------------------------------------------------------------------
open(unit=8, file='head.txt', status='unknown')
read(8, *) nx, ny, nz, nbridges, factor, maxsteps, &
           xreg, yreg, zreg, step, tol, &
           qx, qy, qz, q1, q2, &
           twistFlag_int, RK4flag_int, scottFlag_int, csFlag_int, &
           point0, point1, point2
 close(8)
!----------------------------------------------------------------------------
twistFlag= twistFlag_int .eq. 1
RK4flag=RK4flag_int .eq. 1
if (.not. RK4Flag) step=0.25
scottFlag= scottFlag_int .eq. 1
 csFlag= csFlag_int .eq. 1 
vflag=(qx .ne. 1) .and. (qy .ne. 1) .and. (qz .ne. 1) .and. ( .not. csflag)
q0flag=(zreg(0) .eq. 0.0) .and. (zreg(1) .eq. 0.0) .and. ( .not. csflag)

if (q0flag) cut_coordinate=0.0
if (csflag) then
	length01=norm2(point1-point0)
	length02=norm2(point2-point0)

	 !ev*: unit basis vector of the cut plane
	ev1=(point1-point0)/length01
	ev2=(point2-point0)/length02
	ev3(0)=dble(ev1(1))*ev2(2)-dble(ev1(2))*ev2(1)
	ev3(1)=dble(ev1(2))*ev2(0)-dble(ev1(0))*ev2(2)
	ev3(2)=dble(ev1(0))*ev2(1)-dble(ev1(1))*ev2(0)
else
	!CS_ND: normal direction of the cross section
	if(vflag) CS_ND=2	
	
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

!----------------------------------------------------------------------------
! for RKF45
a21=1./4.
a31=3./32.;      a32=9./32.
a41=1932./2197.; a42=-7200./2197.;  a43= 7296./2197.
a51=439./216.;   a52= -8.;          a53=  3680./513.;    a54=-845./4104.
a61= -8./27.;    a62= 2.;           a63= -3544./2565.;   a64= 1859./4104.; a65= -11./40.
b1 = 16./135.;   b3 = 6656./12825.; b4 = 28561./56430.;  b5 = -9./50.;     b6 = 2./55. 
 ce1 = 1./360.;  ce3= -128./4275.; ce4 = -2197./75240.; ce5 = 1./50.;     ce6 = 2./55.
!----------------------------------------------------------------------------

NaN=transfer(2143289344, 1.0)
q1tq2=q1*q2
q2m1=q2-1
q1m1=q1-1
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

!maxsteps=nint(4*(nx+ny+nz)/step)*10
mmaxsteps=-maxsteps

CALL OMP_set_num_threads(nbridges)

allocate(     q(0:q1-1, 0:q2-1))
allocate(  qtmp(0:q1-1, 0:q2-1))
allocate(   bnr(0:q1-1, 0:q2-1))
allocate(length(0:q1-1, 0:q2-1))
allocate(reF(0:2, 0:q1-1, 0:q2-1))
allocate(reboundary(-2:q1+1, -2:q2+1))

if (twistFlag) allocate(twistF(0:q1-1, 0:q2-1))
if (scottFlag) allocate(q_perp(0:q1-1, 0:q2-1))
!----------------------------------------------------------------------------
! read Bx, By, Bz
allocate(Bfield(0:2, 0:nxm1, 0:nym1, 0:nzm1))
allocate(Bfield_tmp(0:nxm1, 0:nym1, 0:nzm1, 0:2))
open(unit=8, file='b3d.bin', Access="stream")
read(8) Bfield_tmp
!$OMP PARALLEL DO  PRIVATE(s), schedule(DYNAMIC)
do s=0,2
	Bfield(s,:,:,:)=Bfield_tmp(:,:,:,s) 
enddo
!$OMP END PARALLEL DO
 close(8)
deallocate(Bfield_tmp)
!----------------------------------------------------------------------------
!allocate(gradUVBfield(0:2,0:2, 0:nxm1, 0:nym1, 0:nzm1))
!!$OMP PARALLEL DO  PRIVATE(k), schedule(DYNAMIC) 
!do k=0,nzm1
!	call grad_UVB_calculate(k)
!enddo	
!!$OMP END PARALLEL DO
!----------------------------------------------------------------------------
if(twistFlag)then
	allocate(curlB(0:2, 0:nxm1, 0:nym1, 0:nzm1))
	!$OMP PARALLEL DO  PRIVATE(k), schedule(DYNAMIC) 
	do k=0,nzm1
		call curlB_calculate(k)
	enddo	
	!$OMP END PARALLEL DO
endif
!----------------------------------------------------------------------------
end subroutine initialize
