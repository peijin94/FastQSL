module qfactor_common
	integer:: qx, qy, qz, nbridges, q1, q2, q1m1, q2m1, q1tq2, Normal_index
	integer(1), allocatable:: rsboundary(:,:), reboundary(:,:), rboundary_tmp(:, :), sign2d(:,:)
	real:: xreg(0:1), yreg(0:1), zreg(0:1), cut_coordinate, delta, delta1, &
	point0(0:2), point1(0:2), point2(0:2), ev1(0:2), ev2(0:2), ev3(0:2)
	real, allocatable:: rsF(:, :, :), reF(:, :, :), bnr(:, :), &
	q_perp(:,:), q(:, :), qtmp(:, :), qtmp2(:, :), length(:, :), twist(:, :)
	logical:: twistFlag, vflag, q0flag, cflag, csflag, scottFlag
	logical, allocatable:: tangent_Flag(:, :)
end module qfactor_common


module trace_common
	integer:: nx, ny, nz, nxm1, nym1, nzm1, r0max(0:2), mmaxsteps, maxsteps
	real:: xmax, ymax, zmax, xmin, ymin, zmin, pmin(0:2), pmax(0:2), step, & 
	min_step, min_foot_step, tol, min_incline, NaN
	real(8), parameter:: pi=3.141592653589793D0
	logical:: RK4flag, grad3DFlag, curlB3DFlag
end module trace_common


module field_common
	real, allocatable:: Bfield(:, :, :, :), CurlB(:,:, :, :), &
	grad_unit_vec_Bfield(:, :, :, :, :)
end module field_common


module rkf45_common
	real:: c2,c3,c4,c5,c6, &
	a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65, &
	b1,b3,b4,b5,b6, &
	ce1,ce3,ce4,ce5,ce6;
end module rkf45_common


subroutine round_weight(vp, round, weight)
use trace_common
implicit none
real:: w(0:1,0:2), weight(0:1,0:1,0:1), vp(0:2)
integer:: round(0:1,0:2), i, j, k
!----------------------------------------------------------------------------
round(0,:)=floor(vp)
w(1,:)=vp-round(0,:)

do i=0,2
	if ( .not. (vp(i) .ge. 0.0)) then 
	! this way can avoid the crash of vp(i) .eq. NaN (caused by B=0), comparing to vp(i) .lt. 0.0
		round(0,i)=0
		w(1,i)=0.0		
	else if ( vp(i) .ge. pmax(i)) then
		round(0,i)=r0max(i)
		w(1,i)=1.0	
	endif
enddo 
 
round(1,:)=round(0,:)+1
w(0,:)=1.0-w(1,:)
forall(i=0:1,j=0:1,k=0:1) weight(i,j,k)=w(i,0)*w(j,1)*w(k,2)

end subroutine round_weight


!! trilinear interpolation
subroutine interpolateB(vp, bp)
use trace_common
use field_common
implicit none
real:: vp(0:2), bp(0:2), weight(0:1,0:1,0:1), Bcell(0:2,0:1,0:1,0:1)
integer:: round(0:1,0:2), i, j, k, s
!----------------------------------------------------------------------------

call round_weight(vp, round, weight)

! this line makes the code faster if compiled by ifort (gfortran will skip this line), 
! maybe due to raising them to the cache of CPU
Bcell=Bfield(:, round(:,0), round(:,1), round(:,2))

! this way is faster than the later one if compiled by gfortran
! the efficiencies of two ways are same for ifort
bp=0.0
do s=0,2
do k=0,1
do j=0,1
do i=0,1
	bp(s)=bp(s)+weight(i,j,k)*Bfield(s, round(i,0), round(j,1), round(k,2))
enddo
enddo
enddo
enddo

!forall(i=0:2) Bp(i)=sum(weight*Bfield(i, round(:,0), round(:,1), round(:,2)))
end subroutine interpolateB


subroutine interpolate_unit_vec_B(vp, unit_vec_bp)
!Unit vector B
implicit none
real::vp(0:2), unit_vec_bp(0:2), bp(0:2)
!---------------------------------------------------------------------------
call interpolateB(vp, bp)
unit_vec_bp= bp/norm2(bp)
end subroutine interpolate_unit_vec_B


subroutine interpolateAlpha(vp, unit_vec_bp, alpha)
use trace_common
use field_common
implicit none
real:: weight(0:1,0:1,0:1), vp(0:2), bp(0:2), unit_vec_bp(0:2), CurlBp(0:2), alpha, curlB_cell(0:2,0:1,0:1,0:1)
integer:: round(0:1,0:2), i, j, k
!----------------------------------------------------------------------------
call round_weight(vp, round, weight)
forall(i=0:2) Bp(i)=sum(weight*Bfield(i, round(:,0), round(:,1), round(:,2)))

if (curlB3DFlag) then
	forall(i=0:2) CurlBp(i)=sum(weight*curlB(i, round(:,0), round(:,1), round(:,2)))
else
	do k=0,1
	do j=0,1
	do i=0,1
		if (weight(i,j,k) .ne. 0.0) then
			call curlB_grid(round(i,0), round(j,1), round(k,2), curlB_cell(:,i,j,k))
		else
			!avoid NaN 
			curlB_cell(:,i,j,k)=0.0
		endif
	enddo
	enddo
	enddo
	forall(i=0:2) CurlBp(i)=sum(weight*curlB_cell(i,:,:,:))
endif

unit_vec_bp= bp/norm2(bp)
alpha=sum(DBLE(curlbp)*bp)/(sum(DBLE(bp)*bp))
end subroutine interpolateAlpha


!d vp/dt=Bp/Bt
subroutine RK4(dt, vp0, vp1, alpha, alphaFlag)
use field_common
implicit none
real ::  dt, alpha, vp0(0:2), vp1(0:2), &
k1(0:2), k2(0:2), k3(0:2), k4(0:2)
logical:: alphaFlag
!----------------------------------------------------------------------------
if (alphaFlag) then 
	call interpolateAlpha(vp0, k1, alpha)
else
	call interpolate_unit_vec_B(vp0, k1)
endif

call interpolate_unit_vec_B(vp0+dt*1./3.*k1, k2)
call interpolate_unit_vec_B(vp0+dt*(-1./3.*k1+k2), k3)
call interpolate_unit_vec_B(vp0+dt*(k1-k2+k3), k4)

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
call interpolateB(vp0+dt*1./3.*k1, bp)
k2=bp/bp(b_dim)
k2(b_dim)=1.0
call interpolateB(vp0+dt*(-1./3.*k1+k2), bp)
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
real:: k1(0:2), k2(0:2), k3(0:2), k4(0:2), k5(0:2), k6(0:2), vp0(0:2), vp1(0:2)
real:: dt, dt0, dt1, alpha, error, min_error, tol_this, dvp(0:2), scale_dt
logical:: continue_flag, alphaFlag
integer:: rb, rb_index
!----------------------------------------------------------------------------
continue_flag=.true.

if (alphaFlag) then 
	call interpolateAlpha(vp0, k1, alpha)
else
	call interpolate_unit_vec_B(vp0, k1)
endif

do while ( continue_flag ) 

	call interpolate_unit_vec_B(vp0+dt*a21*k1, k2)
	call interpolate_unit_vec_B(vp0+dt*(a31*k1+ a32*k2),k3)   
	call interpolate_unit_vec_B(vp0+dt*(a41*k1+ a42*k2+ a43*k3), k4)   
	call interpolate_unit_vec_B(vp0+dt*(a51*k1+ a52*k2+ a53*k3+ a54*k4),k5)   
	call interpolate_unit_vec_B(vp0+dt*(a61*k1+ a62*k2+ a63*k3+ a64*k4+ a65*k5),k6)

	vp1 = vp0 + dt*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6)
	dvp = dt*(ce1*k1+ce3*k3+ce4*k4+ce5*k5+ce6*k6)
	error = norm2(dvp)
	continue_flag =.false.
	if (abs(dt) .gt. min_step) then
		
		if (.not.(all(pmin<=vp1 .and. vp1<=pmax)))  then	
			call vp_rboundary(vp1, rb, rb_index)
			
			if (rb .ne. 7) then
				if (abs(dt) .ge. 2.*norm2(vp1-vp0))  then
					dt=dt/2.
				else
					if( mod(rb, 2) .eq. 1) then
					 	dt0=pmin(rb_index)- vp0(rb_index)
					 	dt1= vp1(rb_index)-pmin(rb_index)
					else
					 	dt0=pmax(rb_index)- vp0(rb_index)
					 	dt1= vp1(rb_index)-pmax(rb_index)
					endif			
					dt=dt*dt0/(dt0+dt1)*0.95
				endif
				
				if (abs(dt) .lt. min_step) then	
					dt=sign(min_step, dt)			
					continue_flag=.false.
				else
					continue_flag=.true.
				endif
				
			else
				dt=dt*0.618
				if (abs(dt) .lt. min_step) dt=sign(min_step, dt)
				continue_flag=.true.
			endif
			
			cycle

		endif
		
		if  (error .gt. tol_this) then
			dt=dt*0.618
			if (abs(dt) .lt. min_step) dt=sign(min_step, dt)
			continue_flag=.true.
			cycle
		endif
	endif
	
enddo

!(0.618)^-5 \approx 11.0932
min_error=(tol_this/11.0932)/((100./abs(dt))**5.)
if (error .lt. min_error) then
	if (error .ne. 0.0) dt=sign(100., dt)	
else
	scale_dt=((tol_this/error)**0.2)*0.618
	dt=dt*scale_dt
endif

if (abs(dt) .lt. min_step) dt=sign(min_step, dt)

end subroutine RKF45


subroutine vp_rboundary(vp, rb, rb_index)
use trace_common
implicit none
real:: vp(0:2), xp, yp, zp
integer:: rb, rb_index
!----------------------------------------------------------------------------
if (all(pmin<=vp .and. vp<=pmax)) then
	rb=0
	return
endif

xp=vp(0)
yp=vp(1)
zp=vp(2)
  
if ((xp .ge. xmin) .and. (xp .le. xmax) .and. &
    (yp .ge. ymin) .and. (yp .le. ymax) .and. &
    (zp .lt. zmin)) then 
      rb=1
      rb_index=2
      return
endif

if ((xp .ge. xmin) .and. (xp .le. xmax) .and. &
    (yp .ge. ymin) .and. (yp .le. ymax) .and. &
    (zp .gt. zmax)) then 
      rb=2
      rb_index=2
      return
endif

if ((xp .ge. xmin) .and. (xp .le. xmax) .and. &
    (yp .lt. ymin) .and. &
    (zp .ge. zmin) .and. (zp .le. zmax)) then 
      rb=3
      rb_index=1
      return
endif

if ((xp .ge. xmin) .and. (xp .le. xmax) .and. &
    (yp .gt. ymax) .and. &
    (zp .ge. zmin) .and. (zp .le. zmax)) then 
      rb=4
      rb_index=1
      return
endif

if ((xp .lt. xmin) .and. &
    (yp .ge. ymin) .and. (yp .le. ymax) .and. &
    (zp .ge. zmin) .and. (zp .le. zmax)) then 
      rb=5
      rb_index=0
      return
endif

if ((xp .gt. xmax) .and. &
    (yp .ge. ymin) .and. (yp .le. ymax) .and. &
    (zp .ge. zmin) .and. (zp .le. zmax)) then 
      rb=6
      rb_index=0
      return
endif

rb=7

end subroutine vp_rboundary


subroutine correct_foot(vp, vp1, sign_dt, rb)
use trace_common
implicit none
real:: dt, dt0, dt1, alpha
real:: vp(0:2), vp0(0:2), vp1(0:2), vp_orig(0:2), vp1_orig(0:2), vp_tmp(0:2)
integer:: sign_dt, rb, rb_index, maxsteps2, it, pmin_mark(0:2), pmax_mark(0:2)
!---------------------------------------------------------------------------
if (any((vp .eq. pmin) .or. (vp .eq. pmax))) then
	pmin_mark=0
	pmax_mark=0
	vp1=vp
	where(vp .eq. pmin) pmin_mark=1
	where(vp .eq. pmax) pmax_mark=1
	if (sum(pmin_mark+pmax_mark) .eq. 1) then
		if (pmin_mark(2) .eq. 1) rb=1
		if (pmax_mark(2) .eq. 1) rb=2
		if (pmin_mark(1) .eq. 1) rb=3
		if (pmax_mark(1) .eq. 1) rb=4
		if (pmin_mark(0) .eq. 1) rb=5
		if (pmax_mark(0) .eq. 1) rb=6	
	else 
		rb=7
	endif
	return
endif

call vp_rboundary(vp1, rb, rb_index)
if (rb .eq. 0) return

vp0=vp
vp1_orig=vp1
vp_orig =vp
if (rb .ne. 7) then 
	if( mod(rb, 2) .eq. 1) then
		dt0=pmin(rb_index)- vp0(rb_index)
		dt1= vp1(rb_index)-pmin(rb_index)
	else
		dt0=pmax(rb_index)- vp0(rb_index)
		dt1= vp1(rb_index)-pmin(rb_index)
	endif

	if  (RK4flag) then
		dt=    step*abs(dt0/(dt0+dt1))*sign_dt*0.95
	else
		dt=min_step*abs(dt0/(dt0+dt1))*sign_dt*0.95
	endif
	
	if (abs(dt) .ge. min_foot_step) then
		call RK4(dt, vp0, vp, alpha, .false.)
		do while( .not.( all(pmin<=vp .and. vp<=pmax)) .and. (abs(dt) .ge. min_foot_step) )
			dt= dt*0.9
			if (abs(dt) .ge. min_foot_step)  call RK4(dt, vp0, vp, alpha, .false.)
		enddo
	endif
endif

dt=min_foot_step*sign_dt
vp_tmp=vp
it=0
maxsteps2=norm2(vp1_orig-vp_orig)/abs(dt)*2

do while( all(pmin<=vp_tmp .and. vp_tmp<=pmax))
	call RK4(dt, vp_tmp, vp1, alpha, .false.)
	it=it+1	
	if (it .ge. maxsteps2) then
		vp0=vp_orig
		vp1=vp1_orig
		exit
	endif
	vp0=vp_tmp
	vp_tmp=vp1
end do

call vp_rboundary(vp1, rb, rb_index)
if (rb .eq. 7) return
	
if( mod(rb, 2) .eq. 1) then
	dt0=pmin(rb_index)- vp0(rb_index)
	dt1= vp1(rb_index)-pmin(rb_index)
else
	dt0=pmax(rb_index)- vp0(rb_index)
 	dt1= vp1(rb_index)-pmax(rb_index)
endif


if (abs(dt0+dt1) .le. 0.05*norm2(vp0-vp1)) then 	
	if((dt0+dt1) .ne. 0.0) then 
		vp1=(vp0*dt1+vp1*dt0)/(dt0+dt1) 
	else
		vp1=vp0
	endif
else
	call RK4_Boundary(dt0, vp0, vp1, rb_index)
endif

end subroutine correct_foot


subroutine trace_bline(vp0, rs, re, rbs, rbe, line_length, twist0, twistFlag, incline)
!+
! NAME :
!   trace_bline
!   
!PURPOSE:
!     Uses RK4 or RKF45 solver to integrate a line of force for a 3D magnetic field

!INPUTS:
!     vp0        - the position to start the line at
!     twistFlag  - compute the twist or not
!     incline    - abs(bn/ norm2(bp)) at the cross section, for adjusting step or tol

!OUTPUTS:
!     rs           - pixel coordinates of the start point of field lines
!     re           - pixel coordinates of the end point of field lines
!     rbs,rbe      - categorize field lines based on where they thread the boundary of the field cube
!     line_length  -length of the field line
!     twist0       - twist number of the field line
!-
use trace_common
implicit none
real:: dt, vp0(0:2), vp1(0:2), vp(0:2), vp_tmp(0:2), rs(0:2), re(0:2), bp(0:2), unit_vec_bp(0:2), &
twist0, line_length, alpha, alpha0, dL, dL0, dtwist, tol_this, step_this, incline, incline_this
integer:: it, rb, rbs, rbe, sign_dt
logical:: twistFlag, z0Flag
!----------------------------------------------------------------------------
twist0=0.
line_length=0.
!----------------------------------------------------------------------------
z0flag= vp0(2) .eq. zmin

if (incline .le. min_incline)  then
	incline_this = min_incline
else
	incline_this = incline
endif

if (RK4flag) then 
	step_this=step*incline_this
	if (step_this .le. min_step) step_this=min_step
else
	step_this=min_step
	tol_this=tol*(incline_this**1.5)
endif

do sign_dt=-1,1,2
	if (z0flag) then
		
		call interpolateB(vp0, bp)
		
		if (bp(2)*sign_dt .le. 0) then
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
	
	vp=vp0
	it=0
	dt=step_this*sign_dt
	dL=0.
	do while( all(pmin<=vp .and. vp<=pmax) .and. it < maxsteps)
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
				twist0=twist0+dtwist
			endif
			
			alpha0=alpha
		endif
		
		it=it+1	
		vp_tmp=vp
		
		vp=vp1	
		
		
	end do
	
	call correct_foot(vp_tmp, vp1, sign_dt, rb)
	
	
	dL=norm2(vp1-vp_tmp)	
	line_length=line_length+dL
	if (twistflag) then 
		call interpolateAlpha(vp1, unit_vec_bp, alpha)
		dtwist=(alpha0+alpha)/2.*dL
		twist0=twist0+dtwist
	endif
	
	if (sign_dt .eq. -1) then
		rs =vp1
		rbs=rb
	else
		re =vp1
		rbe=rb
	endif
enddo
if (twistflag) twist0=twist0/(4.0*pi)
END subroutine trace_bline


subroutine curlB_grid(i, j, k, CurlBp)
use trace_common
use field_common
implicit none
integer:: i, j, k
real:: dBxdy, dBxdz, dBydx, dBydz, dBzdx, dBzdy, CurlBp(0:2)
!----------------------------------------------------------------------------
if (i .eq. 0) then
	dBydx =-1.5*Bfield(1,0,j,k)   +2.0*Bfield(1,1,j,k)     -0.5*Bfield(1,2,j,k)
	dBzdx =-1.5*Bfield(2,0,j,k)   +2.0*Bfield(2,1,j,k)     -0.5*Bfield(2,2,j,k)
else if (i .eq. nxm1) then
	dBydx = 1.5*Bfield(1,nxm1,j,k)-2.0*Bfield(1,nxm1-1,j,k)+0.5*Bfield(1,nxm1-2,j,k)
	dBzdx = 1.5*Bfield(2,nxm1,j,k)-2.0*Bfield(2,nxm1-1,j,k)+0.5*Bfield(2,nxm1-2,j,k)
else
	dBydx = (Bfield(1,i+1,j,k)-Bfield(1,i-1,j,k))*0.5
	dBzdx = (Bfield(2,i+1,j,k)-Bfield(2,i-1,j,k))*0.5
endif
!----------------------------------------------------------------------------
if (j .eq. 0) then
	dBxdy =-1.5*Bfield(0,i,0,k)   +2.0*Bfield(0,i,1,k)     -0.5*Bfield(0,i,2,k)
	dBzdy =-1.5*Bfield(2,i,0,k)   +2.0*Bfield(2,i,1,k)     -0.5*Bfield(2,i,2,k)
else if (j .eq. nym1) then
	dBxdy = 1.5*Bfield(0,i,nym1,k)-2.0*Bfield(0,i,nym1-1,k)+0.5*Bfield(0,i,nym1-2,k)
	dBzdy = 1.5*Bfield(2,i,nym1,k)-2.0*Bfield(2,i,nym1-1,k)+0.5*Bfield(2,i,nym1-2,k)
else
	dBxdy = (Bfield(0,i,j+1,k)-Bfield(0,i,j-1,k))*0.5
	dBzdy = (Bfield(2,i,j+1,k)-Bfield(2,i,j-1,k))*0.5
endif
!----------------------------------------------------------------------------
if (k .eq. 0) then
	dBydz =-1.5*Bfield(1,i,j,0)   +2.0*Bfield(1,i,j,1)     -0.5*Bfield(1,i,j,2)
	dBxdz =-1.5*Bfield(0,i,j,0)   +2.0*Bfield(0,i,j,1)     -0.5*Bfield(0,i,j,2)
else if (k .eq. nzm1) then
	dBydz = 1.5*Bfield(1,i,j,nzm1)-2.0*Bfield(1,i,j,nzm1-1)+0.5*Bfield(1,i,j,nzm1-2)
	dBxdz = 1.5*Bfield(0,i,j,nzm1)-2.0*Bfield(0,i,j,nzm1-1)+0.5*Bfield(0,i,j,nzm1-2)
else
	dBydz = (Bfield(1,i,j,k+1)-Bfield(1,i,j,k-1))*0.5
	dBxdz = (Bfield(0,i,j,k+1)-Bfield(0,i,j,k-1))*0.5
endif
!----------------------------------------------------------------------------
 curlBp(0) = dbzdy - dbydz
 curlBp(1) = dbxdz - dbzdx
 curlBp(2) = dbydx - dbxdy

END subroutine curlB_grid


subroutine logical2int(logical_in, interger_out)
implicit none
logical:: logical_in
integer:: interger_out
!----------------------------------------------------------------------------
if (logical_in) then
	interger_out=1
else
	interger_out=0
endif

end subroutine logical2int


subroutine initialize()
use qfactor_common
use trace_common
use field_common
use rkf45_common
implicit none
integer:: i, j, k, s, twistFlag_int, csFlag_int, RK4flag_int, scottFlag_int, &
q0flag_int, vflag_int, cflag_int, reclen
real, allocatable:: Bfield_tmp(:, :, :, :)
!----------------------------------------------------------------------------
open(unit=8, file='head.txt', status='old')
read(8, *) nx, ny, nz, nbridges, delta, maxsteps, &
           xreg, yreg, zreg, step, tol, &
           twistFlag_int, RK4flag_int, scottFlag_int, csFlag_int
 close(8)
!----------------------------------------------------------------------------
twistFlag= twistFlag_int .eq. 1
  RK4Flag=   RK4flag_int .eq. 1
scottFlag= scottFlag_int .eq. 1
   csFlag=    csFlag_int .eq. 1 

   q0Flag=(zreg(0) .eq. 0) .and. ( zreg(1) .eq. 0) .and. (.not. csflag)    
    vFlag=(xreg(1) .ne. xreg(0)) .and. (yreg(1) .ne. yreg(0)) .and. (zreg(1) .ne. zreg(0)) .and. (.not. csflag)
    cFlag=(.not. vflag) .and. (.not. q0flag)


if (csflag) then
	Normal_index =-1
	point0=[xreg(0),yreg(0),zreg(0)]
	point1=[xreg(1),yreg(1),zreg(0)]
	point2=[xreg(0),yreg(0),zreg(1)]
	
	q1=norm2(point1-point0)/delta+1
	q2=norm2(point2-point0)/delta+1
	qx=0; qy=0; qz=0
	!ev*: elementary vector of the cut plane
	ev1   =(point1-point0)/norm2(point1-point0)
	ev2   =(point2-point0)/norm2(point2-point0)
	ev3(0)=dble(ev1(1))*ev2(2)-dble(ev1(2))*ev2(1)
	ev3(1)=dble(ev1(2))*ev2(0)-dble(ev1(0))*ev2(2)
	ev3(2)=dble(ev1(0))*ev2(1)-dble(ev1(1))*ev2(0)
	ev3=ev3/norm2(ev3)
else
	qx=(xreg(1)-xreg(0))/delta+1
	qy=(yreg(1)-yreg(0))/delta+1
	qz=(zreg(1)-zreg(0))/delta+1
	if(vflag)then
		q1=qx; q2=qy; Normal_index=2
	endif
	if (qx .eq. 1) then 
		q1=qy; q2=qz; Normal_index=0; cut_coordinate=xreg(0)
	endif	
	if (qy .eq. 1) then 
		q1=qx; q2=qz; Normal_index=1; cut_coordinate=yreg(0)
	endif
	if (qz .eq. 1) then 
		q1=qx; q2=qy; Normal_index=2; cut_coordinate=zreg(0)
	endif
endif

call logical2int(q0flag, q0flag_int)
call logical2int( cflag,  cflag_int)
call logical2int( vflag,  vflag_int)

! tell IDL these informations
open(unit=8, file='tail.txt', status='replace')
write(8, *) qx, qy, qz, q1, q2
write(8, *) q0flag_int, cflag_int, vflag_int
 close(8)
!----------------------------------------------------------------------------
! for RKF45
if (.not. RK4flag) then
	a21=   1./4.
	a31=   3./32.;   a32=    9./32.
	a41=1932./2197.; a42=-7200./2197.; a43=  7296./2197.
	a51= 439./216.;  a52=-8.;          a53=  3680./513.;   a54= -845./4104.
	a61=  -8./27.;   a62= 2.;          a63= -3544./2565.;  a64= 1859./4104.; a65=-11./40.
	b1 =  16./135.;   b3= 6656./12825.; b4= 28561./56430.;  b5=   -9./50.;    b6=  2./55. 
	ce1=   1./360.;  ce3= -128./4275.; ce4= -2197./75240.; ce5=    1./50.;   ce6=  2./55.
endif
!----------------------------------------------------------------------------
q1tq2=q1*q2
q2m1 =q2-1
q1m1 =q1-1
nxm1 =nx-1
nym1 =ny-1
nzm1 =nz-1
xmin =0
ymin =0
zmin =0
pmin =[xmin,ymin,zmin]
xmax =nxm1
ymax =nym1
zmax =nzm1
pmax =[xmax,ymax,zmax]
r0max=[nxm1,nym1,nzm1]-1
NaN=transfer(2143289344, 1.0)

!maxsteps    =nint(4*(nx+ny+nz)/step)*10
mmaxsteps    =-maxsteps
min_incline  =0.05
min_step     =1.0/4.0
min_foot_step=1.0/8.0
delta1	=delta/2.0
!----------------------------------------------------------------------------
! read Bx, By, Bz

allocate(Bfield(0:2, 0:nxm1, 0:nym1, 0:nzm1))
allocate(Bfield_tmp(0:nxm1, 0:nym1, 0:nzm1, 0:2))
open(unit=8, file='b3d.bin', Access='stream', status='old')
read(8) Bfield_tmp
 close(8)
 
 
inquire(iolength=reclen) 1.0 ! if reclen .eq. 4, compiled by gfortran; if reclen .eq. 1, compiled by ifort;
if (reclen .eq. 4) then 
	! the setting of 4 threads achieves the best performance here for gfortran
	CALL OMP_set_num_threads(minval([4,nbridges])) 
else
	CALL OMP_set_num_threads(nbridges)
endif

!$OMP PARALLEL DO  PRIVATE(k), schedule(DYNAMIC)
do k=0, nzm1
	forall(s=0:2) Bfield(s,:,:,k)=Bfield_tmp(:,:,k,s)
enddo
!$OMP END PARALLEL DO
deallocate(Bfield_tmp)

if (reclen .eq. 4) CALL OMP_set_num_threads(nbridges)
!----------------------------------------------------------------------------
!Even for method 3 of Pariat_2012_A&A_541_A78, we deal most tangent part with Scott_2017_ApJ_848_117.
!But this part only takes a small proportion, it is unnecessary to allocate grad_unit_vec_Bfield if only deal with one cross section
grad3DFlag= scottFlag .or. vFlag 

if (grad3DFlag) then
	allocate(grad_unit_vec_Bfield(0:2,0:2, 0:nxm1, 0:nym1, 0:nzm1))
	!$OMP PARALLEL DO  PRIVATE(i, j, k), schedule(DYNAMIC) 
	do k=0, nzm1
	do j=0, nym1
	do i=0, nxm1
		call grad_unit_vec_B_grid(i, j, k, grad_unit_vec_Bfield(:,:,i,j,k))
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO
endif
!----------------------------------------------------------------------------
 curlB3DFlag= twistFlag .and. vFlag ! only allocate this array for 3D case
 
if(curlB3DFlag)then
	allocate(curlB(0:2, 0:nxm1, 0:nym1, 0:nzm1))

	!$OMP PARALLEL DO  PRIVATE(i, j, k), schedule(DYNAMIC) 
	do k=0, nzm1
	do j=0, nym1
	do i=0, nxm1
		call curlB_grid(i, j, k, curlB(:,i,j,k))
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO	
endif
!----------------------------------------------------------------------------
allocate(     q(0:q1m1, 0:q2m1))
allocate(  qtmp(0:q1m1, 0:q2m1))
allocate(   bnr(0:q1m1, 0:q2m1))
allocate(length(0:q1m1, 0:q2m1))
allocate(reF(0:2, 0:q1m1, 0:q2m1))
allocate(reboundary(-2:q1+1, -2:q2+1))
allocate(rboundary_tmp(0:q1m1, 0:q2m1))

if (twistFlag) allocate(twist(0:q1m1, 0:q2m1))
if (scottFlag) allocate(q_perp(0:q1m1, 0:q2m1))

if (vflag .or. cflag) then
	allocate(rsboundary(-2:q1+1, -2:q2+1))
	allocate(rsF(0:2, 0:q1m1, 0:q2m1))
	allocate(tangent_Flag(-1:q1, -1:q2))
endif

end subroutine initialize
