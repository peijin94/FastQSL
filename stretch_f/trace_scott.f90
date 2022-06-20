!for the second way of calculating grad unit vec B, please uncomment the subroutine of interpolate_grad_unit_vec_B2 or interpolate_grad_unit_vec_B3, 
!and change line of "call interpolate_grad_unit_vec_B(vp, unit_vec_bp, grad_unit_vec_B)" to 
!"call interpolate_grad_unit_vec_B2(vp, unit_vec_bp, grad_unit_vec_B)" or "call interpolate_grad_unit_vec_B3(vp, unit_vec_bp, grad_unit_vec_B)"

!subroutine  interpolate_grad_unit_vec_B2(vp, unit_vec_bp, grad_unit_vec_B)
!!!the way of https://github.com/Kai-E-Yang/QSL
!implicit none
!real:: vp(0:2),vp1(0:2),vp2(0:2), unit_vec_bp(0:2), unit_vec_bp1(0:2), unit_vec_bp2(0:2), grad_unit_vec_B(0:2,0:2)
!integer:: i
!!----------------------------------------------------------------------------
!do i=0,2
!	vp1=vp
!	vp2=vp
!	vp1(i)=vp(i)-0.001
!	vp2(i)=vp(i)+0.001
!	call interpolate_unit_vec_B(vp1, unit_vec_bp1)
!	call interpolate_unit_vec_B(vp2, unit_vec_bp2)
!	grad_unit_vec_B(i,0:2)=(unit_vec_bp2-unit_vec_bp1)/0.002
!enddo
!call interpolate_unit_vec_B(vp, unit_vec_bp)
!end subroutine interpolate_grad_unit_vec_B2


!subroutine interpolate_grad_unit_vec_B3(vp, unit_vec_bp, grad_unit_vec_B)
!!the way of https://bitbucket.org/tassev/qsl_squasher/src/hg/
!use trace_common
!use field_common
!implicit none
!real:: w(0:1,0:2), wd(0:1), weight(0:1,0:1,0:1), vp(0:2), bp(0:2), unit_vec_bp(0:2), &
!grad_unit_vec_B(0:2,0:2), unit_vec_B_cell(0:2,0:1,0:1,0:1)
!integer:: round(0:1,0:2), i, j, k
!!----------------------------------------------------------------------------
!round(0,:)=floor(vp)
!w(1,:)=vp-round(0,:)

!do i=0,2
!	if ( .not. (vp(i) .ge. 0.0)) then
!		round(0,i)=0	
!		w(1,i)=0.0		
!	else if ( vp(i) .ge. pmax(i)) then
!		round(0,i)=r0max(i)
!		w(1,i)=1.0
!	endif
!enddo

!round(1,:)=round(0,:)+1
!w(0,:)=1.0-w(1,:)

!forall(i=0:1,j=0:1,k=0:1) weight(i,j,k)=w(i,0)*w(j,1)*w(k,2)
!forall(i=0:2) Bp(i)=sum(weight*Bfield(i, round(:,0), round(:,1), round(:,2)))
!unit_vec_bp=bp/norm2(bp)

!wd(0)=-1.0
!wd(1)= 1.0

!forall(i=0:1,j=0:1,k=0:1) unit_vec_B_cell(:,i,j,k)= &
!Bfield(:, round(i,0), round(j,1), round(k,2))/norm2(Bfield(:, round(i,0), round(j,1), round(k,2)))

!forall(i=0:1,j=0:1,k=0:1) weight(i,j,k)=wd(i)*w(j,1)*w(k,2)
!forall(i=0:2)  grad_unit_vec_B(0,i)=sum(weight*unit_vec_B_cell(i,:,:,:))

!forall(i=0:1,j=0:1,k=0:1) weight(i,j,k)=w(i,0)*wd(j)*w(k,2)
!forall(i=0:2)  grad_unit_vec_B(1,i)=sum(weight*unit_vec_B_cell(i,:,:,:))

!forall(i=0:1,j=0:1,k=0:1) weight(i,j,k)=w(i,0)*w(j,1)*wd(k)
!forall(i=0:2)  grad_unit_vec_B(2,i)=sum(weight*unit_vec_B_cell(i,:,:,:))

!end subroutine interpolate_grad_unit_vec_B3


subroutine interpolate_grad_unit_vec_B(vp, unit_vec_bp, grad_unit_vec_B)
use trace_common
use field_common
implicit none
real:: weight(0:1,0:1,0:1), vp(0:2), bp(0:2), unit_vec_bp(0:2), grad_unit_vec_B(0:2,0:2), grad_unit_vec_B_cell(0:2,0:2,0:1,0:1,0:1)
integer:: round(0:1,0:2), i, j, k
!----------------------------------------------------------------------------
call round_weight(vp, round, weight)
forall(i=0:2) Bp(i)=sum(weight*Bfield(i, round(:,0), round(:,1), round(:,2)))
unit_vec_bp= bp/norm2(bp)

if (grad3DFlag) then
	forall(i=0:2, j=0:2) grad_unit_vec_B(i,j)=sum(weight*grad_unit_vec_Bfield(i, j, round(:,0), round(:,1), round(:,2)))
else
!outputs are identical as the upper way, but don't requiring grad_unit_vec_Bfield.
!the efficiency is 1/4.79(gfortran) or 1/2.21(ifort) times of the upper
	do k=0,1
	do j=0,1
	do i=0,1
		if (weight(i,j,k) .ne. 0.0) then
			call grad_unit_vec_B_grid(round(i,0), round(j,1), round(k,2), grad_unit_vec_B_cell(:,:,i,j,k))
		else
			!avoid NaN
			grad_unit_vec_B_cell(:,:,i,j,k)=0.0
		endif
	enddo
	enddo
	enddo
	forall(i=0:2, j=0:2) grad_unit_vec_B(i,j)=sum(weight*grad_unit_vec_B_cell(i, j, :, :, :))
endif

end subroutine interpolate_grad_unit_vec_B


subroutine f_scott(vector9, vector9_k)
use trace_common
implicit none
real::  vector9(0:8), vector9_k(0:8), grad_unit_vec_B(0:2,0:2), vp(0:2), unit_vec_bp(0:2)
integer:: i
!----------------------------------------------------------------------------
vp=vector9(0:2)
call interpolate_grad_unit_vec_B(vp, unit_vec_bp, grad_unit_vec_B)
vector9_k(0:2)= unit_vec_bp
forall(i=0:2) vector9_k(3+i)= dot_product(vector9(3:5), grad_unit_vec_B(0:2,i))
forall(i=0:2) vector9_k(6+i)= dot_product(vector9(6:8), grad_unit_vec_B(0:2,i))
end subroutine f_scott


subroutine f_scott_boundary(vector9, vector9_k, b_dim)
implicit none
real::  vector9(0:8), vector9_k(0:8)
integer:: b_dim
!----------------------------------------------------------------------------
call f_scott(vector9, vector9_k)
vector9_k=vector9_k/vector9_k(b_dim)
vector9_k(b_dim)=1.0
end subroutine f_scott_boundary


subroutine RK4_scott(dt, vector9, vector9_1)
implicit none
real:: dt, vector9(0:8), vector9_1(0:8), k1(0:8), k2(0:8), k3(0:8), k4(0:8)
!----------------------------------------------------------------------------
call f_scott(vector9, k1)
call f_scott(vector9+dt*1./3.*k1, k2)
call f_scott(vector9+dt*(-1./3.*k1+k2), k3)
call f_scott(vector9+dt*(k1-k2+k3), k4)
vector9_1=vector9+dt/8.0*(k1+3.*k2+3.*k3+k4)
end subroutine RK4_scott


subroutine RK4_scott_boundary(dt, vector9, vector9_1, b_dim)
implicit none
real:: dt, vector9(0:8),vector9_1(0:8), k1(0:8), k2(0:8), k3(0:8), k4(0:8)
integer:: b_dim
!----------------------------------------------------------------------------
call f_scott_boundary(vector9, k1, b_dim)
call f_scott_boundary(vector9+dt*1./3.*k1, k2, b_dim)
call f_scott_boundary(vector9+dt*(-1./3.*k1+k2), k3, b_dim)
call f_scott_boundary(vector9+dt*(k1-k2+k3), k4, b_dim)
vector9_1=vector9+dt/8.0*(k1+3.*k2+3.*k3+k4)
end subroutine RK4_scott_boundary


subroutine RKF45_scott(dt, vector9, vector9_1)
use trace_common
use rkf45_common
implicit none
real:: vector9(0:8),vector9_1(0:8), k1(0:8), k2(0:8), k3(0:8), k4(0:8), k5(0:8), k6(0:8)
real:: dt, dt0, dt1, error, min_error, dvp(0:2), scale_dt, vp0(0:2), vp1(0:2)
logical:: continue_flag
integer:: rb, rb_index
!----------------------------------------------------------------------------

continue_flag=.true.
vp0=vector9(0:2)
call f_scott(vector9, k1)

do while ( continue_flag ) 

	call f_scott(vector9+dt*a21*k1, k2)
	call f_scott(vector9+dt*(a31*k1+ a32*k2), k3)   
	call f_scott(vector9+dt*(a41*k1+ a42*k2+ a43*k3), k4)   
	call f_scott(vector9+dt*(a51*k1+ a52*k2+ a53*k3+ a54*k4), k5)   
	call f_scott(vector9+dt*(a61*k1+ a62*k2+ a63*k3+ a64*k4+ a65*k5), k6)

	vector9_1 = vector9+dt*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6)
	
	dvp=dt*(ce1*k1(0:2)+ce3*k3(0:2)+ce4*k4(0:2)+ce5*k5(0:2)+ce6*k6(0:2))
	error = norm2(dvp)
	
	continue_flag =.false.
	if (abs(dt) .gt. min_step) then
		vp1=vector9_1(0:2)
		if (.not.(all(pmin<=vp1 .and. vp1<=pmax) ))  then
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
				if (abs(dt) .lt. min_step) dt=sign(min_step,dt)
				continue_flag=.true.
			endif
			
			cycle
		endif
		
		if  (error .gt. tol) then
			dt=dt*0.618
			if (abs(dt) .lt. min_step) dt=sign(min_step,dt)
			continue_flag=.true.
			cycle
		endif
	endif
	
enddo

!(0.618)^-5 \approx 11.0932
min_error=(tol/11.0932)/((100./abs(dt))**5.)
if (error .lt. min_error) then
	if (error .ne. 0.0) dt=sign(100.,dt)
else
	scale_dt=((tol/error)**0.2)*0.618
	dt=dt*scale_dt
endif

if (abs(dt) .lt. min_step) dt=sign(min_step,dt)

end subroutine RKF45_scott


subroutine correct_foot_scott(vector9, vector9_1, sign_dt, rb)
use trace_common
implicit none
real:: dt, dt0, dt1, vp(0:2), vp0(0:2), vp1(0:2)
real:: vector9(0:8), vector9_0(0:8), vector9_1(0:8), vector9_orig(0:8), vector9_1_orig(0:8)
integer:: sign_dt, rb, rb_index, maxsteps2, it, pmin_mark(0:2), pmax_mark(0:2)
!----------------------------------------------------------------------------
vp=vector9(0:2)
if (any((vp .eq. pmin) .or. (vp .eq. pmax))) then
	pmin_mark=0
	pmax_mark=0
	vector9_1=vector9
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
endif



vp1=vector9_1(0:2)
call vp_rboundary(vp1, rb, rb_index)
if (rb .eq. 0) return
vector9_orig  =vector9
vector9_1_orig=vector9_1
vector9_0=vector9

if (rb .ne. 7) then 
	vp0=vector9_0(0:2)
	vp1=vector9_1(0:2)
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
		call RK4_scott(dt, vector9_0, vector9)
		vp=vector9(0:2)
		do while( .not.( all(pmin<=vp .and. vp<=pmax)) .and. (abs(dt) .ge. min_foot_step) )
			dt= dt*0.9
			if (abs(dt) .ge. min_foot_step) call RK4_scott(dt, vector9_0, vector9)
			vp=vector9(0:2)
		enddo
	endif
endif

dt=min_foot_step*sign_dt
vp=vector9(0:2)
it=0
maxsteps2=norm2(vector9_1(0:2)-vector9(0:2))/abs(dt)*2

do while(all(pmin<=vp .and. vp<=pmax))
	call RK4_scott(dt, vector9, vector9_1) 
	it=it+1
	if (it .ge. maxsteps2) then
		vector9_0=vector9_orig
		vector9_1=vector9_1_orig
		exit
	endif	
	vector9_0=vector9
	vector9  =vector9_1
	vp       =vector9(0:2)
end do

vp =vector9_0(0:2)
vp1=vector9_1(0:2)

call vp_rboundary(vp1, rb, rb_index)
if (rb .eq. 7) return

if( mod(rb, 2) .eq. 1) then
	 dt0=pmin(rb_index)-  vp(rb_index)
	 dt1= vp1(rb_index)-pmin(rb_index)
else
	 dt0=pmax(rb_index)-  vp(rb_index)
	 dt1= vp1(rb_index)-pmax(rb_index)
endif

if (abs(dt0+dt1) .le. 0.05*norm2(vp-vp1)) then
	if(abs(dt0+dt1) .ne. 0.0) then
		vector9_1=(vector9_0*dt1+vector9_1*dt0)/(dt0+dt1)
	else
		vector9_1=vector9_0
	endif
else
	call RK4_scott_Boundary(dt0, vector9_0, vector9_1, rb_index)
endif
end subroutine correct_foot_scott


subroutine grad_unit_vec_B_grid(i, j, k, grad_unit_vec_B)
use trace_common
use field_common
implicit none
integer:: i, j, k
real::grad_unit_vec_B(0:2, 0:2)
!----------------------------------------------------------------------------
if (i .eq. 0) then
	grad_unit_vec_B(0,:) = &
	-1.5*Bfield(:,0,j,k)/norm2(Bfield(:,0,j,k))&
	+2.0*Bfield(:,1,j,k)/norm2(Bfield(:,1,j,k))&
	-0.5*Bfield(:,2,j,k)/norm2(Bfield(:,2,j,k))
else if (i .eq. nxm1) then	
	grad_unit_vec_B(0,:) = &
	 1.5*Bfield(:,i  ,j,k)/norm2(Bfield(:,i  ,j,k))&
	-2.0*Bfield(:,i-1,j,k)/norm2(Bfield(:,i-1,j,k))&
	+0.5*Bfield(:,i-2,j,k)/norm2(Bfield(:,i-2,j,k))
else
	grad_unit_vec_B(0,:) = &
	(Bfield(:,i+1,j,k)/Norm2(Bfield(:,i+1,j,k))-Bfield(:,i-1,j,k)/Norm2(Bfield(:,i-1,j,k)))*0.5
endif
!----------------------------------------------------------------------------
if (j .eq. 0) then
	grad_unit_vec_B(1,:) = &
	-1.5*Bfield(:,i,0,k)/norm2(Bfield(:,i,0,k))&
	+2.0*Bfield(:,i,1,k)/norm2(Bfield(:,i,1,k))&
	-0.5*Bfield(:,i,2,k)/norm2(Bfield(:,i,2,k))
else if (j .eq. nym1) then	
	grad_unit_vec_B(1,:) = &
	 1.5*Bfield(:,i,j  ,k)/norm2(Bfield(:,i,j  ,k))&
	-2.0*Bfield(:,i,j-1,k)/norm2(Bfield(:,i,j-1,k))&
	+0.5*Bfield(:,i,j-2,k)/norm2(Bfield(:,i,j-2,k))	
else
	grad_unit_vec_B(1,:) = &
	(Bfield(:,i,j+1,k)/Norm2(Bfield(:,i,j+1,k))-Bfield(:,i,j-1,k)/Norm2(Bfield(:,i,j-1,k)))*0.5
endif
!----------------------------------------------------------------------------
if (k .eq. 0) then
	grad_unit_vec_B(2,:) = &
	-1.5*Bfield(:,i,j,0)/norm2(Bfield(:,i,j,0))&
	+2.0*Bfield(:,i,j,1)/norm2(Bfield(:,i,j,1))&
	-0.5*Bfield(:,i,j,2)/norm2(Bfield(:,i,j,2))
else if (k .eq. nzm1) then
	grad_unit_vec_B(2,:) = &
	 1.5*Bfield(:,i,j,k  )/norm2(Bfield(:,i,j,k  ))&
	-2.0*Bfield(:,i,j,k-1)/norm2(Bfield(:,i,j,k-1))&
	+0.5*Bfield(:,i,j,k-2)/norm2(Bfield(:,i,j,k-2))
else
	grad_unit_vec_B(2,:) = &
	(Bfield(:,i,j,k+1)/Norm2(Bfield(:,i,j,k+1))-Bfield(:,i,j,k-1)/Norm2(Bfield(:,i,j,k-1)))*0.5
endif

END subroutine grad_unit_vec_B_grid


!Scott_2017_ApJ_848_117
subroutine trace_scott(vp0, q0, q_perp0, rs, re, rbs, rbe, line_length, twist0, twistFlag)
use trace_common
implicit none
real :: dt, vp(0:2), vp0(0:2), q0, bp(0:2), b0square, Bn_s, Bn_e, q_perp0
real :: unit_vec_bp(0:2), rs(0:2),re(0:2), line_length, twist0, alpha, alpha0, dL, dL0, dtwist
real :: vector9(0:8), vector9_0(0:8), vector9_1(0:8), vector9_s(0:8), vector9_e(0:8), &
u0(0:2), us(0:2), ue(0:2), v0(0:2), vs(0:2), ve(0:2), b0(0:2), bs(0:2), be(0:2), &
vs1(0:2), ve1(0:2), us1(0:2), ue1(0:2)
integer:: it, rb, rbe, rbs, e_index, s_index, sign_dt, maxdim, index1, index2
logical:: z0flag, twistFlag
!----------------------------------------------------------------------------
twist0     =0.0
line_length=0.0
!----------------------------------------------------------------------------
z0flag= vp0(2) .eq. zmin
call interpolateB(vp0, b0)

b0square=dot_product(b0, b0)

maxdim=sum(maxloc(abs(b0)))-1
index1=mod(maxdim+1,3)
index2=mod(maxdim+2,3)
v0(maxdim)= b0(index1)
v0(index1)=-b0(maxdim)
v0(index2)=0.
v0        =v0/norm2(v0)

u0(0)=dble(b0(1))*v0(2)-dble(b0(2))*v0(1)
u0(1)=dble(b0(2))*v0(0)-dble(b0(0))*v0(2)
u0(2)=dble(b0(0))*v0(1)-dble(b0(1))*v0(0)
u0   =u0/norm2(u0)
!----------------------------------------------------------------------------

do sign_dt=-1, 1, 2	
	vector9(0:2)=vp0
	vector9(3:5)=u0
	vector9(6:8)=v0
	
	if (z0flag) then
		if ( b0(2)*sign_dt .le. 0) then 
			if (sign_dt .eq. -1) then
				vector9_s=vector9
				rbs      =1
			else
				vector9_e=vector9
				rbe      =1
			endif		
			cycle
		endif
	endif
	
	
	it=0
	dL=0.
	vp=vp0
	if (RK4flag) then
		dt=    step*sign_dt
	else
		dt=min_step*sign_dt
	endif
	
	
	do while( all(pmin<=vp .and. vp<=pmax) .and. abs(it) < maxsteps)

		line_length=line_length+dL		
		
		if (RK4flag) then  
		 	call   RK4_scott(dt, vector9, vector9_1) 
		else			
			call RKF45_scott(dt, vector9, vector9_1) 
		endif
		
		dL0=dL
		dL =norm2(vector9_1(0:2)-vector9(0:2))
		
		if (twistflag) then 
			call interpolateAlpha(vp, unit_vec_bp, alpha)
			if (it .ne. 0) then
				dtwist=(alpha0+alpha)/2.*dL0
				twist0=twist0+dtwist
			endif
			
			alpha0=alpha
		endif
				
		it       =it+sign_dt	
		vector9_0=vector9
		vector9  =vector9_1
		vp       =vector9_1(0:2)
	end do
	
	call correct_foot_scott(vector9_0, vector9_1, sign_dt, rb)
	
	dL=norm2(vector9_1(0:2)-vector9_0(0:2))
	line_length=line_length+dL
	
	if (twistflag) then
		vp=vector9_1(0:2)
		call interpolateAlpha(vp, unit_vec_bp, alpha)
		dtwist=(alpha0+alpha)/2.*dL
		twist0=twist0+dtwist
	endif
	
	if (sign_dt .eq. -1) then
		vector9_s=vector9_1
		rbs      =rb
	else
		vector9_e=vector9_1
		rbe      =rb
	endif
enddo


if (twistflag) twist0=twist0/(4.0*pi)

if ((rbs .eq. 0) .or. (rbe .eq. 0) .or. (rbs .eq. 7) .or. (rbe .eq. 7)) then 
	q0=NaN
	q_perp0=NaN
	return 
endif

rs=vector9_s(0:2)
call interpolateB(rs, bs)
s_index=(6-rbs)/2
Bn_s=bs(s_index)
us=vector9_s(3:5)
vs=vector9_s(6:8)

re=vector9_e(0:2)
call interpolateB(re, be)
e_index=(6-rbe)/2
Bn_e=be(e_index)
ue=vector9_e(3:5)
ve=vector9_e(6:8)


us1=us-us(s_index)/bs(s_index)*bs
vs1=vs-vs(s_index)/bs(s_index)*bs
ue1=ue-ue(e_index)/be(e_index)*be
ve1=ve-ve(e_index)/be(e_index)*be


q0      =abs( dot_product(Ue1,Ue1)*dot_product(Vs1,Vs1)  &
        +     dot_product(Us1,Us1)*dot_product(Ve1,Ve1)  &
        - 2.0*dot_product(Ue1,Ve1)*dot_product(Us1,Vs1))/&
     ( b0square / abs(Bn_s*Bn_e))

ue1=ue-dot_product(ue,be)/norm2(be)*(be/norm2(be))
ve1=ve-dot_product(ve,be)/norm2(be)*(be/norm2(be))
us1=us-dot_product(us,bs)/norm2(bs)*(bs/norm2(bs))
vs1=vs-dot_product(vs,bs)/norm2(bs)*(bs/norm2(bs))

q_perp0 =abs( dot_product(Ue1,Ue1)*dot_product(Vs1,Vs1)  &
        +     dot_product(Us1,Us1)*dot_product(Ve1,Ve1)  &
        - 2.0*dot_product(Ue1,Ve1)*dot_product(Us1,Vs1))/&
     ( b0square / (norm2(bs)*norm2(be)))
 
end subroutine trace_scott
