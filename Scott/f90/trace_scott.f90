subroutine modify_foot_scott(vector9, vector9_1, sign_dt, rb)
use trace_common
implicit none
real :: dt, dt1, vp(0:2), vp1(0:2)
real :: vector9(0:8), vector9_0(0:8), vector9_1(0:8), vector9_orig(0:8), vector9_1_orig(0:8)
integer:: sign_dt, rb, rb_index, maxsteps2, it
!----------------------------------------------------------------------------
vp1=vector9_1(0:2)
call vp_rboundary(vp1, rb, rb_index)
if (rb .eq. 0) return

vector9_orig  =vector9
vector9_1_orig=vector9_1
dt=1.0/16* sign_dt
vp=vector9(0:2)
it=0
maxsteps2=norm2(vector9_1(0:2)-vector9(0:2))/abs(dt)*2
do while( minval(vp-pmin)>=0 .and. maxval(vp-pmax)<=0)
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
vp=vector9_0(0:2)
vp1=vector9_1(0:2)

call vp_rboundary(vp1, rb, rb_index)
if (rb .eq. 7) return

if( mod(rb, 2) .eq. 1) then
	 dt =pmin(rb_index)-  vp(rb_index)
	 dt1= vp1(rb_index)-pmin(rb_index)
else
	 dt =pmax(rb_index)-  vp(rb_index)
	 dt1= vp1(rb_index)-pmax(rb_index)
endif

if((dt+dt1) .ne. 0.0) vector9_1=(vector9_0*dt1+vector9_1*dt)/(dt+dt1)
end subroutine modify_foot_scott


!subroutine graduvb(vp, grad_UVB)
!implicit none
!real:: vp(0:2),vp1(0:2),vp2(0:2), uvbp1(0:2), uvbp2(0:2), grad_UVB(0:2,0:2)
!integer:: i
!!----------------------------------------------------------------------------
!do i=0,2
!	vp1=vp
!	vp2=vp
!	vp1(i)=vp(i)-0.001
!	vp2(i)=vp(i)+0.001
!	call interpolateUVB(vp1, uvbp1)
!	call interpolateUVB(vp2, uvbp2)
!	grad_UVB(i,0:2)=(uvbp2-uvbp1)/0.002
!enddo
!end subroutine graduvb


subroutine interpolate_graduvb(vp, uvbp, grad_UVB)
use trace_common
use field_common
implicit none
real::w(0:1,0:2), weigh(0:1,0:1,0:1), vp(0:2), bp(0:2), uvbp(0:2), grad_UVB(0:2,0:2)
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
forall(i=0:2,j=0:2) grad_UVB(i,j)=sum(gradUVBfield(i, j, round(:,0) ,round(:,1), round(:,2))*weigh)
uvbp=bp/norm2(bp)

end subroutine interpolate_graduvb


subroutine f_scott(vector9, vector9_k)
implicit none
real ::  vector9(0:8),vector9_k(0:8), grad_UVB(0:2,0:2), vp(0:2), uvbp(0:2)
integer:: i
!----------------------------------------------------------------------------
vp=vector9(0:2)
!call interpolateUVB(vp, uvbp)
!call graduvb(vp, grad_UVB)
call interpolate_graduvb(vp, uvbp, grad_UVB)
vector9_k(0:2)= uvbp
forall(i=0:2) vector9_k(3+i)= dot_product(vector9(3:5),grad_UVB(0:2,i))
forall(i=0:2) vector9_k(6+i)= dot_product(vector9(6:8),grad_UVB(0:2,i))
end subroutine f_scott


subroutine RK4_scott(dt, vector9, vector9_1)
implicit none
real :: dt, vector9(0:8),vector9_1(0:8), k1(0:8), k2(0:8), k3(0:8), k4(0:8)
!----------------------------------------------------------------------------
call f_scott(vector9, k1)
call f_scott(vector9+dt*1./3.*k1, k2)
call f_scott(vector9+dt*(-1./3.*k1+k2), k3)
call f_scott(vector9+dt*(k1-k2+k3), k4)
vector9_1=vector9+dt/8.0*(k1+3.*k2+3.*k3+k4)
end subroutine RK4_scott


subroutine RKF45_scott(dt, vector9, vector9_1)
use trace_common
use rkf45_common
implicit none
real:: vector9(0:8),vector9_1(0:8), k1(0:8), k2(0:8), k3(0:8), k4(0:8), k5(0:8), k6(0:8),vp0(0:2),vp1(0:2)
real:: dt, dt0, dt1, step_error, min_step_error
real:: dvp(0:2), scale_dt
logical:: continue_flag
integer:: rb, rb_index
!----------------------------------------------------------------------------

continue_flag=.true.
vp0=vector9(0:2)
call f_scott(vector9, k1)

do while ( continue_flag ) 

	call f_scott(vector9+dt*a21*k1, k2)
	call f_scott(vector9+dt*(a31*k1+ a32*k2),k3)   
	call f_scott(vector9+dt*(a41*k1+ a42*k2+ a43*k3), k4)   
	call f_scott(vector9+dt*(a51*k1+ a52*k2+ a53*k3+ a54*k4),k5)   
	call f_scott(vector9+dt*(a61*k1+ a62*k2+ a63*k3+ a64*k4+ a65*k5),k6)

	vector9_1 = vector9 + dt*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6)
	
	dvp=dt*(ce1*k1(0:2)+ce3*k3(0:2)+ce4*k4(0:2)+ce5*k5(0:2)+ce6*k6(0:2))
	step_error = norm2(dvp)
	
	continue_flag =.false.
	if (abs(dt) .gt. step) then
		vp1=vector9_1(0:2)
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
		
		if  (step_error .gt. tol) then
			dt=dt*0.618
			if (abs(dt)  .lt. step) dt=sign(step,dt)
			continue_flag=.true.
			cycle
		endif
	endif
	
enddo

!11.0932=(0.618)^-5
min_step_error=(tol/11.0932)/((100./abs(dt))**5.)
if (step_error .lt. min_step_error) then
	dt=sign(100.,dt)
else
	scale_dt=((tol/step_error)**0.2)*0.618
	dt=dt*scale_dt
endif

if (abs(dt)  .lt. step) dt=sign(step,dt)

end subroutine RKF45_scott


subroutine grad_UVB_calculate(k)
use trace_common
use field_common
implicit none
integer:: i,j,k
!----------------------------------------------------------------------------
do j=0,nym1
do i=0,nxm1
!----------------------------------------------------------------------------
	if (i .eq. 0) then
		gradUVBfield(0,:,i,j,k) = &
		-3.0*Bfield(:,0,j,k)/ norm2(Bfield(:,0,j,k))&
		+4.0*Bfield(:,1,j,k)/ norm2(Bfield(:,1,j,k))&
		-    Bfield(:,2,j,k)/ norm2(Bfield(:,2,j,k))
	else if (i .eq. nxm1) then	
		gradUVBfield(0,:,i,j,k) = &
		 3.0*Bfield(:,nxm1  ,j,k)/ norm2(Bfield(:,nxm1  ,j,k))&
		-4.0*Bfield(:,nxm1-1,j,k)/ norm2(Bfield(:,nxm1-1,j,k))&
		+    Bfield(:,nxm1-2,j,k)/ norm2(Bfield(:,nxm1-2,j,k))
	else
		gradUVBfield(0,:,i,j,k) = Bfield(:,i+1,j,k)/Norm2(Bfield(:,i+1,j,k))-Bfield(:,i-1,j,k) /Norm2(Bfield(:,i-1,j,k))
	endif
!----------------------------------------------------------------------------
	if (j .eq. 0) then
		gradUVBfield(1,:,i,j,k) = &
		-3.0*Bfield(:,i,0,k)/ norm2(Bfield(:,i,0,k))&
		+4.0*Bfield(:,i,1,k)/ norm2(Bfield(:,i,1,k))&
		-    Bfield(:,i,2,k)/ norm2(Bfield(:,i,2,k))
	else if (j .eq. nym1) then	
		gradUVBfield(1,:,i,j,k) = &
		 3.0*Bfield(:,i,nym1  ,k)/ norm2(Bfield(:,i,nym1  ,k))&
		-4.0*Bfield(:,i,nym1-1,k)/ norm2(Bfield(:,i,nym1-1,k))&
		+    Bfield(:,i,nym1-2,k)/ norm2(Bfield(:,i,nym1-2,k))	
	else
		gradUVBfield(1,:,i,j,k) = Bfield(:,i,j+1,k)/Norm2(Bfield(:,i,j+1,k))-Bfield(:,i,j-1,k) /Norm2(Bfield(:,i,j-1,k))
	endif
!----------------------------------------------------------------------------
	if (k .eq. 0) then
		gradUVBfield(2,:,i,j,k) = &
		-3.0*Bfield(:,i,j,0)/ norm2(Bfield(:,i,j,0))&
		+4.0*Bfield(:,i,j,1)/ norm2(Bfield(:,i,j,1))&
		-    Bfield(:,i,j,2)/ norm2(Bfield(:,i,j,2))
	else if (k .eq. nzm1) then	
		gradUVBfield(2,:,i,j,k) = &
		 3.0*Bfield(:,i,j,nym1  )/ norm2(Bfield(:,i,j,nym1  ))&
		-4.0*Bfield(:,i,j,nym1-1)/ norm2(Bfield(:,i,j,nym1-1))&
		+    Bfield(:,i,j,nym1-2)/ norm2(Bfield(:,i,j,nym1-2))
	else
		gradUVBfield(2,:,i,j,k) = Bfield(:,i,j,k+1)/Norm2(Bfield(:,i,j,k+1))-Bfield(:,i,j,k-1) /Norm2(Bfield(:,i,j,k-1))
	endif
!----------------------------------------------------------------------------

	gradUVBfield(:,:,i,j,k)=gradUVBfield(:,:,i,j,k)/2.
enddo
enddo

END subroutine grad_UVB_calculate


!Scott_2017_ApJ_848_117
subroutine trace_scott(vp0, q0, q_perp0, rs, re, rbs, rbe, line_length, twist, twistFlag)
use trace_common
implicit none
real :: dt, vp(0:2), vp0(0:2), q0, bp(0:2), Nb0, Nbs, Nbe, q_perp0
real :: incline, uvbp(0:2), rs(0:2),re(0:2), line_length, twist, alpha, alpha0, dL, dL0, dtwist
integer:: it, rb, rbe, rbs, e_ind, s_ind, sign_dt
real:: vector9(0:8), vector9_0(0:8), vector9_1(0:8), vector9_s(0:8), vector9_e(0:8),&
u0(0:2), us(0:2), ue(0:2), v0(0:2), vs(0:2), ve(0:2), b0(0:2), bs(0:2), be(0:2), &
vs1(0:2), ve1(0:2), us1(0:2), ue1(0:2)
logical:: z0flag, twistFlag
!----------------------------------------------------------------------------
q0=0.0
q_perp0=0.0
twist=0.
line_length=0.
!----------------------------------------------------------------------------
z0flag= vp0(2) .eq. 0.0
call interpolateB(vp0, b0)

Nb0=Norm2(b0)

v0=[b0(1),-b0(0),0.]
v0=v0 - dot_product(v0,b0)*b0
v0=v0/norm2(v0)

u0(0)=dble(b0(1))*v0(2)-dble(b0(2))*v0(1)
u0(1)=dble(b0(2))*v0(0)-dble(b0(0))*v0(2)
u0(2)=dble(b0(0))*v0(1)-dble(b0(1))*v0(0)
u0=u0/norm2(u0)


do sign_dt=-1,1,2	
	vector9(0:2)=vp0
	vector9(3:5)=u0
	vector9(6:8)=v0
	
	if (z0flag) then
		if ( b0(2)*sign_dt .le. 0) then 
			if (sign_dt .eq. -1) then
				vector9_s=vector9
				rbs=1
			else
				vector9_e=vector9
				rbe=1
			endif		
			cycle
		endif
	endif
	
	
	it=0
	dt=step*sign_dt
	vp= vp0
	dL=0.
	do while( minval(vp-pmin)>=0 .and. maxval(vp-pmax)<=0 .and. abs(it) < maxsteps)

		line_length=line_length+dL
		
		
		if (RK4flag) then  
		 	call RK4_scott(dt, vector9, vector9_1) 
		else			
			call RKF45_scott(dt, vector9, vector9_1) 
		endif
		
		dL0=dL
		dL=norm2(vector9_1(0:2)-vector9(0:2))
		
		if (twistflag) then 
			call interpolateAlpha(vp, uvbp, alpha, twistflag)
			if (it .ne. 0) then
				dtwist=(alpha0+alpha)/2.*dL0
				twist=twist+dtwist
			endif
			
			alpha0=alpha
		endif
				
		it=it+sign_dt	
		vector9_0=vector9
		vector9=vector9_1
		vp=vector9_1(0:2)
	end do
	
	call modify_foot_scott(vector9_0, vector9_1, sign_dt, rb)
	
	dL=norm2(vector9_1(0:2)-vector9_0(0:2))	
	line_length=line_length+dL
	
	if (twistflag) then
		vp=vector9_1(0:2)
		call interpolateAlpha(vp, uvbp, alpha, twistflag)
		dtwist=(alpha0+alpha)/2.*dL
		twist=twist+dtwist
	endif
	
	if (sign_dt .eq. -1) then
		vector9_s=vector9_1
		rbs=rb
	else
		vector9_e=vector9_1
		rbe=rb
	endif
enddo

if (twistflag) twist=twist/(4.0*pi)




if ((rbs .eq.0) .or. (rbe .eq.0) .or. (rbs .eq.7) .or. (rbe .eq.7)) return 

rs=vector9_s(0:2)
call interpolateB(rs, bs)
s_ind=(6-rbs)/2
Nbs=bs(s_ind)
us=vector9_s(3:5)
vs=vector9_s(6:8)

re=vector9_e(0:2)
call interpolateB(re, be)
e_ind=(6-rbe)/2
Nbe=be(e_ind)
ue=vector9_e(3:5)
ve=vector9_e(6:8)


us1=us-us(s_ind)/bs(s_ind)*bs
vs1=vs-vs(s_ind)/bs(s_ind)*bs
ue1=ue-ue(e_ind)/be(e_ind)*be
ve1=ve-ve(e_ind)/be(e_ind)*be


q0 =abs( dot_product(Ue1,Ue1)*dot_product(Vs1,Vs1)  &
   +     dot_product(Us1,Us1)*dot_product(Ve1,Ve1)  &
   - 2.0*dot_product(Ue1,Ve1)*dot_product(Us1,Vs1))/&
     (abs( (Nb0*Nb0) / (Nbs*Nbe)))

ue1=ue-dot_product(ue,be)/norm2(be)*(be/norm2(be))
ve1=ve-dot_product(ve,be)/norm2(be)*(be/norm2(be))
us1=us-dot_product(us,bs)/norm2(bs)*(bs/norm2(bs))
vs1=vs-dot_product(vs,bs)/norm2(bs)*(bs/norm2(bs))

q_perp0 =abs( dot_product(Ue1,Ue1)*dot_product(Vs1,Vs1)  &
        +     dot_product(Us1,Us1)*dot_product(Ve1,Ve1)  &
        - 2.0*dot_product(Ue1,Ve1)*dot_product(Us1,Vs1))/&
     ( (Nb0*Nb0) / (norm2(bs)*norm2(be)))
     
end subroutine trace_scott
