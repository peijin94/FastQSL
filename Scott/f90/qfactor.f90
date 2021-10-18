! ifort -o qfactor.x qfactor.f90 -fopenmp -mcmodel=large -O3  -xHost -ipo
! gfortran -o qfactor.x qfactor.f90 -fopenmp -mcmodel=large -O3
! the efficiency of ifort (version 2021.2.0) is slightly (~1.3 times) faster than gfortran (version 9.3.0)

include  'trace_bline.f90'
include  'trace_scott.f90'

subroutine qfactor0_bridge(i,j)
use qfactor_common
implicit none
real:: bzp, bp(0:2), vp(0:2), rs(0:2), re(0:2), tw, line_length, q0, q_perp0
integer:: i, j, i12, rbs, rbe, s_ind, e_ind
!----------------------------------------------------------------------------

call ij2vp(i,j,vp)

if (scottFlag) then 
	call trace_scott(vp, q0, q_perp0, rs, re, rbs, rbe, line_length, tw, twistFlag)
	q(i,j)=q0
	q_perp(i,j)=q_perp0
else
	call trace_bline(vp, rs, re, rbs, rbe, line_length, tw, twistFlag)
endif


length(i, j)=line_length
if (twistFlag) twistF(i, j)=tw  

call interpolateB(vp, bp)
bzp=bp(2)
  

if( bzp>0.0) then 
	reboundary(i, j)=rbe
	reF(:, i, j)=re	
	if ( (rbe .ne. 0) .and.(rbe .ne. 7)) then 
		call interpolateB(re, bp)
		e_ind=(6-rbe)/2
		bnr(i, j)=abs(DBLE(bzp)/bp(e_ind))
	endif
else 
	reboundary(i, j)=rbs
	reF(:, i, j)=rs
	if ( (rbs .ne. 0) .and.(rbs .ne. 7)) then 	
		call interpolateB(rs, bp)
		s_ind=(6-rbs)/2
		bnr(i, j)=abs(DBLE(bzp)/bp(s_ind))
	endif
endif

if (bzp .eq. 0.0) reboundary(i, j)=0


if(.not.vflag) then
	i12=i+j*q1
	if (mod(i12+1, 10*nbridges*q1) .eq. 0) call show_time(float(i12+1)/(q1tq2)*100)
endif

END subroutine qfactor0_bridge


subroutine qfactor0_calculate(i,j)
use qfactor_common
implicit none
integer :: i12, i, j, rb, rbs, rbe, m_ind, index1, index2, sign_index1, sign_index2
logical:: bkey1,bkey2, bkey11, bkey12, bkey21, bkey22, margin_flag 
real:: nxx,nxy,nyx,nyy, q0, vp(0:2), q_perp0
real:: rs(0:2), re(0:2), line_length, bnr0, twist
!----------------------------------------------------------------------------
! calculate the Q-factor	
rb=reboundary(i,j)

if ((rb .eq. 0).or.(rb .eq. 7) .or. (bnr(i,j) .eq. 0.0)) then
	q(i,j)=NaN
	return
endif

m_ind=(6-rb)/2
index1=mod(m_ind+1,3)
index2=mod(m_ind+2,3)


margin_flag= (i .eq. 0) .or. (j .eq. 0) .or. (i .eq. q1-1) .or.(j .eq. q2-1)

bkey1= (.not. margin_flag) .and. ( rb .eq. reboundary(i+1,j)).and.( rb .eq. reboundary(i-1,j)) 
bkey2= (.not. margin_flag) .and. ( rb .eq. reboundary(i,j+1)).and.( rb .eq. reboundary(i,j-1))

bkey11=(i+2 .le. q1-1) .and. ( rb .eq. reboundary(i+1,j)) .and. ( rb .eq. reboundary(i+2,j))
bkey12=(i-2 .ge.    0) .and. ( rb .eq. reboundary(i-1,j)) .and. ( rb .eq. reboundary(i-2,j)) .and. (.not. bkey11)
bkey21=(j+2 .le. q2-1) .and. ( rb .eq. reboundary(i,j+1)) .and. ( rb .eq. reboundary(i,j+2))
bkey22=(j-2 .ge.    0) .and. ( rb .eq. reboundary(i,j-1)) .and. ( rb .eq. reboundary(i,j-2)) .and. (.not. bkey21)

if (bkey11) sign_index1= 1
if (bkey12) sign_index1=-1
if (bkey21) sign_index2= 1
if (bkey22) sign_index2=-1

if (bkey1) then
	nxx=reF(index1, i+1,j)-reF(index1, i-1,j)
	nyx=reF(index2, i+1,j)-reF(index2, i-1,j)
else if (bkey11 .or. bkey12) then
	nxx=(-3*reF(index1, i,j)+4*reF(index1, i+sign_index1,j)-reF(index1, i+2*sign_index1,j))*sign_index1	
	nyx=(-3*reF(index2, i,j)+4*reF(index2, i+sign_index1,j)-reF(index2, i+2*sign_index1,j))*sign_index1
endif


if (bkey2) then
	nxy=reF(index1, i,j+1)-reF(index1, i,j-1)
	nyy=reF(index2, i,j+1)-reF(index2, i,j-1)
else if (bkey21 .or.bkey22) then
	nxy=(-3*reF(index1, i,j)+4*reF(index1, i,j+sign_index2)-reF(index1, i,j+2*sign_index2))*sign_index2	
	nyy=(-3*reF(index2, i,j)+4*reF(index2, i,j+sign_index2)-reF(index2, i,j+2*sign_index2))*sign_index2
endif

if ((bkey1 .or. bkey11 .or. bkey12) .and.&
    (bkey2 .or. bkey21 .or. bkey22)) then
	q(i,j)=(nxx*nxx+nxy*nxy+nyx*nyx+nyy*nyy) * factor**2.0 / 4.0 / bnr(i,j)
else	

 	vp=[float(i)*delta+xreg(0), float(j)*delta+yreg(0), 0.0]
	call trace_scott(vp, q0, q_perp0, rs, re, rbs, rbe, line_length, twist, .false.)
	q(i,j)=q0
endif
end subroutine qfactor0_calculate

subroutine qfactor0()
use qfactor_common
implicit none
integer :: i,j
!----------------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE(i,j), schedule(DYNAMIC) 
	DO j= 0, q2m1
	DO i= 0, q1m1 
		call  qfactor0_bridge(i,j)
	enddo	
	enddo
!$OMP END PARALLEL DO


if (scottFlag) return

!$OMP PARALLEL DO PRIVATE(i,j), schedule(DYNAMIC) 
	DO j= 0, q2m1
	DO i= 0, q1m1 
		call  qfactor0_calculate(i,j)
	enddo	
	enddo
!$OMP END PARALLEL DO

end subroutine qfactor0


subroutine qcs_bridge(i,j)
use qfactor_common
use trace_common
implicit none
real:: bp(0:2), vp(0:2), br(0:2), rs(0:2), re(0:2), &
	tw, bn, be, bs, bnt_tmp, line_length, &
	vpa1(0:2), rsa1(0:2), rea1(0:2), vpa2(0:2), rsa2(0:2), rea2(0:2), &
	vpb1(0:2), rsb1(0:2), reb1(0:2), vpb2(0:2), rsb2(0:2), reb2(0:2), &
	nxx, nxy, nyx, nyy, sxx, sxy, syx, syy, exx, exy, eyx, eyy, delta1, q0, q_perp0	
integer:: i, j, i12, rbs, rbe,&
	rbsa1, rbea1, rbsa2, rbea2, rbsb1, rbeb1, rbsb2, rbeb2, &
	maxdim, index1, index2, pmin_mark(0:2), pmax_mark(0:2)
logical::surface_flag, edge_flag, bkey
integer::s_ind,s_index1,s_index2,e_ind,e_index1,e_index2, index_a,index_b
!----------------------------------------------------------------------------

call ij2vp(i,j,vp)

if (scottFlag) then 
	call trace_scott(vp, q0, q_perp0, rs, re, rbs, rbe, line_length, tw, twistFlag)
	q(i,j)=q0
	q_perp(i,j)=q_perp0
else
	call trace_bline(vp, rs, re, rbs, rbe, line_length, tw, twistFlag)
endif

rsF(:, i, j)=rs
reF(:, i, j)=re
rsboundary( i, j)=rbs
reboundary( i, j)=rbe

length(i, j)=line_length
if (twistFlag) twistF(i, j)=tw


if(.not.vflag) then
	i12=i+j*q1
	if (mod(i12+1, 10*nbridges*q1) .eq. 0) call show_time(float(i12+1)/(q1tq2)*100)
endif


if(scottFlag) return

if ((rbe .eq. 0) .or.(rbs .eq. 0) .or. (rbe .eq. 7) .or.(rbs .eq. 7))  return
e_ind=(6-rbe)/2
s_ind=(6-rbs)/2
call interpolateB(re, br)
be=br(e_ind)
call interpolateB(rs, br)
bs=br(s_ind)
!----------------------------------------------------------------------------

! deal with field lines touching the cut plane
pmin_mark=0
pmax_mark=0

where(vp .lt. (pmin+delta) ) pmin_mark=1
where(vp .gt. (pmax-delta) ) pmax_mark=1

select case( sum(pmin_mark+pmax_mark))
	case(0)
		surface_flag=.false.
		edge_flag=.false.
	case(1)
		surface_flag=.true.		
		edge_flag=.false.
		maxdim=sum(maxloc(pmin_mark+pmax_mark))-1
	case(2:3)
		surface_flag=.true.
		edge_flag=.true.
		tflag(i, j)=1! unreal, but for preventing bug
		q(i,j)=nan
end select	

call interpolateB(vp, bp)

if (csflag) then 
	bn=sum(DBLE(bp)*ev3)
else
	bn=bp(CS_ND)
endif

bnr(i, j)=abs(DBLE(bs)*be/(DBLE(bn)**2.))

if (((abs(bn)*5. .le. norm2(bp)) .or. surface_flag) .and. ( .not.  edge_flag )) then
!use the plane quasi-perp to the field line
	tflag(i, j)=1
	
	if (.not. surface_flag) maxdim=sum(maxloc(abs(bp)))-1
	
	index_a=mod(maxdim+1,3)
	index_b=mod(maxdim+2,3)
	delta1=delta/1.
	bn=bp(maxdim)
	bnt_tmp=abs(bs*be/(dble(bn)**2.))	
	
	vpa1=vp
	vpa1(index_a)=vp(index_a)+delta1
	call trace_bline(vpa1, rsa1, rea1, rbsa1, rbea1, line_length, tw, .false.)
	vpa2=vp
	vpa2(index_a)=vp(index_a)-delta1
	call trace_bline(vpa2, rsa2, rea2, rbsa2, rbea2, line_length, tw, .false.)
	vpb1=vp
	vpb1(index_b)=vp(index_b)+delta1
	call trace_bline(vpb1, rsb1, reb1, rbsb1, rbeb1, line_length, tw, .false.)
	vpb2=vp
	vpb2(index_b)=vp(index_b)-delta1
	call trace_bline(vpb2, rsb2, reb2, rbsb2, rbeb2, line_length, tw, .false.)

	bkey= ( rbs .eq. rbsa1).and.( rbs .eq. rbsa2) .and. &
	      ( rbs .eq. rbsb1).and.( rbs .eq. rbsb2) .and. &
	      ( rbe .eq. rbea1).and.( rbe .eq. rbea2) .and. &
	      ( rbe .eq. rbeb1).and.( rbe .eq. rbeb2)

	s_index1=mod(s_ind+1,3)
	s_index2=mod(s_ind+2,3)

	e_index1=mod(e_ind+1,3)
	e_index2=mod(e_ind+2,3)

	if ( bkey ) then
		sxx = rsa1(s_index1)-rsa2(s_index1)
		syx = rsa1(s_index2)-rsa2(s_index2)
		exx = rea1(e_index1)-rea2(e_index1)
		eyx = rea1(e_index2)-rea2(e_index2)
		sxy = rsb1(s_index1)-rsb2(s_index1)
		syy = rsb1(s_index2)-rsb2(s_index2)
		exy = reb1(e_index1)-reb2(e_index1)
		eyy = reb1(e_index2)-reb2(e_index2)
	
		nxx =  exx*syy - exy*syx
		nxy = -exx*sxy + exy*sxx
		nyx =  eyx*syy - eyy*syx
		nyy = -eyx*sxy + eyy*sxx	
		qtmp(i,j) = (nxx*nxx + nxy*nxy + nyx*nyx + nyy*nyy) * bnt_tmp /(delta1**4.)/ 16.
	else
		call trace_scott(vp, q0, q_perp0, rs, re, rbs, rbe, line_length, tw, .false.)
		qtmp(i,j)=q0
	endif

endif

END subroutine qcs_bridge


subroutine qcs_calculate(i,j)
use qfactor_common
implicit none
integer:: i12, i, j, rbs, rbe
real:: sxx, sxy, syx, syy, exx, exy, eyx, eyy, nxx, nxy, nyx, nyy, vp(0:2), q0, q_perp0
integer:: s_ind, s_index1, s_index2, sign_s_index1, sign_s_index2, e_ind, e_index1, e_index2, sign_e_index1, sign_e_index2
logical:: bkeys1, bkeys2, bkeys11, bkeys12, bkeys21, bkeys22, bkeye1,bkeye2, bkeye11, bkeye12, bkeye21, bkeye22
logical:: margin_flag1, margin_flag2
real:: rs(0:2), re(0:2), line_length, twist
!----------------------------------------------------------------------------

rbs=rsboundary(i,j)
rbe=reboundary(i,j)

if((rbs .eq. 0) .or. (rbe .eq. 0).or. (rbs .eq. 7) .or. (rbe .eq. 7)) then 
	q(i,j)=NaN
	return
endif

if (tflag(i,j) .eq. 1) then 
	q(i,j)=qtmp(i, j)
	return
endif

margin_flag1= (i .eq. 0) .or. (i .eq. q1-1) 
margin_flag2= (j .eq. 0) .or. (j .eq. q2-1)

bkeys1 = (.not. margin_flag1) .and. ( rbs .eq. rsboundary(i+1,j)).and.( rbs .eq. rsboundary(i-1,j)) 
bkeys2 = (.not. margin_flag2) .and. ( rbs .eq. rsboundary(i,j+1)).and.( rbs .eq. rsboundary(i,j-1))
bkeye1 = (.not. margin_flag1) .and. ( rbe .eq. reboundary(i+1,j)).and.( rbe .eq. reboundary(i-1,j)) 
bkeye2 = (.not. margin_flag2) .and. ( rbe .eq. reboundary(i,j+1)).and.( rbe .eq. reboundary(i,j-1))


bkeye11=(i+2 .le. q1-1) .and. ( rbe .eq. reboundary(i+1,j)) .and. ( rbe .eq. reboundary(i+2,j))
bkeye12=(i-2 .ge.    0) .and. ( rbe .eq. reboundary(i-1,j)) .and. ( rbe .eq. reboundary(i-2,j)) .and. (.not. bkeye11)
bkeye21=(j+2 .le. q2-1) .and. ( rbe .eq. reboundary(i,j+1)) .and. ( rbe .eq. reboundary(i,j+2))
bkeye22=(j-2 .ge.    0) .and. ( rbe .eq. reboundary(i,j-1)) .and. ( rbe .eq. reboundary(i,j-2)) .and. (.not. bkeye21)



bkeys11=(i+2 .le. q1-1) .and. ( rbs .eq. rsboundary(i+1,j)) .and. ( rbs .eq. rsboundary(i+2,j))
bkeys12=(i-2 .ge.    0) .and. ( rbs .eq. rsboundary(i-1,j)) .and. ( rbs .eq. rsboundary(i-2,j)) .and. (.not. bkeys11)
bkeys21=(j+2 .le. q2-1) .and. ( rbs .eq. rsboundary(i,j+1)) .and. ( rbs .eq. rsboundary(i,j+2))
bkeys22=(j-2 .ge.    0) .and. ( rbs .eq. rsboundary(i,j-1)) .and. ( rbs .eq. rsboundary(i,j-2)) .and. (.not. bkeys21)

s_ind=(6-rbs)/2
s_index1=mod(s_ind+1,3)
s_index2=mod(s_ind+2,3)

e_ind=(6-rbe)/2
e_index1=mod(e_ind+1,3)
e_index2=mod(e_ind+2,3)


if (bkeys11) sign_s_index1= 1
if (bkeys12) sign_s_index1=-1
if (bkeys21) sign_s_index2= 1
if (bkeys22) sign_s_index2=-1

if (bkeye11) sign_e_index1= 1
if (bkeye12) sign_e_index1=-1
if (bkeye21) sign_e_index2= 1
if (bkeye22) sign_e_index2=-1


if (bkeys1) then
	sxx = rsF(s_index1, i+1, j)-rsF(s_index1, i-1, j)
	syx = rsF(s_index2, i+1, j)-rsF(s_index2, i-1, j)
else if ((bkeys11) .or. (bkeys12)) then 
	sxx = (-3*rsF(s_index1, i, j)+4*rsF(s_index1, i+sign_s_index1, j)-rsF(s_index1, i+2*sign_s_index1, j))*sign_s_index1
	syx = (-3*rsF(s_index2, i, j)+4*rsF(s_index2, i+sign_s_index1, j)-rsF(s_index2, i+2*sign_s_index1, j))*sign_s_index1
endif

if (bkeys2) then
	sxy = rsF(s_index1, i, j+1)-rsF(s_index1, i, j-1)
	syy = rsF(s_index2, i, j+1)-rsF(s_index2, i, j-1)
else if ((bkeys21) .or. (bkeys22)) then 
	sxy = (-3*rsF(s_index1, i, j)+4*rsF(s_index1, i, j+sign_s_index2)-rsF(s_index1, i, j+2*sign_s_index2))*sign_s_index2
	syy = (-3*rsF(s_index2, i, j)+4*rsF(s_index2, i, j+sign_s_index2)-rsF(s_index2, i, j+2*sign_s_index2))*sign_s_index2
endif


if (bkeye1) then
	exx = reF(e_index1, i+1, j)-reF(e_index1, i-1, j)
	eyx = reF(e_index2, i+1, j)-reF(e_index2, i-1, j)
else if ((bkeye11) .or. (bkeye12)) then 
	exx = (-3*reF(e_index1, i, j)+4*reF(e_index1, i+sign_e_index1, j)-reF(e_index1, i+2*sign_e_index1, j))*sign_e_index1
	eyx = (-3*reF(e_index2, i, j)+4*reF(e_index2, i+sign_e_index1, j)-reF(e_index2, i+2*sign_e_index1, j))*sign_e_index1
endif

if (bkeye2) then
	exy = reF(e_index1, i, j+1)-reF(e_index1, i, j-1)
	eyy = reF(e_index2, i, j+1)-reF(e_index2, i, j-1)
else if ((bkeye21) .or. (bkeye22)) then 
	exy = (-3*reF(e_index1, i, j)+4*reF(e_index1, i, j+sign_e_index2)-reF(e_index1, i, j+2*sign_e_index2))*sign_e_index2
	eyy = (-3*reF(e_index2, i, j)+4*reF(e_index2, i, j+sign_e_index2)-reF(e_index2, i, j+2*sign_e_index2))*sign_e_index2
endif



if ((bkeys1 .or.  bkeys11 .or. bkeys12) .and. &
    (bkeys2 .or.  bkeys21 .or. bkeys22) .and. &
    (bkeye1 .or.  bkeye11 .or. bkeye12) .and. &
    (bkeye2 .or.  bkeye21 .or. bkeye22) ) then
    
	nxx =  exx*syy - exy*syx
	nxy = -exx*sxy + exy*sxx
	nyx =  eyx*syy - eyy*syx
	nyy = -eyx*sxy + eyy*sxx

	q(i,j) = (nxx*nxx + nxy*nxy + nyx*nyx + nyy*nyy) * bnr(i,j) / (delta**4.)/16.
else
	
	call ij2vp(i,j,vp)
	call trace_scott(vp, q0, q_perp0, rs, re, rbs, rbe, line_length, twist, .false.)
	q(i,j)=q0
endif

end subroutine qcs_calculate


subroutine ij2vp(i, j, vp)
use qfactor_common
implicit none
integer:: i, j
real:: vp(0:2)
!----------------------------------------------------------------------------
if (csflag) then 
	vp=point0+dble(i*delta)*ev1+dble(j*delta)*ev2
else
	select case(CS_ND)
		case(0) 				
			vp=[cut_coordinate, float(i)*delta+yreg(0), float(j)*delta+zreg(0)]
		case(1) 
			vp=[float(i)*delta+xreg(0), cut_coordinate, float(j)*delta+zreg(0)]
		case(2) 
			vp=[float(i)*delta+xreg(0), float(j)*delta+yreg(0), cut_coordinate]
	end select	
endif
end subroutine ij2vp

subroutine qcs()
use qfactor_common
implicit none
integer :: i,j
!----------------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE(i,j), schedule(DYNAMIC) 
	DO j= 0, q2m1
	DO i= 0, q1m1 
		call  qcs_bridge(i,j)
	enddo	
	enddo
!$OMP END PARALLEL DO

if (scottFlag) return

!$OMP PARALLEL DO PRIVATE(i,j), schedule(DYNAMIC) 
	DO j= 0, q2m1
	DO i= 0, q1m1 
		call  qcs_calculate(i,j)
	enddo	
	enddo
!$OMP END PARALLEL DO

end subroutine qcs


subroutine show_time(percent)
real::percent
character(len=30):: time_str(3)
integer:: times(8)
!----------------------------------------------------------------------------
call date_and_time( time_str(1) , time_str(2),  time_str(3), times )		 
print 600, percent, times(5), times(6), times(7)
600 format( '         ', F6.2, '%        ' ,I2.2, ':', I2.2, ':', I2.2)
end subroutine show_time


program qfactor
use qfactor_common
use trace_common
use field_common
implicit none
integer :: k, i, j
real :: bzp, bp(0:2), vp(0:2)
integer(1),allocatable::rboundary3d(:, :, :)
real,allocatable::q3d(:,:,:), q_perp3d(:,:,:), twist3d(:,:,:)
!----------------------------------------------------------------------------

call initialize()

print*, '  _____________________________________'
print*, '        schedule         time'
call show_time(0.0)
!----------------------------------------------------------------------------
if (q0flag) then
!z=0 	
	call qfactor0()
	call show_time(100.0)
	
	open(unit=8, file='qfactor0.bin', access='stream')
	write(8) q, reboundary(0:q1-1, 0:q2-1), length, Bnr
	close(8)
	
	if (twistFlag) then
		open(unit=8, file='twist.bin', access='stream')
		write(8) twistF
		close(8)
	endif
	
	if(scottFlag) then	
		open(unit=8, file='q_perp.bin', access='stream')
		write(8) q_perp
		close(8)		
	endif
	
	open(unit=8, file='reF.bin', access='stream')
	write(8)reF
	close(8)
else
	allocate(rsboundary(-2:q1+1, -2:q2+1))
	allocate(rsF(0:2, 0:q1-1, 0:q2-1))
	allocate(tflag(0:q1-1, 0:q2-1))
!----------------------------------------------------------------------------
	if(vflag)then
	!q3d
		allocate(rboundary3d(0:qx-1, 0:qy-1, 0:qz-1))		
		allocate(q3d(0:qx-1, 0:qy-1, 0:qz-1))
		if (scottFlag) allocate(q_perp3d(0:qx-1, 0:qy-1, 0:qz-1))
		if (twistFlag) allocate(twist3d(0:qx-1,0:qy-1, 0:qz-1))

		call qfactor0()
		call show_time(1.0/qz*100)
		q3d(:, :, 0)=q
		if (scottFlag) q_perp3d(:, :, 0)=q_perp
		if (twistFlag)  twist3d(:, :, 0)=twistF
		
		do j=0, qy-1
		do i=0, qx-1		
			vp=[xreg(0)+delta*i, yreg(0)+delta*j, 0.0]
			call interpolateB(vp, bp)
			bzp=bp(2)   
			if( bzp>0.0) then 
				rboundary3d(i, j, 0)=1+8*reboundary(i, j)
			else 
				rboundary3d(i, j, 0)=8+reboundary(i, j)
			endif				
		enddo
		enddo		

		
		do k=1,qz-1			
			cut_coordinate=zreg(0)+k*delta
			tflag=0
			call qcs()						
			q3d(:, :, k)=q
			if (scottFlag) q_perp3d(:, :, k)=q_perp
			rboundary3d(:, :, k)=rsboundary(0:qx-1,0:qy-1)+8*reboundary(0:qx-1,0:qy-1)
			if (twistFlag) twist3d(:, :, k)=twistF
			call show_time(float(k+1)/qz*100)
		enddo
		
		open(unit=8, file='q3d.bin', access='stream')
		write(8) q3d, rboundary3d
		close(8)
		deallocate(q3d, rboundary3d)
		
		if(scottFlag) then	
			open(unit=8, file='q_perp3d.bin', access='stream')
			write(8) q_perp3d
			close(8)		
		endif
		
		if (twistFlag) then
			open(unit=8, file='twist3d.bin', access='stream')
			write(8) twist3d
			close(8)
			deallocate(twist3d)
		endif
!----------------------------------------------------------------------------	
	else
	!qcs
		tflag=0
		call qcs()
		call show_time(100.0)
		
		open(unit=8, file='qcs.bin', access='stream')
		write(8) q, rsboundary(0:q1-1, 0:q2-1), reboundary(0:q1-1, 0:q2-1), length
	 	close(8)
	 	
		if (twistFlag) then
			open(unit=8, file='twist.bin', access='stream')
			write(8) twistF
			close(8)
		endif
	 	
		if(scottFlag) then	
			open(unit=8, file='q_perp.bin', access='stream')
			write(8) q_perp
			close(8)		
		endif
	
		open(unit=8, file='reF.bin', access='stream')	
		write(8) reF
		close(8)
		open(unit=8, file='rsF.bin', access='stream')	
		write(8) rsF
		close(8)
	endif
	deallocate(rsF, rsboundary, tflag)
endif
!----------------------------------------------------------------------------
deallocate(gradUVBfield, Bfield, q, qtmp, reF, bnr, reboundary)
if (twistFlag) deallocate(twistF, curlB)
if (scottFlag) deallocate(q_perp)
end program qfactor
