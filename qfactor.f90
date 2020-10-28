include  'trace_bline.f90'

subroutine qfactor0_bridge(i12)
use qfactor_common
implicit none
real:: bzp, bp(0:2), vp(0:2), vp1(0:2), br(0:2), rs(0:2), re(0:2), tw, &
	x1r(0:2), x2r(0:2), y1r(0:2), y2r(0:2), line_length
integer:: i1, i2, i12, rbs, rbe, times(8)
 character(len=30):: time_str(3)
logical::bkey
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
i1=mod(i12, q1+2)-1
i2=i12/(q1+2)-1

tw=0.0
vp=[float(i1)*delta+xreg(0), float(i2)*delta+yreg(0), 0.0]


if (.not. (minval(vp-pmin)>=0 .and. maxval(vp-pmax)<=0 )) then
	reF(:, i1, i2)=NAN
	reboundary( i1, i2)=0
	return   
endif

call trace_bline(vp, rs, re, rbs, rbe, tw, line_length)

call interpolateB(vp, bp)
bzp=bp(2)
   
if (( (i1 .eq. -1) .or. (i2 .eq. -1) .or. (i1 .eq. q1) .or.(i2 .eq. q2))) then

	if( bzp>0.0) then 
		reboundary(i1, i2)=rbe
		reF(:, i1, i2)=re
	else 
		reboundary(i1, i2)=rbs
		reF(:, i1, i2)=rs
	endif	
	if (bzp .eq. 0.0) 	reboundary(i1, i2)=0
	
	return
endif

length(i1, i2)=line_length
if (twistFlag) twistF(i1, i2)=tw  
   
   
if( bzp>0.0) then 
	reboundary(i1, i2)=rbe
	reF(:, i1, i2)=re
		
	call interpolateB(re, br)
	select case(rbe)
		case(0) 				
				bnr(i1, i2)=NaN 
		case(1:2) 
				bnr(i1, i2)=abs(DBLE(bzp)/br(2))
		case(3:4) 
				bnr(i1, i2)=abs(DBLE(bzp)/br(1))
		case(5:6) 
				bnr(i1, i2)=abs(DBLE(bzp)/br(0))
		end select		
	else 
	reboundary(i1, i2)=rbs
	r+F(:, i1, i2)=rs
		
	call interpolateB(rs, br)
	select case(rbs)
		case(0) 
				bnr(i1, i2)=NaN 
		case(1:2) 
				bnr(i1, i2)=abs(DBLE(bzp)/br(2))
		case(3:4) 
				bnr(i1, i2)=abs(DBLE(bzp)/br(1))
		case(5:6) 
				bnr(i1, i2)=abs(DBLE(bzp)/br(0))
	end select		
endif

if (bzp .eq. 0.0) 	reboundary(i1, i2)=0
	

if(.not.vqflag) then
	if (mod(i12+1, 10*nbridges*q1) .eq. 0) then
		call date_and_time( time_str(1) , time_str(2),  time_str(3), times )		 
		print 400, float(i12+1)/(q1tq2)*100, times(5), times(6), times(7)
		400 format( '         ', F6.2, '%        ' I2.2, ':', I2.2, ':', I2.2)
	endif
endif

END


!;###########################################################
subroutine qfactor0_calculate(i12)
use qfactor_common
implicit none
integer :: i12, i, j, rb
logical:: bkey

! calculate the Q-factor	
i=mod(i12, q1)
j=i12/q1


rb=reboundary(i,j)
bkey= ( rb .eq. reboundary(i+1,j)).and.( rb .eq. reboundary(i-1,j)) &
.and.( rb .eq. reboundary(i,j+1)).and.( rb .eq. reboundary(i,j-1))

if (.NOT. bkey) then 
	q(i,j)=NAN
else
	if (bnr(i,j) .ne. 0.0) then
		select case(rb)
			case(0) 
				q(i,j) =NAN
			case(1:2) 					
				q(i,j)=((reF(0, i+1,j)-reF(0, i-1,j))**2 + (reF(1, i+1,j)-reF(1, i-1,j))**2 +&
				(reF(0, i,j+1)-reF(0, i,j-1))**2 + (reF(1, i,j+1)-reF(1, i,j-1))**2) * factor**2.0 / 4.0 / bnr(i,j)
			case(3:4) 
				q(i,j)=((reF(0, i+1,j)-reF(0, i-1,j))**2 + (reF(2, i+1,j)-reF(2, i-1,j))**2 + &
				(reF(0, i,j+1)-reF(0, i,j-1))**2 + (reF(2, i,j+1)-reF(2, i,j-1))**2) * factor**2.0 / 4.0/ bnr(i,j)   
			case(5:6) 
				q(i,j)= ((reF(2, i+1,j)-reF(2, i-1,j))**2 + (reF(1, i+1,j)-reF(1, i-1,j))**2 + &
				(reF(2, i,j+1)-reF(2, i,j-1))**2 + (reF(1, i,j+1)-reF(1, i,j-1))**2) * factor**2.0 / 4.0 / bnr(i,j)
		end select		
	else 
		q(i,j)=NAN
	endif
endif
     
end
!###########################################################



subroutine qfactor0()
use qfactor_common
implicit none
integer :: i, i_end
 !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

i_end=(q1+2)*(q2+2)-1
!$OMP PARALLEL DO SHARED(i_end), PRIVATE(i), schedule(DYNAMIC) 
	DO i= 0, i_end
		call  qfactor0_bridge(i)
	enddo	
!$OMP END PARALLEL DO

i_end=q1tq2-1
!$OMP PARALLEL DO SHARED(i_end), PRIVATE(i), schedule(DYNAMIC) 
	DO i= 0, i_end
		call  qfactor0_calculate(i)
	enddo
!$OMP END PARALLEL DO
end
!###########################################################

subroutine qcs_bridge(i12)
use qfactor_common
implicit none
 character(len=30):: time_str(3)
real:: bp(0:2), vp(0:2), br(0:2), rs(0:2), re(0:2), &
	tw, bn, be, bs, bnt_tmp, line_length, &
	vpa1(0:2), rsa1(0:2), rea1(0:2), vpa2(0:2), rsa2(0:2), rea2(0:2), &
	vpb1(0:2), rsb1(0:2), reb1(0:2), vpb2(0:2), rsb2(0:2), reb2(0:2), &
	nxx, nxy, nyx, nyy, sxx, sxy, syx, syy, exx, exy, eyx, eyy
integer:: i1, i2, i12, j, rbs, rbe, times(8), &
	rbsa1, rbea1, rbsa2, rbea2, rbsb1, rbeb1, rbsb2, rbeb2, &
	maxdim, index1, index2, pmin_mark(0:2), pmax_mark(0:2)
logical::bkeys, bkeye, surface_flag, edge_flag
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
i1=mod(i12, q1+2)-1
i2=i12/(q1+2)-1

tw=0.0

if (csflag) then 
	vp=point0+dble(i1*delta)*ev1+dble(i2*delta)*ev2
else
	select case(CS_ND)
		case(0) 				
			vp=[cut_coordinate, float(i1)*delta+yreg(0), float(i2)*delta+zreg(0)]
		case(1) 
			vp=[float(i1)*delta+xreg(0), cut_coordinate, float(i2)*delta+zreg(0)]
		case(2) 
			vp=[float(i1)*delta+xreg(0), float(i2)*delta+yreg(0), cut_coordinate]
	end select	
endif


if (.not. (minval(vp-pmin)>=0 .and. maxval(vp-pmax)<=0 )) then
	reF(:, i1, i2)=NAN
	reboundary(i1, i2)=0
	rsF(:, i1, i2)=NAN
	rsboundary(i1, i2)=0
	return      
endif





call trace_bline(vp, rs, re, rbs, rbe, tw, line_length)
rsF(:, i1, i2)=rs
reF(:, i1, i2)=re
rsboundary( i1, i2)=rbs
reboundary( i1, i2)=rbe





if (( (i1 .eq. -1) .or. (i2 .eq. -1) .or. (i1 .eq. q1) .or.(i2 .eq. q2))) return

length(i1, i2)=line_length
if (twistFlag) twistF(i1, i2)=tw


call interpolateB(re, br)
select case(rbe)
	case(0) 				
		be=NAN  
	case(1:2) 
		be=br(2)
	case(3:4) 
		be=br(1)
	case(5:6) 
		be=br(0)
end select		
		
call interpolateB(rs, br)
select case(rbs)
	case(0) 				
		bs=NAN  
	case(1:2) 
		bs=br(2)
	case(3:4) 
		bs=br(1)
	case(5:6) 
		bs=br(0)
end select  
 !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
! deal with field lines touching the cut plane
pmin_mark=0
pmax_mark=0
where(vp .le. (pmin+delta/2) ) pmin_mark=1
where(vp .ge. (pmax-delta/2) ) pmax_mark=1

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
		tflag(i1, i2)=1! no real, but for prevent bug
		q(i1,i2)=nan
end select	



call interpolateB(vp, bp)

if (csflag) then 
	bn=sum(DBLE(bp)*ev3)
else
	bn=bp(CS_ND)
endif

bnr(i1, i2)=abs(DBLE(bs)*be/(DBLE(bn)**2.))


if (((abs(bn)*50 .le. sqrt(sum(DBLE(bp)*bp))) .or. surface_flag) .and. ( .not.  edge_flag )) then	
	! thse points are calculated with Method 3 in Pariat and Demoulin (2012) 
	
	
	tflag(i1, i2)=1
	if (.not. surface_flag) maxdim=sum(maxloc(abs(bp)))-1
	
	bn=bp(maxdim)
	bnt_tmp=abs(bs*be/(dble(bn)**2.))
	index1=mod(maxdim+1,3)
	index2=mod(index1+1,3)
	
	vpa1=vp
	vpa1(index1)=vp(index1)+delta
	call trace_bline(vpa1, rsa1, rea1, rbsa1, rbea1, tw, line_length)
	vpa2=vp
	vpa2(index1)=vp(index1)-delta
	call trace_bline(vpa2, rsa2, rea2, rbsa2, rbea2, tw, line_length)
	vpb1=vp
	vpb1(index2)=vp(index2)+delta
	call trace_bline(vpb1, rsb1, reb1, rbsb1, rbeb1, tw, line_length)
	vpb2=vp
	vpb2(index2)=vp(index2)-delta
	call trace_bline(vpb2, rsb2, reb2, rbsb2, rbeb2, tw, line_length)

     bkeys = ( rbs .eq. rbsa1) .and.( rbs .eq. rbsa2) .and.( rbs .eq. rbsb1) .and. ( rbs .eq. rbsb2) 
     bkeye = ( rbe .eq. rbea1) .and.( rbe .eq. rbea2) .and.( rbe .eq. rbeb1) .and. ( rbe .eq. rbeb2) 

	if ((.not. bkeys) .or. (.not. bkeye)) then
		 qtmp(i1, i2)=NAN
	else
         ! launch point(x+,y+)     
		select case(rbs)
			case(0)
				 qtmp(i1, i2)=NAN
			case(1:2) 
				sxx = rsa1(0)-rsa2(0)
				sxy = rsb1(0)-rsb2(0)
				syx = rsa1(1)-rsa2(1)
				syy = rsb1(1)-rsb2(1)
			case(3:4) 
				sxx = rsa1(2)-rsa2(2)
				sxy = rsb1(2)-rsb2(2)
				syx = rsa1(0)-rsa2(0)
				syy = rsb1(0)-rsb2(0)		
			case(5:6) 		
				sxx = rsa1(1)-rsa2(1)
				sxy = rsb1(1)-rsb2(1)
				syx = rsa1(2)-rsa2(2)
				syy = rsb1(2)-rsb2(2)
			end select				
		
					 ! target point (x-,y-)
		select case(rbe)
			case(0) 				
				qtmp(i1, i2)=NAN
			case(1:2) 
				exx = rea1(0)-rea2(0)
				exy = reb1(0)-reb2(0)
				eyx = rea1(1)-rea2(1)
				eyy = reb1(1)-reb2(1)
			case(3:4) 
				exx = rea1(2)-rea2(2)
				exy = reb1(2)-reb2(2)
				eyx = rea1(0)-rea2(0)
				eyy = reb1(0)-reb2(0)		
			case(5:6) 		
				exx = rea1(1)-rea2(1)
				exy = reb1(1)-reb2(1)
				eyx = rea1(2)-rea2(2)
				eyy = reb1(2)-reb2(2)				
		end select		
		         
        ! em=((exx,exy),(eyx,eyy))
        ! nm=em##sm
      
		     nxx =  exx*syy - exy*syx
		     nxy = -exx*sxy + exy*sxx
		     nyx =  eyx*syy - eyy*syx
		     nyy = -eyx*sxy + eyy*sxx
		     
		     ! Eq.(21-22) in Pariat & Demoulin (2012)
		     qtmp(i1, i2)= (nxx*nxx + nxy*nxy + nyx*nyx + nyy*nyy) * bnt_tmp /(delta**4.)/16. 
		     
	endif
	
endif
  
if(.not.vqflag) then
	if (mod(i12+1, 10*nbridges*q1) .eq. 0) then
		call date_and_time( time_str(1) , time_str(2),  time_str(3), times )		 
		print 500, float(i12+1)/(q1tq2)*100, times(5), times(6), times(7)
		500 format( '         ', F6.2, '%        ' I2.2, ':', I2.2, ':', I2.2)   
	endif
endif

END
!###########################################################
subroutine qcs_calculate(i12)
use qfactor_common
implicit none
integer :: i12, i, j, rbs, rbe
logical:: bkeys, bkeye
real:: sxx, sxy, syx, syy, exx, exy, eyx, eyy, nxx, nxy, nyx, nyy

!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
! calculate the Q-factor
i=mod(i12, q1)
j=i12/q1

	if (tflag(i,j) .eq. 1) then 
		q(i,j)=qtmp(i, j)
		return
	endif

     
	rbs=rsboundary(i,j)
	bkeys = ( rbs .eq. rsboundary(i+1,j)).and.( rbs .eq. rsboundary(i-1,j)).and.&
	( rbs .eq. rsboundary(i,j+1)).and.( rbs .eq. rsboundary(i,j-1)) 
		   
	rbe=reboundary(i,j)
	bkeye = ( rbe .eq. reboundary(i+1,j)).and.( rbe .eq. reboundary(i-1,j)).and.&
	( rbe .eq. reboundary(i,j+1)).and.( rbe .eq. reboundary(i,j-1))
                      
	if ((.not. bkeys) .or. (.not. bkeye)) then
		q(i,j)=NAN
	else
         ! launch point(x+,y+)     
         select case(rbs)
			case(0)
			case(1:2) 
				sxx = rsF(0, i+1, j)-rsF(0, i-1, j)
				sxy = rsF(0, i, j+1)-rsF(0, i, j-1)
				syx = rsF(1, i+1, j)-rsF(1, i-1, j)
				syy = rsF(1, i, j+1)-rsF(1, i, j-1)
			case(3:4) 
				sxx = rsF(2, i+1, j)-rsF(2, i-1, j)
				sxy = rsF(2, i, j+1)-rsF(2, i, j-1)
				syx = rsF(0, i+1, j)-rsF(0, i-1, j)
				syy = rsF(0, i, j+1)-rsF(0, i, j-1)
			case(5:6) 
				sxx = rsF(1, i+1, j)-rsF(1, i-1, j)
				sxy = rsF(1, i, j+1)-rsF(1, i, j-1)
				syx = rsF(2, i+1, j)-rsF(2, i-1, j)
				syy = rsF(2, i, j+1)-rsF(2, i, j-1)
		end select		
		
		! target point (x-,y-)
		select case(rbe)
			case(0)
			case(1:2) 
				exx = reF(0, i+1, j)-reF(0, i-1, j)
				exy = reF(0, i, j+1)-reF(0, i, j-1)
				eyx = reF(1, i+1, j)-reF(1, i-1, j)
				eyy = reF(1, i, j+1)-reF(1, i, j-1)
			case(3:4) 
				exx = reF(2, i+1, j)-reF(2, i-1, j)
				exy = reF(2, i, j+1)-reF(2, i, j-1)
				eyx = reF(0, i+1, j)-reF(0, i-1, j)
				eyy = reF(0, i, j+1)-reF(0, i, j-1)
			case(5:6) 
				exx = reF(1, i+1, j)-reF(1, i-1, j)
				exy = reF(1, i, j+1)-reF(1, i, j-1)
				eyx = reF(2, i+1, j)-reF(2, i-1, j)
				eyy = reF(2, i, j+1)-reF(2, i, j-1)
		end select
		         
        ! em=((exx,exy),(eyx,eyy))
        ! nm=em##sm
      
         nxx =  exx*syy - exy*syx
         nxy = -exx*sxy + exy*sxx
         nyx =  eyx*syy - eyy*syx
         nyy = -eyx*sxy + eyy*sxx
         
         ! Eq.(21-22) in Pariat & Demoulin (2012)
         q(i,j) = (nxx*nxx + nxy*nxy + nyx*nyx + nyy*nyy) * bnr(i,j) * (factor**4.)/16. 
         
         if((rbs .eq. 0) .or. (rbe .eq. 0)) then 
         	q(i,j)=nan
		 endif
		 
    endif    

end

!###########################################################
subroutine qcs()
use qfactor_common
implicit none
integer :: i,i_end
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

i_end=(q1+2)*(q2+2)-1
!$OMP PARALLEL DO SHARED(i_end), PRIVATE(i), schedule(DYNAMIC) 
	DO i= 0, i_end
		call  qcs_bridge(i)
	enddo	
!$OMP END PARALLEL DO

i_end=q1tq2-1
!$OMP PARALLEL DO SHARED(i_end), PRIVATE(i), schedule(DYNAMIC) 
	DO i= 0, i_end
		call  qcs_calculate(i)
	enddo	
!$OMP END PARALLEL DO
end

!###########################################################
program qfactor
use qfactor_common
implicit none
integer :: i, i1, i2, rec1, times(8)
integer,allocatable::rboundary3d(:, :, :)
integer,allocatable::rboundary(:, :)
 character(len=30):: binfile, time_str(3)
real,allocatable::q3d(:,:,:), twist3d(:,:,:)
real :: bzp, bp(0:2), vp(0:2)
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
 call initialize()
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
if (q0flag) then
!z=0
	call qfactor0()
	binfile='qfactor0.bin'
	allocate(rboundary(0:q1-1, 0:q2-1))
	rboundary=reboundary(0:q1-1, 0:q2-1)
	open(unit=8, file=binfile	, access="direct",recl=q1*q2)
	write(8, rec=1) q
	write(8, rec=2) rboundary
	deallocate(rboundary)
	write(8, rec=3) length
	write(8, rec=4) Bnr
	if (twistFlag) write(8, rec=5) twistF
	close(8, status='KEEP')

else
	allocate(rsboundary(-1:q1, -1:q2))
	allocate(rsF(0:2, -1:q1, -1:q2))
	allocate(tflag(0:q1-1, 0:q2-1))
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
	if(vqflag)then
	!q3d
	allocate(rboundary3d(0:qx-1, 0:qy-1, 0:qz-1))
		call qfactor0()		
		
		allocate(q3d(0:qx-1, 0:qy-1, 0:qz-1))
		if (twistFlag) then 
			allocate(twist3d(0:qx-1,0:qy-1, 0:qz-1));  twist3d(:, :, 0)=twistF
		endif
		q3d(:, :, 0)=q
		do i1=0, q1
			do i2=0, q2
				vp=[xreg(0)+delta*i1, yreg(0)+delta*i2, 0.0]
				call interpolateB(vp, bp)
				bzp=bp(2)   
				if( bzp>0.0) then 
					rboundary3d(i1, i2, 0)=1+7*reboundary(i1, i2)
				else 
					rboundary3d(i1, i2, 0)=7+reboundary(i1, i2)
				endif				
			enddo
		enddo
		
		call date_and_time( time_str(1) , time_str(2),  time_str(3), times )		 
		print 600, 1.0/(qz)*100, times(5), times(6), times(7)
		600 format( '         ', F6.2, '%        ' I2.2, ':', I2.2, ':', I2.2) 
		
		do i=1,qz-1	
			cut_coordinate=zreg(0)+i*delta
			tflag=0
			call qcs()
						
			q3d(:, :, i)=q
			rboundary3d(:, :, i)=rsboundary(0:q1-1, 0:q2-1)+7*reboundary(0:q1-1, 0:q2-1)
			if (twistFlag) twist3d(:, :, i)=twistF
		
			call date_and_time( time_str(1) , time_str(2),  time_str(3), times ) 
			print 600, float(i+1)/(qz)*100, times(5), times(6), times(7)	  
			
		enddo
	
		binfile='q3d.bin'
		open(unit=8, file=binfile, access="direct",recl=qx*qy)
		do rec1=0, qz-1 
		  write(8,rec=1+rec1) q3d(:, :,rec1)
		  write(8,rec=1+rec1+qz) rboundary3d(:, :,rec1)
		enddo
		close(8, status='KEEP')
		deallocate(q3d, rboundary3d)
	
		if (twistFlag) then
			binfile='twist3d.bin'
			open(unit=8, file=binfile, access="direct",recl=qx*qy)
			do rec1=0, qz-1 
			  write(8,rec=1+rec1) twist3d(:, :,rec1)
			enddo
			close(8, status='KEEP')
			deallocate(twist3d)
		endif
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
	else
	!qcs
		tflag=0
		call qcs()
		binfile='qcs.bin'
		allocate(rboundary(0:q1-1, 0:q2-1))		
		open(unit=8, file=binfile, access="direct", recl=q1*q2)
		write(8, rec=1) q
		rboundary=rsboundary(0:q1-1, 0:q2-1)
		write(8, rec=2) rboundary
		rboundary=reboundary(0:q1-1, 0:q2-1)
		write(8, rec=3) rboundary
		write(8, rec=4) length
		deallocate(rboundary)
		if (twistFlag) write(8, rec=5) twistF
	 	close(8, status='KEEP')	 	
	endif
	deallocate(rsF, rsboundary, tflag)
endif
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
deallocate(Bfield, q, qtmp, reF, bnr, reboundary)
if (twistFlag) deallocate(twistF, curlB)
end

