subroutine  interpolate_grad_unit_vec_B2(vp, unit_vec_bp, grad_unit_vec_B)
!!the way of https://github.com/Kai-E-Yang/QSL
implicit none
real:: vp(0:2),vp1(0:2),vp2(0:2), unit_vec_bp(0:2), unit_vec_bp1(0:2), unit_vec_bp2(0:2), grad_unit_vec_B(0:2,0:2)
integer:: i
!----------------------------------------------------------------------------
do i=0,2
	vp1=vp
	vp2=vp
	vp1(i)=vp(i)-0.001
	vp2(i)=vp(i)+0.001
	call interpolate_unit_vec_B(vp1, unit_vec_bp1)
	call interpolate_unit_vec_B(vp2, unit_vec_bp2)
	grad_unit_vec_B(i,0:2)=(unit_vec_bp2-unit_vec_bp1)/0.002
enddo
call interpolate_unit_vec_B(vp, unit_vec_bp)
end subroutine interpolate_grad_unit_vec_B2


subroutine interpolate_grad_unit_vec_B3(vp, unit_vec_bp, grad_unit_vec_B)
!the way of https://bitbucket.org/tassev/qsl_squasher/src/hg/
use trace_common
use field_common
implicit none
real:: w(0:1,0:2), dw(0:1), weight(0:1,0:1,0:1), vp(0:2), bp(0:2), unit_vec_bp(0:2), &
grad_unit_vec_B(0:2,0:2), unit_vec_B_cell(0:2,0:1,0:1,0:1)
integer:: round(0:1,0:2), i, j, k
!----------------------------------------------------------------------------
round(0,:)=floor(vp)
w(1,:)=vp-round(0,:)

do i=0,2
	if ( .not. (vp(i) .ge. 0.0)) then
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
forall(i=0:2) Bp(i)=sum(weight*Bfield(i, round(:,0), round(:,1), round(:,2)))
unit_vec_bp=bp/norm2(bp)

dw(0)=-1.0
dw(1)= 1.0

forall(i=0:1,j=0:1,k=0:1) unit_vec_B_cell(:,i,j,k)= &
Bfield(:, round(i,0), round(j,1), round(k,2))/norm2(Bfield(:, round(i,0), round(j,1), round(k,2)))

forall(i=0:1,j=0:1,k=0:1) weight(i,j,k)=dw(i)*w(j,1)*w(k,2)
forall(i=0:2)  grad_unit_vec_B(0,i)=sum(weight*unit_vec_B_cell(i,:,:,:))

forall(i=0:1,j=0:1,k=0:1) weight(i,j,k)=w(i,0)*dw(j)*w(k,2)
forall(i=0:2)  grad_unit_vec_B(1,i)=sum(weight*unit_vec_B_cell(i,:,:,:))

forall(i=0:1,j=0:1,k=0:1) weight(i,j,k)=w(i,0)*w(j,1)*dw(k)
forall(i=0:2)  grad_unit_vec_B(2,i)=sum(weight*unit_vec_B_cell(i,:,:,:))

end subroutine interpolate_grad_unit_vec_B3



subroutine interpolate_grad_unit_vec_B3_stretch(vp, unit_vec_bp, grad_unit_vec_B)
!the way of https://bitbucket.org/tassev/qsl_squasher/src/hg/
use trace_common
use field_common
implicit none
real:: w(0:1,0:2), dw(0:1), weight(0:1,0:1,0:1), vp(0:2), bp(0:2), unit_vec_bp(0:2), &
grad_unit_vec_B(0:2,0:2), unit_vec_B_cell(0:2,0:1,0:1,0:1), vpBound
integer:: round(0:1,0:2), i, j, k, index_i, index_j, index_k
!----------------------------------------------------------------------------
call vp_index(vp, vpBound, index_i, index_j, index_k)

round(0,:)=[index_i, index_j, index_k]
round(1,:)=round(0,:)+1

w(0,0)=(xa(index_i+1)-vpBound(0))/dxa(index_i)
w(0,1)=(ya(index_j+1)-vpBound(1))/dya(index_j)
w(0,2)=(za(index_k+1)-vpBound(2))/dza(index_k)

w(1,:)=1.0-w(0,:)

forall(i=0:1,j=0:1,k=0:1) weight(i,j,k)=w(i,0)*w(j,1)*w(k,2)
forall(i=0:2) Bp(i)=sum(weight*Bfield(i, round(:,0), round(:,1), round(:,2)))
unit_vec_bp=bp/norm2(bp)

forall(i=0:1,j=0:1,k=0:1) unit_vec_B_cell(:,i,j,k)= &
Bfield(:, round(i,0), round(j,1), round(k,2))/norm2(Bfield(:, round(i,0), round(j,1), round(k,2)))

dw(0)=-1.0
dw(1)= 1.0

forall(i=0:1,j=0:1,k=0:1) weight(i,j,k)=dw(i)/dxa(index_i)*w(j,1)*w(k,2)
forall(i=0:2)  grad_unit_vec_B(0,i)=sum(weight*unit_vec_B_cell(i,:,:,:))

forall(i=0:1,j=0:1,k=0:1) weight(i,j,k)=w(i,0)*dw(j)/dya(index_j)*w(k,2)
forall(i=0:2)  grad_unit_vec_B(1,i)=sum(weight*unit_vec_B_cell(i,:,:,:))

forall(i=0:1,j=0:1,k=0:1) weight(i,j,k)=w(i,0)*w(j,1)*dw(k)/dza(index_k)
forall(i=0:2)  grad_unit_vec_B(2,i)=sum(weight*unit_vec_B_cell(i,:,:,:))

end subroutine interpolate_grad_unit_vec_B3_stretch


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
