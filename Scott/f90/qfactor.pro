PRO qfactor, bx, by, bz, xreg=xreg, yreg=yreg, zreg=zreg, factor=factor, maxsteps=maxsteps, $
                         twistFlag=twistFlag, csFlag=csFlag, $
                         nbridges=nbridges,  step=step, odir=odir, fstr=fstr, tol=tol, $
                         RK4Flag=RK4Flag, scottFlag=scottFlag, no_preview=no_preview
;+
; PURPOSE:
;   Calculate the squashing factor Q at the photosphere or a cross section
;   or a box volume, given a 3D magnetic field with uniform grids
; 
;   For details see:
;     Liu et al. (2016, ApJ)  
;
; INPUTS               
;   Bx, By, Bz:      3D magnetic field 
;   
; OPTIONAL INPUTS  
;   xreg, yreg, zreg:  pixel coordinates of the field region, in arrays of two elements
;                      default is to include all available pixels		               
;   	             --- If (xreg[0] ne xreg[1]) AND (yreg[0] ne yreg[1]) AND (zreg[0] ne zreg[1]) AND NOT csFlag
;   	                    caculate Q in a box volume
;   	                    invoke vFlag
;  		             --- If (xreg[0] eq xreg[1]) OR (yreg[0] eq yreg[1]) or (zreg[0] eq zreg[1])
;  		                    calculate Q in a cross section perpendicular to X or Y or Z axis
;  		                    invoke cFlag 
;  		             --- If zreg=[0, 0], calculate Q at the photosphere
;  		                    invoke z0Flag 
;  		             --- If csFlag is set, see below
;  		                    invoke cFlag 
;      
;   factor:     to bloat up the original resolution; i.e., grid spacing = 1/factor; default is 4 
;   
;   maxsteps:   maxium steps for stracing a field line at one direction; default is 40*(nx+ny+nz)/step
;
;   nbridges:   number of processors to engage; default is 8
;
;   twistFlag:  to calulate twist number Tw; see Liu et al. (2016, ApJ); default is 0
;   
;   csFlag:     to calulate Q in a cross section defined by three points; default is 0
;               point0=[xreg[0],yreg[0],zreg[0]] ; origin
;               point1=[xreg[1],yreg[1],zreg[0]] ; point0 -> point1, first axis
;               point2=[xreg[0],yreg[0],zreg[1]] ; point0 -> point2, second axis
;				   
;   step:       step size in tracing field lines (for RK4); default=1/4
;		
;   odir:       directory to save the results
;		
;   fstr:       filename to save the results 
;
;   tol:        tolerance of a step in RKF45; default is 10.^(-4)
;
;   RK4Flag:    to trace bline by RK4; default is 0B (RKF45)
;
;   scottFlag:  caculate Q and Q_perp by the method of Scott_2017_ApJ_848_117; default is 0B (method 3 of Pariat_2012_A&A_541_A78)
; 
;   no_preview: don't preduce PNG images for preview; default is 0B
;
; OUTPUTS:
;   q/qcs/q3d:   squashing factor   
;   slogq:       sign(Bz) x log_10(Q), calulated with all the field lines
;   slogq_orig:  only include field lines with both footpoints on the photoshere
;   q_perp:      q_perp in Titov 2007
;   twist:       twist number = \int \mu_0 J \cdot B /4\pi B^2 ds
;   rboundary:   nature of the ends of field lines, see 'subroutine vp_rboundary' in trace_bline.f90
;                0 - inside
;                1 - end at zmin
;                2 - end at zmax
;                3 - end at ymin
;                4 - end at ymax
;                5 - end at xmin
;                6 - end at xmax
;                7 - others 
;  rF, rsF, reF: coordinates of mapping points
; 
; PROCEDURES USED:
; ------IDL procedures used; those included in SolarSoft are not specified
;   check_slash()
;   doppler_color
;   doppler_color_mix            
;   sign2d()
;      
; ------FORTRAN procedures used
;   qfactor.f90
;   trace_bline.f90
;   trace_scott.f90
;
; ------COMPILATION: 
;      ifort -o qfactor.x qfactor.f90 -fopenmp -mcmodel=large -O3 -xHost -ipo
;   gfortran -o qfactor.x qfactor.f90 -fopenmp -mcmodel=large -O3
;
;  the efficiency of ifort (version 2021.2.0) is slightly (~1.3 times) faster than gfortran (version 9.3.0)
;      
;      
; MODIFICATION HISTORY:
;   Developed by R. Liu, J. Chen and Peijing, Zhang @ USTC
;   
;   Jun 30,2014 R. Liu @ USTC, IDL edition 
;   Jul  2,2014 R. Liu, introduce nchunks > nbridges to utilize idle child processes
;   Apr 21,2015 R. Liu and J. Chen, deal with field lines pass through the boundary other than bottom
;   Apr 29,2015 R. Liu and J. Chen, further optimization on shared memory
;   Apr 30,2015 R. Liu, correct the size of qmap
;   Apr 27,2015 R. Liu, qcs
;   Jun 15,2015 J. Chen, Fortran Edition, modify foot points with RK4_Boundary
;   Jul  8,2015 J. Chen, q3d
;   Oct 29,2015 J. Chen, deal with field lines touching the cut plane: use the plane quasi-perp to the field line;
;	  
;   Nov 1, 2015 J. Chen,
;		(1) fuse qcs and qfactor in qfactor.f90;  
;		(2) the cross section can be paralleled to the photosphere;
;		(3) set tmp_dir;
;		(4) further optimization on thread task allocation method;
;	  
;   Nov 4, 2015 J. Chen, introduce rboundary3d=byte(rsboundary3d+8*reboundary3d)
;		So if rboundary3d[i, j, k] eq 9B, the definition of q3d[i, j, k] is based on the bottom surface (Titov 2003).
;		
;   Jun 22,2016 J. Chen, add the map of field line length
;   Feb 15,2017 J. Chen, accelerated by 10%
;
;   Oct 30,2017 J. Chen,
;       add the map of Bnr=abs(Bn_local/Bn_target) at the photosphere
;       correct the output of Q: Modify -1,0 to NAN
;
;   Aug 28,2018 J. Chen, supplement Q at maginal points
;   May  1,2021 Peijing Zhang and J. Chen, trace field line with RKF45; RK4 is remainded, modify classic RK4 to the RK4 with 3/8-rule
;   May  1,2021 J. Chen, adapted to gfortran compiler
;   Jun  1,2021 J. Chen, supplement Q at the point where surronding 4 points have different mapping planes;
;                       forcibly convert the input Bx, By, Bz to float arrays  (Real(4) in Fortran)
;   Jun  5,2021 J. Chen, deal part of the points where surronding 4 points have different mapping planes with the method of Scott_2017_ApJ_848_117
;   Jun 10,2021 J. Chen, provide keyword maxsteps at input
;   Jun 11,2021 J. Chen, add the coordinates of mapping points to '*.sav' data
;   Jun 13,2021 J. Chen, mark 0 for inside and 7 for others in 'subroutine vp_rboundary'
;   Jun  5,2021 J. Chen, switch the order of indexes of Bfield in trace_bline.f90 for a higher efficiency (>1.6 times faster)
;   Jul  9,2021 J. Chen, provide the option with the method of Scott_2017_ApJ_848_117, and q_perp as an output
;   Jul 24,2021 J. Chen, improve I/O between IDL and Fortran
;
;   This software is provided without any warranty. Permission to use,
;   copy, modify. Distributing modified or unmodified copies is granted,
;   provided this disclaimer and information are included unchanged.
;-
sbx=size(Bx)
sby=size(By)
sbz=size(Bz)

if sbx[0] ne 3 or sby[0] ne 3 or sbz[0] ne 3 then message,'Bx, By and Bz must be 3D arrays!'
if sbx[1] ne sby[1] or sbx[1] ne sbz[1] or $
   sbx[2] ne sby[2] or sbx[2] ne sbz[2] or $
   sbx[3] ne sby[3] or sbx[3] ne sbz[3] then message,'Bx, By and Bz must have the same dimensions!' 
nx=sbz[1] & ny=sbz[2] & nz=sbz[3]
b2dz=Bz[*,*,0]

if ~keyword_set(xreg) then xreg=[0, nx-1]
if ~keyword_set(yreg) then yreg=[0, ny-1]
if ~keyword_set(zreg) then zreg=[0,    0]
if n_elements(xreg) ne 2 or n_elements(yreg) ne 2 or n_elements(zreg) ne 2 then $
   message,'xreg, yreg and zreg must be 2-element arrays!'

if keyword_set(factor) then factor=long(factor) else factor=4L
factor_str='f'+string(factor,format='(i02)')
if keyword_set(step) then step=step else step=1./4.

qx=abs(xreg[1]-xreg[0])*factor+1
qy=abs(yreg[1]-yreg[0])*factor+1
qz=abs(zreg[1]-zreg[0])*factor+1

; calculate Q in a cross section parallel to the photosphere
zFlag = xreg[0] ne xreg[1] AND yreg[0] ne yreg[1] AND zreg[0] eq zreg[1] 
if zFlag then begin 
	q1=qx & q2=qy
endif
if zFlag and zreg[0] eq 0 then z0Flag=1B else z0Flag=0B

; calculate Q in a cross section perpendicular to X- or Y-axis 
xFlag = xreg[0] eq xreg[1] AND yreg[0] ne yreg[1] AND zreg[0] ne zreg[1] 
if xFlag then begin 
	q1=qy & q2=qz
endif

yFlag = xreg[0] ne xreg[1] AND yreg[0] eq yreg[1] AND zreg[0] ne zreg[1]
if yFlag then begin 
	q1=qx & q2=qz
endif

; calculate Q in a cross section perpendicular to the photosphere
if keyword_set(csFlag) then csFlag=1B else csFlag=0B
cFlag = xFlag OR yFlag OR zFlag OR csFlag AND NOT z0Flag
 

if csFlag then begin
	point0=[xreg[0],yreg[0],zreg[0]]
	point1=[xreg[1],yreg[1],zreg[0]]
	point2=[xreg[0],yreg[0],zreg[1]]
	if 	((total(point0 eq point1) eq 3 ) or (total(point0 eq point2) eq 3 ) or (total(point1 eq point2) eq 3 )) then $
      message,'Something is wrong with the cross section .......'
	q1=sqrt(total(Double(point1-point0)*(point1-point0)))*factor+1
	q2=sqrt(total(Double(point2-point0)*(point2-point0)))*factor+1
endif else begin
;	point* not use here but need input to head.txt
	point0=fltarr(3) & point1=fltarr(3) & point2=fltarr(3)	
endelse


if cFlag then begin
	if (q1 eq 1 OR q2 eq 1) then message,'Something is wrong with the cross section .......'
end

; caculate Q in a box volume 
vFlag = xreg[0] ne xreg[1] AND yreg[0] ne yreg[1] AND zreg[0] ne zreg[1] AND NOT csFlag
if vFlag then begin
	q1=qx & q2=qy
	if qx eq 1 OR qy eq 1 OR qz eq 1 then message,'Something is wrong with the box volume .......'
endif

if vflag then head_str='q3d_' else head_str='qcs_'
cut_str='' 
if xFlag then cut_str='_x'+string(xreg[0],format='(I0)')
if yFlag then cut_str='_y'+string(yreg[0],format='(I0)')
if zFlag then cut_str='_z'+string(zreg[0],format='(I0)')
if ~keyword_set(fstr) then fstr = head_str + factor_str + cut_str
if  keyword_set(odir) then odir=check_slash(odir) else odir= check_slash(curdir())+'qfactor/'
if   ~dir_exist(odir) then file_mkdir,odir
print,'Results to be saved in ', odir+fstr+'.sav' 

; the directory for transmission between Fortran and IDL
;tmp_dir= check_slash(odir+'tmp')
tmp_dir='/dev/shm/tmp/'  ; use RAM
if not file_exist(tmp_dir)  then file_mkdir, tmp_dir

if ~keyword_set(nbridges)   then nbridges=8
max_threads=!CPU.HW_NCPU
nbridges=nbridges < max_threads
if  keyword_set(twistFlag)  then twistFlag =1B else twistFlag =0B
if  keyword_set(RK4Flag)    then RK4Flag   =1B else RK4Flag   =0B
if  keyword_set(scottFlag)  then scottFlag =1B else scottFlag =0B
if  keyword_set(no_preview) then no_preview=1B else no_preview=0B
if ~keyword_set(tol)        then tol=10.0^(-4.)
if ~keyword_set(maxsteps)   then maxsteps=long(4*(nx+ny+nz)/step)
;###############################################################

dummy=file_search(tmp_dir+'*.txt',count=nf)
if nf ne 0 then spawn,'rm '+tmp_dir+'*.txt'
dummy=file_search(tmp_dir+'*.bin',count=nf)
if nf ne 0 then spawn,'rm '+tmp_dir+'*.bin'
;----------------------------------------------------------------------------------------------
; mark the area for calculation on magnetogram
if (~no_preview) then begin
	cur_device=!D.name
	SET_PLOT, 'Z' 
	loadct,0
	DEVICE, SET_RESOLUTION=[nx,ny]

	if max(abs(b2dz)) ge 1000 then begin
		tv,bytscl(b2dz,min=-1000,max=1000,/nan)
	endif else begin
		tv,bytscl(b2dz,min=-max(abs(b2dz))/2.,max=max(abs(b2dz))/2.,/nan)
	endelse

	if csflag then begin
		plots,[xreg[0],xreg[1]],[yreg[0],yreg[1]],/dev
	endif else begin
		plots,[xreg[0],xreg[1],xreg[1],xreg[0],xreg[0]],[yreg[0],yreg[0],yreg[1],yreg[1],yreg[0]],/dev
	endelse
	im=tvread()
	write_png,odir+fstr+'_bz.png',im
	set_plot, cur_device
endif
;----------------------------------------------------------------------------------------------
;  transmit data to Fortran
get_lun,unit
openw,  unit, tmp_dir+'head.txt'
printf, unit, nx, ny, nz, nbridges, factor, long(maxsteps)
printf, unit, float(xreg), float(yreg), float(zreg),  step, tol
printf, unit, long(qx), long(qy), long(qz), long(q1), long(q2)
printf, unit, long(twistFlag), long(RK4flag), long(scottFlag), long(csflag)
printf, unit, float(point0), float(point1), float(point2)
close,  unit

openw,  unit, tmp_dir+'b3d.bin'
writeu, unit, float(temporary(Bx)), float(temporary(By)), float(temporary(Bz))
close,  unit
;----------------------------------------------------------------------------------------------
; calculate in Fortran
cdir=curdir()
CD, tmp_dir
tnow=systime(1)
spawn,'qfactor.x' ; if not known by the system, specify the path
tend=systime(1)
cd, cdir
tcalc=tend-tnow

if (tcalc ge 3600) then begin
	time_elapsed=string(tcalc/3600.0,format='(f0.2)')+' hours'
endif else begin 
	if (tcalc ge 60) then begin 
		time_elapsed=string(tcalc/60.0,format='(f0.2)')+' minutes'
	endif else begin 
		time_elapsed=string(tcalc,format='(f0.2)')+' seconds'
	endelse
endelse
print, time_elapsed+' elapsed for calculating qfactor' 
; ################################### retrieving results ###################################################### 
; load color table
loadct,0
tvlct,r_0,g_0,b_0,/get
doppler_color
tvlct,r_doppler,g_doppler,b_doppler,/get
doppler_color_mix
tvlct,r_doppler_mix,g_doppler_mix,b_doppler_mix,/get

; Q at the photosphere ----------------------------------------------------------------------------------------------
IF z0Flag THEN BEGIN
  ; get data From Fortran
	q=fltarr(q1,q2)
	rboundary=lonarr(q1,q2)
  	length=fltarr(q1,q2)
	Bnr=fltarr(q1,q2)
	openr, unit, tmp_dir+'qfactor0.bin'
	readu, unit, q, rboundary, length, Bnr
	close, unit
	
  ; save results and slogQ map
  ; only deal with field lines with both footpoints anchored at the bottom
 	sign_congrid_b2dz=sign2d(congrid(b2dz[xreg[0]:xreg[1],yreg[0]:yreg[1]], qx, qy, /interp, /minus_one))
  	ss=where(rboundary ne 1,complement=nss)
	qtmp=q
	if(ss[0] ne -1) then qtmp[ss]=!values.F_NAN  
  	slogq_orig=sign_congrid_b2dz*alog10(qtmp>1)
  	
  ; include all field lines  
  	slogq=sign_congrid_b2dz*alog10(q>1.)
  	
  	if (~no_preview) then begin
		im=bytscl(slogq_orig,min=-5,max=5,/nan)
		write_png,odir+fstr+'_slogq_orig.png',im,r_doppler,g_doppler,b_doppler
		
		im=bytscl(length,min=0,max=sqrt(nx^2.0+ny^2.0+nz^2.0),/nan)  
		write_png,odir+fstr+'_length.png',im,r_0,g_0,b_0
		
		im=bytscl(alog10(Bnr),min=-2,max=2,/nan)
		write_png,odir+fstr+'_lg(Bnr).png',im,r_doppler,g_doppler,b_doppler		
		
		im_slogq_mix=bytarr(q1,q2)
		if( ss[0] ne -1) then im_slogq_mix[ ss]=bytscl(slogq[ ss],min=-5,max=5,/nan,top=127)+128B  ; green-white-yellow
		if(nss[0] ne -1) then im_slogq_mix[nss]=bytscl(slogq[nss],min=-5,max=5,/nan,top=127)    ; blue-white-red
		write_png,odir+fstr+'_slogq_mix.png',im_slogq_mix, r_doppler_mix,g_doppler_mix,b_doppler_mix
	endif
		
    ; q_perp map
	if scottFlag then begin
		q_perp=fltarr(q1,q2)
		openr, unit, tmp_dir+'q_perp.bin'
		readu, unit, q_perp
		close, unit
		
		
		qtmp=q_perp
		if(ss[0] ne -1) then qtmp[ss]=!values.F_NAN
		slogq_perp_orig=sign_congrid_b2dz*alog10(qtmp>1)
		slogq_perp=sign_congrid_b2dz*alog10(q_perp>1.)
		
		if (~no_preview) then begin
			im=bytscl(slogq_perp_orig,min=-5,max=5,/nan)
			write_png,odir+fstr+'_slogq_perp_orig.png',im,r_doppler,g_doppler,b_doppler	
			if( ss[0] ne -1) then im_slogq_mix[ ss]=bytscl(slogq_perp[ ss],min=-5,max=5,/nan,top=127)+128B  ; green-white-yellow
			if(nss[0] ne -1) then im_slogq_mix[nss]=bytscl(slogq_perp[nss],min=-5,max=5,/nan,top=127)    ; blue-white-red
			write_png,odir+fstr+'_slogq_perp.png',im_slogq_mix, r_doppler_mix,g_doppler_mix,b_doppler_mix
		endif
	endif
  
  
  ; twist map
	if twistFlag then begin
		twist=fltarr(qx,qy)
		openr, unit, tmp_dir+'twist.bin'
		readu, unit, twist
		close, unit
		if (~no_preview) then begin
			im=bytscl(twist,min=-5,max=5,/nan)
			write_png,odir+fstr+'_twist.png',im,r_doppler,g_doppler,b_doppler
		endif
	endif

	rF=fltarr(3,q1,q2)
	openr, unit, tmp_dir+'reF.bin'
	readu, unit, rF
	close, unit
	
	if scottFlag then begin
	 	if twistFlag then save,filename=odir+fstr+'.sav', slogq, slogq_orig, q, length,  Bnr, rboundary, xreg, yreg, zreg, factor, rF, q_perp, slogq_perp, slogq_perp_orig, twist $ 
		             else save,filename=odir+fstr+'.sav', slogq, slogq_orig, q, length,  Bnr, rboundary, xreg, yreg, zreg, factor, rF, q_perp, slogq_perp, slogq_perp_orig
	endif else begin
		if twistFlag then save,filename=odir+fstr+'.sav', slogq, slogq_orig, q, length,  Bnr, rboundary, xreg, yreg, zreg, factor, rF, twist $ 
		             else save,filename=odir+fstr+'.sav', slogq, slogq_orig, q, length,  Bnr, rboundary, xreg, yreg, zreg, factor, rF
	endelse
ENDIF 

; Q in the 3D box volume ----------------------------------------------------------------------------------------------
IF vflag THEN BEGIN
  ;get data From Fortran
	q3d=fltarr(qx,qy,qz)
	rboundary3d=bytarr(qx,qy,qz)	
	openr,unit, tmp_dir+'q3d.bin'
	readu, unit, q3d, rboundary3d
	close, unit
  ;	rboundary3d=rsboundary3d+8*reboundary3d  (It has been calculated in Fortran)
  ;	So if rboundary3d[i, j, k] eq 9B, the definition of q3d[i, j, k]  is based on the bottom surface

	if scottFlag then begin
		q_perp3d=fltarr(qx,qx,qz)
		openr, unit, tmp_dir+'q_perp3d.bin'
		readu, unit, q_perp3d
		close, unit
	endif
  

	if twistFlag then begin
		twist3d=fltarr(qx,qy,qz)
		openr,unit, tmp_dir+'twist3d.bin'
		readu, unit, twist3d
		close, unit
	endif
  ; save resutls 
	if scottFlag then begin
		if twistFlag then save,filename=odir+fstr+'.sav', q3d, rboundary3d, xreg, yreg, zreg, factor, q_perp3d, twist3d $
		             else save,filename=odir+fstr+'.sav', q3d, rboundary3d, xreg, yreg, zreg, factor, q_perp3d
	endif else begin
		if twistFlag then save,filename=odir+fstr+'.sav', q3d, rboundary3d, xreg, yreg, zreg, factor, twist3d $
		             else save,filename=odir+fstr+'.sav', q3d, rboundary3d, xreg, yreg, zreg, factor	
	endelse
			         
	sign_congrid_b2dz=sign2d(congrid(b2dz[xreg[0]:xreg[1],yreg[0]:yreg[1]], qx, qy, /interp, /minus_one))
  
	if (zreg[0] eq 0.0) then begin
		qtmp=q3d[*,*,0]
		rb=rboundary3d[*,*,0]
		ss=where(rb ne 9,complement=nss)
		if(ss[0] ne -1) then qtmp[ss]=!values.F_NAN
		slogq_orig=sign_congrid_b2dz*alog10(qtmp>1)
		if (~no_preview) then begin
			im=bytscl(slogq_orig,min=-5,max=5,/nan)
			write_png,odir+fstr+'_z0_slogq_orig.png',im,r_doppler,g_doppler,b_doppler
		endif
	endif  
ENDIF 

; Q at the cross section ----------------------------------------------------------------------------------------------
IF cFlag THEN BEGIN

	qcs=fltarr(q1,q2)
	rsboundary=lonarr(q1,q2)		
	reboundary=lonarr(q1,q2)
	length=fltarr(q1,q2)
	openr, unit, tmp_dir+'qcs.bin'
	readu, unit, qcs, rsboundary, reboundary, length
	close, unit
	
	logq=alog10(qcs>1.)
	qcs_orig=qcs
	ss1=where(rsboundary ne 1,complement=nss1)
	ss2=where(reboundary ne 1,complement=nss1)
	if(ss1[0] ne -1) then qcs_orig[ss1]=!values.F_NAN
	if(ss2[0] ne -1) then qcs_orig[ss2]=!values.F_NAN

	
	if (~no_preview) then begin
		im=bytscl(length,min=0,max=sqrt(nx^2.0+ny^2.0+nz^2.0),/nan)
		write_png,odir+fstr+'_length.png',im,r_0,g_0,b_0
		im=bytscl(logq,min=0,max=5,/nan)
		write_png,odir+fstr+'_logq.png',im,r_0,g_0,b_0
		im=bytscl(alog10(qcs_orig>1.),min=0,max=5,/nan)
		write_png,odir+fstr+'_logq_orig.png',im,r_0,g_0,b_0
	endif

	if scottFlag then begin
		q_perp=fltarr(q1,q2)
		openr, unit, tmp_dir+'q_perp.bin'
		readu, unit, q_perp
		close, unit
		
		q_perp_orig=q_perp
		if(ss1[0] ne -1) then q_perp_orig[ss1]=!values.F_NAN
		if(ss2[0] ne -1) then q_perp_orig[ss2]=!values.F_NAN
		slogq_perp_orig=alog10(q_perp_orig>1)
		slogq_perp=alog10(q_perp>1.)
		
		
		if (~no_preview) then begin
			im=bytscl(slogq_perp_orig,min=0,max=5,/nan)
			write_png,odir+fstr+'_logq_perp_orig.png',im, r_0, g_0, b_0
			im_slogq_perp=bytscl(slogq_perp,min=0,max=5,/nan)
			write_png,odir+fstr+'_logq_perp.png',im_slogq_perp, r_0, g_0, b_0
		endif
	endif

	if twistFlag then begin	
		twist=fltarr(q1,q2)
		openr, unit, tmp_dir+'twist.bin'
		readu, unit, twist
		close, unit
		if (~no_preview) then begin
			im=bytscl(twist,min=-5,max=5,/nan)
			write_png,odir+fstr+'_twist.png',im,r_doppler,g_doppler,b_doppler
		endif
	endif	   

	reF=fltarr(3,q1,q2)
	openr, unit, tmp_dir+'reF.bin'
	readu, unit, reF
	close, unit
	rsF=fltarr(3,q1,q2)
	openr, unit, tmp_dir+'reF.bin'
	readu, unit, rsF
	close, unit
	
	
	if scottFlag then begin
		if twistFlag then save,filename=odir+fstr+'.sav', qcs, qcs_orig, length, rsboundary, reboundary, xreg, yreg, zreg, factor, csFlag, rsF, reF, q_perp, q_perp_orig, twist $
		             else save,filename=odir+fstr+'.sav', qcs, qcs_orig, length, rsboundary, reboundary, xreg, yreg, zreg, factor, csFlag, rsF, reF, q_perp, q_perp_orig	
	endif else begin
		if twistFlag then save,filename=odir+fstr+'.sav', qcs, qcs_orig, length, rsboundary, reboundary, xreg, yreg, zreg, factor, csFlag, rsF, reF, twist $
		             else save,filename=odir+fstr+'.sav', qcs, qcs_orig, length, rsboundary, reboundary, xreg, yreg, zreg, factor, csFlag, rsF, reF
	endelse
	
ENDIF
;----------------------------------------------------------------------------------------------	
; hourse keeping
spawn, 'rm -r '+tmp_dir
free_lun, unit
print,'----------------------Done----------------------'
END

