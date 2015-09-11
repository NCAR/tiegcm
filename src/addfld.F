!
      module addfld_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Save fields to secondary history (B. Foster: November, 2005)
! This module replaces the old addfsech routines.
!
! 7/13/06 btf: Using fsechist%data associated status to determine 
!              firstcall, rather than unames.
!
      use params_module
      implicit none
      contains
!-----------------------------------------------------------------------
      subroutine addfld(name,long_name,units,
     |  f,dname1,lb1,ub1,dname2,lb2,ub2,idx)
!
! Purpose:
!   Save user-provided data f to secondary history field structure.
!
!   Basically, the user provides a 2d array in any dimension combination
!   (geo or mag) containing all or part of the global domain, along with 
!   the field name, long_name, and units. This data is saved in either
!   2d (lon,lat), or 3d (lon,lat,lev) data arrays (fsechist(i)%data) for 
!   the final history variable at the current model time.
!
!   If all 3 spatial dimensions are desired on the history, idx is the 
!   index of the 3rd dimension (i.e., the dimension not specified by 
!   dname1 and dname2). In this case, addfld is called from inside a
!   loop over the 3rd dimension. If only 2 spatial dimensions are desired 
!   (lon,lat) then idx can be zero. However, only the levs dimension can 
!   be excluded (i.e., idx can be zero only if saving a horizontal slice, 
!   see table below)
!
!   This routine saves the data and field information in the derived type 
!   fsechist(n) (see fields.F). The final variable on the history will be 
!   either (lon,lat,lev,time), or (lon,lat,time), geo or mag (see subs 
!   def_fsech and wrfsech in nchist.F).
!
! On Input:
!   dname1     Name (char string) of 1st dimension of f (see table below)
!   lb1,ub1    Lower and upper integer bounds of 1st dimension of f
!              (must be within full domain dimension sizes, see table below)
!   dname2     Name (char string) of 2nd dimension of f (see table below)
!   lb2,ub2    Lower and upper integer bounds of 2nd dimension of f
!              (must be within full domain dimension sizes, see table below)
!   f(lb1:ub1,lb2:ub2) 
!              Full or sub domain data to be saved for secondary history
!   idx        Index of 3rd dimension (or 0 if writing horizontal slice only)
!              (idx can be 0 only if providing f(lon,lat) or f(lat,lon))
!   name       Name of field (should be in namelist SECFLDS) (<= 16 chars)
!   long_name  Optional long name of the field (<= 80 chars)
!   units      Optional units of the field (<= 16 chars)
!
! On Output:
!   fsechist(n) is defined (derived type fields_sech, see fields.F). Data 
!   component has been allocated, and data from f (full or sub-domain) 
!   has been saved to fsechist(n)%data(lon,lat,lev)
!
! Caveat: If k > 0, midpoints is assumed.
!
! Note 1: nlevp1==nilevp1==nmlevp1==nimlevp1
! Note 2: For geo fields, only i=3:nlonp2 are written to histories (see sub wrfsech)
! Note 3: If the input lon dimension bounds are the full domain dimension
!         (e.g., lb1,ub1==1,nlonp4), then it is assumed that only the master
!         task is making the call for the current field. In this case (task0_only), 
!         the data are not gathered to the master task when writing the history
!         (see sub mp_gather2root_fsech in mpi.F)
!
! Following is a table of valid combinations for dname1 and dname2, with the
!   corresponding full-domain dimension sizes. Final fields on the secondary 
!   histories will be one of two types: f(lon,lat,lev,time) (when idx > 0), 
!   or f(lon,lat,time) (when idx==0).
!
! dname1 dname2 idx  dim sizes         Description
!                    (full domain)     (of input f)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! lev    lon    j    nlevp1,nlonp4     geo lat slice at j (midpoint levs)
! ilev   lon    j    nilevp1,nlonp4    geo lat slice at j (interface levs)
! lon    lev    j    nlonp4,nlevp1     geo lat slice at j (midpoint levs)
! lon    ilev   j    nlonp4,nilevp1    geo lat slice at j (interface levs)
!
! lev    lat    i    nlevp1,nlat       geo lon slice at i (midpoint levs)
! ilev   lat    i    nilevp1,nlat      geo lon slice at i (interface levs)
! lat    lev    i    nlat,nlevp1       geo lon slice at i (midpoint levs)
! lat    ilev   i    nlat,nilevp1      geo lon slice at i (interface levs)
!
! lon    lat    k    nlonp4,nlat       geo horizontal slice at k
! lat    lon    k    nlat,nlonp4       geo horizontal slice at k
! lon    lat    0    nlonp4,nlat       geo horizontal slice (2d)
! lat    lon    0    nlat,nlonp4       geo horizontal slice (2d)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! mlev   mlon   j    nmlevp1,nmlonp1   mag lat slice at j (midpoint levs)
! imlev  mlon   j    nimlevp1,nmlonp1  mag lat slice at j (interface levs)
! mlon   mlev   j    nmlonp1,nmlevp1   mag lat slice at j (midpoint levs)
! mlon   imlev  j    nmlonp1,nimlevp1  mag lat slice at j (interface levs)
!
! mlev   mlat   i    nmlevp1,nmlat     mag lon slice at i (midpoint levs)
! imlev  mlat   i    nimlevp1,nmlat    mag lon slice at i (interface levs)
! mlat   mlev   i    nmlat,nmlevp1     mag lon slice at i (midpoint levs)
! mlat   imlev  i    nmlat,nimlevp1    mag lon slice at i (interface levs)
!
! mlon   mlat   k    nmlonp1,nmlat     mag horizontal slice at k
! mlat   mlon   k    nmlat,nmlonp1     mag horizontal slice at k
! mlon   mlat   0    nmlonp1,nmlat     mag horizontal slice (2d)
! mlat   mlon   0    nmlat,nmlonp1     mag horizontal slice (2d)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! Example 1: save (lev,lon) array at j for 3d+time field on the history:
!   real :: f(lev0:lev1,lon0:lon1) ! field on mpi task subdomains
!   do j=lat0,lat1
!     do i=lon0,lon1
!       do k=lev0,lev1
!         f(k,i) = [expression, probably involving k,i,j indices in other arrays]
!       enddo
!     enddo
!     call addfld('MY_3D_DIAG','My long name','my units',
!       f,'lev',lev0,lev1,'lon',lon0,lon1,j)
!   enddo
!
! Example 2: save (lon,lat) array for 2d+time field on the history:
!   real :: f(lon0:lon1,lat0:lat1) ! field on mpi task subdomains
!   do j=lat0,lat1
!     do i=lon0,lon1
!       f(i,j) = [expression, probably involving i,j indices in other arrays]
!     enddo
!   enddo
!   call addfld('MY_2D_DIAG','My long name','my units',
!     f,'lon',lon0,lon1,'lat',lat0,lat1,0)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! Use modules:
      use fields_module,only: fsechist,f4d,nf4d,shortname_len
      use hist_module,only: isechist,modeltime
      use init_module,only: istep
      use mpi_module,only: lat0,lat1,lon0,lon1
      implicit none
!
! Args:
      character(len=*),intent(in) :: dname1,dname2,name,long_name,units
      integer,intent(in) :: lb1,ub1,lb2,ub2,idx
      real,intent(in) :: f(lb1:ub1,lb2:ub2)
!
! Local:
      integer :: ier,i,j,k,ii,ix,ixf4d,iseries_p,iseries_s
      integer :: idim1,idim2,idim3,ihdim1,ihdim2,ihdim3
      logical :: time2write,wrprim,wrsech,newseries_prim,
     |  newseries_sech,found
      logical :: mag,geo,firstcall
      real :: fmin,fmax
      character(len=8) :: dname3
      integer :: istep1=-1
!
! Valid dimension names on history (dname1 and dname2 must be one of these)
      integer,parameter :: nvalid=8
      character(len=8) :: valid_dnames(nvalid) =
     |  (/'lev     ','ilev    ','lon     ','lat     ',  ! geo dim names
     |    'mlev    ','imlev   ','mlon    ','mlat    '/) ! mag dim names
!
! Dimensions of global fsechist(ix)%data corresponding to dim names:
      integer :: valid_dims(nvalid) =
     |  (/nlevp1    ,nilevp1   ,nlonp4    ,nlat      ,  ! geo dim sizes
     |    nmlevp1   ,nimlevp1  ,nmlonp1   ,nmlat     /) ! mag dim sizes
!
! External:
      integer,external :: strloc
      logical,external :: wrhist
!
! Return silently if not writing secondary histories:
      if (isechist==0) return
!
! Save only if we are writing a secondary history this step:
      time2write = wrhist(istep,modeltime,
     |  wrprim, newseries_prim, iseries_p,
     |  wrsech, newseries_sech, iseries_s) 
      if (.not.wrsech) return
      if (istep1 < 0) istep1 = istep
!
! Check for name in fsechist(mxfsech) (namelist SECFLDS)
! If field name was not included in namelist SECFLDS, then the data
! pointer is not allocated, and the field is not defined on secondary
! history. See checks for data allocation (.not.associated(fsechist(ix)%data))
! in nchist.F and mpi.F.
!
      ix = strloc(fsechist%short_name,mxfsech,name) ! mxfsech in params.F
      if (ix <= 0) then ! field was not requested by SECFLDS in namelist read file
!
! Return silently to avoid cluttering stdout with warnings.
! (often addfld calls are left open in the code, but the field was not
!  requested SECFLDS in the namelist read file)
!
        return
      endif
!
      firstcall = .false.
      if (.not.associated(fsechist(ix)%data)) firstcall = .true.
!
! Return with warning if field requested is already saved as 4d field 
! (the 4d field will be written to the secondary history by nchist module). 
! That is, the user added a call to addfld in which the field name (1st arg) 
! is already in use as a 4d field (e.g., call addfld('TN'...))
!
      ixf4d = strloc(f4d%short_name,nf4d,name) 
      if (ixf4d > 0) return
!
! Call is for a diagnostic secondary history field (2d or 3d, mag or geo).
!     write(6,"('addfld: ix=',i3,' name=',a)") ix,name
!
! Dim types 1 and 2 must be different:
      if ((index(dname1,'lev') > 0 .and. index(dname2,'lev') > 0).or.
     |    (index(dname1,'lon') > 0 .and. index(dname2,'lon') > 0).or.
     |    (index(dname1,'lat') > 0 .and. index(dname2,'lat') > 0)) then
        write(6,"('>>> addfld: redundant dimension types: dname1,2=',
     |    2a)") dname1,dname2
        call shutdown('addfld')
      endif
!
! Set idim1,idim2 (also validate dname1,2):
      idim1 = 0 ; idim2 = 0
      i = strloc(valid_dnames,nvalid,dname1)
      if (i > 0) idim1 = valid_dims(i)
      i = strloc(valid_dnames,nvalid,dname2)
      if (i > 0) idim2 = valid_dims(i)
      if (idim1==0) then
        write(6,"(/,'>>> addfld: invalid dname1 = ',a,' for field ',a)") 
     |    dname1,name
        write(6,"('Valid dim names: ',8a8)") valid_dnames
        call shutdown('addfld')
      endif
      if (idim2==0) then
        write(6,"(/,'>>> addfld: invalid dname2 = ',a,' for field ',a)") 
     |    dname2,name
        write(6,"('Valid dim names: ',8a8)") valid_dnames
        call shutdown('addfld')
      endif
!
! Validate dimension bounds:
      if (lb1 < 1 .or. ub1 > idim1) then
        write(6,"(/,'>>> addfld bad bounds: field=',a,
     |    ' dname1=',a,' idim1=',i4,' lb1=',i4,' ub1=',i4)") 
     |    trim(name),dname1,idim1,lb1,ub1
        call shutdown('addfld')
      endif
      if (lb2 < 1 .or. ub2 > idim2) then
        write(6,"(/,'>>> addfld bad bounds: field=',a,
     |    ' dname2=',a,' idim2=',i4,' lb2=',i4,' ub2=',i4)") 
     |    trim(name),dname2,idim2,lb2,ub2
        call shutdown('addfld')
      endif
!
! If idx==0, it must refer to lev dimension:
      if (idx==0.and.
     |  ((index(dname1,'lev') > 0.and.index(dname2,'lon') > 0).or.
     |   (index(dname1,'lev') > 0.and.index(dname2,'lat') > 0).or.
     |   (index(dname1,'lon') > 0.and.index(dname2,'lev') > 0).or.
     |   (index(dname1,'lat') > 0.and.index(dname2,'lev') > 0))) then
        write(6,"(/,'>>> addfld: idx can be zero only when saving',
     |    ' f(lon,lat) or f(lat,lon): dname1,2=',2a)") dname1,dname2
        call shutdown('addfld')
      endif
!
! If longitude dimension is global, assume task0_only, i.e., assume
! field is presented by master task only -- in this case gather2root
! must not be called on this field (sub mp_gather2root_sech in mpi.F)
!
      fsechist(ix)%task0_only = .false.
      if ((index(dname1,'lon') > 0 .and. lb1==1.and.ub1==idim1) .or.
     |    (index(dname2,'lon') > 0 .and. lb2==1.and.ub2==idim2))
     |  fsechist(ix)%task0_only = .true.
!
! Mag or geo:
      geo = .true. ; mag = .false.
      i = strloc(valid_dnames(5:8),4,dname1)
      if (i > 0) then
        geo = .false.
        mag = .true.
      endif
!
! Dims must be all mag or all geo:
! (valid_dnames(1:4) are geo, valid_dnames(5:8) are mag)
!
      if (geo) then
        i = strloc(valid_dnames(5:8),4,dname1)
        ii = strloc(valid_dnames(5:8),4,dname2)
        if (i > 0.or.ii > 0) then
          write(6,"(/,'>>> addfld: cannot mix geo and mag dimensions:',
     |      ' dname1,2=',2a)") dname1,dname2
        endif
      else
        i  = strloc(valid_dnames(1:4),4,dname1)
        ii = strloc(valid_dnames(1:4),4,dname2)
        if (i > 0.or.ii > 0) then
          write(6,"(/,'>>> addfld: cannot mix geo and mag dimensions:',
     |      ' dname1,2=',2a)") dname1,dname2
        endif
      endif
!
! Set full-domain secondary output array dimensions from dname1,2:
! (geographic or magnetic, lev,lat,or lon)
      fsechist(ix)%dimnames(:) = '        '
      fsechist(ix)%dimsizes(:) = 0

      if (index(dname1,'lon') > 0) fsechist(ix)%dimsizes(1) = idim1
      if (index(dname1,'lat') > 0) fsechist(ix)%dimsizes(2) = idim1
      if (index(dname1,'lev') > 0) fsechist(ix)%dimsizes(3) = idim1
      if (index(dname2,'lon') > 0) fsechist(ix)%dimsizes(1) = idim2
      if (index(dname2,'lat') > 0) fsechist(ix)%dimsizes(2) = idim2
      if (index(dname2,'lev') > 0) fsechist(ix)%dimsizes(3) = idim2
!
      if (index(dname1,'lon') > 0) fsechist(ix)%dimnames(1) = dname1
      if (index(dname1,'lat') > 0) fsechist(ix)%dimnames(2) = dname1
      if (index(dname1,'lev') > 0) fsechist(ix)%dimnames(3) = dname1
      if (index(dname2,'lon') > 0) fsechist(ix)%dimnames(1) = dname2
      if (index(dname2,'lat') > 0) fsechist(ix)%dimnames(2) = dname2
      if (index(dname2,'lev') > 0) fsechist(ix)%dimnames(3) = dname2
!
! Set third dimension of full-domain output array according to idx:
!
! 11/05 btf: If the 3rd dimension is lev, we don't know
!   if idx > 0 refers to midpoints or interfaces. This will 
!   happen only when addfld is called from a k-loop and idx > 0.
!
      idim3 = 0                                ! 2d array
      fsechist(ix)%ndims = 2
! 3d field:
      if (idx > 0) then                        ! 3d array
        fsechist(ix)%ndims = 3
        if (index(dname1,'lev')==0.and.index(dname2,'lev')==0) then
          idim3 = nlevp1 ; if (mag) idim3 = nmlevp1 ! but could be nilevp1 or nimlevp1
        elseif (index(dname1,'lat')==0.and.index(dname2,'lat')==0) then
          idim3 = nlat ; if (mag) idim3 = nmlat
        elseif (index(dname1,'lon')==0.and.index(dname2,'lon')==0) then
          idim3 = nlonp4 ; if (mag) idim3 = nmlonp1
        endif
        if (idim3==0) then
          write(6,"('>>> addfld: cannot determine idim3: dname1,2=',a,',
     |      ',a)") trim(dname1),trim(dname2)
          call shutdown('addfld')
        endif
        dname3 = '        '
        do i=1,nvalid
          if (idim3==valid_dims(i)) dname3 = valid_dnames(i)
        enddo
        do i=1,3
          if (len_trim(fsechist(ix)%dimnames(i))==0) then
            fsechist(ix)%dimsizes(i) = idim3
            fsechist(ix)%dimnames(i) = dname3
          endif
        enddo
        if (idx > idim3) then
          write(6,"(/,'>>> addfld: field ',a,': bad idx=',i4,
     |      ' dname3=',a,' idim3=',i4,' (idx should not be > idim3)')")
     |      trim(name),idx,dname3,idim3
          call shutdown('addfld')
        endif
      endif
!
      do i=1,fsechist(ix)%ndims
        if (len_trim(fsechist(ix)%dimnames(i))==0.or.
     |    fsechist(ix)%dimsizes(i)==0) then
          write(6,"(/,'>>> addfld: could not determine dimension ',i2,
     |      ': field ',a,': dimname=',a,' dimsize=',i3)") i,
     |      trim(name),fsechist(ix)%dimnames(i),fsechist(ix)%dimsizes(i)
          call shutdown('addfld')
        endif
      enddo
!
!     write(6,"('addfld: field ',a,' fsechist(ix)%dimnames=',a,',',a,
!    |  ',',a,' dimsizes=',3i4)") name,trim(fsechist(ix)%dimnames(1)),
!    |  trim(fsechist(ix)%dimnames(2)),trim(fsechist(ix)%dimnames(3)),
!    |  fsechist(ix)%dimsizes
!
! Allocate fsechist(ix)%data, either 2d (lon,lat) or 3d (lon,lat,lev):
! (field is mag or geo)
!
! 3d (i,j,k):
      if (.not.associated(fsechist(ix)%data)) then
        ihdim1=fsechist(ix)%dimsizes(1)
        ihdim2=fsechist(ix)%dimsizes(2)
        ihdim3=fsechist(ix)%dimsizes(3)
        if (idx > 0) then
          allocate(fsechist(ix)%data(ihdim1,ihdim2,ihdim3),stat=ier)
          if (ier /= 0) then
            write(6,"(/,'>>> Error allocating 3d (ihdim1=',i3,
     |        ',ihdim2=',i3,',ihdim3=',i3,') for field ',a)") 
     |        ihdim1,ihdim2,ihdim3,trim(fsechist(ix)%short_name)
            call shutdown('addfld')
          else
            write(6,"(/,'Allocated 3d sech field ',a,'(',a,'=',i3,
     |        ',',a,'=',i3,',',a,'=',i3,')')") 
     |        trim(fsechist(ix)%short_name),
     |        trim(fsechist(ix)%dimnames(1)),ihdim1,
     |        trim(fsechist(ix)%dimnames(2)),ihdim2,
     |        trim(fsechist(ix)%dimnames(3)),ihdim3
          endif
!
! 2d (i,j):
        else ! idx==0 (no lev dimension)
          allocate(fsechist(ix)%data(ihdim1,ihdim2,1),stat=ier)
          if (ier /= 0) then
            write(6,"(/,'>>> Error allocating 2d (ihdim1=',i3,
     |        ' ihdim2=',i3,' for field ',a)") ihdim1,ihdim2,
     |        trim(fsechist(ix)%short_name)
            call shutdown('addfld')
          else
            write(6,"(/,'Allocated 2d sech field ',a,'(',a,'=',i3,
     |        ',',a,'=',i3,')')") trim(fsechist(ix)%short_name),
     |        trim(fsechist(ix)%dimnames(1)),ihdim1,
     |        trim(fsechist(ix)%dimnames(2)),ihdim2
          endif
        endif
        fsechist(ix)%data = spval ! init whole array to missing value
      endif ! this field data not associated
!
! Assign data to fsech(ix)%data(i,j,k) from f, according to dim names:
!
! f(lev,lon), idx=j
      if (index(dname1,'lev') > 0.and.index(dname2,'lon') > 0) then
        do k=lb1,ub1
          do i=lb2,ub2
            fsechist(ix)%data(i,idx,k) = f(k,i)
          enddo

!         if (trim(fsechist(ix)%short_name)=='RP')
!    |    write(6,"('addfld: idx=',i3,' k=',i3,' name=',a,
!    |      ' data(:,idx,k)=',/,(6e12.4))") idx,k,
!    |      fsechist(ix)%short_name,fsechist(ix)%data(:,idx,k)

        enddo
!
! f(lev,lat), idx=i
      elseif (index(dname1,'lev') > 0.and.index(dname2,'lat') > 0) then 
        do k=lb1,ub1
          do j=lb2,ub2
            fsechist(ix)%data(idx,j,k) = f(k,j)
          enddo
        enddo
!
! f(lon,lev), idx=j
      elseif (index(dname1,'lon') > 0.and.index(dname2,'lev') > 0) then
        do k=lb2,ub2
          do i=lb1,ub1
            fsechist(ix)%data(i,idx,k) = f(i,k)
          enddo
        enddo
!
! f(lon,lat), idx = k or 0
      elseif (index(dname1,'lon') > 0.and.index(dname2,'lat') > 0) then 
        if (idx==0) then
          do j=lb2,ub2
            do i=lb1,ub1
              fsechist(ix)%data(i,j,1) = f(i,j)
            enddo
          enddo
        else  ! 3d
          do j=lb2,ub2
            do i=lb1,ub1
              fsechist(ix)%data(i,j,idx) = f(i,j)
            enddo
          enddo
        endif
!
! f(lat,lev), idx = i
      elseif (index(dname1,'lat') > 0.and.index(dname2,'lev') > 0) then
        do i=lb2,ub2
          do j=lb1,ub1  
            fsechist(ix)%data(idx,j,k) = f(j,k)
          enddo
        enddo
!
! f(lat,lon), idx = k or 0
      elseif (index(dname1,'lat') > 0.and.index(dname2,'lon') > 0) then
        if (idx==0) then
          do j=lb1,ub1
            do i=lb2,ub2
              fsechist(ix)%data(i,j,1) = f(j,i)
            enddo
          enddo
        else
          do j=lb1,ub1
            do i=lb2,ub2
              fsechist(ix)%data(i,j,idx) = f(j,i)
            enddo
          enddo
        endif
      endif
!
! Set other structure components:
      fsechist(ix)%long_name  = long_name
      fsechist(ix)%units      = units
      fsechist(ix)%geo = geo
      fsechist(ix)%mag = mag
!
! Report to stdout:
      if (firstcall) then
        write(6,"(/,'Initialized diagnostic secondary history field ',
     |    a,' (ix=',i3,'):')") trim(fsechist(ix)%short_name),ix
        write(6,"('  short_name = ',a)") trim(fsechist(ix)%short_name)
        write(6,"('  long_name  = ',a)") trim(fsechist(ix)%long_name)
        write(6,"('  units      = ',a)") trim(fsechist(ix)%units) 
        write(6,"('  geo        = ',l1)")  fsechist(ix)%geo
        write(6,"('  mag        = ',l1)")  fsechist(ix)%mag
        write(6,"('  dimnames   = ',3a)")  fsechist(ix)%dimnames
        write(6,"('  dimsizes   = ',3i4)") fsechist(ix)%dimsizes
        write(6,"('  ndims      = ', i2)") fsechist(ix)%ndims
        write(6,"('  task0_only = ', l5)") fsechist(ix)%task0_only
      endif ! first call for this field
!
! Print non-spval min,max, i.e., portion contributed by current task:
!     if (fsechist(ix)%ndims==3) then
!       call fminmaxspv(fsechist(ix)%data(:,:,:),
!    |    fsechist(ix)%dimsizes(1)*fsechist(ix)%dimsizes(2)*
!    |    fsechist(ix)%dimsizes(3),fmin,fmax,spval)
!       write(6,"('  3-d min,max (non-spval portion so far)= ',
!    |    2e12.4)") fmin,fmax
!     else
!       call fminmaxspv(fsechist(ix)%data(:,:,1),
!    |    fsechist(ix)%dimsizes(1)*fsechist(ix)%dimsizes(2),
!    |    fmin,fmax,spval)
!       write(6,"('  2-d min,max (non-spval portion so far)= ',
!    |    2e12.4)") fmin,fmax
!     endif
      end subroutine addfld
!-----------------------------------------------------------------------
      end module addfld_module
