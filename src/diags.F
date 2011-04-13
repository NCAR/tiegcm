!
      module diags_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use addfld_module,only: addfld
      implicit none
!
! March, 2011 btf: add "sanctioned" diagnostics.
! Field name lengths (s.a., fields.F):
      integer,parameter :: 
     |  longname_len  = 80, ! length of field long name 
     |  shortname_len = 16, ! length of field short name 
     |  units_len     = 16  ! length of field units attribute
!
! number of diagnostics (see sub init_diags):
      integer,parameter :: ndiag = 14
      type diag_type
        character(len=longname_len)  :: long_name
        character(len=shortname_len) :: short_name
        character(len=units_len)     :: units
        character(len=shortname_len) :: levels ! "lev", "ilev", or "none"
      end type diag_type
      type(diag_type) :: diags(ndiag)

      contains
!-----------------------------------------------------------------------
      subroutine init_diags(iprint)
      implicit none
!
! Args:
      integer,intent(in) :: iprint ! if > 0, print diags list to stdout
!
! Init to empty strings:
!
      integer :: n
      do n=1,ndiag
        diags(n)%long_name  = ' '
        diags(n)%short_name = ' '
        diags(n)%units      = ' '
        diags(n)%levels     = ' '
      enddo
!
! Fields order is alphabetical by short_name:
!
      n = 1
      diags(n)%short_name = 'COOL_EXP'
      diags(n)%long_name  = 'Total Explicit Cooling'
      diags(n)%units      = 'erg/K/s'  ! Table 4.3 of Model Description
      diags(n)%levels     = 'lev'
!
      n = n+1
      diags(n)%short_name = 'COOL_IMP'
      diags(n)%long_name  = 'Total Implicit Cooling'
      diags(n)%units      = '1/s'       ! Table 4.3 of Model Description
      diags(n)%levels     = 'lev'
!
      n = n+1
      diags(n)%short_name = 'DEN'
      diags(n)%long_name  = 'Total Density'
      diags(n)%units      = 'g/cm3'
      diags(n)%levels     = 'ilev'
!
      n = n+1
      diags(n)%short_name = 'HEATING'
      diags(n)%long_name  = 'Total Heating'
      diags(n)%units      = 'erg/g/s'
      diags(n)%levels     = 'lev'
!
      n = n+1
      diags(n)%short_name = 'HMF2'
      diags(n)%long_name  = 'Height of F2' 
      diags(n)%units      = 'km'
      diags(n)%levels     = 'none' ! hmf2 is 2d lon x lat
!
      n = n+1
      diags(n)%short_name = 'NMF2'
      diags(n)%long_name  = 'Density at HMF2' 
      diags(n)%units      = 'cm3'
      diags(n)%levels     = 'none' ! nmf2 is 2d lon x lat
!
! Units of joule heating are erg/g/s (qjoule.F and p.141 of MD, 
! however Table 4.3 of the MD says qjoule input to dt.F is in erg/K/s)
! I am going w/ erg/g/s to be compatable w/ total heating.
      n = n+1
      diags(n)%short_name = 'QJI_TI'
      diags(n)%long_name  = 'Joule Heating from TI' 
      diags(n)%units      = 'erg/g/s'   
      diags(n)%levels     = 'lev'
!
      n = n+1
      diags(n)%short_name = 'QJI_TN'
      diags(n)%long_name  = 'Joule Heating from TN'
      diags(n)%units      = 'erg/g/s'
      diags(n)%levels     = 'lev'
!
      n = n+1
      diags(n)%short_name = 'SIGMA_HAL'
      diags(n)%long_name  = 'Hall Conductivity'
      diags(n)%units      = 'S/m'
      diags(n)%levels     = 'lev'
!
      n = n+1
      diags(n)%short_name = 'SIGMA_PED'
      diags(n)%long_name  = 'Pedersen Conductivity'
      diags(n)%units      = 'S/m'
      diags(n)%levels     = 'lev'
!
      n = n+1
      diags(n)%short_name = 'TEC'
      diags(n)%long_name  = 'Total Electron Content'
      diags(n)%units      = 'TEC'
      diags(n)%levels     = 'none' ! 2d lon x lat
!
      n = n+1
      diags(n)%short_name = 'UI_ExB'
      diags(n)%long_name  = 'Zonal Ion Drift (ExB)'
      diags(n)%units      = 'm/s'
      diags(n)%levels     = 'ilev'
!
      n = n+1
      diags(n)%short_name = 'VI_ExB'
      diags(n)%long_name  = 'Meridional Ion Drift (ExB)'
      diags(n)%units      = 'm/s'
      diags(n)%levels     = 'ilev'
!
      n = n+1
      diags(n)%short_name = 'WI_ExB'
      diags(n)%long_name  = 'Vertical Ion Drift (ExB)'
      diags(n)%units      = 'm/s'
      diags(n)%levels     = 'ilev'

      if (n /= ndiag) then
        write(6,"('>>> init_diags: n=',i4,' ndiag=',i4)")
     |    n,ndiag
        call shutdown('init_diags ndiag')
      endif
!
! Print table to stdout:
      if (iprint > 0) then
        write(6,"(/,72('-'))")
        write(6,"('List of Available Diagnostic Fields:')")
        write(6,"('Shortnames may be added to namelist SECFLDS',/)")
        write(6,"('Field   Shortname       Units           Levels  ',
     |    'Longname')")
        do n=1,ndiag
          write(6,"(i4,4x,a,a,a,a)") n,diags(n)%short_name,
     |      diags(n)%units,diags(n)%levels(1:8),trim(diags(n)%long_name)
        enddo
        write(6,"(72('-'),/)")
      endif
      end subroutine init_diags
!-----------------------------------------------------------------------
      subroutine mkdiag_DEN(name,rho,lev0,lev1,lon0,lon1,lat)
!
! Save diagnostic DEN (Total Density). This routine is called from 
! dt.F, and simply passes rho to addfld at the current latitude.
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in) :: rho(lev0:lev1,lon0:lon1)
!
! Local:
      integer :: ix
!
! Check that field name is a diagnostic, and was requested:
      ix = checkf(name) ; if (ix==0) return
!
! Save to secondary history. Note rho is calculated on interface
! levels, hence 'ilev' in the addfld call:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,rho,'ilev',lev0,lev1,'lon',lon0,lon1,lat)
      end subroutine mkdiag_DEN
!-----------------------------------------------------------------------
      subroutine mkdiag_QJI_TI(name,qji_ti,lev0,lev1,lon0,lon1,lat)
!
! From qjoule.F, levels at midpoints.
!
! Args
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in) :: qji_ti(lev0:lev1,lon0-2:lon1+2)
!
! Local:
      integer :: ix
!
! Check that field name is a diagnostic, and was requested:
      ix = checkf(name) ; if (ix==0) return
!
! Save on secondary history:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,qji_ti(lev0:lev1-1,lon0:lon1),diags(ix)%levels,
     |  lev0,lev1-1,'lon',lon0,lon1,lat)
      end subroutine mkdiag_QJI_TI
!-----------------------------------------------------------------------
      subroutine mkdiag_QJI_TN(name,qji_tn,lev0,lev1,lon0,lon1,lat)
!
! From qjoule.F, levels at midpoints.
!
! Args
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in) :: qji_tn(lev0:lev1,lon0-2:lon1+2)
!
! Local:
      integer :: ix
!
! Check that field name is a diagnostic, and was requested:
      ix = checkf(name) ; if (ix==0) return
!
! Save on secondary history:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,qji_tn(lev0:lev1-1,lon0:lon1),diags(ix)%levels,
     |  lev0,lev1-1,'lon',lon0,lon1,lat)

      end subroutine mkdiag_QJI_TN
!-----------------------------------------------------------------------
      subroutine mkdiag_TEC(name,tec,z,electrons,lev0,lev1,lon0,lon1,
     |  lat)
!
! Calculate 2d global height-integrated diagnostic TEC (Total Electron Content)
! This is called from sub elden, elden.F at latitude lat.
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: z,electrons
      real,intent(out) :: tec(lon0:lon1)
!
! Local:
      integer :: i,k,ix
!
! Check that field name is a diagnostic, and was requested:
      ix = checkf(name) ; if (ix==0) return
!
! Integrate electron content in height at current latitude:
      tec(:) = 0.
      do i=lon0,lon1
        do k=lev0,lev1-1
          tec(i) = tec(i)+(z(k+1,i)-z(k,i))*electrons(k,i)
        enddo
      enddo
!     write(6,"('mkdiag_TEC call addfld: lat=',i4,' ',a,' min,max=',
!    |  2e12.4)") lat,name,minval(tec),maxval(tec)
!
! Save on secondary histories:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,tec(:),'lon',lon0,lon1,'lat',lat,lat,0)
      end subroutine mkdiag_TEC
!-----------------------------------------------------------------------
      subroutine mkdiag_UI(name,ui,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Save 3d diagnostic UI_ExB. This is called from ionvel.F. 
! This routine simply passes ui from ionvel to addfld.
! Levels at interfaces. Incoming units are cm/s, but is put
! on the history in m/s.
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,intent(in) :: ui(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2)
!
! Local:
      integer :: lat,ix
!
! Check that field name is a diagnostic, and was requested:
      ix = checkf(name) ; if (ix==0) return
!
! Save on secondary histories (convert from cm/s to m/s):
      do lat=lat0,lat1
        call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |    diags(ix)%units, ui(:,lon0:lon1,lat)/100.,diags(ix)%levels,
     |    lev0,lev1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      end subroutine mkdiag_UI
!-----------------------------------------------------------------------
      subroutine mkdiag_VI(name,vi,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Save 3d diagnostic VI_ExB. This is called from ionvel.F. 
! This routine simply passes vi from ionvel to addfld.
! Levels at interfaces. Incoming units are cm/s, but is put
! on the history in m/s.
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,intent(in) :: vi(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2)
!
! Local:
      integer :: lat,ix
!
! Check that field name is a diagnostic, and was requested:
      ix = checkf(name) ; if (ix==0) return
!
! Save on secondary histories (convert from cm/s to m/s):
      do lat=lat0,lat1
        call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |    diags(ix)%units, vi(:,lon0:lon1,lat)/100.,diags(ix)%levels,
     |    lev0,lev1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      end subroutine mkdiag_VI
!-----------------------------------------------------------------------
      subroutine mkdiag_WI(name,wi,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Save 3d diagnostic WI_ExB. This is called from ionvel.F. 
! This routine simply passes wi from ionvel to addfld.
! Levels at interfaces. Incoming units are cm/s, but is put
! on the history in m/s.
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,intent(in) :: wi(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2)
!
! Local:
      integer :: lat,ix
!
! Check that field name is a diagnostic, and was requested:
      ix = checkf(name) ; if (ix==0) return
!
! Save on secondary histories (convert from cm/s to m/s):
      do lat=lat0,lat1
        call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |    diags(ix)%units, wi(:,lon0:lon1,lat)/100.,diags(ix)%levels,
     |    lev0,lev1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      end subroutine mkdiag_WI
!-----------------------------------------------------------------------
      subroutine mkdiag_HEAT(name,total_heat,lev0,lev1,lon0,lon1,lat)
!
! Total heating from dt.F. Levels on midpoints.
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in) :: total_heat(lev0:lev1,lon0:lon1)
!
! Local:
      integer :: ix
!
! Check that field name is a diagnostic, and was requested:
      ix = checkf(name) ; if (ix==0) return
!
! Save on secondary histories:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,total_heat(lev0:lev1-1,:),diags(ix)%levels,
     |  lev0,lev1-1,'lon',lon0,lon1,lat)
      end subroutine mkdiag_HEAT
!-----------------------------------------------------------------------
      subroutine mkdiag_COOLEXP(name,cool_exp,lev0,lev1,lon0,lon1,lat)
!
! Explicit cooling from newton.F, levels="lev" (midpoints)
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in) :: cool_exp(lev0:lev1,lon0-2:lon1+2)
!
! Local:
      integer :: ix
!
! Check that field name is a diagnostic, and was requested:
      ix = checkf(name) ; if (ix==0) return
!
! Save on secondary histories:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,cool_exp(:,lon0:lon1),diags(ix)%levels,
     |  lev0,lev1,'lon',lon0,lon1,lat)
      end subroutine mkdiag_COOLEXP
!-----------------------------------------------------------------------
      subroutine mkdiag_COOLIMP(name,cool_imp,lev0,lev1,lon0,lon1,lat)
!
! Implicit cooling from newton.F, levels="lev" (midpoints)
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in) :: cool_imp(lev0:lev1,lon0-2:lon1+2)
!
! Local:
      integer :: ix
!
! Check that field name is a diagnostic, and was requested:
      ix = checkf(name) ; if (ix==0) return
!
! Save on secondary histories:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,cool_imp(:,lon0:lon1),diags(ix)%levels,
     |  lev0,lev1,'lon',lon0,lon1,lat)
      end subroutine mkdiag_COOLIMP
!-----------------------------------------------------------------------
      subroutine mkdiag_SIGMAPED(name,sigmaped,lev0,lev1,lon0,lon1,lat)
!
! From lamdas.F. Levels at midpoints (altho lamdas in lamdas.F are
! at interfaces)
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in) :: sigmaped(lev0:lev1,lon0-2:lon1+2)
!
! Local:
      integer :: ix
!
! Check that field name is a diagnostic, and was requested:
      ix = checkf(name) ; if (ix==0) return
!
! Save on secondary histories:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,sigmaped(lev0:lev1-1,lon0:lon1),
     |  diags(ix)%levels,lev0,lev1-1,'lon',lon0,lon1,lat)
      end subroutine mkdiag_SIGMAPED
!-----------------------------------------------------------------------
      subroutine mkdiag_SIGMAHAL(name,sigmahal,lev0,lev1,lon0,lon1,lat)
!
! From lamdas.F. Levels at midpoints (altho lamdas themselves are at
! interfaces)
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in) :: sigmahal(lev0:lev1,lon0-2:lon1+2)
!
! Local:
      integer :: ix
!
! Check that field name is a diagnostic, and was requested:
      ix = checkf(name) ; if (ix==0) return
!
! Save on secondary histories:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,sigmahal(lev0:lev1-1,lon0:lon1),
     |  diags(ix)%levels,lev0,lev1-1,'lon',lon0,lon1,lat)
      end subroutine mkdiag_SIGMAHAL
!-----------------------------------------------------------------------
      subroutine mkdiag_HNMF2(name,ht,ne,lev0,lev1,lon0,lon1,lat)
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in),dimension(lev0:lev1,lon0:lon1) :: ht,ne
!
! Local:
      integer :: ix
      real,dimension(lon0:lon1) :: hmf2,nmf2
!
! Check that field name is a diagnostic, and was requested:
      ix = checkf(name) ; if (ix==0) return
!
! Calculate hmf2,nmf2 from ne at current latitude:
      call hnmf2(ht,ne,hmf2,nmf2,lev0,lev1,lon0,lon1,lat) 
!
! Save 2d field on secondary histories:
      if (trim(name)=='HMF2') then
        hmf2 = hmf2*1.e-5 ! cm to km
        call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |    diags(ix)%units,hmf2,'lon',lon0,lon1,'lat',lat,lat,0)
      elseif (trim(name)=='NMF2') then
        call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |    diags(ix)%units,nmf2,'lon',lon0,lon1,'lat',lat,lat,0)
      endif
      end subroutine mkdiag_HNMF2
!-----------------------------------------------------------------------
      integer function checkf(name)
      use input_module,only: secflds,mxfsech
!
! Args:
      character(len=*),intent(in) :: name
!
! External:
      integer,external :: ixfindc ! util.F
!
! Check for field in namelist SECFLDS -- if not found, return silently
      checkf = ixfindc(secflds,mxfsech,name)
      if (checkf==0) return
!
! Get index to diags -- if not found, return:
      checkf = ixfindc(diags(:)%short_name,ndiag,name)
      if (checkf==0) then
        write(6,"('>>> Cannot find field ',a,' in list of ',
     |    'available diagnostics -- returning')") name
      endif
      return
      end function checkf
!-----------------------------------------------------------------------
      subroutine hnmf2(ht,ne,hmf2,nmf2,lev0,lev1,lon0,lon1,lat)
!
! Find the NmF2 and hmF2 based on the electron profile.
! btf: This is taken from an idl procedure provided by Liying.
!
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in),dimension(lev0:lev1,lon0:lon1) :: ht,ne
      real,intent(out),dimension(lon0:lon1) :: hmf2,nmf2
!
! Local:
      integer :: k,kx,i
      real :: h(3),p(3),h12,h22,h32,deltx,atx,ax,btx,bx,ctx,cx
!
      do i=lon0,lon1
        kx = 0
        kloop: do k = lev1-1,2,-1
          if (ne(k,i) >= ne(k-1,i).and.ne(k,i) >= ne(k+1,i)) then
            kx = k
            exit kloop
          endif
        enddo kloop
        if (kx==0) then
          write(6,"('>>> hnmf2: lat=',i4,' could not find kx -- ',
     |      'setting hmf2=nmf2=0')") lat
          hmf2(:) = 0. ; nmf2(:) = 0.
          return
        endif
        h = (/ht(kx-1,i),ht(kx,i),ht(kx+1,i)/)
        p = (/ne(kx-1,i),ne(kx,i),ne(kx+1,i)/)
        h12 = h(1)*h(1)
        h22 = h(2)*h(2)
        h32 = h(3)*h(3)
        deltx=h12*h(2)+h22*h(3)+h32*h(1)-h32*h(2)-h12*h(3)-h22*h(1)
        atx=p(1)*h(2)+p(2)*h(3)+p(3)*h(1)-h(2)*p(3)-h(3)*p(1)-h(1)*p(2)
        ax=atx/deltx
        btx=h12*p(2)+h22*p(3)+h32*p(1)-h32*p(2)-h12*p(3)-h22*p(1)
        bx=btx/deltx
        ctx=h12*h(2)*p(3)+h22*h(3)*p(1)+h32*h(1)*p(2)-h32*h(2)*p(1)-
     |    h12*h(3)*p(2)-h22*h(1)*p(3)
        cx=ctx/deltx
        hmf2(i)=-(bx/(2.*ax))
        nmf2(i)=-((bx*bx-4.*ax*cx)/(4.*ax))
      enddo ! i=lon0,lon1
!     write(6,"('hnmf2: lat=',i4,' hmf2 min,max=',2e12.4,
!    |  ' nmf2 min,max=',2e12.4)") lat,minval(hmf2),maxval(hmf2),
!    |  minval(nmf2),maxval(nmf2)
      end subroutine hnmf2
!-----------------------------------------------------------------------
!
! This is Liying's original IDL procedure hnmf2.pro. This was converted
! to fortran in subroutine hnmf2 above.
!
! PRO hnmf2, X, Y, XOUX, YOUX
! ;
! ; Find the NmF2 and hmF2 based on the electron profile
! ;
! ; input
! ;     X: height
! ;     Y: Ne profile
! ;
! ; output
! ;     XOUX: hmF2
! ;     YOUX: NmF2
! ;
! nz=n_elements(x)
! for k=nz-2,1,-1 do begin
!   if ((y[k] ge y[k-1]) and (y[k] ge y[k+1])) then break
! endfor
! h=[x[k-1],x[k],x[k+1]]
! p=[y[k-1],y[k],y[k+1]]
! 
! h12=h(0)*h(0)
! h22=h(1)*h(1)
! h32=h(2)*h(2)
! DELTX=h12*h(1)+h22*h(2)+h32*h(0)-h32*h(1)-h12*h(2)-h22*h(0)
! ATX=p(0)*h(1)+p(1)*h(2)+p(2)*h(0)-h(1)*p(2)-h(2)*p(0)-h(0)*p(1)
! AX=ATX/DELTX
! BTX=h12*p(1)+h22*p(2)+h32*p(0)-h32*p(1)-h12*p(2)-h22*p(0)
! BX=BTX/DELTX
! CTX=h12*h(1)*p(2)+h22*h(2)*p(0)+h32*h(0)*p(1)-h32*h(1)*p(0)-h12*h(2)*p(1)-h22*h(0)*p(2)
! CX=CTX/DELTX
! XOUX=-(BX/(2.*AX))
! YOUX=-((BX*BX-4.*AX*CX)/(4.*AX))
!
! RETURN
! END
!-----------------------------------------------------------------------
      end module diags_module
