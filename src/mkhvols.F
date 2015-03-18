!
      module mk_hvols
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Enable history mss path file name expansion, e.g.:
! "/USER/file1","to","/USER/file5","by","1" becomes
! "/USER/file1","/USER/file2",...,"/USER/file5" 
! See function mkhvols. This is called from input.
!
      contains
!-----------------------------------------------------------------------
      logical function isdigit(s)
      implicit none
      character(len=1),intent(in) :: s
      character(len=1) :: digits(10)=
     +  (/'0','1','2','3','4','5','6','7','8','9'/)
      integer :: i
      isdigit = .false.
      do i=1,10
        if (s==digits(i)) then
          isdigit = .true.
          exit
        endif
      enddo
      end function isdigit
!-----------------------------------------------------------------------
      function ndollars(s)
      character(len=*),intent(in) :: s
      integer :: i
      ndollars = 0
      do i=1,len(s)
        if (s(i:i)=='$') ndollars = ndollars+1 
      enddo
      end function ndollars
!-----------------------------------------------------------------------
      integer function decode_int(str,sout)
      implicit none
!
! Return integer decoded from digits in string str.
! Only use digits from end of str. Also return string
! sout, with integers
!
      character(len=*),intent(in) :: str
      character(len=*),intent(out) :: sout
      character(len=len(str)) :: snum,s
      character(len=8) :: format
      integer :: i,ipos,slen
!
      decode_int = -1
      s = adjustl(str)
      slen = len_trim(s)
      if (slen>9999) then
        write(6,"('>>> decode_int: string too long:',
     +    ' slen=',i10)") slen 
        return
      endif
      ipos = 0
      snum = ' '
      sout = trim(s)
      do i=slen,1,-1
        if (.not.isdigit(s(i:i))) then
          if (ipos > 0) then
            exit
          else
            cycle
          endif
        endif
        ipos = ipos+1
        snum(slen-ipos+1:slen-ipos+1) = s(i:i)
        sout(i:i) = '$'
      enddo
      if (ipos==0) return
      snum = adjustl(snum)
      write(format,"('(i',i4,')')") len_trim(snum)
!     read(trim(snum),format,err=100) decode_int ! ibm does not like this
      read(snum,format,err=100) decode_int
      return
 100  write(6,"('>>> decode_int: error making integer from ',
     +  'string: ',a)") trim(snum)
      decode_int = -1
      end function decode_int
!-----------------------------------------------------------------------
      subroutine encode_str(str,num,ndigits)
      implicit none
!
! Given str containing n characters '$', return
! same string, but with the '$' replaced by num.
!
      character(len=*),intent(inout) :: str
      character(len=len(str)) :: newstr
      integer,intent(in) :: num,ndigits
      integer :: i,ipos,ndollar
      character(len=16) :: format,numstr
!
      if (ndigits<=0) return
      write(format,"('(i',i4,'.',i4,')')") ndigits,ndigits
      format = trim(adjustl(format))
      ipos = index(str,'$')
      ndollar = ndollars(str)
      write(numstr,format,err=100) num
      call replace_substr(str,str(ipos:ipos+ndollar-1),numstr,newstr)
      str = newstr
      return
 100  write(6,"('>>> encode_str: error writing integer ',
     +  'num=',i12,' to string.')") num
      return
      end subroutine encode_str
!-----------------------------------------------------------------------
      subroutine replace_substr(s,subold,subnew,sout)
      character(len=*),intent(in) :: s,subold,subnew
      character(len=*),intent(out) :: sout
      integer :: i,ii,ipos,lonsubnew,lensubold,lens
!
! Init sout to input s:
      sout = s
!
! Get position of first char of subold in s:
      ipos = index(s,subold)
      if (ipos <= 0) then
        write(6,"('>>> WARNING replace_substr: could not find subold ',
     |    a,' in string ',a)") trim(subold),trim(s)
        return
      endif
!
! Check subnew:
      lensubnew = len_trim(subnew)
      if (lensubnew <= 0) then
        write(6,"('>>> WARNING replace_substr: zero length subnew.')")
        return
      endif
!
! Check if enough space in sout:
      lensubold = len_trim(subold)
      lens = len_trim(s)
      if (lens-lensubold+lensubnew > len(sout)) then
        write(6,"('>>> WARNING replace_substr: sout not long enough:',
     |    ' lens-lensubold+lensubnew=',i4,' len(sout)=',i4)")
     |    lens-lensubold+lensubnew,len(sout)
        return
      endif
!
! Replace subold with subnew in sout:
      sout(ipos:ipos+lensubnew-1) = trim(subnew)
!
! Complete sout with remainder of s:
      ii = ipos+lensubnew
      do i=ipos+lensubold,lens
        sout(ii:ii) = s(i:i)  
        ii = ii+1
      enddo
      end subroutine replace_substr
!-----------------------------------------------------------------------
      integer function mkhvols(vols_in,vols_out,mxout)
      implicit none
!
! Given array of strings vols_in as mss paths, return 
!   list of vols_out, which are same as vols_in unless
!   vols_in(2) = 'to', in which case vols_in is, e.g.:
!   "vol5","to","vol10","by","n". In this case expand
!   to "vol5","vol6",...,"vol10" (use decode_int)
! Return number of defined elements in vols_out
!
! Args:
      character(len=*),intent(in) :: vols_in(:)
      integer,intent(in) :: mxout
      character(len=*),intent(out) :: vols_out(mxout)
!
! Locals:
      character(len=240) :: vol_first,vol_last,template_first,
     +  template_last,template_dum
      integer :: nvin,i,n_first,n_last,n_delta,ipos0,ipos1,
     |  ndol_first,ndigits,ndig
!
! Check input sizes:
      mkhvols = 0
      nvin = size(vols_in)
      if (nvin > mxout) then
        write(6,"('>>> mkhvols: too many volume inputs =',
     +    i3,' (must be <= mxout=',i3)") nvin,mxout
        nvin = mxout
      elseif (nvin<=0) then
        write(6,"('>>> mkhvols: need vols_in')")
        return
      endif
      vols_out = ' '
      if (nvin==1) then
        vols_out(1) = vols_in(1)
        return
      endif
!
! No expansion -> just echo non-blank elements of vols_in to vols_out:
!
      if (trim(vols_in(2))/='to'.and.trim(vols_in(2))/='TO') then
        do i=1,nvin
          if (len_trim(vols_in(i)) > 0) then
            mkhvols = mkhvols+1
            vols_out(mkhvols) = vols_in(i)
          endif  
        enddo
        return
      endif
!
! Do expansion of "vols1","to","vols2","by","n"
!
      write(6,"('Doing expansion of mss histvol paths:')")
      write(6,"(a,' ',a,' ',a,' ',a,' ',a)")
     +  trim(vols_in(1)),trim(vols_in(2)),trim(vols_in(3)),
     +  trim(vols_in(4)),trim(vols_in(5))
!
      if (trim(vols_in(4))/='by'.and.trim(vols_in(4))/='BY') then
        write(6,"('>>> mkhvols: if 2nd vols is ""',a,'"" then ',
     +    '4th vol must be ""by"" or ""BY""')") trim(vols_in(2))
        return
      endif
!
! vol_first is first volume, n_first is number from first volume:
      vol_first = vols_in(1)
      vol_first = trim(vol_first)
      n_first = decode_int(vol_first,template_first)
      if (n_first==-1) then
        write(6,"('>>> mkhvols: could not get number from ',
     +    'first volume: vols_in(1)=',a)") vols_in(1)
        return
      endif
      ndol_first = ndollars(template_first)
!
! vol_last is last volume, n_last is number from last volume:
      vol_last  = vols_in(3)
      vol_last = trim(vol_last)
      n_last = decode_int(vol_last,template_last)
      if (n_last==-1) then
        write(6,"('>>> mkhvols: could not get number from ',
     +    'last volume: vols_in(3)=',a)") vols_in(3)
        return
      endif
!
! n_delta is delta integer between volumes:
      n_delta = decode_int(vols_in(5),template_dum)
      if (n_delta==-1) then
        write(6,"('>>> mkhvols: could not get delta from ',
     +    'vols_in(5)=',a)") vols_in(5)
        return 
      endif
!
! Check delta:
      if (n_first > n_last) then
        write(6,"('>>> mkhvols: n_first must be <= n_last:',
     +    ' n_first=',i5,' n_last=',i5)") n_first,n_last
        return
      endif 
      if (n_delta > n_last-n_first) then
        write(6,"('>>> mkhvols: n_delta must be <= n_last-n_first:',
     +    ' n_delta=',i5,' n_first=',i5,' n_last=',i5)")
     +    n_delta,n_first,n_last
      endif
!
! Warn if first template /= last template (i.e., root names are
! different, not including number of dollars). (Not fatal -- will 
! use 1st template in this case)
!
      ipos0 = index(template_first,'$')
      ipos1 = index(template_last,'$')
      if (template_first(1:ipos0-1)/=template_last(1:ipos1-1)) then
        write(6,"('>>> WARNING mkhvols: root names of first and ',
     +    'last mss files are different:',/'  template_first=',a,
     +    ' template_last=',a,'  (will use first)')")
     +    trim(template_first),trim(template_last) 
      endif
      ndigits = ndol_first
!
! Loop from number in first vol to number in last vol, by delta,
! building vols_out w/ encode_str.
!
      do i=n_first,n_last,n_delta
        ndig = int(log10(real(i)))+1	! number of digits in i
        if (ndig > ndigits) ndigits = ndigits+(ndig-ndigits)
        if (mkhvols+1 > mxout) then
          write(6,"('>>> WARNING mkhvols: too many vols -- ',
     +      'stopping at mxhvols=',i3)") mkhvols
          exit
        endif
        mkhvols = mkhvols+1
        vols_out(mkhvols) = template_first
        call encode_str(vols_out(mkhvols),i,ndigits)

!       write(6,"('mkhvols: i=',i3,' ndigits=',i2,' mkhvols=',i3,
!    |    ' template=',a,' vols_out=',a)") i,ndigits,mkhvols,
!    |    trim(template_first),trim(vols_out(mkhvols))

      enddo
!     call shutdown('mkhvols')
      end function mkhvols
      end module mk_hvols
