!
      subroutine vtsetup
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Define user states and activities for vampir tracing: 
!
      implicit none
      integer :: ier,code
      character(len=16) :: state,activity
#ifdef VT
      code = 100 ; state = 'mp_bndlats' ; activity = 'Communication'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 101 ; state = 'mp_bndlons' ; activity = 'Communication'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 102 ; state = 'mp_bndlons_f3d' ; activity = 'Communication'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
!     code = 103 ; state = 'mp_pole_lat' ; activity = 'Communication'
      code = 103 ; state = 'mp_polelat'  ; activity = 'Communication'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 104 ; state = 'mp_gatherlons' ; activity = 'Communication'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 105 ; state = 'mp_scatterlons' ; activity = 'Communication'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 106 ; state='mp_gatherlons_f3d' ; activity='Communication'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 107 ; state='mp_scatterlons_f3d' ; activity='Communication'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 108 ; state = 'mp_periodic_f4d' ; activity='Communication'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 109 ; state = 'mp_periodic_f3d' ; activity='Communication'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 110 ; state = 'mp_dynpot' ; activity='Communication'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 111 ; state = 'mp_bndlats_kmh' ; activity='Communication'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 112 ; state = 'mp_bndlons_kmh' ; activity='Communication'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 113 ; state = 'oplus' ; activity='ModelCode'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 114 ; state = 'minor' ; activity='ModelCode'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 115 ; state = 'advance' ; activity='ModelCode'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 116 ; state = 'addiag' ; activity='ModelCode'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 117 ; state = 'efield' ; activity='ModelCode'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 118 ; state = 'qrj' ; activity='ModelCode'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 119 ; state = 'hdif1' ; activity='ModelCode'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 120 ; state = 'hdif2' ; activity='ModelCode'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 121 ; state = 'ionvel' ; activity='ModelCode'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
! State "ions" includes subs qrj, qinite, xray, and aurora.
! (qrj also has its own separate state, see above)
      code = 122 ; state = 'ions' ; activity='ModelCode'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
! Filtering activity from subs filter_minor (minor.F), filter_op (oplus.F),
! and filter_w (swdot.F):
      code = 123 ; state = 'filter_minor' ; activity='Filtering'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
      code = 124 ; state = 'filter_op' ; activity='Filtering'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
      code = 125 ; state = 'filter_w' ; activity='Filtering'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
! Resume "ModelCode" activity:
      code = 126 ; state = 'settei' ; activity='ModelCode'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 127 ; state = 'timestep' ; activity='ModelCode'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 128 ; state = 'qinite' ; activity='ModelCode'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 129 ; state = 'xray' ; activity='ModelCode'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)
!
      code = 130 ; state = 'aurora' ; activity='ModelCode'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)

      code = 131 ; state = 'filter_tn' ; activity='Filtering'
      call vtsymdef(code, state, activity, ier)
      write(6,"('vttrace: Code: ',i3,' State: ',a,' Activity: ',a)")
     |  code,trim(state),trim(activity)

#endif
      end subroutine vtsetup
