module ionden_module

  use params_module,only: rp

  implicit none

  contains
!-----------------------------------------------------------------------
  elemental subroutine calculate_prodloss( &
    xnmbar,o2,o1,n2,n2d,n4s,no, &
    ne,op,o2p,np,n2p,op2p,op2d, &
    rk1,rk2,rk3,rk19,rk20,rk25,ra1,ra2,ra3,beta9, &
    qop,qo2p,qnp,qn2p,qnop,qop2p,qop2d, &
    op_prod,o2p_prod,np_prod,n2p_prod,nop_prod,op2p_prod,op2d_prod, &
    op_loss,o2p_loss,np_loss,n2p_loss,nop_loss,op2p_loss,op2d_loss)
! calculate production and loss rates of ions

    use cons_module,only:rmass_o2,rmass_o1, &
      rmass_n2,rmass_n2d,rmass_n4s,rmass_no
    use chemrates_module,only:rk4,rk5,rk6,rk7,rk8,rk9, &
      rk10,rk16,rk17,rk18,rk21,rk22,rk23,rk24,rk26,rk27

    real(kind=rp),intent(in) :: &
      xnmbar,o2,o1,n2,n2d,n4s,no, &
      ne,op,o2p,np,n2p,op2p,op2d, &
      rk1,rk2,rk3,rk19,rk20,rk25,ra1,ra2,ra3,beta9, &
      qop,qo2p,qnp,qn2p,qnop,qop2p,qop2d
    real(kind=rp),intent(out) :: &
      op_prod,o2p_prod,np_prod,n2p_prod,nop_prod,op2p_prod,op2d_prod, &
      op_loss,o2p_loss,np_loss,n2p_loss,nop_loss,op2p_loss,op2d_loss

    real(kind=rp),parameter :: o1min = 1e6_rp
    real(kind=rp) :: o2_cm3,o1_cm3,n2_cm3,n2d_cm3,n4s_cm3,no_cm3

    o2_cm3 = xnmbar*o2/rmass_o2
    o1_cm3 = xnmbar*o1/rmass_o1
    n2_cm3 = xnmbar*n2/rmass_n2
    n2d_cm3 = xnmbar*n2d/rmass_n2d
    n4s_cm3 = xnmbar*n4s/rmass_n4s
    no_cm3 = xnmbar*no/rmass_no

    op_prod = qop + rk8*o1_cm3*np + &
      rk18*o1_cm3*op2p + rk19*op2p*ne + rk21*op2p + &
      rk24*o1_cm3*op2d + rk25*op2d*ne + rk27*op2d
    op_loss = rk1*o2_cm3 + rk2*n2_cm3 + rk10*n2d_cm3

    o2p_prod = qo2p + &
      rk1*o2_cm3*op  + rk6 *o2_cm3*np + &
      rk9*o2_cm3*n2p + rk26*o2_cm3*op2d
    o2p_loss = rk4*n4s_cm3 + rk5*no_cm3 + ra2*ne

! 2024/01 Haonan Wu: With upward extension, around Z=10 near the winter pole,
! np_loss becomes too small, leading to N+ explosion.
! Cap minimum O number density here to prevent N+ explosion.
    np_prod = qnp + rk10*n2d_cm3*op + rk17*n2_cm3*op2p
    np_loss = (rk6+rk7)*o2_cm3 + rk8*max(o1_cm3,o1min)

    n2p_prod = qn2p + rk16*n2_cm3*op2p + rk23*n2_cm3*op2d
    n2p_loss = rk3*o1_cm3 + rk9*o2_cm3 + ra3*ne

    nop_prod = qnop + &
      rk2*n2_cm3 *op  + rk3*o1_cm3*n2p + &
      rk4*n4s_cm3*o2p + rk5*no_cm3*o2p + &
      rk7*o2_cm3 *np  + beta9*no_cm3
    nop_loss = ra1*ne

    op2p_prod = qop2p
    op2p_loss = (rk16+rk17)*n2_cm3 + rk18*o1_cm3 + (rk19+rk20)*ne + rk21 + rk22

    op2d_prod = qop2d + rk20*op2p*ne + rk22*op2p
    op2d_loss = rk23*n2_cm3 + rk24*o1_cm3 + rk25*ne + rk26*o2_cm3 + rk27

  endsubroutine calculate_prodloss
!-----------------------------------------------------------------------
  elemental function ionden_local(step,ionden,prod,loss) result(ionden_upd)
! solve for dn/dt = Q - L*n (no transport)

    real(kind=rp),intent(in) :: step,ionden,prod,loss
    real(kind=rp) :: ionden_upd

    ionden_upd = ionden*exp(-loss*step) + prod/loss*(1-exp(-loss*step))

  endfunction ionden_local
!-----------------------------------------------------------------------
endmodule ionden_module
