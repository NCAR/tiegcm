module emission_module
! volume emission rates for O(1D) at 6300 A and N(2D) at 5200 A
! based on Rees, M. H., and Roble, R. G. (1975)
! Observations and theory of the formation of stable auroral red arcs
! Reviews of Geophysics, 13(1), 201–242, doi:10.1029/RG013i001p00201

  use params_module,only:rp

  implicit none

  contains
!-----------------------------------------------------------------------
  elemental function ver6300(o1_cm3,n2_cm3,ne,te) result(ver)

    real(kind=rp),intent(in) :: o1_cm3,n2_cm3,ne,te
    real(kind=rp) :: ver

    real(kind=rp),parameter :: &
      AD = 9.1e-3_rp, & ! Einstein level of O(1D) [s^-1]
      Al = 6.9e-3_rp, & ! transition coef for O(1D) [s^-1]
      SD = 5e-11_rp, &  ! rate coef of quenching species N2 [cm^3/s]
      e1D = 1.96_rp     ! energy required to excite O(1D) [eV]
    real(kind=rp) :: Le_O1D,eta_1D

! Eq (61) [eV/cm^3/s]
    Le_O1D = 1.07e-10_rp*ne*o1_cm3*sqrt(te)*exp(-2.27e4_rp/te)* &
      (0.406_rp + 0.357e-4_rp*te - &
       (0.333_rp + 0.183e-4_rp*te)*exp(-1.37e4_rp/te) - &
       (0.456_rp + 0.174e-4_rp*te)*exp(-2.97e4_rp/te))

! Eq (64) [cm^-3/s]
    eta_1D = Le_O1D/e1D

! Eq (65) [cm^-3/s]
    ver = Al*eta_1D/(AD+SD*n2_cm3)

  endfunction ver6300
!-----------------------------------------------------------------------
  elemental function ver5200(o2_cm3,n2d_cm3,ne,te) result(ver)

    real(kind=rp),intent(in) :: o2_cm3,n2d_cm3,ne,te
    real(kind=rp) :: ver

    real(kind=rp),parameter :: &
      AD = 1.06e-5_rp, & ! Einstein level of N(2D) [s^-1]
      e2D = 2.37_rp      ! energy required to excite N(2D) [eV]
    real(kind=rp) :: omega,Le_N2D,eta_2D,SD

! Eq (4) in Henry, R. J. W., and Williams, R. E. (1968)
! Collision strengths and photoionization cross sections for nitrogen, oxygen, and neon
! Astronomical Society of the Pacific, 80(477), 669, doi:10.1086/128709

! 2025/09 Haonan Wu:
! Rees and Roble (1975) refers to Table 1 in Henry and Williams (1968),
! while Table 1 is the temperature dependence of gamma, the energy integral of omega.
! The energy dependence of omega is given in Figure 1 without an explicit expression.
! From Figure 1, N(1,3) does not seem to change much with energy,
! therefore I will simply use a constant omega of 0.5 for the calculation
    omega = 0.5_rp

! Eq (63) [eV/cm^3/s]
    Le_N2D = 3.4e-18_rp/sqrt(te)*exp(-2.74e4_rp/te)*omega*ne*n2d_cm3

! Eq (66) [cm^-3/s]
    eta_2D = Le_N2D/e2D

! Eq (67) [s^-1]
    SD = 1.0e-12_rp*o2_cm3 + 8.54e-7_rp/sqrt(te)*omega*ne

! Eq (68) [cm^-3/s]
    ver = eta_2D*AD/(AD+SD)

  endfunction ver5200
!-----------------------------------------------------------------------
endmodule emission_module
