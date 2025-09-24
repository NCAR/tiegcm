module emission_rates

    implicit none
    contains
    subroutine ver6300(ne,o1,te,n2, &
        lev0,lev1,lon0,lon1,lat0,lat1)

        ! volume emission rate for O(1D) at 630.0 nm
        ! based on eq (65) of Rees and Roble (1975, Rev. of Geophys. and Space Phys.).
        AD = 9.1e-3 ! [/s], Einstein level of O1D.
        Al = 6.9e-3 ! [/c], transition coef for O1D.
        SD = 5e-11  ! [cm^3/s], rate coef of N2, the quenching species.
        e1D = 1.96  ! [eV], energy required to excite O(1D).
        Le_O1D = 1.07e-10*NE.*O.*sqrt(Te).*exp(-2.27e4./Te) &
           .*( 0.406 + 0.357e-4.*Te - (0.333+0.183e-4.*Te).*exp(-1.37e4./Te) - (0.456+0.174e-4.*Te).*exp(-2.97e4./Te) ) ! [eV/cm^3/s]
        eta_1D = Le_O1D./e1D ! [/cm^3/s]
        ver = Al./AD.*eta_1D./(1+SD./AD.*N2) ! [Photons/cm^3/s]

    end subroutine ver6300

end module emission_rates
