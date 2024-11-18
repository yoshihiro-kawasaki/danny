module danny_dust_m_

    use physical_constants
    use danny_variables_m_
    implicit none

contains

subroutine setDustProperty

    use danny_variables_m_
    implicit none

    integer i

    dust_to_gas_mass_ratio_ = 1.0d-2
    dust_internal_density_ = 3.0d0
    dust_grain_radius_ = 1.0d-5
    dust_grain_cross_section_ = PI * dust_grain_radius_**2.0d0
    dust_grain_mass_ = (4.0d0*PI/3.0) * dust_grain_radius_**3.0d0 * dust_internal_density_
    dust_number_denisty_ = (dust_to_gas_mass_ratio_* 1.36*ATOMIC_MASS_UNIT*number_density_) / dust_grain_mass_
    dust_abundances_ = dust_number_denisty_ / number_density_
    dust_temperature_ = temperature_

    binding_sites_per_a_dust_grain_ = 4.0*PI*dust_grain_radius_**2.0d0 * BINDING_SITES_DENSITY_

    do i = 1, number_of_species_
        if (species_name_(i)(1:1) == 'G') then
            species_mass_(i) = dust_grain_mass_
        end if 
    end do

end subroutine setDustProperty

end module danny_dust_m_