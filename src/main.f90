program main

    use danny_variables_m_
    use danny_read_m_
    use danny_dust_m_
    use danny_m_

    implicit none

    character(STRING_BUF_) :: input_file = "input/input"
    character(STRING_BUF_) :: output_file = "output.txt"

    number_density_ = 2.0d5
    temperature_    = 10.0d0
    ionization_rate_ = 1.3d-17

    call readInputFile(input_file)

    print *, species_file_
    print *, binding_energy_file_
    print *, reaction_file_
    print *, abundances_file_
    print *, element_file_
    print *, activation_energy_file_
    print *, enthalpy_file_
    print *, number_of_gaseous_species_
    print *, number_of_species_
    print *, number_of_reactions_
    print *, number_of_elements_
    print *, number_of_gaseous_reactions_
    print *, is_three_phase_reaction_

    call allocateArray

    call readElementFile

    call readSpeciesFile

    call readAbundancesFile

    call readReactionFile

    call readBindingEnergyFile

    call readSurfaceActivationEnergyFile

    call readEnthalpyOfFormationFile

    call setDustProperty

    call setBranchingRatio

    call calculateRateCoefficient

    ! call checkReadReaction

    ! call allocateLsodesWorkArray

    ! call checkReactionRateCoefficient


    call timeIntegration2(output_file)

    call checkChemicalState
    
end program main