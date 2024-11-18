module danny_variables_m_

    use physical_constants
    implicit none 

    integer          :: number_of_elements_
    integer          :: number_of_species_          ! total species 
    integer          :: number_of_species_plus_one_
    integer          :: number_of_gaseous_species_
    integer          :: number_of_grains_species_           ! dust grrain species (..G-2, G-1, G0, G+1, G+2, ..)
    integer          :: number_of_reactions_
    integer          :: number_of_gaseous_reactions_
    double precision :: number_density_
    double precision :: temperature_ 
    double precision :: ionization_rate_

    integer, parameter     :: STRING_BUF_ = 100
    integer, parameter     :: SPECIES_STRING_BUF_ = 15
    character(STRING_BUF_)  :: element_file_
    character(STRING_BUF_)  :: species_file_
    character(STRING_BUF_)  :: binding_energy_file_
    character(STRING_BUF_)  :: reaction_file_
    character(STRING_BUF_)  :: abundances_file_
    character(STRING_BUF_)  :: activation_energy_file_
    character(STRING_BUF_)  :: enthalpy_file_
    

    character(SPECIES_STRING_BUF_), allocatable   :: species_name_(:)
    integer, allocatable                          :: number_of_constituent_element_of_species_(:, :)
    character(SPECIES_STRING_BUF_), allocatable   :: element_name_(:)
    double precision, allocatable                 :: elemental_mass_(:)
    double precision, allocatable                 :: species_abundances_(:)
    double precision, allocatable                 :: species_mass_(:)
    double precision, allocatable                 :: species_charge_(:)
    double precision, allocatable                 :: hall_beta_(:)

    integer, allocatable          :: type_of_reaction_(:)
    integer, allocatable          :: index_of_reaction_species_(:,:)
    character(SPECIES_STRING_BUF_), allocatable :: reaction_species_name_(:, :)
    integer, allocatable          :: number_of_reaction_species_(:)
    integer, allocatable          :: index_of_relevant_reaction_of_species_(:,:)

    double precision, allocatable :: branching_ratio_(:)
    double precision, allocatable :: reaction_desorption_probability_(:)
    double precision, allocatable :: reaction_desorption_probability_bare_(:)

    integer, allocatable          :: index_dust_grains_(:)
    double precision, allocatable :: value_for_reaction_(:, :)
    double precision, allocatable :: reaction_rate_coefficient_(:)
    double precision, allocatable :: lower_temperature_limit_(:)
    double precision, allocatable :: upper_temperature_limit_(:)
    integer, allocatable          :: formula_id_of_gas_reaction_(:)
    integer, allocatable          :: reaction_id_(:)
    integer, allocatable          :: type_id_start_(:)
    integer, allocatable          :: type_id_end_(:)

    double precision, allocatable :: binding_energy_(:)      ! H2O ice
    double precision, allocatable :: binding_energy_bare_(:)
    double precision, allocatable :: diffusion_barrier_(:)
    double precision, allocatable :: diffusion_barrier_bare_(:)
    double precision, allocatable :: vibrational_frequency_(:)
    double precision, allocatable :: vibrational_frequency_bare_(:)
    double precision, allocatable :: enthalpy_of_formation_(:)


    double precision :: dust_to_gas_mass_ratio_
    double precision :: dust_internal_density_
    double precision :: dust_grain_radius_
    double precision :: dust_grain_mass_ 
    double precision :: dust_grain_cross_section_ 
    double precision :: dust_number_denisty_
    double precision :: dust_abundances_
    double precision :: dust_temperature_

    double precision :: magnetic_field_

    double precision, allocatable :: x_species_(:)
    double precision xdust_

    double precision :: diff_m_H2O_

    double precision :: coverage_of_H2O_on_dust_surface_   ! coverage fraction of H2O ice
    double precision :: coverage_of_bare_on_dust_surface_  ! coverage fraction of bare (silicate)  

    double precision :: relative_tolerance_ = 1.0d-5
    double precision :: absolute_tolerance_ = 1.0d-15
    double precision :: integration_time_ = 1.0d6 * SOLAR_YEAR

    logical :: is_thermal_ionization_ = .false.
    logical :: is_three_phase_reaction_
    logical :: is_chemical_desorption_

    double precision :: total_abundance_of_dust_surface_species_
    double precision :: total_abundance_of_dust_mantle_species_
    double precision :: number_of_surface_layers_
    double precision :: number_of_mantle_layers_
    double precision :: number_of_total_layers_
    double precision :: number_of_surface_active_layer_ = 4.0d0
    double precision :: total_desorption_rate_
    double precision :: total_accretion_rate_

    ! double precision sigmaO
    ! double precision sigmaH
    ! double precision sigmaP
    ! double precision etaO 
    ! double precision etaH 
    ! double precision etaA

    ! double precision, allocatable :: sigmaO_i(:)
    ! double precision, allocatable :: sigmaH_i(:)
    ! double precision, allocatable :: sigmaP_i(:)

    integer :: index_electron_ = 0
    integer :: index_H_ = 0
    integer :: index_H2_ = 0
    integer :: index_sH_ = 0
    integer :: index_sH2_ = 0
    integer :: index_mH_ = 0
    integer :: index_mH2_ = 0
    integer :: index_sH2O_ = 0

    double precision :: binding_sites_per_a_dust_grain_ ! ダスト１つのbinding siteの数 = BINDING_SITES_DENSITY * (ダストの表面積4pi a^2)


    ! integer, parameter :: NOT_FOUND_SPECIES_ = -1
    integer :: NOT_FOUND_SPECIES_ = -1

    ! constants for reactions
    double precision, parameter :: WIDTH_OF_BARRIER_      = 1.0d-8 ! 1 Å = 1d-8 cm
    double precision, parameter :: TRANSFERRED_ENERGY_    = 2.0d-3 * ELECTRON_VOLT
    double precision, parameter :: BINDING_SITES_DENSITY_ = 1.0d15
    double precision, parameter :: POLARIZABILITY_        = 1.0d0
    ! double precision, parameter :: DUST_INTERNAL_DENSITY = 2.0
    double precision, parameter :: GAS_MOLECULAR_WEIGHT_  = 2.34d0
    double precision, parameter :: GAS_MOLECULAR_MASS_    = GAS_MOLECULAR_WEIGHT_ * PROTON_MASS
    double precision, parameter :: POTASSIUM_ABUNDANCE_   = 3.04d-7

    double precision, parameter :: VISUAL_EXTINCTION_ = 1.0d1

    double precision, parameter :: MINIMUM_RATE_COEFFICIENT_ = 1.0d-99
    double precision, parameter :: MINIMUM_ABUNDANCES_       = 1.0d-40

    character(SPECIES_STRING_BUF_), parameter :: NO_SPECIES_ = "*"

    ! reaction_type ID
    integer, parameter :: TYPE_ID_1_  = 1
    integer, parameter :: TYPE_ID_2_  = 2
    integer, parameter :: TYPE_ID_3_  = 3
    integer, parameter :: TYPE_ID_4_  = 4
    integer, parameter :: TYPE_ID_5_  = 5
    integer, parameter :: TYPE_ID_6_  = 6
    integer, parameter :: TYPE_ID_7_  = 7
    integer, parameter :: TYPE_ID_8_  = 8
    integer, parameter :: TYPE_ID_10_ = 10
    integer, parameter :: TYPE_ID_11_ = 11

    integer, parameter :: TYPE_ID_GRAIN_CHARGED_SPECIS_      = 15
    integer, parameter :: TYPE_ID_ACCRETION_GAS_SPECIES_     = 16
    integer, parameter :: TYPE_ID_THERMAL_DESORPTION_        = 17
    integer, parameter :: TYPE_ID_CR_DESORPTION_             = 18
    integer, parameter :: TYPE_ID_PHOTO_DESORPTION_UV_       = 19
    integer, parameter :: TYPE_ID_PHOTO_DESORPTION_UV_CR_    = 20
    integer, parameter :: TYPE_ID_GRAIN_SURFACE_REACTION_    = 21
    integer, parameter :: TYPE_ID_GRAIN_MANTLE_REACTION_     = 22
    integer, parameter :: TYPE_ID_PHOTO_DISS_CR_GRAINS_      = 23
    integer, parameter :: TYPE_ID_PHOTO_DISS_UV_GRAINS_      = 24
    integer, parameter :: TYPE_ID_SURAFCE_TO_MANTLE_         = 25
    integer, parameter :: TYPE_ID_MANTLE_TO_SURFACE_         = 26

    integer, parameter :: MAX_NUMBER_OF_TYPE_ID_             = 26

    ! variables for LSODES
    integer :: lrw_, liw_
    integer, allocatable :: iwork_(:)
    double precision, allocatable :: rwork_(:)

contains

subroutine allocateArray

    implicit none

    number_of_species_plus_one_ = number_of_species_ + 1
    NOT_FOUND_SPECIES_ = number_of_species_plus_one_

    ! allcatae
    allocate(species_name_(number_of_species_))
    species_name_(1:number_of_species_) = ''

    allocate(number_of_constituent_element_of_species_(number_of_elements_, number_of_species_))
    number_of_constituent_element_of_species_(1:number_of_elements_, 1:number_of_species_) = 0

    allocate(species_abundances_(number_of_species_))
    species_abundances_(1:number_of_species_) = MINIMUM_ABUNDANCES_ !< 0.0d0 

    allocate(species_mass_(number_of_species_))
    species_mass_(1:number_of_species_) = 0.0d0 

    allocate(species_charge_(number_of_species_))
    species_charge_(1:number_of_species_) = 0.0d0 

    allocate(element_name_(number_of_elements_))
    element_name_(1:number_of_elements_) = ''

    allocate(elemental_mass_(number_of_elements_))
    elemental_mass_(1:number_of_elements_) = 0.0d0

    allocate(type_of_reaction_(number_of_reactions_))
    type_of_reaction_(1:number_of_reactions_) = 0

    allocate(index_of_reaction_species_(8, number_of_reactions_))
    index_of_reaction_species_(1:8, 1:number_of_reactions_) = 0

    allocate(reaction_species_name_(8, number_of_reactions_))
    reaction_species_name_(:,:) = ""

    allocate(number_of_reaction_species_(number_of_species_))
    number_of_reaction_species_(:) = 0

    allocate(branching_ratio_(number_of_reactions_))
    branching_ratio_(:) = 0.0d0

    allocate(value_for_reaction_(5, number_of_reactions_))
    value_for_reaction_(1:5, 1:number_of_reactions_) = 0.0d0

    allocate(reaction_rate_coefficient_(number_of_reactions_))
    reaction_rate_coefficient_(1:number_of_reactions_) = 0.0d0

    allocate(lower_temperature_limit_(number_of_reactions_))
    lower_temperature_limit_(1:number_of_reactions_) = 0.0d0 
    
    allocate(upper_temperature_limit_(number_of_reactions_))
    upper_temperature_limit_(1:number_of_reactions_) = 0.0d0 

    allocate(formula_id_of_gas_reaction_(number_of_reactions_))
    formula_id_of_gas_reaction_(1:number_of_reactions_) = 0

    allocate(reaction_id_(number_of_reactions_))
    reaction_id_(1:number_of_reactions_) = 0

    allocate(type_id_start_(MAX_NUMBER_OF_TYPE_ID_))
    type_id_start_(1:MAX_NUMBER_OF_TYPE_ID_) = 0

    allocate(type_id_end_(MAX_NUMBER_OF_TYPE_ID_))
    type_id_end_(1:MAX_NUMBER_OF_TYPE_ID_) = 0

    allocate(x_species_(number_of_species_))
    x_species_(1:number_of_species_) = 0.0d0 

    allocate(binding_energy_(number_of_species_))
    binding_energy_(1:number_of_species_) = 0.0d0

    allocate(binding_energy_bare_(number_of_species_))
    binding_energy_bare_(1:number_of_species_) = 0.0d0 

    allocate(diffusion_barrier_(number_of_species_))
    diffusion_barrier_(1:number_of_species_) = 0.0d0 

    allocate(diffusion_barrier_bare_(number_of_species_))
    diffusion_barrier_bare_(1:number_of_species_) = 0.0d0 

    allocate(vibrational_frequency_(number_of_species_))
    vibrational_frequency_(1:number_of_species_) = 0.0d0

    allocate(vibrational_frequency_bare_(number_of_species_))
    vibrational_frequency_bare_(1:number_of_species_) = 0.0d0

    allocate(enthalpy_of_formation_(number_of_species_))
    enthalpy_of_formation_(:) = 1.0d0

    ! allocate(hall_beta(number_of_species))
    ! hall_beta(1:number_of_species) = 0.0d0
    ! allocate(sigmaO_i(number_of_species))
    ! allocate(sigmaH_i(number_of_species))
    ! allocate(sigmaP_i(number_of_species))

end subroutine allocateArray

end module danny_variables_m_