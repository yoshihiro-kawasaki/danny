
module danny_m_

    use physical_constants
    use danny_variables_m_

    implicit none

contains

! ######################################################################################################
!  Calculate reaction rate coefficient
! ######################################################################################################

subroutine setBranchingRatio 

    use danny_variables_m_
    implicit none

    integer i, j, k, nbranch, idx_r1, idx_r2, idx_p1, idx_p2, idx_p3, idx_p1_2, idx_p2_2, idx_p3_2, idx_p, n_atom
    double precision reaction_enthalpy, exothermicity, prob_reac_bare, prob_reac_ice, pcd_bare, pcd_ice, eps_m
    integer idx_p_array(3)

    do i = 1, number_of_reactions_

        branching_ratio_(i) = 1.0d0

        ! ダスト表面反応
        if ((type_of_reaction_(i) == TYPE_ID_GRAIN_SURFACE_REACTION_) .or. &
            (type_of_reaction_(i) == TYPE_ID_GRAIN_MANTLE_REACTION_)) then 

            nbranch = 0
            idx_r1 = index_of_reaction_species_(1, i)
            idx_r2 = index_of_reaction_species_(2, i)
            idx_p1 = index_of_reaction_species_(4, i)
            idx_p2 = index_of_reaction_species_(5, i)
            idx_p3 = index_of_reaction_species_(6, i)

            !分岐する反応を数える(反応物が全てsurface(s) o rmantle(m)だけをカウント)
            do j = number_of_gaseous_reactions_, number_of_reactions_
                if (((idx_r1 == index_of_reaction_species_(1, j)) .and. (idx_r2 == index_of_reaction_species_(2, j))) .or. &
                    ((idx_r1 == index_of_reaction_species_(2, j)) .and. (idx_r2 == index_of_reaction_species_(1, j)))) then

                    idx_p1_2 = index_of_reaction_species_(4, j)
                    idx_p2_2 = index_of_reaction_species_(5, j)
                    idx_p3_2 = index_of_reaction_species_(6, j)

                    if ((idx_p2_2 == NOT_FOUND_SPECIES_) .and. (idx_p3_2 == NOT_FOUND_SPECIES_)) then

                        if ((species_name_(idx_p1_2)(1:1) == 's') .or. (species_name_(idx_p1_2)(1:1) == 'm')) then 
                            nbranch = nbranch + 1
                        end if

                    else if (idx_p3_2 == NOT_FOUND_SPECIES_) then 

                        if (((species_name_(idx_p1_2)(1:1) == 's') .or. (species_name_(idx_p1_2)(1:1) == 'm')) .and. &
                            ((species_name_(idx_p2_2)(1:1) == 's') .or. (species_name_(idx_p2_2)(1:1) == 'm'))) then 
                                nbranch = nbranch + 1
                        end if

                    else

                        if (((species_name_(idx_p1_2)(1:1) == 's') .or. (species_name_(idx_p1_2)(1:1) == 'm')) .and. &
                            ((species_name_(idx_p2_2)(1:1) == 's') .or. (species_name_(idx_p2_2)(1:1) == 'm')) .and. &
                            ((species_name_(idx_p3_2)(1:1) == 's') .or. (species_name_(idx_p3_2)(1:1) == 'm'))) then 
                                nbranch = nbranch + 1
                        end if

                    endif
                end if
            end do

            ! branchin_ratioの計算(単純に分岐数nbranchで割る)
            if (nbranch == 0) then 
                branching_ratio_(i) = 0.0d0 
            else 
                branching_ratio_(i) = branching_ratio_(i) / dble(nbranch)
            end if

            ! 同じ化学種どうしの反応の場合、2で割る
            if (idx_r1 == idx_r2) then
                branching_ratio_(i) = branching_ratio_(i) * 0.5d0
            end if

            ! 化学脱着の確率の計算
            if (is_chemical_desorption_ .and. (type_of_reaction_(i) == TYPE_ID_GRAIN_SURFACE_REACTION_)) then
                ! sA + sB → sC + sD
                ! sA + sB →  C +  D
                ! 反応熱が生じたら、その熱がダスト表面からの脱着に使われることがある。化学脱着する確率はbranching ratioに含める。

                ! 反応エンタルピーの計算 (単位は kJ/mol) : 反応エンタルピー = (生成物のエンタルピーの和) - (反応物のエンタルピーの和)
                ! 反応エンタルピー < 0 : 発熱反応, 反応エンタルピー > 0 : 吸熱反応
                reaction_enthalpy = enthalpy_of_formation_(idx_p1) - (enthalpy_of_formation_(idx_r1) &
                    + enthalpy_of_formation_(idx_r2))
                if (idx_p2 .ne. NOT_FOUND_SPECIES_) reaction_enthalpy = reaction_enthalpy + enthalpy_of_formation_(idx_p2)
                if (idx_p3 .ne. NOT_FOUND_SPECIES_) reaction_enthalpy = reaction_enthalpy + enthalpy_of_formation_(idx_p3)
                exothermicity = - reaction_enthalpy ! 発熱(吸熱)量 (発熱を正とするためマイナスかけてる)

                ! 単位変換 kJからK (ボルツマン定数で割る)
                exothermicity = exothermicity*1.0d3 / (BOLTZMANN_CONSTANT*1.0d-7) ! 1.0d-7はボルツマン定数をcgsからSI単位系にするときに出る数
                ! 単位変換 molから１回の反応
                exothermicity = exothermicity / AVOGADRO_CONSTANT

                if ((enthalpy_of_formation_(idx_r1) <= -9999.0d0) .or. (enthalpy_of_formation_(idx_r2) <= -9999.0d0) .or.&
                    (enthalpy_of_formation_(idx_p1) <= -9999.0d0)) then
                    exothermicity = -1.0
                end if
                if ((idx_p2 .ne. NOT_FOUND_SPECIES_) .and. (enthalpy_of_formation_(idx_p2) <= -9999.0d0)) then
                    exothermicity = -1.0
                end if
                if ((idx_p3 .ne. NOT_FOUND_SPECIES_) .and. (enthalpy_of_formation_(idx_p3) <= -9999.0d0)) then
                    exothermicity = -1.0
                end if

                if (exothermicity < 0.0d0) then

                    ! 吸熱反応の場合は脱着しない
                    if ((species_name_(idx_p1)(1:1) == 's') .and. &
                       ((idx_p2 .ne. NOT_FOUND_SPECIES_) .and. (species_name_(idx_p2)(1:1) == 's')) .and. &
                       ((idx_p3 .ne. NOT_FOUND_SPECIES_) .and. (species_name_(idx_p3)(1:1) == 's'))) then 
                        prob_reac_bare = 1.0d0
                        prob_reac_ice  = 1.0d0
                    else if ((species_name_(idx_p1)(1:1) == 's') .and. &
                        ((idx_p2 .ne. NOT_FOUND_SPECIES_) .and. (species_name_(idx_p2)(1:1) == 's')) .and. &
                        (idx_p3 == NOT_FOUND_SPECIES_)) then 
                        prob_reac_bare = 1.0d0 
                        prob_reac_ice = 1.0d0 
                    else if ((species_name_(idx_p1)(1:1) == 's') .and. &
                        (idx_p2 == NOT_FOUND_SPECIES_) .and. &
                        (idx_p3 == NOT_FOUND_SPECIES_)) then 
                        prob_reac_bare = 1.0d0 
                        prob_reac_ice  = 1.0d0 
                    else
                        prob_reac_bare = 0.0d0
                        prob_reac_ice  = 0.0d0
                    end if

                else

                    ! 各生成物の脱着率 Minissale et al. 2016
                    idx_p_array(1) = idx_p1
                    idx_p_array(2) = idx_p2
                    idx_p_array(3) = idx_p3
                    n_atom = 0
                    do j = 1, 3
                        idx_p = idx_p_array(j)
                        if (idx_p == NOT_FOUND_SPECIES_) cycle
                        do k = 1, number_of_elements_
                            n_atom = n_atom + number_of_constituent_element_of_species_(k, idx_p)
                        end do
                    end do
                    n_atom = 3 * n_atom

                    prob_reac_bare = 1.0d0
                    prob_reac_ice  = 1.0d0
                    do j = 1, 3 
                        idx_p = idx_p_array(j)
                        if (idx_p == NOT_FOUND_SPECIES_) cycle
                        eps_m = (120.0d0*ATOMIC_MASS_UNIT - species_mass_(idx_p))**2.0d0 / &
                            (120.0d0*ATOMIC_MASS_UNIT + species_mass_(idx_p))**2.0d0

                        pcd_bare = exp(-binding_energy_(idx_p) / (eps_m*exothermicity/dble(n_atom))) ! 表面がbare(silicate)のとき化学脱着する確率
                        pcd_ice  = pcd_bare * 0.1d0 ! 表面がH2Oiceのとき化学脱着する確率

                        if (species_name_(idx_p)(1:1) == 's') then
                            ! 今考えている化学種(idx_p)が化学脱着しない場合の確率は 1 - Pcd_(bare, ice)
                            prob_reac_bare = prob_reac_bare * (1.0d0 - pcd_bare)
                            prob_reac_ice  = prob_reac_ice  * (1.0d0 - pcd_ice)
                        else
                            ! 今考えている化学種(idx_p)が化学脱着する場合の確率
                            prob_reac_bare = prob_reac_bare * pcd_bare
                            prob_reac_ice  = prob_reac_ice * pcd_ice
                        end if

                    end do

                    ! special case (Cazaux et al. 2016)
                    ! sO + sH 
                    if (((species_name_(idx_r1) == 'sO') .and. (species_name_(idx_r2) == 'sH')) .or. &
                        ((species_name_(idx_r1) == 'sH') .and. (species_name_(idx_r2) == 'sO')) ) then 
                        pcd_ice = 0.25d0 
                        if (species_name_(idx_p1)(1:1) == 's') then 
                            prob_reac_ice = 1.0 - pcd_ice
                        else 
                            prob_reac_ice = pcd_ice
                        end if
                    end if

                    ! sOH + sH
                    if (((species_name_(idx_r1) == 'sOH') .and. (species_name_(idx_r2) == 'sH')) .or. &
                        ((species_name_(idx_r1) == 'sH') .and. (species_name_(idx_r2) == 'sOH')) ) then 
                        pcd_ice = 0.30d0 
                        if (species_name_(idx_p1)(1:1) == 's') then 
                            prob_reac_ice = 1.0 - pcd_ice
                        else 
                            prob_reac_ice = pcd_ice
                        end if
                    end if

                    ! sN + sN
                    if ((species_name_(idx_r1) == 'sN') .and. (species_name_(idx_r2) == 'sN')) then 
                        pcd_ice = 0.50d0 
                        if (species_name_(idx_p1)(1:1) == 's') then 
                            prob_reac_ice = 1.0 - pcd_ice
                        else 
                            prob_reac_ice = pcd_ice
                        end if
                    end if
                    
                    ! sH + sH
                    if ((species_name_(idx_r1) == 'sH') .and. (species_name_(idx_r2) == 'sH')) then 
                        pcd_ice = 1.0d0
                        pcd_bare = 1.0d0
                        if (species_name_(idx_p1)(1:1) == 's') then 
                            prob_reac_ice = 1.0d0 - pcd_ice
                            prob_reac_bare = 1.0d0 - pcd_bare
                        else 
                            prob_reac_ice = pcd_ice
                            prob_reac_bare = pcd_bare
                        end if
                    end if

                end if

                value_for_reaction_(3, i) = branching_ratio_(i) * prob_reac_bare
                value_for_reaction_(4, i) = branching_ratio_(i) * prob_reac_ice

            else 

                value_for_reaction_(3, i) = branching_ratio_(i)
                value_for_reaction_(4, i) = branching_ratio_(i)
            
            end if ! ! is_chemical_desorption

        end if ! end ダスト表面反応

    end do ! do i = 1, number_of_reactions

end subroutine setBranchingRatio


subroutine calculateRateCoefficient

    use danny_variables_m_
    implicit none 

    integer i 

    ! 1. Dissociation or ionization of species due to direct collision with cosmic-ray particles.
    do i = type_id_start_(TYPE_ID_1_), type_id_end_(TYPE_ID_1_)
        reaction_rate_coefficient_(i) = value_for_reaction_(1, i) * ionization_rate_
    end do

    ! ! 2. Dissociation or ionization of species due to UV photons emitted following H2 excitation.
    ! do i = type_id_start_(TYPE_ID_2_), type_id_end_(TYPE_ID_2_)
    !     ! (1/(1-omega) = 2, omega = 0.5)
    !     reaction_rate_coefficient_(i) = value_for_reaction_(1, i) * ionization_rate_ * 2.0d0 / number_density_
    ! end do

    ! 3. Dissociation or ionization of neutral species by UV photons with a standard interstellar UV field.
    do i = type_id_start_(TYPE_ID_3_), type_id_end_(TYPE_ID_3_)
        reaction_rate_coefficient_(i) = value_for_reaction_(1, i) * exp(-value_for_reaction_(3, i)*VISUAL_EXTINCTION_)
    end do

    ! 4-8. Bimolecular reactions includes all chemical reactions between two species.
    call calculateGasPhaseReactionRateCoefficient

    do i = type_id_start_(TYPE_ID_10_), type_id_end_(TYPE_ID_10_)
    end do

    write(*,*) dust_abundances_

    call calculateDustChargedGasReactionRateCoefficient

    call calculateThermalDesorptionRateCoefficient

    call calculateCosmicRayDesorptionRateCoefficient

    ! do i = type_id_start_(TYPE_ID_PHOTO_DISS_CR_GRAINS_), type_id_end_(TYPE_ID_PHOTO_DISS_CR_GRAINS_)
    !     ! (1/(1-omega) = 2, omega = 0.5)
    !     reaction_rate_coefficient_(i) = value_for_reaction_(1, i) * ionization_rate_ * 2.0d0 / number_density_
    ! end do

    do i = type_id_start_(TYPE_ID_PHOTO_DISS_UV_GRAINS_), type_id_end_(TYPE_ID_PHOTO_DISS_UV_GRAINS_)
        reaction_rate_coefficient_(i) = value_for_reaction_(1, i) * exp(-value_for_reaction_(3, i)*VISUAL_EXTINCTION_)
    end do

end subroutine calculateRateCoefficient


subroutine calculateRateCoefficientDependAbundances(abundances)

    use danny_variables_m_
    implicit none

    double precision, intent(in) :: abundances(number_of_species_)

    integer i 
    double precision x_sH2O, dsite

    total_abundance_of_dust_surface_species_ = 0.0d0 
    total_abundance_of_dust_mantle_species_  = 0.0d0 
    do i = number_of_gaseous_species_, number_of_species_
        if (species_name_(i)(1:1) == 's') then 
            total_abundance_of_dust_surface_species_ = total_abundance_of_dust_surface_species_ + abundances(i)
        end if
        if (species_name_(i)(1:1) == 'm') then 
            total_abundance_of_dust_mantle_species_ = total_abundance_of_dust_mantle_species_ + abundances(i)
        end if
    end do
    dsite = 1.0 / (binding_sites_per_a_dust_grain_ * dust_abundances_)
    number_of_surface_layers_ = total_abundance_of_dust_surface_species_ * dsite 
    number_of_mantle_layers_ = total_abundance_of_dust_mantle_species_ * dsite
    number_of_total_layers_ = number_of_surface_layers_ + number_of_mantle_layers_
    x_sH2O = abundances(index_sH2O_)

    ! The fraction of the adsorption site covered by H2O ice and bare sites
    coverage_of_H2O_on_dust_surface_ = x_sH2O  / (dust_abundances_ * binding_sites_per_a_dust_grain_)
    coverage_of_H2O_on_dust_surface_ = min(coverage_of_H2O_on_dust_surface_, 1.0d0)
    coverage_of_bare_on_dust_surface_ = 1.0d0 - coverage_of_H2O_on_dust_surface_
    ! coverage_of_H2O_on_dust_surface_ = min(number_of_surface_layers_, 1.0d0)
    ! coverage_of_bare_on_dust_surface_ = 1.0d0 - coverage_of_H2O_on_dust_surface_

    ! 2. Dissociation or ionization of species due to UV photons emitted following H2 excitation.
    do i = type_id_start_(TYPE_ID_2_), type_id_end_(TYPE_ID_2_)
        ! (1/(1-omega) = 2, omega = 0.5)
        reaction_rate_coefficient_(i) = value_for_reaction_(1, i) * ionization_rate_ * 2.0d0 * abundances(index_H_)
    end do

    do i = type_id_start_(TYPE_ID_PHOTO_DISS_CR_GRAINS_), type_id_end_(TYPE_ID_PHOTO_DISS_CR_GRAINS_)
        ! (1/(1-omega) = 2, omega = 0.5)
        reaction_rate_coefficient_(i) = value_for_reaction_(1, i) * ionization_rate_ * 2.0d0 * abundances(index_H_)
    end do

    call calculateAccretionGasOnGrainsRateCoefficient

    call calculatePhotoDesorptionByUV

    call calculatePhotoDesorptionByUVCR

    call calculateDustSurfaceReactionRateCoefficient

    if (is_three_phase_reaction_) then 

        call calculateDustMantleReactionRateCoefficient

    end if

end subroutine calculateRateCoefficientDependAbundances


subroutine checkReactionRateCoefficient

    use danny_variables_m_
    implicit none 

    character(STRING_BUF_) :: rate_check_file = "check_reaction_rate.txt"
    integer :: file = 14
    integer i, idx1, idx2, idx3, idx4, idx5, idx6, idx7, idx8

    call calculateRateCoefficient

    open(file, file = rate_check_file, status='replace')

    do i = 1, number_of_reactions_

        idx1 = index_of_reaction_species_(1, i)
        idx2 = index_of_reaction_species_(2, i)
        idx3 = index_of_reaction_species_(3, i)
        idx4 = index_of_reaction_species_(4, i)
        idx5 = index_of_reaction_species_(5, i)
        idx6 = index_of_reaction_species_(6, i)
        idx7 = index_of_reaction_species_(7, i)
        idx8 = index_of_reaction_species_(8, i)
        write(file, *) reaction_species_name_(1:8, i), type_of_reaction_(i), reaction_rate_coefficient_(i)

    end do

    close(file)

end subroutine


subroutine calculateGasPhaseReactionRateCoefficient

    use danny_variables_m_
    implicit none

    integer i, j, k, n
    double precision T, T_300, Tinv, sqrtT, inv300
    double precision distmin(10), distmax(10)
    integer indice(10)

    T = temperature_
    T_300 = temperature_ / 300.0d0
    Tinv = 1.0 / temperature_
    sqrtT = sqrt(temperature_)
    inv300 = 1.0d0 / 300.0d0

    j = 1
    k = 1
    indice(:) = 0
    distmin(:) = 9999.0d0
    distmax(:) = 9999.0d0 

    do i = type_id_start_(TYPE_ID_4_), type_id_end_(TYPE_ID_8_)

        ! Modified Arrhenius
        if (formula_id_of_gas_reaction_(i) == 3) then

            if (T < lower_temperature_limit_(i)) then

                reaction_rate_coefficient_(i) = value_for_reaction_(1, i) * &
                    (lower_temperature_limit_(i)*inv300)**value_for_reaction_(2, i) * &
                    exp(-value_for_reaction_(3, i)/lower_temperature_limit_(i))

            else if (T > upper_temperature_limit_(i)) then

                reaction_rate_coefficient_(i) = value_for_reaction_(1, i) * &
                    (upper_temperature_limit_(i)*inv300)**value_for_reaction_(2, i) * &
                    exp(-value_for_reaction_(3, i)/upper_temperature_limit_(i))

            else

                reaction_rate_coefficient_(i) = value_for_reaction_(1, i) * &
                    T_300**(value_for_reaction_(2, i)) * &
                    exp(-value_for_reaction_(3, i)*Tinv)

            end if

            ! Check for the presence of several rate coefficients present in the network for the same reaction
            if (reaction_id_(i+1) == reaction_id_(i)) then

                indice(j) = i 
                distmin(j) = lower_temperature_limit_(i) - T 
                distmax(j) = T - upper_temperature_limit_(i)
                j = j + 1

            end if

            if ((reaction_id_(i+1) .ne. reaction_id_(i)) .and. (j .ne. 1)) then 

                indice(j) = i 
                distmin(j) = lower_temperature_limit_(i) - T 
                distmax(j) = T - upper_temperature_limit_(i)

                do k = 1, j 
                    n = indice(k)
                    if (T < lower_temperature_limit_(n)) reaction_rate_coefficient_(n) = 0.0d0 
                    if (T > upper_temperature_limit_(n)) reaction_rate_coefficient_(n) = 0.0d0 
                end do

                if (maxval(reaction_rate_coefficient_(indice(1:j))) < 1.0d-99) then

                    if (minval(abs(distmin)) < minval(abs(distmax))) then
                        n = indice(minloc(abs(distmin), dim=1))
                        reaction_rate_coefficient_(n) = value_for_reaction_(1,n) * &
                            (lower_temperature_limit_(n)*inv300)**value_for_reaction_(2, n) * &
                            exp(-value_for_reaction_(3, n)/lower_temperature_limit_(n))
                    else 
                        n = indice(minloc(abs(distmax), dim=1))
                        reaction_rate_coefficient_(n) = value_for_reaction_(1,n) * &
                            (upper_temperature_limit_(n)*inv300)**value_for_reaction_(2, n) * &
                            exp(-value_for_reaction_(3, n)/upper_temperature_limit_(n))
                    end if

                end if

                j = 1
                indice(:) = 0
                distmin(:) = 9999.0d0 
                distmax(:) = 9999.0d0 

            end if

        end if ! Modified Arrhenius

        ! ionpol1
        if (formula_id_of_gas_reaction_(i) == 4) then

            if (T < lower_temperature_limit_(i)) then

                reaction_rate_coefficient_(i) = value_for_reaction_(1, i) * value_for_reaction_(2, i) * &
                    (0.62d0 + 0.4767d0*value_for_reaction_(3, i)*sqrt(300.0d0/lower_temperature_limit_(i)))

            else if (T > upper_temperature_limit_(i)) then

                reaction_rate_coefficient_(i) = value_for_reaction_(1, i) * value_for_reaction_(2, i) * &
                    (0.62d0 + 0.4767d0*value_for_reaction_(3, i)*sqrt(300.0d0/upper_temperature_limit_(i)))

            else 

                reaction_rate_coefficient_(i) = value_for_reaction_(1, i) * value_for_reaction_(2, i) * &
                    (0.62d0 + 0.4767d0*value_for_reaction_(3, i)*sqrt(300.0d0*Tinv))

            end if

            ! Check for the presence of several rate coefficients present in the network for the same reaction
            if (reaction_id_(i+1) == reaction_id_(i)) then
                indice(j) = i 
                distmin(j) = lower_temperature_limit_(i) - T 
                distmax(j) = T - upper_temperature_limit_(i)
                j = j + 1
            end if

            if ((reaction_id_(i+1) .ne. reaction_id_(i)) .and. (j .ne. 1)) then 

                indice(j) = i 
                distmin(j) = lower_temperature_limit_(i) - T 
                distmax(j) = T - upper_temperature_limit_(i)

                do k = 1, j 
                    n = indice(k)
                    if (T < lower_temperature_limit_(n)) reaction_rate_coefficient_(n) = 0.0d0 
                    if (T > upper_temperature_limit_(n)) reaction_rate_coefficient_(n) = 0.0d0 
                end do

                if (maxval(reaction_rate_coefficient_(indice(1:j))) < 1.0d-99) then

                    if (minval(abs(distmin)) < minval(abs(distmax))) then
                        n = indice(minloc(abs(distmin), dim=1))
                        reaction_rate_coefficient_(n) = value_for_reaction_(1,n) * &
                            (lower_temperature_limit_(n)*inv300)**value_for_reaction_(2, n) * &
                            exp(-value_for_reaction_(3, n)/lower_temperature_limit_(n))
                    else 
                        n = indice(minloc(abs(distmax), dim=1))
                        reaction_rate_coefficient_(n) = value_for_reaction_(1,n) * &
                            (upper_temperature_limit_(n)*inv300)**value_for_reaction_(2, n) * &
                            exp(-value_for_reaction_(3, n)/upper_temperature_limit_(n))
                    end if

                end if

                j = 1
                indice(:) = 0
                distmin(:) = 9999.0d0 
                distmax(:) = 9999.0d0 

            end if

        end if ! ionpol1

        ! ionpol2
        if (formula_id_of_gas_reaction_(i) == 5) then

            if (T < lower_temperature_limit_(i)) then 

                reaction_rate_coefficient_(i) = value_for_reaction_(1, i) * value_for_reaction_(2, i) * &
                    (1.0d0 + 0.0967*value_for_reaction_(3, i)*sqrt(300.0d0/lower_temperature_limit_(i)) + &
                    + (value_for_reaction_(3, i))**2.0d0 * 300.0d0/(10.526d0*lower_temperature_limit_(i)))

            else if (T > upper_temperature_limit_(i)) then 

                reaction_rate_coefficient_(i) = value_for_reaction_(1, i) * value_for_reaction_(2, i) * &
                    (1.0d0 + 0.0967*value_for_reaction_(3, i)*sqrt(300.0d0/upper_temperature_limit_(i)) + &
                    + (value_for_reaction_(3, i))**2.0d0 * 300.0d0/(10.526d0*upper_temperature_limit_(i)))

            else

                reaction_rate_coefficient_(i) = value_for_reaction_(1, i) * value_for_reaction_(2, i) * &
                    (1.0d0 + 0.0967*value_for_reaction_(3, i)*sqrt(300.0d0*Tinv) + &
                    + (value_for_reaction_(3, i))**2.0d0 * 300.0d0*Tinv/(10.526d0))

            end if

            ! Check for the presence of several rate coefficients present in the network for the same reaction
            if (reaction_id_(i+1) == reaction_id_(i)) then
                indice(j) = i 
                distmin(j) = lower_temperature_limit_(i) - T 
                distmax(j) = T - upper_temperature_limit_(i)
                j = j + 1
            end if

            if ((reaction_id_(i+1) .ne. reaction_id_(i)) .and. (j .ne. 1)) then 

                indice(j) = i 
                distmin(j) = lower_temperature_limit_(i) - T 
                distmax(j) = T - upper_temperature_limit_(i)

                do k = 1, j 
                    n = indice(k)
                    if (T < lower_temperature_limit_(n)) reaction_rate_coefficient_(n) = 0.0d0 
                    if (T > upper_temperature_limit_(n)) reaction_rate_coefficient_(n) = 0.0d0 
                end do

                if (maxval(reaction_rate_coefficient_(indice(1:j))) < 1.0d-99) then

                    if (minval(abs(distmin)) < minval(abs(distmax))) then
                        n = indice(minloc(abs(distmin), dim=1))
                        reaction_rate_coefficient_(n) = value_for_reaction_(1,n) * &
                            (lower_temperature_limit_(n)*inv300)**value_for_reaction_(2, n) * &
                            exp(-value_for_reaction_(3, n)/lower_temperature_limit_(n))
                    else 
                        n = indice(minloc(abs(distmax), dim=1))
                        reaction_rate_coefficient_(n) = value_for_reaction_(1,n) * &
                            (upper_temperature_limit_(n)*inv300)**value_for_reaction_(2, n) * &
                            exp(-value_for_reaction_(3, n)/upper_temperature_limit_(n))
                    end if

                end if

                j = 1
                indice(:) = 0
                distmin(:) = 9999.0d0 
                distmax(:) = 9999.0d0 

            end if 

        end if ! ionpol2

    end do

end subroutine calculateGasPhaseReactionRateCoefficient


subroutine calculateDustChargedGasReactionRateCoefficient

    use danny_variables_m_
    implicit none 
    integer i, ts, te, idx1, idx2
    double precision nu, tau_coef, tau, theta_nu, vth_coef, vth, si 

    vth_coef = sqrt(8.0d0*BOLTZMANN_CONSTANT*temperature_/(PI))
    tau_coef = dust_grain_radius_*BOLTZMANN_CONSTANT*temperature_ / (CHARGE_UNIT**2.0d0)

    ts = type_id_start_(TYPE_ID_GRAIN_CHARGED_SPECIS_)
    te = type_id_end_(TYPE_ID_GRAIN_CHARGED_SPECIS_)

    do i = ts, te 

        idx1 = index_of_reaction_species_(1, i)
        idx2 = index_of_reaction_species_(2, i)

        if (idx1 == index_electron_) then 
            si = 0.60d0
        else
            si = 1.0d0
        end if 

        vth = vth_coef / sqrt(species_mass_(idx1))
        nu = species_charge_(idx2) / species_charge_(idx1)
        tau = tau_coef / (species_charge_(idx1)**2.0d0)
        if (nu <= 0.0d0) then
            theta_nu = 0.0d0
        else
            theta_nu = nu / (1.0 + nu**(-0.50d0))
        end if

        if (nu == 0.0d0) then
            reaction_rate_coefficient_(i) = si*vth*dust_grain_cross_section_ * (1.0d0 + sqrt(PI/(2.0d0*tau)))
        else if (nu < 0.0d0) then
            reaction_rate_coefficient_(i) = si*vth*dust_grain_cross_section_ * &
                (1.0d0 - nu/tau) * (1.0d0 + sqrt(2.0 / (tau - 2.0d0*nu)))
        else 
            reaction_rate_coefficient_(i) = si*vth*dust_grain_cross_section_ * &
                (1.0d0 + (4.0d0*tau + 3.0d0*nu)**(-0.5d0))**2.0d0 * exp(-theta_nu/tau)
        end if

    end do

end subroutine calculateDustChargedGasReactionRateCoefficient


subroutine calculateAccretionGasOnGrainsRateCoefficient

    use danny_variables_m_
    implicit none 

    integer i, ts, te, idx1
    double precision vth_coef, si, T, si_ice, si_bare

    T = temperature_
    vth_coef = sqrt(8.0d0*BOLTZMANN_CONSTANT*temperature_/(PI))

    ts = type_id_start_(TYPE_ID_ACCRETION_GAS_SPECIES_)
    te = type_id_end_(TYPE_ID_ACCRETION_GAS_SPECIES_)

    do i = ts, te 

        idx1 = index_of_reaction_species_(1, i)

        ! sticking probability for H or H2 (Chaabouni et al. 2012)
        if (idx1 == index_H_) then
            ! ice
            si_ice = (1.0d0 + 2.5d0*T/52.0d0) / (1.0d0 + T/52.0d0)**(2.5d0)
            ! silicate
            si_bare = (1.0d0 + 2.5d0*T/25.0d0) / (1.0d0 + T/25.0d0)**(2.5d0)
            si = coverage_of_H2O_on_dust_surface_ * si_ice + coverage_of_bare_on_dust_surface_ * si_bare
        else if (idx1 == index_H2_) then
            ! ice 
            si_ice = (1.0d0 + 2.5d0*T/87.0d0) / (1.0d0 + T/87.0d0)**(2.5d0)
            ! silicate
            si_bare = (1.0d0 + 2.5d0*T/56.0d0) / (1.0d0 + T/56.0d0)**(2.5d0)
            si = coverage_of_H2O_on_dust_surface_ * si_ice + coverage_of_bare_on_dust_surface_ * si_bare
        else
            si = 1.0d0 
        end if

        reaction_rate_coefficient_(i) = si * (vth_coef / sqrt(species_mass_(idx1))) * dust_grain_cross_section_ * &
            dust_number_denisty_

    end do

end subroutine calculateAccretionGasOnGrainsRateCoefficient


subroutine calculateThermalDesorptionRateCoefficient

    use danny_variables_m_
    implicit none 

    integer i, ts, te, idx1
    double precision Tinv

    Tinv = 1.0d0 / temperature_

    ts = type_id_start_(TYPE_ID_THERMAL_DESORPTION_)
    te = type_id_end_(TYPE_ID_THERMAL_DESORPTION_)

    do i = ts, te 

        idx1 = index_of_reaction_species_(1, i)
        reaction_rate_coefficient_(i) = vibrational_frequency_(idx1) * exp(-binding_energy_(idx1)*Tinv) 
        if (reaction_rate_coefficient_(i) < MINIMUM_RATE_COEFFICIENT_) reaction_rate_coefficient_(i) = 0.0d0

    end do

end subroutine calculateThermalDesorptionRateCoefficient


subroutine calculateCosmicRayDesorptionRateCoefficient

    ! Hasegawa & Herbst 1993a
    use danny_variables_m_
    implicit none 
    integer i, ts, te, idx1
    double precision Tinv, f70K 

    Tinv = 1.0d0 / 70.0d0
    f70K = 3.16d-19

    ts = type_id_start_(TYPE_ID_CR_DESORPTION_)
    te = type_id_end_(TYPE_ID_CR_DESORPTION_)

    do i = ts, te

        idx1 = index_of_reaction_species_(1, i)
        reaction_rate_coefficient_(i) = f70K * vibrational_frequency_(idx1) * exp(-binding_energy_(idx1)*Tinv)
        if (reaction_rate_coefficient_(i) < MINIMUM_RATE_COEFFICIENT_) reaction_rate_coefficient_(i) = 0.0d0

    end do

end subroutine calculateCosmicRayDesorptionRateCoefficient


subroutine calculatePhotoDesorptionByUV

    ! Ruaud et al. 2016
    use danny_variables_m_
    implicit none 
    integer i, ts, te 
    double precision S_UV, F_UV, Ypd

    S_UV = 1.0d0 
    F_UV = 1.0d8 
    Ypd = 1.0d-4 

    ts = type_id_start_(TYPE_ID_PHOTO_DESORPTION_UV_)
    te = type_id_end_(TYPE_ID_PHOTO_DESORPTION_UV_)

    do i = ts, te 

        reaction_rate_coefficient_(i) = F_UV * S_UV * exp(-2.0d0*VISUAL_EXTINCTION_) * Ypd / BINDING_SITES_DENSITY_
        if (number_of_total_layers_ >=  number_of_surface_active_layer_) then
            reaction_rate_coefficient_(i) = reaction_rate_coefficient_(i) * number_of_surface_active_layer_ &
                / number_of_total_layers_
        end if

    end do

end subroutine calculatePhotoDesorptionByUV


subroutine calculatePhotoDesorptionByUVCR

    ! Ruaud et al. 2016
    use danny_variables_m_
    implicit none 
    integer i, ts, te 
    double precision S_UV_CR, F_UV_CR, Ypd

    S_UV_CR = 1.0d0 
    F_UV_CR = 1.0d4 
    Ypd = 1.0d-4 

    ts = type_id_start_(TYPE_ID_PHOTO_DESORPTION_UV_CR_)
    te = type_id_end_(TYPE_ID_PHOTO_DESORPTION_UV_CR_)

    do i = ts, te 

        reaction_rate_coefficient_(i) = F_UV_CR * S_UV_CR * Ypd / BINDING_SITES_DENSITY_
        if (number_of_total_layers_ >=  number_of_surface_active_layer_) then
            reaction_rate_coefficient_(i) = reaction_rate_coefficient_(i) * number_of_surface_active_layer_ &
                / number_of_total_layers_ 
        end if

    end do

end subroutine calculatePhotoDesorptionByUVCR


subroutine calculateDustSurfaceReactionRateCoefficient

    ! Hasegawa et al. 1992, Furuya et al. 2015, Ruaud et al. 2016
    use danny_variables_m_
    implicit none 

    integer i, ts, te, idx1, idx2
    double precision khop1, khop2, Tinv, kappa, dsite, numax, kappa0, ktunn1, ktunn2, act_Tinv
    double precision pcd_ice, pcd_bare, Pcd

    ts = type_id_start_(TYPE_ID_GRAIN_SURFACE_REACTION_)
    te = type_id_end_(TYPE_ID_GRAIN_SURFACE_REACTION_)

    Tinv = 1.0d0 / temperature_
    dsite = 1.0d0 / (binding_sites_per_a_dust_grain_ * dust_number_denisty_)

    do i = ts, te 

        idx1 = index_of_reaction_species_(1, i)
        idx2 = index_of_reaction_species_(2, i)

        khop1 = vibrational_frequency_(idx1) * exp(-diffusion_barrier_(idx1)*Tinv)
        khop2 = vibrational_frequency_(idx2) * exp(-diffusion_barrier_(idx2)*Tinv)

        if ((idx1 == index_sH2_) .and. (idx2 == index_sH2_)) then
            reaction_rate_coefficient_(i) = khop1 * dsite
            cycle
        end if

        ! tunnelling for H or H2 (Hasegawaet al. 1992)
        ! if (idx1 == index_sH_) then
        !     ktunn1 = vibrational_frequency_(idx1) * exp(-2.0d0*WIDTH_OF_BARRIER_/DIRAC_CONSTANT* &
        !         sqrt(2.0d0*species_mass_(idx1)*BOLTZMANN_CONSTANT*diffusion_barrier_(idx1)))
        !     if (ktunn1 > khop1) then
        !         khop1 = ktunn1
        !     end if
        ! end if

        ! if (idx2 == index_sH_) then
        !     ktunn2 = vibrational_frequency_(idx2) * exp(-2.0d0*WIDTH_OF_BARRIER_/DIRAC_CONSTANT* &
        !         sqrt(2.0d0*species_mass_(idx2)*BOLTZMANN_CONSTANT*diffusion_barrier_(idx2)))
        !     if (ktunn2 > khop2) then
        !         khop2 = ktunn2
        !     end if
        ! end if

        ! if (idx1 == index_sH2_) then 
        !     ktunn1 = vibrational_frequency_(idx1) * exp(-2.0d0*WIDTH_OF_BARRIER_/DIRAC_CONSTANT* &
        !         sqrt(2.0d0*species_mass_(idx1)*BOLTZMANN_CONSTANT*diffusion_barrier_(idx1)))
        !     if (ktunn1 > khop1) then
        !         khop1 = ktunn1
        !     end if
        ! end if

        ! if (idx2 == index_sH2_) then 
        !     ktunn2 = vibrational_frequency_(idx2) * exp(-2.0d0*WIDTH_OF_BARRIER_/DIRAC_CONSTANT* &
        !         sqrt(2.0d0*species_mass_(idx2)*BOLTZMANN_CONSTANT*diffusion_barrier_(idx2)))
        !     if (ktunn2 > khop2) then
        !         khop2 = ktunn2
        !     end if
        ! end if

        kappa = 1.0d0
        ! activation barrier (see danny_read_m readSurfaceActivationEnergyFile)
        ! value_for_reaction(1, i) : activation energy [T]
        ! value_for_reaction(2, i) : tunnelling effect probabilty
        if (value_for_reaction_(1, i) > 0.0d0) then 
            act_Tinv = value_for_reaction_(1, i)*Tinv
            if (act_Tinv > value_for_reaction_(2, i)) then
                act_Tinv = value_for_reaction_(2, i)
            end if
            kappa0 = exp(-act_Tinv)
            numax  = max(vibrational_frequency_(idx1), vibrational_frequency_(idx2))
            kappa = numax*kappa0 / (numax*kappa0 + khop1 + khop2)
        end if

        pcd_bare = value_for_reaction_(3, i)
        pcd_ice = value_for_reaction_(4, i)
        pcd = coverage_of_bare_on_dust_surface_*pcd_bare + coverage_of_H2O_on_dust_surface_*pcd_ice
        ! if ((idx1 == index_sH2_) .and. (idx2 == index_sH2_)) pcd = 1.0d0

        reaction_rate_coefficient_(i) = kappa * pcd * (khop1 + khop2) * dsite

    end do

end subroutine calculateDustSurfaceReactionRateCoefficient


subroutine calculateDustMantleReactionRateCoefficient

    use danny_variables_m_
    implicit none 

    integer i, ts, te, idx1, idx2
    double precision khop1, khop2, Tinv, kappa, dsite, numax, kappa0, ktunn1, ktunn2, act_Tinv

    ts = type_id_start_(TYPE_ID_GRAIN_MANTLE_REACTION_)
    te = type_id_end_(TYPE_ID_GRAIN_MANTLE_REACTION_)

    Tinv = 1.0d0 / temperature_
    dsite = 1.0d0 / (binding_sites_per_a_dust_grain_*dust_number_denisty_)

    do i = ts, te 

        idx1 = index_of_reaction_species_(1, i)
        idx2 = index_of_reaction_species_(2, i)

        khop1 = vibrational_frequency_(idx1) * exp(-diffusion_barrier_(idx1)*Tinv)
        khop2 = vibrational_frequency_(idx2) * exp(-diffusion_barrier_(idx2)*Tinv)

        ! tunnelling for H or H2 (Hasegawaet al. 1992)
        ! if (idx1 == index_mH_) then
        !     ktunn1 = vibrational_frequency_(idx1) * exp(-2.0d0*WIDTH_OF_BARRIER_/DIRAC_CONSTANT* &
        !         sqrt(2.0d0*species_mass_(idx1)*BOLTZMANN_CONSTANT*diffusion_barrier_(idx1)))
        !     if (ktunn1 > khop1) then
        !         khop1 = ktunn1
        !     end if
        ! end if

        ! if (idx2 == index_mH_) then
        !     ktunn2 = vibrational_frequency_(idx2) * exp(-2.0d0*WIDTH_OF_BARRIER_/DIRAC_CONSTANT* &
        !         sqrt(2.0d0*species_mass_(idx2)*BOLTZMANN_CONSTANT*diffusion_barrier_(idx2)))
        !     if (ktunn2 > khop2) then
        !         khop2 = ktunn2
        !     end if
        ! end if

        ! if (idx1 == index_mH2_) then 
        !     ktunn1 = vibrational_frequency_(idx1) * exp(-2.0d0*WIDTH_OF_BARRIER_/DIRAC_CONSTANT* &
        !         sqrt(2.0d0*species_mass_(idx1)*BOLTZMANN_CONSTANT*diffusion_barrier_(idx1)))
        !     if (ktunn1 > khop1) then
        !         khop1 = ktunn1
        !     end if
        ! end if

        ! if (idx2 == index_mH2_) then 
        !     ktunn2 = vibrational_frequency_(idx2) * exp(-2.0d0*WIDTH_OF_BARRIER_/DIRAC_CONSTANT* &
        !         sqrt(2.0d0*species_mass_(idx2)*BOLTZMANN_CONSTANT*diffusion_barrier_(idx2)))
        !     if (ktunn2 > khop2) then
        !         khop2 = ktunn2
        !     end if
        ! end if

        kappa = 1.0d0
        ! activation barrier (see danny_read_m readSurfaceActivationEnergyFile)
        ! value_for_reaction(1, i) : activation energy [T]
        ! value_for_reaction(2, i) : tunnelling effect probabilty
        if (value_for_reaction_(1, i) > 0.0d0) then 
            act_Tinv = value_for_reaction_(1, i)*Tinv
            if (act_Tinv > value_for_reaction_(2, i)) then
                act_Tinv = value_for_reaction_(2, i)
            end if
            kappa0 = exp(-act_Tinv)
            numax  = max(vibrational_frequency_(idx1), vibrational_frequency_(idx2))
            kappa = numax*kappa0 / (numax*kappa0 + khop1 + khop2)
        end if

        reaction_rate_coefficient_(i) = branching_ratio_(i) * kappa * (khop1 + khop2) * dsite

        if (number_of_mantle_layers_ > 1.0d0) then
            reaction_rate_coefficient_(i) = reaction_rate_coefficient_(i) / number_of_mantle_layers_
        end if 


        
    end do

end subroutine calculateDustMantleReactionRateCoefficient


subroutine calculateDustSurfaceAndMantleSwappingRateCoefficient(species_abundances)

    ! Ruaud et al. 2016
    use danny_variables_m_
    implicit none 
    double precision, intent(in) :: species_abundances(number_of_species_)

    integer i, ts, te, idx1
    double precision alpha_loss, alpha_gain
    double precision swap_mantle_to_surface, swap_surface_to_mantle, sum_swap_mantle_to_surface
    double precision y, Tinv, r1, r2

    swap_mantle_to_surface = 0.0d0 
    swap_surface_to_mantle = 0.0d0 
    sum_swap_mantle_to_surface = 0.0d0 
    Tinv = 1.0d0 / temperature_

    ! mantle to surface
    ts = type_id_start_(TYPE_ID_MANTLE_TO_SURFACE_)
    te = type_id_end_(TYPE_ID_MANTLE_TO_SURFACE_)

    alpha_loss = total_abundance_of_dust_mantle_species_ / total_abundance_of_dust_surface_species_
    if (alpha_loss >= 1.0d0)  alpha_loss = 1.0d0 
    r1 = - alpha_loss * total_desorption_rate_ / total_abundance_of_dust_mantle_species_

    ! write(*,*) "total_desorption = ", total_desorption_rate_
    ! write(*,*) "total_accretion  = ", total_accretion_rate_

    do i = ts, te  

        idx1 = index_of_reaction_species_(1, i)
        y = species_abundances(idx1)
        if (y > MINIMUM_ABUNDANCES_) then
            if (number_of_mantle_layers_ < 1.0d0) then 
                swap_mantle_to_surface = vibrational_frequency_(idx1) * exp(-binding_energy_(idx1)*Tinv)
            else 
                swap_mantle_to_surface = vibrational_frequency_(idx1) * exp(-binding_energy_(idx1)*Tinv) &
                    / number_of_mantle_layers_
            end if 
            sum_swap_mantle_to_surface = sum_swap_mantle_to_surface + swap_mantle_to_surface * y
        else
            swap_mantle_to_surface = 0.0d0 
        end if

        r2 = swap_mantle_to_surface

        reaction_rate_coefficient_(i) = r1 + r2

    end do

    ! surface to mantle
    ts = type_id_start_(TYPE_ID_SURAFCE_TO_MANTLE_)
    te = type_id_end_(TYPE_ID_SURAFCE_TO_MANTLE_)

    alpha_gain = number_of_surface_layers_ / number_of_surface_active_layer_ !< 0.5 : beta = 2 active layer
    r1 = alpha_gain * total_accretion_rate_ / total_abundance_of_dust_surface_species_

    do i = ts, te 

        idx1 = index_of_reaction_species_(1, i)
        y = species_abundances(idx1)
        if (y > MINIMUM_ABUNDANCES_) then 
            swap_surface_to_mantle = sum_swap_mantle_to_surface / total_abundance_of_dust_surface_species_
        else
            swap_surface_to_mantle = 0.0d0 
        end if 

        r2 = swap_surface_to_mantle

        reaction_rate_coefficient_(i) = r1 + r2 

    end do

end subroutine calculateDustSurfaceAndMantleSwappingRateCoefficient


! ######################################################################################################
!  Integration
! ######################################################################################################


subroutine timeIntegration(output_file)

    use danny_variables_m_
    implicit none 

    character(STRING_BUF_), intent(in) :: output_file
    integer IOPT,ISTATE, ITASK, ITOL, MF
    double precision T, TOUT, T_END, RTOL, ATOL(number_of_species_)
    integer nstep, i, j
    double precision tstep

    double precision x(number_of_species_)

    ITOL   = 2
    RTOL   = 1.0e-5 ! reltol
    ATOL   = 1.0e-15 !abstol
    ITASK  = 1
    ISTATE = 1
    IOPT   = 1

    call allocateLsodesWorkArray

    T     = 0.0d0
    TOUT  = 1.0d0 * SOLAR_YEAR * 1.0d-1
    T_END = 1.0e6* SOLAR_YEAR
    nstep = 100
    tstep = 10.0**(log10(T_END/TOUT)/(nstep - 1.0))

    open(14, file=output_file, status='replace')
    write(14, *) species_name_(:)

    call setInitialCondition(x)

    write(14, *) T, x(:)

    MF = 121

    atol(:) = 1.0e-15

    call setWorkArrays

    do i = 1, nstep

        ISTATE = 1
        call setWorkArrays

        do j = 1, number_of_species_
            atol(j) = max(1.0d-20, 1.0d-15*x(j))
        end do

        call DLSODES(ordinaryDifferentialEquation, number_of_species_, x, T, TOUT, ITOL, RTOL, atol, ITASK, ISTATE, &
            IOPT, rwork_, lrw_, iwork_, liw_, jacobianJth, MF)

        do j = 1, number_of_species_
            if (x(j) < MINIMUM_ABUNDANCES_) x(j) = MINIMUM_ABUNDANCES_
        end do

        write(14, *) T/SOLAR_YEAR, x(:)
        write(*, *) T/SOLAR_YEAR, number_of_surface_layers_, number_of_mantle_layers_

        if (ISTATE .lt. 0) ISTATE = 1

        TOUT = TOUT * tstep

    end do

    close(14)
    species_abundances_(:) = x(:)

end subroutine timeIntegration


subroutine timeIntegration2(output_file)

    use danny_variables_m_
    implicit none

    ! external ordinary_differential_equation, jacobian
    character(STRING_BUF_) output_file
    double precision T, TOUT, T_END, dt
    integer nstep, i
    double precision tstep

    double precision x(number_of_species_)

    call allocateLsodesWorkArray

    T     = 0.0d0
    TOUT  = 1.0d0 * SOLAR_YEAR * 1.0d-1
    T_END = 1.0e6 * SOLAR_YEAR
    nstep = 100
    tstep = 10.0**(log10(T_END/TOUT)/(nstep - 1.0))

    call setInitialCondition(x)

    open(14, file=output_file, status='replace')
    write(14, *) species_name_(:)
    write(14, *) T, x(:)

    call setWorkArrays

    do i = 1, nstep

        dt = TOUT - T

        call timeIntegrationOneStep(dt, x)

        T = T + dt

        write(14, *) TOUT/SOLAR_YEAR, x(:)
        write(*, *) T/SOLAR_YEAR, number_of_surface_layers_, number_of_mantle_layers_
        
        TOUT = TOUT * tstep

    end do

    close(14)
    species_abundances_(:) = x(:)

end subroutine timeIntegration2



subroutine timeIntegrationOneStep(dt, temp_abundance)

    implicit none 
    double precision, intent(in) :: dt
    double precision, intent(inout) :: temp_abundance(number_of_species_)

    integer IOPT,ISTATE, ITASK, ITOL, MF
    double precision t, tout, RTOL, ATOL(number_of_species_)
    integer i

    ITOL   = 2
    RTOL   = 1.0e-5 ! reltol
    ATOL   = 1.0e-15 !abstol
    ITASK  = 1
    ISTATE = 1
    IOPT   = 1

    t      = 0.0d0
    tout   = dt

    MF = 121
    ! MF = 222

    atol(:) = 1.0e-15

    ! call set_work_aarys

    do while (t < tout)

        ISTATE = 1

        do i  = 1, number_of_species_
            atol(i) = max(1.0d-20, 1.0d-15*temp_abundance(i))
        end do

        call setWorkArrays

        call DLSODES(ordinaryDifferentialEquation, number_of_species_, temp_abundance, T, TOUT, ITOL, RTOL, atol, ITASK, ISTATE, &
            IOPT, rwork_, lrw_, iwork_, liw_, jacobianJth, MF)

        ! call DLSODES(ordinaryDifferentialEquation, number_of_species_, temp_abundance, T, TOUT, ITOL, RTOL, atol, ITASK, ISTATE, &
        !     IOPT, rwork_, lrw_, iwork_, liw_, dummy_jac, MF)

        do i = 1, number_of_species_
            if (temp_abundance(i) < MINIMUM_ABUNDANCES_) temp_abundance(i) = MINIMUM_ABUNDANCES_
        end do

        if (ISTATE .ne. 2) then 
            write(*,*) "ISTATE = ", ISTATE
            ! ISTATE = 1
        end if

    end do

end subroutine timeIntegrationOneStep

! subroutine timeIntegrationOneStep(dt, temp_abundance)

!     implicit none 
!     double precision, intent(in) :: dt
!     double precision, intent(inout) :: temp_abundance(number_of_species_)

!     integer IOPT,ISTATE, ITASK, ITOL, MF
!     double precision t, tout, RTOL, ATOL(number_of_species_), x(number_of_species_)
!     double precision t_temp
!     integer i

!     ITOL   = 2
!     RTOL   = 1.0e-5 ! reltol
!     ATOL   = 1.0e-15 !abstol
!     ITASK  = 1
!     ISTATE = 1
!     IOPT   = 1

!     t      = 0.0d0
!     tout   = dt

!     MF = 121
!     ! MF = 222

!     atol(:) = 1.0e-15

!     ! call set_work_aarys

!     do while (t < tout)

!         ISTATE = 1

!         do while (ISTATE .eq. 2)

!             x(:) = temp_abundance(:)

!             do i  = 1, number_of_species_
!                 atol(i) = max(1.0d-20, 1.0d-15*x(i))
!             end do

!             call setWorkArrays

!             call DLSODES(ordinaryDifferentialEquation, number_of_species_, x, T, TOUT, ITOL, RTOL, atol, ITASK, ISTATE, &
!                 IOPT, rwork_, lrw_, iwork_, liw_, jacobianJth, MF)

!             do i = 1, number_of_species_
!                 if (temp_abundance(i) < MINIMUM_ABUNDANCES_) temp_abundance(i) = MINIMUM_ABUNDANCES_
!             end do

!             if (ISTATE .ne. 2) then 
!                 write(*,*) "ISTATE = ", ISTATE
!                 ! ISTATE = 1
!                 RTOL = RTOL*1.1d0
!             end if

!         end do

!         temp_abundance(:) = x(:)

!     end do

! end subroutine timeIntegrationOneStep

subroutine setInitialCondition(x)

    use danny_variables_m_
    implicit none
    
    double precision :: electron_abund = 0.0d0

    double precision, intent(out) :: x(number_of_species_)
    integer i 

    x(:) = species_abundances_(:)

    do i = 1, number_of_species_
        if ((species_charge_(i) .ne. 0.0d0) .and. (species_name_(i) .ne. "e-")) then
            electron_abund = electron_abund + species_charge_(i) * x(i)
        end if
        if ((species_name_(i)(1:1) == 'G') .and. (species_charge_(i) == 0.0d0))  then
            x(i) = dust_abundances_
        end if
    end do

    x(index_electron_) = electron_abund

end subroutine setInitialCondition


subroutine checkChemicalState

    use danny_variables_m_
    implicit none

    integer i 
    double precision total_charge, cation_density, anion_density
    double precision total_dust_number_density_result, mean_dust_charge, error_dust

    total_charge = 0.0d0
    cation_density = 0.0d0
    anion_density = 0.0d0
    total_dust_number_density_result = 0.0d0
    mean_dust_charge = 0.0d0 

    do i = 1, number_of_species_
        total_charge = total_charge + species_charge_(i) * species_abundances_(i)
        if (species_charge_(i) > 0.0d0) then 
            cation_density = cation_density + species_charge_(i) * species_abundances_(i)
        else if (species_charge_(i) < 0.0d0) then 
            anion_density = anion_density + species_charge_(i) * species_abundances_(i)
        end if 

        if (species_name_(i)(1:1) .eq. "G") then 
            total_dust_number_density_result = total_dust_number_density_result + species_abundances_(i)
            mean_dust_charge = mean_dust_charge + species_charge_(i) * species_abundances_(i)
        end if

    end do

    if (total_dust_number_density_result == 0.0d0) then
        mean_dust_charge = 0.0d0 
    else 
        mean_dust_charge = mean_dust_charge / total_dust_number_density_result
    end if 
    total_dust_number_density_result = total_dust_number_density_result * number_density_
    error_dust = abs((total_dust_number_density_result - dust_number_denisty_) / dust_number_denisty_)
    
    write(*,*) "initial total dust number density = ", dust_number_denisty_
    write(*,*) "finish total dust number density  = ", total_dust_number_density_result
    write(*,*) "error dust number density         = ", error_dust
    write(*,*) "mean dust charge                  = ", mean_dust_charge
    write(*,*) "total charge                      = ", total_charge
    write(*,*) "total cation number density       = ", cation_density
    write(*,*) "total antion number density       = ", anion_density
    write(*,*) "electron number density           = ", species_abundances_(index_electron_)

end subroutine checkChemicalState


subroutine allocateLsodesWorkArray

    use danny_variables_m_
    implicit none

    integer i, j
    integer count_nonzeros, max_nonzeros
    double precision x(number_of_species_), PDJ(number_of_species_), dxdt(number_of_species_)
    double precision t
    double precision IAN(number_of_species_), JAN(number_of_species_)

    ! ヤコビアンの各列毎で0ではない成分をカウントし、その最大値を取得
    max_nonzeros = 0
    x(:) = 1.0d-5
    ! call setInitialCondition(x)
    dxdt(:) = 0.0d0
    t = 0.0d0

    call calculateRateCoefficient
    call calculateRateCoefficientDependAbundances(x)
    call ordinaryDifferentialEquation(number_of_species_, 0.0d0, x, dxdt)

    ! stop
    ! if (is_three_phase_reaction_) call calculateDustSurfaceAndMantleSwappingRateCoefficient(x)

    do j = 1, number_of_species_

        call jacobianJth(number_of_species_, t, x, j, IAN, JAN, PDJ)

        count_nonzeros = 0
        do i = 1, number_of_species_
            if (PDJ(i) .eq. 0.0d0) then
                count_nonzeros = count_nonzeros + 1
            end if
        end do

        if (count_nonzeros > max_nonzeros) then
            max_nonzeros = count_nonzeros
        end if

    end do

    ! lsodesで使用するwork arraysの初期化
    liw_ = 30
    ! liw = 31 + 3*max_nonzeros*number_of_species + 21*number_of_species
    lrw_ = 20 + 2*max_nonzeros*number_of_species_ + 15*number_of_species_

    allocate(iwork_(liw_))
    allocate(rwork_(lrw_))

    iwork_(:) = 0
    rwork_(:) = 0.0d0 

end subroutine allocateLsodesWorkArray


subroutine setWorkArrays

    use danny_variables_m_
    implicit none 

    iwork_(:) = 0
    rwork_(:) = 0.0d0 

    rwork_(6) = 3.154d14
    iwork_(6) = 3000

end subroutine setWorkArrays


subroutine ordinaryDifferentialEquation(neq, t, x, dxdt)

    use danny_variables_m_
    implicit none 

    integer neq
    double precision, intent(in) :: t
    double precision, intent(in) :: x(neq)
    double precision, intent(out) :: dxdt(neq)

    integer i, ir1, ir2, ir3, ip1, ip2, ip3, ip4, ip5, type_id
    double precision :: k = 0.0d0
    double precision dxdt_gain(neq), dxdt_loss(neq)
    double precision total_dxdt_gain, total_dxdt_loss

    dxdt(1:neq) = 0.0d0
    dxdt_gain(1:neq) = 0.0d0
    dxdt_loss(1:neq) = 0.0d0

    call calculateRateCoefficientDependAbundances(x)

    do i = 1, number_of_reactions_

        type_id = type_of_reaction_(i)
        if ((type_id == TYPE_ID_SURAFCE_TO_MANTLE_) .or. (type_id == TYPE_ID_MANTLE_TO_SURFACE_)) cycle

        ir1 = index_of_reaction_species_(1, i)
        ir2 = index_of_reaction_species_(2, i)
        ir3 = index_of_reaction_species_(3, i)

        ip1 = index_of_reaction_species_(4, i)
        ip2 = index_of_reaction_species_(5, i)
        ip3 = index_of_reaction_species_(6, i)
        ip4 = index_of_reaction_species_(7, i)
        ip5 = index_of_reaction_species_(8, i)

        if (ir2 == NOT_FOUND_SPECIES_) then
            k = reaction_rate_coefficient_(i) * x(ir1)
        else 
            if (ir3 == NOT_FOUND_SPECIES_) then
                k = reaction_rate_coefficient_(i) * x(ir1) * x(ir2) * number_density_
            else
                k = reaction_rate_coefficient_(i) * x(ir1) * x(ir2) * x(ir3) * number_density_ * number_density_
            end if
        end if

        dxdt(ir1) = dxdt(ir1) - k
        if (ir2 .ne. NOT_FOUND_SPECIES_) dxdt(ir2) = dxdt(ir2) - k
        if (ir3 .ne. NOT_FOUND_SPECIES_) dxdt(ir3) = dxdt(ir3) - k

        dxdt(ip1) = dxdt(ip1) + k
        if (ip2 .ne. NOT_FOUND_SPECIES_) dxdt(ip2) = dxdt(ip2) + k
        if (ip3 .ne. NOT_FOUND_SPECIES_) dxdt(ip3) = dxdt(ip3) + k
        if (ip4 .ne. NOT_FOUND_SPECIES_) dxdt(ip4) = dxdt(ip4) + k
        if (ip5 .ne. NOT_FOUND_SPECIES_) dxdt(ip5) = dxdt(ip5) + k

        dxdt_loss(ir1) = dxdt_loss(ir1) - k 
        if (ir2 .ne. NOT_FOUND_SPECIES_) dxdt_loss(ir2) = dxdt_loss(ir2) - k
        if (ir3 .ne. NOT_FOUND_SPECIES_) dxdt_loss(ir3) = dxdt_loss(ir3) - k 

        dxdt_gain(ip1) = dxdt_gain(ip1) + k
        if (ip2 .ne. NOT_FOUND_SPECIES_) dxdt_gain(ip2) = dxdt_gain(ip2) + k
        if (ip3 .ne. NOT_FOUND_SPECIES_) dxdt_gain(ip3) = dxdt_gain(ip3) + k
        if (ip4 .ne. NOT_FOUND_SPECIES_) dxdt_gain(ip4) = dxdt_gain(ip4) + k
        if (ip5 .ne. NOT_FOUND_SPECIES_) dxdt_gain(ip5) = dxdt_gain(ip5) + k


    end do

    if (is_three_phase_reaction_) then
        total_dxdt_gain = 0.0d0 
        total_dxdt_loss = 0.0d0 
        do i = number_of_gaseous_species_, number_of_species_
            if (species_name_(i)(1:1) == 's') then
                total_dxdt_gain = total_dxdt_gain + dxdt_gain(i)
                total_dxdt_loss = total_dxdt_loss + dxdt_loss(i)
            end if
        end do

        total_desorption_rate_ = total_dxdt_loss
        total_accretion_rate_ = total_dxdt_gain
        call calculateDustSurfaceAndMantleSwappingRateCoefficient(x)

        do i = type_id_start_(TYPE_ID_SURAFCE_TO_MANTLE_), type_id_end_(TYPE_ID_MANTLE_TO_SURFACE_)

            ir1 = index_of_reaction_species_(1, i)
            ip1 = index_of_reaction_species_(4, i)

            k = reaction_rate_coefficient_(i) * x(ir1)
            dxdt(ir1) = dxdt(ir1) - k 
            dxdt(ip1) = dxdt(ip1) + k 

        end do

    end if

end subroutine ordinaryDifferentialEquation

subroutine jacobianJth(neq, t, x, j, IAN, JAN, PDJ)

    use danny_variables_m_
    implicit none

    integer, intent(in) :: neq, j 
    double precision, intent(in) :: t 
    double precision, intent(in) :: IAN(neq), JAN(neq)
    double precision, intent(in) :: x(neq)
    double precision, intent(out) :: PDJ(neq)

    integer i, ir1, ir2, ir3, ip1, ip2, ip3, ip4, ip5, reaction_idx
    double precision :: k = 0.0d0

    do i = 1, number_of_reaction_species_(j)

        reaction_idx = index_of_relevant_reaction_of_species_(i, j)

        ir1 = index_of_reaction_species_(1, reaction_idx)
        ir2 = index_of_reaction_species_(2, reaction_idx)
        ir3 = index_of_reaction_species_(3, reaction_idx)

        ip1 = index_of_reaction_species_(4, reaction_idx)
        ip2 = index_of_reaction_species_(5, reaction_idx)
        ip3 = index_of_reaction_species_(6, reaction_idx)
        ip4 = index_of_reaction_species_(7, reaction_idx)
        ip5 = index_of_reaction_species_(8, reaction_idx)

        if (ir2 == NOT_FOUND_SPECIES_) then 

            if (ir1 == j) then
                k = reaction_rate_coefficient_(reaction_idx)
                PDJ(ir1) = PDJ(ir1) - k 
                PDJ(ip1) = PDJ(ip1) + k
                if (ip2 .ne. NOT_FOUND_SPECIES_) PDJ(ip2) = PDJ(ip2) + k
                if (ip3 .ne. NOT_FOUND_SPECIES_) PDJ(ip3) = PDJ(ip3) + k
                if (ip4 .ne. NOT_FOUND_SPECIES_) PDJ(ip4) = PDJ(ip4) + k
                if (ip5 .ne. NOT_FOUND_SPECIES_) PDJ(ip5) = PDJ(ip5) + k
            end if

        else if (ir3 == NOT_FOUND_SPECIES_) then 

            if (ir1 == j) then
                k = reaction_rate_coefficient_(reaction_idx) * x(ir2) * number_density_
                PDJ(ir1) = PDJ(ir1) - k
                PDJ(ir2) = PDJ(ir2) - k
                PDJ(ip1) = PDJ(ip1) + k
                if (ip2 .ne. NOT_FOUND_SPECIES_) PDJ(ip2) = PDJ(ip2) + k
                if (ip3 .ne. NOT_FOUND_SPECIES_) PDJ(ip3) = PDJ(ip3) + k
                if (ip4 .ne. NOT_FOUND_SPECIES_) PDJ(ip4) = PDJ(ip4) + k
                if (ip5 .ne. NOT_FOUND_SPECIES_) PDJ(ip5) = PDJ(ip5) + k
            end if

            if (ir2 == j) then
                k = reaction_rate_coefficient_(reaction_idx) * x(ir1) * number_density_
                PDJ(ir1) = PDJ(ir1) - k
                PDJ(ir2) = PDJ(ir2) - k
                PDJ(ip1) = PDJ(ip1) + k
                if (ip2 .ne. NOT_FOUND_SPECIES_) PDJ(ip2) = PDJ(ip2) + k
                if (ip3 .ne. NOT_FOUND_SPECIES_) PDJ(ip3) = PDJ(ip3) + k
                if (ip4 .ne. NOT_FOUND_SPECIES_) PDJ(ip4) = PDJ(ip4) + k
                if (ip5 .ne. NOT_FOUND_SPECIES_) PDJ(ip5) = PDJ(ip5) + k
            end if

        else 

            if (ir1 == j) then
                k = reaction_rate_coefficient_(reaction_idx) * x(ir2) * x(ir3) * number_density_**2.0d0
                PDJ(ir1) = PDJ(ir1) - k
                PDJ(ir2) = PDJ(ir2) - k
                PDJ(ir3) = PDJ(ir3) - k
                PDJ(ip1) = PDJ(ip1) + k
                if (ip2 .ne. NOT_FOUND_SPECIES_) PDJ(ip2) = PDJ(ip2) + k
                if (ip3 .ne. NOT_FOUND_SPECIES_) PDJ(ip3) = PDJ(ip3) + k
                if (ip4 .ne. NOT_FOUND_SPECIES_) PDJ(ip4) = PDJ(ip4) + k
                if (ip5 .ne. NOT_FOUND_SPECIES_) PDJ(ip5) = PDJ(ip5) + k
            end if

            if (ir2 == j) then
                k = reaction_rate_coefficient_(reaction_idx) * x(ir1) * x(ir3) * number_density_**2.0d0
                PDJ(ir1) = PDJ(ir1) - k
                PDJ(ir2) = PDJ(ir2) - k
                PDJ(ir3) = PDJ(ir3) - k
                PDJ(ip1) = PDJ(ip1) + k
                if (ip2 .ne. NOT_FOUND_SPECIES_) PDJ(ip2) = PDJ(ip2) + k
                if (ip3 .ne. NOT_FOUND_SPECIES_) PDJ(ip3) = PDJ(ip3) + k
                if (ip4 .ne. NOT_FOUND_SPECIES_) PDJ(ip4) = PDJ(ip4) + k
                if (ip5 .ne. NOT_FOUND_SPECIES_) PDJ(ip5) = PDJ(ip5) + k
            end if

            if (ir3 == j) then
                k = reaction_rate_coefficient_(reaction_idx) * x(ir1) * x(ir2) * number_density_**2.0d0
                PDJ(ir1) = PDJ(ir1) - k
                PDJ(ir2) = PDJ(ir2) - k
                PDJ(ir3) = PDJ(ir3) - k
                PDJ(ip1) = PDJ(ip1) + k
                if (ip2 .ne. NOT_FOUND_SPECIES_) PDJ(ip2) = PDJ(ip2) + k
                if (ip3 .ne. NOT_FOUND_SPECIES_) PDJ(ip3) = PDJ(ip3) + k
                if (ip4 .ne. NOT_FOUND_SPECIES_) PDJ(ip4) = PDJ(ip4) + k
                if (ip5 .ne. NOT_FOUND_SPECIES_) PDJ(ip5) = PDJ(ip5) + k
            end if

        end if

    end do

end subroutine jacobianJth

end module danny_m_