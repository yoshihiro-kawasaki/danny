module danny_read_m_

    use danny_variables_m_
    implicit none 

contains

subroutine splitString(string, split_list, split_size)

    ! 文字列stringを ":" のところで分割して、分割したものをsplit_listに、分割
    ! 数をsplit_sizeにそれぞれ格納して返す

    implicit none 
    character(STRING_BUF_), intent(in) :: string
    character(STRING_BUF_), intent(inout) :: split_list(20)
    integer split_size, i, s, ind, j

    ind = 1
    j = 1

    do i = 1, STRING_BUF_
        if (string(i:i) == ':') then
            s = i - 1
            split_list(j) = trim( adjustl( string(ind:s) ) )
            ind = i + 1
            j = j + 1
        end if
    end do

    split_list(j) = trim( adjustl( string(ind:) ) )
    split_size = j

end subroutine splitString


subroutine readInputFile(input_file)

    ! use global_variables
    implicit none
    character(STRING_BUF_), intent(in) :: input_file
    character(STRING_BUF_) line
    character(STRING_BUF_) line_split(20)
    integer :: file_input = 10
    integer n, err, n_row, split_size, ibool

    open(file_input, file = input_file, status='old')

    ! input file の行数を調べる
    n = 0
    count_row : do
        read(file_input, '(a)', iostat=err) line
        if (err == 0) then
            n = n + 1
        else if (err > 0) then
            exit
        end if 
    end do count_row
    n_row = n
    rewind(file_input)

    ! input file 読み込み
    do n = 1, n_row 

        read(file_input, '(a)') line
        call splitString(line, line_split, split_size)

        if (line_split(1) == "species_file") then 

            species_file_ = line_split(2)

        else if (line_split(1) == "binding_energy_file") then

             binding_energy_file_ = line_split(2)

        else if (line_split(1) == "reaction_file") then

            reaction_file_ = line_split(2)

        else if (line_split(1) == "abundance_file") then

            abundances_file_ = line_split(2)

        else if (line_split(1) == "element_file") then 

            element_file_ = line_split(2)

        else if (line_split(1) == "activation_energy_file") then 
            
            activation_energy_file_ = line_split(2)

        else if (line_split(1) == "enthalpy_file") then 

            enthalpy_file_ = line_split(2)

        else if (line_split(1) == "number_of_gaseous_species") then 

            read(line_split(2), *) number_of_gaseous_species_

        else if (line_split(1) == "number_of_species") then
            
            read(line_split(2), *) number_of_species_

        else if (line_split(1) == "number_of_elements") then 

            read(line_split(2), *) number_of_elements_

        else if (line_split(1) == "number_of_reactions") then

            read(line_split(2), *) number_of_reactions_

        else if (line_split(1) == "number_of_gas_phase_reactions") then

            read(line_split(2), *) number_of_gaseous_reactions_

        else if (line_split(1) == "is_three_phase_reaction") then 
            
            read(line_split(2), *) ibool 
            if (ibool == 1) then
                is_three_phase_reaction_ = .true.
            else if (ibool == 0) then
                is_three_phase_reaction_ = .false.
            else 
                !pass
            endif

        else if (line_split(1) == "is_chemical_desorption") then 

            read(line_split(2), *) ibool
            if (ibool == 1) then
                is_chemical_desorption_ = .true.
            else if (ibool == 0)  then 
                is_chemical_desorption_ = .false.
            else 
                !pass
            end if

        else
            ! pass

        end if 

    end do

    close(file_input)

end subroutine readInputFile


subroutine findIndexSpecies(spe_name, index_of_species)

    ! sname(化学種)のindexを取得する
    ! snameがspecies_name(:)に含めれてなければNOT_FOUND_SPECIESを返す

    implicit none

    character(SPECIES_STRING_BUF_), intent(in) :: spe_name
    integer, intent(out) :: index_of_species

    integer i

    do i = 1, number_of_species_
        if (species_name_(i) .eq. spe_name) then
            index_of_species = i  
            return
        end if  
    end do

    index_of_species = NOT_FOUND_SPECIES_

    return

end subroutine findIndexSpecies


subroutine readElementFile 

    ! 元素ファイルの読み込み

    implicit none 

    integer :: file = 10
    integer i

    open(file, file = element_file_, status='old')

    read(file, *)
    do i = 1, number_of_elements_
        read(file, *) element_name_(i), elemental_mass_(i)
    end do

    close(file)

end subroutine readElementFile


subroutine readSpeciesFile 

    implicit none 

    integer :: file = 10
    integer i

    open(file, file = species_file_, status='old')
    read(file, *)
    do i = 1, number_of_species_
        read(file, *) species_name_(i), species_mass_(i), species_charge_(i), &
            number_of_constituent_element_of_species_(1:number_of_elements_, i)
        species_mass_(i) = species_mass_(i) * ATOMIC_MASS_UNIT
        if (species_name_(i) .eq. "e-") then
            species_mass_(i) = ELECTRON_MASS
        end if 
    end do

    ! Hなどの特定の化学種のインデックスの取得
    do i = 1, number_of_species_
        if (species_name_(i) .eq. "H") then
            index_H_ = i 
        end if 
        if (species_name_(i) == "H2") then
            index_H2_ = i 
        end if 
        if (species_name_(i) .eq. "e-") then
            index_electron_ = i 
        end if 
        if (species_name_(i) == "sH") then 
            index_sH_ = i 
        end if
        if (species_name_(i) == "sH2") then
            index_sH2_ = i 
        end if 
        if (species_name_(i) == "mH") then 
            index_mH_ = i 
        end if 
        if (species_name_(i) == "mH2") then 
            index_mH2_ = i 
        end if
        if (species_name_(i) .eq. "sH2O") then 
            index_sH2O_ = i 
        end if
    end do

end subroutine readSpeciesFile


subroutine readAbundancesFile

    implicit none

    integer :: file = 10
    integer i, err, na, j
    character(20) dummy, sname
    double precision abund

    open(file, file = abundances_file_, status='old')

    ! ファイルの行数のカウント
    na = 0
    count_abundances : do 
        read(file, *, iostat=err) dummy
        if (err == 0) then
            na = na + 1
        else if (err > 0) then
            exit
        end if
    end do count_abundances

    rewind(file)

    ! ファイル読み込み
    read(file, *)
    do i = 1, na-1 
        read(file, *) sname, dummy, abund 
        call findIndexSpecies(sname, j)
        if (j .ne. NOT_FOUND_SPECIES_) then
            species_abundances_(j) = abund 
        end if
    end do

    close(file)

end subroutine readAbundancesFile


subroutine readReactionFile 

    implicit none

    integer :: file = 10
    integer i, j, type_id, reac_id, formula_id, idstart, idend
    integer idx1, idx2, idx3, max_nj, count
    integer use_species_for_reaction(number_of_reactions_, number_of_species_)
    character(SPECIES_STRING_BUF_) dummy, r1, r2, r3, p1, p2, p3, p4, p5
    double precision tmin, tmax

    open(file, file = reaction_file_, status='old')

    do i = 1, number_of_reactions_

        read(file, *) r1, r2, r3, p1, p2, p3, p4, p5, value_for_reaction_(1:5, i), &
            dummy, type_id, tmin, tmax, formula_id, reac_id 

        reaction_species_name_(1, i) = r1
        reaction_species_name_(2, i) = r2
        reaction_species_name_(3, i) = r3
        reaction_species_name_(4, i) = p1
        reaction_species_name_(5, i) = p2
        reaction_species_name_(6, i) = p3 
        reaction_species_name_(7, i) = p4 
        reaction_species_name_(8, i) = p5

        call findIndexSpecies(r1, index_of_reaction_species_(1, i))
        call findIndexSpecies(r2, index_of_reaction_species_(2, i))
        call findIndexSpecies(r3, index_of_reaction_species_(3, i))
        call findIndexSpecies(p1, index_of_reaction_species_(4, i))
        call findIndexSpecies(p2, index_of_reaction_species_(5, i))
        call findIndexSpecies(p3, index_of_reaction_species_(6, i))
        call findIndexSpecies(p4, index_of_reaction_species_(7, i))
        call findIndexSpecies(p5, index_of_reaction_species_(8, i))

        ! CRPの反応はH2の存在量に依存するため、反応率方程式を計算するときのために2番目の反応物のindexをH2にする
        ! if ((type_id == TYPE_ID_2_) .or. (type_id == TYPE_ID_PHOTO_DISS_CR_GRAINS_)) then
        !     if (index_of_reaction_species_(2, i) == NOT_FOUND_SPECIES_) then
        !         index_of_reaction_species_(2, i) = index_H2_
        !     end if
        ! end if

        type_of_reaction_(i) = type_id 
        lower_temperature_limit_(i) = tmin
        upper_temperature_limit_(i) = tmax
        formula_id_of_gas_reaction_(i) = formula_id
        reaction_id_(i) = reac_id

    end do

    close(file)


    ! 各reaction type idのはじめと終わりのindexを取得
    idstart = type_of_reaction_(1)
    type_id_start_(idstart) = 1 !type_of_reaction(1)
    do i = 2, number_of_reactions_
        if (type_of_reaction_(i) .eq. type_of_reaction_(i-1)) then
            cycle
        else 
            type_id_start_(type_of_reaction_(i)) = i 
            type_id_end_(type_of_reaction_(i-1)) = i-1
        end if
    end do
    idend = type_of_reaction_(number_of_reactions_)
    type_id_end_(idend) = number_of_reactions_


    ! 各化学種が関連する反応を取得
    do j = 1, number_of_species_

        count = 1
        do i = 1, number_of_reactions_

            idx1 = index_of_reaction_species_(1, i)
            idx2 = index_of_reaction_species_(2, i)
            idx3 = index_of_reaction_species_(3, i)

            if ((idx1 == j) .or. (idx2 == j) .or. (idx3 == j)) then 
                use_species_for_reaction(count, j) = i 
                count = count + 1
            end if

        end do
        number_of_reaction_species_(j) = count - 1

    end do

    max_nj = maxval(number_of_reaction_species_)
    allocate(index_of_relevant_reaction_of_species_(max_nj, number_of_species_))
    index_of_relevant_reaction_of_species_(:,:) = use_species_for_reaction(1:max_nj, 1:number_of_species_)

end subroutine readReactionFile

subroutine readBindingEnergyFile 

    implicit none

    integer :: file = 10
    integer i, err, na, ig, is, im
    character(20) dummy, sname, s_sname, m_sname
    double precision bd, mass, bd_bare

    open(file, file = binding_energy_file_, status='old')

    ! ファイルの行数のカウント
    na = 0
    count_abundances : do 
        read(file, *, iostat=err) dummy
        if (err == 0) then
            na = na + 1
        else if (err > 0) then
            exit
        end if
    end do count_abundances

    rewind(file)

    ! ファイル読み込み
    do i = 1, na

        read(file, *) sname, mass, bd, bd_bare
        ! '//'は文字列を結合する演算子
        s_sname = 's' // sname
        m_sname = 'm' // sname

        call findIndexSpecies(sname, ig)
        if (ig .ne. NOT_FOUND_SPECIES_) then 
            binding_energy_(ig) = bd
        end if

        ! surface species
        call findIndexSpecies(s_sname, is)
        if (is .ne. NOT_FOUND_SPECIES_) then 
            binding_energy_(is) = bd
            diffusion_barrier_(is) = 0.4d0 * bd
            ! Hのdifusion barrier
            if (s_sname .eq. 'sH') diffusion_barrier_(is) = 230.0d0
            if (s_sname .eq. 'sH2') diffusion_barrier_(is) = 220.0d0 
        end if

        ! mantle species
        call findIndexSpecies(m_sname, im)
        if (im .ne. NOT_FOUND_SPECIES_) then 
            binding_energy_(im) = bd
            diffusion_barrier_(im) = 0.8d0 * bd
            if (sname .eq. 'mH2O') diff_m_H2O_ = diffusion_barrier_(im) 
        end if

    end do  

    close(file)

    ! H, H2, C, N, Oを除くmantle種で、Ediff_m < Ediff_m(H2O)となるものは、Ediff_m = Ediff_m(H2O)とする。
    do i = number_of_gaseous_species_+1, number_of_species_
        if ((species_name_(i)(1:1) .eq. 'm') .and. (binding_energy_(i) .ne. 0.0) & 
            .and. (diffusion_barrier_(i) < diff_m_H2O_)) then

                if ((species_name_(i) .ne. 'mH') .and. (species_name_(i) .ne. 'mH2') .and. &
                    (species_name_(i) .ne. 'mC') .and. (species_name_(i) .ne. 'mN') .and. &
                    (species_name_(i) .ne. 'mO')) then

                        diffusion_barrier_(i) = diff_m_H2O_
                        !diffusion_barrier_bare(i) = diff_m_H2O

                end if

        end if
    end do


    ! calculate vibrational frequency
    do i = number_of_gaseous_species_+1, number_of_species_
        if (binding_energy_(i) .ne. 0.0) then
            vibrational_frequency_(i) = sqrt(2.0*BINDING_SITES_DENSITY_*BOLTZMANN_CONSTANT* &
                binding_energy_(i) / (PI*PI*species_mass_(i)))
        end if
    end do

end subroutine readBindingEnergyFile


subroutine readSurfaceActivationEnergyFile 

    implicit none 

    integer :: file = 10
    integer i, err, na, j, ir1, ir2, ir3, ip1, ip2, ip3, ip4, ip5
    character(SPECIES_STRING_BUF_) dummy, r1, r2, r3, p1, p2, p3, p4, p5
    double precision e_ac, reduced_mass

    open(file, file = activation_energy_file_, status='old')

    ! ファイルの行数のカウント
    na = 0
    count_activation : do 
        read(file, *, iostat=err) dummy
        if (err == 0) then
            na = na + 1
        else if (err > 0) then
            exit
        end if
    end do count_activation

    rewind(file)

    do i = 1, na 

        read(file, *) r1, r2, r3, p1, p2, p3, p4, p5, e_ac
        call findIndexSpecies(r1, ir1)
        call findIndexSpecies(r2, ir2)
        call findIndexSpecies(r3, ir3)
        call findIndexSpecies(p1, ip1)
        call findIndexSpecies(p2, ip2)
        call findIndexSpecies(p3, ip3)
        call findIndexSpecies(p4, ip4)
        call findIndexSpecies(p5, ip5)

        do j = number_of_gaseous_reactions_, number_of_reactions_
            if ((ir1 == index_of_reaction_species_(1, j)) .and. ir2 == index_of_reaction_species_(2, j) .and. &
                (ir3 == index_of_reaction_species_(3, j)) .and. ip1 == index_of_reaction_species_(4, j) .and. &
                (ip2 == index_of_reaction_species_(5, j)) .and. ip3 == index_of_reaction_species_(6, j) .and. &
                (ip4 == index_of_reaction_species_(7, j)) .and. ip5 == index_of_reaction_species_(8, j)) then 
                    
                    ! count = count + 1
                    ! write(*,*) count, r1, r2, r3, p1, p2, p3, p4, p5, e_ac
                    value_for_reaction_(1, j) = e_ac
                    reduced_mass = species_mass_(ir1) * species_mass_(ir2) / (species_mass_(ir1) + species_mass_(ir2))
                    value_for_reaction_(2, j) = 2.0d0*WIDTH_OF_BARRIER_/DIRAC_CONSTANT * &
                        sqrt(2.0d0*reduced_mass*BOLTZMANN_CONSTANT*e_ac)
            end if 
        end do

    end do

end subroutine readSurfaceActivationEnergyFile


subroutine readEnthalpyOfFormationFile

    implicit none

    integer :: file = 10
    integer i, err, na, idx
    character(SPECIES_STRING_BUF_) dummy, sname, s_sname, m_sname
    double precision enthalpy

    open(file, file = enthalpy_file_, status='old')

    na = 0
    count_activation : do 
        read(file, *, iostat=err) dummy
        if (err == 0) then
            na = na + 1
        else if (err > 0) then
            exit
        end if
    end do count_activation

    rewind(file)

    do i = 1, na 
        read(file, *) sname, enthalpy
        s_sname = 's' // sname
        m_sname = 'm' // sname
        call findIndexSpecies(sname, idx)
        if (idx .ne. NOT_FOUND_SPECIES_) enthalpy_of_formation_(idx) = enthalpy
        call findIndexSpecies(s_sname, idx)
        if (idx .ne. NOT_FOUND_SPECIES_) enthalpy_of_formation_(idx) = enthalpy
        call findIndexSpecies(m_sname, idx)
        if (idx .ne. NOT_FOUND_SPECIES_) enthalpy_of_formation_(idx) = enthalpy
    end do

end subroutine


! #################################################################################################
! #################################################################################################


subroutine checkReadSpecies

    use danny_variables_m_
    implicit none 

    character(STRING_BUF_) :: filename = "check_read_species.txt"
    integer :: file = 14
    integer i

    open(file, file = filename, status='replace')

    do i = 1, number_of_species_
        write(file, *) species_name_(i), species_mass_(i), species_charge_(i)
    end do

    close(file)

end subroutine checkReadSpecies


subroutine checkReadReaction

    use danny_variables_m_
    implicit none 

    character(STRING_BUF_) :: filename = "check_read_reaction.txt"
    integer :: file = 14
    integer i, ir1, ir2, ir3, ip1, ip2, ip3, ip4, ip5
    character(SPECIES_STRING_BUF_) r1, r2, r3, p1, p2, p3, p4, p5

    open(file, file = filename, status='replace')

    ! do i = 1, number_of_reactions_
    do i = type_id_start_(TYPE_ID_GRAIN_SURFACE_REACTION_), type_id_end_(TYPE_ID_GRAIN_SURFACE_REACTION_)

        ir1 = index_of_reaction_species_(1, i)
        ir2 = index_of_reaction_species_(2, i)
        ir3 = index_of_reaction_species_(3, i)
        ip1 = index_of_reaction_species_(4, i)
        ip2 = index_of_reaction_species_(5, i)
        ip3 = index_of_reaction_species_(6, i)
        ip4 = index_of_reaction_species_(7, i)
        ip5 = index_of_reaction_species_(8, i)

        r1 = species_name_(ir1)

        if (ir2 .eq. NOT_FOUND_SPECIES_) then
            if (type_of_reaction_(i) .eq. TYPE_ID_1_) then 
                r2 = "CR"
            else if (type_of_reaction_(i) .eq. TYPE_ID_2_) then
                r2 = "CRP"
            else if (type_of_reaction_(i) .eq. TYPE_ID_3_) then
                r2 = "Photon"
            else
                r2 = NO_SPECIES_
            end if
        
        else 
            r2 = species_name_(ir2)
        end if

        if (ir3 .eq. NOT_FOUND_SPECIES_) then 
            r3 = NO_SPECIES_
        else 
            r3 = species_name_(ir3)
        end if
        
        if (ip1 .eq. NOT_FOUND_SPECIES_) then 
            p1 = NO_SPECIES_
        else 
            p1 = species_name_(ip1)
        end if

        if (ip2 .eq. NOT_FOUND_SPECIES_) then
            p2 = NO_SPECIES_
        else 
            p2 = species_name_(ip2)
        end if

        if (ip3 .eq. NOT_FOUND_SPECIES_) then 
            p3 = NO_SPECIES_
        else 
            p3 = species_name_(ip3)
        end if

        if (ip4 .eq. NOT_FOUND_SPECIES_) then 
            p4 = NO_SPECIES_
        else 
            p4 = species_name_(ip4)
        end if

        if (ip5 .eq. NOT_FOUND_SPECIES_) then 
            p5 = NO_SPECIES_
        else 
            p5 = species_name_(ip5)
        end if

        write(file, *) r1, r2, r3, p1, p2, p3, p4, p5, branching_ratio_(i), &
            value_for_reaction_(3, i), value_for_reaction_(4, i)
        ! write(*,*) r1, r2, r3, p1, p2, p3, p4, p5
        ! write(*, *) ir1, ir2, ir3, ip1
    end do

    close(file) 

end subroutine checkReadReaction

end module danny_read_m_