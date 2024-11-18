import numpy as np 

input_dir = "input/"
species_file = input_dir + "species.dat"
reaction_file = input_dir + "reactions.dat"

grain_charge_min = -1
grain_charge_max = 1

isGrainChargedGasSpeciesReaction = True # not include grain surface species
isGrainReaction = True                  # include grain surface species
isGrainSurfaceReaction = True
isChemicalDesorption = True
isH2Desorption = True                   # ダスト表面のH2が脱着する反応
isThreePhaseReaction = False            # gas, grain surface, grain mantleの3相を考慮するかどうか

###########################################################################################################
###########################################################################################################

## filename
gas_species_file1    = "chemical_data/gas_species.dat"
gas_reactions_file1  = "chemical_data/kida.uva.2014.dat"
# gas_species_file1    = "chemical_data2/gas_species.dat"
# gas_reactions_file1  = "chemical_data2/gas_reaction.dat"
binding_energy_file  = "chemical_data/binding_energy.dat"
# binding_energy_file  = "Ruaud2016/binding_energy.dat"
grain_reaction_file  = "chemical_data/kida_surface_reaction.dat"
# grain_reaction_file  = "Ruaud2016/surface_reaction.dat"
abundance_file       = "chemical_data/abundances.dat"
# element_file         = "chemical_data/element.dat"
element_file         = "chemical_data/element.dat"
activate_energy_file = "chemical_data/surface_activation_energy.dat"
enthalpy_file        = "chemical_data/enthalpy_of_formation.dat"

###########################################################################################################
###########################################################################################################

## type_id
type_id_grain_charged_gas_species_reaction = 15
type_id_accretion_gas_species = 16
type_id_thermal_desorption = 17
type_id_cosmic_ray_desorption = 18
type_id_photo_desorption_by_UV = 19
type_id_photo_desorption_by_UV_CR = 20
type_id_grain_surface_reaction = 21
type_id_grain_mantle_reaction  = 22
type_id_photo_dissociation_by_CR_on_grains = 23
type_id_photo_dissociation_by_UV_on_grains = 24
type_id_surface_to_mantle = 25
type_id_mantle_to_surface = 26

###########################################################################################################
###########################################################################################################

element_file_input = input_dir + "element.dat"
file = open(element_file_input, mode="w")
element_list = []
count_element = 0
with open(element_file) as f:
    for line in f:
        file.write(line)
        if not (line[0] == "#"):
            count_element += 1
            line = line.split(" ")
            line = [x for x in line if not x == ""]
            line[-1] = line[-1].replace("\n", "")
            line[-1] = float(line[-1])
            element_list.append(line)

file.close()

###########################################################################################################
###########################################################################################################


## read and write species file

gas_species_list1 = np.loadtxt(gas_species_file1, str).tolist()
gas_species_list_all = gas_species_list1
gas_species_list = []
gas_species = []
gas_species_mass = []
charged_gas_species = []
dust_grain_species = []

count_gas_species = 0
count_surface_species = 0
count_mantle_species = 0
count_grain_species = 0

file = open(species_file, mode="w")

def set_species_write(species_list):
    """
    species_listをファイルに書き込みための文字列に直す
    """
    nlen = len(species_list[0])
    if (len(species_list[1]) == 1):
        blanks = 12 - nlen
    elif (len(species_list[1]) == 2):
        blanks = 14 - nlen - len(species_list[1])
    else:
        blanks = 15 - nlen - len(species_list[1])
    s = species_list[0] + " "*blanks

    # mass
    blanks = 4 - len(species_list[1])
    s += " "*blanks + species_list[1]

    # charge
    blanks = 8 - len(species_list[2])
    s += " "*blanks + species_list[2]

    for i in range(13):
        blanks = 5 - len(species_list[i+3])
        s += " "*blanks + species_list[i+3]
    s += "\n"
    return s

## input species fileの一番最初の行への記述
file.write("#name       mass  charge")
for ii in range(len(element_list)):
    blanks = 5 - len(element_list[ii][0])
    file.write(" "*blanks + element_list[ii][0])
file.write("\n")


# GRAINをリストから除く
for spe in gas_species_list_all:
    if spe[0][:5] == "GRAIN":
        pass
    else:
        gas_species_list.append(spe)

## gas phase species
for spe in gas_species_list:
    ## calculate species mass
    mass = 0
    ii = 0
    for si in spe[2:]:
        mass += float(si) * element_list[ii][1]
        ii += 1
    spe.insert(1, str(int(mass)))
    s = set_species_write(spe)
    count_gas_species += 1
    file.write(s)
    gas_species.append(spe[0])
    gas_species_mass.append([spe[0], spe[1]])
    ## charged species
    if not (spe[2] == "0"):
        charged_gas_species.append([spe[0], spe[1], spe[2]]) ## name mass charge


grain_surface_species = []
if isGrainReaction:

    # grain surface species
    grain_species_list = np.loadtxt(binding_energy_file, str).tolist()
    temp_grain_species = np.loadtxt(binding_energy_file, str).T[0].tolist()
    for spe in temp_grain_species:
        if spe in gas_species:
            ind = gas_species.index(spe)
            slist = gas_species_list[ind].copy()
            if (slist[0] == "GRAIN-"):
                print(ind, spe, slist)
            slist[0] = "s" + slist[0]
            #slist.insert(1, gas_species_mass[ind][1])
            s = set_species_write(slist)
            file.write(s)
            count_surface_species += 1
            grain_surface_species.append(spe)
        else:
            print(spe)
            pass


    if isThreePhaseReaction:
        ## grain mantle species
        for spe in grain_surface_species:
            if spe in gas_species:
                ind = gas_species.index(spe)
                slist = gas_species_list[ind].copy()
                slist[0] = "m" + slist[0]
                #slist.insert(1, gas_species_mass[ind][1])
                s = set_species_write(slist)
                file.write(s)
                count_mantle_species += 1
            else:
                print(spe)
                pass

## dust grain species
if (isGrainChargedGasSpeciesReaction or isGrainReaction):
    for c in range(grain_charge_min, grain_charge_max+1):
        if c <= 0:
            gname = "G" + "c" + str(c)
        else:
            gname = "Gc+" + str(c)
        glist = [gname, "0", str(c)]
        for i in range(13):
            glist.append("0")
        s = set_species_write(glist)
        file.write(s)
        count_grain_species += 1
        dust_grain_species.append([gname, str(c)])


file.close()

total_species = count_gas_species + count_surface_species + count_mantle_species + count_grain_species
#print(count_gas_species, count_surface_species, count_mantle_species, count_grain_species, total_species)
print("gas species     = ", count_gas_species)
print("surface species = ", count_surface_species)
print("mantle species  = ", count_mantle_species)
print("grain species   = ", count_grain_species)
print("total species   = ", total_species)

###########################################################################################################
###########################################################################################################

## read and write gas rhase reaction file
number_of_gas_reaction_type = 12
gas_reaction_type = [[] for i in range(number_of_gas_reaction_type)]

reaction_CPR = []          # Photo-processes (Photo-dissociation) induced by cosmic rays (secondary photons)
reaction_PHOTON = []       # Photo-processes (Photo-dissociation) by UV photons with a standard interstellar UV field
total_reactions = 0
total_grain_reactions = 0

def position_check(alist):
    '''
    alist : 化学反応の関する化学種や反応係数に用いる値が格納されたリスト

    alistにおいて、反応・生成に関わる化学種の最後の位置+1の値を返す
    '''
    nlen = len(alist)
    for i in range(nlen):
        #print(alist[i])
        try:
            b = float(alist[i])
            return i
        except:
            continue
    return -1

file = open(reaction_file, mode="w")

def set_reaction_write(reaction_list):
    '''
    reaction_listをファイルに書き出すための文字列に直す。
    '''
    s = reaction_list[0]
    for i in range(1, 8):
        blanks = 12 - len(reaction_list[i-1])
        s += " "*blanks + reaction_list[i]
    for i in range(13):
        # if i == 0:
        #     blanks = 5
        # else:
        #     blanks = 2
        blanks = 12 - len(reaction_list[i+8])
        s += " "*blanks + reaction_list[i+8]
    #s.replace("\n", "")
    #s += "\n"
    return s


#### gas phase reactions ####
count_gas_phase_reactions = 0

ast = ["*", "CR", "CRP", "Photon"]

with open(gas_reactions_file1) as f:
    for line in f:
        # １行読み込み
        line = line.split(" ")
        line = [x for x in line if not x == ""]
        if (line[0][0] == "!") or (line[0][0] == "#"):
            continue
        nn = position_check(line) 
        line.insert(2, "*")
        for i in range(7-nn):
            line.insert(i+nn+1, "*")
        reac = line[:8]

        # 化学反応において、gas_speciesに含まれていない化学種が存在sする場合は含めないようにする
        temp = 0
        for r in reac:
            if not ((r in gas_species) or (r in ast)):
                temp = 1
                break
        if (temp == 1):
            #print(line)
            continue

        # type_id毎のlistに格納
        type_id = int(line[-7])
        gas_reaction_type[type_id].append(line)

# reaction fileへのの書き込み
for id in range(1, number_of_gas_reaction_type):
    num_reac = len(gas_reaction_type[id])
    if (num_reac == 0):
        continue
    else:
        for rlist in gas_reaction_type[id]:
            rlist[-1] = rlist[-1].replace("\n", "")
            rlist[-1] += "\n"
            s = set_reaction_write(rlist)
            file.write(s)
            count_gas_phase_reactions += 1

# grain surface speciesも考える場合に必要なPhoto-process reactionsの取得
if isGrainReaction:
    # Photo-processes reactions 
    for rlist in gas_reaction_type[2]:
        reac = rlist[:8]
        tempi = 0
        for s in reac:
            if (s == "*") or (s == "CRP") or (s in grain_surface_species):
                pass
            else:
                tempi += 1
                break
        if tempi == 0:
            reaction_CPR.append(rlist)

    for rlist in gas_reaction_type[3]:
        reac = rlist[:8]
        tempi = 0
        for s in reac:
            if (s == "*") or (s == "Photon") or (s in grain_surface_species):
                pass
            else:
                tempi += 1
        if tempi == 0:
            reaction_PHOTON.append(rlist)

total_reactions += count_gas_phase_reactions
# print(count_gas_phase_reactions, total_reactions)
print("gas phase reactions = ", count_gas_phase_reactions)

###################################################
#### grains and (charged) gas species reaction ####
###################################################
count_grain_charged_gas_reactions = 0
if isGrainChargedGasSpeciesReaction:
    len_dust_grain_species = len(dust_grain_species)
    for spe in charged_gas_species:
        sname = spe[0]
        sch = float(spe[2])
        if (sname == "e-"):
            for i in range(1, len_dust_grain_species):
                slist = [sname, dust_grain_species[i][0], "*", dust_grain_species[i-1][0], "*", "*", "*", "*", spe[1], spe[2], dust_grain_species[i][1]] # mass charge_gas charge_dust
                for i in range(9): # 12
                    slist.append("0")
                slist.append("0\n")
                slist[-7] = str(type_id_grain_charged_gas_species_reaction)
                s = set_reaction_write(slist)
                file.write(s)
                count_grain_charged_gas_reactions += 1
        elif (sname[:-1] in gas_species):
            if (sch < 0):
                for i in range(1, len_dust_grain_species):
                    slist = [sname, dust_grain_species[i][0], "*", dust_grain_species[i-1][0], sname[:-1], "*", "*", "*", spe[1], spe[2], dust_grain_species[i][1]]
                    for i in range(9):
                        slist.append("0")
                    slist.append("0\n")
                    slist[-7] = str(type_id_grain_charged_gas_species_reaction)
                    s = set_reaction_write(slist)
                    file.write(s)
                    count_grain_charged_gas_reactions += 1
            elif (sch > 0):
                for i in range(0, len_dust_grain_species-1):
                    slist = [sname, dust_grain_species[i][0], "*", dust_grain_species[i+1][0], sname[:-1], "*", "*", "*", spe[1], spe[2], dust_grain_species[i][1]]
                    for i in range(9):
                        slist.append("0")
                    slist.append("0\n")
                    slist[-7] = str(type_id_grain_charged_gas_species_reaction)
                    s = set_reaction_write(slist)
                    file.write(s)
                    count_grain_charged_gas_reactions += 1
        else:
            # 今後実装する
            pass 
            #print(sname)


    total_grain_reactions += count_grain_charged_gas_reactions
    total_reactions += count_grain_charged_gas_reactions
    # print(count_grain_charged_gas_reactions, total_grain_reactions, total_reactions)
    print("grains and charged gas species reaction = ", count_grain_charged_gas_reactions)

###################################################
####          Grain Surface Reaction           ####
###################################################
if isGrainReaction:
    #### accretion (neutral species) ####
    count_accretion = 0
    for spelist in grain_species_list:
        spe = spelist[0]
        gspe = "s" + spelist[0]
        mass = spelist[1] # species mass (amu)
        be = spelist[2]  # binding energy
        try:
            ind = gas_species.index(spe)
        except:
            continue
        slist = [spe, "*", "*", gspe, "*", "*", "*", "*", mass]
        for i in range(11):
            slist.append("0")
        slist.append("0\n")
        slist[-7] = str(type_id_accretion_gas_species)
        s = set_reaction_write(slist)
        file.write(s)
        count_accretion += 1

    total_grain_reactions += count_accretion
    total_reactions += count_accretion
    # print(count_accretion, total_grain_reactions, total_reactions)
    print("accretion on grains = ", count_accretion)

    #### thermal desorption ####
    count_thermal_desorption = 0
    for spe in grain_surface_species:
        if spe in gas_species:
            ind1 = gas_species.index(spe)
            ind2 = grain_surface_species.index(spe)
            mass = gas_species_mass[ind1][1]
            bd = grain_species_list[ind2][2]
            gspe = "s" + spe
            # slist = [gspe, "*", "*", spe, "*", "*", "*", "*", mass, bd]
            slist = [gspe, "*", "*", spe, "*", "*", "*", "*", "0", "0"]
            for i in range(10):
                slist.append("0")
            slist.append("0\n")
            slist[-7] = str(type_id_thermal_desorption)
            s = set_reaction_write(slist)
            file.write(s)
            count_thermal_desorption += 1

    print("tharmal desorption  = ", count_thermal_desorption)
    total_grain_reactions += count_thermal_desorption

    #### cosmic-ray desoption ####
    count_CR_desorption = 0
    for spe in grain_surface_species:
        if spe in gas_species:
            ind1 = gas_species.index(spe)
            ind2 = grain_surface_species.index(spe)
            mass = gas_species_mass[ind1][1]
            bd = grain_species_list[ind2][2]
            gspe = "s" + spe
            # slist = [gspe, "*", "*", spe, "*", "*", "*", "*", mass, bd]
            slist = [gspe, "*", "*", spe, "*", "*", "*", "*", "0", "0"]
            for i in range(10):
                slist.append("0")
            slist.append("0\n")
            slist[-7] = str(type_id_cosmic_ray_desorption)
            s = set_reaction_write(slist)
            file.write(s)
            count_CR_desorption += 1

    print("cosmiac rays desorption = ", count_CR_desorption)
    total_grain_reactions += count_CR_desorption

    #### photo desoption by UV ####
    count_photo_desorption_by_UV = 0
    for spe in grain_surface_species:
        gspe = "s" + spe
        slist = [gspe, "*", "*", spe, "*", "*", "*", "*"]
        for i in range(12):
            slist.append("0")
        slist.append("0\n")
        slist[-7] = str(type_id_photo_desorption_by_UV)
        s = set_reaction_write(slist)
        file.write(s)
        count_photo_desorption_by_UV += 1

    print("photo desorption by UV = ", count_photo_desorption_by_UV)
    total_grain_reactions += count_photo_desorption_by_UV

    #### photo desoption by UV_CR ####
    count_photo_desorption_by_UV_CR = 0
    for spe in grain_surface_species:
        gspe = "s" + spe
        slist = [gspe, "*", "*", spe, "*", "*", "*", "*"]
        for i in range(12):
            slist.append("0")
        slist.append("0\n")
        slist[-7] = str(type_id_photo_desorption_by_UV_CR)
        s = set_reaction_write(slist)
        file.write(s)
        count_photo_desorption_by_UV_CR += 1

    print("photo desorption by UV_CR = ", count_photo_desorption_by_UV_CR)
    total_grain_reactions += count_photo_desorption_by_UV_CR

    if isGrainSurfaceReaction:
        #### grain surface reaction (two body) ####
        count_grain_surface_reactions = 0
        linep = []
        with open(grain_reaction_file) as f:
            for line in f:
                line = line.split(" ")
                line = [x for x in line if not x == ""]
                if (line[0][0] == "!") or (line[0][0] == "#"):
                    continue

                ## "*"の挿入
                nn = position_check(line)
                line.insert(2, "*")
                for i in range(7-nn):
                    line.insert(i+nn+1, "*")
                ## 重複する反応を除く
                if not (count_grain_surface_reactions == 0):
                    if line[:8] == linep[:8]:
                        continue
                linep = line.copy()

                ## 反応の化学種がgrain_speciesに含まれているか確認
                spe = line[:8] # "s"をつける前の化学種反応種のlist
                int_temp = 0
                for i in range(len(line[:8])):
                    if line[i] == "*":
                        pass
                    elif line[i] in grain_surface_species:
                        line[i] = "s" + line[i]
                        pass
                    else:
                        int_temp += 1
                        continue
                if int_temp != 0:
                    continue

                ## make reaction list
                idx1 = grain_surface_species.index(spe[0])
                idx2 = grain_surface_species.index(spe[1])
                mass1 = grain_species_list[idx1][1]
                bd1 = grain_species_list[idx1][2]
                mass2 = grain_species_list[idx2][1]
                bd2 = grain_species_list[idx2][2]
                # slist = [line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], mass1, bd1, mass2, bd2]
                slist = [line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], "0", "0", "0", "0"]
                for i in range(8):
                    slist.append("0")
                slist.append("0\n")
                slist[-7] = str(type_id_grain_surface_reaction)

                ## ファイル書き込み
                s = set_reaction_write(slist)
                file.write(s)
                count_grain_surface_reactions += 1
                
                # 化学脱着
                if isChemicalDesorption:
                    if (slist[4] == "*") and (slist[5] == "*"):
                        # sA + sB → C
                        slist_cd = slist.copy()
                        slist_cd[3] = slist_cd[3][1:]
                        s = set_reaction_write(slist_cd)
                        file.write(s)
                        count_grain_surface_reactions += 1

                    elif slist[5] == "*":
                        # sA + sB → sC + D
                        slist_cd = slist.copy()
                        slist_cd[4] = slist_cd[4][1:]
                        s = set_reaction_write(slist_cd)
                        file.write(s)
                        count_grain_surface_reactions += 1
                        # sA + sB → C + sD
                        slist_cd = slist.copy()
                        slist_cd[3] = slist_cd[3][1:]
                        s = set_reaction_write(slist_cd)
                        file.write(s)
                        count_grain_surface_reactions += 1

                        # sA + sB → C + D
                        slist_cd[4] = slist_cd[4][1:]
                        s = set_reaction_write(slist_cd)
                        file.write(s)
                        count_grain_surface_reactions += 1

                    else:
                        # sA + sB → sC + sD + E
                        slist_cd = slist.copy()
                        slist_cd[5] = slist_cd[5][1:]
                        s = set_reaction_write(slist_cd)
                        file.write(s)
                        count_grain_surface_reactions += 1
                        # sA + sB → sC + D + sE
                        slist_cd = slist.copy()
                        slist_cd[4] = slist_cd[4][1:]
                        s = set_reaction_write(slist_cd)
                        file.write(s)
                        count_grain_surface_reactions += 1
                        # sA + sB → sC + D +  E
                        slist_cd[5] = slist_cd[5][1:]
                        s = set_reaction_write(slist_cd)
                        file.write(s)
                        count_grain_surface_reactions += 1
                        # sA + sB →  C + sD + sE
                        slist_cd = slist.copy()
                        slist_cd[3] = slist_cd[3][1:]
                        s = set_reaction_write(slist_cd)
                        file.write(s)
                        count_grain_surface_reactions += 1
                        # sA + sB →  C + sD +  E
                        slist_cd[5] = slist_cd[5][1:]
                        s = set_reaction_write(slist_cd)
                        file.write(s)
                        count_grain_surface_reactions += 1
                        # sA + sB →  C +  D + sE
                        slist_cd = slist.copy()
                        slist_cd[3] = slist_cd[3][1:]
                        slist_cd[4] = slist_cd[4][1:]
                        s = set_reaction_write(slist_cd)
                        file.write(s)
                        count_grain_surface_reactions += 1
                        # sA + sB →  C +  D + E
                        slist_cd[5] = slist_cd[5][1:]
                        s = set_reaction_write(slist_cd)
                        file.write(s)
                        count_grain_surface_reactions += 1

            if isH2Desorption:
                # ダスト表面のH2が脱着する反応
                # sH2 + sH2 → sH2 + H2
                slist = ["sH2", "sH2", "*", "sH2", "H2", "*", "*", "*", "0", "0", "0", "0"]
                for i in range(8):
                    slist.append("0")
                slist.append("0\n")
                slist[-7] = str(type_id_grain_surface_reaction)

                ## ファイル書き込み
                s = set_reaction_write(slist)
                file.write(s)
                count_grain_surface_reactions += 1

        print("grains surface reactions = ", count_grain_surface_reactions)
        total_grain_reactions += count_grain_surface_reactions

    if isThreePhaseReaction:
        #### grain mantle reaction (two body) : same grain surface reaction ####
        count_grain_mantle_reactions = 0
        linep = []
        with open(grain_reaction_file) as f:
            for line in f:
                line = line.split(" ")
                line = [x for x in line if not x == ""]
                if (line[0][0] == "!") or (line[0][0] == "#"):
                    continue

                ## "*"の挿入
                nn = position_check(line)
                line.insert(2, "*")
                for i in range(7-nn):
                    line.insert(i+nn+1, "*")
                ## 重複する反応を除く
                if not (count_grain_mantle_reactions == 0):
                    if line[:8] == linep[:8]:
                        continue
                linep = line.copy()

                ## 反応の化学種がgrain_speciesに含まれているか確認
                spe = line[:8]
                int_temp = 0
                for i in range(len(line[:8])):
                    if line[i] == "*":
                        pass
                    elif line[i] in grain_surface_species:
                        line[i] = "m" + line[i]
                        pass
                    else:
                        int_temp += 1
                        continue
                if int_temp != 0:
                    continue

                ## make reaction list
                idx1 = grain_surface_species.index(spe[0])
                idx2 = grain_surface_species.index(spe[1])
                mass1 = grain_species_list[idx1][1]
                bd1 = grain_species_list[idx1][2]
                mass2 = grain_species_list[idx2][1]
                bd2 = grain_species_list[idx2][2]
                # slist = [line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], mass1, bd1, mass2, bd2]
                slist = [line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], "0", "0", "0", "0"]
                for i in range(8):
                    slist.append("0")
                slist.append("0\n")
                slist[-7] = str(type_id_grain_mantle_reaction)

                ## ファイル書き込み
                s = set_reaction_write(slist)
                file.write(s)
                count_grain_mantle_reactions += 1

        print("grain mantle reactions = ", count_grain_mantle_reactions)
        total_grain_reactions += count_grain_mantle_reactions

    if isGrainSurfaceReaction:
        #### Photo-processes (Photo-dissociation) induced by cosmic rays (secondary photons) on grains ####
        count_photodiss_by_CR_on_grains = 0
        # surface
        for re in reaction_CPR:
            r = re.copy()
            for i in range(8):
                if (r[i] == "*") or (r[i] == "CRP"):
                    pass
                else:
                    r[i] = "s" + r[i]
            r[-7] = str(type_id_photo_dissociation_by_CR_on_grains)
            s = set_reaction_write(r)
            file.write(s)
            count_photodiss_by_CR_on_grains += 1

        if isThreePhaseReaction:
            # mantle
            for re in reaction_CPR:
                r = re.copy()
                for i in range(8):
                    if (r[i] == "*") or (r[i] == "CRP"):
                        pass
                    else:
                        r[i] = "m" + r[i]
                r[-7] = str(type_id_photo_dissociation_by_CR_on_grains)
                s = set_reaction_write(r)
                file.write(s)
                count_photodiss_by_CR_on_grains += 1

        print("photdiss by CRs on grains = ", count_photodiss_by_CR_on_grains)
        total_grain_reactions += count_photodiss_by_CR_on_grains

    if isGrainSurfaceReaction:
        #### Photo-processes (Photo-dissociation) by UV photons with a standard interstellar UV field on grains ####
        count_photodiss_by_UR_on_grains = 0
        #surface
        for re in reaction_PHOTON:
            r = re.copy()
            for i in range(8):
                if (r[i] == "*") or (r[i] == "Photon"):
                    pass
                else:
                    r[i] = "s" + r[i]
            r[-7] = str(type_id_photo_dissociation_by_UV_on_grains)
            s = set_reaction_write(r)
            file.write(s)
            count_photodiss_by_UR_on_grains += 1
        
        if isThreePhaseReaction:
            #mantle species
            for re in reaction_PHOTON:
                r = re.copy()
                for i in range(8):
                    if (r[i] == "*") or (r[i] == "Photon"):
                        pass
                    else:
                        r[i] = "m" + r[i]
                r[-7] = str(type_id_photo_dissociation_by_UV_on_grains)
                s = set_reaction_write(r)
                file.write(s)
                count_photodiss_by_UR_on_grains += 1

        print("photdiss by UR on grains = ", count_photodiss_by_UR_on_grains)
        total_grain_reactions += count_photodiss_by_UR_on_grains

    if (isGrainSurfaceReaction and isThreePhaseReaction):
        #### surface to mantle ####
        count_surface_to_mantle = 0
        for spe in grain_surface_species:
            gspe = "s" + spe
            mspe = "m" + spe
            slist = [gspe, "*", "*", mspe, "*", "*", "*", "*"]
            for i in range(12):
                slist.append("0")
            slist.append("0\n")
            slist[-7] = str(type_id_surface_to_mantle)
            s = set_reaction_write(slist)
            file.write(s)
            count_surface_to_mantle += 1

        print("surface to mantle = ", count_surface_to_mantle)
        total_grain_reactions += count_surface_to_mantle

        #### mantle to surface ####
        count_mantle_to_surface = 0
        for spe in grain_surface_species:
            gspe = "s" + spe
            mspe = "m" + spe
            slist = [mspe, "*", "*", gspe, "*", "*", "*", "*"]
            for i in range(12):
                slist.append("0")
            slist.append("0\n")
            slist[-7] = str(type_id_mantle_to_surface)
            s = set_reaction_write(slist)
            file.write(s)
            count_mantle_to_surface += 1

        print("mantle to surface = ", count_mantle_to_surface)
        total_grain_reactions += count_mantle_to_surface

file.close()

########################################################################################
########################################################################################

file = open("input/binding_energy.dat", mode="w")

count_binding_energy = 0
with open(binding_energy_file) as f:
    for line in f:
        file.write(line)
        if not (line[0] == "#"):
            count_binding_energy += 1

file.close()

########################################################################################
########################################################################################

binding_energy_file_out = "input/binding_energy.dat"
file = open(binding_energy_file_out, mode="w")

count_binding_energy = 0
with open(binding_energy_file) as f:
    for line in f:
        if not (line[0] == "#"):
            file.write(line)
            count_binding_energy += 1

file.close()

########################################################################################
########################################################################################

abundance_file_input = input_dir + "abundances.dat"
file = open(abundance_file_input, mode="w")
count_abundances = 0
with open(abundance_file) as f:
    for line in f:
        file.write(line)
        if not (line[0] == "#"):
            count_abundances += 1

file.close()

########################################################################################
########################################################################################


## surface activate energy
activate_energy_file_input = input_dir + "surface_activation_energy.dat"
if isGrainSurfaceReaction:
    #file = open(activate_energy_file_input, mode="w")
    with open(activate_energy_file_input, mode='w') as fout:
        with open(activate_energy_file) as f:
            for line in f:
                sline = line.split()
                temp = 0
                for i in range(8):
                    if sline[i] == "*":
                        pass
                    else:
                        if sline[i] in grain_surface_species:
                            sline[i] = "s" + sline[i]
                        else:
                            temp += 1
                if temp == 0:
                    print("{:<15}".format(sline[0]), "{:<15}".format(sline[1]), "{:<15}".format(sline[2]), "{:<15}".format(sline[3]),
                        "{:<15}".format(sline[4]), "{:<15}".format(sline[5]), "{:<15}".format(sline[6]), "{:<15}".format(sline[7]), 
                        "{:<15}".format(sline[8]), end="\n",file=fout)
            
        if isThreePhaseReaction:
            with open(activate_energy_file) as f:
                for line in f:
                    sline = line.split()
                    temp = 0
                    for i in range(8):
                        if sline[i] == "*":
                            pass
                        else:
                            if sline[i] in grain_surface_species:
                                sline[i] = "m" + sline[i]
                            else:
                                temp += 1
                    if temp == 0:
                        print("{:<15}".format(sline[0]), "{:<15}".format(sline[1]), "{:<15}".format(sline[2]), "{:<15}".format(sline[3]),
                            "{:<15}".format(sline[4]), "{:<15}".format(sline[5]), "{:<15}".format(sline[6]), "{:<15}".format(sline[7]), 
                            "{:<15}".format(sline[8]), end="\n",file=fout)

        

# file.close()

########################################################################################
########################################################################################

## enthalpy of formation
enthalpy_file_input = input_dir + "enthalpy_of_formation.dat"
with open(enthalpy_file_input, mode="w") as file:
    with open(enthalpy_file, mode="r") as f:
        for line in f:
            if line[0] == "#":
                continue
            else:
                line = line.replace("\n", "")
                print(line, end="\n", file=file)

########################################################################################
########################################################################################

total_reactions = count_gas_phase_reactions + total_grain_reactions
# print(count_gas_phase_reactions, total_grain_reactions, total_reactions)
print("total grain surface reactions = ", total_grain_reactions)
print("total reactions = ", total_reactions)

threephase = 0
if isThreePhaseReaction:
    threephase = 1
else:
    threephase = 0

chemical_desorption = 0
if isChemicalDesorption:
    chemical_desorption = 1

file = open("input/input", mode="w")
file.write("species_file                  : " + species_file + "\n")
file.write("binding_energy_file           : " + binding_energy_file_out + "\n")
file.write("reaction_file                 : " + reaction_file + "\n")
file.write("abundance_file                : " + abundance_file_input + "\n")
file.write("element_file                  : " + element_file_input + "\n")
file.write("activation_energy_file        : " + activate_energy_file_input + "\n")
file.write("enthalpy_file                 : " + enthalpy_file_input + "\n")
file.write("number_of_gaseous_species     : " + str(count_gas_species) + "\n")
file.write("number_of_species             : " + str(total_species) + "\n")
file.write("number_of_elements            : " + str(count_element) + "\n")
file.write("number_of_reactions           : " + str(total_reactions) + "\n")
file.write("number_of_gas_phase_reactions : " + str(count_gas_phase_reactions) + "\n")
file.write("is_chemical_desorption        : " + str(chemical_desorption) + "\n")
file.write("is_three_phase_reaction       : " + str(threephase) + "\n")

file.close()