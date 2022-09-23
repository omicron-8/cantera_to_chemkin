import os

def convertMech(solution, mechOut, thermOut=None, tranOut=None):
    """
    Write chemkin mechanism files from a cantera solution. 
    """
    mech_file_name = os.path.join(mechOut)
    f = open(mech_file_name, 'w+', newline='\n')

    calorie    = 4184.0        # calories per kiloJoule
    atm        = 101325.0      # Pascal per atmosphere
    debye      = 3.33564e-30   # 1 debye = d coulomb-meters
    kB         = 1.380649e-23  # Boltzmann constant m2*kg/(s2*K)
    com_line   = ('!'+ "-"*72 + '\n')
    
    def removeChar(string, char_to_replace):
        for char in char_to_replace:
            string = string.replace(char, "")
        return string

    def get_arrhenius(rate_object, order):
        A = '{:<11.3E}'.format(rate_object.pre_exponential_factor*10.0 ** (3*(order-1)) )
        n = '{:<9.3f}'.format(rate_object.temperature_exponent)
        E = '{:<9.2f}'.format(rate_object.activation_energy / calorie)
        return A + n + E


    # Write elements and species
    f.write(com_line)
    f.write('ELEMENTS\n')
    for element in solution.element_names:
        f.write(element+' ')
    f.write('\nEND\n')
    f.write('SPECIES\n')
    string      = ''
    # Use space equal to longest name + 2
    max_len     = max([len(sp_name) for sp_name in solution.species_names])+2
    line_length = 0
    for sp_name in solution.species_names:
        line_length += max_len
        if line_length >= 65:
            string     += '\n'
            line_length = 0
        string += sp_name + ' '*(max_len - len(sp_name))
    f.write(string)
    f.write('\nEND\n')

    # Redirect output to thermo file
    if thermOut != None:
        thermo_file_name = os.path.join(thermOut)
        f_save = f
        f = open(thermo_file_name, 'w+', newline='\n')
    else:
        f.write(com_line)
    
    f.write('THERMO ALL\n' + '   300.000  1000.000  5000.000\n')

    # Write thermo data
    for sp_index in range(solution.n_species):
        species     = solution.species(sp_index)
        name        = str(solution.species(sp_index).name)
        nasa_coeffs = solution.species(sp_index).thermo.coeffs
        t_low   = '{:.3f}'.format(species.thermo.min_temp)
        t_max   = '{:.3f}'.format(species.thermo.max_temp)
        t_mid   = '{:.3f}'.format(species.thermo.coeffs[0])
        t_range = str(t_low) + '  ' + str(t_max) + '  ' + t_mid
        species_comp = ''
        for atom in species.composition:
            species_comp += '{:<4}'.format(atom)
            species_comp += str(int(species.composition[atom]))
        
        if solution.phase_of_matter == 'gas':
            species_phase = 'G'
        elif solution.phase_of_matter == 'liquid':
            species_phase = 'L'
        else:
            raise Exception("Unknown phase: " + solution.phase_of_matter)
        
        f.write( ('{:<18}'.format(name) + 
                  '{:<6}'.format('    ') +
                  '{:<20}'.format(species_comp) +
                  '{:<4}'.format(species_phase) +
                  '{:<31}'.format(t_range) + '1\n' ) )
        f.write( ('{:>15.8e}'.format(nasa_coeffs[1]) +
                  '{:>15.8e}'.format(nasa_coeffs[2]) +
                  '{:>15.8e}'.format(nasa_coeffs[3]) +
                  '{:>15.8e}'.format(nasa_coeffs[4]) +
                  '{:>15.8e}'.format(nasa_coeffs[5]) + '    2\n') )
        f.write( ('{:>15.8e}'.format(nasa_coeffs[6]) +
                  '{:>15.8e}'.format(nasa_coeffs[7]) +
                  '{:>15.8e}'.format(nasa_coeffs[8]) +
                  '{:>15.8e}'.format(nasa_coeffs[9]) +
                  '{:>15.8e}'.format(nasa_coeffs[10]) + '    3\n') )
        f.write( ('{:>15.8e}'.format(nasa_coeffs[11]) +
                  '{:>15.8e}'.format(nasa_coeffs[12]) +
                  '{:>15.8e}'.format(nasa_coeffs[13]) +
                  '{:>15.8e}'.format(nasa_coeffs[14]) + '                   4\n') )
        
    f.write('END\n')

    # Redirect output the main file
    if thermOut != None:
        f.close()
        f = f_save
        print("Wrote CK thermo file to '" + thermo_file_name + "'.")

    # Write reactions
    f.write(com_line)
    f.write('REACTIONS\n')
    for reac_index in range(solution.n_reactions):
        equation_string = str(solution.reaction_equation(reac_index))
        equation_string = equation_string.replace(' ', '')
        equation_object = solution.reaction(reac_index)
        equation_type = type(equation_object).__name__

        if equation_type == 'ElementaryReaction':
            order = sum(equation_object.reactants.values())
            if len(equation_object.orders) > 0:
                order = sum(equation_object.orders.values())
            arrhenius = get_arrhenius(equation_object.rate, order)
            f.write('{:<41}'.format(equation_string) + arrhenius + '\n')
            for key in equation_object.orders:
                # Explicit orders
                f.write('   FORD / '+key+' '+str(equation_object.orders[key])+' /\n')
        
        elif equation_type == 'ThreeBodyReaction':
            order = sum(equation_object.reactants.values()) + 1
            arrhenius = get_arrhenius(equation_object.rate, order)
            f.write('{:<41}'.format(equation_string) + arrhenius + '\n')
            efficiencies = str(equation_object.efficiencies)
            efficiencies = efficiencies.replace(':','/')
            efficiencies = removeChar(efficiencies,' {}\'')
            efficiencies = efficiencies.replace(',','/ ')
            if len(efficiencies) > 0:
                f.write('   ' + efficiencies + '\n')

        elif equation_type == 'FalloffReaction':
            coeff_sum = sum(equation_object.reactants.values())
            arr_high = get_arrhenius(equation_object.high_rate, coeff_sum)
            f.write('{:<41}'.format(equation_string) +arr_high + '\n')
            arr_low = get_arrhenius(equation_object.low_rate, coeff_sum + 1)
            f.write('{:<41}'.format('   LOW  /') + arr_low + ' /\n')
            j = equation_object.falloff.parameters
            # If optional Arrhenius data included:
            if equation_object.falloff.type == 'Lindemann':
                pass
            elif equation_object.falloff.type == 'Troe':
                if j[3] == 0:
                    f.write('   TROE /' +
                            '  ' + '{:.3E}'.format(j[0]) +
                            '  ' + '{:.3E}'.format(j[1]) +
                            '  ' + '{:.3E}'.format(j[2]) + ' /\n')
                else:
                    f.write('   TROE /' +
                            '  ' + '{:.3E}'.format(j[0]) +
                            '  ' + '{:.3E}'.format(j[1]) +
                            '  ' + '{:.3E}'.format(j[2]) +
                            '  ' + '{:.3E}'.format(j[3]) +' /\n')
            elif equation_object.falloff.type == 'SRI':
                f.write('   SRI /' +
                        '  ' + '{:.3E}'.format(j[0]) +
                        '  ' + '{:.3E}'.format(j[1]) +
                        '  ' + '{:.3E}'.format(j[2]) + ' /\n')
            else:
                raise Exception("Unknown falloff type: " + equation_object.falloff.type)
            efficiencies = str(equation_object.efficiencies)
            efficiencies = efficiencies.replace(':','/')
            efficiencies = removeChar(efficiencies,' {}\'')
            efficiencies = efficiencies.replace(',','/ ')
            if len(efficiencies) > 0:
                f.write('   ' + efficiencies + '\n')
        elif equation_type == 'PlogReaction':
            f.write('{:<41}'.format(equation_string) +
                    '{:>9}'.format(0) +
                    '{:>9}'.format(0) +
                    '{:>11}'.format(0) + '\n')
            order = sum(equation_object.reactants.values())
            for rate in equation_object.rates:
                pressure  = '{:.3E}'.format(rate[0] / atm)
                arrhenius = get_arrhenius(rate[1], order)
                f.write('{:<29}'.format('   PLOG / ') +
                        '{:<12}'.format(pressure) + arrhenius + ' /\n')
        else:
            raise Exception("Unknown reaction type: " + equation_type)
        
        # Dupluicate option
        if equation_object.duplicate is True:
            f.write(' DUPLICATE' +'\n')
    f.write('END')
    f.close()

    # Write transport file
    if tranOut != None:
        tran_file_name = os.path.join(tranOut)
        f = open(tran_file_name, 'w+', newline='\n')
        for sp_index in range(solution.n_species):
            species     = solution.species(sp_index)
            name        = str(solution.species(sp_index).name)
            if species.transport.geometry == 'atom':
                geom = 0
            elif species.transport.geometry == 'linear':
                geom = 1
            elif species.transport.geometry == 'nonlinear':
                geom = 2
            else:
                raise Exception("Unknown molecule geometry: " + species.transport.geometry)
            f.write('{:<17}'.format(name) +
                    '{:>3}'.format(geom) +
                    '{:>10.3f}'.format(species.transport.well_depth/kB) + # m / kB
                    '{:>10.3f}'.format(species.transport.diameter*1e10) + # Angstrom
                    '{:>10.3f}'.format(species.transport.dipole/debye)  + # debye  
                    '{:>10.3f}'.format(species.transport.polarizability*1e30) + # cubic Angstrom
                    '{:>10.3f}'.format(species.transport.rotational_relaxation) +
                    '\n')        
        f.close()
        print("Wrote CK transport file to '" + tran_file_name + "'.")
    
    print("Wrote CK mechanism file to '" + mech_file_name + "'.")
    return
