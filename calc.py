#!/usr/bin/python
#
# ========================


__all__ = ['calc_mat_props']

import os
# -------------------------------------------------------------------------------
# Import modules
# -------------------------------------------------------------------------------
from ast import arg
import csv
from xml.dom.minidom import Element
from numpy import mean, arange, digitize, array, log10
from numpy import round as rnd
from copy import deepcopy
from classes import vtk_data
import sys
from bisect import bisect
from operator import itemgetter
import time

import data_import

# -------------------------------------------------------------------------------
# Functions for calculating material data
# -------------------------------------------------------------------------------
def calc_mat_props(parts, param, vtk, seperate): #parts = List with Objekts from Elementtype bsp: [quad_test, lin_hex]
    """ Find material properties of each part """
    #print("vtk :" , vtk)
    #print("parts", parts)
    # first check that elements are all within CT volume
    _check_elements_in_CT(parts, vtk)

    # calculate the material properties
    for p in range(len(parts)):
        if parts[p].ignore != True:
            print(parts[p].ele_type)
            #print(parts[p].elindex)
            print(parts[p].numtype)
            parts[p],count_hu, cdb_indx_u50 = _assign_mat_props(parts[p], param, vtk,seperate)
    print("Anzahl parts:",len(parts))

    #print("HU_integ :",parts[0].elements[0].HU_integ)
    #group modulus values
    print(count_hu)
    part,indi,indi2 = _refine_materials(parts,param,count_hu,seperate, cdb_indx_u50)

    #indiz. from Elements under and over 50
    ini = [indi,indi2]
    #print("ini",ini)

    return part,ini

def _assign_mat_props(part, param, vtk, seperate):
    """ Find material properties of each element and add density and moduli to class part """

    # define equations from parameters file
    e_hu, rhoQCT_equ, rhoAsh_equ, modulus_equ = _define_equations(param)


    # V3
    if param['integration'] == 'E':

        # find modulus for each element

        #modulis = []
        #density = []
        # for e in part.elements:
            #modulis.append([e.integral(param['intSteps'], vtk, rhoQCT_equ, rhoAsh_equ, modulus_equ)][0])
            #density.append([e.integral(param['intSteps'], vtk, rhoQCT_equ, rhoAsh_equ, modulus_equ)][1])

        if seperate == True:

            diff = [] #decides if under 50 or not , 0 == under 50
            vol = []
            hu_count = 0
            hu_value = []
            hu_value2 = []
            #moduli_gesamt = []
            moduli2 = []
            for r in part.elements:
                #list with result from Hu integration returns HU und BMD
                result = list(r.integral_hu(vtk,rhoQCT_equ))
                #print(r.indx)
                if result[0] <= 50.:
                    diff.append(0)
                    hu_count+=1
                    hu_value.append(result[0])

                else:
                    diff.append(1)

                vol.append(r.volume)
                hu_value2.append(result[0])

            #volume under 50 over 50  ges
            R = _get_fraction(vol, diff)
            #R= 0.2
            write_R(R)
            print("R :",R)
            #print(diff)
            print("count_hu :", hu_count)

            cdb_indx_u50 = [i for i in range(len(diff)) if diff[i] == 0]# starts with zero
            print(cdb_indx_u50)


            count_2 = 0
            for k in part.elements:
                result3 = k.integral2(vtk,diff,count_2,R,cdb_indx_u50,e_hu,rhoQCT_equ,rhoAsh_equ,modulus_equ)
                moduli2.append(result3)
                count_2+=1

        #list with e, den and hu
        #like bonemat, old version
        moduli = [e.integral(vtk, rhoQCT_equ, rhoAsh_equ, modulus_equ) for e in part.elements] #moduli returns list with tupel with len elements with value of E, den, Hu integratet

        #print(moduli)
        #print(moduli2)

        #zum testen, ob  seperate fubnktioniert (anzahl der ersten EInträge)
        #hu_count = 2
        modulis = [] #oeigin moduli like bonemat
        density = []
        density2 = []
        modulis2 = [] #new moduki with new equation
        volume1 = []
        volume2 = []

        #for d in moduli:
        for d in moduli:
            modulis.append(d[0])
            density.append(d[1])

        for l in moduli2:
            modulis2.append(l[0])
            density2.append(l[1])


        for e in part.elements:
            volume2.append(e.volume)
            #print(e.volume)

        #print(modulis)
        #print(modulis2)

        print(max(modulis))
        print(max(modulis2))
        print(min(modulis))
        print(min(modulis2))


        e_u = []
        for i in range(len(diff)):
            if diff[i] == 0:

                e_u.append(modulis2[i])


        print(hu_value2)
        print(hu_value)
        print(e_u)
        print(modulis2)
        #print(max(e_u)+","+min(e_u))
        print("diff")
        print("E_new", sum(e_u)/len(e_u))



        # save
        part = _save_modulus_density(part, modulis, density, modulis2, density2) #modulis 2 = new moduli


        return part, hu_count, cdb_indx_u50

def _save_modulus_density(part, modulus, density, modulus2,density2):
    """ Save modulus data in part class """


    part.moduli = modulus
    part.density = density

    part.moduli2 = modulus2
    part.density2 = density2

    csv_moduli(modulus, density)

    return part

def csv_moduli(moduli, density):
    factor = 1000 # because density is in g/cm^3 to write in mg/cm^3 or in kg/m^3, its rho qct
    dens = [i * factor for i in density]
    mod = sorted(moduli, reverse=True)
    den = sorted(dens, reverse = True)

    #mod = moduli
    #den = density
    #h = hu

    d_max = max(den)
    d_min = min(den)

    m_max = max(mod)
    m_min = min(mod)
    m_mean = sum(mod)/len(mod)



    with open('moduli_txt.txt', 'w') as f:
        f.write("maximum :" + str(round(m_max,8)))
        f.write("\nminimum :" + str(round(m_min,8)))
        f.write("\nmean :" + str(round(m_mean,8)))
        for i in mod:
            f.write('\n')
            f.write(str(round(i,8)))

    with open('density_txt.txt', 'w') as d:
        d.write("maximum :" + str(round(d_max,8)))
        d.write("\nminimum :" + str(round(d_min,8)))
        for j in den:
            d.write('\n')
            d.write(str(round(j,8)))

def _check_elements_in_CT(parts, vtk):
    for p in parts: #p in parts = Object from rpart = part_data from fuction _get_part_data_r, parts = list with Objekt parts
        node_data = array(_get_node_data(p))#transform list in a array with shape i x y [i= kootenanzahl, y = 4 (anazhl inhalte) )
        #print(node_data)
        #print(node_data.T)
        #print(node_data.T.shape)

        i, x, y, z = node_data.T #KS ist von CT und CDB gleich, da von STL (Amira,CT) übernommen in Catia uns Ansys

        #print(i)
        #print(x)
        #print("vtkx:", vtk.x)
        #print("vtkx:", vtk.x.shape)
        #print(y)
        #print("vtky:", vtk.y)
        #print("vtky:", vtk.y.shape)
        #print(z)
        #print("vtkz:", vtk.z)
        #print("vtkz:", vtk.z.shape)
        if min(x) < min(vtk.x):
            n_indx = i[x.tolist().index(min(x))]
            raise ValueError("Error: Node " + repr(int(n_indx)) + " has an x-coordinate of: " + repr(
                min(x)) + " which is outside the CT volume\n"
                          "Dataset minimum x value is: " + repr(min(vtk.x)))
        elif min(y) < min(vtk.y):
            n_indx = i[y.tolist().index(min(y))]
            raise ValueError("Error: Node " + repr(int(n_indx)) + " has a y-coordinate of: " + repr(
                min(y)) + " which is outside the CT volume\n"
                          "Dataset minimum y value is: " + repr(min(vtk.y)))
        elif min(z) < min(vtk.z):
            n_indx = i[z.tolist().index(min(z))]
            raise ValueError("Error: Node " + repr(int(n_indx)) + " has a z-coordinate of: " + repr(
                min(z)) + " which is outside the CT volume\n"
                          "Dataset minimum z value is: " + repr(min(vtk.z)))
        elif max(x) > max(vtk.x):
            n_indx = i[x.tolist().index(max(x))]
            raise ValueError("Error: Node " + repr(int(n_indx)) + " has an x-coordinate of: " + repr(
                max(x)) + " which is outside the CT volume\n"
                          "Dataset maximum x value is: " + repr(max(vtk.x)))
        elif max(y) > max(vtk.y):
            n_indx = i[y.tolist().index(max(y))]
            raise ValueError("Error: Node " + repr(int(n_indx)) + " has a y-coordinate of: " + repr(
                max(y)) + " which is outside the CT volume\n"
                          "Dataset maximum y value is: " + repr(max(vtk.y)))
        elif max(z) > max(vtk.z):
            n_indx = i[z.tolist().index(max(z))]
            raise ValueError("Error: Node " + repr(int(n_indx)) + " has a z-coordinate of: " + repr(
                max(z)) + " which is outside the CT volume\n"
                          "Dataset maximum z value is: " + repr(max(vtk.z)))

def _get_node_data(part):
    """ Identifies node data using part class """

    # join node index with co-ordinate data
    nodes = {}

    for e in part.elements:
        for n in range(len(e.nodes)):
            nodes[e.nodes[n]] = e.pts[n]



    # create array of data
    node_data = []

    for n in list(nodes.keys()):
        if len(nodes[n]) == 1:
            nodes[n].append(0.0)
            nodes[n].append(0.0)
        elif len(nodes[n]) == 2:
            nodes[n].append(0.0)

        node_data.append([int(n), nodes[n][0], nodes[n][1], nodes[n][2]])

    #print(sorted(node_data))
    #print(len(node_data))
    #print(type(node_data))

    return sorted(node_data) #list with i list einträgen(i = Elementanzahl) [[1,2,3,4], [....],[....]] node index and coordinates from node

def _get_fraction(vol,dif):
    vol1 = []
    vol2 = []
    for i in range(len(dif)):
        #print(i)
        if dif[i] == 0:
            vol1.append(vol[i])
        else:
            vol2.append(vol[i])

    v1 = sum(vol1) #under 50
    v2 = sum(vol2)#over 50
    v3 = sum(vol)#all volume

    with open ("V.txt", "w") as v:
        v.write("under 50 :" + str(v1) + "\n")
        v.write("over 50 :" + str(v2) + "\n")
        v.write("all:" + str(v3))

    print("all volume :",v3)
    print("under 50 :",v1)

    R = 1 - (v2/v3)
    return R


# -------------------------------------------------------------------------------
# Functions for calculating modulus values
# -------------------------------------------------------------------------------
def _define_equations(param): # wie funktionen in parametern/variablen speichern ???, dauert lamge, da für jedes Element abgefragt wird
    """ From parameters inputs define equation lambda functions """

    # define calibration calculation
    def rhoQCT(HU): #returns in /cm^3
        #hu = []
        #hu.append(HU)
        #print("HU_Wert :", hu)
        #huu = sorted(HU, reverse=True)
        #print(huu)
        #if param['HUValues'].lower() == 'yes':
        #    if HU <= 300:
        #        res = param['rhoQCTar1'] + (param['rhoQCTbr1'] * HU)
        #        if res < 0:
        #            return 0.000001
        #        else:
        #            return res
        #    elif HU >= 301 and HU <= 700:
        #        res = param['rhoQCTar2'] + (param['rhoQCTbr2'] * HU)
        #        if res < 0:
        #            return 0.000001
        #        else:
        #            return res
        #    elif HU >= 700 and HU <= 1000:
        #        res = param['rhoQCTar3'] + (param['rhoQCTbr3'] * HU)
        #        if res < 0:
        #            return 0.000001
        #        else:
        #            return res


        #if param['HUValues'].lower() == 'no': # geht für jedes element durch , braucht lange zeit
            #print(param['HUValues'].lower ())
        os.chdir(sys.path[0]) #öffnet es für jedes element, dauert super lange, in die andere Schleife, wo aufgerufen wird machen !!!!!
        with open("HU_phantoms.txt", "r") as f:
            objekt = f.read()
            if objekt == "None":

                res = param['rhoQCTa'] + (param['rhoQCTb'] * HU)
                if res < 0:
                    return 0.000001
                else:
                    return res
            else:

                coef = objekt.split(";")
                coef[1] = float(coef[1])
                coef[0] = float(coef[0])

                result = coef[0] + (coef[1] * HU)
                #print("calibration: on ")
                if result < 0:
                #print("HU smaleer than 0, density, too ")
                    return 0.000001

                else:
                    return result #in g/cm^3


    def hu_e(q,R):
        """equation follows from fat """
        E = (10**(-8.54 + (4.05 * (log10(400 * (q / (q + 1.4 * R))))))*1000)

        if E < 0:
            return 0.000001
        elif E == None:
            return 0.000001
        else:
            return E


    def rhoAsh(Q):
        res = param['rhoAsha11'] + (param['rhoAshb11'] * Q)
        if res < 0:
            return 0.000001
        else:
            return res

    # define modulus calculation in MPA
    if param['numEparam'] == 'single':
        def modulus(Ash):
            res = param['Eb11'] * (Ash ** param['Ec11'])
            if res < 0:
                return 0.000001
            else:
                return res
            # IR
    elif param['numEparam'] == 'double':
        #in that case, Ea11 is 0 , Eb11/ Ec11 default value for all
        def modulus(Ash):
            if Ash < param['Ethresh1']:
                res = param['Eb1'] * (Ash ** param['Ec1'])
                if res < 0:
                    res = 0.000001

            elif Ash >= param['Ethresh1']:
                res = param['Eb11'] * (Ash ** param['Ec11'])
                if res < 0:
                    res = 0.000001
            return res



    #elif param['numEparam'] == 'triple':
        #def modulus(Ash):
        #    if Ash < param['Ethresh1']:
        #        res = param['Ea1'] + (param['Eb1'] * (Ash ** param['Ec1']))
        #    elif (Ash >= param['Ethresh1']) & (Ash <= param['Ethresh2']):
        #        res = param['Ea2'] + (param['Eb2'] * (Ash ** param['Ec2']))
        #    elif Ash > param['Ethresh2']:
        #        res = param['Ea3'] + (param['Eb3'] * (Ash ** param['Ec3']))
        #    else:
        #        raise ValueError("Error: modulus for density value: " + repr(Ash) + " cannot be calculated")
        #    if res < 0:
        #        return 0.000001
        #    else:
        #        return res
    else:
        IOError("Error: " + param['numEparam'] + " is not a valid input for numCTparam.  Must be 'single' or 'triple'")
    return hu_e, rhoQCT, rhoAsh, modulus


# -------------------------------------------------------------------------------
# Functions for grouping modulus values
# -------------------------------------------------------------------------------
def _refine_materials(parts, param,count_hu,seperate, cdb_indx):
    """ Group the materials into bins separated by the gapValue parameter """


    #moduli = _get_all_modulus_values(parts) # einzelne moduli der Elemente, kein gab value
    #print("parts_moduli: ", parts[0].moduli, parts[1].moduli)
    # limit the moduli values for each part
    for p in range(len(parts)):
        if parts[p].ignore != True:
            parts[p].moduli, parts[p].moduli2,parts[p].sep_mat,indi,indi2  = _limit_num_materials(parts[p].moduli,parts[p].moduli2,parts[p].density2,
                                                                            param['gapValue'],
                                                                            param['minVal'],
                                                                            param['groupingDensity'], count_hu,seperate,cdb_indx)


    #print(indi,indi2)
    return parts,indi,indi2

def _get_all_modulus_values(parts):
    """ Create list of all modulus values for all parts in model """

    moduli = []
    for p in parts:
        if p.ignore != True:
            moduli.extend(p.moduli)

    return moduli

def _get_mod_intervals(moduli, max_materials):
    """ Find the refined moduli values """

    min_mod = float(min(moduli))
    max_mod = float(max(moduli))
    if len(set(moduli)) < max_materials:
        mod_interval = 0.0
    else:
        mod_interval = (max_mod - min_mod) / float(max_materials)

    return mod_interval

def _limit_num_materials(moduli, moduli2, density2, gapValue, minVal, groupingDensity,count_hu,seperate,cdb_indx):
    """ Groups the moduli into bins and calculates the max/mean for each """
    if gapValue == 0:
        return moduli
    else:
        #warn user if there are a lot of bins to calculate
        binLength = (max(moduli) - min(moduli)) / gapValue #
        if binLength > 10000:
            print('\n    WARNING:')
            print('    You have specified a very small gap size relative to your maximum modulus')
            print(('    So your modulus "bin" size is: ' + repr(round(binLength))))
            print('    This will take so long to calculate it may crash your computer ***')
            answ = eval(input('    Would you like to continue (y/n)'))
            if (answ == 'n') | (answ == 'N') | (answ == 'no') | (answ == 'No') | (answ == 'NO'):
                sys.exit()
            else:
                print('\n')

        # calculate the modulus values
        #bins = arange(max(moduli), minVal - gapValue, -gapValue).tolist()
        #indices_h, sorted_moduli_h = list(zip(*sorted(enumerate(hu), key=itemgetter(1))))
        indices, sorted_moduli = list(zip(*sorted(enumerate(moduli), key=itemgetter(1))))

        print(count_hu)

        """
        enumerate: adding iteration to list/ tuple  from [a,b,c] --> [(0,a),(1,b),(2,c)]
        sorted/itemgetter(): sot list with key, itemgetter with the value you want to sort from tuple
        zip with star: create iteration (list or tuple...) from two objects, with stare return multiple results
        """
        #element indices and sorted list of all moduli, not grouped (read like red the cdb, for element indx)
        #print("indices :", indices)
        #print("sorted_ moduli :", sorted_moduli)



        # work way through list adding modified modulus values
        new_moduli = [minVal] * len(moduli)
        while len(sorted_moduli) > 0:
            b = bisect(sorted_moduli, sorted_moduli[-1] - gapValue)  # position of next value outside gab

            if groupingDensity == 'max':
                val = sorted_moduli[-1]
            elif groupingDensity == 'mean':
                val = sum(sorted_moduli[b:]) / len(sorted_moduli[b:])
                # print(val)
            else:
                raise IOError("Error: groupingDensity should be 'max' or 'mean' in parameters file")
            if val < minVal:
                val = minVal
            for i in indices[b:]:
                new_moduli[i] = val
            sorted_moduli = sorted_moduli[:b]
            indices = indices[:b]


        if seperate == True:


            indices5 = [i for i in list(range(len(moduli2))) if i not in cdb_indx]
            moduli5 = [moduli2[i] for i in indices5]
            low_moduli5 = [moduli2[i] for i in cdb_indx]
            #print(moduli2)
            #print(moduli5)
            #print(low_moduli5)

            indices55, sorted_moduli5 = list(zip(*sorted(enumerate(moduli5), key=itemgetter(1))))
            new_moduli5 = [minVal] * len(moduli5)
            while len(sorted_moduli5) > 0:
                b = bisect(sorted_moduli5, sorted_moduli5[-1] - gapValue)  # position of next value outside gab

                if groupingDensity == 'max':
                    val = sorted_moduli[-1]
                elif groupingDensity == 'mean':
                    val = sum(sorted_moduli5[b:]) / len(sorted_moduli5[b:])
                    # print(val)
                else:
                    raise IOError("Error: groupingDensity should be 'max' or 'mean' in parameters file")
                if val < minVal:
                    val = minVal
                for i in indices55[b:]:
                    new_moduli5[i] = val
                sorted_moduli5 = sorted_moduli5[:b]
                indices55 = indices55[:b]

            value_limit5 = sum(low_moduli5)/len(low_moduli5)
            #print(low_moduli5)
            #print(len(low_moduli5))
            #print(value_limit5)

            for i in range(len(cdb_indx)):
                new_moduli5.insert(cdb_indx[i], value_limit5)
            #print((sorted(new_moduli5)))




            indices2, sorted_moduli2 = list(zip(*sorted(enumerate(moduli2), key=itemgetter(1))))
            #print(indices2)
            #print(sorted_moduli2)

            sorted_dens = sorted(density2)
            new_moduli2 = [minVal] * len(moduli)
            indices3 = indices2[:]
            sorted_moduli3 = sorted_moduli2[:]
            while len(sorted_moduli2) > 0:
                c = bisect(sorted_moduli2, sorted_moduli2[-1] - gapValue)  # position of next value outside gab

                if groupingDensity == 'max':
                    val2 = sorted_moduli[-1]
                elif groupingDensity == 'mean':
                    val2 = sum(sorted_moduli2[c:]) / len(sorted_moduli2[c:])
                    # print(val)
                else:
                    raise IOError("Error: groupingDensity should be 'max' or 'mean' in parameters file")
                if val2 < minVal:
                    val2 = minVal
                for i in indices2[c:]:
                    new_moduli2[i] = val2
                sorted_moduli2 = sorted_moduli2[:c]
                indices2 = indices2[:c]


            indi = list(indices3[:count_hu])
            Indi = [i+1 for i in indi]
            indi2 = list(indices3[count_hu:])
            Indi2 = [i+1 for i in indi2]
            value_limit = (sum(sorted_moduli3[:count_hu])) / count_hu
            value_limit_dens = (sum(sorted_dens[:count_hu])/count_hu)*0.000000001
            #print(value_limit_dens)

            val_lim = [value_limit,value_limit_dens]
            val_lim5 = [value_limit5, value_limit_dens]
            print("value_limit",val_lim5)

            for i in indi:
                new_moduli2[i] = value_limit

        #print((sorted(new_moduli)))
        #print((sorted(new_moduli2)))
        #print("indi", Indi)
        #print(len(indi))
        #print("indi2", Indi2)
        #print(len(Indi2))

        return new_moduli,new_moduli5,val_lim5,cdb_indx,indices5 #indi2 = over 50

        #return new_moduli,new_moduli2,val_lim,Indi,Indi2 #indi2 = over 50

def write(new_moduli):
    new_moduli_sortet = sorted(new_moduli, reverse=True)
    set_moduli = set(new_moduli_sortet)
    set_moduli_list = sorted(list(set_moduli), reverse=True)

    count = 1
    p = []
    s = new_moduli_sortet[:]
    s.append(0)
    #print(s)
    for i in range(len(s) - 1):
        if s[i] == s[i + 1]:
            count += 1
        elif s[i] != s[i + 1]:
            p.append(count)
            count = 1
    #print(p)
    #print(new_moduli_sortet)

    with open('new_moduli_txt.txt', 'w') as w:
        w.write("max :" + str(max(new_moduli_sortet)))
        w.write("\nmin :" + str(min(new_moduli_sortet)))
        w.write("\nmean :" + str(sum(new_moduli_sortet) / len(new_moduli_sortet)))
        w.write("\nlen :" + str(len(new_moduli_sortet)))
        for l in new_moduli_sortet:
            w.write("\n")
            w.write(str(round(l)))

    with open("set_moduli_txt.txt", 'w') as s:
        s.write("max :" + str(max(set_moduli_list)))
        s.write("\nmin :" + str(min(set_moduli_list)))
        s.write("\nmean :" + str(sum(set_moduli_list) / len(set_moduli_list)))
        s.write("\nlen: " + str(len(set_moduli_list)))
        for l in set_moduli_list:
            s.write("\n")
            s.write(str(l))

    with open("len_elements_pro_bin.txt", "w") as j:
        for i in p:
            j.write(str(i))
            j.write("\n")


def write_R(R):

    with open("R.txt", "w") as i:
        i.write(str(R))




#############################################################################################################


#vtk = data_import.import_ct_data(r"C:\Users\8iries\Desktop\Praktikum\Programm\series-01")
#param = data_import.import_parameters_xml(r"C:\Users\8iries\Desktop\Praktikum\Programm\params.html")
#parts = data_import.import_mesh(r"C:\Users\8iries\Desktop\Praktikum\Programm\mesh_P5.cdb", param)
#calc_mat_props(parts, param, vtk)


#############################################################################################################





#vtk = data_import.import_ct_data(r"C:\Users\8iries\Desktop\Praktikum\Programm\CT_dataset")
#param = data_import.import_parameters_xml(r"C:\Users\8iries\Desktop\Praktikum\Programm\params.xml")
#parts = data_import.import_mesh(r"C:\Users\8iries\Desktop\Praktikum\Programm\ct2.cdb", param)

#calc_mat_props(parts, param, vtk)

