#!/usr/bin/python
#
# py_bonemat_ansys - data output
# ===============================
#
# Created by Abdul Hanan

__all__ = ['output_ansys_inp']

import os
import sys
# -------------------------------------------------------------------------------
# Import modules
# -------------------------------------------------------------------------------
import re
from general import _read_text
from data_import import _get_lines, _flat_List, _get_marks
from datetime import datetime
import shutil

# -------------------------------------------------------------------------------
# Functions to output new ansys file
# -------------------------------------------------------------------------------

def output_ansys_inp(parts, mesh_file, param, ansys_v,seperate,indi):
    #print(type(parts))
    #print(len(parts))
    #print(parts[0].elements)

    #normal integration without seperation V3
    indexo, densities, modindex, modules = modulus_values(parts)

    # Reading ANSYS File
    files = _read_text(mesh_file)

    # Getting all Elements
    oelements = get_allElements(files, indexo)

    # Getting Nodes
    onodes, nend = get_nodes(files)

    # Getting heads and marks
    head1, head2, Nhead, head3, mark1, mark2, eend, end = get_texts(files, ansys_v)

    # Getting Element Types for Block
    etype = get_etype(parts)

    #print(oelements)
    #print(type(oelements))
    # Forming material block ansys19
    materials = ansys19(param, modindex, modules, densities)


    if seperate == True:


        indexo2, densities2, modindex2, modules2, seperate_mat,sep_mat_one = modulus_values2(parts, seperate, indi)

        # get elements
        oelements_2, oelements_3, oelements_4,nodeinx_sep2,nodeinx_sep3 = get_allElements2(files, indexo2, indi, sep_mat_one, seperate)

        #get nodes
        onodes2, onodes3, nend = get_nodes2(files, nodeinx_sep2,nodeinx_sep3)

        #materials
        materials2, materials3, materials4 = ansys192(param, modindex2, modules2, densities2, seperate_mat)

        # Getting Element Types for Block
        etype = get_etype(parts)

        # Getting heads and marks
        head1, head2, Nhead, head3, mark1, mark2, eend, end = get_texts(files, ansys_v)

        """Head überarbeiten, len Elements"""
        Nhead_2,Ehead_2,Nhead_3,Ehead_3 = Nhead_sep(Nhead, head3, indi,nodeinx_sep2,nodeinx_sep3)

        #under 50
        complete2 = head1 + head2 + etype + Nhead_2 + mark1 + onodes2 + nend + Ehead_2 + mark2 + oelements_2 + eend + materials2 + "\n" + end
        #over 50
        complete3 = head1 + head2 + etype + Nhead_3 + mark1 + onodes3 + nend + Ehead_3 + mark2 + oelements_3 + eend + materials3 + "\n" + end

        complete4 = head1 + head2 + etype +Nhead + mark1 + onodes + nend + head3 + mark2 + \
                    oelements_3 + oelements_4 +eend + materials3 + materials4 +"\n" + end



    #head hinzu gfügt = Nodeblock
    complete = head1 + head2 + etype +Nhead + mark1 + onodes + nend + head3 + mark2 + oelements + eend + materials + "\n" + end


    cdb = [complete,complete2,complete3,complete4]

    write_cdb(cdb)

    #objekt = open ("pathnew.txt", 'r') # get path from txt from GUI input
    #filepath = objekt.read()
    #print(filepath)

    #objekt_2 = open("name.txt", 'r') # get path from txt from GUI input
    #name = objekt_2.read()
    #name = [name + ".cdb"]
    #print(name)

    #with open(name[0], 'w') as oupf: # verschieben der Daten in angebenen Ordner
    #    oupf.write(complete)

    #shutil.move(name[0], filepath + "/" + name[0]) #for run with test run

def write_cdb(cdb):
    with open ("name.txt") as i:
        n = i.read()

    with open ("pathnew.txt") as f:
        filepath = f.read()

    name_0 = [[str(n) +".cdb"],[str(n)+ "_HU_u_50" + ".cdb"],[str(n)+ "_HU_o_50" + ".cdb"],[str(n)+"_HU_sep_one"+".cdb"]]
    name = [filepath+"/"+i[0] for i in name_0]

    print(name)

    #os.chdir(filepath)

    #for i in range(len(cdb)):
    with open(name[0],'w') as out:
        out.write(cdb[0])
    with open(name[1],'w') as out:
        out.write(cdb[1])
    with open(name[2],'w') as out:
        out.write(cdb[2])
    with open(name[3],'w') as out:
        out.write(cdb[3])


    #with open(name[0], 'w') as oupf:  # verschieben der Daten in angebenen Ordner
    #    oupf.write(complete)

    #shutil.move(name[0], filepath + "/" + name[0])  # for run with test run

def Nhead_sep(nhead,head3,indi,node_inx2,node_inx3):
    #print(nhead)
    #print(etype)
    #print(onodes)
    #print(olements)
    #print(nend)

    #print(len(node_inx2))
    #print(len(node_inx3))
    #print(len(indi[0]))
    #print(len(indi[1]))

    #Nodes
    n2 = "      "+str(len(node_inx2))
    #print(len(node_inx2))
    words = [i.strip() for i in nhead.split(",")]
    words[3] = n2
    words[4] = n2
    nhead_2 = ",".join(words)+"\n"

    #Elements
    N2 = "     "+str(len(indi[0]))
    words_e = [i.strip() for i in head3.split(",")]
    words_e[3] = N2
    words_e[4] = N2
    e_Head3_2=",".join(words_e)+"\n"

    # Nodes
    n3 = "      " + str(len(node_inx3))
    #print(len(node_inx3))
    words3 = [i.strip() for i in nhead.split(",")]
    words3[3] = n3
    words3[4] = n3
    nhead_3 = ",".join(words3) + "\n"

    # Elements
    N3 = "     " + str(len(indi[1]))
    words_e3 = [i.strip() for i in head3.split(",")]
    words_e3[3] = N3
    words_e3[4] = N3
    e_Head3_3 = ",".join(words_e3) + "\n"

    return nhead_2,e_Head3_2,nhead_3,e_Head3_3

def get_nodes2(files,sep_n2,sep_n3):

    file = files[:-4]
    nodes, nend = _get_onodes(file)
    # print(nodes)
    nodes = [n.replace("+00", "+0") for n in nodes]
    nodes = [n.replace("-00", "-0") for n in nodes]
    node = [n[1:9] + n[10:18] + n[19:27] + "  " + n[27:47] + "  " + n[47:67] + "  " + n[67:] + "\n" for n in nodes]

    nodes_final2 = [node[i - 1] for i in sep_n2]
    onodes_sep2 = ''.join(e for e in nodes_final2)

    nodes_final3 = [node[k - 1] for k in sep_n3]
    onodes_sep3 = ''.join(l for l in nodes_final3)

    nend = nend + '-1\n'


    return onodes_sep2,onodes_sep3,nend

def get_nodes(files):
    # Getting Nodes
    file = files[:-4]
    nodes, nend = _get_onodes(file)
    #print(nodes)
    nodes = [n.replace("+00", "+0") for n in nodes]
    nodes = [n.replace("-00", "-0") for n in nodes]
    node = [n[1:9] + n[10:18] + n[19:27] + "  " + n[27:47] + "  " + n[47:67] + "  " + n[67:] + "\n" for n in nodes]
    onodes = ''.join(e for e in node)
    nend = nend + '-1\n'



    #print(onodes, nend)
    return onodes, nend

def get_texts(files, ansys_v):
    # Heading For ANSYS

    datetimes = str(datetime.now())

    #head1 = "!! Generated by PyMat " + datetimes + "\n\n/PREP7\n"
    head1 = "/COM,ANSYS RELEASE " + str(ansys_v) + "       " + datetimes + "\n/PREP7" + "\n/NOPR\n"

    # Heading for Nodes
    pattern = "ENDIF\n(.*?)\n"

    head2 = re.findall(pattern, files)
    head2 = ''.join(head2) + "\n"

    #pattern2 = "ET,        1,187\n(.*?)\n" #for NBLOCK....
    #head = re.findall(pattern2, files)
    #head = head[0] + "\n"

    Nhead = _get_Nblock()

    # Get node mark
    # x = _get_marks (files)
    # mark1 = "("+x[0]+")\n"
    mark1 = "(3i8,6e22.13)\n"

    # Get elemeng mark
    # mark2 = "("+x[1]+")\n"
    mark2 = "(19i8)\n"

    # Element ending statement
    eend = "      -1\n"

    # Heading for Elements
    pattern2 = "-1,\n(.*?)\n"
    head3 = re.findall(pattern2, files) # files = cdb lines
    head3 = ''.join(head3) + "\n"

    end = "/GO\nFINISH"

    return head1, head2, Nhead, head3, mark1, mark2, eend, end

def _get_Nblock():
    #"""get the NBLOCK head """
    os.chdir(sys.path[0])
    with open("pathcdb.txt", "r") as f:
        pathcdb = f.read()
    with open(pathcdb, "r") as f:
        for row in f:
            if "NBLOCK" in row:
                nline = row
                break
    return nline

def ansys192(param, modindex2, modules2, densities2,mat_sep):

    materials_sep2 = ""
    materials_sep3 = ""

    #for HU under 50
    dens = [mat_sep[1]]
    modul = [mat_sep[0]]

    for v in range(len(list(modul))):
        common_line = "MPTEMP,R5.0, 1, 1,  0.00000000    ,"
        modulus = "MPDATA,R5.0, 1,EX," + ",  1" + ", 1, " + str(modul[v]) + "    ,"
        poison = "MPDATA,R5.0, 1,NUXY," + ",  1" + ", 1, " + str("{0:.8f}".format(param)) + "    ,"
        density = "MPDATA,R5.0, 1,DENS," + ",  1" + ", 1, " + str(dens[v]) + "    ,"
        materials_sep2 += common_line + "\n" + modulus + "\n" + common_line + "\n" + poison + "\n" + common_line + "\n" + density + "\n"

    #for HU over 50
    for v in range(len(modules2)):
        common_line = "MPTEMP,R5.0, 1, 1,  0.00000000    ,"
        modulus = "MPDATA,R5.0, 1,EX," + modindex2[v] + ", 1, " + str(modules2[v]) + "    ,"
        poison = "MPDATA,R5.0, 1,NUXY," + modindex2[v] + ", 1, " + str("{0:.8f}".format(param)) + "    ,"
        density = "MPDATA,R5.0, 1,DENS," + modindex2[v] + ", 1, " + str(densities2[v]) + "    ,"
        materials_sep3 += common_line + "\n" + modulus + "\n" + common_line + "\n" + poison + "\n" + common_line + "\n" + density + "\n"

    #for under 50 in cdb booth
    modindex4 = str(len(modules2)+1).rjust(4)
    common_line2 = "MPTEMP,R5.0, 1, 1,  0.00000000    ,"
    modulus2 = "MPDATA,R5.0, 1,EX," + modindex4 + ", 1, " + str(modul[0]) + "    ,"
    poison2 = "MPDATA,R5.0, 1,NUXY," + modindex4+ ", 1, " + str("{0:.8f}".format(param)) + "    ,"
    density2 = "MPDATA,R5.0, 1,DENS," + modindex4+ ", 1, " + str(dens[0]) + "    ,"
    materials_sep4 = common_line2 + "\n" + modulus2 + "\n" + common_line2 + "\n" + poison2 + "\n" + common_line2 + "\n" + density2 + "\n"




    return materials_sep2, materials_sep3,materials_sep4

def ansys19(param, modindex, modules, densities):
    # Materials for Ansys 2019
    materials = ""
    for v in range(len(modules)):
        common_line = "MPTEMP,R5.0, 1, 1,  0.00000000    ,"
        modulus = "MPDATA,R5.0, 1,EX," + modindex[v] + ", 1, " + str(modules[v]) + "    ,"
        poison = "MPDATA,R5.0, 1,NUXY," + modindex[v] + ", 1, " + str("{0:.8f}".format(param)) + "    ,"
        density = "MPDATA,R5.0, 1,DENS," + modindex[v] + ", 1, " + str(densities[v]) + "    ,"
        materials += common_line + "\n" + modulus + "\n" + common_line + "\n" + poison + "\n" + common_line + "\n" + density + "\n"

    return materials

def get_etype(parts):
    typ = []
    for num in range(len(parts)):
        types = parts[num].numtype #types = dic. {1 : 187}
        typ.append(types)
        #print(typ)
        etype = ""
    for i in typ:
        for k, v in i.items():
            etype += 'ET,' + str(k) + ',' + str(v) + '\n'
            # etype = ET, 1, 178
            #print("etype :", etype)
    return etype

def modulus_values(parts):

    # Limited Modulus Values Calculation and arrangement
    set_data = get_arranged_modulus(parts)

    module = list(set_data.values())
    modules = [round(c, 8) for c in module] # rounded modules from set in list

    # Get modulus indexes for MPBLOCK
    modindex = get_mod_index(set_data)

    # Modulus Index PLacement
    # Fliping Keys to Values
    indexes = {y: x for x, y in set_data.items()}

    t = []
    xy = [0] * len(parts)

    for p in range(len(parts)): #anzahl elements


        for i in parts[p].moduli: #i = value E -Modul von Element, gruppiert....
            if i in list(indexes.keys()): # list indexes.keys = value of E, Anzahl Materialgruppen
                t.append(indexes[i])
        xy[p] = t
        t = []

    indexo = []
    for s in xy:
        for t in s:
            indexo.append(t.rjust(8))

    # Limitng Denisty Values according to Modulus Values and Density coversion

    xny = [i.strip() for i in indexo]
    dmod = [int(r) for r in xny]
    density = _get_all_density_values(parts)
    zipped = zip(dmod, density)
    we = {}
    for t in zipped:
        key, val = t
        we.setdefault(key, []).append(val)
    vals = []
    for names, values in we.items():
        vals.append("{avg}".format(avg=sum(values) / len(values)))
    zipper = zip(we.keys(), vals)
    valss = dict(zipper)
    sd = sorted(valss.items())
    densit = []
    for k, v in sd:
        densit.append(v)
    ints = [round(float(s), 8) for s in densit] #
    #intes = [(s * 0.5084) * 0.000000001 + 0.000000001 for s in ints] # what correlation ????????????
    intes = [round(s * 0.000000001, 16) for s in ints] #density in t/mm^3 for ansys import to match with MPA for E , * e-9

    densities = [str(e) for e in intes]


    return indexo, densities, modindex, modules

def modulus_values2(parts, seperate, indi):


    modul = parts[0].sep_mat[0]
    densi = parts[0].sep_mat[1]
    seperate_mat = [modul, densi]
    all_moduli2 = []

    for p in parts:
        if p.ignore != True:
            all_moduli2.extend(p.moduli2)
        mat_values2 = sorted(list(set(all_moduli2)), reverse=True)
        print(mat_values2)
        mat_values2.remove(seperate_mat[0])
        #mat_values2.pop() #last remove because its in the HU under 50 cdb
        set_data = dict([(repr(m + 1), mat_values2[m]) for m in range(len(mat_values2))])

    # Limited Modulus Values Calculation and arrangement
    module = list(set_data.values())
    modules = [round(c, 8) for c in module]  # rounded modules from set in list
    sep_mat_one = str(len(module)+1).rjust(8)
    print(modules)
    print(sep_mat_one)

    # Get modulus indexes for MPBLOCK
    modindex = get_mod_index(set_data)
    print("modindex :",modindex)

    # Modulus Index PLacement
    # Fliping Keys to Values
    indexes = {y: x for x, y in set_data.items()}

    t = []
    xy = [0] * len(parts)

    for p in range(len(parts)):  # anzahl elements
        for i in parts[p].moduli2:  # i = value E -Modul von Element, gruppiert....
            if i in list(indexes.keys()):  # list indexes.keys = value of E, Anzahl Materialgruppen
                t.append(indexes[i])
        xy[p] = t
        t = []

    indexo = []
    for s in xy:
        for t in s:
            indexo.append(t.rjust(8))
    #print(indexo)
    # Limitng Denisty Values according to Modulus Values and Density coversion

    xny = [i.strip() for i in indexo]
    dmod = [int(r) for r in xny]
    density = []
    for p in parts:
        if p.ignore != True:
            density.extend(p.density)
    zipped = zip(dmod, density)
    we = {}
    for t in zipped:
        key, val = t
        we.setdefault(key, []).append(val)
    vals = []
    for names, values in we.items():
        vals.append("{avg}".format(avg=sum(values) / len(values)))
    zipper = zip(we.keys(), vals)
    valss = dict(zipper)
    sd = sorted(valss.items())
    densit = []
    for k, v in sd:
        densit.append(v)
    ints = [round(float(s), 8) for s in densit]  #
    # intes = [(s * 0.5084) * 0.000000001 + 0.000000001 for s in ints] # what correlation ????????????
    intes = [round(s * 0.000000001, 16) for s in
             ints]  # density in t/mm^3 for ansys import to match with MPA for E , * e-9

    densities = [str(e) for e in intes]

    return indexo, densities, modindex, modules, seperate_mat, sep_mat_one

def get_allElements2(files, indexo2,indi,sep_mat_one,seperate):


    elements2,elements3,elements4, node_sep2, node_sep3 = _get_oelem2(files, indexo2, indi,sep_mat_one, seperate)

    elements = [elements2,elements3,elements4] # elements 2 == HU under 50

    typelements =[_get_typesele(i) for i in elements]

    list2 = [0] * len(typelements[0])
    list3 = [0] * len(typelements[1])
    list4 = list2[:]

    for l in typelements[0]:
        list2[typelements[0].index(l)] = _get_arrange(l)
    values2 = []
    for i in list2:
        values2.append(_get_ofelem(i))
    oele2_sep = ''.join(e for e in values2)

    for l in typelements[1]:
        list3[typelements[1].index(l)] = _get_arrange(l)
    values3 = []
    for i in list3:
        values3.append(_get_ofelem(i))
    oele3_sep = ''.join(e for e in values3)

    for l in typelements[2]:
        list4[typelements[2].index(l)] = _get_arrange(l)
    values4 = []
    for i in list4:
        values4.append(_get_ofelem(i))
    oele4_sep = ''.join(e for e in values4)


    return oele2_sep, oele3_sep, oele4_sep, node_sep2,node_sep3

def get_allElements(files, indexo):

    elements = _get_oelem(files, indexo)
    typelements = _get_typesele(elements)
    #print(typelements)
    list3 = [0] * len(typelements)
    for l in typelements:
        list3[typelements.index(l)] = _get_arrange(l)
    #print(list3)
    values = []
    for i in list3:
        values.append(_get_ofelem(i))
    oelements = ''.join(e for e in values)

    return oelements

def _get_all_density_values(parts):
    """ Create list of all density values for all parts in model """

    density = []
    for p in parts:
        if p.ignore != True:
            density.extend(p.density)
    #print(sorted(density, reverse=True))



    return density

def _get_typesele(elements):
    tens = []
    eights = []
    sixes = []
    fours = []
    for i in elements:
        if len(i) == 168:
            tens.append(i)
        elif len(i) == 152:
            eights.append(i)
        elif len(i) == 136:
            sixes.append(i)
        elif len(i) == 120:
            fours.append(i)

    list1 = [tens] + [eights] + [sixes] + [fours]
    list2 = [x for x in list1 if x != []]

    return list2

def get_mod_index(set_data):
    modsindex = list(set_data.keys())
    #print(modsindex)
    modsindexes = [str(b) for b in modsindex] # str for writing txt ?
    ind = ["     " + c for c in modsindexes]
    return ind

def get_arranged_modulus(parts):
    #print("parts :" ,parts)
    #print(parts[0].moduli)
    all_moduli = []

    for p in parts:
        if p.ignore != True:
            #print(p.ignore)
            all_moduli.extend(p.moduli)# hängt part moduli an
        set_data = _create_material_elesets(all_moduli)

        #print("set Data :", set_data)
        #print(set_data)
    return set_data

def _create_material_elesets(moduli): #Moduli List eder einzelnen E Modul für jedes Element
    mat_values = sorted(list(set(moduli)), reverse=True)


    #print("mat_value",mat_values)
    mat_sets = dict([(repr(m + 1), mat_values[m]) for m in range(len(mat_values))])

    #print("mat_set",mat_sets)
    return mat_sets

def _get_arrange(oelement):
    d = [] #oelement = liste mit strings , strings = jewei,s eine Zeile aus Elementsblock
    for t in oelement:
        d.append(int(t[8:16].strip()))
    #print(d)
    arranged = [x for _, x in sorted(zip(d, oelement))]

    #e = []
    #for i in oelement:
    #    e.append(i)
    #arranged = e

    return arranged

def _get_ofelem(elements):
    elo = []
    for t in elements:
        if len(t) == 168:
            one = t[:152]
            two = t[152:168]
            elo.append(one + '\n')
            elo.append(two + '\n')
        else:
            elo.append(t + '\n')
    elements = ''.join(e for e in elo)
    return elements

def _get_onodes(lines):
    x = _get_marks(lines)
    pens = _get_lines(x[0], '-1,', lines)
    line = pens.replace(')\n', '')
    n_line = line.split('\n')
    last = n_line[-1]
    n_lines = n_line[:-1]


    return n_lines, last

def _get_oelem2(lines,indexes2,ind,sep_mat_one,seperate):
    #print(indexes2)
    x = _get_marks(lines)
    elems = _get_lines(x[1], '     -1\n', lines)
    elem = elems.replace(')\n', '')
    # print(elem)
    elem = elem.rstrip()
    # print(elem)
    all_splits = elem.split('\n')
    ret = []
    for i in range(0, len(all_splits), 2):
        ret.append(all_splits[i])
        try:
            ret[-1] += "\n" + all_splits[i + 1]
        except:
            pass
    quad = []
    for t in ret:
        if len(t) == 211:
            quad.append(t.replace('\n', ''))
        else:
            quad.append(t)
    ret = []
    for i in range(0, len(all_splits), 2):
        ret.append(all_splits[i])
        try:
            ret[-1] += "\n" + all_splits[i + 1]
        except:
            pass
    quad = []
    for t in ret:
        if len(t) == 211:
            quad.append(t.replace('\n', ''))
        else:
            quad.append(t.split('\n'))
    finals = _flat_List(quad)

    # seperate 50 HU and extra cdb
    if seperate == True:
        ele_final2 = []
        ele_final3 = []
        ele_final4 = []

        #print(sep_mat_one)

        finals2 = [finals[i] for i in ind[0]]
        finals3 = [finals[k] for k in ind[1]]

        #change indx in finals3
        indx = sorted([i[102:110] for i in finals3])



        nodes2 = []
        nodes3 = []
        nump = [[112, 120], [122, 130], [132, 140], [142, 150], [152, 160], [162, 170], [172, 180], [182, 190],
                [192, 200], [202, 210]]
        for j in finals2:
            for k in nump:
                nodes2.append(int(j[k[0]:k[1]].strip()))

        for f in finals3:
            for o in nump:
                nodes3.append(int(f[o[0]:o[1]].strip()))
        Nodes2 = sorted(list(set(nodes2)))
        Nodes3 = sorted(list(set(nodes3)))

        elindex2 = [e[2:10] for e in finals2]
        elindex3 = [e[2:10] for e in finals3]

        b = []
        for i in range(
                len(finals2)):  # len finals = 24(anzahl Elemente) ; len finals [i] = 210 ; x = indexes[i] + elindex[i]...... wird zu ....


            a = elindex2[i] + elindex2[i] + elindex2[i] + finals2[i][42:50] + finals2[i][52:60] + finals2[
                                                                                                                   i][
                                                                                                               62:70] + \
                finals2[i][72:80] + finals2[i][82:90] + finals2[i][92:100] + finals2[i][102:110] + finals2[i][112:120] + \
                finals2[i][122:130] + finals2[i][132:140] + finals2[i][142:150] + finals2[i][152:160] + finals2[i][
                                                                                                        162:170] + \
                finals2[i][172:180] + finals2[i][182:190] + finals2[i][192:200] + finals2[i][202:210]
            b.append(a)
            x = '       1' + b[i]
            ele_final2.append(x)
            y = sep_mat_one+b[i]
            ele_final4.append(y)



        s = [ind[1].index(i) for i in sorted(ind[1])]
        for i in range(len(finals3)):  # len finals = 24(anzahl Elemente) ; len finals [i] = 210 ; x = indexes[i] + elindex[i]...... wird zu ....

            x = indexes2[i] + elindex3[s[i]] + elindex3[s[i]] + elindex3[s[i]]+\
                finals3[s[i]][42:50] + finals3[s[i]][52:60] + finals3[s[i]][62:70] + \
                finals3[s[i]][72:80] + finals3[s[i]][82:90] + finals3[s[i]][92:100] +\
                indx[i] + finals3[s[i]][112:120] + finals3[s[i]][122:130] + finals3[s[i]][132:140] +\
                finals3[s[i]][142:150] + finals3[s[i]][152:160] + finals3[s[i]][162:170] + \
                finals3[s[i]][172:180] + finals3[s[i]][182:190] + finals3[s[i]][192:200] + finals3[s[i]][202:210]
            ele_final3.append(x)



    else:
        ele_final2 = []
        Nodes2 = []
        ele_final3 = []
        Nodes3 = []

    return ele_final2,ele_final3,ele_final4,Nodes2,Nodes3

def _get_oelem(lines, indexes):

    x = _get_marks(lines)
    elems = _get_lines(x[1], '     -1\n', lines)

    elem = elems.replace(')\n', '')
    #print(elem)
    elem = elem.rstrip()
    #print(elem)
    all_splits = elem.split('\n')
    ret = []
    for i in range(0, len(all_splits), 2):
        ret.append(all_splits[i])
        try:
            ret[-1] += "\n" + all_splits[i + 1]
        except:
            pass
    quad = []
    for t in ret:
        if len(t) == 211:
            quad.append(t.replace('\n', ''))
        else:
            quad.append(t)
    ret = []
    for i in range(0, len(all_splits), 2):
        ret.append(all_splits[i])
        try:
            ret[-1] += "\n" + all_splits[i + 1]
        except:
            pass
    quad = []
    for t in ret:
        if len(t) == 211:
            quad.append(t.replace('\n', ''))
        else:
            quad.append(t.split('\n'))
    finals = _flat_List(quad)

    #seperate 50 HU and extra cdb


    elindex = [e[2:10] for e in finals] #index for material
    ele = []

    for i in range(len(finals)):  #len finals = 24(anzahl Elemente) ; len finals [i] = 210 ; x = indexes[i] + elindex[i]...... wird zu ....

        if len(finals[i]) == 210:
            x = indexes[i] + elindex[i] + elindex[i] + elindex[i]+ finals[i][42:50] + finals[i][52:60] + finals[i][62:70] + \
                finals[i][72:80] + finals[i][82:90] + finals[i][92:100] + finals[i][102:110] + finals[i][112:120] + \
                finals[i][122:130] + finals[i][132:140] + finals[i][142:150] + finals[i][152:160] + finals[i][162:170] + \
                finals[i][172:180] + finals[i][182:190] + finals[i][192:200] + finals[i][202:210]

            ele.append(x)
        elif len(finals[i]) == 190:
            x = indexes[i] + elindex[i] +elindex[i] + elindex[i] + finals[i][42:50] + finals[i][52:60] + finals[i][62:70] + \
                finals[i][72:80] + finals[i][82:90] + finals[i][92:100] + finals[i][102:110] + finals[i][112:120] + \
                finals[i][122:130] + finals[i][132:140] + finals[i][142:150] + finals[i][152:160] + finals[i][162:170] + \
                finals[i][172:180] + finals[i][182:190]
            ele.append(x)
        elif len(finals[i]) == 170:
            x = indexes[i] + elindex[i] + '       1       1' + finals[i][42:50] + finals[i][52:60] + finals[i][62:70] + \
                finals[i][72:80] + finals[i][82:90] + finals[i][92:100] + finals[i][102:110] + finals[i][112:120] + \
                finals[i][122:130] + finals[i][132:140] + finals[i][142:150] + finals[i][152:160] + finals[i][162:170]
            ele.append(x)
        elif len(finals[i]) == 150:
            x = indexes[i] + elindex[i] + '       1       1' + finals[i][42:50] + finals[i][52:60] + finals[i][62:70] + \
                finals[i][72:80] + finals[i][82:90] + finals[i][92:100] + finals[i][102:110] + finals[i][112:120] + \
                finals[i][122:130] + finals[i][132:140] + finals[i][142:150]
            ele.append(x)
    #return finals
    #print(ele)
    return ele

