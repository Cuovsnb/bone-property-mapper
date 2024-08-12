#!/usr/bin/python
#

# ===============================
#

__all__ = ['import_parameters_xml', 'import_mesh', 'import_ct_data']

# -------------------------------------------------------------------------------
# Import modules
# -------------------------------------------------------------------------------
import os
import math
import numpy as np
import sys
import pydicom
from numpy import mean, array, concatenate, linspace, arange, size
from numpy import zeros, floor, diff, prod, matrix
import re
from general import _read_text, _remove_spaces, _refine_spaces
from general import _remove_eol_r
from general import _rot_matrix
from classes import linear_tet, quad_tet, linear_wedge, linear_hex
from classes import vtk_data, part, rpart, spart
#from lxml import etree

# -------------------------------------------------------------------------------
# Functions for importing parameters in XML
# -------------------------------------------------------------------------------
def import_parameters_xml(fle, gui_param):
    """Importing paramater data from XML"""
    # Parse all elements of XML file, IR change from XML to GUI
    #tree = etree.parse("Configuration_file.xml")
    #root = tree.getroot()
    #maps = [j for k in root.iter() for j in k.items()]
    # Storing and Mapping Data as Dictionary
    param = {}
    #b = []
    #a = []
    #c = []
    #intervals = []

    #for d in maps:
    #    if d[0] == 'b':
    #        b.append(d[1])
    #    elif d[0] == 'a':
    #        a.append(d[1])
    #    elif d[0] == 'c':
    #        c.append(d[1])
    #    elif d[0] == 'IntervalsType':
    #        intervals.append(d[1])
    #    elif d[0] == 'ROCalibrationCorrectionIsActive':
    #        param['calibrationCorrect'] = d[1]
    #    elif d[0] == 'ROIntercept':
    #        param['rhoQCTa'] = d[1]
    #    elif d[0] == 'ROSlope':
    #        param['rhoQCTb'] = d[1]
    #    elif d[0] == 'MinElasticity':
    #        param['minVal'] = d[1]
    #    elif d[0] == 'CalculationModality':
    #        param['integration'] = d[1]
    #    elif d[0] == 'StepsNumber':
    #        param['intSteps'] = d[1]
    #    elif d[0] == 'GapValue':
            #param['gapValue'] = d[1] IR
    #        param['gapValue'] = str(gui_param["gap"])
    #    elif d[0] == 'DensitySelection':
    #        param['groupingDensity'] = d[1]
    #    elif d[0] == 'PoissonRatio':
    #        param['poisson'] = d[1]

    #    # ------Possible Threshold Values------ #
    #    elif d[0] == 'RhoQCT1':
    #        param['rhoThresh1'] = d[1]
    #    elif d[0] == 'RhoQCT2':
    #        param['rhoThresh2'] = d[1]
    #    elif d[0] == 'RhoAsh1':
    #        param['Ethresh1'] = d[1]
    #    elif d[0] == 'RhoAsh2':
    #        param['Ethresh2'] = d[1]

    #param['rhoAshb1'] = b[1]
    #param['rhoAshb2'] = b[2]
    #param['rhoAshb3'] = b[3]
    #param['Eb1'] = b[5]
    #param['Eb2'] = b[6]
    #param['Eb3'] = b[7]
    #param['rhoAsha1'] = a[1]
    #param['rhoAsha2'] = a[2]
    #param['rhoAsha3'] = a[3]
    #param['Ea11'] = a[4]
    #param['Eb11'] = b[4]
    #param['Ec11'] = c[0]
    #param['Ea1'] = a[5]
    #param['Ea2'] = a[6]
    #param['Ea3'] = a[7]
    #param['Ec1'] = a[1]
    #param['Ec2'] = a[2]
    #param['Ec3'] = a[3]
    #param['numCTparam'] = intervals[0]
    #param['numEparam'] = intervals[0] #geändert auf 0 statt 1
    #param['rhoAsha11'] = a[0]
    #param['rhoAshb11'] = b[0]

    print(param)
    #IR
    param['rhoQCTa'] = "-0.003935729"
    param['rhoQCTb'] = "0.000791701"
    param["ansys"] = gui_param["ansys_version"]
    param["poisson"] = str(gui_param["poisson"])
    param["Ethresh1"] = str(gui_param["treshhold_value"])
    param['rhoAsha11'] = str(gui_param['b'])
    param['rhoAshb11'] = str(gui_param['a'])
    param['Eb11'] = str(gui_param['c_cort'])
    param['Ec11'] = str(gui_param['d_cort'])
    param['numEparam'] = str(gui_param["treshhold_check"])
    param['Eb1'] = str(gui_param['c_trab'])
    param['Ec1'] = str(gui_param['d_trab'])
    param['gapValue'] = str(gui_param["gap"])
    param['minVal'] = '1e-06'
    param['groupingDensity'] = 'Mean'
    param['integration'] = 'E'
    #param['intSteps'] = '4'
    #param['HUValues'] ='no'

    print(param)

    #poissonratio = input('      Assign custom poisson ratio? [Yes/No]: ')
    #if poissonratio.lower() == 'yes':
    #    param['poisson'] = input('      Enter new poisson value: ')
    #    print('      New poisson ratio value: ' + param['poisson'])
    #elif poissonratio == 'No' or 'no':
    #    print('      Poisson ratio default value: ' + param['poisson'])

    # Cleaning Values
    for k, v in param.items():
        if param[k].lower() == 'true':
            param[k] = True
        elif param[k].lower() == 'false':
            param[k] = False
        elif param[k].lower() == 'none':
            param[k] = None
        else:
            try:
                if type(v) == str:
                    param[k] = float(v)
            except:
                if k != 'integration' and k != 'ignore':
                    param[k] = v.lower()

    #HUparam = input('      Calibrate as per HU range?   [Yes/No]: ')
    #if HUparam.lower() == 'yes':
        #param['HUValues'] = 'yes'
        #print('      Default Calibration: Off')
        #param['rhoQCTar1'] = float(input('      Enter rhoQCTa for HU range 1-300: '))
        #param['rhoQCTbr1'] = float(input('      Enter rhoQCTb for HU range 1-300: '))
        #param['rhoQCTar2'] = float(input('\n      Enter rhoQCTa for range 301-700: '))
        #param['rhoQCTbr2'] = float(input('      Enter rhoQCTb for range 301-700: '))
        #param['rhoQCTar3'] = float(input('\n      Enter rhoQCTa for range 701-1000: '))
        #param['rhoQCTbr3'] = float(input('      Enter rhoQCTb for range 701-1000: '))
        #print('\n')
    #else:
    param['HUValues'] = 'no'
    print('      Default Caliberation: On')

    return param
    #return param as dic

# -------------------------------------------------------------------------------
# Functions for importing ansys realistic cdb mesh data
# -------------------------------------------------------------------------------
def import_mesh(fle, param):
    """ Reads in mesh data """
    print("fle:",fle, param)
    filetype = _what_mesh_filetype(fle)
    if filetype == '.inp':
        mesh = _import_abq_mesh(fle, param)
    elif filetype == '.ntr':
        mesh = _import_ntr_mesh(fle)
    elif filetype == '.neutral':
        mesh = _import_ntr_mesh(fle)
    elif filetype == '.cdb':
        mesh = _import_ans_mesh(fle)


    else:
        raise IOError('Mesh Import Error: Unknown filetype')
    #print(mesh[0].elements[0].indx)
    #print(mesh[0].elements[0].pts)
    #print(mesh[0].elements[100].nodes)

    #for i in range(len(mesh[0].elindex)):
        #a = mesh[0].elements[i].pts
        #print(a)

    #for i in range(10):
        #a = mesh[0].elements[i].pts
        #print(a)
   # print(mesh[0].elements.nodes)

    return mesh

def _what_mesh_filetype(fle):
    """ Determines mesh file format """
    lines = _read_text(fle)
    if '*Heading' in lines:
        return '.inp'
    elif '1G' in lines:
        return '.ntr'
    elif '.ntr' in fle:
        return '.ntr'
    elif '.neutral' in fle:
        return '.ntr'
    elif '.cdb' in fle:
        return '.cdb'
    else:
        return None

def _import_ans_mesh(lines):
    lines = _read_text(lines)
    #print(lines)
    # Get element types as different parts
    count = 0

    elements = _get_relem(lines)
    #print("elements:",elements)
    #print(elements[153])
    #print(len(elements))
    #print(elements[0])
    #print(type(elements[0]))


    part = []

    for e in elements:
        #print(e)
        #print(len(e))

        #print(len(e))
        #print(e)
        #print(e[0])
        if len(e) == 210:
            part.append('quad_tet')
        elif len(e) == 190:
            part.append('linear_hex')
        elif len(e) == 170:
            part.append('linear_wedge')
        elif len(e) == 150:
            part.append('linear_tet')
        #else:
           # print("no ele type avaible ")

    part = list(dict.fromkeys(part))
    part_data = [0] * len(part) #'list with objekts', objekt = for each elementtype

    for p in part:

        part_data[part.index(p)] = _get_part_data_r(lines, p, elements)


    #print(part_data[0].elements.pts) # wie auf elements zugreifen ? Objekt in Objekt ?
    #print(quad_tet.nodes)
    #print(type(part_data[0].elements))
    #f = 0
    #for i in part_data[0].elements:
        #f +=1
        #print(i.indx)
        #print(i.pts)
        #print(i.nodes)
        #if f == 2:
            #break
    return part_data

def _get_qt_elem(elements):
    eles = [i for i in elements if len(i) == 210]
    x = []
    t = []
    for i in eles:
        x.append(i[100:])
        t.append(i[10:21])
    ele = [" ".join(v.split()).split(" ") for v in x]
    #print("t:",t)
    return ele,t

def _get_lh_elem(elements):
    eles = [i for i in elements if len(i) == 190]
    x = []
    t = []
    for i in eles:
        x.append(i[100:])
        t.append(i[10:21])
    ele = [" ".join(v.split()).split(" ") for v in x]
    #print(t)
    return ele,t

def _get_lw_elem(elements):
    eles = [i for i in elements if len(i) == 170]
    x = []
    t = []
    for i in eles:
        x.append(i[100:])
        t.append(i[10:21])
    ele = [" ".join(v.split()).split(" ") for v in x]
    return ele,t

def _get_lt_elem(elements):
    eles = [i for i in elements if len(i) == 150]
    x = []
    t = []
    for i in eles:
        x.append(i[100:])
        t.append(i[10:21])
    ele = [" ".join(v.split()).split(" ") for v in x]
    return ele,t

def _get_ele_index(element):
    #elindexs = [e[10:21] for e in elements] #from materialindex to elindex
    #print("elindexes :", elindexs)
    #print(element)

    #elindex = [e.strip() for e in elindexs]

    elindex = [e.strip() for e in element]
    return elindex

def _get_part_data_r(lines, p, elements):
    nodes = _get_rnodes(lines) #elements noch alle elemente , nicht getrennt und zugeordnet , p = [quad_test, lin_tet] (z.B)

    if p == 'quad_tet':
        element = _get_qt_elem(elements)
    elif p == 'linear_hex':
        element = _get_lh_elem(elements)

    elif p == 'linear_wedge':
        element = _get_lw_elem(elements)
    elif p == 'linear_tet':
        element = _get_lt_elem(elements)

    elindex = _get_ele_index(element[1]) #get the right indxfrom selectet elements with len

    eletype = p

    numtype = _get_numtype(element[0], elindex)


    rparts = _create_rpart(p, element[0], eletype, nodes, elindex, numtype)

    return rparts

def _get_numtype(elements, elindex):
    ziped = list(zip(elindex, elements))
    #print("zip :", ziped)
    ss = {}
    for l in ziped:
        key, val = l
        ss.setdefault(key, []).append(val)
    typed = {}
    for k, v in ss.items():
        for e in v:
            if len(e) == 11:
                typed[k] = '187'
            elif len(e) == 9:
                typed[k] = '185'
            elif len(e) == 7:
                typed[k] = '183'
            elif len(e) == 5:
                typed[k] = '181'

    return typed

def _create_rpart(partname, elements, eletype, nodes, elindex, numtype, transform=[[0., 0., 0]], ignore=False):
    """ Creates part class from input data """
    # create the part
    #elements = all elements
    new_part = rpart(partname, eletype, elindex, numtype, transform, ignore)

    # add elements to part
    for e in elements:

        pts = [nodes[int(n)] for n in e[1:]]

        if eletype == 'quad_tet':
            ele = quad_tet(int(e[0]), pts, e[1:])
        elif eletype == 'linear_tet':
            ele = linear_tet(int(e[0]), pts, e[1:])
        elif eletype == 'linear_wedge':
            ele = linear_wedge(int(e[0]), pts, e[1:])
        elif eletype == 'linear_hex':
            ele = linear_hex(int(e[0]), pts, e[1:])
            #le = linear_hex(int(eleindex, nodes, nodeindex)

        new_part.add_relement(ele)
        #print(new_part.add_relement(ele))   #wieso None ???

    return new_part

# Getting Starting and Ending points for node & element extraction
def _get_marks(lines):
    if 'CMBLOCK' in lines:
        line = _get_lines('NBLOCK', 'CMBLOCK', lines)
        x = re.findall("\((.+)\)", line)
    else:
        x = re.findall("\((.+)\)", lines)

    return x

# Nodes extraction
def _get_rnodes(lines):
    #print("Halo")
    x = _get_marks(lines)
    pens = _get_lines(x[0], '-1,', lines)
    line = pens.replace(')\n', '')
    n_lines = line.split('\n')
    n_lines = n_lines[:-1]

    # Getting indexes
    index_node = []
    for i in n_lines:
        index_node.append(i[:9])
    node_indexes = []
    for t in range(len(index_node)):
        node_indexes.append(index_node[t].strip())
    node_indexes = map(int, node_indexes)
    node_indexes = list(node_indexes)

    # Getting coordinates
    x = []
    for i in n_lines:
        x.append(i[27:])
    no = []
    n = 21
    for l in x:
        no.append([l[index: index + n].strip() for index in range(0, len(l), n)])
    dimensions = [[float(i) for i in t] for t in no]


    # Conversion to dictionary
    zips = zip(node_indexes, dimensions)
    final = dict(zips)

    return final

def _get_relem(lines):

    x = _get_marks(lines)
    elem = _get_lines(x[1], '     -1\n', lines)
    #print(elem)
    elem = elem.replace(')\n', '')
    elem = elem.rstrip()
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
    rets = []
    for i in range(0, len(all_splits), 2):
        rets.append(all_splits[i])
        try:
            rets[-1] += "\n" + all_splits[i + 1]
        except:
            pass
    quad = []
    for t in ret:
        if len(t) == 211:
            quad.append(t.replace('\n', ''))
        else:
            quad.append(t.split('\n'))
    finals = _flat_List(quad)
    #print(finals[0])
    #print(finals[1])
    #print(finals[2])
    return finals

def _flat_List(nestedList):
    ### Converts a nested list to a flat list ###
    flatList = []
    # Iterate over all the elements in given list
    for elem in nestedList:
        # Check if type of element is list
        if isinstance(elem, list):
            # Extend the flat list by adding contents of this element (list)
            flatList.extend(_flat_List(elem))
        else:
            # Append the elemengt to the list
            flatList.append(elem)
    #print(flatList)
    return flatList

def _get_lines(startstring, endstring, lines):
    """ Identifies smallest matching string which meets start/end criteria """
    #print(startstring, endstring)
    # create pattern
    start = re.compile(startstring, re.IGNORECASE)
    end = re.compile(endstring, re.IGNORECASE)
    #print(start,end)

    # locate all incidences of start and end strings (case insensitive)
    start_i, end_i = [], []
    i, j = 0, 0
    if startstring == '':
        start_i = [0]
    else:
        while start.search(lines, i) != None:
            start_i.append(start.search(lines, i).end())
            i = start_i[-1] + 1
    if endstring == '$':
        end_i = [len(lines) - 1]
    else:
        while end.search(lines, j) != None:
            end_i.append(end.search(lines, j).start())
            j = end_i[-1] + 1

    # check start and end strings have been found
    if (start_i == []) | (end_i == []):
        return ''

    # determine minimum string
    comb = []
    [[comb.append((i, j)) for i in start_i] for j in end_i]
    d = [x[1] - x[0] for x in comb]
    pos = [i for i in d if i > 0]
    #print(pos)

    # check string has length
    if pos == []:
        return ''

    # return the lines
    startstop = comb[d.index(min(pos))]
    #print(startstop)
    final_lines = lines[startstop[0]:startstop[1]]
    #print(final_lines)

   # for i in final_lines:
        #print(i)



    return final_lines

# -------------------------------------------------------------------------------
# Functions for importing CT data
# -------------------------------------------------------------------------------
from numpy import int16, arange

def import_ct_data(ct_data): #path
    """ Reads CT dataset from specified file or directory """

    if ".vtk" in ct_data:
        data = _import_vtk_ct_data(ct_data)
    else:
        data = _import_dicoms(ct_data)

    return data

def _import_dicoms(ct_data):
    """ Imports dicom data and rearranges voxels to vtk lookup """


    dicom_order, dicom_data = _import_dicom_ct_data(ct_data)
    #print("dicom Data :",dicom_data)
    #print("dicom order :",dicom_order)

    lookup = _import_dicom_lookup(dicom_order, dicom_data, ct_data) # viel zeit, berechnet rescale slope and intercep to calculate HU
    X, Y, Z = _calc_dicom_xyz(dicom_data)

    #print("XYZ :", X,Y,Z)

    dicom = vtk_data(X, Y, Z, lookup) #creates Objekt from VTK class

    #print(dicom.lookup)
    #print("dicom dimen :", dicom.dimen)
    return dicom

def _import_dicom_ct_data(dicom_dir): #path
    """ Returns dicom CT image data and z order of files """
    #dicom_dir = path ct
    # Check DICOM dataset for any error
    _check_files(dicom_dir)

    # find files in folder
    dicom_files = os.listdir(dicom_dir)

    # sort by z position
    dicom_order = _sort_dicoms(dicom_files, dicom_dir)

    # gather data
    dicom_data = _gather_dicom_data(dicom_order, dicom_dir)


    # check validity of image metadata
    if _check_dicom_data(dicom_data) is not None:
        raise IOError("Dicom Import Error: " + _check_dicom_data(dicom_data))

    # return CT voxel data
    return dicom_order, dicom_data

def _import_dicom_lookup(dicom_order, dicom_data, ct_data):
    """ Iterates through sorted dicom slices and stores lookup data """
    #print(dicom_data)
    lookup = array([])
    #print(dicom_data, dicom_order)
    for s in dicom_order:

        dic = _read_dicom(s, ct_data)
        pixel_array = dic.pixel_array

        #np.set_printoptions(threshold=np.inf)
        #print(pixel_array)
        #print(pixel_array.shape)

        #print("Intercept, slope",dicom_data[14][0], dicom_data[15][0])

        # rescale intensity if required
        if (dicom_data[14][0] != 1.0) or (dicom_data[15][0] != 0.0):           #rescale slope and rescale intercept
            pixel_array = (pixel_array * dicom_data[14][0]) + dicom_data[15][0] # calculats HU
        # add data to lookup array
        lookup = concatenate((lookup, pixel_array.flatten())) # array with HU for aech voxel

    return lookup

def _calc_dicom_xyz(dicom_data):
    """ Calculates range of X, Y, and Z voxel co-ordinates from dicom information """
    X = arange(dicom_data[5][0], dicom_data[5][0] + (dicom_data[2][0] * dicom_data[0][0]), dicom_data[2][0])
    Y = arange(dicom_data[6][0], dicom_data[6][0] + (dicom_data[3][0] * dicom_data[1][0]), dicom_data[3][0])
    Z = arange(dicom_data[7][0], dicom_data[7][0] + (dicom_data[4][0] * len(dicom_data[0])), dicom_data[4][0])

    #print(X)
    #print(Z)
    return [X, Y, Z]

def _convert_dicom_ct_data(ct_data):
    """ Converts Dicom scans to a VTK ascii file format """

    dicom_order, dicom_data = _import_dicom_ct_data(ct_data)
    lookup = _import_dicom_lookup(dicom_order, dicom_data, ct_data)
    _write_vtk_header(dicom_data, ct_data + ".vtk")
    _write_vtk_lookup(lookup, ct_data)

def _check_files(dir):

    for root, dirs, files in os.walk(dir, topdown=False):
        for name in files:
            # Remove hidden files in directory
            if name.startswith('.'):
                os.remove(os.path.join(root, name))
            # Add .dcm extension if does not have in directory
            elif '.' not in name:
                os.rename(os.path.join(root, name), os.path.join(root, name) + '.dcm')
            # Give error of directory do not contain .dcm files
            elif '.dcm' not in name:
                print('Invalid DICOM dataset!')

def _gather_dicom_data(dicom_order, path):
    """ Collect image data for the dicoms in the folder """

    data = [[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []]

    #get SLice Thickness
    os.chdir(path)  # because without: permission denined, new packages (gdcm, pyiplo.....)
    s_s = []
    for p in range(2):
        d = pydicom.dcmread(path + "/" + dicom_order[p])
        s_s.append(d)
    ss = abs(round(float(s_s[0].ImagePositionPatient[2])-float(s_s[1].ImagePositionPatient[2]),4))

    for d in dicom_order:

        dic = _read_dicom(d, path)

        data[0].append(dic.Rows)
        data[1].append(dic.Columns)
        data[2].append(float(dic.PixelSpacing[0]))
        data[3].append(float(dic.PixelSpacing[1]))
        data[4].append(float(ss))
        #if hasattr(dic, 'SliceThickness'):
        #    data[4].append(float(dic.SliceThickness))
        #elif hasattr(dic, 'SpacingBetweenSlices'):
        #    data[4].append(float(dic.SpacingBetweenSlices))
        #else:
        #    data[4].append(float(dic.PixelSpacing[0]))
        # if dic == 'SliceThickness':
        #     data[4].append(float(dic.SliceThickness))
        # elif dic == 'SpacingBetweenSlices':
        #     data[4].append(float(dic.SpacingBetweenSlices))
        # elif dic != 'SpacingBetweenSlices' and dic != 'SliceThickness':
        #     data[4].append(float(dic.PixelSpacing[0]))
        data[5].append(float(dic.ImagePositionPatient[0]))
        data[6].append(float(dic.ImagePositionPatient[1]))
        data[7].append(float(dic.ImagePositionPatient[2]))
        data[8].append(float(dic.ImageOrientationPatient[0]))
        data[9].append(float(dic.ImageOrientationPatient[1]))
        data[10].append(float(dic.ImageOrientationPatient[2]))
        data[11].append(float(dic.ImageOrientationPatient[3]))
        data[12].append(float(dic.ImageOrientationPatient[4]))
        data[13].append(float(dic.ImageOrientationPatient[5]))
        data[14].append(float(dic.RescaleSlope))
        data[15].append(float(dic.RescaleIntercept))

    return data

def _check_dicom_data(dicoms):

    """ Check that the dicom images are valid """
    if len(dicoms[0]) == 1:
        return None

    else:
        if not len(set(dicoms[0])) == 1:
            return "Dicom rows different sizes"
        elif not len(set(dicoms[1])) == 1:
            return "Dicom columns different sizes"
        elif not len(set(dicoms[2])) == 1:
            return "Dicom images have different pixel spacings"
        elif not len(set(dicoms[3])) == 1:
            return "Dicom images have different pixel spacings"
        elif not len(set(dicoms[4])) == 1:
            return "Dicom images have different thicknesses"
        elif not len(set(dicoms[5])) == 1:
            return "Dicom images have different origins"
        elif not len(set(dicoms[6])) == 1:
            return "Dicom images have different origins"
        elif not len(set([round(c, 1) for c in diff(dicoms[7])])) == 1:
            return "Dicom slices are not sequential"
        elif not len(set(dicoms[8])) == 1:
            return "Dicom images have different row cosine orientation"
        elif not len(set(dicoms[9])) == 1:
            return "Dicom images have different row cosine orientation"
        elif not len(set(dicoms[10])) == 1:
            return "Dicom images have different row cosine orientation"
        elif not len(set(dicoms[11])) == 1:
            return "Dicom images have different col cosine orientation"
        elif not len(set(dicoms[12])) == 1:
            return "Dicom images have different col cosine orientation"
        elif not len(set(dicoms[13])) == 1:
            return "Dicom images have different col cosine orientation"
        elif [dicoms[8][0], dicoms[9][0], dicoms[10][0], dicoms[11][0], dicoms[12][0], dicoms[13][0]] != [1.0, 0.0, 0.0,
                                                                                                          0.0, 1.0,
                                                                                                          0.0]:
            print([dicoms[8][0], dicoms[9][0], dicoms[10][0], dicoms[11][0], dicoms[12][0], dicoms[13][0]])

            return "ImageOrientation parameter must be [1,0,0,0,1,0]"

def _read_dicom(fle, path): #fle = filename/picture name ; path = name of file bzw. Ordnername
    """ Reads dicom file information """
    #print(fle)
    #d = pydicom.dcmread(os.path.join(path, fle), force=True)
    #change from IR because BG Dicom are other
    os.chdir(path)#because without: permission denined, new packages (gdcm, pyiplo.....)
    d = pydicom.dcmread(path + "/" + fle)

    return d

def _sort_dicoms(dicom_files, path):
    """ Iterates through dicom files and returns sorted z order """

    z = [pydicom.dcmread(os.path.join(path, d), force=True).ImagePositionPatient[2] for d in dicom_files] #dcmread to file_read
    dicom_order = list(zip(*sorted(zip(z, dicom_files))))[1]
    return dicom_order

def _dicom_xyz_data(dicom_data):
    """ Find X,Y,Z voxel grid increments from dicom data """

    rows = dicom_data[0][0]
    columns = dicom_data[1][0]
    slices = dicom_data[7]
    xspacing, yspacing = [dicom_data[2][0], dicom_data[3][0]]
    slice_thickness = dicom_data[4][0]

    # calculate X-coordinates
    xstart = dicom_data[5][0]
    xstop = xstart + (columns * yspacing)
    X = arange(xstart, xstop, yspacing)

    # calculate Y-coordinates
    ystart = dicom_data[6][0]
    ystop = ystart + (rows * xspacing)
    Y = arange(ystart, ystop, xspacing)

    # calculate Z-coordinates
    zstart = min(slices)
    zstop = zstart + (slice_thickness * len(dicom_data[0]))
    Z = arange(zstart, zstop, slice_thickness)

    return X, Y, Z

def _import_vtk_ct_data(ct_data):
    """ Creates python array data from vtk file """

    # read in text
    lines = _read_text(ct_data)
    lines = _remove_eol_r(lines)
    header = lines.split('LOOKUP')[0]

    # read in X data
    xlines = _refine_vtk_lines(header, 'X_COORDINATES \d+', 'double', 'Y_COORDINATES')
    X = [float(x) for x in xlines]

    # read in Y data
    ylines = _refine_vtk_lines(header, 'Y_COORDINATES \d+', 'double', 'Z_COORDINATES')
    Y = [float(y) for y in ylines]

    # read in Z data
    zlines = _refine_vtk_lines(header, 'Z_COORDINATES \d+', 'double', 'CELL_DATA')
    Z = [float(z) for z in zlines]

    # read in lookup data
    lookup_lines = _refine_vtk_lines(lines, 'LOOKUP_TABLE', 'default', '')
    lookup = [int(l) for l in lookup_lines]

    # create vtk class
    vtk = vtk_data(X, Y, Z, lookup)

    # return data

    return vtk

def _refine_vtk_lines(lines, key1, key2, key3):
    """ Find lines bewteen key words and remove end-of-line characters """

    # find lines
    p = re.compile(key1 + ' ' + key2 + '\n([-.\d\n ]+)' + key3)
    if len(p.findall(lines)) > 0:
        lines = p.findall(lines)[0]
    else:
        p = re.compile(key1 + ' ' + 'float' + '\n([-.\d\n ]+)' + key3)
        if len(p.findall(lines)) > 0:
            lines = p.findall(lines)[0]
        else:
            raise ValueError("Error reading VTK header, unrecognised format")
    # find numbers
    p2 = re.compile(r'[-.\d]+')

    return p2.findall(lines)

def _write_vtk_header(dicom_data, fle):
    """ Writes VTK data header """

    x, y, z = _dicom_xyz_data(dicom_data)
    with open(fle, 'w') as oupf:
        dimen = [len(x), len(y), len(z)]
        oupf.write('# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET RECTILINEAR_GRID\n')
        oupf.write('DIMENSIONS ' + repr(dimen[0]) + ' ' + repr(dimen[1]) + ' ' + repr(dimen[2]) + '\n')
        oupf.write('X_COORDINATES ' + repr(dimen[0]) + ' double\n')
        _create_vtk_coord_str(x, 9, oupf)
        oupf.write('Y_COORDINATES ' + repr(dimen[1]) + ' double\n')
        _create_vtk_coord_str(y, 9, oupf)
        oupf.write('Z_COORDINATES ' + repr(dimen[2]) + ' float\n')
        _create_vtk_coord_str(z, 9, oupf)
        oupf.write('CELL_DATA ' + repr(prod([d - 1 for d in dimen])) + '\n')
        oupf.write('POINT_DATA ' + repr(prod(dimen)) + '\n')
        oupf.write('SCALARS scalars short\n')
        oupf.write('LOOKUP_TABLE default\n')

def _write_vtk_lookup(lookup, path):
    """ Writes VTK lookup data """

    fle = path + ".vtk"
    count = 0
    for d in lookup:
        with open(fle, 'a') as oupf:
            if count < 8:
                oupf.write(repr(int(d)) + ' ')
                count += 1
            else:
                oupf.write(repr(int(d)) + '\n')
                count = 0
    with open(fle, 'a') as oupf:
        oupf.write('\n')

def _create_vtk_coord_str(coords, max_num, oupf, perform_round=True):
    """ Create string for numbers (rounded if specified in chunks of max_num and write to file """

    # perform rounding to 4 decimal places (default)
    if perform_round:
        coords = ["%0.4f" % (c,) for c in coords]
    else:
        coords = [repr(c) for c in coords]
    #print(coords)
    #print(len(coords))
    #print(type(len(coords)))
    #print(max_num)
    #print(type(max_num))
    #print(len(coords) / max_num)
    ## create string

    for c in range(math.ceil((len(coords) / max_num))):
        coords_str = " ".join(coords[c * max_num: (c + 1) * max_num]) + '\n'
        oupf.write(coords_str)
    if (len(coords) % max_num) > 0:
        coords_str = " ".join(coords[-(len(coords) % max_num):])
        oupf.write(coords_str)
        oupf.write('\n')

    return

def test_output(lines):
    """ Check lines information has been imported """

    if lines is None:
        raise ValueError("Error reading input file, check format")
    if lines == []:
        raise ValueError("Error reading input file, check format")

#################################################################################################################
#Test Files Pelvis

#
#import_mesh(r"C:\Users\8iries\Desktop\Praktikum\Programm\mesh_P5.cdb")
#import_ct_data(r"C:\Users\8iries\Desktop\Praktikum\Programm\series-01")
#import_parameters_xml(r"C:\Users\8iries\Desktop\Praktikum\Programm\params.xml")

##################################################################################################################
#Test File LCP Plate

#vtk = import_ct_data(r"C:\Users\8iries\Desktop\Praktikum\Test_Programm\3B4C66CC")
#param = import_parameters_xml(r"C:\Users\8iries\Desktop\Praktikum\Programm\params.xml")
#parts = import_mesh(r"C:\Users\8iries\Desktop\Praktikum\Test_Programm\mesh_body.cdb")

###kleiner tetraeder
#import_mesh(r"C:\Users\8iries\Desktop\Praktikum\Test_Programm\Mesh_tet\mesh_tet_com.cdb")
#import_ct_data(r"C:\Users\8iries\Desktop\Praktikum\Test_Programm\test_lcp")

#_import_ans_mesh(r'C:\Users\8iries\Desktop\Praktikum\Test_Programm\mesh_body.cdb')
#_convert_dicom_ct_data(r"C:\Users\8iries\Desktop\Praktikum\Test_Programm\Realistic\Ct_dataset")
#grid.plot()


##################################################################################################################
#Test_Data_Bonemat


#vtk = import_ct_data(r"C:\Users\8iries\Desktop\Praktikum\Programm\CT_dataset")
#param = import_parameters_xml(r"C:\Users\8iries\Desktop\Praktikum\Programm\params.xml")
#parts = import_mesh(r"C:\Users\8iries\Desktop\Praktikum\Programm\ct2.cdb", param)


#Mesh = open(r'pathcdb.txt', 'r')  # get path from txt from GUI input
#parts = Mesh.read()

#Param = open(r'pathxml.txt', 'r')
#param = Param.read()


#import_mesh(parts,param)


#C:/Users/8iries/Desktop/Ba/FB188694
#C:/Users/8iries/Desktop/Ba/Gelantine_Block_files/2xGelan_Block_grob_hex.cdb
#param = "p "
#parts = import_mesh(r"C:/Users/8iries/Desktop/Ba/Gelantine_Block_files/2xGelan_Block_grob_hex.cdb", param)