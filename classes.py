#!/usr/bin/python
#
# py_bonemat_abaqus - class definitions
# =====================================
#
# Created by Elise Pegg, University of Bath

__all__ = ['vtk_data', 'linear_tet', 'quad_tet', 'linear_wedge', 'linear_hex', 'part', 'rpart']

# -------------------------------------------------------------------------------
# Import modules
# -------------------------------------------------------------------------------
#from concurrent.futures import process
from numpy.linalg import det
from numpy import mean, arange, matrix, array
from copy import deepcopy
from bisect import bisect_left, bisect_right
from operator import itemgetter
#from itertools import product
#import numpy as np
#from multiprocessing import Pool
#import pyvista

# -------------------------------------------------------------------------------
# Part class abaqus mesh
# -------------------------------------------------------------------------------



class part:
    """ Part class which represents a mesh"""

    __slots__ = ("name", "elements", "ele_name", "ele_type", "moduli", "transform", "ignore")

    def __init__(self, name, ele_name, ele_type, transform=[[0., 0., 0]], ignore=False):
        self.name = name
        self.ele_type = ele_type
        self.ele_name = ele_name
        self.elements = []
        self.moduli = []
        self.transform = transform
        self.ignore = ignore

    def add_element(self, ele):
        # add element to part
        self.elements.append(ele)


# -------------------------------------------------------------------------------
# Part class simple cdb mesh
# -------------------------------------------------------------------------------
class spart:
    """ Part class which represents a realistic mesh"""

    __slots__ = ("elements", "ele_type", "moduli", "transform", "ignore")

    def __init__(self, ele_type, transform=[[0., 0., 0]], ignore=False):
        self.ele_type = ele_type
        self.elements = []
        self.moduli = []
        self.transform = transform
        self.ignore = ignore

    def add_relement(self, ele):
        # add element to part
        self.elements.append(ele)


# -------------------------------------------------------------------------------
# Part class realistic cdb mesh
# -------------------------------------------------------------------------------

class rpart:
    """ Part class which represents a realistic mesh"""

    __slots__ = ("name", "elements", "ele_type", "elindex", "numtype", "moduli","moduli2","sep_mat", "transform", "ignore", "density", "density2")

    def __init__(self, name, ele_type, elindex, numtype, transform=[[0., 0., 0]], ignore=False):
        self.name = name
        self.ele_type = ele_type
        self.elindex = elindex #eig. Partindex....
        self.numtype = numtype
        self.elements = []
        self.moduli = []
        self.moduli2 = []
        self.sep_mat = []
        self.transform = transform
        self.ignore = ignore
        self.density = []
        self.density2 = []


    def add_relement(self, ele):
        # add element to part
        self.elements.append(ele)

# -------------------------------------------------------------------------------
# CT data class
# -------------------------------------------------------------------------------


class vtk_data:
    """ Vtk ct class """

    __slots__ = ("x", "y", "z", "dimen", "vox", "lookup", "HU", "E", "E_int", "E_m")

    def __init__(self, x, y, z, lookup): #x,y,z coordinates from ct data/ grid (first oone Image Position..mittelpunkt der Voxel) , loo9kup = array with HU data for each voxel
        self.x = x
        self.y = y
        self.z = z
        self.dimen = [len(self.x), len(self.y), len(self.z)] #dimension/ Anzahl aller Koordinaten je achse
        self.vox = set([])
        self.lookup = array(lookup) #HU
        #IR
        self.HU = []
        self.E = []
        self.E_int = []
        self.E_m = []

    def get_voxels(self, tet):
        # define variables
        voxels = []
        in_element = False

        # check element is within voxel range
        tet_pts = list(zip(*[tet.pts[0], tet.pts[1], tet.pts[2], tet.pts[3]]))
        in_vox = True
        if (min(tet_pts[0]) < self.x[0]):
            raise ValueError('Mesh Error: Element(s) outside of CT data range. X coord min=' + repr(min(tet_pts[0])))
        if (max(tet_pts[0]) > self.x[-1]):
            raise ValueError('Mesh Error: Element(s) outside of CT data range. X coord max=' + repr(max(tet_pts[0])))
        if (min(tet_pts[1]) < self.y[0]):
            raise ValueError('Mesh Error: Element(s) outside of CT data range. Y coord min=' + repr(min(tet_pts[1])))
        if (max(tet_pts[1]) > self.y[-1]):
            raise ValueError('Mesh Error: Element(s) outside of CT data range. Y coord max=' + repr(max(tet_pts[1])))
        if (min(tet_pts[2]) < self.z[0]):
            raise ValueError('Mesh Error: Element(s) outside of CT data range. Z coord min=' + repr(min(tet_pts[2])))
        if (max(tet_pts[2]) > self.z[-1]):
            raise ValueError('Mesh Error: Element(s) outside of CT data range. Z coord max=' + repr(max(tet_pts[3])))

        # find grid data
        [[xstart, xend], [ystart, yend], [zstart, zend]] = self.find_grid(tet_pts)

        # test if element is smaller than voxel
        if (xstart == xend - 1):
            if (ystart == yend - 1):
                if (zstart == zend - 1):
                    in_element = True
                    ele_centroid = mean(tet_pts, 0).tolist()
                    x = [bisect_right(self.x, ele_centroid[0]) - 1, bisect_left(self.x, ele_centroid[0])]
                    y = [bisect_right(self.y, ele_centroid[1]) - 1, bisect_left(self.y, ele_centroid[1])]
                    z = [bisect_right(self.z, ele_centroid[2]) - 1, bisect_left(self.z, ele_centroid[2])]
                    eight_vox = _xyz_comb(x, y, z)
                    for v in eight_vox:
                        voxel_lookup = _calc_lookup(v[0], v[1], v[2], self.dimen)
                        self.vox.add(voxel_lookup)
                        voxels.append(voxel_lookup)

        # for each voxel test if in tet, if is return lookup
        if in_element == False:
            for x in range(xstart, xend + 1):
                for y in range(ystart, yend + 1):
                    for z in range(zstart, zend + 1):
                        if _calc_lookup(x, y, z, self.dimen) not in self.vox:
                            if tet.in_tet([self.x[x], self.y[y], self.z[z]]):
                                voxel_lookup = _calc_lookup(x, y, z, self.dimen)
                                self.vox.add(voxel_lookup)
                                voxels.append(voxel_lookup)
        return voxels

    def find_grid(self, pts, box=False): # pts = liste mit coordinates from scaled nodes, self.x... Coordinates from ct,
        # calculate grid, pts = calculatet x from class qud_tet
        #print(pts)
        #print(self.x)
        #print(type(self.x))
        #print(pts[0])
        #print(type(pts[0]))

        xstart = bisect_right(self.x, min(pts[0])) - 1
        xend = bisect_left(self.x, max(pts[0]))
        ystart = bisect_right(self.y, min(pts[1])) - 1
        yend = bisect_left(self.y, max(pts[1]))
        zstart = bisect_right(self.z, min(pts[2])) - 1
        zend = bisect_left(self.z, max(pts[2]))
        # if need start and end to be different (box = true) add one to end
        if box:
            if xstart == xend:
                if xend == (len(self.x) - 1):
                    xstart -= 1
                else:
                    xend += 1
            if ystart == yend:
                if yend == (len(self.y) - 1):
                    ystart -= 1
                else:
                    yend += 1
            if zstart == zend:
                if zend == (len(self.z) - 1):
                    zstart -= 1
                else:
                    zend += 1

        #print("start-end-CT: ", [[xstart, xend], [ystart, yend], [zstart, zend]]) #Index der vo selx. und self.y ... index dr Liste von Ursprungs CT Koordinaten
        return [[xstart, xend], [ystart, yend], [zstart, zend]]

    def interpolateScalar_hu(self,xyz,rhoQCT):

        [xi, yi, zi] = self.find_grid([[xyz[0]], [xyz[1]], [xyz[2]]],
                                      True)  # index values (find grid (x,y,z) with x,y,z from x = calculatetd coordinates (normal coordinates)

        box = [[xi[0], yi[0], zi[0]],
               [xi[0], yi[0], zi[1]],
               [xi[0], yi[1], zi[0]],
               [xi[0], yi[1], zi[1]],
               [xi[1], yi[0], zi[0]],
               [xi[1], yi[0], zi[1]],
               [xi[1], yi[1], zi[0]],
               [xi[1], yi[1], zi[1]]]
        # print("bx: ",box)
        # box return the indexes of the array from self.x , self.y, self.z ....., (orig. CT-Coordinates)

        # define the origin/Ursprung , x_i [0] = x_start, ... einVoxel der 8, mit den Startdaten
        origin_indx = 0
        origin = [self.x[xi[origin_indx]],
                  self.y[yi[origin_indx]],
                  self.z[zi[origin_indx]]]
        # print("origin: ",origin)

        # for each corner of the box, find the scalar value (HU ???)
        # b[i] = coordinates of corner of box, dimension, anzahl der koordinaten , self.lookup = array with HU Value
        # self.lookup[index of array]
        # print("self.lookup :", self.lookup, self.lookup.shape)
        scalar = [self.lookup[_calc_lookup(b[0], b[1], b[2], self.dimen)] for b in
                  box]  # scalar = HU = List of HU, self.lookup[position/indx in array for box corners]
        # print(self.dimen)
        # value of array to HU
        #self.HU.append(scalar)


        density = [rhoQCT(q) for q in scalar]

        differences = [self.x[xi[1]] - origin[0],
                       self.y[yi[1]] - origin[1],
                       self.z[zi[1]] - origin[2]]

        # calculate the relative position (%) of the co-ordinate within the box, xyz[n] = scaled kartesisch coordinates of node
        rel_xyz = [(xyz[n] - origin[n]) / differences[n] for n in [0, 1, 2]]

        # interpolate density
        cc00 = (rel_xyz[0] * (density[4] - density[0])) + density[0]
        cc01 = (rel_xyz[0] * (density[5] - density[1])) + density[1]
        cc11 = (rel_xyz[0] * (density[7] - density[3])) + density[3]
        cc10 = (rel_xyz[0] * (density[6] - density[2])) + density[2]
        cc0 = (rel_xyz[1] * (cc10 - cc00)) + cc00
        cc1 = (rel_xyz[1] * (cc11 - cc01)) + cc01
        cc = (rel_xyz[2] * (cc1 - cc0)) + cc0

        #interpolste HU
        ccc00 = (rel_xyz[0] * (scalar[4] - scalar[0])) + scalar[0]
        ccc01 = (rel_xyz[0] * (scalar[5] - scalar[1])) + scalar[1]
        ccc11 = (rel_xyz[0] * (scalar[7] - scalar[3])) + scalar[3]
        ccc10 = (rel_xyz[0] * (scalar[6] - scalar[2])) + scalar[2]
        ccc0 = (rel_xyz[1] * (ccc10 - ccc00)) + ccc00
        ccc1 = (rel_xyz[1] * (ccc11 - ccc01)) + ccc01
        ccc = (rel_xyz[2] * (ccc1 - ccc0)) + ccc0


        return ccc,cc

    def interpolateScalar2(self, xyz,diff,count,R, cdb_indx, e_hu, modulus, rhoAsh, rhoQCT): #x,y,z = coordinatef of scaled node coordinates 8kartensicsch)
        # calculate CT box surrounding the point ##in ct grid
        #print("HU :", self.lookup)
        #x-start_x-end index values for box (auch für y,z)
        [xi, yi, zi] = self.find_grid([[xyz[0]], [xyz[1]], [xyz[2]]], True)  # index values (find grid (x,y,z) with x,y,z from x = calculatetd coordinates (normal coordinates)

        box = [[xi[0], yi[0], zi[0]],
               [xi[0], yi[0], zi[1]],
               [xi[0], yi[1], zi[0]],
               [xi[0], yi[1], zi[1]],
               [xi[1], yi[0], zi[0]],
               [xi[1], yi[0], zi[1]],
               [xi[1], yi[1], zi[0]],
               [xi[1], yi[1], zi[1]]]
        #print("bx: ",box)
        # box return the indexes of the array from self.x , self.y, self.z ....., (orig. CT-Coordinates)

        # define the origin/Ursprung , x_i [0] = x_start, ... einVoxel der 8, mit den Startdaten
        origin_indx = 0
        origin = [self.x[xi[origin_indx]],
                  self.y[yi[origin_indx]],
                  self.z[zi[origin_indx]]]

        # calculate the dimensions of the box, length of the axis
        differences = [self.x[xi[1]] - origin[0],
                       self.y[yi[1]] - origin[1],
                       self.z[zi[1]] - origin[2]]

        # calculate the relative position (%) of the co-ordinate within the box, xyz[n] = scaled kartesisch coordinates of node
        rel_xyz = [(xyz[n] - origin[n]) / differences[n] for n in [0, 1, 2]]


        # for each corner of the box, find the scalar value (HU ???)
        #b[i] = coordinates of corner of box, dimension, anzahl der koordinaten , self.lookup = array with HU Value

        scalar = [self.lookup[_calc_lookup(b[0], b[1], b[2], self.dimen)] for b in box] #scalar = HU = List of HU, self.lookup[position/indx in array for box corners]
        #print(scalar)

        C00 = (rel_xyz[0] * (scalar[4] - scalar[0])) + scalar[0]
        C01 = (rel_xyz[0] * (scalar[5] - scalar[1])) + scalar[1]
        C11 = (rel_xyz[0] * (scalar[7] - scalar[3])) + scalar[3]
        C10 = (rel_xyz[0] * (scalar[6] - scalar[2])) + scalar[2]
        C0 = (rel_xyz[1] * (C10 - C00)) + C00
        C1 = (rel_xyz[1] * (C11 - C01)) + C01
        C = (rel_xyz[2] * (C1 - C0)) + C0



        # apply modulus calculation (only changes value if V3) IR: for 8 box corners a value
        density = [rhoQCT(q) for q in scalar] # first density rhoQct


        #if diff[count] == 0:
            #scalars = [e_hu(p,R) for p in density]
            #if count <= 10:
                #print("U50")
                #print("HU ",C, scalar)
                #print("E ", scalars)


        if C <= 50:
            scalars = [e_hu(p,R) for p in density]
        else:
            scalars = [modulus(rhoAsh(rhoQCT(s))) for s in scalar]  #E-Modul, 8 values for one node (box values of grid)


        #print(scalars)
        #self.E.append(scalars)


        #print("scalae E :", scalars)
        #e_m = sum(scalars)/len(scalars)
        #self.E_m.append(e_m)
        #print("E_m :", e_m)

        #print(scalars)

        # interpolate scalar E
        c00 = (rel_xyz[0] * (scalars[4] - scalars[0])) + scalars[0]
        c01 = (rel_xyz[0] * (scalars[5] - scalars[1])) + scalars[1]
        c11 = (rel_xyz[0] * (scalars[7] - scalars[3])) + scalars[3]
        c10 = (rel_xyz[0] * (scalars[6] - scalars[2])) + scalars[2]
        c0 = (rel_xyz[1] * (c10 - c00)) + c00
        c1 = (rel_xyz[1] * (c11 - c01)) + c01
        c = (rel_xyz[2] * (c1 - c0)) + c0

        #if diff[count] == 0:
            #scalars = [e_hu(p,R) for p in density]
            #if count <= 30:
             #   print("U50")
             #   print("HU ",C, scalar)
             #   print("E ", scalars)
             #   print("c", c)
        #print(c)
        #if diff[count] == 0:
        #    print("E :",c)
        #with open ("c.txt", "w") as f:
            #f.write(str(c))

        self.E_int.append(c) # 10 Werte pro Element, mit interpoliertem wert von ct-box um den Eckknoten
        #print("interpoliert E :", c)

        # interpolate density
        cc00 = (rel_xyz[0] * (density[4] - density[0])) + density[0]
        cc01 = (rel_xyz[0] * (density[5] - density[1])) + density[1]
        cc11 = (rel_xyz[0] * (density[7] - density[3])) + density[3]
        cc10 = (rel_xyz[0] * (density[6] - density[2])) + density[2]
        cc0 = (rel_xyz[1] * (cc10 - cc00)) + cc00
        cc1 = (rel_xyz[1] * (cc11 - cc01)) + cc01
        cc = (rel_xyz[2] * (cc1 - cc0)) + cc0

        #interpolte HU
        ccc00 = (rel_xyz[0] * (scalar[4] - scalar[0])) + scalar[0]
        ccc01 = (rel_xyz[0] * (scalar[5] - scalar[1])) + scalar[1]
        ccc11 = (rel_xyz[0] * (scalar[7] - scalar[3])) + scalar[3]
        ccc10 = (rel_xyz[0] * (scalar[6] - scalar[2])) + scalar[2]
        ccc0 = (rel_xyz[1] * (ccc10 - ccc00)) + ccc00
        ccc1 = (rel_xyz[1] * (ccc11 - ccc01)) + ccc01
        ccc = (rel_xyz[2] * (ccc1 - ccc0)) + ccc0

        #print(ccc)
        #print(c)
        return c, cc, ccc # scales, density,hu ein Wert eines Punktes eines Elements aus  10 interpoliertwn punkten Punkten

    def interpolateScalar(self, xyz, modulus, rhoAsh, rhoQCT): #x,y,z = coordinatef of scaled node coordinates 8kartensicsch)
        # calculate CT box surrounding the point ##in ct grid
        #print("HU :", self.lookup)
        #x-start_x-end index values for box (auch für y,z)
        [xi, yi, zi] = self.find_grid([[xyz[0]], [xyz[1]], [xyz[2]]], True)  # index values (find grid (x,y,z) with x,y,z from x = calculatetd coordinates (normal coordinates)

        box = [[xi[0], yi[0], zi[0]],
               [xi[0], yi[0], zi[1]],
               [xi[0], yi[1], zi[0]],
               [xi[0], yi[1], zi[1]],
               [xi[1], yi[0], zi[0]],
               [xi[1], yi[0], zi[1]],
               [xi[1], yi[1], zi[0]],
               [xi[1], yi[1], zi[1]]]
        #print("bx: ",box)
        # box return the indexes of the array from self.x , self.y, self.z ....., (orig. CT-Coordinates)

        # define the origin/Ursprung , x_i [0] = x_start, ... einVoxel der 8, mit den Startdaten
        origin_indx = 0
        origin = [self.x[xi[origin_indx]],
                  self.y[yi[origin_indx]],
                  self.z[zi[origin_indx]]]

        # calculate the dimensions of the box, length of the axis
        differences = [self.x[xi[1]] - origin[0],
                       self.y[yi[1]] - origin[1],
                       self.z[zi[1]] - origin[2]]

        # calculate the relative position (%) of the co-ordinate within the box, xyz[n] = scaled kartesisch coordinates of node
        rel_xyz = [(xyz[n] - origin[n]) / differences[n] for n in [0, 1, 2]]


        # for each corner of the box, find the scalar value (HU ???)
        #b[i] = coordinates of corner of box, dimension, anzahl der koordinaten , self.lookup = array with HU Value

        scalar = [self.lookup[_calc_lookup(b[0], b[1], b[2], self.dimen)] for b in box] #scalar = HU = List of HU, self.lookup[position/indx in array for box corners]



        # apply modulus calculation (only changes value if V3) IR: for 8 box corners a value
        density = [rhoQCT(q) for q in scalar] # first density rhoQct



        scalars = [modulus(rhoAsh(rhoQCT(s))) for s in scalar]  # E-Modul, 8 values for one node (box values of grid)
        #self.E.append(scalars)


        #print("scalae E :", scalars)
        #e_m = sum(scalars)/len(scalars)
        #self.E_m.append(e_m)
        #print("E_m :", e_m)

        #print(scalars)

        # interpolate scalar E
        c00 = (rel_xyz[0] * (scalars[4] - scalars[0])) + scalars[0]
        c01 = (rel_xyz[0] * (scalars[5] - scalars[1])) + scalars[1]
        c11 = (rel_xyz[0] * (scalars[7] - scalars[3])) + scalars[3]
        c10 = (rel_xyz[0] * (scalars[6] - scalars[2])) + scalars[2]
        c0 = (rel_xyz[1] * (c10 - c00)) + c00
        c1 = (rel_xyz[1] * (c11 - c01)) + c01
        c = (rel_xyz[2] * (c1 - c0)) + c0

        self.E_int.append(c) # 10 Werte pro Element, mit interpoliertem wert von ct-box um den Eckknoten
        #print("interpoliert E :", c)

        # interpolate density
        cc00 = (rel_xyz[0] * (density[4] - density[0])) + density[0]
        cc01 = (rel_xyz[0] * (density[5] - density[1])) + density[1]
        cc11 = (rel_xyz[0] * (density[7] - density[3])) + density[3]
        cc10 = (rel_xyz[0] * (density[6] - density[2])) + density[2]
        cc0 = (rel_xyz[1] * (cc10 - cc00)) + cc00
        cc1 = (rel_xyz[1] * (cc11 - cc01)) + cc01
        cc = (rel_xyz[2] * (cc1 - cc0)) + cc0

        #interpolte HU
        ccc00 = (rel_xyz[0] * (scalar[4] - scalar[0])) + scalar[0]
        ccc01 = (rel_xyz[0] * (scalar[5] - scalar[1])) + scalar[1]
        ccc11 = (rel_xyz[0] * (scalar[7] - scalar[3])) + scalar[3]
        ccc10 = (rel_xyz[0] * (scalar[6] - scalar[2])) + scalar[2]
        ccc0 = (rel_xyz[1] * (ccc10 - ccc00)) + ccc00
        ccc1 = (rel_xyz[1] * (ccc11 - ccc01)) + ccc01
        ccc = (rel_xyz[2] * (ccc1 - ccc0)) + ccc0


        return c, cc, ccc # scales, density,hu ein Wert eines Punktes eines Elements aus  10 interpoliertwn punkten Punkten


# -------------------------------------------------------------------------------
# Miscellaneous functions
# -------------------------------------------------------------------------------
def _calc_lookup(x, y, z, dimen): # x,y,z = from box the x,y,z value , = index of list from ct grid, calculates position for list from scalar
    #print("xHallo :", x)
    #print("y :", y)
    #print("z :", z)
    #print("dim :", dimen)
    #print(x + (y * dimen[0]) + (z * dimen[0] * dimen[1]))
    return x + (y * dimen[0]) + (z * dimen[0] * dimen[1])


def _xyz_comb(x, y, z):
    """ Find voxels within x, y, z ranges """
    combinations = []
    for i in x:
        for j in y:
            for k in z:
                combinations.append([i, j, k])
    return combinations


# -------------------------------------------------------------------------------
# Element classes
# -------------------------------------------------------------------------------
# linear tetrahedron (C3D4)
class linear_tet:
    """ Linear tetrahedral element class """

    __slots__ = ("indx", "pts", "nodes", "volume")

    def __init__(self, indx, pts, nodes):
        self.indx = indx
        self.pts = pts
        self.nodes = nodes
        self.volume = 0

    def jacobian(self):
        # create base matrix from vertices for in_tet test
        return [[1, 1, 1, 1],
                [self.pts[0][0], self.pts[1][0], self.pts[2][0], self.pts[3][0]],
                [self.pts[0][1], self.pts[1][1], self.pts[2][1], self.pts[3][1]],
                [self.pts[0][2], self.pts[1][2], self.pts[2][2], self.pts[3][2]]]

    def in_tet(self, pt):
        # test if point is within tet
        test = True
        for n in [0, 1, 2, 3]:
            tmp_mtx = deepcopy(self.jacobian())
            tmp_mtx[1][n] = pt[0]
            tmp_mtx[2][n] = pt[1]
            tmp_mtx[3][n] = pt[2]
            if det(tmp_mtx) < 0:
                test = False
        return test

    def integral(self, numSteps, vtk, rhoQCT=lambda a: a, rhoAsh=lambda b: b, modulus=lambda c: c):
        # calculate integral co-ordinates for the element
        step = 1.0 / numSteps
        volume = 0
        integral = 0
        count = 0
        denint = 0 #IR
        for t in arange(step / 2., 1, step):
            for s in arange(step / 2., 1 - t, step):
                for r in arange(step / 2., 1 - s - t, step):
                    count += 1
                    # find jacobian
                    detJ = det(self.jacobian())
                    # add up volume
                    volume += (detJ / 6.)
                    # calculate shape functions
                    w = [1 - r - s - t,
                         r,
                         s,
                         t]
                    # find co-ordinate for each iteration using shape functions
                    x = [0, 0, 0]
                    for n in [0, 1, 2]:
                        for i in range(4):
                            x[n] += w[i] * self.pts[i][n]
                    # for each co-ordinate, find corresponding CT co-ordinate

                    scales, den, HU = vtk.interpolateScalar(x, modulus, rhoAsh, rhoQCT)
                    integral += (scales) * (detJ / 6.)
                    denint += (den) * (detJ / 6.)
                    # integral += (vtk.interpolateScalar(x, modulus, rhoAsh, rhoQCT)) * (detJ / 6.)

        self.volume = volume / count
        integ = integral / volume
        density = denint / volume
        return integ, density


# Quadratic tetrahedron (C3D10)
class quad_tet:
    """ Quadratic tetrahedral element class """

    __slots__ = ("indx", "pts", "nodes", "volume", "E_integ", "HU_integ","hu_50")

    def __init__(self, indx, pts, nodes):
        #print(indx)
        self.indx = indx
        self.pts = pts
        self.nodes = nodes
        self.volume = []
        self.E_integ = []
        self.HU_integ = []
        self.hu_50 = " "

    def jacobian(self, l, r, s, t):
        # create jacobian for integration methods #IR for volume later
        p = self.pts
        J = [[((p[0][i] * ((4 * l) - 1)) + (p[4][i] * 4 * r) + (p[6][i] * 4 * s) + (p[7][i] * 4 * t)),  # {dN/dl}
              ((p[1][i] * ((4 * r) - 1)) + (p[4][i] * 4 * l) + (p[5][i] * 4 * s) + (p[8][i] * 4 * t)),  # {dN/dr}
              ((p[2][i] * ((4 * s) - 1)) + (p[5][i] * 4 * r) + (p[6][i] * 4 * l) + (p[9][i] * 4 * t)),  # {dN/ds}
              ((p[3][i] * ((4 * t) - 1)) + (p[7][i] * 4 * l) + (p[8][i] * 4 * r) + (p[9][i] * 4 * s))] for i in
             [0, 1, 2]]  # {dN/dt}


        # make the matrix square
        J = [[1., 1., 1., 1.]] + J

        return J

    def integral_hu(self,vtk, rhoQCT=lambda b: b):

        p = self.pts
        step = 1.0 / 4
        volume = 0
        integral = 0
        denint = 0
        huint = 0
        count2 = 0
        count = 0
        for t in arange(step / 2., 1,
                        step):  # local coordinates bzw. Formfunktionen Einheitstetraeder, for numerical integration
            for s in arange(step / 2., 1 - t, step):
                for r in arange(step / 2., 1 - s - t, step):
                    count += 1
                    l = 1 - r - s - t
                    # print("tsrl :", t, s, r, l)
                    # find jacobian
                    detJ = det(self.jacobian(l, r, s, t))
                    # add up volume
                    volume += (detJ / 6.)
                    # calculate shape functions #ANsatzfunktionen des quad. tet elements
                    w = [(2 * l - 1) * l,
                         (2 * r - 1) * r,
                         (2 * s - 1) * s,
                         (2 * t - 1) * t,
                         4 * l * r,
                         4 * r * s,
                         4 * l * s,
                         4 * l * t,
                         4 * r * t,
                         4 * s * t]
                    # find co-ordinate for each iteration using shape functions, x = coordinates of scaled element nodes
                    x = [0, 0, 0]  # zero for adding the values later
                    co = [0, 1, 2]  # for x,y,z from p
                    for n in co:
                        for i in range(10):  # for 10 nodes from tet_quad
                            x[n] += w[i] * p[i][
                                n]  # coordinates is sum of 10 times iteration shapefunction 1 * nodei , all for x,y,z (iterates co)
                            # x of element of one node = summ of his shape function * coordinate of x of node
                            # transformation der Integrationspubkte zu globalen x,y,z koordinaten  ?

                    # print("coordinates x :", x)
                    # print(x)
                    # for each co-ordinate, find corresponding CT co-ordinate# 10 mal durchgeführt für jeden skalierten Knoten
                    hu,dens = vtk.interpolateScalar_hu(x,rhoQCT)  # for 10 tomes numeric iteration, 10 times new box, and HU
                    integral += (hu) * (detJ / 6.)
                    denint += dens *(detJ/6.)
                    #print(detJ)

        self.volume = volume / count
        HU = integral / volume
        density = denint/volume
        #print(HU)
        return HU, density, HU

    def integral(self, vtk, rhoQCT=lambda a: a, rhoAsh=lambda b: b, modulus=lambda c: c):
        # calculate integral co-ordinates for the element
        p = self.pts
        step = 1.0 / 4
        volume = 0
        integral = 0
        denint = 0
        huint = 0
        count2 = 0
        count = 0

        for t in arange(step / 2., 1, step):            #local coordinates bzw. Formfunktionen Einheitstetraeder, for numerical integration
            for s in arange(step / 2., 1 - t, step):
                for r in arange(step / 2., 1 - s - t, step):
                    count += 1
                    l = 1 - r - s - t
                    #print("tsrl :", t, s, r, l)
                    # find jacobian
                    detJ = det(self.jacobian(l, r, s, t))
                    # add up volume
                    volume += (detJ / 6.)
                    # calculate shape functions #ANsatzfunktionen des quad. tet elements
                    w = [(2 * l - 1) * l,
                         (2 * r - 1) * r,
                         (2 * s - 1) * s,
                         (2 * t - 1) * t,
                         4 * l * r,
                         4 * r * s,
                         4 * l * s,
                         4 * l * t,
                         4 * r * t,
                         4 * s * t]
                    # find co-ordinate for each iteration using shape functions, x = coordinates of scaled element nodes
                    x = [0, 0, 0] # zero for adding the values later
                    co = [0, 1, 2]  #for x,y,z from p
                    for n in co:
                        for i in range(10): # for 10 nodes from tet_quad
                            x[n] += w[i] * p[i][n]    # coordinates is sum of 10 times iteration shapefunction 1 * nodei , all for x,y,z (iterates co)
                                                    # x of element of one node = summ of his shape function * coordinate of x of node
                                                    # transformation der Integrationspubkte zu globalen x,y,z koordinaten  ?


                    #print("coordinates x :", x)
                    #print(x)
                    # for each co-ordinate, find corresponding CT co-ordinate# 10 mal durchgeführt für jeden skalierten Knoten
                    #HU = vtk.interpolateScalar_hu(x,e_hu)

                    scales, den, HU = vtk.interpolateScalar(x,modulus, rhoAsh, rhoQCT) # for 10 tomes numeric iteration, 10 times new box, and HU


                    #print("s:h :",count, " :",scales,HU)
                    integral += (scales) * (detJ / 6.) #'e Modul, V3'
                    denint += (den) * (detJ / 6.)
                    huint += (HU) * (detJ/6.)


        #count2 +=1
        #self.volume = volume/count
        density = denint / volume
        integ = integral / volume
        hounsfiled = huint/volume

        return integ, density, hounsfiled

    def integral2(self, vtk, diff, count_2, R, cdb_indx_u50, e_hu=lambda d: d, rhoQCT=lambda a: a, rhoAsh=lambda b: b, modulus=lambda c: c):
        # calculate integral co-ordinates for the element
        p = self.pts
        step = 1.0 / 4
        volume = 0
        integral = 0
        denint = 0
        huint = 0
        count = 0

        for t in arange(step / 2., 1,
                        step):  # local coordinates bzw. Formfunktionen Einheitstetraeder, for numerical integration
            for s in arange(step / 2., 1 - t, step):
                for r in arange(step / 2., 1 - s - t, step):
                    count += 1
                    l = 1 - r - s - t
                    # print("tsrl :", t, s, r, l)
                    # find jacobian
                    detJ = det(self.jacobian(l, r, s, t))
                    # add up volume
                    volume += (detJ / 6.)
                    # calculate shape functions #ANsatzfunktionen des quad. tet elements
                    w = [(2 * l - 1) * l,
                         (2 * r - 1) * r,
                         (2 * s - 1) * s,
                         (2 * t - 1) * t,
                         4 * l * r,
                         4 * r * s,
                         4 * l * s,
                         4 * l * t,
                         4 * r * t,
                         4 * s * t]
                    # find co-ordinate for each iteration using shape functions, x = coordinates of scaled element nodes
                    x = [0, 0, 0]  # zero for adding the values later
                    co = [0, 1, 2]  # for x,y,z from p
                    for n in co:
                        for i in range(10):  # for 10 nodes from tet_quad
                            x[n] += w[i] * p[i][n]  # coordinates is sum of 10 times iteration shapefunction 1 * nodei , all for x,y,z (iterates co)
                            # x of element of one node = summ of his shape function * coordinate of x of node
                            # transformation der Integrationspubkte zu globalen x,y,z koordinaten  ?

                    # print("coordinates x :", x)
                    # print(x)
                    # for each co-ordinate, find corresponding CT co-ordinate# 10 mal durchgeführt für jeden skalierten Knoten
                    # HU = vtk.interpolateScalar_hu(x,e_hu)

                    scales, den, HU = vtk.interpolateScalar2(x, diff, count_2, R,cdb_indx_u50, e_hu, modulus, rhoAsh,
                                                            rhoQCT)  # for 10 tomes numeric iteration, 10 times new box, and HU

                    # print("s:h :",count, " :",scales,HU)
                    integral += (scales) * (detJ / 6.)  # 'e Modul, V3'
                    denint += (den) * (detJ / 6.)
                    huint += (HU) * (detJ / 6.)

        #self.volume = volume / count
        density = denint / volume
        integ = integral / volume
        hounsfiled = huint / volume
        #print(hounsfiled)
        #if diff[count_2] == 0:
        #    print(hounsfiled, integ)
        #for i in cdb_indx_u50:
        #        if count_2 == i:
        #            print(hounsfiled, integ)

        return integ, density, hounsfiled

#Wedge element

class linear_wedge:
    """ Linear wedge element class """

    __slots__ = ("indx", "pts", "nodes", "volume")

    def __init__(self, indx, pts, nodes):
        self.indx = indx
        self.pts = pts
        self.nodes = nodes
        self.volume = []

    def jacobian(self, r, s, t):
        # create jacobian for integration methods
        p = self.pts
        J = [[(p[1][i] - p[0][i]) * (1 - t) + (p[4][i] - p[3][i]) * t,
              (p[2][i] - p[0][i]) * (1 - t) + (p[5][i] - p[3][i]) * t,
              (p[4][i] - p[1][i]) * r + (p[5][i] - p[2][i]) * s + (p[3][i] - p[0][i]) * (1 - r - s)] for i in [0, 1, 2]]

        return J

    def integral(self, numSteps, vtk, rhoQCT=lambda a: a, rhoAsh=lambda b: b, modulus=lambda c: c):
        # calculate integral co-ordinates for the element
        p = self.pts
        step = 1.0 / numSteps
        volume = 0
        integral = 0
        denint = 0
        count = 0
        # iterate through tetrahedral volume, where:
        #    l = natural co-ordinate 1
        #    r = natural co-ordinate 2
        #    s = natural co-ordinate 3
        #    t = natural co-ordinate 4
        for t in arange(step / 2., 1., step):
            for s in arange(step / 2., 1., step):
                for r in arange(step / 2., 1. - s, step):
                    count += 1
                    # find jacobian
                    detJ = det(self.jacobian(r, s, t))
                    # add up volume
                    volume += detJ / 2.
                    # calculate shape functions
                    w = [(1 - r - s) * (1 - t),
                         s * (1 - t),
                         r * (1 - t),
                         (1 - r - s) * t,
                         s * t,
                         r * t]

                    # find co-ordinate for each iteration using shape functions
                    x = [0, 0, 0]
                    for n in [0, 1, 2]:
                        for i in range(6):
                            x[n] += w[i] * p[i][n]
                    # for each co-ordinate, find corresponding CT co-ordinate
                    scales, den = vtk.interpolateScalar(x, modulus, rhoAsh, rhoQCT)
                    integral += scales * (detJ / 2.)
                    # integral += (vtk.interpolateScalar(x, modulus, rhoAsh, rhoQCT)) * (detJ / 2.)
                    denint += (den) * (detJ / 2.)

        self.volume = volume / count
        density = denint / volume
        integ = integral / volume

        return integ, density

    # Hexahedral element


class linear_hex:
    """ Hexahedral element class """

    # __slots__ = ("indx", "pts", "nodes", "jacobian")

    def __init__(self, indx, pts, nodes):
        self.indx = indx
        self.pts = pts
        self.nodes = nodes
        self.volume = []

    def jacobian(self, r, s, t):
        # create jacobian for integration methods
        p = self.pts
        J = [[(p[1][i] - p[0][i]) * (1 - s) * (1 - t) + (p[2][i] - p[3][i]) * (1 + s) * (1 - t) + (
                    p[5][i] - p[4][i]) * (1 - s) * (1 + t) + (p[6][i] - p[7][i]) * (1 + s) * (1 + t),
              (p[3][i] - p[0][i]) * (1 - r) * (1 - t) + (p[2][i] - p[1][i]) * (1 + r) * (1 - t) + (
                          p[7][i] - p[4][i]) * (1 - r) * (1 + t) + (p[6][i] - p[5][i]) * (1 + r) * (1 + t),
              (p[4][i] - p[0][i]) * (1 - r) * (1 - s) + (p[5][i] - p[1][i]) * (1 + r) * (1 - s) + (
                          p[6][i] - p[2][i]) * (1 + r) * (1 + s) + (p[7][i] - p[3][i]) * (1 - r) * (1 + s)] for i in
             [0, 1, 2]]

        return J

    def integral(self, vtk, rhoQCT=lambda a: a, rhoAsh=lambda b: b, modulus=lambda c: c):
        # calculate integral co-ordinates for the element
        p = self.pts
        step = 1.0 / 2 #change from num step
        volume = 0
        integral = 0
        denint = 0
        count = 0
        for t in arange(-1 + step, 1, step * 2):
            for s in arange(-1 + step, 1, step * 2):
                for r in arange(-1 + step, 1, step * 2):
                    count += 1
                    # find jacobian
                    detJ = det(self.jacobian(r, s, t))
                    # add up volume
                    volume += (detJ / 8.)
                    # calculate shape functions
                    w = [((1 - r) * (1 - s) * (1 - t)) / 8.,
                         ((1 + r) * (1 - s) * (1 - t)) / 8.,
                         ((1 + r) * (1 + s) * (1 - t)) / 8.,
                         ((1 - r) * (1 + s) * (1 - t)) / 8.,
                         ((1 - r) * (1 - s) * (1 + t)) / 8.,
                         ((1 + r) * (1 - s) * (1 + t)) / 8.,
                         ((1 + r) * (1 + s) * (1 + t)) / 8.,
                         ((1 - r) * (1 + s) * (1 + t)) / 8., ]
                    # find co-ordinate for each iteration using shape functions
                    x = [0, 0, 0]
                    for n in [0, 1, 2]:
                        for i in range(8):
                            x[n] += w[i] * p[i][n]
                    # for each co-ordinate, find corresponding CT co-ordinate
                    scales, den = vtk.interpolateScalar(x, modulus, rhoAsh, rhoQCT)
                    integral += scales * (detJ / 8.)
                    denint += den * (detJ / 8.)

                    # integral += (vtk.interpolateScalar(x, modulus, rhoAsh, rhoQCT)) * (detJ / 8.)

        self.volume = volume / (count * 8)

        density = denint / volume
        integ = integral / volume

        return integ, density


#Test File

#vtk = import_ct_data(r"C:\Users\8iries\Desktop\Praktikum\Programm\CT_dataset")
#param = import_parameters_xml(r"C:\Users\8iries\Desktop\Praktikum\Programm\params.xml")
#parts = import_mesh(r"C:\Users\8iries\Desktop\Praktikum\Programm\mesh_CT-dataset2.cdb", param)

