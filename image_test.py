import os
import sys

import SimpleITK as sitk
import pyvista as pv
import numpy as np
from src.threeD.loaddicomfile import load_dcm_info
import pydicom
import ansys.mapdl.reader as reader


#my_theme = themes.DarkTheme
#my_theme.color = "white"
#np.set_printoptions(threshold=sys.maxsize)
#path = r'G:\Ries\Backup\Praktikum\Programm\Patienten\002_Patient\Bild-Daten\CT\EBDD08BC\9333C43C'


def main():




    #p = pv.Plotter(shape=(1,2))
    #p.set_background("black")


    #path to ct
    os.chdir(sys.path[0])
    with open('pathct.txt', 'r') as f:
        path = f.read()

    #read dcm pictures and sort them with z value
    pic = []
    imgs = os.listdir(path)
    def dicom(dicom_files,path):
        z = [pydicom.dcmread(os.path.join(path,d), force=True).ImagePositionPatient[2] for d in dicom_files]
        dicom_order = list(zip(*sorted(zip(z, dicom_files))))[1]
        return dicom_order
    sorted_dicom = dicom(imgs, path)

    #read each dicom and pixel data and save them in list
    pixel_list = []
    for i in sorted_dicom:
        v = pydicom.dcmread(path + "/" +i)
        pic.append(v)

    for j in range(len(pic)):
        v = pic[j].pixel_array
        pixel_list.append(v)

    num_slices = len(pixel_list)

    #get the slice thickeness from z value diefferenz
    ss = abs(round(float(pic[0].ImagePositionPatient[2])-float(pic[1].ImagePositionPatient[2]),4))
    pixel_array = np.array(pixel_list)
    pixel_array_T = pixel_array.T

    #create uniform grid for CT Data, gave them cell geometry of CT
    grid = pv.UniformGrid()

    x = pic[0].ImagePositionPatient[0]
    y = pic[0].ImagePositionPatient[1]
    z = pic[0].ImagePositionPatient[2]
    slice_info = [z, ss, num_slices]

    grid.origin = (float(x), float(y), float(z))
    grid.spacing = (float(pic[0].PixelSpacing[0]),float(pic[0].PixelSpacing[0]), ss )
    grid.dimensions = np.array(pixel_array_T.shape)+1
    grid.cell_data["HU"] = pixel_array_T.flatten(order="F")#

    #create treshhold mesh, ct data for z_slice and volume rendering
    ct_treshhold = grid.threshold(value=(0,1000))
    ct_treshold_2 = grid.threshold(value=(-500,2000))
    #create orthogonal slices
    #slices = grid.slice_orthogonal(x=0,y=0,z=0)
    #dargs = dict(cmap = "gist_gray")


    #cdb file
    with open('pathcdb.txt', 'r') as f:
        path_cdb = f.read()
    archive_cdb = reader.Archive(path_cdb)  # storage cdb as unstructured Grid
    grid_cdb = archive_cdb._parse_vtk(force_linear=True)
    #cdb_grid_slice = grid_cdb.slice_orthogonal(x=0, y=0, z=0)

    #coordinate axis for plotting zero point
    X = np.array([[0,0,0], [200,0,0]])
    Y = np.array([[0,0,0],[0,200,0]])
    Z = np.array([[0,0,0], [0,0,200]])

    #p.subplot(0,0)
    #p.add_lines(X, color = "red", width=10)
    #p.add_lines(Y, color = "yellow", width=10)
    #p.add_lines(Z, color = "green", width=10)
    #add mesh to plotter
    #p.add_mesh_threshold(ct_treshhold, cmap = "gist_gray", opacity = 0.9)
    #p.add_mesh(grid_cdb, color="red", opacity = 0.9)

    #p.subplot(0,1)
    #p.add_mesh_slice(grid_cdb, style = "wireframe", generate_triangles = True, assign_to_axis = 0, normal_rotation= False)
    #p.add_mesh_slice(grid, assign_to_axis = 0, normal_rotation=False)
    #adding slices

    #p.add_mesh(slices, cmap="gist_gray")
    #p.add_mesh(cdb_grid_slice, color = "red", style ="wireframe")
    #p.add_camera_orientation_widget()

    #p.add_axes(box=True)
    #p.add_axes_at_origin()
    #p.add_bounding_box(color='grey', corner_factor=0.5, line_width=True, opacity=1.0)
    #p.show()
    return ct_treshhold,X,Y,Z,grid_cdb,ct_treshold_2, slice_info



if __name__ == '__main__':
    main()