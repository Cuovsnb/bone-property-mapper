#!/usr/bin/python
#
# py_bonemat_ansys - run
# ==========================
#
# Created by Abdul Hanan

__all__ = ['run']

# -------------------------------------------------------------------------------
# Import modules
# -------------------------------------------------------------------------------
import sys
import os
import general, data_import, calc, data_output
from version import __version__
import time
from datetime import date


# -------------------------------------------------------------------------------
# Define run program
# -------------------------------------------------------------------------------
def run(argv0, argv1, argv2, gui_params, seperate, argv3=None):
    # ---------------------------------------------------------------------------
    # create welcome screen
    # ---------------------------------------------------------------------------
    global time_2
    t = time.time()
    d = date.today()
    month = d.strftime("%B")[:3]
    year = d.isocalendar()[0]
    print(("""
    ************** MaterailMapping """ + __version__ + """ ************
    *** Immanuel Ries, Fraunhofer IMWS,   """ + month + """ """ + repr(year) + """ ***
    **************************************************
    """))
    # ---------------------------------------------------------------------------
    # Check input arguments: Simple OR Realistic
    # ---------------------------------------------------------------------------
    if argv3 is not None:
        print('    Image Dataset Type: Simple')
        argv = [argv0, argv1, argv2, argv3]
        if not general.check_argv_second(argv):
            sys.exit(1)
        param_file = argv[0]
        ct_data = argv[1]
        nodes = argv[2]
        elements = argv[3]

        # ---------------------------------------------------------------------------
        # Import parameters
        # ---------------------------------------------------------------------------
        print(("    Importing parameters file: " + param_file))
        if ".txt" in param_file:
            param = data_import.import_parameters(param_file)
        elif ".xml" in param_file:
            param = data_import.import_parameters_xml(param_file)
        # ---------------------------------------------------------------------------
        # Import Abaqus input file data
        # ---------------------------------------------------------------------------
        #print(("    Importing mesh files: " + elements + " and " + nodes))
        parts = data_import.import_ansys_simple(elements, nodes)
        #print(parts)
        # ---------------------------------------------------------------------------
        # Import CT data
        # ---------------------------------------------------------------------------
        #print(("    Importing CT data: " + ct_data))
        vtk_data = data_import.import_ct_data(ct_data)
        # ---------------------------------------------------------------------------
        # Determine material properties for elements within each part
        # ---------------------------------------------------------------------------
        # print("    Calculating material properties")
        # parts = calc.calc_mat_props(parts, param, vtk_data)

        # #---------------------------------------------------------------------------
        # # Write data to new abaqus input file
        # #---------------------------------------------------------------------------
        # print("    Writing material data to new abaqus input file:")
        # print(("\t" + mesh_file[:-4] + "MAT.inp"))
        # data_output.output_abq_inp(parts, mesh_file, param['poisson'])
        # #---------------------------------------------------------------------------
        # # End
        # #---------------------------------------------------------------------------
        # print("""
        # **   !!! Bone material assignment complete !!!   **
        # ***************************************************
        # """)
        # tt = time.time() - t
        # print((" Elapsed time: " + repr(tt)))
        # os._exit(0)

    else:
        print('    Image Dataset Type: Realistic')
        argv = [argv0, argv1, argv2, gui_params]
        print(argv)
        if not general.check_argv(argv):
            sys.exit(1)
        param_file = argv[0]
        ct_data = argv[1]
        mesh_file = argv[2] #path to files
        gui_param = argv [3]


        # ---------------------------------------------------------------------------
        # Import parameters
        # ---------------------------------------------------------------------------
        print(("    Importing parameters file: " + param_file))
        if ".txt" in param_file:
            param = data_import.import_parameters(param_file)
        #elif ".xml" in param_file:
        else:
            param = data_import.import_parameters_xml(param_file, gui_param)
        # ---------------------------------------------------------------------------
        # Import Abaqus input file data ##IR: or Ansys Mesh
        # ---------------------------------------------------------------------------
        print(("    Importing mesh file: " + mesh_file))
        parts = data_import.import_mesh(mesh_file, param)
        #print("Parts1 :", parts)

        # ---------------------------------------------------------------------------
        # Import CT data ##vtk verwirrend, eigentlich ct data importiert
        # ---------------------------------------------------------------------------
        print(("    Importing CT data: " + ct_data))
        vtk_data = data_import.import_ct_data(ct_data)
        # ---------------------------------------------------------------------------
        # Determine material properties for elements within each part
        # ---------------------------------------------------------------------------
        print("    Calculating material properties")
        parts,indi = calc.calc_mat_props(parts, param, vtk_data,seperate) #parts immernoch in Lise

        modulinew = []
        #for i in indi[0]: #indi fängt bei 0 an
        #    modulinew.append(parts[0].moduli2[i])
        print(modulinew)

        write_data(vtk_data.HU,vtk_data.E,vtk_data.E_int)
        # e und e int zwischenschritte Werte (10 Pro Element bzw. mehhr beie_intp)

        # #---------------------------------------------------------------------------
        # # Write data to new abaqus input file
        # #---------------------------------------------------------------------------
        print("    Writing material data to new CDB input file:")
        print(("\t" + "ANSYS_Assigned.cdb"))
        # data_output_ansys.output_ansys_inp(parts, mesh_file, param['poisson'],AVG)
        print(gui_param)
        data_output.output_ansys_inp(parts,mesh_file, param['poisson'], gui_param["ansys_version"],seperate,indi)
        # #---------------------------------------------------------------------------
        # # End
        # #---------------------------------------------------------------------------
        print("""
        **   !!! Materials assigned successfuly !!!   **
        ***************************************************
        """)
        tt = time.time() - t
        time_2 = str(round(float(tt),3))
        print((" Elapsed time: " + repr(tt)))
        print(time_2)
        with open("time.txt", "w") as t:
            t.write(time_2)

        #os._exit(0)

def write_data(Hu, E, E_int):

    os.chdir(sys.path[0])

    #hu_sorted = sorted(Hu)
    #E_sorted = sorted(E)
    #Eint_sorted = sorted(E_int)

    with open("HU.txt", "w") as f:
        for i in Hu:
            f.write(str(i))
            f.write("\n")

    with open("E.txt", "w") as p:
        for j in E:
            p.write(str(j))
            p.write("\n")

    with open ("E_int", "w") as o:
        for k in E_int:
            o.write(str(k))
            o.write("\n")


