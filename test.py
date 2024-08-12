from run import run
#For Simple: Paramaters File | CT-Data | Nodes | Elements
#For Realistic: Paramaters File | CT-Data | Mesh

def run_programm_from_test(gui_calc_params):
    print(gui_calc_params)
    """simple test files"""
    #True or None
    seperate = True


    with open (r'pathcdb.txt', 'r') as M:
        mesh = M.read()
    #Mesh = open(r'pathcdb.txt', 'r')  # get path from txt from GUI input
    #mesh = Mesh.read()
    print(mesh)
    param = "p"

    with open (r'pathct.txt', 'r') as P:
        ct = P.read()
    #CT = open(r'pathct.txt', 'r')
    #ct = CT.read()
    print(ct)
    run(param, ct, mesh, gui_calc_params,seperate)

'''run programm without GUI'''
#g = gui_calc_params , the params you would set in the gui
#g = {'ansys_version': '', 'treshhold_check': 'single', 'a': 0.87719, 'b': 0.07895, 'c_cort': 14664.0, 'd_cort': 1.49, 'treshhold_value': 0.5, 'c_trab': 14664.0, 'd_trab': 1.49, 'gap': 200.0, 'poisson': 0.3}
#run_programm_from_test(g)

#pylibjpeg
#gdcm
#pylibjpeg-libjpeg


#param = r"C:\Users\8iries\Desktop\Praktikum\Programm\params.xml"
#ct = r"C:\Users\8iries\Desktop\Praktikum\Programm\CT_dataset"
#mesh = r"C:\Users\8iries\Desktop\Praktikum\Programm\ct2.cdb"

""" own pelvis BONE"""
#param = r"C:\Users\8iries\Desktop\Praktikum\Programm\params.xml"
#ct = r"C:\Users\8iries\Desktop\Praktikum\Programm\Patienten\002_Patient\Bild-Daten\CT\EBDD08BC\9333C43C"
#mesh = r"C:\Users\8iries\Desktop\Praktikum\Programm\Patienten\002_Patient\Patient_02_without_Contacts.cdb"


"""Bonemat_Testfile_Realistic"""
#ct = r"C:\Users\8iries\Desktop\Praktikum\Programm\Realistic\CT_dataset"
#param = r"C:\Users\8iries\Desktop\Praktikum\Test_Programm\params.xml"
#mesh = r"C:\Users\8iries\Desktop\Praktikum\Programm\Realistic\Mesh.cdb"

#run('Ansys_Testing_Files/Simple/Configuration_file.xml','Ansys_Testing_Files/Simple/Volume.vtk','Ansys_Testing_Files/Simple/Nodes.lis','Ansys_Testing_Files/Simple/Elements.lis')
#run(param,ct,mesh)
#run('Abaqus Example Files/example_parameters.txt','Abaqus Example Files/example_ct_data.vtk','Abaqus Example Files/example_abaqus_mesh.inp')




