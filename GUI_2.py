import sys
# Setting the Qt bindings for QtPy
import os

#os.environ["QT_API"] = "pyqt5"

import test
from qtpy.QtWidgets import QDialog, QApplication, QFileDialog, QTreeWidgetItem
import xml.etree.ElementTree as et
from qtpy import QtWidgets
from qtpy import uic
import ansys.mapdl.reader as reader
import pyvista as pv
from pyvistaqt import QtInteractor, MainWindow
import image_test
import praph
from src.threeD import threeD_module
import numpy as np
import pydicom
from src.threeD.loaddicomfile import load_dcm_info



class MyMainWindow(MainWindow):

    def __init__(self, parent=None, show=True):
        QtWidgets.QMainWindow.__init__(self, parent)
        uic.loadUi("mainwindow.ui", self)
        # create the frame
        # setzt neues V Layout um pyvista Widget in frame anzupassen
        vlayout = QtWidgets.QVBoxLayout()
        # add the pyvista interactor object
        # setzt plot in frame widget
        self.plotter = QtInteractor(self.frame)
        self.plotter_2 = QtInteractor(self.frame_2)

        # setzt in layout das widget
        vlayout.addWidget(self.plotter)
        vlayout.addWidget(self.plotter_2)
        # setzt in frame das layout mit dem pyvista plot
        self.frame.setLayout(vlayout)
        self.frame_2.setLayout(vlayout)
        self.pathct = []
        self.pathcdb = []
        self.pathnew = []
        self.calibration_decicion = None

        #für neues QWindow  für CT_slice
        self.refresh3D = True
        self.ui_3D = None
        self.justclose = False
        self.maindirectory = os.getcwd()
        self.first3D = True

        # interaktionen mit Buttons
        self.pushcdb.clicked.connect(self.browsefiles_cdb)
        self.showcdb.clicked.connect(self.plot)
        self.pushdic.clicked.connect(self.browsefiles_ct)
        self.pushnew.clicked.connect(self.browsefiles_new)
        self.showdic.clicked.connect(self.showdicom)
        self.showcdblonly.clicked.connect(self.showcdbonly)
        self.pushsto.clicked.connect(self.run)
        self.calibration.clicked.connect(self.cali)
        self.radioButton.clicked[bool].connect(self.grid_in_slice)

        #self.show()

        if show:
            self.show()

    def cali(self):
        self.calib = praph.Calibration()
        self.calib.show()

    def browsefiles_cdb(self):
        global path_cdb
        path_cdb = QFileDialog.getOpenFileName(self, "Open CDB", "C:\Benutzer", "CDB File(*.cdb)")
        self.linecdb.setText(path_cdb[0])
        self.pathcdb = self.linecdb.text()
        with open('pathcdb.txt', 'w') as f:
            for line in self.pathcdb:
                f.write(line)

    def browsefiles_ct(self):
        global path_ct
        path_ct = QFileDialog.getExistingDirectory(self, "Open CT-File", "C:\Benutzer")
        self.linedic.setText(path_ct)
        self.pathct = self.linedic.text()
        print(sys.path)
        with open('pathct.txt', 'w') as f:
            for line in self.pathct:
                f.write(line)

    def browsefiles_new(self):

        global path_new
        path_new = QFileDialog.getExistingDirectory(self, "Search Storage", "C:\Benutzer")
        name = self.linenew.displayText()
        self.linesto.setText(path_new)
        self.pathnew = self.linesto.displayText()

        with open ('pathnew.txt', 'w') as f:
            for line in self.pathnew:
                f.write(line)

        with open ('name.txt', 'w') as i:
            for line in name:
                i.write(line)

    def create_mesh_z(self, value):
        """

        :type value: object
        """
        res = int(value)
        sphere = self.grid_ct.slice(normal="z", origin=(0, 0, res))
        mesh = self.grid_cdb.slice(normal="z", origin=(0, 0, res))
        self.plotter_2.add_mesh(sphere, name="sphere", cmap="bone")
        self.plotter_2.add_mesh(mesh, name="mesh", color="red", style="wireframe")
        return

    def grid_in_slice(self, down):
        if down:
            self.pic = self.plotter_2.add_mesh(self.grid_cdb, color="red", style="wireframe")
        else:
            self.plotter_2.remove_actor(self.pic)

    def plot(self):
        '''plotting 3d ct and cdb'''


        self.plotter.set_background("black")
        self.plotter_2.set_background("black")

        #get ct and cdb data for 3d plot
        self.dataset_ct, x,y,z,self.grid_cdb,self.grid_ct, slice_info = image_test.main() #slice_info 0 [startpunkt, slice,thickness, anzahl slices]
        #print(path_ct)



        #add volume and mesh to the plotter widget

        self.plotter.add_mesh_threshold(self.dataset_ct, cmap = "gist_gray", opacity = 1 )
        #self.plotter.add_volume(ct_data, cmap="gist_gray", opacity= "geom", opacity_unit_distance=28) # ad to the plotter widget the volume etc...
        self.plotter.add_mesh(self.grid_cdb, color="red", opacity = 1)
        self.plotter.add_lines(x, color="red", width=10)
        self.plotter.add_lines(y, color="yellow", width=10)
        self.plotter.add_lines(z, color="green", width=10)
        self.plotter.add_camera_orientation_widget()
        #self.plotter.add_axes(box=True)
        #self.plotter.add_axes_at_origin()
        #point = np.array([0.0, 0.0, 0.0])
        #labels = ['Point 0, 0, 0']
        #self.plotter.add_point_labels(point, labels, show_points=True, point_size=20, point_color="r", render_points_as_spheres=True, font_size=22)



        self.plotter_2.add_slider_widget(
            self.create_mesh_z,
            [int(slice_info[0]), int(slice_info[0] + (slice_info[1]*slice_info[2]))],
            title="Z",
            title_opacity=0.5,
            title_color="green",
            pointa=(0.03, 0.9),
            pointb=(0.4, 0.9),
            # fmt="%0.9f",
            title_height=0.05,
            color = "white",
        )
        #if self.check_grid.isChecked == True:
            #self.plotter_2.add_mesh(self.grid_cdb, color="red", style="wireframe")
        #self.plotter_2.add_mesh(self.slices_ct, cmap="gist_gray")
        #self.plotter_2.add_mesh(self.slices_cdb,  color = "red", style ="wireframe")
        self.plotter_2.add_bounding_box(color='grey', corner_factor=0.5, line_width=True, opacity=1.0)

        self.plotter_2.show()
        self.plotter.show()

    def showdicom(self, path):

        os.chdir(self.maindirectory)
        if self.refresh3D:
            if self.first3D:
                self.ui_3D = threeD_module.CthreeD()
                self.first3D = False
            self.ui_3D.show()
            self.refresh3D = False
            self.ui_3D.rejected.connect(self.closeCT)
        elif self.justclose:
            pass
        else:
            if not self.threeDaction.isChecked():
                self.threeDaction.setChecked(True)
            self.ui_3D.raise_()

        #self.ui_3D = threeD_module.CthreeD()
        #self.ui_3D.show()
        #self.ui_3D.rejected.connect(self.close3d)

    def showcdbonly(self):
        filename = path_cdb[0]
        archive = reader.Archive(filename)
        grid = archive._parse_vtk(force_linear=True)
        grid.plot(color='r', show_edges=True)

        #self.okButton.clicked.connect(self.accept)
        #self.cancelButton.clicked.connect(self.reject)

    def closeCT(self):
        os.chdir(self.maindirectory)
        self.justclose = True
        #self.threeDaction.setChecked(False)
        self.ui_3D.__init__()
        self.refresh3D = True
        self.justclose = False

    def params_gui(self):
        gui_calc_params = {}

        if self.checkBox_192.isChecked() == True:
            gui_calc_params["ansys_version"] = "19.2"
        elif self.checkBox_21r1.isChecked() == True:
            gui_calc_params["ansys_version"] = "2021 R1"
        elif self.checkBox_22r2.isChecked() == True:
            gui_calc_params["ansys_version"] = "2022 R2"
        else:
            gui_calc_params["ansys_version"] = (str(self.line_ansys.displayText()))

        if self.checkBox_tc.isChecked == True:
            gui_calc_params["treshhold_check"] = 'double'
        else:
            gui_calc_params["treshhold_check"] = 'single'

        gui_calc_params["a"] = float(self.lineEdit_a.displayText())
        gui_calc_params["b"] = float(self.lineEdit_b.displayText())
        gui_calc_params["c_cort"] = float(self.lineEdit_cc.displayText())
        gui_calc_params["d_cort"] = float(self.lineEdit_dc.displayText())
        gui_calc_params["treshhold_value"] = float(self.lineEdit_th.displayText())
        gui_calc_params["c_trab"] = float(self.lineEdit_ct.displayText())
        gui_calc_params["d_trab"] = float(self.lineEdit_dt.displayText())
        gui_calc_params["gap"] = float(self.lineEdit_gap.displayText())
        gui_calc_params["poisson"] = float(self.lineEdit_p.displayText())
        print(gui_calc_params)
        return gui_calc_params

    def run(self):

        # get value into running grogramm, if calibration not clicked, it runs with default calibration
        if self.checkBox.isChecked() != True:
            print("Pfad: Starttt HU ?")
            print(os.getcwd())
            with open('HU_phantoms.txt', 'w') as f:
                #for line in str(self.hu):
                #    f.write(line)
                #f.write('\n')
                #for i in str(rho_data1):
                #    f.write(i)
                f.write("None")

        gui_calc_params = self.params_gui()
        value = 1
        self.show_popup(value)

        test.run_programm_from_test(gui_calc_params)

        value = 0
        self.show_popup(value)

        #main.run()

    def show_popup(self, value):
        if value == 1:
            QtWidgets.QMessageBox.information(None, "Info", "Start Calculation", QtWidgets.QMessageBox.Ok)
        else:
            with open("time.txt" , "r") as i:
                dauer = i.read()
            QtWidgets.QMessageBox.information(None, "Info", "Finish Calculation with time: " + dauer +"s = " +
                                              str(round((float(dauer)/60),3)) + "min = " +
                                              str(round((float(dauer)/(60*60)),3)) + "h " + ", in " +
                                              str(self.pathnew), QtWidgets.QMessageBox.Ok)

def main():
    app = QtWidgets.QApplication(sys.argv)
    window = MyMainWindow()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()



#app = QtWidgets.QApplication(sys.argv)
#window = MainWindow()
#try:
    #app.exec()
#except AttributeError:
    #app.exec_()
#window.show(