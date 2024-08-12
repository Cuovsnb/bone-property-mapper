import os
import sys
import numpy as np
from PyQt5 import QtWidgets
import pyqtgraph as pg
import pydicom as pyd
import itertools
import matplotlib.pyplot as plt


class Calibration(QtWidgets.QMainWindow):

    def __init__(self):
        super().__init__()

        # Basic UI layout

        self.statusbar = QtWidgets.QStatusBar(self)
        self.setStatusBar(self.statusbar)
        self.glw = pg.GraphicsLayoutWidget()
        self.setCentralWidget(self.glw)
        self.pushbutton = QtWidgets.QPushButton(self)
        self.pushbutton2 = QtWidgets.QPushButton(self)
        self.rho = QtWidgets.QLineEdit(self)
        self.func = QtWidgets.QLineEdit(self)
        self.slice = None
        self.path = None
        self.x = [] # coordinates picked
        self.y = [] # coordinates picked
        self.order = []
        self.hu = []
        self.hu_values = []


        #set UI
        self.pushbutton2.move(100,0)
        self.pushbutton.setText("push circle")
        self.pushbutton2.setText("save HU")
        self.rho.move(201,0)
        self.rho.setText("BMD in g/cm^3")
        self.func.move(302,0)
        self.func.setFixedSize(300,30)

        # Make image plot
        self.p1 = self.glw.addPlot()
        self.p1.getViewBox().setAspectLocked()

        # Draw axes and ticks above image/data
        [ self.p1.getAxis(ax).setZValue(10) for ax in self.p1.axes ]
        self.data = self.dicom()

        #self.data = np.random.rand(120, 100)
        self.img = pg.ImageItem(self.data)
        self.p1.addItem(self.img)
        self.circle = pg.TargetItem(size=70, symbol="o")
        self.circle.getViewBox()
        self.p1.addItem(self.circle)
        # Centre axis ticks on pixel
        self.img.setPos(0.5, 0.5)

        # Swap commented lines to choose between hover or click events
        self.p1.scene().sigMouseMoved.connect(self.mouseMovedEvent)
        #self.p1.scene().sigMouseClicked.connect(self.mouseClickedEvent)
        self.pushbutton.clicked.connect(self.pixel)
        self.pushbutton2.clicked.connect(self.save)


    def save(self):

        # linear regression
        # get data from input
        rho_data = self.rho.displayText()
        rho_data1 = rho_data.split(",")
        rho_data2 = [float(slot) for slot in rho_data1]  # cast to float in List

        self.rho_array = np.array(rho_data2)
        self.hu_array = np.array(self.hu)


        if len(self.hu) == len(rho_data1):

            #get coeffitientren
            b = self.coef(self.hu_array, self.rho_array) # b_o = n, b_1 = x for mx +n
            print("Estimated coefficients:\nb_0 = {}  \
                      \nb_1 = {}".format(b[0], b[1]))

            #plot line and points
            self.plot_line(self.hu_array, self.rho_array, b)

            self.func.setText("y = " + str(round(b[0],3)) + "+ m * " + str(round(b[1],3)) + "  R^2 = " + str(self.r))


            #b = list(b)
            b = str(b[0]) + ";" + str(b[1])
            with open('HU_phantoms.txt', 'w') as f:
                #for line in str(self.hu):
                #    f.write(line)
                #f.write('\n')
                #for i in str(rho_data1):
                #    f.write(i)
                #f.write(str(b[0]) + "," + str(b[1])) b[1] = m
                f.write(b)
            print(self.hu_array)
            with open("hu_data_praph.txt", 'w') as hu:

                for i in self.hu_values:
                    hu.write(str(i)+"\n")



            #b[1] = m #probieren nur test , muss später weg, geth schneller als un GUI_2
            with open("HU_phantoms.txt", "r") as f:
                objekt = f.read()
                print(objekt)
                if objekt != "None":
                    coef = objekt.split(";")
                    coef[1] = float(coef[1])
                    coef[0] = float(coef[0])

                    print(coef)
                    result = coef[0] + (coef[1] * 300)
                    if result < 0:
                        return 0.000001


                    print(result)

    def plot_line(self,x,y,b):
        # plotting the actual points as scatter plot
        plt.scatter(x, y, color="m",
                    marker="o", s=30)

        # predicted response vector
        y_pred = b[0] + b[1] * x

        # plotting the regression line
        plt.plot(x, y_pred, color="g" )

        # putting labels
        plt.xlabel('HU Value')
        plt.ylabel('Density')


        # function to show plot
        plt.show()

    def coef(self,x,y):

        n = np.size(x)

        # mean of x and y vector
        m_x = np.mean(x)
        m_y = np.mean(y)

        # calculating cross-deviation and deviation about x
        SS_xy = np.sum(y * x) - n * m_y * m_x
        SS_xx = np.sum(x * x) - n * m_x * m_x

        # calculating regression coefficients
        b_1 = SS_xy / SS_xx
        b_0 = m_y - b_1 * m_x

        #R^2 = erklärte/ gesamte Abweichung

        y_mean = (sum(y) / (len(y)))
        y_reg = [b_0 + (b_1 * i) for i in x]
        erkl = sum([(i - y_mean) ** 2 for i in y_reg])
        ges = sum([(i - y_mean) ** 2 for i in y])


        self.r = round(erkl/ges, 3)

        return (b_0, b_1)

    def pixel(self):

        #middlepoint of circle
        point = self.circle.pos()

        #picked pixel (amount of pixel side)
        self.a = 7

        #start, end of Fläche
        x_start = int(point[0] - self.a +1)
        x_end = int(point[0] + self.a +1)
        y_start = int(point[1] - self.a +1)
        y_end = int(point[1] + self.a +1)


        #range of x,y coordinates
        x = list(range(x_start, x_end, 1))
        y = list(range(y_start, y_end,1))

        #coordinates of fläche, erzeugt alle koordinatenpunkte von listen x,y
        self.points = list(itertools.product(x,y))

        #coordinates in own list
        X=[]
        Y=[]
        for i in self.points:
            X.append(i[0])
            Y.append(i[1])

        #plot coordinates
        line = pg.ScatterPlotItem()
        line.addPoints(X, Y)
        self.p1.addItem(line)

        #index of slices for HU
        slice_number = []
        for i in [self.slice, self.slice*2, self.slice*3, self.slice*4, self.slice*5, self.slice*6]:
            slice_number.append(i)

        #print(self.points)
        #new_points = set(self.points)
        #print(new_points)
        #self.points = list(new_points)
        #print(self.points)

        #get pixel value of all 6 slices, value_arra = array with all pixel values from all slices
        for i in slice_number:
            plan = pyd.read_file(self.path + "/" + self.order[i])
            image_2d = plan.pixel_array.astype(float)
            array = np.array(image_2d)  # create array
            array2 = np.rot90(array, k=3)
            value_array = []
            for l in range(len(self.points)):
                value_array.append(int(array2[self.points[l][0]][self.points[l][1]]))




        #value of array to HU
        hu_array = []
        for i in value_array:
            j = (i * int(self.slope)) + int(self.intercept)
            hu_array.append(j)
            Hu_array = sorted(hu_array)

        self.hu_values.append(Hu_array)



        #median hu score from each circle
        hu_score = sum(hu_array)/len(hu_array)
        print(hu_score)

        #add to list
        self.hu.append(hu_score)
        print(self.hu)

    def dicom(self):
        #get path
        os.chdir(sys.path[0])
        objekt = open(r'pathct.txt', 'r')  # QFileDialog.getExistingDirectory(self, 'choose dicom directory')
        self.path = objekt.read()  # path tot ct

        #sort dycoms by Z coordinate
        dicom_files = os.listdir((self.path))
        z = [pyd.dcmread(os.path.join(self.path, d), force=True).ImagePositionPatient[2] for d in dicom_files]
        dicom_order = list(zip(*sorted(zip(z, dicom_files))))[1]
        self.order = dicom_order

        #get intercept and slope for HU
        dcm_data = pyd.dcmread(self.path + "/" + dicom_order[0])
        self.slope = str(dcm_data.RescaleSlope) # most 1
        self.intercept = str(dcm_data.RescaleIntercept)

        #max and slice nuber (3 slices for pixel HU)
        max = len(os.listdir(self.path))
        num = int(max * 0.5)
        self.slice = int(max/7)

        #read dicom pixel data
        plan = pyd.read_file(self.path + "/" + self.order[num])# picture in middle of dataset, num
        image_2d = plan.pixel_array.astype(float)
        array = np.array(image_2d) #create array
        array2 = np.rot90(array, k = 3) # rotate picture

        return array2

    def mouseClickedEvent(self, event):
        self.mouseMovedEvent(event.pos())
        print(self.circle.pos())

    def mouseMovedEvent(self, pos):
        # Check if event is inside image, and convert from screen/pixels to image xy indicies
        if self.p1.sceneBoundingRect().contains(pos):
            mousePoint = self.p1.getViewBox().mapSceneToView(pos)
            x_i = round(mousePoint.x())
            y_i = round(mousePoint.y())

            #if x_i <= 512 and x_i >=0 and y_i <= 512 and y_i >=0 :
            #    self.x.append(x_i)
            #    self.y.append(y_i)
            #    self.points.append(tuple([x_i, y_i]))

            #line = pg.ScatterPlotItem()
            #line.addPoints(self.x, self.y)
            #self.p1.addItem(line)



            if x_i >= 0 and x_i <= self.data.shape[0] and y_i >= 0 and y_i <= self.data.shape[1]:
                self.statusbar.showMessage("({}, {}) = {:0.2f}".format(x_i, y_i, self.data[x_i -1, y_i-1]))
                return

        self.statusbar.clearMessage()

def main():
    import sys
    app = QtWidgets.QApplication(sys.argv)
    mainwindow = Calibration()
    mainwindow.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()