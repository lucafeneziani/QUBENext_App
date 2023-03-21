import numpy as np
from PyQt5 import QtGui
from PyQt5.QtGui import QPixmap, QFont
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QMainWindow, QDesktopWidget, QPushButton, QLabel, QFileDialog, QMessageBox #, QVBoxLayout, QLineEdit, QWidget, QGridLayout
import pyqtgraph as pg
import os

import functions
import configuration_params
from constants import CLINICAL_RANGE_PERC


########################################################################################################################################

class QApp(QMainWindow):
    def __init__(self):
        super(QApp, self).__init__()
        self.GUIstruct()
    
    def GUIstruct(self):
        width = QDesktopWidget().screenGeometry().width()
        height = QDesktopWidget().screenGeometry().height()
        self.setGeometry(0,0,round(width),round(height))
        self.setWindowTitle("QApp")
        self.setStyleSheet('background: whitesmoke')

        #########################################################################
        # ENABLE/DISABLE VARIABLES

        self.mlic_enable     = False
        self.bortfeld_enable = False
        self.calibZ_enable   = False
        self.xpos_enable     = False
        self.ypos_enable     = False
        self.calibX_enable   = False
        self.calibY_enable   = False
        self.isflipped       = False

        #########################################################################
        # GUI STYLE
        
        button_style = 'background-color: None'
        enable_style = 'background-color: None'
        textstyle = QFont()
        textstyle.setBold(True)
        textstyle2 = QFont(None, 15)
        textstyle2.setBold(True)
        self.pen_data = pg.mkPen(color=(255, 0, 0), width=3)
        self.pen_fit = pg.mkPen(color=(0, 0, 255), width=3, style=Qt.DashLine)

        #########################################################################
        # CONFIGURATION PARAMS

        self.detector_params  = configuration_params.configs

        #########################################################################
        # TOP

        # Bar
        self.barlabel = QLabel(self)
        self.pixmap = QPixmap('/Users/lucafeneziani/Desktop/QUBENext_App/images/bar.png')
        self.barlabel.setPixmap(self.pixmap)
        self.barlabel.resize(round(width), round(0.08*height))
        self.barlabel.move(round(0.0*width), round(0.0*height))

        # De.Tec.Tor Logo
        self.logolabel = QLabel(self)
        self.pixmap = QPixmap('/Users/lucafeneziani/Desktop/QUBENext_App/images/DeTecTor.png')
        self.pixmap = self.pixmap.scaled(round(0.25*width), round(0.1*height),Qt.KeepAspectRatio,Qt.SmoothTransformation)
        self.logolabel.setPixmap(self.pixmap)
        self.logolabel.resize(round(0.27*width), round(0.08*height))
        self.logolabel.move(round(0.73*width), round(0.0*height))
        self.logolabel.setStyleSheet('background-color: None')

        # Device logo
        self.logolabel = QLabel(self)
        self.pixmap = QPixmap('/Users/lucafeneziani/Desktop/QUBENext_App/images/detectors/logo_QUBENext.png')
        #self.pixmap = QPixmap('./images/detectors/logo_'+self.detector_params['detector_type']+'.png')
        self.pixmap = self.pixmap.scaled(round(0.15*width), round(0.1*height),Qt.KeepAspectRatio,Qt.SmoothTransformation)
        self.logolabel.setPixmap(self.pixmap)
        self.logolabel.resize(round(0.27*width), round(0.08*height))
        self.logolabel.move(round(0.05*width), round(0.0*height))
        self.logolabel.setStyleSheet('background-color: None')

        #########################################################################
        # LEFT SIDE
        
        # Z Plot
        self.ZPlot = pg.PlotWidget(self)
        self.ZPlot.setLabel('left','counts')
        self.ZPlot.setLabel('bottom','channels')
        self.ZPlot.setTitle('Profile Z total counts')
        self.ZPlot.resize(round(0.65*width), round(0.25*height))
        self.ZPlot.move(round(0.02*width), round(0.1*height))
        self.ZPlot.setBackground('w')

        # Shaw Z Raw Data button
        self.shawZraw = QPushButton(self)
        self.shawZraw.setText("Raw Data")
        self.shawZraw.resize(round(0.07*width), round(0.05*height))
        self.shawZraw.move(round(0.67*width), round(0.1*height))
        self.shawZraw.clicked.connect(self.Shaw_Z_Data)
        self.shawZraw.setStyleSheet(button_style)
        self.shawZraw.setEnabled(False)

        # Shaw Z Analyzed Data button
        self.shawZfit = QPushButton(self)
        self.shawZfit.setText("Fit Data")
        self.shawZfit.resize(round(0.07*width), round(0.05*height))
        self.shawZfit.move(round(0.67*width), round(0.15*height))
        self.shawZfit.clicked.connect(self.Shaw_Z_Fit)
        self.shawZfit.setStyleSheet(button_style)
        self.shawZfit.setEnabled(False)

        # Reset Z Plot button
        self.resetZplot = QPushButton(self)
        self.resetZplot.setText("Reset axes")
        self.resetZplot.resize(round(0.07*width), round(0.05*height))
        self.resetZplot.move(round(0.67*width), round(0.3*height))
        self.resetZplot.clicked.connect(lambda: self.ZPlot.getPlotItem().enableAutoRange())
        self.resetZplot.setStyleSheet(button_style)
        self.resetZplot.setEnabled(False)

        # X Plot
        self.XPlot = pg.PlotWidget(self)
        self.XPlot.setLabel('left','counts')
        self.XPlot.setLabel('bottom','channels')
        self.XPlot.setTitle('Profile X total counts')
        self.XPlot.resize(round(0.28*width), round(0.25*height))
        self.XPlot.move(round(0.02*width), round(0.4*height))
        self.XPlot.setBackground('w')

        # Shaw X Raw Data button
        self.shawXraw = QPushButton(self)
        self.shawXraw.setText("Raw Data")
        self.shawXraw.resize(round(0.07*width), round(0.05*height))
        self.shawXraw.move(round(0.3*width), round(0.4*height))
        self.shawXraw.clicked.connect(self.Shaw_X_Data)
        self.shawXraw.setStyleSheet(button_style)
        self.shawXraw.setEnabled(False)

        # Shaw X Analyzed Data button
        self.shawXfit = QPushButton(self)
        self.shawXfit.setText("Fit Data")
        self.shawXfit.resize(round(0.07*width), round(0.05*height))
        self.shawXfit.move(round(0.3*width), round(0.45*height))
        self.shawXfit.clicked.connect(self.Shaw_X_Fit)
        self.shawXfit.setStyleSheet(button_style)
        self.shawXfit.setEnabled(False)

        # Reset X Plot button
        self.resetXplot = QPushButton(self)
        self.resetXplot.setText("Reset axes")
        self.resetXplot.resize(round(0.07*width), round(0.05*height))
        self.resetXplot.move(round(0.3*width), round(0.6*height))
        self.resetXplot.clicked.connect(lambda: self.XPlot.getPlotItem().enableAutoRange())
        self.resetXplot.setStyleSheet(button_style)
        self.resetXplot.setEnabled(False)

        # Y Plot
        self.YPlot = pg.PlotWidget(self)
        self.YPlot.setLabel('left','counts')
        self.YPlot.setLabel('bottom','channels')
        self.YPlot.setTitle('Profile Y total counts')
        self.YPlot.resize(round(0.28*width), round(0.25*height))
        self.YPlot.move(round(0.39*width), round(0.4*height))
        self.YPlot.setBackground('w')
        
        # Shaw Y Raw Data button
        self.shawYraw = QPushButton(self)
        self.shawYraw.setText("Raw Data")
        self.shawYraw.resize(round(0.07*width), round(0.05*height))
        self.shawYraw.move(round(0.67*width), round(0.4*height))
        self.shawYraw.clicked.connect(self.Shaw_Y_Data)
        self.shawYraw.setStyleSheet(button_style)
        self.shawYraw.setEnabled(False)

        # Shaw Y Analyzed Data button
        self.shawYfit = QPushButton(self)
        self.shawYfit.setText("Fit Data")
        self.shawYfit.resize(round(0.07*width), round(0.05*height))
        self.shawYfit.move(round(0.67*width), round(0.45*height))
        self.shawYfit.clicked.connect(self.Shaw_Y_Fit)
        self.shawYfit.setStyleSheet(button_style)
        self.shawYfit.setEnabled(False)

        # Reset Y Plot button
        self.resetYplot = QPushButton(self)
        self.resetYplot.setText("Reset axes")
        self.resetYplot.resize(round(0.07*width), round(0.05*height))
        self.resetYplot.move(round(0.67*width), round(0.6*height))
        self.resetYplot.clicked.connect(lambda: self.YPlot.getPlotItem().enableAutoRange())
        self.resetYplot.setStyleSheet(button_style)
        self.resetYplot.setEnabled(False)
        
        #########################################################################
        # RESULTS LABEL

        self.labelResults = QLabel(self)
        self.labelResults.move(round(0.02*width), round(0.7*height))
        self.labelResults.resize(round(0.75*width), round(0.25*height))
        self.labelResults.setStyleSheet('background-color: None; border: 2px solid blue; border-radius: 10px')
        #
        self.labelResults = QLabel(self)
        self.labelResults.move(round(0.02*width), round(0.71*height))
        self.labelResults.resize(round(0.75*width), round(0.05*height))
        self.labelResults.setText('\t\tProfile Z\t\t\t\t\t Profile X\t\t\t\t\t  Profile Y')
        self.labelResults.setFont(textstyle2)
        self.labelResults.setStyleSheet('background-color: None')
        
        # Profile Z
        self.labelResultsZt = QLabel(self)
        self.labelResultsZt.move(round(0.04*width), round(0.7*height))
        self.labelResultsZt.resize(round(0.15*width), round(0.25*height))
        self.labelResultsZt.setText('\n\n\nPeak position:\n\nPeak-plateau ratio:\n\nClinical range (R{:d}):\n\nPeak width (@{:d}%):'.format(int(CLINICAL_RANGE_PERC*100),int(CLINICAL_RANGE_PERC*100)))
        self.labelResultsZt.setFont(textstyle)
        self.labelResultsZt.setStyleSheet('background-color: None')
        #
        self.labelResultsZ = QLabel(self)
        self.labelResultsZ.move(round(0.16*width), round(0.7*height))
        self.labelResultsZ.resize(round(0.10*width), round(0.25*height))
        self.labelResultsZ.setText('\n\n\n--\n\n--\n\n--\n\n--')
        self.labelResultsZ.setFont(textstyle)
        self.labelResultsZ.setStyleSheet('background-color: None')
        
        # Profile X
        self.labelResultsXt = QLabel(self)
        self.labelResultsXt.move(round(0.29*width), round(0.7*height))
        self.labelResultsXt.resize(round(0.15*width), round(0.25*height))
        self.labelResultsXt.setText('\n\n\nMean:\n\nSigma:\n\nFWHM:\n\nWidth (@90%):')
        self.labelResultsXt.setFont(textstyle)
        self.labelResultsXt.setStyleSheet("background-color: None")
        #
        self.labelResultsX = QLabel(self)
        self.labelResultsX.move(round(0.41*width), round(0.7*height))
        self.labelResultsX.resize(round(0.10*width), round(0.25*height))
        self.labelResultsX.setText('\n\n\n--\n\n--\n\n--\n\n--')
        self.labelResultsX.setFont(textstyle)
        self.labelResultsX.setStyleSheet("background-color: None")
        
        # Profile Y
        self.labelResultsYt = QLabel(self)
        self.labelResultsYt.move(round(0.54*width), round(0.7*height))
        self.labelResultsYt.resize(round(0.15*width), round(0.25*height))
        self.labelResultsYt.setText('\n\n\nMean:\n\nSigma:\n\nFWHM:\n\nWidth (@90%):')
        self.labelResultsYt.setFont(textstyle)
        self.labelResultsYt.setStyleSheet("background-color: None")
        #
        self.labelResultsY = QLabel(self)
        self.labelResultsY.move(round(0.66*width), round(0.7*height))
        self.labelResultsY.resize(round(0.10*width), round(0.25*height))
        self.labelResultsY.setText('\n\n\n--\n\n--\n\n--\n\n--')
        self.labelResultsY.setFont(textstyle)
        self.labelResultsY.setStyleSheet("background-color: None")

        #########################################################################
        # RIGHT SIDE

        # Z File load button
        self.loadZdata = QPushButton(self)
        self.loadZdata.setText("Load Z Data")
        self.loadZdata.resize(round(0.2*width), round(0.05*height))
        self.loadZdata.move(round(0.79*width), round(0.1*height))
        self.loadZdata.clicked.connect(self.Load_Z_Data)
        self.loadZdata.setStyleSheet(button_style)

        # Z File label
        self.labelZfile = QLabel(self) 
        self.labelZfile.setText("File IN Z:")
        self.labelZfile.resize(round(0.045*width), round(0.03*height))
        self.labelZfile.move(round(0.8*width), round(0.148*height))
        self.labelZfile.setStyleSheet('border: 1px solid gray; border: None')
        #
        self.labelZfile = QLabel(self) 
        self.labelZfile.setText('')
        self.labelZfile.resize(round(0.135*width), round(0.03*height))
        self.labelZfile.move(round(0.845*width), round(0.148*height))
        self.labelZfile.setStyleSheet('background-color: lightgray; border: None')

        # Z File reverse button
        self.reverseZdata = QPushButton(self)
        self.reverseZdata.setText("Reverse Z Data")
        self.reverseZdata.resize(round(0.2*width), round(0.05*height))
        self.reverseZdata.move(round(0.79*width), round(0.18*height))
        self.reverseZdata.clicked.connect(self.Reverse_Z_Data)
        self.reverseZdata.setStyleSheet(button_style)
        self.reverseZdata.setEnabled(False)

        # Enable Z switch button
        self.enableZ = QPushButton(self)
        self.enableZ.setText("Enable Z")
        self.enableZ.resize(round(0.1*width), round(0.05*height))
        self.enableZ.move(round(0.79*width), round(0.23*height))
        self.enableZ.clicked.connect(self.Enable_Z)
        self.enableZ.setStyleSheet(enable_style)
        self.enableZ.setEnabled(False)

        # Enable Bortfeld fit switch button
        self.enableBortfeld = QPushButton(self)
        self.enableBortfeld.setText("Bortfeld fit")
        self.enableBortfeld.resize(round(0.1*width), round(0.05*height))
        self.enableBortfeld.move(round(0.89*width), round(0.23*height))
        self.enableBortfeld.clicked.connect(self.Enable_Bortfeld)
        self.enableBortfeld.setStyleSheet(enable_style)
        self.enableBortfeld.setEnabled(False)

        # Z Calibration load button
        self.loadZcalib = QPushButton(self)
        self.loadZcalib.setText("Load Z Calib")
        self.loadZcalib.resize(round(0.2*width), round(0.05*height))
        self.loadZcalib.move(round(0.79*width), round(0.28*height))
        self.loadZcalib.clicked.connect(self.Load_Z_Calibration)
        self.loadZcalib.setStyleSheet(button_style)
        self.loadZcalib.setEnabled(False)

        # Z Calib label 
        self.labelZcalib = QLabel(self) 
        self.labelZcalib.setText("Calib Z:")  
        self.labelZcalib.resize(round(0.045*width), round(0.03*height))
        self.labelZcalib.move(round(0.8*width), round(0.328*height))
        self.labelZcalib.setStyleSheet('border: 1px solid gray; border: None')
        #
        self.labelZcalib = QLabel(self) 
        self.labelZcalib.setText('')
        self.labelZcalib.resize(round(0.135*width), round(0.03*height))
        self.labelZcalib.move(round(0.845*width), round(0.328*height))
        self.labelZcalib.setStyleSheet('background-color: lightgray; None')

        # Enable Z Calib switch button
        self.enableZcalib = QPushButton(self)
        self.enableZcalib.setText("Apply calib Z")
        self.enableZcalib.resize(round(0.1*width), round(0.05*height))
        self.enableZcalib.move(round(0.79*width), round(0.36*height))
        self.enableZcalib.clicked.connect(self.Apply_Z_Calib)
        self.enableZcalib.setStyleSheet(enable_style)
        self.enableZcalib.setEnabled(False)

        # X File load button
        self.loadXdata = QPushButton(self)
        self.loadXdata.setText("Load X Data")
        self.loadXdata.resize(round(0.2*width), round(0.05*height))
        self.loadXdata.move(round(0.79*width), round(0.4*height))
        self.loadXdata.clicked.connect(self.Load_X_Data)
        self.loadXdata.setStyleSheet(button_style)

        # X File label 
        self.labelXfile = QLabel(self) 
        self.labelXfile.setText("File IN X:")  
        self.labelXfile.resize(round(0.045*width), round(0.03*height))
        self.labelXfile.move(round(0.8*width), round(0.448*height))
        self.labelXfile.setStyleSheet('border: 1px solid gray; border: None')
        #
        self.labelXfile = QLabel(self) 
        self.labelXfile.setText('')
        self.labelXfile.resize(round(0.135*width), round(0.03*height))
        self.labelXfile.move(round(0.845*width), round(0.448*height))
        self.labelXfile.setStyleSheet('background-color: lightgray; border: None')

        # Enable X switch button
        self.enableX = QPushButton(self)
        self.enableX.setText("Enable X")
        self.enableX.resize(round(0.1*width), round(0.05*height))
        self.enableX.move(round(0.79*width), round(0.48*height))
        self.enableX.clicked.connect(self.Enable_X)
        self.enableX.setStyleSheet(enable_style)
        self.enableX.setEnabled(False)

        # Y File load button
        self.loadYdata = QPushButton(self)
        self.loadYdata.setText("Load Y Data")
        self.loadYdata.resize(round(0.2*width), round(0.05*height))
        self.loadYdata.move(round(0.79*width), round(0.53*height))
        self.loadYdata.clicked.connect(self.Load_Y_Data)
        self.loadYdata.setStyleSheet(button_style)

        # Y File label 
        self.labelYfile = QLabel(self) 
        self.labelYfile.setText("File IN Y:")  
        self.labelYfile.resize(round(0.045*width), round(0.03*height))
        self.labelYfile.move(round(0.8*width), round(0.578*height))
        self.labelYfile.setStyleSheet('border: 1px solid gray; border: None')
        #
        self.labelYfile = QLabel(self) 
        self.labelYfile.setText('')
        self.labelYfile.resize(round(0.135*width), round(0.03*height))
        self.labelYfile.move(round(0.845*width), round(0.578*height))
        self.labelYfile.setStyleSheet('background-color: lightgray; border: None')

        # Enable Y switch button
        self.enableY = QPushButton(self)
        self.enableY.setText("Enable Y")
        self.enableY.resize(round(0.1*width), round(0.05*height))
        self.enableY.move(round(0.79*width), round(0.61*height))
        self.enableY.clicked.connect(self.Enable_Y)
        self.enableY.setStyleSheet(enable_style)
        self.enableY.setEnabled(False)

        # X/Y Calibration load button
        self.loadXYcalib = QPushButton(self)
        self.loadXYcalib.setText("Load X/Y Calib")
        self.loadXYcalib.resize(round(0.2*width), round(0.05*height))
        self.loadXYcalib.move(round(0.79*width), round(0.66*height))
        self.loadXYcalib.clicked.connect(self.Load_XY_Calibration)
        self.loadXYcalib.setStyleSheet(button_style)
        self.loadXYcalib.setEnabled(False)

        # X/Y Calib label 
        self.labelXYcalib = QLabel(self) 
        self.labelXYcalib.setText("Calib X/Y:")  
        self.labelXYcalib.resize(round(0.045*width), round(0.03*height))
        self.labelXYcalib.move(round(0.8*width), round(0.708*height))
        self.labelXYcalib.setStyleSheet('border: 1px solid gray; border: None')
        #
        self.labelXYcalib = QLabel(self) 
        self.labelXYcalib.setText('')
        self.labelXYcalib.resize(round(0.135*width), round(0.03*height))
        self.labelXYcalib.move(round(0.845*width), round(0.708*height))
        self.labelXYcalib.setStyleSheet('background-color: lightgray; border: None')

        # Enable X Calib switch button
        self.enableXcalib = QPushButton(self)
        self.enableXcalib.setText("Apply calib X")
        self.enableXcalib.resize(round(0.1*width), round(0.05*height))
        self.enableXcalib.move(round(0.79*width), round(0.74*height))
        self.enableXcalib.clicked.connect(self.Apply_X_Calib)
        self.enableXcalib.setStyleSheet(enable_style)
        self.enableXcalib.setEnabled(False)

        # Enable Y Calib switch button
        self.enableYcalib = QPushButton(self)
        self.enableYcalib.setText("Apply calib Y")
        self.enableYcalib.resize(round(0.1*width), round(0.05*height))
        self.enableYcalib.move(round(0.89*width), round(0.74*height))
        self.enableYcalib.clicked.connect(self.Apply_Y_Calib)
        self.enableYcalib.setStyleSheet(enable_style)
        self.enableYcalib.setEnabled(False)

        # Reset All button
        self.resetall = QPushButton(self)
        self.resetall.setText("Reset all")
        self.resetall.resize(round(0.2*width), round(0.05*height))
        self.resetall.move(round(0.79*width), round(0.8*height))
        self.resetall.clicked.connect(self.Reset_all)
        self.resetall.setStyleSheet('background-color: lightgray; color: red')

        # Analyze button
        self.analyze = QPushButton(self)
        self.analyze.setText("ANALYZE")
        self.analyze.setFont(textstyle)
        self.analyze.resize(round(0.2*width), round(0.1*height))
        self.analyze.move(round(0.79*width), round(0.85*height))
        self.analyze.clicked.connect(self.Analyze)
        self.analyze.setStyleSheet('background-color: lightgreen')
        self.analyze.setEnabled(False)

        #########################################################################

        

    ########################################################################################################################################
    def Load_Z_Data(self):
            
        try:
            '''
            data = np.loadtxt('./data/mlic/20161019_011009_QUBE.dat', dtype=str, delimiter = '\t')
            '''
            file = QFileDialog.getOpenFileName(self, os.getcwd())[0]
            data = np.loadtxt(file, dtype=str, delimiter = '\t')
            self.labelZfile.setText(file.split('/')[-1])
            
            self.Z_Time = data[-1][0].astype(float)
            self.Z_Data = data[-1][1::].astype(float)

            self.Z_data_x = range(len(self.Z_Data))
            self.Z_data_y = self.Z_Data
            self.ZPlot.clear()
            self.Zraw = self.ZPlot.plot(self.Z_data_x, self.Z_data_y, pen = self.pen_data)
            self.ZPlot.setLabel('bottom','channels')
            self.ZPlot.getPlotItem().enableAutoRange()

            self.shawZraw.setEnabled(True)
            self.reverseZdata.setEnabled(True)
            self.enableZ.setEnabled(True)
            self.enableBortfeld.setEnabled(True)
            self.loadZcalib.setEnabled(True)
            self.analyze.setEnabled(True)
            self.resetZplot.setEnabled(True)

        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Warning")
            msg.setInformativeText('Unrecognized data format - Check the upload')
            msg.exec_()


    def Load_X_Data(self):

        try:
            '''
            data = np.loadtxt('./data/strip/20151127_005129_miniQ-X.dat', dtype=str, delimiter = '\t')
            '''
            file = QFileDialog.getOpenFileName(self, os.getcwd())[0]
            data = np.loadtxt(file, dtype=str, delimiter = '\t')
            self.labelXfile.setText(file.split('/')[-1])
            
            self.X_Time = data[-1][0].astype(float)
            if self.detector_params['strip_integral'] == 'yes':
                self.X_Data = data[-1][1:-1].astype(float)
                self.X_integral = data[-1][-1].astype(float)
            else:
                self.X_Data = data[-1][1::].astype(float)

            self.X_data_x = range(len(self.X_Data))
            self.X_data_y = self.X_Data
            self.XPlot.clear()
            self.Xraw = self.XPlot.plot(self.X_data_x, self.X_data_y, pen = self.pen_data)
            self.XPlot.setLabel('bottom','channels')
            self.XPlot.getPlotItem().enableAutoRange()
            
            self.shawXraw.setEnabled(True)
            self.enableX.setEnabled(True)
            self.loadXYcalib.setEnabled(True)
            self.analyze.setEnabled(True)
            self.resetXplot.setEnabled(True)
            self.Xdataplot = True
        
        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Warning")
            msg.setInformativeText('Unrecognized data format - Check the upload')
            msg.exec_()

  

    def Load_Y_Data(self):

        try:
            '''
            data = np.loadtxt('./data/strip/20151127_005129_miniQ-Y.dat', dtype=str, delimiter = '\t')
            '''
            file = QFileDialog.getOpenFileName(self, os.getcwd())[0]
            data = np.loadtxt(file, dtype=str, delimiter = '\t')
            self.labelYfile.setText(file.split('/')[-1])
        
            self.Y_Time = data[-1][0].astype(float)
            if self.detector_params['strip_integral'] == 'yes':
                self.Y_Data = data[-1][1:-1].astype(float)
                self.Y_integral = data[-1][-1].astype(float)
            else:
                self.Y_Data = data[-1][1::].astype(float)

            self.Y_data_x = range(len(self.Y_Data))
            self.Y_data_y = self.Y_Data
            self.YPlot.clear()
            self.Yraw = self.YPlot.plot(self.Y_data_x, self.Y_data_y, pen = self.pen_data)
            self.YPlot.setLabel('bottom','channels')
            self.YPlot.getPlotItem().enableAutoRange()
            
            self.shawYraw.setEnabled(True)
            self.enableY.setEnabled(True)
            self.loadXYcalib.setEnabled(True)
            self.analyze.setEnabled(True)
            self.resetYplot.setEnabled(True)
            self.Ydataplot = True
        
        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Warning")
            msg.setInformativeText('Unrecognized data format - Check the upload')
            msg.exec_()

    def Reverse_Z_Data(self):
        self.Z_data_y = np.flip(self.Z_data_y)
        if self.isflipped:
            self.isflipped = False
        else:
            self.isflipped = True

        self.ZPlot.removeItem(self.Zraw)
        self.Zraw = self.ZPlot.plot(self.Z_data_x, self.Z_data_y, pen = self.pen_data)

        return
    

    def Enable_Z(self):
        if self.mlic_enable:
            self.mlic_enable = False
            self.enableZ.setStyleSheet('background-color: None; color: None')
        else:
            self.mlic_enable = True
            self.enableZ.setStyleSheet('background-color: None; color: green')
        return


    def Enable_Bortfeld(self):
        if self.bortfeld_enable:
            self.bortfeld_enable = False
            self.enableBortfeld.setStyleSheet('background-color: None; color: None')
        else:
            self.bortfeld_enable = True
            self.enableBortfeld.setStyleSheet('background-color: None; color: green')
        return
    

    def Load_Z_Calibration(self):
        
        try:
            file = QFileDialog.getOpenFileName(self, os.getcwd())[0]
            data = np.loadtxt(file, dtype=str, delimiter = '\t')
            self.labelZcalib.setText(file.split('/')[-1])

            self.calibZ_vector_orig = np.transpose(data)[1][1::].astype(float)
            self.calibZ_vector_flip = np.flip(self.calibZ_vector_orig)

            self.enableZcalib.setEnabled(True)

        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Warning")
            msg.setInformativeText('Unrecognized data format - Check the upload')
            msg.exec_()

        return
    

    def Load_XY_Calibration(self):

        try:
            file = QFileDialog.getOpenFileName(self, os.getcwd())[0]
            data = np.loadtxt(file, dtype=str, delimiter = '\t')
            self.labelXYcalib.setText(file.split('/')[-1])

            self.calibX_vector = np.transpose(data)[1][1::].astype(float)
            self.calibY_vector = np.transpose(data)[2][1::].astype(float)

            self.enableXcalib.setEnabled(True)
            self.enableYcalib.setEnabled(True)
        
        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Warning")
            msg.setInformativeText('Unrecognized data format - Check the upload')
            msg.exec_()

        return
    

    def Apply_Z_Calib(self):
        try:
            if self.isflipped:
                self.calibZ_vector = self.calibZ_vector_flip
            else:
                self.calibZ_vector = self.calibZ_vector_orig

            if self.calibZ_enable:
                self.Z_data_y = self.Z_data_y / self.calibZ_vector
                self.ZPlot.removeItem(self.Zraw)
                self.Zraw = self.ZPlot.plot(self.Z_data_x, self.Z_data_y, pen = self.pen_data)
                self.calibZ_enable = False
                self.enableZcalib.setStyleSheet('background-color: None; color: None')
            else:
                self.Z_data_y = self.Z_data_y * self.calibZ_vector
                self.ZPlot.removeItem(self.Zraw)
                self.Zraw = self.ZPlot.plot(self.Z_data_x, self.Z_data_y, pen = self.pen_data)
                self.calibZ_enable = True
                self.enableZcalib.setStyleSheet('background-color: None; color: green')

        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Warning")
            msg.setInformativeText('Error during calibration applying')
            msg.exec_()

        return
    

    def Apply_X_Calib(self):

        try:
            if self.calibX_enable:
                self.X_data_y = self.X_data_y / self.calibX_vector
                self.XPlot.removeItem(self.Xraw)
                self.Xraw = self.XPlot.plot(self.X_data_x, self.X_data_y, pen = self.pen_data)
                self.calibX_enable = False
                self.enableXcalib.setStyleSheet('background-color: None; color: None')
            else:
                self.X_data_y = self.X_data_y * self.calibX_vector
                self.XPlot.removeItem(self.Xraw)
                self.Xraw = self.XPlot.plot(self.X_data_x, self.X_data_y, pen = self.pen_data)
                self.calibX_enable = True
                self.enableXcalib.setStyleSheet('background-color: None; color: green')

        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Warning")
            msg.setInformativeText('Error during calibration applying')
            msg.exec_()
        return
    

    def Apply_Y_Calib(self):

        try:
            if self.calibY_enable:
                self.Y_data_y = self.Y_data_y / self.calibY_vector
                self.YPlot.removeItem(self.Yraw)
                self.Yraw = self.YPlot.plot(self.Y_data_x, self.Y_data_y, pen = self.pen_data)
                self.calibY_enable = False
                self.enableYcalib.setStyleSheet('background-color: None; color: None')
            else:
                self.Y_data_y = self.Y_data_y * self.calibY_vector
                self.YPlot.removeItem(self.Yraw)
                self.Yraw = self.YPlot.plot(self.Y_data_x, self.Y_data_y, pen = self.pen_data)
                self.calibY_enable = True
                self.enableYcalib.setStyleSheet('background-color: None; color: green')
        
        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Warning")
            msg.setInformativeText('Error during calibration applying')
            msg.exec_()

        return
    

    def Enable_X(self):
        if self.xpos_enable:
            self.xpos_enable = False
            self.enableX.setStyleSheet('background-color: None; color: None')
        else:
            self.xpos_enable = True
            self.enableX.setStyleSheet('background-color: None; color: green')
        return
    

    def Enable_Y(self):
        if self.ypos_enable:
            self.ypos_enable = False
            self.enableY.setStyleSheet('background-color: None; color: None')
        else:
            self.ypos_enable = True
            self.enableY.setStyleSheet('background-color: None; color: green')
        return
    

    def Shaw_Z_Data(self):
        if self.Zraw.isVisible():
            self.Zraw.hide()
        else:
            self.Zraw.show()
        return
    

    def Shaw_Z_Fit(self):
        if self.Zfit.isVisible():
            self.Zfit.hide()
        else:
            self.Zfit.show()
        return

    
    def Shaw_X_Data(self):
        if self.Xraw.isVisible():
            self.Xraw.hide()
        else:
            self.Xraw.show()
        return
    

    def Shaw_X_Fit(self):
        if self.Xfit.isVisible():
            self.Xfit.hide()
        else:
            self.Xfit.show()
        return

    def Shaw_Y_Data(self):
        if self.Yraw.isVisible():
            self.Yraw.hide()
        else:
            self.Yraw.show()
        return
    

    def Shaw_Y_Fit(self):
        if self.Yfit.isVisible():
            self.Yfit.hide()
        else:
            self.Yfit.show()
        return
    
    
    def Reset_all(self):

        # Plot
        self.ZPlot.clear()
        self.XPlot.clear()
        self.YPlot.clear()
        self.ZPlot.setLabel('bottom','channels')
        self.XPlot.setLabel('bottom','channels')
        self.YPlot.setLabel('bottom','channels')

        # Data
        self.Z_Data = []
        self.X_Data = []
        self.Y_Data = []

        # Variables
        self.mlic_enable     = False
        self.bortfeld_enable = False
        self.calibZ_enable   = False
        self.xpos_enable     = False
        self.ypos_enable     = False
        self.calibX_enable   = False
        self.calibY_enable   = False
        self.isflipped       = False

        # Label
        self.labelZfile.setText('')
        self.labelXfile.setText('')
        self.labelYfile.setText('')
        self.labelZcalib.setText('')
        self.labelXYcalib.setText('')
        self.labelResultsZ.setText('\n\n\n--\n\n--\n\n--\n\n--')
        self.labelResultsX.setText('\n\n\n--\n\n--\n\n--\n\n--')
        self.labelResultsY.setText('\n\n\n--\n\n--\n\n--\n\n--')
        self.labelResults.setText('')
        self.enableZ.setStyleSheet('background-color: None; color: None')
        self.enableX.setStyleSheet('background-color: None; color: None')
        self.enableY.setStyleSheet('background-color: None; color: None')
        self.enableBortfeld.setStyleSheet('background-color: None; color: None')
        self.enableZcalib.setStyleSheet('background-color: None; color: None')
        self.enableXcalib.setStyleSheet('background-color: None; color: None')
        self.enableYcalib.setStyleSheet('background-color: None; color: None')

        # button
        self.shawZraw.setEnabled(False)
        self.shawZfit.setEnabled(False)
        self.shawXraw.setEnabled(False)
        self.shawXfit.setEnabled(False)
        self.shawYraw.setEnabled(False)
        self.shawYfit.setEnabled(False)
        self.reverseZdata.setEnabled(False)
        self.enableZ.setEnabled(False)
        self.enableBortfeld.setEnabled(False)
        self.loadZcalib.setEnabled(False)
        self.enableZcalib.setEnabled(False)
        self.enableX.setEnabled(False)
        self.enableY.setEnabled(False)
        self.loadXYcalib.setEnabled(False)
        self.enableXcalib.setEnabled(False)
        self.enableYcalib.setEnabled(False)
        self.analyze.setEnabled(False)
        self.resetZplot.setEnabled(False)
        self.resetXplot.setEnabled(False)
        self.resetYplot.setEnabled(False)
        
        return
    

    def Check_direction(self):
        left_edge  = np.mean(self.Z_data_y[0:10])
        right_edge = np.mean(self.Z_data_y[len(self.Z_data_y)-10::])
        
        if left_edge < right_edge:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
            msg.setText("Warning")
            msg.setInformativeText('The Z data seems to be reversed, do you want to continue?')
            ret = msg.exec_()

            if ret == QMessageBox.Ok:
                return 'ok'
            else:
                return 'no'

    
    def Analyze(self):

        self.labelResultsZ.setText('\n\n\n--\n\n--\n\n--\n\n--')
        self.labelResultsX.setText('\n\n\n--\n\n--\n\n--\n\n--')
        self.labelResultsY.setText('\n\n\n--\n\n--\n\n--\n\n--')
        
        try:
            self.Zfit.hide()
        except:
            pass
        try:
            self.Xfit.hide()
        except:
            pass
        try:
            self.Yfit.hide()
        except:
            pass
            

        # MLIC Analysis (Z Profile)
        if self.mlic_enable:

            if self.Check_direction() == 'no':
                return

            self.Zres = functions.mlic_analysis(self.Z_data_y, self.bortfeld_enable)
            self.Z_data_x = self.Zres['coordinates_raw']
            self.Z_data_y = self.Zres['raw_data']
            self.Z_fit_x = self.Zres['coordinates_fit']
            self.Z_fit_y = self.Zres['fit_data']

            self.ZPlot.clear()
            self.ZPlot.setLabel('bottom','depth [cm w.e.]')
            self.Zraw = self.ZPlot.plot(self.Z_data_x, self.Z_data_y, pen = self.pen_data)
            self.Zfit = self.ZPlot.plot(self.Z_fit_x, self.Z_fit_y, pen = self.pen_fit)
            self.ZPlot.getPlotItem().enableAutoRange()
            self.shawZfit.setEnabled(True)
            self.labelResultsZ.setText('\n\n\n{:.2f}\t{}\n\n{:.2f}\t{}\n\n{:.2f}\t{}\n\n{:.2f}\t{}'.format(self.Zres['peak_pos']['value'],self.Zres['peak_pos']['unit'],
                                                                                           self.Zres['pp_ratio']['value'],self.Zres['pp_ratio']['unit'],
                                                                                           self.Zres['cl_range']['value'],self.Zres['cl_range']['unit'],
                                                                                           self.Zres['peak_width']['value'],self.Zres['peak_width']['unit']))
        
        # STRIP Analysis (X Profile)
        if self.xpos_enable:

            self.Xres = functions.strip_analysis(self.X_data_y)
            self.X_data_x = self.Xres['coordinates_raw']
            self.X_data_y = self.Xres['raw_data']
            self.X_fit_x = self.Xres['coordinates_fit']
            self.X_fit_y = self.Xres['fit_data']

            self.XPlot.clear()
            self.XPlot.setLabel('bottom','[mm]')
            self.Xraw = self.XPlot.plot(self.X_data_x, self.X_data_y, pen = self.pen_data)
            self.Xfit = self.XPlot.plot(self.X_fit_x, self.X_fit_y, pen = self.pen_fit)
            self.XPlot.getPlotItem().enableAutoRange()
            self.shawXfit.setEnabled(True)
            self.labelResultsX.setText('\n\n\n{:.2f}\t{}\n\n{:.2f}\t{}\n\n{:.2f}\t{}\n\n{:.2f}\t{}'.format(self.Xres['mean']['value'],self.Xres['mean']['unit'],
                                                                                           self.Xres['sigma']['value'],self.Xres['sigma']['unit'],
                                                                                           self.Xres['fwhm']['value'],self.Xres['fwhm']['unit'],
                                                                                           self.Xres['peak_width']['value'],self.Xres['peak_width']['unit']))
        
        # STRIP Analysis (Y Profile)
        if self.ypos_enable:

            self.Yres = functions.strip_analysis(self.Y_data_y)
            self.Y_data_x = self.Yres['coordinates_raw']
            self.Y_data_y = self.Yres['raw_data']
            self.Y_fit_x = self.Yres['coordinates_fit']
            self.Y_fit_y = self.Yres['fit_data']

            self.YPlot.clear()
            self.YPlot.setLabel('bottom','[mm]')
            self.Yraw = self.YPlot.plot(self.Y_data_x, self.Y_data_y, pen = self.pen_data)
            self.Yfit = self.YPlot.plot(self.Y_fit_x, self.Y_fit_y, pen = self.pen_fit)
            self.YPlot.getPlotItem().enableAutoRange()
            self.shawYfit.setEnabled(True)
            self.labelResultsY.setText('\n\n\n{:.2f}\t{}\n\n{:.2f}\t{}\n\n{:.2f}\t{}\n\n{:.2f}\t{}'.format(self.Yres['mean']['value'],self.Yres['mean']['unit'],
                                                                                           self.Yres['sigma']['value'],self.Yres['sigma']['unit'],
                                                                                           self.Yres['fwhm']['value'],self.Yres['fwhm']['unit'],
                                                                                           self.Yres['peak_width']['value'],self.Yres['peak_width']['unit']))
        
        if self.mlic_enable == False and self.xpos_enable == False and self.ypos_enable == False:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Warning")
            msg.setInformativeText('Nothing to analyze:\nload files and press enable button')
            msg.exec_()
        
        return