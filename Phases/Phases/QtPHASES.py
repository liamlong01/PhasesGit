#!/usr/bin/python3
# -*- coding: utf-8 -*-


#!/usr/bin/env python



from __future__ import unicode_literals #p lib



import sys #python library
import os #python libraries

import mplCanvas #calls mplcanvas

from PyQt5 import QtGui, QtCore, QtWidgets #p lib for gui codes

from Mesh.ConditionEditor.MeshEditor import MeshEditor #calling mesheditor class from mesheditor file in mesh/conditioneditor folders
from Mesh.Generation.MeshGenerator import MeshGenerator
import Mesh.Mesh as Mesh #calling mesh file from mesh folder

progname = os.path.basename(sys.argv[0])
progversion = "0.1"



from PhasesController import PhasesController #calling class from file
from PlotController import PlotController


global app_path #file mgt
app_path = os.path.dirname(os.path.realpath('__file__'))

class fileInputDialog(QtWidgets.QDialog):

    """This is a dialog window that allows the user to specify files which will be used to define 
    a mesh, its boundary conditions, its initial conditions, some parameters 
    and an optional output text file for storing output data. 
    Once the files have been specified the newFiles_Specified signal will be emited"""
    
    #signals
    newFiles_Specified = QtCore.pyqtSignal(dict) 
    

    def __init__(self,parent, files):
        """
        Initializes the dialog widgets and layout. Two parameters:
            parent -    a reference to the QWidget parent
            files -     a dictionary that must coontain the exact keys: 'Project', 'Mesh', 'Boundary Conditions', 'Initial Conditions', and 'Output'
                        The items must also be strings but they can be anything. 
                        The items will be used as the inital text for the appropriate file input text box.
            """
        
        QtWidgets.QDialog.__init__(self, parent = parent)

        self.setWhatsThis("""PHASES requires three sets of data to start a simulation.
These sets of data can be loaded here through the use of text files in csv format.
The mesh file descibes the geometry of the mesh.
The boundary conditions file specifies the boundaries.
The initial conditions specifies the conditions at the start of the simulation.
Optionally an output text file can be specified to store the calculated values.
For detailed information regarding file formats pleases refer to ____.docx """)
        
        self.setWindowTitle('File Input Dialog')
        
        self.fileInputs =[]
        self.oldfiles = files
        

        self.layout = QtWidgets.QFormLayout(self)
        
        # a temporary dictionary is created so that if
        # the user hits cancel or quit previous settings will not be affected
        self.tempfiles = files.copy()

        self.fileInputs.append(self.addSelector('Project'))
        self.fileInputs[0].fileAdded.connect(self.addProject)

        self.fileInputs.append(self.addSelector('Mesh'))
        self.fileInputs.append(self.addSelector('Initial Conditions'))
        self.fileInputs.append(self.addSelector('Boundary Conditions'))
        self.fileInputs.append(self.addSelector('Output'))           
            
        #TODO: these buttons should be changed to a button box for cross platform readability
        OkButton = QtWidgets.QPushButton('Ok', parent = self)
        OkButton.clicked.connect(self.accept)
        CancelButton = QtWidgets.QPushButton('Cancel', parent = self)
        CancelButton.clicked.connect(self.reject)
        Buttons = QtWidgets.QHBoxLayout()
        Buttons.addWidget(OkButton)
        Buttons.addWidget(CancelButton)

        self.layout.addRow('', Buttons)
    
    def addSelector(self, fileType):
        """
        Initializes a single file selector.
        Dialog requires one for each file to be specified.
        """
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
        fileIn = file_selector(self, fileType, self.tempfiles[fileType])
        fileIn.fileAdded.connect(self.addTempFile)
        self.layout.addRow("%s File:" % fileType, fileIn)
        return fileIn

    def execute(self):
       
        """
        Executes the dialog. Emits and returns specified files if ok.
        Does nothing if user quits dialog or hits cancel.
        """
        if  self.exec(): #ok button clicked
            for file in self.tempfiles:
                print(self.tempfiles[file])
                if not os.path.exists(self.tempfiles[file]):
                    print('error here')
                    errBox = QtWidgets.QErrorMessage(parent = self)
                    errBox.showMessage('File path %s  does not exist' %file)
                    return None
            self.newFiles_Specified.emit(self.tempfiles)
            return self.tempfiles
        else:
            return None #ok button clicked

    def addProject(self, filetype, file): # definetext files path
        """
        When a project text file is added by the user this will
        automatically populate the other files as specified by the 
        last few lines of the project folder.
        
            filetype - should be 'Project' as that is emitted by the file_input class.
                        otherwise this function does nothing
                        
            file - a string of the directory of the chosen project file
        """
        
        if filetype != 'Project':
            return
  
        # open and read file
        prj = open(file,'r')
        Files =[]
        for line in prj.readlines()[-7:]:   # the last 7 lines should contain all file directories
            Files.append(line.strip('\n'))  # also getting rid of pesky '\n' at end of lines
       

        # finding the files 
        self.addTempFile('Mesh',Files[1])
        self.addTempFile('Initial Conditions',Files[2])
        self.addTempFile('Boundary Conditions',Files[3])
        self.addTempFile('Output',Files[4])

        # adding directories to file_input textboxes
        for i in range(4):
            self.fileInputs[i+1].textBox.setText(Files[i+1])
            
        #close file    
        prj.close()


    def addTempFile(self, file_name, file):
        """
        Method to add files to a temporary dictionary for them.
        This method acts as a Qt slot so that the files can be changed
        when appropriate signals are emitted.
        
            file_name - string describing the type of file
            file - a string of the directory of the file
        """
        self.tempfiles[file_name] = file


class file_selector(QtWidgets.QWidget):

    """ 
    A widget for selecting files.
    It contains a label, a text input 
    and a button to open a file selection dialog.
    
    fileAdded is a siganl that will be emitted when a file dialog is used to specify a file.
    The first string is the type of file, and the second string is the file directory.
    
    
    """
    fileAdded = QtCore.pyqtSignal(str, str)

    def __init__(self, parent,file_name,initText):
        """
        parent - reference to parent QWidget
        file_name - name of the file,to be used as the label
        initText - initial text to be placed in the text box.file
        """
        QtWidgets.QWidget.__init__(self, parent=parent)

        layout = QtWidgets.QHBoxLayout(self)
        self.name = file_name
        self.textBox = QtWidgets.QLineEdit(self)
        self.textBox.setText(initText)
        layout.addWidget(self.textBox)

        fileDialog = QtWidgets.QFileDialog(parent)
        fileDialog.fileSelected.connect(self.textBox.setText)
        fileDialog.fileSelected.connect(self.sendSig)

        button = QtWidgets.QPushButton('...', parent = self)
        button.clicked.connect(lambda: fileDialog.exec())
        layout.addWidget(button)

    def sendSig(self, file):
        self.fileAdded.emit(self.name, file)


class ApplicationWindow(QtWidgets.QMainWindow):
    """This is the main window of the entire application. 
    It initializes all the components, controls the layout, 
    and maps the Qt siganls and slots that are passed between different widgets"""

    def __init__(self):
 
        # initialization of th QMainWindow class and attributes
        QtWidgets.QMainWindow.__init__(self)
    
       
       
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        
     
        self.setWindowTitle("application main window")
        
        
        # initializing the main widget up as a tab widget
        # this will allow different widgets to be active at the same time
        # and they will be structured in a typical format with tabs at the top
        self.main_widget = QtWidgets.QTabWidget(self)
        
        
        
        # these properties will be initalized later
        self.mesh = None
        

        # Setup of the matplotlib plot area
        self.plot_widget = QtWidgets.QWidget(self)
        l = QtWidgets.QVBoxLayout(self.plot_widget)
       
        self.mainPlot = mplCanvas.DynamicCanvas(self.plot_widget, dpi=100)
       
        
        l.addWidget(self.mainPlot)

        # adding the plot area to the main widget
        self.main_widget.addTab(self.plot_widget, 'Plot')
        
        # Adds the main widget to the center of the window
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)
        
        #initializes the widget for setting mesh geometry, inital conditions, and boundary conditions
        self.meshEdit = MeshEditor()
        
        self.meshGenerate = MeshGenerator()
        self.meshGenerate.meshGenerated.connect(self.setMesh)
        
        # tabChanger will control some parameter as different tabs are used in the main widget
        self.main_widget.currentChanged.connect(self.tabChanger) 
        
        

        # setting up the plot controls
        self.plotcontrols = PlotController(self, self.mainPlot)
        self.plotcontrols.plotSettingsUpdated.connect(self.updatePlot) # allows the graph to be updated when settings are changed
        
        
        # This adds the plotcontrol widgets to a dock window as a container
        # This allows the plotcontrols to be convieniently placed to the side of the application window
        self.controlsDock = QtWidgets.QDockWidget('Plot Controls', parent = self)
        self.controlsDock.setWidget(self.plotcontrols)
        self.dockControls()
        
        #setting up controls for running phases
        self.phasesControls = PhasesController(self)
        self.phasesControls.executed.connect(self.phasesExec)
        
        
        # Also adding the phasescontroller to a dockWidget
        # This is set to hide by default and showld be shown by calling phasesDock.show()
        # once a valid mesh has been created or loaded. This functionality is implemented
        # in the dockPhases method within this class
        self.phasesDock = QtWidgets.QDockWidget('Parameters', parent = self)
        scrollarea = QtWidgets.QScrollArea(parent = self)
        scrollarea.setWidget(self.phasesControls)
        self.phasesDock.setWidget(scrollarea)
        self.phasesDock.hide()


        # Setting up the menu bar
        # This is a great way of giving the user access to many different functions
        # in a small area with a simple layout
        self.file_menu = QtWidgets.QMenu('&File', self)
        self.file_menu.addAction('&Quit', self.fileQuit,
                                 QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)
       
        self.mesh_menu = QtWidgets.QMenu('&Mesh', self)
        self.mesh_menu.addAction('&Load Mesh From File', self.loadMeshFromFile)
        self.mesh_menu.addAction('Generate Mesh', self.openMeshGenerate)
        self.mesh_menu.addSeparator()
        self.mesh_menu.addAction('Edit Mesh Conditions', self.openMeshEdit)
        self.mesh_menu.addSeparator()
        self.mesh_menu.addAction('Save Mesh To File', self.saveMeshToFile)
        
        

        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.mesh_menu)

        self.help_menu = QtWidgets.QMenu('&Help', self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.help_menu)
        self.help_menu.addAction('&About', self.about)

        self.plot_menu = QtWidgets.QMenu('&Plot',self)
        self.plot_menu.addAction('Plot Settings', self.dockControls)
        self.plot_menu.addAction('Execute Phases', self.phasesExec)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.plot_menu)

        # a container to the datafiles which may be used to create the mesh
        self.dataFiles ={
                         'Project': '',
                         'Mesh': '',
                         'Output':'',
                         'Boundary Conditions':'',
                         'Initial Conditions':'' 
                         }
                         
    # This function is used to pass updated plotSettings to the matplotlib widget
    # The plotController does not have a reference to the mesh and the data assoxiated with it
    # so the signal needs to be passed here first so the mesh can be passed along 
    # with the plot settings          
    def updatePlot(self, plotType, timeStep, doNodes):
        
        """This method is used to pass updated plotSettings to the matplotlib widget
    The plotController does not have a reference to the mesh and the data assoxiated with it
    so the signal needs to be passed here first so the mesh can be passed along 
    with the plot settings.
    
            plotType -  a member of the mplCanvas.PlotType instance
            timeStep -  an integer specifiying the timeStep which is to be graphed
    
    
    """
        if self.mesh is not None:
            self.mainPlot.update_figure(plotType, timeStep, self.mesh, doNodes)

    def dockControls(self):
        """A method for docking the plotcontrols to the left side of the window.
        If they are already docked there this should do nothing"""
        
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea,self.controlsDock)
        
        
    def dockPhases(self,parameters):
        """
        A method for docking the phasescontrols to the right side of the window.
        If they are already docked there this should do nothing.
        
            parameters - specifies initial parameters which will populate
                            the phasesController param input text boxes initally
        """
        self.phasesControls.setAllParams(parameters)
        self.phasesDock.show()
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea,self.phasesDock)

    def phasesExec(self, params = None):
        """
        This is what tells the mesh to send itself to phases, the C++ library 
        which will then perform all of the number crunching. The results of the simulation will be stored in the mesh class
        After the simulation is completed, the graph widget is reset.
        
            steps - the number of time steps which are to be simulated.
                    If not specified defaults the the global variable stepsDefault.
                    
            params - the parameters to run the simulation with. If not specified, will use the parameters already contained in mesh.
        """
        
        
        if params is not None:
            self.mesh.params = params
        
        self.mesh.sendToAdda()
  
        self.plotcontrols.reset(self.mesh.params['tstp'])
        

    def about(self):
        """Short dialog explaining what the application is"""
        QtWidgets.QMessageBox.about(self, "About", """This is a prototype GUI Application for PHASES""")

    #TODO: Implement file-checking to make sure the file formatting is correct
    def loadMeshFromFile(self):
        """Opens a dialog prompting the user to specify project files that specify a mesh and its initial conditions.
        A mesh will be created immediately form the specified data files.
        File validation has not yet been implemented so improperly formatted files will likely 
        crash the application or case significant errors.
        """
        
        dialog = fileInputDialog(self, self.dataFiles)
       
      
        newFiles = dialog.execute()
        if newFiles is None:
            return
            
        self.plotcontrols.reset()     
        self.dataFiles = newFiles    
           
        self.setMesh(Mesh.loadMeshFromFile(self.dataFiles))
        
    def setMesh(self,mesh):   # displays the mesh parameters to the user which they can then edit
        self.mesh = mesh
        self.dockPhases(self.mesh.params)
        self.meshEdit.setMesh(self.mesh)
        
    def saveMeshToFile(self):
        dialog = QtWidgets.QInputDialog(self)
        
        directory,ok = dialog.getText(self, 'Input Directory', 'Directory:')
        if ok:
            name,ok = dialog.getText(self, 'Name the project', 'Name:')
            if ok:
                self.mesh.params=self.phasesControls.getAllParams()
                Mesh.saveMeshToFile(self.mesh,directory,name)
   
    def openMeshEdit(self):
        """A method to open the mesh generation widget in a new tab 
        within the main application widget.
        """
        self.main_widget.addTab(self.meshEdit, 'Initial Conditions')
            
        self.main_widget.setCurrentWidget(self.meshEdit)
                 
        
    def openMeshGenerate(self):
        """A method to open the mesh generation widget in a new tab 
        within the main application widget.
        """
        self.main_widget.addTab(self.meshGenerate, 'Mesh Generation')
            
        self.main_widget.setCurrentWidget(self.meshGenerate)
        
    def tabChanger(self, index):
        """
        Stops the graph widget from animating while other widgets are active.
        The animation causes some slow down which is annoying if the user is trying to do other things
        
            index - the index of the tab that now has focus
            
        """
        if index != 0:
            self.plotcontrols.canAnimate = False
           # self.meshEdit.drawMesh(self.mesh)
        else:
            self.plotcontrols.canAnimate = True
            
    def message(self, msg, title = 'Error!'):
        
        """A method for creating popup boxes to display error messages.
        
            msg - string which will be displayed to the user.
            title - title of popup, defaults to 'Error!'
              
        """
        
        mb = QtWidgets.QMessageBox(parent = self)
        mb.setText(msg)
        mb.setWindowTitle(title)
        mb.exec()
        
    def fileQuit(self):
        self.close()

    def closeEvent(self,ce):
        self.fileQuit()


icon = os.path.join(os.path.dirname(__file__),'Resources/Icon.png')

qApp = QtWidgets.QApplication(sys.argv)

aw = ApplicationWindow()

aw.setWindowTitle("PHASES")

aw.setGeometry(100,100,1180,820)


qApp.setWindowIcon(QtGui.QIcon(icon))


aw.show()


sys.exit(qApp.exec_())


