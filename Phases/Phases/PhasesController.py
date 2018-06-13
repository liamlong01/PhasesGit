# -*- coding: utf-8 -*-

from PyQt5 import  QtCore, QtWidgets

# this dictionary contains all of the tooltips for each parameter used by phases  
# these will be used by the PhasesController class to describe them to the user          
paramDescriptions = {
       'se': """ There are 3 main solution blocks in cntrl, C(concentration equation), T (temperature) and U/V/P (velocities/ pressures),
se indicates which block is solved, numbered 1,2,3 or 4 (all together), 
Options: se = 1(C), 2(T), 3(UVP), 12(CT), 13(CUVP), 23(TUVP), 4(CTUVP)""",

       'ao': 'An addition capability of droplet flows was added later, ao = 2 means that the droplet flow solver is added, use ao = 1 (without droplets)',
       'tk': 'Number of times cycling through the temperature equation loop until a convergence tolerance is reached',
       'tvk': 'Same as tk, but instead the number of time through the temperature/velocity loop',
       'tstp': 'Number of time steps',
       'df':'Size of timestep',
       'lr':'Domain height',
       'wr':'Domain width',
       'tol': 'Residual tolerance',
       'ur':'Boundary velocity reference',
       'pcs' : 'Concentrarion of solidus',
       'pcl': 'Concentration of liquidus',
       'tmlt': 'Melting point',
       'tsol' : 'Temperature of solidus',
       'bt': 'Thermal expansivity',
       'bs': 'Sollutal expansivity',
       'tmin': 'Minimum temperature',
       'tmax': 'Maximum temperature',
       'vsc': 'Kinematic viscosity',
       'prs' : 'Solid density',
       'prl' : 'Fluid density',
       'pks' : 'Solid conductivity',
       'pkl' : 'Liquid Conductivity',
       'pds0' : 'Solid diffusivity',
       'pdl0' : 'Fluid diffusivity',
       'phs' : 'Solid specific heat',
       'phl' : 'Fluid specific heat',
       'pl' : 'Latent heat of fusion',
        'k' : 'Thermal Conductivity'
    }

class PhasesController (QtWidgets.QWidget):
    
    """
    This class is used to control the parameters that are used during a phases simulation.
    It contains a QLineEdit widget for each parameter. Once the exucute button is pressed, 
    a dialog will pop up asking for the number of timesteps, and the executed signal will be emitted.
    This signal will contain an integer that was inputted through the dialog by the user, 
    and a dictionary of all the currently specified parameters
    
    """
    
    executed = QtCore.pyqtSignal(dict)
    
       
    def __init__(self, parent, initParams = None): #if initParams is none need to populate with all params equal to none
        QtWidgets.QWidget.__init__(self, parent=parent)
        
        #dictionaries are unordered so to have a consistent an order of parameters needs to be defined
        self.orderOfParams = ['se','ao','tk','tvk','tstp','df','lr','wr','tol','ur','pcs','pcl','tmlt','tsol','bt','bs','tmin','tmax','vsc','prs','prl','pds0','pdl0','pks','pkl','phs','phl','pl', 'k']
        
        self.layout = QtWidgets.QFormLayout(self)
        self.inputs = {}
        for key in self.orderOfParams:
            if initParams is not None:
               pInput = paramInput(key,str(initParams[key]))
            else:
               pInput = paramInput(key,'')
            #pInput.paramChanged.connect(self.setParam)
            self.inputs[key] = pInput
            
            self.layout.addRow('%s:' %key, pInput)
        
        
        testButton = QtWidgets.QPushButton('PHASES', self)
        testButton.clicked.connect(self.test) 
        self.layout.addRow('Execute:', testButton)

        self.setLayout(self.layout)
        
        
    def getParam(self,paramToGet):
        return float(self.inputs[paramToGet].text())
        
    def getAllParams(self):
        params = {}
        for key in self.inputs.keys():
            params[key] = self.getParam(key)
        return params
            
        
    def setParam(self,param, toWhat):
        self.inputs[param].setText(str(toWhat))
       

    def setAllParams(self,params):
        for key in params.keys():
            self.setParam(key,params[key])
           
            

    def test(self):
       # dialog = QtWidgets.QInputDialog(self)
        #steps, ok= dialog.getInt(self, 'Input Number of Time Steps', 'Time Steps:')
        #if ok:
        self.executed.emit(self.getAllParams())

class paramInput(QtWidgets.QLineEdit):

    paramChanged = QtCore.pyqtSignal(str,float)
    
    """
    An extension of the QLineEdit class to add some additional functionality for convenience.
    Takes two strings as input to the __init__ method.
    
        param - a description of the param the instance will represent
        initVal - the value to place in the textbox initially.
    
    
    When text is changed, emits a signal called paramChnaged which contains  a string describing
    the param that was changed, and a float of which value it was changed to. The description comes from
    the value param passed in the __init__ method.
    
    
    """

    def __init__(self, param, initVal):
        QtWidgets.QLineEdit.__init__(self,initVal)
        self.textChanged.connect(self.paramEdited)
        self.setToolTip(paramDescriptions[param])
        self.param = param

    def paramEdited(self,toWhat):
        """
        The method that is called whenever the string in the text box changes.
        Only allows strings that represent floats
        Or the empty string, so that editing values is smooth
        """
        
        try:
            newParam = float(toWhat)
        except ValueError: 
            # the param has been changed to something that is not a valid float
            # the only valid string for this is the empty string 
            # if it is anything else we set it to the empty string
            if toWhat != '':       
               self.setText('') 
            return
            
        self.paramChanged.emit(self.param, newParam)