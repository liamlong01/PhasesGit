# -*- coding: utf-8 -*-

from PyQt5 import QtCore, QtWidgets


from Phases.mplCanvas import PlotType

stepDefault = 5

class PlotController(QtWidgets.QWidget):
    
    """
    A widget that allow the user to change to change the plot settings
    Current options include selecting the type of plot, which timestep to plot, and 
    animating the plot with a given time interval bewtween each step.
    
    The signal plotSettingsUpdated will be emitted when the plot is to be updated.
    It will specify the type of plot as an instance of the mplCanvas.PlotType class,
    and int specifying the required timestep.
    
    The signal will be emitted only when the user clicks the 'Apply' button, 
    unless the graph is currently being animated in which case it well be emitted every interval
    as given by the user.
    
    Additionally once animated is selected, animation will only begin after 'Apply' is clicked.
    
    Also important to note is that this class does not care about the implementation of the graphWidget
    so long as it has property called 'plotTypes' that is list of instances of the class PlotType.
    This class can then be easily reimplemented should the graph widget be required to change to extend its
    functionality to things like 3d graphing.
    
    """
    
    plotSettingsUpdated = QtCore.pyqtSignal(PlotType, int, bool)

    def __init__(self, parent, graphWidget):
        """
        parent - a reference to the parent widget
        graphWidget - a reference to the widget for plotting
        """
        QtWidgets.QWidget.__init__(self, parent=parent)

        self.maxIndex = 5 #default value to later by a dynamic thing

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.animation)
 
         #setup for the plot controls
 
        grid = QtWidgets.QGridLayout(self)
            
        self.plotSelect = QtWidgets.QComboBox(self)
        self.graph = graphWidget
        labels = []
        for plot in graphWidget.plotTypes:
            labels.append(plot.label)
        
        self.plotSelect.addItems(labels)
        self.plotSelect.currentIndexChanged.connect(self.plotChanged)
        grid.addWidget(self.plotSelect,4,0,1,2)

        label = QtWidgets.QLabel('Step #:', self)
        self.stepCounter  = QtWidgets.QSpinBox(self)
        self.stepCounter.setRange(1,self.maxIndex)
        grid.addWidget(label,0,0)
        grid.addWidget(self.stepCounter,0,1)

        label = QtWidgets.QLabel('Animate:', self)
        self.isAnimated = QtWidgets.QCheckBox(self)
        grid.addWidget(label,1,0)
        grid.addWidget(self.isAnimated,1,1)
      
        apply = QtWidgets.QPushButton('Apply', parent = self)
        apply.clicked.connect(self.update_controls)
        grid.addWidget(apply,7,1)
       

        label = QtWidgets.QLabel('Animation Interval (ms):', self)
        self.interval  = QtWidgets.QSpinBox(self)
        self.interval.setRange(0,3000)
        self.interval.setValue(1000)
        grid.addWidget(label,2,0)
        grid.addWidget(self.interval,2,1)
        
        label = QtWidgets.QLabel('Plot Nodes :', self)
        self.isNodes = QtWidgets.QCheckBox(self)
        grid.addWidget(label,3,0)
        grid.addWidget(self.isNodes,3,1)
        
        
        grid.addWidget(self.graph.plotTypes[0].settings(),5,0,2,2)

        grid.setRowStretch(8,1) 
       
        self.setLayout(grid)
        
        
        self.canAnimate = True
        self.hide()

    def plotChanged(self, plotIndex):
        """
        When the plot select combobox slection is changed this method is called.
        
        This method will remove the old plot settings widget and replace with a new
        one as specified by the selected PlotType in the graph widget.
        
        Each PlotType should have a method called settings which returns the required widget.
        That settings method is called by this method
        """
        
        newWidget = self.graph.plotTypes[plotIndex].settings()
        oldWidget = self.layout().itemAtPosition(5,0).widget()
        
        self.layout().removeWidget(oldWidget)
        oldWidget.deleteLater()
        del oldWidget
        
        self.layout().addWidget(newWidget,5,0,2,2)

         
    def updateGraph(self):
        for plot in self.graph.plotTypes:
            if self.plotSelect.currentText() == plot.label:
                currentPlot = plot
                break
        else:
            currentPlot = None
        self.plotSettingsUpdated.emit(currentPlot, self.stepCounter.value(), self.isNodes.checkState())
        
    def update_controls(self):

        #Error checking before anything is changed
        plotType = self.plotSelect.currentText()
        if plotType == "Select Type of Plot":
           mb = QtWidgets.QMessageBox(parent = self)
           mb.setText("Graph type has not been selected.")
           mb.exec()
           return None     #exits function call

        try:
            self.updateGraph() 
        except FileNotFoundError:
           mb = QtWidgets.QMessageBox(parent = self)
           mb.setText("Error: Mesh or Output File Not Specified")
           mb.exec()
           return None 
        
        if self.isAnimated.checkState():
            self.timer.start(self.interval.value())
        else:
            self.timer.stop()

    def animation(self):
        
        if not self.isAnimated.checkState():
                self.timer.stop()
                return
                
        if self.canAnimate:

            current = self.stepCounter.value()

            if current == self.maxIndex:
                self.stepCounter.setValue(1)
            else:
                self.stepCounter.setValue(current + 1)

            self.updateGraph()
            
            

    def setMaxIndex(self, toWhat):
        self.maxIndex = toWhat
        self.stepCounter.setRange(1,self.maxIndex)

    def reset(self, steps = stepDefault):
        self.plotSelect.setCurrentIndex(0)
        self.stepCounter.setValue(1)
        self.setMaxIndex(steps)
        self.timer.stop()
        self.isAnimated.setCheckState(0)
        self.updateGraph()
         

