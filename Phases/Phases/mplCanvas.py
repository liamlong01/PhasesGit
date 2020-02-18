#tThis class performs the necessary initilizations to embed a matplotlib graph in qt

from matplotlib.backends import qt_compat


from matplotlib import cm

from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import style

from PyQt5 import QtWidgets, QtCore

#import graphs


style.use('ggplot')

class MplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        
        self.axes = self.fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        #self.axes.hold(False)

        self.compute_initial_figure()

        #
        FigureCanvas.__init__(self, self.fig)
        
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self):
        pass


class MyStaticMplCanvas(MplCanvas):
    """Simple canvas with a sine plot."""

    def compute_initial_figure(self):
        t = arange(0.0, 3.0, 0.01)
        s = sin(2*pi*t)
        self.axes.plot(t, s)


class DynamicCanvas(MplCanvas):
    """A canvas that can updated from calls through the methond update_figure"""

    def __init__(self, parent, *args, **kwargs):
        MplCanvas.__init__(self, parent, *args, **kwargs)
     

        #Default settings
        self.index = 1
        self.initializePlotTypes()
        
    def initializePlotTypes(self):
        self.plotTypes = []
        self.plotTypes.append(Default(self))
        self.plotTypes.append(tempContour(self))
        self.plotTypes.append(isothermalLines(self))
        self.plotTypes.append(fliqContour(self))
        self.plotTypes.append(pressureContour(self))
        self.plotTypes.append(entropyContour(self))
        self.plotTypes.append(velVectors(self))
       
        self.plotTypes[-1].addContourPlots(self.plotTypes[1:-1])
            
        
    def update_figure(self, plotType,timeStep, mesh, doNodes):

        self.fig.clear()
        self.axes = self.fig.add_subplot(111)

        
     
        X, Y = mesh.getXY()

        
        plotType.plotFunction(mesh,timeStep)

        if doNodes:
            #X, Y = mesh.getXY()
            self.axes.scatter(X,Y)
        
        self.draw()

        try:     
            fig.tight_layout()
        except:
            pass

        
class PlotType(object):
    """
    An abstarct class that can be used by other classes to implement plot functions.
    
        label - a string label to describe the type of plot
        canvas - the associated graph widget
        plotFunction - the function to be used to plot
    """

    def __init__(self, canvas):
        self.label = None
        self.canvas = canvas
    def plotFunction(self,*args,**kwargs):
        pass
    def settings(self):
        """
        This method returns a QWidget that can allow user inputs that control the formatting of the graph.
        This method is not required to be overridden, an empty widget will be returned if it is not.
        
        Other places in the application are responsible for calling this method and formatting the widget
        in the application window.
        """
        return QtWidgets.QWidget()
    
    
"""
The following classes inherit from PlotType. They can read data from a mesh type object
and output it to an mplCanvas. The formatting of the graphs is contained in the graphs.py file.
"""
    
class Default(PlotType):
    
    def __init__(self,canvas):
        super().__init__(canvas)
        self.label = 'Select Type of Plot'
        
class ColorMapPlot(PlotType):
    """
    Implements the ability to select colormaps for any plot that requires them.
    
    Should be inherited by any PlotType which requires the user to select a colormap
    """
    
    def __init__(self,canvas):
        super().__init__(canvas)
        self.label = 'Select Type of Plot'
        self.colormaps = {'Bone':cm.bone_r, 'Hot':cm.hot_r,'Jet':cm.jet,'Blues':cm.Blues, 'Blue-Purple':cm.BuPu, 'Cool-Warm':cm.coolwarm}
        self.reverseColormaps = {'Bone':cm.bone, 'Hot':cm.hot,'Jet':cm.jet_r,'Blues':cm.Blues_r, 'Blue-Purple':cm.BuPu_r, 'Cool-Warm':cm.coolwarm_r}
        self.cm = 'Jet'
        
    def settings(self):
        widget = QtWidgets.QGroupBox('ColorMap Plot Settings')
        layout = QtWidgets.QFormLayout()
        
        self.cmSelect = QtWidgets.QComboBox()
     
        self.cmSelect.addItems(self.colormaps.keys())
        layout.addRow('Set Colormap :',self.cmSelect)
        self.cmSelect.setCurrentText('Jet')
        
        self.reverse = QtWidgets.QCheckBox()
        layout.addRow('Reverse Colormap :',self.reverse)
        
        widget.setLayout(layout)
        return widget
        
    def getCM(self):
        if not self.reverse.isChecked():
            return self.colormaps[self.cmSelect.currentText()] 
        else:
           return self.reverseColormaps[self.cmSelect.currentText()]
    
 
class ContourPlot(ColorMapPlot):
    
    """
    A class that implements the settings method of the PlotType class for contour plots.
    
    Provides generic settings such as the colormap selection and setting max/min limits.
    
    This should be considered an abstract class that can be inherited to then implement additional
    plotting functionality, namely the PlotType function plotFunction.
    """
    
    def __init__(self,canvas, label):
        super().__init__(canvas)
        self.label = label
        self.numContours = 5
         
    def settings(self):
        widget = super().settings()
        layout = widget.layout()
        widget.setTitle('Contour Plot Settings')
        
        contourSelect = QtWidgets.QSpinBox()
        contourSelect.setMinimum(1)
        contourSelect.setMaximum(500)
     
        contourSelect.setValue(self.numContours)
        contourSelect.valueChanged.connect(self.setNumContours)
        
        layout.addRow('Number of Contours :', contourSelect)
       
       # this is a widget that can be used to add on a velocity vector field
       #  although a dditional functionality for this would need to be added similar to
       # that which is used for adding contours to vector fields
       
       # addVectorField = QtWidgets.QComboBox()
       # addVectorField.addItem('')
       # addVectorField.addItem('Velocity Vector Field')
        
        #layout.addRow('Add Vector Field :', addVectorField)
  
        return widget


    def setNumContours(self, toWhat):
        self.numContours = toWhat
        
    
class VectorField(ColorMapPlot):
    """
    This should be considered an abstract class that can be inherited to then implement additional
    plotting functionality, namely the PlotType function plotFunction.
    """
    
    def __init__(self,canvas,label):
        super().__init__(canvas)
        self.label = label
        self.scale = 9
        self.headlength = 6
        
    def settings(self):
        self.widget = super().settings()
        layout = self.widget.layout()
        self.widget.setTitle('Vector Field Settings')
        
        scaleSelect = QtWidgets.QSlider()
        scaleSelect.setMinimum(1)
        scaleSelect.setMaximum(100)
        scaleSelect.setOrientation(QtCore.Qt.Horizontal)
        scaleSelect.setValue(125 - self.scale*5)
        scaleSelect.valueChanged.connect(self.setScale)
        
        headLength = QtWidgets.QSlider()
        headLength.setMinimum(1)
        headLength.setMaximum(10)
        headLength.setOrientation(QtCore.Qt.Horizontal)
        headLength.setValue(self.headlength)
        headLength.valueChanged.connect(self.setHeadlength)
        
        self.addContour = QtWidgets.QComboBox()
        self.addContour.addItem('No Contour')
        self.addContour.addItems(self.contours.keys())
        self.addContour.currentTextChanged.connect(self.setContourWidget)
        
      
        
        layout.addRow('Vector Length :',scaleSelect)
        layout.addRow('Vector Head Length :',headLength)
        layout.addRow('Add Contour Plot :',self.addContour)
                
        self.contourRow = layout.rowCount() # gets the index number of the next row
                                            # this will be used to place a contour settings widget
                                            # this widget will need be changable thus we need the row number
                                            # doing it this way  means we do not need to no anything about
                                            # widgets added to the layout in super classes
          
        blank = QtWidgets.QWidget()                                    
        layout.addRow(blank)  # a placeholder widget which may later be replaced by contour setti ngs widget
        
        return self.widget
        
    def setScale(self, scale):
        self.scale = (125-scale)/5.0

    def setHeadlength(self, length):
        self.headlength = length
    
    def addContourPlots(self, contourplots):
        self.contours = {}
        for plot in contourplots:
            self.contours[plot.label] = plot
    def setContourWidget(self,label):
        layout = self.widget.layout()
        newWidget = self.contours[label].settings()
        oldWidget = layout.itemAt(self.contourRow,QtWidgets.QFormLayout.SpanningRole).widget()
        
        layout.removeWidget(oldWidget)
        oldWidget.deleteLater()
        del oldWidget
    
        layout.setWidget(self.contourRow, QtWidgets.QFormLayout.SpanningRole,newWidget)
        
    def plotContour (self,mesh,timeStep):
         key = self.addContour.currentText()
         if key  != 'No Contour':
             contour = self.contours[key]
             contour.plotFunction(mesh,timeStep)
            
        
    
class tempContour(ContourPlot):
    
    def __init__(self,canvas, label = 'Temperature Contours'):
        super().__init__(canvas, label)
        self.units = {'Kelvin' : (self.kelvin,'K'), 'Celsius': (self.celsius,'C'), 'Fahrenheit': (self.fahrenheit,'F')}
        
    def plotFunction(self, mesh,timeStep):
        X,Y = mesh.getXY()
        T = mesh.getTemperature(timeStep)
        unitConverter,unit = self.getUnits()
        T= unitConverter(T)
        V = [i*(unitConverter(mesh.tempRange[1])-unitConverter(mesh.tempRange[0]))/self.numContours+unitConverter(mesh.tempRange[0]) 
                for i in range(self.numContours+1)]
        self.canvas.axes.set_title(' Temperature Contours')
        cont=self.canvas.axes.contourf(X,Y,T, V, cmap = self.getCM())
        
        divider=make_axes_locatable(self.canvas.axes)
        clbr = self.canvas.fig.colorbar(mappable=cont)
        clbr.set_label(r'Temperature Reference $^\circ %s$' % unit) 
        
    def settings(self):
        widget = super().settings() # gets the widget from the superclass Contour
                                    # this is done so that an additional widget can be added to control the units of the plot
        widget.setTitle('Temperature Contour Settings')
        
        layout = widget.layout()

        self.unitSelect = QtWidgets.QComboBox()
        self.unitSelect.addItems(self.units)
        self.unitSelect.setCurrentText('Celsius')
        
        layout.addRow('Select Units :', self.unitSelect)
    
        return widget
        
    def getUnits(self):
        return self.units[self.unitSelect.currentText()]
                          
    def kelvin(self, degInKelvin):
        return degInKelvin
        
    def celsius(self, degInKelvin):
        return degInKelvin - 273.15
        
    def fahrenheit(self, degInKelvin):
        degInC = self.celsius(degInKelvin)
        
        return degInC*(1.8) + 32
        
class isothermalLines(tempContour):
    
    def __init__(self,canvas, label = 'Isothermal Lines'):
        super().__init__(canvas, label = label)
  
        
    def plotFunction(self, mesh,timeStep):
        X,Y = mesh.getXY()
        T = mesh.getTemperature(timeStep)
        unitConverter,unit = self.getUnits()
        T= unitConverter(T)
        V = [i*(unitConverter(mesh.tempRange[1])-unitConverter(mesh.tempRange[0]))/self.numContours+unitConverter(mesh.tempRange[0]) 
                for i in range(self.numContours+1)]
        self.canvas.axes.set_title(' Temperature Contours')
        cont=self.canvas.axes.contour(X,Y,T, V, cmap = self.getCM())
        self.canvas.axes.clabel(cont,inline=1) 
         
    
   
        
class fliqContour(ContourPlot):
    
    def __init__(self,canvas):
        super().__init__(canvas, 'Liquid Fraction Contour')
        
        
    def plotFunction(self, mesh,timeStep):
        X,Y = mesh.getXY()
        FLIQ = mesh.getFLIQ(timeStep)
        V = [i*(mesh.fliqRange[1]-mesh.fliqRange[0])/self.numContours + mesh.fliqRange[0]
                for i in range(self.numContours+1)]
                    
        self.canvas.axes.set_title('Liquid Fraction Contours')
        cont=self.canvas.axes.contourf(X,Y,FLIQ, V, cmap = self.getCM())
       
        divider=make_axes_locatable(self.canvas.axes)
        clbr = self.canvas.fig.colorbar(mappable=cont)
        clbr.set_label(r'Liquid Fraction Reference') 
       
        
        
class pressureContour(ContourPlot):
    
    def __init__(self,canvas):
        super().__init__(canvas,'Pressure Contour')
  
        
    def plotFunction(self, mesh,timeStep):
        X,Y = mesh.getXY()
        P = mesh.getPressure(timeStep)
        V = [i*(mesh.pressureRange[1]-mesh.pressureRange[0])/self.numContours + mesh.pressureRange[0]
                for i in range(self.numContours+1)]
                     
        self.canvas.axes.set_title('Pressure Contours')
        cont=self.canvas.axes.contourf(X,Y,P, V, cmap = self.getCM())
        
        divider=make_axes_locatable(self.canvas.axes)
        clbr = self.canvas.fig.colorbar(mappable=cont)
        clbr.set_label(r'Pressure Reference')


class entropyContour(ContourPlot):

    def __init__(self, canvas):
        super().__init__(canvas, 'Entropy Contour')

    def plotFunction(self, mesh, timeStep):
        X, Y = mesh.getXY()
        E = mesh.getEntropy(timeStep)


        self.canvas.axes.set_title('Entropy Contours')
        cont = self.canvas.axes.contourf(X, Y, E, cmap=self.getCM())

        divider = make_axes_locatable(self.canvas.axes)
        clbr = self.canvas.fig.colorbar(mappable=cont)
        clbr.set_label(r'Entropy Reference')

class velVectors(VectorField):
    
    def __init__(self,canvas):
        super().__init__(canvas, 'Velocity Vector Field')
        
    def plotFunction(self, mesh,timeStep):
        self.plotContour(mesh,timeStep)
        
        X,Y = mesh.getXY()
        U,V = mesh.getUV(timeStep)
        
        
        np.seterr(divide='ignore',invalid='ignore')
        magn=np.sqrt((U**2)+(V**2))
        self.canvas.axes.set_ylabel('Hot Wall',fontsize=16)
        #axes.set_xlim([0.0,0.01]) only for mesh2
        #axes.set_ylim([0.0,0.01])
        self.canvas.axes.set_title(r'Vector Field in a Multiphase Flow Simulation')
        #axes.text(6.465454E-02,.0416875,'Solid Region',fontsize=12)
        scale = mesh.velMagRange[1]/np.amax(Y)*self.scale
        quiv = self.canvas.axes.quiver(X,Y,U,V,magn,cmap=self.getCM(), headlength=self.headlength, pivot = 'middle', clim = mesh.velMagRange,units ='y',scale = scale) 
        
        divider=make_axes_locatable(self.canvas.axes)
        #cax=divider.append_axes('right','3%',pad='3%')
        clbr=self.canvas.fig.colorbar(mappable = quiv)
        clbr.set_label(r"Magnitude Reference ($m/s$)")
        
        

    
            
        