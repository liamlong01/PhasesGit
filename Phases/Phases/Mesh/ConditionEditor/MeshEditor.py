# -*- coding: utf-8 -*-

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

from . import Tools #. means import from same directory


class MeshEditor(QWidget):
    
    displayToolChanged = pyqtSignal()


    def __init__(self):
        QWidget.__init__(self)

        layout = QVBoxLayout(self)

        self.toolbar=QToolBar(self)
        layout.addWidget(self.toolbar)
        self.canvas = MeshDisplay(parent=self)

        layout.addWidget(self.canvas)
        
        self.setupToolbar()
        
        self.mesh = None
        
    def setupToolbar(self):
        
        self.singleSelect = Tools.SingleSelect(self.toolbar) 
        self.singleSelect.triggered.connect(lambda: self.setSelectTool(self.singleSelect))
        self.displayToolChanged.connect(self.singleSelect.clearSelections)
        self.toolbar.addAction(self.singleSelect)
        
        self.boxSelect = Tools.BoxSelect(self.toolbar)  
        self.boxSelect.triggered.connect(lambda: self.setSelectTool(self.boxSelect))
        self.displayToolChanged.connect(self.boxSelect.clearSelections)
        self.boxSelect.selectionRequest.connect(lambda func: [func(x) for x in self.canvas.displayTool.selectableObjects()])
        
        
        self.toolbar.addAction(self.boxSelect)
        
        self.toolbar.addSeparator()
        
        self.bcTool = Tools.BCEditTool(self.toolbar)
        self.bcTool.triggered.connect(lambda: self.setDisplayTool(self.bcTool))
        self.bcTool.setShortcut(Qt.CTRL + Qt.Key_B)
        self.toolbar.addAction(self.bcTool)
        
        self.icTool = Tools.ICEditTool(self.toolbar)
        self.icTool.triggered.connect(lambda: self.setDisplayTool(self.icTool))
        self.bcTool.setShortcut(Qt.CTRL + Qt.Key_M)
        self.toolbar.addAction(self.icTool)
        
        self.toolbar.addSeparator()
        
        position = QLabel(self.toolbar)
        self.canvas.mouseMoved.connect(position.setText)
        self.toolbar.addWidget(position)
        
        
        
    def setDisplayTool(self, tool):
        if self.canvas.displayTool is not None:
            self.canvas.displayTool.switchFrom()
        self.canvas.displayTool = tool.switchTo()

        self.displayToolChanged.emit()
        self.canvas.update()
        
    def setSelectTool(self, tool):
        if self.canvas.selectTool is not None:
            self.canvas.selectTool.switchFrom()
        self.canvas.selectTool = tool.switchTo()
        
        
        
        
    def setMesh(self, mesh): # defines scale and position for drawing mesh
        """
        Sets this widget's mesh and performs the necessary calculations for when the mesh is changed.
        """
        self.mesh = mesh
        
        self.bcTool.boundaries = mesh.boundaries
        self.bcTool.params = mesh.params
        
        self.icTool.elements = mesh.elements
        self.icTool.params = mesh.params
        
        maxX = 0
        minX = 0
        maxY = 0
        minY = 0
        
        for node in self.mesh.nodes:
            if node.x > maxX:
                maxX=node.x
            elif node.x<minX:
                minX=node.x
            if node.y >maxY:
                maxY=node.y
            elif node.y<minY:
                minY=node.y
                
        scaleX = 400/(maxX-minX)
        scaleY = 400/(maxY-minY)
        self.setScale(scaleX,scaleY)
        
        self.canvas.transX = 203 - (maxX-minX)/2*scaleX 
        self.canvas.transY = -203 - (maxY-minY)/2*scaleY    

         
    def setScale(self, scaleX,scaleY):
        self.canvas.scaleX = scaleX
        self.icTool.scaleX = scaleX #also sets scale for all other tools
        
        
        self.canvas.scaleY = scaleY
        self.icTool.scaleY = scaleY #also sets scale for all other tools
        
       
        
                
            
        
class MeshDisplay(QWidget):
    
    
    mouseMoved = pyqtSignal(str)
    
    def __init__(self,parent = None):
        QWidget.__init__(self,parent=parent)
        self.setMouseTracking(True) # this makes it so that the mousemovedEvent is called at any point when the mouse moves
                                    # by default this only happens when the mouse is pressed
      
        self.displayTool = Tools.DisplayTool()
        self.selectTool = Tools.SelectTool()
        
     
        self.currentObject = None
        
        self.scaleX=1
        self.scaleY=1
        self.transX=0
        self.transY=0 
        
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.showContextMenu)
        

        
    def showContextMenu(self, point): 
        
        #Takes the possible right-click menu actions from the current display tool 
        #The user can then select these actions in a dropdown menu that appears
        menu = QMenu('Context Menu', self)
        dAction = QAction("Default Action", menu)
        dAction.triggered.connect(lambda: print(self.currentObject))
        menu.addAction(dAction)
        for actionSet in self.displayTool.contextMenuActions(menu, self.selectTool.selections):
            menu.addSeparator()
            for action in actionSet:
                menu.addAction(action)
            
        menu.exec(self.mapToGlobal(point))
    
    def paintEvent(self,event):
        
        def drawMesh(mesh, painter):
            
            if isinstance(self.displayTool,Tools.BCEditTool):
                #drawing boundaries different colors
                for boundary in mesh.boundaries:
                    color = boundary.conditions['temperature'][0][2]*255
                    pen = QPen(QColor(color,165,0))
                    pen.setWidth(8)
                    painter.setPen(pen) 
                    painter.setBrush(QColor(color,165,0))
                    points = [QPointF(node.x*self.scaleX,node.y*self.scaleY) for node in boundary.nodes]
                    painter.drawLine(points[0],points[1])
                
            for element in mesh.elements:

                if isinstance(self.displayTool,Tools.ICEditTool):
                    #drawing elements as differnt colors
                    color = element.initialConditions['temperature']*255
    
                    painter.setPen(QColor(255,255-color,255-color)) 
                    painter.setBrush(QColor(255,255-color,255-color))
                    
                    points = [QPointF(node.x*self.scaleX,node.y*self.scaleY) for node in element.nodes]
                    painter.drawPolygon(QPolygonF(points))
            
                

                painter.setBrush(Qt.black)
                painter.setPen(Qt.black)
                for node in element.nodes:
                    x = node.x*self.scaleX
                    y = node.y*self.scaleY
                    
                    #draws a node
                    painter.drawEllipse(QRectF(x-3,y-3,6,6))
                    
                    prevNode = element.nodes[element.nodes.index(node)-1]
                    prevX = prevNode.x*self.scaleX
                    prevY = prevNode.y*self.scaleY
                    
                    #draws aline between nodes
                    painter.drawLine(QPointF(prevX,prevY),QPointF(x,y))
                    
        painter = QPainter()
        painter.begin(self)
        
        # flips y-axis so that up is positive
        # down is positive in qt by default
        painter.scale(1,-1)
        painter.translate(self.transX,self.transY)

        if self.parent().mesh is not None:
           
            painter.setPen(Qt.black) 
            painter.setBrush(Qt.black)
            
            drawMesh(self.parent().mesh, painter)
            
            
            self.selectTool.drawStuff(painter,self.currentObject)
            
            if self.currentObject is not None:
                self.displayTool.drawStuff(painter, self.currentObject)
                
            if self.selectTool.selections: #if there is at least one selection
                for selection in self.selectTool.selections:
                    self.displayTool.drawStuff(painter, selection)
                    
            self.selectTool.drawStuff(painter,self.currentObject)
           
        
            
        painter.end()
        
           
    def mouseMoveEvent(self, event):

        
        x,y = self.getRealXY(event.x(),event.y())

        self.mouseMoved.emit("x: %s" % x + " y: %s" % y)
        
        self.currentObject = self.displayTool.mouseMoveEvent(x,y)
        self.selectTool.mouseMoveEvent(x,y,self.currentObject)
        
        
    def mousePressEvent(self, event):
        
        
        x,y = self.getRealXY(event.x(),event.y())
        if event.button() == Qt.LeftButton:
            self.selectTool.mousePressEvent(x,y,self.currentObject)
            
    def mouseReleaseEvent(self, event):
        x,y = self.getRealXY(event.x(),event.y())
        self.selectTool.mouseReleaseEvent(x,y,self.currentObject)
        self.update()
        
            
        

    def getRealXY(self, qx, qy):
        # calaculates the real xy-coordinates based on the widget's internal coordinate system
        
        x = (round(qx)-self.transX)/self.scaleX
        y = (round(qy)+self.transY)/self.scaleY*-1

        return x,y
        
    
        