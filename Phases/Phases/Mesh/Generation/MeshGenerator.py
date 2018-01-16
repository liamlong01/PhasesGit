from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *


from .NodeGenerator import NodeGenerator
from .EdgeGenerator import EdgeGenerator
from .RectGenerator import RectGenerator
#from .MeshCreator import MeshCreator

from math import log

from Phases.Mesh.Mesh import Mesh  

class MeshGenerator(QWidget):

    meshGenerated = pyqtSignal(Mesh)

    def __init__(self, parent = None):
        QWidget.__init__(self, parent)
        
        self.display = GeneratorDisplay(parent = self)
        controlWidget = QTabWidget()

        
        
         #node editor
        nodeEdit = NodeGenerator(self.display,parent = self)
        controlWidget.addTab(nodeEdit,'Nodes')
        
        #edge editor
        edgeEdit = EdgeGenerator(nodeEdit.nodes, self.display,parent = self)
        controlWidget.addTab(edgeEdit,'Edges')

        #mesh editor
        rectCreate = RectGenerator(nodeEdit, edgeEdit, self.display, parent = self)
        rectCreate.meshCreated.connect(self.meshGenerated.emit)
        controlWidget.addTab(rectCreate,'Rectangle')

        
        #mesh editor
        meshCreate = QWidget()
        controlWidget.addTab(meshCreate,'Mesh')
        
        #setting up the main display
        
        nodeEdit.nodeEdited.connect(self.display.nodeEdited)
        nodeEdit.nodeEdited.connect(edgeEdit.nodeEdited)
        nodeEdit.nodeDeleted.connect(edgeEdit.nodeDeleted)
        edgeEdit.edgeEdited.connect(self.display.edgeEdited)
        self.display.mouseMoved.connect(nodeEdit.position.setText)
        
        
        #creating the layout
        layout = QHBoxLayout(self)
        layout.addWidget(controlWidget)
        layout.addWidget(self.display)
        layout.setStretch(1, 4)

   
            

        
class GeneratorDisplay(QWidget):
    
    mouseMoved = pyqtSignal(str)
    
    def __init__(self,parent = None):
        QWidget.__init__(self, parent)
        self.setMouseTracking(True) #track mouse position without mouse button pressed
        
        self.nodes = []
        self.edges=[]

        self.transX = 0
        self.transY = 0
        self.scaleX = 1
        self.scaleY = 1
        
    def paintEvent(self, event):
        
        
        painter = QPainter()
        painter.begin(self)
        
        painter.scale(1,-1)
        painter.translate(self.transX,self.transY)
        
        painter.setBrush(Qt.black)
        painter.setPen(Qt.black)
        
        for node in self.nodes:
            x = node.x*self.scaleX
            y = node.y*self.scaleY
            painter.drawEllipse(QRectF(x-3,y-3,6,6))
            
            #drawing number labels for each node
            #the axes need to be filpped twice so that the text
            #is oriented in the correct direction
            
            ndigits = round(log(self.nodes.index(node)+1))+ 1 #number of digits in the number to be placed
                                                                #if the number has multiple digits it needs to be moved over more
      
            painter.scale(1,-1)
            painter.drawText(x-8*ndigits,-y+2, str(self.nodes.index(node)+1))
            painter.scale(1,-1)
            
        for edge in self.edges:
            edge.draw(painter)
        
        painter.end()
        
    def mouseMoveEvent(self,event):
        
        
        self.mouseMoved.emit('x : %s y : %s'%self.getRealXY(event.x(),event.y()))
        
        
    def nodeEdited(self):
        if not self.nodes:
            return
        
        self.update()
        
        maxX = max([node.x for node in self.nodes])
        minX = min([node.x for node in self.nodes])
        deltX = maxX-minX
        
        maxY = max([node.y for node in self.nodes])
        minY = min([node.y for node in self.nodes])
        deltY=maxY-minY
        
        if not maxX == minX:
            self.scaleX = 600/(deltX)
        
       #     if maxX==0:
        #        deltX = 0.00001
         #   else:
          #      deltX = maxX
        if not maxY==minY:
            self.scaleY = 600/(deltY)
            
       #     if maxY==0:
        #        deltY = 0.00001
         #   else:
          #      deltY = maxY
                  
        
        
        
        
        self.transX = 307- ((deltX)/2+minX)*self.scaleX 
        self.transY = -307 - ((deltY)/2+minY)*self.scaleY  

    def edgeEdited(self,edge):
        self.update()
        
        

    def getRealXY(self, qx, qy):
        
        ndigitsX = round(log(self.scaleX/400,10)) + 2
        ndigitsY = round(log(self.scaleY/400,10)) + 2
        
        x = round((round(qx)-self.transX)/self.scaleX,ndigitsX)
        y = round((round(qy)+self.transY)/self.scaleY*-1,ndigitsY)

        return x,y
        
        
        