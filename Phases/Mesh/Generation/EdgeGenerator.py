# -*- coding: utf-8 -*-

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

from Mesh.Mesh import Node 

from math import atan,pi      

class Edge(object):
    def __init__(self,node1, node2, scaler):
        """
        Each edge contains two nodes.
        Scaler should be a function that can be used to find the scalefactors necessary for painting the nodes on a painting widget.
        """

        self.node1 = node1
        self.node2 = node2
        self.scales = scaler
        self.current = False
        
         
    def scaleNodes(self):
        """
        Takes the edges two nodes, and returns them as two QPointFs appropriately scaled based on the display's scale.
        This is useful for painting and should be called before using nodes to create objects to be used by a QPainter
        """
        
        scaleX,scaleY = self.scales()
        
        point1 =  QPointF(self.node1.x*scaleX,self.node1.y*scaleY) 
        point2 = QPointF(self.node2.x*scaleX,self.node2.y*scaleY)
        
        return point1,point2
        
    def draw(self,painter):
        if self.current:
            painter.setPen(Qt.green)
            painter.setBrush(Qt.green)
            
        else:
            painter.setPen(Qt.black)
            painter.setBrush(Qt.black)
    
    def settings(self):
        return QWidget(),[]
        
    def lerp(self,x1,y1,x2,y2,t):
        """
        linear interpretation between two points
        t - a parameter specifying where to interpolate
                t = 0 is p1
                t=1 is p2
                0<t<1 will give a point on the line in between
        """
        x = (1-t)*x1 + t*x2
        y = (1-t)*y1 + t*y2

        return x,y

        
    def split(self):
        pass
        
        
        
    def __eq__(self, otherEdge):
        if not isinstance(otherEdge,Edge):
            return False

        if self.node1 == otherEdge.node1 and self.node2 == otherEdge.node2:
            return True
        else:
            return False
        
        
        
class StraightLine(Edge):
    
    
    
    def __init__(self, node1, node2,scaler):
        super().__init__(node1,node2,scaler)
        self.detail = 4
        
    def draw(self, painter):
        super().draw(painter)
        point1, point2 = self.scaleNodes()

        
        line = QLineF(point1,point2)
        painter.drawLine(line)
        
        ##drawing extra Nodes
        for t in range(self.detail):
            t = (t+1)/(self.detail+1)
            x,y = self.lerp(point1.x(),point1.y(),point2.x(),point2.y(),t)
            point = QPointF(x,y)
            painter.drawEllipse(QRectF(point.x()-2,point.y()-2,4,4))
            
    def settings(self):
        widget = QWidget()
        detailSlider = QSlider()
        detailSlider.setOrientation(Qt.Horizontal)
        detailSlider.setMinimum(0)
        detailSlider.setMaximum(20)
        
        detailSlider.valueChanged.connect(self.setDetail)
        
        layout = QHBoxLayout(widget)
        layout.addWidget(QLabel('Edge Detail:'))
        layout.addWidget(detailSlider)
        
        signals = [detailSlider.valueChanged]
        
        return widget, signals
        
    def setDetail(self, howMuch):
         self.detail = howMuch
    
    def normal(self):
        pass

    def angle(self, refNode = None):
        if refNode == self.node2:
            #swap the nodes
            temp = self.node1
            self.node1=self.node2
            self.node2 = temp
        
        deltX = self.node1.x-self.node2.x
        deltY =self.node1.y-self.node2.y
        
        angle = atan(deltY/deltY)

        if deltX < 0:
            angle = angle+pi
        
    def split(self):
        
        newNodes = []

        for t in range(self.detail):
            t=(t+1)/(self.detail+1)
            x,y = self.lerp(self.node1.x,self.node1.y,self.node2.x,self.node2.y,t)
            newNodes.append(Node(x,y))
            
        prevNode = self.node1
        newEdges = []
        for node in newNodes:
            line = StraightLine(prevNode,node,self.scales)
            line.detail = 0
            newEdges.append(line)
            prevNode = node
        line = StraightLine(prevNode,self.node2,self.scales)
        line.detail = 0
        newEdges.append(line)

        
        return newEdges,newNodes
            
            
            
class EdgeGenerator(QWidget):
    
    edgeEdited = pyqtSignal(Edge)
    
    def __init__(self, nodes,display,parent = None):
        QWidget.__init__(self, parent)
        
        self.edges = display.edges
        self.nodes = nodes #reference to the nodes that can be used to create edges, should not be edited by this class

        self.edgeTypes = ['Straight Line']
        self.currentEdge = None
        

        ###creating edges
        edgeCreate = QGroupBox('Create Edge')
        layout = QGridLayout(edgeCreate)
        
        self.node1Sel = QSpinBox(self)
        self.node1Sel.setMaximum(0)
        self.node1 = lambda: self.nodes[self.node1Sel.value()-1]
        layout.addWidget(QLabel('Node 1:'), 0,0)
        layout.addWidget(self.node1Sel,0,1)
        
        self.node2Sel = QSpinBox(self)
        self.node2Sel.setMaximum(0)
        self.node2 = lambda: self.nodes[self.node2Sel.value()-1]
        layout.addWidget(QLabel('Node 2:'), 1,0)
        layout.addWidget(self.node2Sel,1,1)
        
        self.edgeType = QComboBox()
        self.edgeType.addItems(self.edgeTypes)
        layout.addWidget(self.edgeType,2,0,1,2)
        
        create = QPushButton('Create Edge')
        create.clicked.connect(self.createEdge)
        layout.addWidget(create,3,1)
        layout.setRowStretch(4,1)
        
        ###editing Edges
        self.edgeEditor = QGroupBox('Edit Edges')
        self.edgeSelect = QSpinBox()
        self.edgeSelect.setMaximum(0)
        self.edgeSelect.valueChanged.connect(self.edgeSelected)
        
        editorLayout = QGridLayout(self.edgeEditor)
        editorLayout.addWidget(QLabel('Select Edge:'), 0,0)
        editorLayout.addWidget(self.edgeSelect,0,1)
        editorLayout.addWidget(QWidget(),1,0)
        
        delButton = QPushButton('Delete Edge')
        delButton.clicked.connect(self.deleteEdge)
        editorLayout.addWidget(delButton,3,1)
        
        
        #Save Edges
        edgeSave = QGroupBox('Save Edges')
        saveLayout = QGridLayout(edgeSave)
        saveButton = QPushButton('Save Edges')
        saveButton.clicked.connect(self.saveEdges)
        saveLayout.addWidget(saveButton)
        
        mainLayout = QVBoxLayout(self)
        mainLayout.addWidget(edgeCreate)
        mainLayout.addWidget(self.edgeEditor)
        mainLayout.addWidget(edgeSave)
        mainLayout.addStretch(3)
        

        self.scaler = lambda: (display.scaleX,display.scaleY)

        
    def createEdge(self):
        if self.node1Sel.value() == self.node2Sel.value() or self.node1Sel.value()==0 or self.node2Sel.value()==0:
            return
            
        if self.edgeType.currentText()=='Straight Line':
            self.addStraightLine(self.node1(),self.node2())
            
      

    def addStraightLine(self,firstnode, secondnode):
        edge = StraightLine(firstnode,secondnode, self.scaler)
            
        if edge in self.edges:
         
            return edge
        
        self.edges.append(edge)
        
        self.edgeSelect.setMaximum(len(self.edges))
        
        self.edgeEdited.emit(edge)

        return edge

        
    def edgeSelected(self, index):
        if index ==0:
            return
            
        if self.currentEdge is not None: 
            self.currentEdge.current = False #the old edge is deselected
        
        
        self.currentEdge = self.edges[index-1] #change to new edge
        self.currentEdge.current = True
        
        layout = self.edgeEditor.layout()
        
        newWidget, signalsToDetect = self.currentEdge.settings()
        oldWidget = layout.itemAtPosition(1,0).widget()
        
        layout.removeWidget(oldWidget)
        oldWidget.deleteLater()
        del oldWidget
        
        layout.addWidget(newWidget,1,0,2,2)
        
        for signal in signalsToDetect:
            signal.connect(lambda:self.edgeEdited.emit(self.currentEdge))
                
    def deleteEdge(self, index):
        edge = self.edges.pop(index-1)
        self.edgeSelect.setMaximum(len(self.edges))
        self.edgeEdited.emit(edge)  
    
        
    def saveEdges(self):
        
        
        oldEdges = self.edges.copy()
        for edge in oldEdges:
            self.edges.pop(self.edges.index(edge))
            newEdges,newNodes = edge.split()
            self.edges.extend(newEdges)
            self.nodes.extend(newNodes)  
            
        self.edgeSelect.setMaximum(len(self.edges))
        
        #finally nodes and edges should be organized counterclockwize
        #starting with the rightmostNode
        firstNode = self.nodes[0]
        for node in self.nodes:
            if node.x > currentNode:
                firstNode = node
         
         #edge with the lowest angle at this node is the first edge        
        if firstNode.edge1.angle(refNode = firstNode) < firstNode.edge2.angle(refNode=firstNode) :
            firstEdge = firstNode.edge1
        else:
            firstEdge = firstNode.edge2
       
        
        nextNode = None
        currentNode = firstNode
        ccwNodes = []
        while nextNode != firstNode:
            ccw.Nodes.append(currentNode)
            
        
            
        
    def nodeDeleted(self,node):
        edgesToDel = []
        for edge in self.edges:
            
            if node == edge.node1 or node == edge.node2:
                edgesToDel.append(edge)
         
        for edge in edgesToDel:
            self.edges.pop(self.edges.index(edge))
            self.edgeEdited.emit(edge)

        
        
    def nodeEdited(self):
        self.node1Sel.setMaximum(len(self.nodes))
        self.node2Sel.setMaximum(len(self.nodes))
        
        


       
   
        


            
        
        
        
        
