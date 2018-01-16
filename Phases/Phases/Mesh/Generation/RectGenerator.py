# -*- coding: utf-8 -*-

from PyQt5.QtCore import *
from Mesh import Mesh
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

import numpy as np

class Rectangle(object):
  

    def __init__(self, height, width, nodes,edges):
        self.nodes = nodes
        self.height = height
        self.width = width
        self.edges = edges
        self.topEdges = [edges[0]]
        self.leftEdges = [edges[1]]
        self.botEdges = [edges[2]]
        self.rightEdges = [edges[3]]
        
    def removeEdge(self ,edge):

        if edge in self.rightEdges:
            self.rightEdges.remove(edge)
        elif edge in self.leftEdges:
            self.leftEdges.remove(edge)
        elif edge in self.topEdges:
            self.topEdges.remove(edge)
        elif edge in self.botEdges:
            self.botEdges.remove(edge)

        self.edges.remove(edge)
        
    def delete(self):
        for node in self.nodes:
            pass



class RectGenerator(QWidget):

    meshCreated = pyqtSignal(Mesh.Mesh)
    
    def __init__(self, nodeEditor, edgeEditor, display, parent = None):
        QWidget.__init__(self, parent)

        self.nodeEditor = nodeEditor
        self.edgeEditor = edgeEditor
        self.display = display

        self.rectangles = []

        layout = QVBoxLayout(self)

        #define element size
        res = QGroupBox('Resolution')
        resLayout = QGridLayout(res)

        dxEdit = QLineEdit()
        resLayout.addWidget(QLabel('dx :'),0,0)
        resLayout.addWidget(dxEdit,0,1)
        self.dx = lambda: float(dxEdit.text())

        dyEdit = QLineEdit()
        resLayout.addWidget(QLabel('dy :'),1,0)
        resLayout.addWidget(dyEdit,1,1)
        self.dy = lambda: float(dyEdit.text())


        ###add new square
        newCreate = QGroupBox('Add new square')
        newLayout = QGridLayout(newCreate)

        self.newhEdit = QLineEdit()
        newLayout.addWidget(QLabel('Height :'),0,0)
        newLayout.addWidget(self.newhEdit,0,1)
        heightNew = lambda: float(self.newhEdit.text())

        self.newwEdit = QLineEdit()
        newLayout.addWidget(QLabel('Width :'),1,0)
        newLayout.addWidget(self.newwEdit,1,1)
        widthNew = lambda: float(self.newwEdit.text())
        
        self.alignment = QComboBox()
        self.alignment.addItems(['Middle','Bottom', 'Top'])
        newLayout.addWidget(self.alignment,2,0,1,2)


        addnew = QPushButton('Add')
        newLayout.addWidget(addnew,3,0)
        addnew.clicked.connect(lambda: self.addNewRectangle(heightNew(),widthNew()))
        newLayout.setRowStretch(4,1)
        
        #generate
        generate = QGroupBox('Generate')
        genLayout = QGridLayout(generate)

        genNodes = QPushButton('Generate')
        genLayout.addWidget(genNodes,0,0)
        genNodes.clicked.connect(self.generateNodes)
     

        layout.addWidget(res)
        layout.addWidget(newCreate)
        layout.addWidget(generate)

        layout.addStretch(3)

    def addRectangle(self, refX,refY, height,width):
        node1 = self.nodeEditor.addNode(refX+width,refY+height)
        node2 = self.nodeEditor.addNode(refX,refY+height)
        node3 = self.nodeEditor.addNode(refX,refY)
        node4 = self.nodeEditor.addNode(refX+width,refY)

        edge1 = self.edgeEditor.addStraightLine(node1,node2)
        edge2 = self.edgeEditor.addStraightLine(node2,node3)
        edge3 = self.edgeEditor.addStraightLine(node3,node4)
        edge4 = self.edgeEditor.addStraightLine(node4,node1)

        nodes = [node1, node2, node3, node4]
        edges = [edge1, edge2, edge3, edge4]
        for edge in edges:
            edge.setDetail(0)

        self.rectangles.append(Rectangle(height, width, nodes,edges))

        if len(self.rectangles) > 1:
            self.deleteRectEdge(self.rectangles[-2], self.rectangles[-2].rightEdges[0])
            self.deleteRectEdge(self.rectangles[-1], self.rectangles[-1].leftEdges[0])

            edge = self.edgeEditor.addStraightLine(self.rectangles[-2].nodes[0],self.rectangles[-1].nodes[1])
            edge.setDetail(0)
            if self.rectangles[-2].nodes[0].y > self.rectangles[-1].nodes[1].y:
                self.rectangles[-2].edges.append(edge)
                self.rectangles[-2].rightEdges.append(edge)
            else:
                self.rectangles[-1].edges.append(edge)
                self.rectangles[-1].leftEdges.append(edge)

            edge = self.edgeEditor.addStraightLine(self.rectangles[-2].nodes[3],self.rectangles[-1].nodes[2])
            edge.setDetail(0)
            if self.rectangles[-2].nodes[3].y < self.rectangles[-1].nodes[2].y:
                self.rectangles[-2].edges.append(edge)
                self.rectangles[-2].rightEdges.append(edge)
            else:
                self.rectangles[-1].edges.append(edge)
                self.rectangles[-1].leftEdges.append(edge)
          

    def addNewRectangle(self,height, width):

        #determine if height and width are compatiple with dx and dy
        
        if  self.alignment.currentText() == 'Middle' and self.rectangles!=[]:
            if not (round(height/4/self.dy(),8)).is_integer():
                mb = QMessageBox(parent = self)
                mb.setText('For Middle alignment height/4 must be evenly divisible by dy')
                mb.setWindowTitle('Dimension error')
                mb.exec()
                return

        if not (round(height/self.dy(),8)).is_integer():
            mb = QMessageBox(parent = self)
            mb.setText('Height must be evenly divisible by dy')
            mb.setWindowTitle('Dimension error')
            mb.exec()
            return
        
        if not (round(width/self.dx(),8)).is_integer():
            mb = QMessageBox(parent = self)
            mb.setText('Width must be evenly divisible by dx')
            mb.setWindowTitle('Dimension error')
            mb.exec()
            return


        #first rectangle starts at 0,0
        if self.rectangles == []:
            refX = 0
            refY = 0
        else:
          
            if self.alignment.currentText() == 'Bottom':
                refX = self.rectangles[-1].nodes[3].x #bottom right node provides reference
                refY = self.rectangles[-1].nodes[3].y #bottom right node provides reference

            elif self.alignment.currentText() == 'Middle':
                refX = self.rectangles[-1].nodes[3].x #average top and bottom right node provides reference
                refY = (self.rectangles[-1].nodes[3].y + self.rectangles[-1].nodes[0].y )/2 - height/2 #average bwetween top and bottom provides reference

            elif self.alignment.currentText() == 'Top':
                refX = self.rectangles[-1].nodes[0].x #topright node provides reference
                refY = self.rectangles[-1].nodes[0].y - height #top right node provides reference

        self.addRectangle(refX,refY,height,width)
      


    def deleteRectEdge(self,rect,edge):
        rect.removeEdge(edge)
        self.edgeEditor.deleteEdge(self.edgeEditor.edges.index(edge)+1)

    def generateNodes(self):
        #first lets shift all rectangles sp that all coordinates are positive
        minX = 0
        minY = 0
        for rect in self.rectangles:
            for node in rect.nodes:
                if node.x < minX:
                    minX = node.x
                if node.y < minY:
                    minY = node.y

        for rect in self.rectangles:
            for node in rect.nodes:
                node.x = node.x - minX
                node.y = node.y - minY



        #now generate mesh
        mesh = Mesh.Mesh()
        for rect in self.rectangles:
            nodes =[]
            #making nodes
            nodex = np.arange(rect.nodes[2].x,round(rect.nodes[0].x+self.dx(),8),self.dx())
            nodey = np.arange(rect.nodes[2].y,round(rect.nodes[0].y+self.dy(),8),self.dy())
            for i in range(len(nodex)):
                for j in range(len(nodey)):
                    node = Mesh.Node(nodex[i],nodey[j])

                    #if node is a corner copy the corner node
                    for rnode in rect.nodes:
                        if node == rnode:
                            node = rnode

                    node = mesh.addNode(node)
                    nodes.append(node)
                    
                    
            #making elements
            shift = len(nodey)
            elements = []
            for i in range(len(nodex)-1):
                for j in range(len(nodey)-1):
                    node1 = nodes[shift*i+j+shift+1]
                    node2 = nodes[shift*i+j+1]
                    node3 = nodes[shift*i+j]
                    node4 = nodes[shift*i+j+shift]

                    elemNodes = [node1, node2, node3, node4]
                    elem = Mesh.Element(elemNodes)
                    elements.append(elem)
                    mesh.addElement(elem)
                    elem.initialConditions['temperature'] = 0
                    elem.initialConditions['concentration'] = 0
                    elem.initialConditions['u-velocity'] = 0
                    elem.initialConditions['v-velocity'] = 0

            #defining boundaries
            for edge in rect.edges:


                if edge in rect.topEdges or edge in rect.botEdges:
                    shift = len(nodey)
                    if edge.node1.x > edge.node2.x:
                        startNode = edge.node2
                        endNode  = edge.node1
                    else:
                        startNode = edge.node1
                        endNode = edge.node2

                    if edge in rect.topEdges:
                        side = 'T2'
                    else:
                        side = 'B3'

                elif edge in rect.leftEdges or edge in rect.rightEdges:
                    shift = 1
                    if edge.node1.y > edge.node2.y:
                        startNode = edge.node2
                        endNode  = edge.node1
                    else:
                        startNode = edge.node1
                        endNode = edge.node2
                    if edge in rect.leftEdges:
                        side = 'L1'
                    else:
                        side = 'R4'


                start = nodes.index(startNode)
                end = nodes.index(endNode)
                edgeNodes = np.arange(start,end+shift,shift)
                if edge in rect.leftEdges or edge in rect.topEdges:
                    edgeNodes = np.flipud(edgeNodes)

                for i in range(len(edgeNodes)-1):
                    boundary = Mesh.Boundary([nodes[edgeNodes[i]], nodes[edgeNodes[i+1]]])
                    
                    
                    


                    for elem in elements:
                        if boundary.nodes[0] in elem.nodes:
                            if boundary.nodes[1] in elem.nodes:
                                boundary.element = elem

                    mesh.addBoundary(boundary)
                    boundary.conditions['temperature'] = [[0,1,0],[0,1,0]]
                    boundary.conditions['concentration'] = [[0,1,0],[0,1,0]]
                    boundary.conditions['u-velocity'] = [[0,1,0],[0,1,0]]
                    boundary.conditions['v-velocity'] = [[0,1,0],[0,1,0]]
                    boundary.orientation = side
        #calculating max nodes in y and max nodes in x for mesh files
        nx = 1
        ny = 0
        for rect in self.rectangles:
            if rect.height/self.dy()+1 > ny:
                ny = round(rect.height/self.dy())+1
            nx = nx + round(rect.width/self.dx())
        
        mesh.nx = nx
        mesh.ny = ny

        mesh.dy = self.dy
        mesh.dx = self.dx
        self.meshCreated.emit(mesh)

                


           
       
        
        
        
        
 
