# -*- coding: utf-8 -*-

from PyQt5.QtCore import * #gui
from Phases.Mesh import Mesh  #keep track of elements, nodes and boundaries
from PyQt5.QtGui import * #gui
from PyQt5.QtWidgets import * #gui

import numpy as np # python libraries for arrays and matrices

class Rectangle(object): #describes a rectangle
  

    def __init__(self, height, width, nodes,edges):
        self.nodes = nodes
        self.height = height
        self.width = width
        self.edges = edges
        self.topEdges = [edges[1]]
        self.leftEdges = [edges[0]]
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



class RectGenerator(QWidget): #defines gui layout rectangle

    meshCreated = pyqtSignal(Mesh.Mesh)
    
    def __init__(self, nodeEditor, edgeEditor, display, parent = None):
        QWidget.__init__(self, parent)

        self.nodeEditor = nodeEditor
        self.edgeEditor = edgeEditor
        self.display = display

        self.rectangles = []

        layout = QVBoxLayout(self)

        # gui layout for element size
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
        #defines rectangle height 
        self.newhEdit = QLineEdit()
        newLayout.addWidget(QLabel('Height :'),0,0)
        newLayout.addWidget(self.newhEdit,0,1)
        heightNew = lambda: float(self.newhEdit.text())
        #defines rectangle width
        self.newwEdit = QLineEdit()
        newLayout.addWidget(QLabel('Width :'),1,0)
        newLayout.addWidget(self.newwEdit,1,1)
        widthNew = lambda: float(self.newwEdit.text())
        #alternates betwween L shape and T shape
        self.alignment = QComboBox()
        self.alignment.addItems(['Middle','Bottom', 'Top'])
        newLayout.addWidget(self.alignment,2,0,1,2)

        #sets up button for adding rectangle
        addnew = QPushButton('Add')
        newLayout.addWidget(addnew,3,0)
        addnew.clicked.connect(lambda: self.addNewRectangle(heightNew(),widthNew()))
        newLayout.setRowStretch(4,1)
        
        # sets up button to generates nodes
        generate = QGroupBox('Generate')
        genLayout = QGridLayout(generate)

        ranMesh = QCheckBox('Randomize Mesh')
        genLayout.addWidget(ranMesh,0,0)
        self.randomMesh = lambda: ranMesh.checkState()

        genNodes = QPushButton('Generate')
        genLayout.addWidget(genNodes,1,0)
        genNodes.clicked.connect(self.generateNodes)


     

        layout.addWidget(res) # defines element size, dx and dy
        layout.addWidget(newCreate) # creates rectangle
        layout.addWidget(generate) #generates mesh

        layout.addStretch(3)

    def addRectangle(self, refX,refY, height,width):
        """
        Adds a rectangle with the specified height and width. refX and refY define the location of the bottom left corner.
        """
        #generating the 4 corners of the rectangle in order topright, topleft, bottom left, bottom right
        node1 = self.nodeEditor.addNode(refX+width,refY+height)
        node2 = self.nodeEditor.addNode(refX,refY+height)
        node3 = self.nodeEditor.addNode(refX,refY)
        node4 = self.nodeEditor.addNode(refX+width,refY)

        #generating the four edges of the rectangle in orde of top,right,bottom,left
        edge2 = self.edgeEditor.addStraightLine(node1,node2)
        edge1 = self.edgeEditor.addStraightLine(node2,node3)
        edge3 = self.edgeEditor.addStraightLine(node3,node4)
        edge4 = self.edgeEditor.addStraightLine(node4,node1)

        nodes = [node1, node2, node3, node4]
        edges = [edge1, edge2, edge3, edge4]
        for edge in edges:
            edge.setDetail(0)
        
        #adding the rectangle
        self.rectangles.append(Rectangle(height, width, nodes,edges))

        #if its not the first rectangle, the rectangles will share edges with previous rectangles
        #this logic determines which edges need to be deleted to merge the two together
        if len(self.rectangles) > 1:

            #deleting right edge of previous rectange and left edge of new recrangle
            self.deleteRectEdge(self.rectangles[-2], self.rectangles[-2].rightEdges[0])
            self.deleteRectEdge(self.rectangles[-1], self.rectangles[-1].leftEdges[0])


            #edge connecting top right of previous rect to top left of new rect
            edge = self.edgeEditor.addStraightLine(self.rectangles[-2].nodes[0],self.rectangles[-1].nodes[1])
            edge.setDetail(0)

            #determines which rectangle the new edge belongs to
            if self.rectangles[-2].nodes[0].y > self.rectangles[-1].nodes[1].y:
                #corner of previous rect is above new rect 
                #this edge belongs to previous rect
                self.rectangles[-2].edges.append(edge)
                self.rectangles[-2].rightEdges.append(edge)
            else:
                #corner of new rect is above prev rect 
                #this edge belongs to new rect
                self.rectangles[-1].edges.append(edge)
                self.rectangles[-1].leftEdges.append(edge)
            
            #edge connecting bottom right of previous rect to bottom left of new rect
            edge = self.edgeEditor.addStraightLine(self.rectangles[-2].nodes[3],self.rectangles[-1].nodes[2])
            edge.setDetail(0)
            if self.rectangles[-2].nodes[3].y < self.rectangles[-1].nodes[2].y:
                #corner of previous rect is below new rect 
                #this edge belongs to previous rect
                self.rectangles[-2].edges.append(edge)
                self.rectangles[-2].rightEdges.append(edge)
            else:
                #corner of previous rect is above new rect 
                #this edge belongs to new rect
                self.rectangles[-1].edges.append(edge)
                self.rectangles[-1].leftEdges.append(edge)
          

    def addNewRectangle(self,height, width):
        """
        Adds a new Rectangle to the mesh boundary definition with the specified height and width. The alignment that has been selected is also taken into account.
        """

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
            #calculate location of bottom left corner for new rectangke based on alignment
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
        """
        Deletes the specified edge from rect.
        """
        rect.removeEdge(edge)
        self.edgeEditor.deleteEdge(self.edgeEditor.edges.index(edge)+1)

    def generateNodes(self):
        """
        Main mesh generation algorithm. Splits the region enclosed by the boundary into equally sized rectangles.
       Generates and organizes all nodes, elements and boundaries. Randomizes mesh if the Randomize Mesh checkbox is selected. 
       Emits the meshCreated signal.
        """
        #first lets shift all rectangles up that all coordinates are positive
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
            #making nodes in evenly spaced grid
            nodex = np.arange(round(rect.nodes[2].x,8),round(rect.nodes[0].x+self.dx(),8),self.dx())
            nodey = np.arange(round(rect.nodes[2].y,8),round(rect.nodes[0].y+self.dy(),8),self.dy())
            for i in range(len(nodex)):

                for j in range(len(nodey)):
                    node = Mesh.Node(nodex[i],nodey[j])

                    #if node is a corner copy the corner node to avoid duplicates
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
                    nodes[edgeNodes[i]].isBoundary = True
                    nodes[edgeNodes[i+1]].isBoundary = True
                    
                    
                    


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

        #make elements more randomized, ie not perfect rectangles
        #irregular quadrilaterals
        if self.randomMesh():
            for node in mesh.nodes:
                node.nudge(self.dx(), self.dy())
        else:
            mesh.dx = self.dx()
            mesh.dy = self.dy()
		#reSkew = True	
		#while reSkew:
		#	reSkew = False
		#	for elem in mesh.elements:
		#		if elem.badSkew():
		#			reSkew True
		#			for node in elem.nodes:
		#				node.nudge()
		#				if not elem.badSkew():
		#					break
		#					
        self.meshCreated.emit(mesh)

                


           
       
        
        
        
        
 