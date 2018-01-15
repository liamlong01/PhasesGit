#!/usr/bin/env python
""""
An old prototype which was being designed to allow the user to draw a mesh directly.
The current mesh generator is a much simpler implementation for the user, and much easier to program as well.


This file never gets imported and nothing here is currently used by any other files ,
but it's kept here in case something becomes useful later on.
for example the NewSpanningTree class may prove to be useful in mesh generation techniques.

""""

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from math import *

class tool():

    """
A base class for defining painting tools.
Does nothing but states all methods that would be 
required by a tool.
    """
   
    def __init__(self):
        pass
        
    def mousePressEvent(self,event,occupiedObject):
        self.drawing = True


    def mouseMoveEvent(self,event,occupiedObject):
        pass

    def mouseReleaseEvent(self, event, occupiedObject):
        self.drawing = False
        

    def switch(self):
        return self

    def isBusy(self):
        return False


class toolBarTool(tool, QAction):

    """ adds funcionality to the tool class to make it compatibble with qt widgets and toolbars """

    def __init__(self, icon, text, parent):
        tool.__init__(self)
        QAction.__init__(self, icon, text, parent)


#class ClearTool(toolBarTool):

#    """a tool for drawing lines on a QtWidget"""

#    def __init__(self, parent):
#        icon = QIcon('Resources/clear.png')
#        text = 'Clear'
#        toolBarTool.__init__(self,icon,text,parent)


#    def paints(self, objects):
#        objects.clear()
        
class moveTool(toolBarTool):

    """
A base class for defining painting tools.
Does nothing but states all methods that would be 
required by a tool
    """
   
    def __init__(self,parent):
        icon = QIcon('Resources/pencil.png')
        text = 'Move'
        toolBarTool.__init__(self,icon,text,parent)
        
        self.moving = False
    
        
    def mousePressEvent(self,event,occupiedObject):
        self.moving = True
        self.object = occupiedObject

    def mouseMoveEvent(self,event,occupiedObject):
        pass

    def mouseReleaseEvent(self, event):
        self.moving = False
        

    def switch(self):
        return self

    def isBusy(self):
        return False      

class LineTool(toolBarTool):

    """a tool for drawing lines on a QtWidget"""

    def __init__(self, parent):
        icon = QIcon('Resources/pencil.png')
        text = 'Line'
        toolBarTool.__init__(self,icon,text,parent)
        
        self.drawingMode = None


    def mousePressEvent(self, point, occupiedObject):

        self.drawingMode = True
        
        
       
        self.originalNode = occupiedObject.createNewNode(point)
       
        self.linetodraw= Line()
        self.linetodraw.setP1(self.originalNode)
        

    def mouseMoveEvent(self,point, occupiedObject):
        if self.drawingMode is None:
            return
        currentNode = occupiedObject.centerNode(point)
        xDist = currentNode.x() - self.originalNode.x()
        yDist = currentNode.y() - self.originalNode.y()

        
        if fabs(xDist) > fabs(yDist):
            self.linetodraw.setP2(Node(self.originalNode.x() + xDist, self.originalNode.y()))
           
            
        else:
            self.linetodraw.setP2(Node(self.originalNode.x(), self.originalNode.y() + yDist))
 
        return self.linetodraw
        

    def mouseReleaseEvent(self, occupiedObject):
        if self.linetodraw.isNull():
            return []

        self.drawingMode = None
        self.linetodraw.setP2(occupiedObject.addNode(self.linetodraw.node2))
        if not self.linetodraw.isNull():
          
            return [self.linetodraw]
 

    def switch(self):
        self.drawingMode = None
        
        return self

    
    def isBusy(self):
        if self.drawingMode is not None:
            return True
        else:
            return False


class RectTool(toolBarTool):

    """a tool for drawing rectangles on a QtWidget"""

    def __init__(self, parent):
        icon = QIcon('Resources/rect.png')
        text='Rect'
        toolBarTool.__init__(self,icon,text,parent)

    def mousePressEvent(self, event):
        print('Woohoo rects!')

  
class Section(QRectF):

    def __init__(self):
        QRectF.__init__(self)
  

    def createNewNode(self, point):
        return Node(point.x(),point.y())

    def centerNode(self, point):
        return Node(point.x(),point.y())

    def draw(self, painter):

        painter.drawRect(self)

    def contextMenuActions(self, menu):
        return []

    def addNode(self, node):
        return node

class Node(Section):        

    def __init__(self,x,y):
        Section.__init__(self)
        self.lines =[]
        self.setCoords(x,y)
        self.Actions = NodeActions(self)
        self.discovered = False

        self.ancestors =[] #for SpanningTree structure
        
    def createNewNode(self, point):
        return self

    def centerNode(self, point):
        return self

    def addNode(self, node):
        if node.center() == self.center():
            return self
        else:
            return node
    
    def setCoords(self,x,y):
        self.setTopLeft(QPointF(x-3,y+3))
        self.setBottomRight(QPointF(x+3,y-3))
        
   
    def setX(self, x):
        self.setTopLeft(QPointF(x-3, self.y()+3))
        self.setBottomRight(QPointF(x+3, self.y()-3))
       

    def setY(self, y):
        self.setTopLeft(QPointF(self.x()-3, y+3))
        self.setBottomRight(QPointF(self.x()+3, y-3))

    def x(self):
        return self.center().x()
    
    def y(self):
        return self.center().y()

    def draw(self,painter):
        tempPen = painter.pen()
        tempBrush = painter.brush()
        if self.discovered:
            painter.setPen(Qt.red)
            painter.setBrush(Qt.red)
        painter.drawEllipse(self)
        painter.setPen(tempPen)
        painter.setBrush(tempBrush)

    def addLine(self, line):
        self.lines.append(line)

    def subLine(self, line):
        self.lines.remove(line)
        if not self.lines:
            del self

    def contextMenuActions(self,menu):
        return self.Actions.actions(menu)

    def lineTo(self, diffNode):
        for line in self.lines:
            if line.otherNode(self) == diffNode:
                return line
        return None

    #for generating SpanningTree structure with Nodes
    def addAncestor(self, node):
        self.ancestors.extend(node.ancestors)
        self.ancestors.append(node)

class NodeActions(QWidget):
 
    def __init__(self,node,parent=None):
        QWidget.__init__(self,parent)
        self.node = node

    def actions(self, menu):
        action = QAction("Print Node Info",menu)
        action.triggered.connect(self.printInfo)
        return [action]   
    
    def printInfo(self):
        print(self.node.lines)
        print(self.node.discovered)


class Line(Section):        


    def __init__(self):
        Section.__init__(self)
        self.node1 = None
        self.node2 = None
        self.Actions = LineActions(self)
        self.Actions.pointsChanged.connect(self.setPoints)

    def setPoints(self, point1, point2):
        self.setP1(point1)
        self.setP2(point2)

    def setP1(self,node):
        if self.node1 is not None:
            self.node1.subLine(self)
        node.addLine(self)
        self.node1  = node
        self.setTopLeft(node.center()+QPointF(-3, 3))
        

    def setP2(self,node):
        if self.node2 is not None:
            self.node2.subLine(self)
        node.addLine(self)
        self.node2 = node
        self.setBottomRight(node.center() +QPointF(3, -3))    

    def otherNode(self, node):
        if node == self.node1:
            return self.node2
        elif node == self.node2:
            return self.node1

    def createNewNode(self, point):
        newNode = self.centerNode(point)
        self.Actions.lineSplit.emit(self, newNode)
        return newNode
        
       
    def centerNode(self,point):
       if self.isVertical():
            point.setX(self.center().x())
       else:
            point.setY(self.center().y())
        
       return Node(point.x(),point.y())

    def addNode(self, node):
        if (self.isVertical() and node.center().x()==self.center().x()) \
            or node.center().y()==self.center().y():
            self.Actions.lineSplit.emit(self, node)
        return node

    def contextMenuActions(self,menu):
        return self.Actions.actions(menu)

    def isVertical(self):
        return (fabs(self.height()) > fabs(self.width()))

    def draw(self,painter):
        thinLine = QLineF(self.node1.center(), self.node2.center())

        painter.drawLine(thinLine)
        self.node1.draw(painter)
        self.node2.draw(painter)

    def Area(self):
        return [self.node1,self.node2, self]

    def isNull(self):
        if self.node1 ==None or self.node2==None:
            return True
    

        
        
class LineActions(QWidget):

    pointsChanged = pyqtSignal(Node,Node)
    lineSplit = pyqtSignal(Line, Node)

    def __init__(self,line,parent=None):
        QWidget.__init__(self,parent)
        self.line = line

    def actions(self, menu):
        action = QAction("Set Line length",menu)
        action.triggered.connect(self.lengthInput)
        return [action]

    def lengthInput(self):
        dialog = QInputDialog(self)
        length, ok= dialog.getDouble(self, 'Input Line Length', 'Length:')
         
        if self.line.isVertical():
            point1 = Node(self.line.center().x(), self.line.center().y() + length/2.0)
            point2 = Node(self.line.center().x(), self.line.center().y() - length/2.0)
        else:
            point1 = Node(self.line.center().x() - length/2.0, self.line.center().y())
            point2 = Node(self.line.center().x() + length/2.0, self.line.center().y())

        self.pointsChanged.emit(point1, point2)
        


class undoableAction():

    def __init__(self, *args, **kwargs):
        pass
    
    def undo(self):
        pass
    def redo(self):
        pass

class LineAdded(undoableAction):

    def __init__(self, parent):
        self.parent=parent

    def undo(self):
        self.line = self.parent.AllLines[-1]
        del self.parent.AllLines[-1]

    def redo(self):
        self.parent.AllLines.append(self.line)
        
class CanvasCleared(undoableAction):

    def __init__(self, parent, lines, nodes):
        self.parent=parent
        self.lines=lines.copy()
        self.nodes = nodes.copy()

    def undo(self):
        self.parent.AllLines = self.lines
        self.parent.AllNodes = self.nodes

    def redo(self):
        self.parent.AllLines = []
        self.parent.AllNodes = []

class Grid(object):

    def __init__(self, distance, topLeft, bottomRight):
        self.distance = distance #distance in pixels should eventually change to be a scale
        width = bottomRight.x() - topLeft.x()
        self.width = width + (distance- width % distance)
        height = bottomRight.y() - topLeft.y() 
        self.height = height + (distance -  height % distance)
        self.topLeft = topLeft #use a point
        self.bottomRight = Node(topLeft.x()+self.width, topLeft.y() + self.height)


    def draw(self, painter):
        
        for y in range(int(self.height/self.distance+1)):
            y = y*self.distance
            line = QLineF(QPointF(self.topLeft.x(), y), QPointF(self.bottomRight.x(), y))
            painter.drawLine(line)
        for x in range(int(self.width/self.distance+1)):
            x = x*self.distance
            line = QLineF(QPointF(x,self.topLeft.y()), QPointF(x, self.bottomRight.y()))
            painter.drawLine(line)


class MeshNode():
    pass
class MeshElement():
    pass




class Mesh(object):

    def __init__(self, elementSize, boundingPoints):
        self.boundaries = {
                    'left': [],
                    'right': [],
                    'top' : [],
                    'bottom' : [],
            }
        
        self.sortEdges(boundingEdges)
        boundingRect = None
       
        self.grid = Grid(elementsize, topReft, bottomRight)


        def sortEdges (self, edges):
            for edge in edges:
                pass

class meshGenerator(QWidget):

    def __init__(self):
        QWidget.__init__(self)

        layout = QVBoxLayout(self)

        self.toolbar=QToolBar(self)
        layout.addWidget(self.toolbar)
        self.canvas = meshPainter()

        layout.addWidget(self.canvas)

        #painting tools
        moveAction = QAction(QIcon('Resources/move.png'),'Move', self.toolbar)
        moveAction.triggered.connect(self.canvas.clear)
        self.toolbar.addAction(moveAction)

        self.rectTool = RectTool(self.toolbar)
        self.rectTool.triggered.connect(self.switchToRect)
        self.toolbar.addAction(self.rectTool)

        self.lineTool = LineTool(self.toolbar)
        self.lineTool.triggered.connect(self.switchToLine)
        self.toolbar.addAction(self.lineTool)

        clearAction = QAction(QIcon('Resources/clear.png'),'Clear', self.toolbar)
        clearAction.triggered.connect(self.canvas.clear)
        clearAction.setShortcut(Qt.CTRL + Qt.Key_P)
        self.toolbar.addAction(clearAction)

        self.toolbar.addSeparator()

        undoAction = QAction(QIcon('Resources/undo.png'),'Undo', self.toolbar)
        undoAction.triggered.connect(self.canvas.undo)
        undoAction.setShortcut(Qt.CTRL + Qt.Key_Z)
        self.toolbar.addAction(undoAction)

        redoAction = QAction(QIcon('Resources/redo.png'),'Redo', self.toolbar)
        redoAction.triggered.connect(self.canvas.redo)
        redoAction.setShortcut(Qt.CTRL + Qt.Key_A)
        self.toolbar.addAction(redoAction)

        self.toolbar.addSeparator()

        meshDrawAction = QAction(QIcon('Resources/DrawMesh.png'),'Make Mesh', self.toolbar)
        meshDrawAction.triggered.connect(self.canvas.clear)
        self.toolbar.addAction(meshDrawAction)

        meshAction = QAction(QIcon('Resources/mesh.png'),'Make Mesh', self.toolbar)
        meshAction.triggered.connect(self.createMesh)
        self.toolbar.addAction(meshAction)


        position = QLabel(self.toolbar)
        self.canvas.mouseMoved.connect(position.setText)
        self.toolbar.addWidget(position)
       
    def switchToRect(self):
        self.canvas.switchTool(self.rectTool)

    def switchToLine(self):
        self.canvas.switchTool(self.lineTool)

    def createMesh(self):
        dialog = QInputDialog(self)
        numNodes, ok= dialog.getInt(self, 'Generate Mesh', 'Number of Elements:')
        cycle = self.canvas.isMeshValid()
        if cycle:
            points =[]
            for node in cycle:
                points.append(node.center())
            boundary = QPolygonF(points)
            area = boundary.boundingRect()

            #   First determine if edges are left, bottom, right, or top
            #   Create a grid with set scale within bounding rectangle
            #   snap lines in grid to edges:
            #       -eg. for a left edge take the line just to the left and snap it over
            #       -if a line is exactly on an edge, use that one
            #   create nodes where grid lines intersect within polygon
            #   create elements from nodes
            
    def setMesh(self, mesh):
        canvas.AllLines = []
        canvas.AllNodes = []
        for element in mesh.elements:
            for node in element:
                canvas.Allnodes.append()

class meshPainter(QWidget):

    mouseMoved = pyqtSignal(str)

    def __init__(self):
    
        QWidget.__init__(self)
        self.currentObject = Section()
        self.setMouseTracking(True)
        
        self.grid  = Grid(25, Node(0,0), Node(1650, 875))

        self.tool= tool()
 
        self.currentShape = None
        self.boundary = None

        self.AllLines = []
        self.AllNodes =[]

        self.previousActions =[]
        self.undoneActions =[]

        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.showContextMenu)


    def switchTool(self, action):
        self.tool = action.switch()
        

    def clear(self):
        
        self.addAction(CanvasCleared(self, self.AllLines, self.AllNodes))
        self.AllLines = []
        self.AllNodes =[]
        self.currentObject = Section()
        self.update()

    def undo(self):
        if self.previousActions:
            self.previousActions[-1].undo()
            self.undoneActions.append(self.previousActions[-1])
            del self.previousActions[-1]
            self.update()

    def redo(self):
        if self.undoneActions:
            self.undoneActions[-1].redo()
            self.previousActions.append(self.undoneActions[-1])
            del self.undoneActions[-1]
            self.update()

    def showContextMenu(self, point): 
     
        menu = QMenu('Context Menu', self)
        dAction = QAction("Default Action", menu)
        
        menu.addAction(dAction)
        menu.addSeparator()

        for action in self.currentObject.contextMenuActions(menu):
            menu.addAction(action)

        menu.exec(self.mapToGlobal(point))

    def paintEvent(self, event):
       
        #if no shapes have been added do nothing
        

        painter = QPainter()
        painter.begin(self)
 
        painter.setPen(Qt.lightGray)
        self.grid.draw(painter)

        

        painter.setPen(Qt.black)
        painter.setBrush(Qt.black)
        if self.currentShape is not None:
            self.currentShape.draw(painter)
         
        if self.boundary is not None:
            painter.drawPolygon(self.boundary)
        
        
        if self.AllLines:
            for line in self.AllLines:
                line.draw(painter)

        if self.currentObject is not None:
            painter.setPen(Qt.green)
            painter.setBrush(Qt.green)
            self.currentObject.draw(painter)

        

        painter.end()
    
    def mousePressEvent(self, event):
        
        if event.button() == Qt.LeftButton:
            self.tool.mousePressEvent(event.pos(), self.currentObject)
        elif event.button() ==Qt.RightButton:
            self.customContextMenuRequested.emit(event.pos())       
    
    def mouseMoveEvent(self, event):
        self.isOccupied(event.pos())
        

        self.currentShape = self.tool.mouseMoveEvent(event.pos(),self.currentObject)
        x = round(event.x())
        y = round(event.y())

        self.mouseMoved.emit("x: %s" % x + " y: %s" % y)
        self.update()
      
    def mouseReleaseEvent(self, event):
        newLines = self.tool.mouseReleaseEvent(self.currentObject)
        if newLines:
            for line in newLines:
                self.addLine(line)
            
        self.update()
    
    def addAction(self, action):
        self.previousActions.append(action)
        self.undoneActions = []

    def isOccupied(self, point):
        for line in self.AllLines:
            for section in line.Area():
                if section.contains(point):
                    self.currentObject=section
                    return True
                else:
                    self.currentObject = Section()

    def addLine(self, line):
        line.Actions.lineSplit.connect(self.splitLine)
        self.AllLines.append(line)
        self.addNode(line.node1)
        self.addNode(line.node2)
        self.addAction(LineAdded(self))

    def addNode(self, node):
        if node not in self.AllNodes:
            self.AllNodes.append(node)  

    def splitLine(self, line, node):
        newLine = Line()
        newLine.setP1(line.node2)
        newLine.setP2(node)
        self.addLine(newLine)
        #self.addAction(LineSplit())
        line.setP2(node)

    def isMeshValid(self):

        for node in self.AllNodes:
            node.discovered=False
                
        tree = NewSpanningTree(self.AllNodes[0])


        allCycles=[]
        for line in self.AllLines:
            if line not in tree.edges and line.node1.discovered:
                #cycle!!!
                path = tree.findPath(line.node1, line.node2)

                path.append(line.node1)     #makes the path complete
               

                edgePath = tree.pathToEdges(path)
                allCycles.append(edgePath)
                #need to find a path in spanning tree from startnode to endnode or vice versa
        
        
        for node in self.AllNodes:
            node.discovered = False
        
        for node in path:
            node.discovered = True

        tree.__del__()

        if len(allCycles) == 1:
            print(allCycles[0])
            return allCycles[0]
        else:
            return False          

    def sizeHint(self):
        return QSize(500, 500)

class NewSpanningTree(object):
    
    def __init__(self, startingNode):
        self.nodeAncestors = {}
        self.edges = []
        self.DFS(startingNode, None)

    def addNode(self, node, prevNode):
        if prevNode  in self.nodeAncestors.keys():
            self.nodeAncestors[node] = [prevNode] + self.nodeAncestors[prevNode]
        else:
            self.nodeAncestors[node]= []
       
            
  
        
        
     # this function creates a spanning tree using a depth-first search algorithm
        # results are stored in nodesVisited and treeEdges.
        #      
        #
        # if any edges connected to a visited node are not part of the spanning tree
        # then a cycle exists within the nodes of the tree
    def DFS(self, node, prevNode):
        node.discovered = True
        nodeFound = False
        for line in node.lines:
                
            nextNode = line.otherNode(node)

            if not nextNode.discovered:
                self.addNode(node, prevNode)
                self.edges.append(nextNode.lineTo(node))
                nodeFound=True
                self.DFS(nextNode, node)
        if not nodeFound:
            self.addNode(node, prevNode)


    #This functions returns the path between two nodes in the SpanningTree structure

    def findPath(self, node1, node2):

        if node2 in self.nodeAncestors[node1]:
            # node1 is an ancestor of node2
            path = self.directPath(node1, node2)

        elif node1 in self.nodeAncestors[node2]:
            # node2 is an ancestor of node1
            path = self.directPath(node2, node1)  
            path.reverse()

        else:
            #the nodes share a common ancestor and the path must be through that ancestor
            for ancestor in self.nodeAncestors[node1]:
                if ancestor in self.nodeAncestors[node2]:
                    #ancestor is in both node1.ancestors and node2.ancestors
                    part1 = self.directPath(node1, ancestor)        # the path from node 1 to common ancestor  
                    part2 = self.directPath(node2, ancestor)        # plus the path from common ancestor to node 2. 
                    path =  part1 + part2[:-1].reverse()            # the [:-1] ensures the common ancestor isn't added twice 
                                                                                                        
                                                                                                        
                    break
      

        return path

    # this function returns the path between two nodes if one is a direct ancestor of the other 
    #   descendant  - a Node further down in the tree
    #   ancestor    - a Node higher up in the tree
    #
    def directPath(self, descendant, ancestor):
        directPath = [descendant]
        while directPath[-1] != ancestor:
            directPath.append(self.nodeAncestors[directPath[-1]][-1])
        return directPath

    def pathToEdges(self, pathNodes):
        pathEdges = []
        for i in range(len(pathNodes)-1):
            line = pathNodes[i].lineTo(pathNodes[i+1])
            if line is not None:
                pathEdges.append(line)
            else: #Error cycle is not valid
                return None
                
        return pathEdges

    

class SpanningTree(object):
    
    def __init__(self, startingNode):
        self.nodes = {}
        self.edges = []
        self.DFS(startingNode, None)

    def addNode(self, node, prevNode):
        self.nodes.append(node)
        if prevNode is not None:
            node.addAncestor(prevNode)
        
        
     # this function creates a spanning tree using a depth-first search algorithm
        # results are stored in nodesVisited and treeEdges.
        #      
        #
        # if any edges connected to a visited node are not part of the spanning tree
        # then a cycle exists within the nodes of the tree
    def DFS(self, node, prevNode):
        node.discovered = True
        nodeFound = False
        for line in node.lines:
                
            nextNode = line.otherNode(node)

            if not nextNode.discovered:
                self.addNode(node, prevNode)
                self.edges.append(nextNode.lineTo(node))
                nodeFound=True
                self.DFS(nextNode, node)
        if not nodeFound:
            self.addNode(node, prevNode)


    #This functions returns the path between two nodes in the SpanningTree structure

    def findPath(self, node1, node2):

        if node2 in node1.ancestors:
            # node1 is an ancestor of node2
            path = self.directPath(node1, node2)

        elif node1 in node2.ancestors:
            # node2 is an ancestor of node1
            path = self.directPath(node2, node1)  
            path.reverse()

        else:
            #the nodes share a common ancestor and the path must be through that ancestor
            for ancestor in node1.ancestors:
                if ancestor in node2.ancestors:
                    #ancestor is in both node1.ancestors and node2.ancestors
                    part1 = self.directPath(node1, ancestor)        # the path from node 1 to common ancestor  
                    part2 = self.directPath(node2, ancestor)        # plus the path from common ancestor to node 2. 
                    path =  part1 + part2[:-1].reverse()            # the [:-1] ensures the common ancestor isn't added twice 
                                                                                                        
                                                                                                        
                    break
      

        return path

    # this function returns the path between two nodes if one is a direct ancestor of the other 
    #   descendant  - a Node further down in the tree
    #   ancestor    - a Node higher up in the tree
    #
    def directPath(self, descendant, ancestor):
        directPath = [descendant]
        while directPath[-1] != ancestor:
            directPath.append(directPath[-1].ancestors[-1])
        return directPath

    def pathToEdges(self, pathNodes):
        pathEdges = []
        for i in range(len(pathNodes)-1):
            line = pathNodes[i].lineTo(pathNodes[i+1])
            if line is not None:
                pathEdges.append(line)
            else: #Error cycle is not valid
                return None
                
        return pathEdges

    def __del__(self):
        for node in self.nodes:
            node.ancestors = []
     
     
        
      
    
           