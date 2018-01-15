# -*- coding: utf-8 -*-
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *

import Mesh.Mesh as Mesh

import os


class Scale(object):
    
    
    def __init__(self,initval):
        self.val  = initval
    
    def __get__ (self,insatnce, owner):
        return self.val
    
    def __set__(self, instance, value):
        self.val = value

class Tool(object):

    """
A base class for defining painting tools.
Does nothing but states all methods that would be 
required by a tool. This makes it easy to interface with each tool through duck-typing.
Each individual tool should inherit from this class to be used within a widget object.
    """
    
    scaleX = Scale(1)   # defined this way the same scaleX and scaleY are used
    scaleY = Scale(1)   # by all instances of the Tool class even if 
                        # they are child classes
    
    def __init__(self,*args, **kwargs):
       pass
        
    def mousePressEvent(self,*args,**kwargs):
        """
        This method states what the tool will do when the mouse is pressed.
        
        """
        self.drawing = True


    def mouseMoveEvent(self,*args,**kwargs):
        """
        This method states what the tool will do when the mouse is moving.
        
        """
        pass

    def mouseReleaseEvent(self, *args,**kwargs):
        """
        This method states what the tool will do when the mouse released.
        """
        self.drawing = False
        

    def switchTo(self):
        """
        this method performs actions that need to be done when thisthe tool is selected.
        Should return a reference to self.
        """
        return self
        
    def switchFrom(self):
        """
        this method performs any actions that should be done  when the tool is deselected
        """
        pass

    def isBusy(self):
        """
        This method returns true if the tool is curretnly performing an action.
        """
        return False
        
    def drawStuff(self,  *args,**kwargs):
       """
       This method is used to draw things on a display by the tool
       
       """
       pass
   
   
   
    def setScaleX(self,new = None):
      
        if new:
            self.scaleX = new
        return self.scaleX
    
    def setScaleY(self,new = None):
        
        if new:
            self.scaleY = new
        
        return self.scaleY
        


class ToolbarTool(Tool, QAction):

    """ 
    Adds funcionality to the Tool class to make it compatable with qt widgets and toolbars.
    This includes inheriting from th Qt QAction class, and an additional method for drawing objects
    with the QPainter class
    
    """

    def __init__(self, icon, text, parent, *args, **kwargs):

        icon = os.path.join(os.path.dirname(__file__),icon)

        Tool.__init__(self,*args, **kwargs)
        QAction.__init__(self, QIcon(icon), text, parent)
        self.setCheckable(True)
        
    def drawStuff(self, painter, *args,**kwargs):
       """
       This method will use the painter which should be an instance of the QPainter class to draw things using the objectToDraw object.
       
       """
       pass
   
    def contextMenuActions(self,  *args,**kwargs):
       """
       This controls what entries are in the context menu that appears when the right mouse button is clicked.
       menu must be an instance of the QMenu object.
       Should return a list of lists of QAction, where each list QActions are to be seperated by menu seperators
      
       
       """
       pass
   
   
    def switchTo(self):
        """
        this method performs actions that need to be done when thisthe tool is selected.
        Should return a reference to self.
        """
        self.setChecked(True)
        return self
        
    def switchFrom(self):
        """
        this method performs any actions that should be done  when the tool is deselected
        """
        self.setChecked(False)
        
        
class DisplayTool(ToolbarTool):
    
    def __init__(self,icon='Icons/Blank.png)', text='', parent=None, *args, **kwargs):
        super().__init__(icon, text, parent, *args, **kwargs)
        
class SelectTool(ToolbarTool):
    
    def __init__(self, icon='Icons/Blank.png)', text='', parent=None, *args, **kwargs):
        super().__init__(icon, text, parent, *args, **kwargs)
        self.selections =[]

    def clearSelections(self):
        self.selections = []
   
class BCEditDialog(QDialog):
        """
        A dialog that allows the user to define new boundary conditions.
        """
        
        def __init__(self,parent, conditions):
            QDialog.__init__(self,parent=parent)
            
            self.layout = QFormLayout(self)
            
            self.oldConditions = conditions
            self.tempConditions = conditions.copy()
            
            
            self.layout.addRow('Boundary Conditions :', QLabel('aT\' + bT = c'))
            
            self.addConditionEdit('temperature')
            self.addConditionEdit('concentration')
            self.addConditionEdit('u-velocity')
            self.addConditionEdit('v-velocity')
            
            #TODO: these buttons should be changed to a button box for cross platform readability
            OkButton = QPushButton('Ok', parent = self)
            OkButton.clicked.connect(self.accept)
            CancelButton = QPushButton('Cancel', parent = self)
            CancelButton.clicked.connect(self.reject)
            Buttons = QHBoxLayout()
            Buttons.addWidget(OkButton)
            Buttons.addWidget(CancelButton)

            self.layout.addRow('', Buttons)
            
            
        def addConditionEdit(self,label):
            
            def coefTextBox(label, index):
            
                textbox = QLineEdit(str(self.tempConditions[label][0][index]),parent = self)
                # set the only possible values to be doubles
                # this also includes the empty string 
                # as a possible value in the text box
                textbox.setValidator(QDoubleValidator(textbox))
                
                changeCond =lambda z: self.addCondition(label,index,z)
                textbox.textChanged.connect(changeCond)
                return textbox
            
            lineLayout = QHBoxLayout()
            
            lineLayout.addWidget(QLabel('a:'))
            lineLayout.addWidget(coefTextBox(label, 0))
            
            lineLayout.addWidget(QLabel('b:'))
            lineLayout.addWidget(coefTextBox(label, 1))
            
            lineLayout.addWidget(QLabel('c:'))
            lineLayout.addWidget(coefTextBox(label, 2))
            
            
            lineWidget = QWidget()
            lineWidget.setLayout(lineLayout)
            
            self.layout.addRow('%s :' %label, lineWidget)
   
            
        def addCondition(self,label,index, toWhat):
            try:
                float(toWhat)
            except ValueError:
                # toWhat is empty string do nothing
                return
                     
            self.tempConditions[label][0][index] = float(toWhat)
            self.tempConditions[label][1][index] = float(toWhat)


        def execute(self):
            if self.exec():
                return self.tempConditions
            else:
                return None
            

            
            
class BCEditTool(DisplayTool):
    """
    The tool used to change boundary conditions.
    Only allows select tools to select boundaries
    """
    
    def __init__(self,parent):
        icon = 'Icons/boundaryConditions.png'
        text = 'Edit Boundary Conditions'
        
        self.boundaries = None
        self.params = None
        
          
        super().__init__(icon,text,parent)
        
        
    def selectableObjects(self):
        return self.boundaries
            
    def mouseMoveEvent(self, x, y):
        
        #this loop will check to see if the current mouse position of x,y is occupied by or near to a boundary object in the mesh
        
        return self.findCurrentBoundary(x,y)
            
        
   
    
    def findCurrentBoundary(self, x, y):
        """
        claculates a line function for each boundary.
        Then returns a boundary if the mouse is currently
        near that line.
        """
        
        if self.boundaries is not None:
            for boundary in self.boundaries:
                xlims = [node.x for node in boundary.nodes]
                ylims = [node.y for node in boundary.nodes]
                
                # Creating a line equation in the form y=mx+b that characterises the boundary
                # This cone be done because we know the nodes at each end which containn their coordinates
                
               
                fuzzFact = 8 # this is used to expand the amount of space taken up by the boundary
                             # so that to select it the user does not need to be pixel perfect 
                             # ie make the line a bit fuzzy
                
                
                if xlims[0]==xlims[1]: #vertical line, m is infinty so this case needs to be handled seperately
                    if y < max(ylims)+fuzzFact/self.scaleY and y >min(ylims)-fuzzFact/self.scaleY:
                         if x < xlims[0]+fuzzFact/self.scaleX and x > xlims[0]-fuzzFact/self.scaleX:
                             return boundary  
      
                else:   
                    m = (ylims[1]-ylims[0])/(xlims[1]-xlims[0])
                    b = ylims[0] - m*xlims[0]
                    linefunc = lambda z: m*z+b
               
                             
                    if x<=max(xlims)+fuzzFact/self.scaleX and x>=min(xlims)-fuzzFact/self.scaleX: 
                        if y < linefunc(x)+fuzzFact/self.scaleY and y > linefunc(x)-fuzzFact/self.scaleY:
                            return boundary
                
                    
        
    def drawStuff(self,painter,boundary):
        """
        Takes a QPainter and uses it to draw the given boundary
        """
        
        points = [QPointF(node.x*self.scaleX,node.y*self.scaleY) for node in boundary.nodes]
        pen = QPen(Qt.yellow)
        pen.setWidth(8)          
        painter.setPen(pen)  
        painter.drawLine(QLineF(points[0],points[1]))
    
    def contextMenuActions(self,menu, boundaries):
        """
        Returns a single QAction that opens a bcEdit dialog when triggered.
        """
        
        if boundaries is not None:
            action = QAction("Change Initial Conditions",menu)
            action.triggered.connect(lambda: self.editBC(boundaries))
            return [[action]]
        else:
            return [[]]
            
    def editBC(self, boundaries):
        
       dialog = BCEditDialog(self.parent(), boundaries[0].conditions)
        
       newConditions = dialog.execute()
       if newConditions is None:
           return
       
       for boundary in boundaries:
           boundary.conditions = newConditions.copy() # copybecause each set of bcs needs to remain independent
        
class ICEditDialog (QDialog):
    """
    A dialog that allows the user to specify initial conditions on mesh elements
    """

    def __init__(self, parent,conditions):
        QDialog.__init__(self,parent=parent)
            
        self.layout = QFormLayout(self)
        
        self.oldConditions = conditions
        self.tempConditions = conditions.copy()

        self.addConditionEdit('temperature')
        self.addConditionEdit('concentration')
        self.addConditionEdit('u-velocity')
        self.addConditionEdit('v-velocity')
        
        #TODO: these buttons should be changed to a button box for cross platform readability
        OkButton = QPushButton('Ok', parent = self)
        OkButton.clicked.connect(self.accept)
        CancelButton = QPushButton('Cancel', parent = self)
        CancelButton.clicked.connect(self.reject)
        Buttons = QHBoxLayout()
        Buttons.addWidget(OkButton)
        Buttons.addWidget(CancelButton)

        self.layout.addRow('', Buttons)  
            
    def addConditionEdit(self,label):
         
        textbox = QLineEdit(str(self.tempConditions[label]),parent = self)
        textbox.setValidator(QDoubleValidator(textbox))
        changeCond =lambda z: self.addCondition(label,z)
        textbox.textChanged.connect(changeCond)
        
        self.layout.addRow('%s :' %label, textbox)
   
        
    def addCondition(self,label,toWhat):
        try:
            float(toWhat)
        except ValueError:
            # toWhat is empty string do nothing
            return
                 
        self.tempConditions[label]= float(toWhat)
        


    def execute(self):
        if self.exec():
            return self.tempConditions
        else:
            return None
            
                  
     
        
class ICEditTool(DisplayTool):
    """
    The tool used to edit initial conditions.
    Only allows select tools to select mesh elements.
    """
    
    def __init__(self,parent):
        icon = 'Icons/initConditions.png'
        text = 'Edit Initial Conditions'
        super().__init__(icon,text,parent)
        
        self.elements = None
        self.params = None
        
    def selectableObjects(self):
        return self.elements
        
    def mouseMoveEvent(self,x,y):
        """
        Returns an element if the mouse is currently inside that element.
        """
        if self.elements is not None:
            for element in self.elements:
                points = [QPointF(node.x,node.y) for node in element.nodes]
                polygon = QPolygonF(points)
                
                if polygon.containsPoint(QPointF(x,y),Qt.OddEvenFill):
                    return element
    
    def drawStuff(self, painter, element):
        """
        Takes a QPainter and uses it to draw the given element.
        """
        points = [QPointF(node.x*self.scaleX,node.y*self.scaleY) for node in element.nodes]
                  
        color = QColor(255,255,0,100) #transparent yellow         
        painter.setPen(color)
        painter.setBrush(color)          
                  
        painter.drawPolygon(QPolygonF(points))
    
    def editIC(self, elements):
        dialog = ICEditDialog(self.parent(), elements[0].initialConditions)
        
        newConditions = dialog.execute()
        if newConditions is None:
           return
           
        for element in elements:
            element.initialConditions = newConditions.copy() #copy becuase each set of ics needs to remain independent
        
 
    
    def subDivide(self, elements):
        """
        divides an element into 4 smaller pieces.
        currently only works with rectangular elements
        """
        for element in elements:
            maxX = element.nodes[0].x
            minX = element.nodes[0].x
            maxY = element.nodes[0].y
            minY = element.nodes[0].y
            for node in element.nodes:
                if node.x > maxX:
                    maxX=node.x
                elif node.x<minX:
                    minX=node.x
                if node.y >maxY:
                    maxY=node.y
                elif node.y<minY:
                    minY=node.y
                    
            x = minX + (maxX-minX)/2.0
            y = minY + (maxY-minY)/2.0

            #initializing new nodes
            centreNode = Mesh.Node(x,y) 
            rightNode = Mesh.Node(maxX,y) 
            leftNode = Mesh.Node(minX,y) 
            bottomNode = Mesh.Node(x,minY) 
            topNode = Mesh.Node(x,maxY) 
            
            #making sure no duplicate nodes are being created
            for node in self.parent().parent().mesh.nodes:
                if node.x == x and node.y == y:
                    centreNode = node
                if node.x == maxX and node.y == y:
                    rightNode = node
                if node.x == minX and node.y == y:
                    leftNode = node
                if node.x == x and node.y == minY:
                    bottomNode = node
                if node.x == x and node.y == maxY:
                    topNode = node
                    
            
            
            newElements = [
                           Mesh.Element([element.nodes[0],topNode,centreNode,rightNode]),
                           Mesh.Element([topNode,element.nodes[1],leftNode,centreNode]),
                           Mesh.Element([centreNode,leftNode,element.nodes[2],bottomNode]),
                           Mesh.Element([rightNode,centreNode,bottomNode, element.nodes[3]])
                          ]
            
            self.parent().parent().mesh.subElement(element)
            for newElement in newElements:
                self.parent().parent().mesh.addElement(newElement)
                newElement.initialConditions = element.initialConditions.copy()
                
            if element.isBoundary:
                for boundary in element.boundaries:
                    if element.nodes[0] in boundary.nodes and element.nodes[1] in boundary.nodes: #top side boundary
                        newBoundary1 = Mesh.Boundary([element.nodes[0],topNode])
                        newBoundary1.element = newElements[0]
                        newElements[0].boundaries.append(newBoundary1)
                        newElements[0].isBoundary = True
                        
                        newBoundary2 = Mesh.Boundary([element.nodes[1],topNode])
                        newBoundary2.element = newElements[1]
                        newElements[1].boundaries.append(newBoundary2)
                        newElements[1].isBoundary = True
                        
                        
                    elif element.nodes[1] in boundary.nodes and element.nodes[2] in boundary.nodes: #left side boundary
                        newBoundary1 = Mesh.Boundary([element.nodes[1],leftNode])
                        newBoundary1.element = newElements[1]
                        newElements[1].boundaries.append(newBoundary1)
                        newElements[1].isBoundary = True
                        
                        newBoundary2 = Mesh.Boundary([element.nodes[2],leftNode])
                        newBoundary2.element = newElements[2]
                        newElements[2].boundaries.append(newBoundary2)
                        newElements[2].isBoundary = True
                        
                    elif element.nodes[2] in boundary.nodes and element.nodes[3] in boundary.nodes: #bottom side boundary
                        newBoundary1 = Mesh.Boundary([element.nodes[2],bottomNode])
                        newBoundary1.element = newElements[2]
                        newElements[2].boundaries.append(newBoundary1)
                        newElements[2].isBoundary = True
                        
                        newBoundary2 = Mesh.Boundary([element.nodes[3],bottomNode])
                        newBoundary2.element = newElements[3]
                        newElements[3].boundaries.append(newBoundary2)
                        newElements[3].isBoundary = True
                        
                    elif element.nodes[3] in boundary.nodes and element.nodes[0] in boundary.nodes: #right side boundary
                        newBoundary1 = Mesh.Boundary([element.nodes[3],rightNode])
                        newBoundary1.element = newElements[3]
                        newElements[3].boundaries.append(newBoundary1)
                        newElements[3].isBoundary = True
                        
                        newBoundary2 = Mesh.Boundary([element.nodes[0],rightNode])
                        newBoundary2.element = newElements[0]
                        newElements[0].boundaries.append(newBoundary2)
                        newElements[0].isBoundary = True
                        
                    
                    
                    newBoundary1.conditions = boundary.conditions.copy()
                    newBoundary2.conditions = boundary.conditions.copy()
                    
                    self.parent().parent().mesh.subBoundary(boundary)
                    self.parent().parent().mesh.addBoundary(newBoundary1)
                    self.parent().parent().mesh.addBoundary(newBoundary2)
                          
    
    
    def contextMenuActions(self, menu, selections):
        """
        Returns two actions.
        One opens an ICEdit dialog for changing initial conditions on selected elements.
        The other subdivides all selected elements into four smaller pieces.
        """
        if selections is not None:
            editAction = QAction("Edit ICs",menu)
            editAction.triggered.connect(lambda: self.editIC(selections))
            
            subdivAction = QAction("Subdivide Element", menu)
            subdivAction.triggered.connect(lambda: self.subDivide(selections))
            return [[editAction,subdivAction]]
        else:
            return [[]]


class SingleSelect(SelectTool):
    """
    Selects objects one at a time.
    Holding Ctrl allows the user to select multiple objects at once
    """
    
    def __init__(self, parent):
        icon = "Icons/DrawMesh.png"
        text = "Single Select"
        super().__init__(icon=icon,text=text,parent=parent)
        
        


    def mousePressEvent(self,x,y,currentObject):
        """
        Selects the current object. If control is currently being pressed
        appends the current object to the list of selections.
        """
        
        if currentObject is not None:
            global qApp #this is a pointer to the manager of the whole application
                        #it is needed to check if keyboard keys are pressed
            
            if qApp.keyboardModifiers() == Qt.ControlModifier:
                if currentObject not in self.selections:
                    self.selections.append(currentObject)
                
            else:
                self.selections = [currentObject]

        
        
        
class BoxSelect(SelectTool):
    """
    Allows the user to draw a box which will then select all objects contained within that box
    """
    
    
    selectionRequest = pyqtSignal(object)
    
    def __init__(self, parent):
        icon = "Icons/rect.png"
        text = "Box Select"
        
        super().__init__(icon=icon,text=text,parent=parent)
        
        self.drawingBox = False
        self.box = None
        
    def mousePressEvent(self,x,y,currentObject):
       
        self.drawingBox = True
        self.selections = []
        self.firstPos = (x,y)
        
    def mouseMoveEvent(self, x, y, currentObject):
        
        if self.drawingBox:
           
            
            # QRectFs are defined by 4 float values,
            # The first two give the minimum x and y coordinates
            # The second two give the width and height of the rectangle
            x1 = min(x, self.firstPos[0])
            y1 = min(y, self.firstPos[1])

            width = abs(x-self.firstPos[0])
            height = abs(y-self.firstPos[1])
            
            self.box = QRectF(x1*self.scaleX,y1*self.scaleY,width*self.scaleX,height*self.scaleY)
            
            self.selections = []
            self.selectionRequest.emit(self.selectionFilter)
            
    def selectionFilter(self, selectable):
        
        """
        Returns true if all the nodes within selectable are within the 
        rectangle created by this tool.
        
        Any selectable object that meets this criteria should be set to selected while this tool 
        is in use.
        """
        
        for node in selectable.nodes:
            if self.box.contains(node.x*self.scaleX,node.y*self.scaleY):
                next
            else:
                return 
        if selectable not in self.selections:
            self.selections.append(selectable)
        return True
        
    
            
            
    def drawStuff(self,painter,*args,**kwargs):
        if self.drawingBox and self.box:
            painter.setPen(QColor(51,104,255)) 
            painter.setBrush(QColor(45,100,245,50))
            painter.drawRect(self.box)
        
        
    def mouseReleaseEvent(self, *args):
        self.drawingBox =False
        self.box = None
        

            
    
        
        
        
        
        
  
        

