# -*- coding: utf-8 -*-
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

from Mesh.Mesh import Node
  

class NodeGenerator(QWidget):
    
    nodeEdited = pyqtSignal()
    nodeDeleted = pyqtSignal(Node)
    
    def __init__(self, display,parent = None):
        QWidget.__init__(self, parent)
        self.nodes = display.nodes
        layout = QVBoxLayout(self)
        
        ### Adding nodes
        addNode = QGroupBox('Add Nodes')
        addLayout = QGridLayout(addNode)
        
        self.xIn = QLineEdit()
        addLayout.addWidget(QLabel('x :'),0,0)
        addLayout.addWidget(self.xIn,0,1)
        x = lambda: float(self.xIn.text())
        
        self.yIn = QLineEdit()
        addLayout.addWidget(QLabel('y :'),1,0)
        addLayout.addWidget(self.yIn,1,1)
        y = lambda: float(self.yIn.text())

        add = QPushButton('Add')
        addLayout.addWidget(add,2,0)
        add.clicked.connect(lambda: self.addNode(x(),y()))
        addLayout.setRowStretch(4,1)
        
        
        #Delete Nodes
        subNode = QGroupBox('Delete Nodes')
        subLayout = QGridLayout(subNode)
        
        self.delNodeSelect = QSpinBox()
        subLayout.addWidget(QLabel('Node #:'),0,0)
        subLayout.addWidget(self.delNodeSelect,0,1)
        
        delButton = QPushButton('Delete Node')
        delButton.clicked.connect(lambda:self.delNode(self.delNodeSelect.value()))
        subLayout.addWidget(delButton,1,1)
        
        subLayout.setRowStretch(2,1)
        
        delAllButton = QPushButton('Clear All Nodes')
        delAllButton.clicked.connect(self.clearNodes)
        subLayout.addWidget(delAllButton,3,1)
        
        #rearrange nodes
        swapNode = QGroupBox('Rearrange Nodes')
        swapLayout = QGridLayout(swapNode)
        
        self.node1 = QSpinBox()
        swapLayout.addWidget(QLabel('Node #1:'),0,0)
        swapLayout.addWidget(self.node1,0,1)
        
        self.node2 = QSpinBox()
        swapLayout.addWidget(QLabel('Node #2:'),1,0)
        swapLayout.addWidget(self.node2,1,1)
        
        swap = QPushButton('Swap Nodes')
        swap.clicked.connect(lambda: self.swapNodes(self.node1.value(),self.node2.value()))
        swapLayout.addWidget(swap,2,1)
        
        
        self.position = QLabel('')
        layout.addWidget(self.position)
        
        layout.addWidget(addNode)
        layout.addWidget(subNode)
        layout.addWidget(swapNode)
        layout.addStretch(2)
        
        self.setIndexRange()
        
        
        
    def addNode(self, x, y):
        self.xIn.setText('')
        self.yIn.setText('')
        
        node = Node(x,y)
        if node in self.nodes:
            return node
            
        self.nodes.append(node)
        self.setIndexRange(maximum = len(self.nodes))
        
        self.nodeEdited.emit()

        return node
        
    def delNode(self,index):
        if index == 0:
            return
            
        node = self.nodes.pop(index-1)
        
        self.setIndexRange(maximum = len(self.nodes))
        self.nodeEdited.emit()
        self.nodeDeleted.emit(node)
        
    def swapNodes(self, index1,index2):
        if index1 == 0 or index2 == 0:
            return
        
        temp = self.nodes[index1-1]

        self.nodes[index1-1] = self.nodes[index2-1]
        self.nodes[index2-1] = temp

        self.nodeEdited.emit()
        
        
    def clearNodes(self):
        for i in range(len(self.nodes)):
            self.delNode(1)
       
    def setIndexRange(self, minimum=0,maximum=0):
        self.delNodeSelect.setMaximum(maximum)
        self.delNodeSelect.setMinimum(minimum)
        
        self.node1.setMaximum(maximum)
        self.node2.setMinimum(minimum)
        
        self.node1.setMaximum(maximum)
        self.node2.setMinimum(minimum)
        
    