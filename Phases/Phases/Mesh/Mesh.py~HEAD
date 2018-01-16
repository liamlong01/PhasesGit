# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 14:30:04 2016

@author: ljl432
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 14:20:40 2016

@author: ljl432
"""
import DataPrep
import PhasesCaller
import numpy as np
from math import sqrt
from copy import deepcopy

def loadMeshFromFile(files):
    """
    This function will initalize and return a Mesh object, using data given from text files 
    specified in files. Files needs to be a dictionary with keys 'Project,'Mesh', 'Boundary Conditions', 'Initial Conditions'
    which are associated with directories to the appropriate data files. this function will use this information to
    create a mesh with the appropriate properties, geometry, boundary conditions, and inital conditions
    
    """
    
    mesh = Mesh()
    mesh.params = DataPrep.sampleF(files['Project'])
    
    anodes,aelem,nbon1,nbon2,nbon3,nbon4,elbon,elems,xy = DataPrep.meshF(files['Mesh'])
    bnodes_surf,ele_surf,concen_b,temp_b,uv_b,vv_b = DataPrep.boundaryF(files['Boundary Conditions'], elbon)
    ic_con,ic_temp,ic_uv,ic_vv = DataPrep.initialconF(files['Initial Conditions'], aelem)
    
    #interpreting and organizing mesh data
    mesh.nnp = anodes
    mesh.nel = aelem
    mesh.nsrf = elbon
    mesh.nx = nbon3
    mesh.ny = nbon4

    
    for coords in xy:   #xy is structured as a list of lists: 
        newNode = (Node(coords[0],coords[1]))
        mesh.addNode(newNode)
        
        
        

    for elem in elems:
        elemNodes=[]
        for node in elem:
            elemNodes.append(mesh.nodes[node-1])
        newElem = Element(elemNodes)
        mesh.addElement(newElem)
       
        
    #interperting and organizing boundary conditions
    for boundary in bnodes_surf:
        index = bnodes_surf.index(boundary)
        refElem = mesh.elements[ele_surf[index]-1]
        refElem.isBoundary = True 
        
        boundary = Boundary([mesh.nodes[boundary[0]-1],mesh.nodes[boundary[1]-1]])
        boundary.element = refElem
        refElem.boundaries.append(boundary)

        mesh.addBoundary(boundary)
        

    def addBoundaryConditions(kind, fromList):
        for i in range(len(fromList)):
            vals = [float(j) for j in fromList[i][:3]]
            mesh.boundaries[int((i/2))].conditions[kind].append(vals)
        
    addBoundaryConditions('concentration',concen_b)
    addBoundaryConditions('temperature',temp_b)
    addBoundaryConditions('u-velocity',uv_b)
    addBoundaryConditions('v-velocity',vv_b)

    for i in range(len(vv_b)):
        mesh.boundaries[int(i/2)].orientation = vv_b[3]
    

    
    #    for i in range(len(concen_b)):
 #       concenVals = concen_b[i][:3]
  #      concenVals[2] = concenVals[2]/(1+1e-32)
   #     
    #    mesh.boundaries[int((i/2))].concentrations.append(concenVals)
        
    
    #interpreting and organizing inital conditions
    def addInitialConditions(kind,fromList,faca,facb):
        for i in range(len(fromList)):
            val = faca*float(fromList[i][1]) +facb
            mesh.elements[i].initialConditions[kind] = val

    addInitialConditions('concentration',ic_con,1,0)
    
    faca = 1/(mesh.params['tmax']-mesh.params['tmin'] +1e-32)
    facb = -mesh.params['tmin']/(mesh.params['tmax']-mesh.params['tmin'] +1e-32)
    addInitialConditions('temperature',ic_temp,faca,facb)
    
    faca = 1/(mesh.params['ur']+1e-32)
    addInitialConditions('u-velocity',ic_uv,faca,0)
    addInitialConditions('v-velocity',ic_vv,faca,0)
    
    return mesh
    
def saveMeshToFile(mesh,directory,name):
   
    ###### write to project file
    
    
   
    prj = open(directory + '\\%s.prj'%name, 'w')
    
    #writing parameters
    orderOfParams = ['se','ao','tk','tvk','tstp','df','lr','wr','tol','ur','pcs','pcl','tmlt','tsol','bt','bs','tmin','tmax','vsc','prs','prl','pds0','pdl0','pks','pkl','phs','phl','pl']
    for param in orderOfParams:
        prj.write('%s\n'%mesh.params[param])
    prj.write('0\n')
    prj.write('1\n')

    #writing file names
    prj.write(directory + '\\%s.prj\n'%name)
    prj.write(directory + '\\mesh_%s.dat\n'%name)
    prj.write(directory + '\\ic_%s.dat\n'%name)
    prj.write(directory + '\\bc_%s.dat\n'%name)
    prj.write(directory + '\\output_%s.dat\n'%name)
    prj.write(directory + '\\batch_%s.dat\n'%name)
    prj.write(directory + '\\custum_%s.dat\n'%name)
    
    prj.close()
    
 
    ###### write to mesh file
    meshFile = open(directory + '\\mesh_%s.dat'%name, 'w')
   
    # first line states the #nodes #elements, dimensions, #boundaries
    firstLine = '%s,%s,%s,%s,%s,%s,%s\n' %(len(mesh.nodes),len(mesh.elements),mesh.nx,mesh.ny,mesh.nx,mesh.ny,len(mesh.boundaries))
    meshFile.write(firstLine)
    
    # writing the nodes in each element
    for element in mesh.elements:
        elemNodes = ''
        for node in element.nodes:
            elemNodes += '%s,'%node.index
        meshFile.write(elemNodes[:-1]+'\n') #removes the last comma and replaces it with a new line before writing
    
    # writing xy coordinates for each node    
    for node in mesh.nodes:
        meshFile.write('%s,%s\n'%(node.x,node.y))
        
    meshFile.close()
    
    
    ##### write to boundary conditions file
    bcFile = open(directory + '\\bc_%s.dat'%name, 'w')
    
    # writing boundary nodes
    for boundary in mesh.boundaries:
        line = '%s,%s\n'%(boundary.nodes[0].index,boundary.nodes[1].index)
        bcFile.write(line)
    
    # writing boundary elements
    for boundary in mesh.boundaries:
        bcFile.write('%s\n'%boundary.element.index)
    
    # writing boundary conditions
    
    def writeBoundary(key):
        for boundary in mesh.boundaries:
            
            line = '%s,%s,%s,%s,%s\n'%(boundary.conditions[key][0][0],boundary.conditions[key][0][1],boundary.conditions[key][0][2], boundary.orientation,boundary.element.index)
            bcFile.write(line)
            
            line = '%s,%s,%s,%s,%s\n'%(boundary.conditions[key][1][0],boundary.conditions[key][1][1],boundary.conditions[key][1][2], boundary.orientation, boundary.element.index)
            bcFile.write(line)        
                     
    writeBoundary('concentration')
    writeBoundary('temperature')
    writeBoundary('u-velocity')
    writeBoundary('v-velocity')
               
    
    bcFile.close()
    
    #write to initial conditions file
    icFile= open(directory + '\\ic_%s.dat'%name, 'w')
    
    def writeInitial(key):
        for element in mesh.elements:
            print(key)
            cond = element.initialConditions[key]
            
           
            icFile.write('%s,%s,%s,%s,%s\n'%(element.index,cond,cond,cond,cond))
    
    
    writeInitial('concentration')
    writeInitial('temperature')
    writeInitial('u-velocity')
    writeInitial('v-velocity')
    
    icFile.close()
    
    #write to output file
    outputFile = open(directory + '\\output_%s.dat'%name, 'w')
    
    
    for i in range(int(mesh.params['tstp'])):
        for node in mesh.nodes:
          
            
            try:
                outputFile.write('%s,%s,%s,%s,%s,%s,%s,%s,%s, %s\n'%(node.index,
                                                             node.temperature[i],
                                                             node.concentration[i],
                                                             mesh.params['prl'],
                                                             node.uVelocity[i],
                                                             node.vVelocity[i],
                                                             node.wVelocity[i],
                                                             1-node.liquidFrac[i],
                                                             node.liquidFrac[i],
                                                            node.entropy[i]))
                
            except IndexError: # the number of timesteps is probably greater then the number of timesteps specified
                next          # all calculated outputs have been written so break the loop
        
    
    outputFile.close()
            

class Mesh(object):
    """
    A class that contains all of the necessary information required to run a phases simulation.
    Has the following properties:
        
        nodes       - a list of references to all nodes as Node objects contained in the mesh
        elements    - a list of references to all elements as Element objects contained in the mesh
        boundaries  - a list of references to all boundaries as Boundary objects contained in the mesh
        
        nnp     - the total number of nodes
        nel     - the total number of elements
        nx      - number of nodes along the x direction
        ny      - number of nodes along the y  direction
        nsrf    - the number of boundary surfaces
        
    And the following methods:
        getXY           - returns two numpy arrays, one for X coordinates, one for Y coordinates of all nodal points in the mesh
        getTemperature  - returns a numpy array giving the temperature of each nodal point at the given timestep
        getUV           - returns two numpy arrays giving the u-velocity and v-velocity of each nodal point at the given timestep
        getFLIQ         - returns a numpy array giving the liquid fraction at each nodal point at the given timestep
        getPressure     - returns a numpy array giving the pressure of each nodal point at the given timestep
        
        setParam        - sets the specified parameter to the given value
        
        sendToAdda      - sends parameters, mesh geometry, boundary conditions, and initial conditions to the Adda function in the PHASES library.
                          Then starts the simulation for the given number of time steps and stores the output data in the mesh Node objects.    
                          This data can then be accessed through the getXY, getTemperature, etc. methods

    """
    
    def __init__(self):
        self.nodes = []
        self.elements = []
        self.boundaries = []
        self.params= { }


        self.nnp = None
        self.nel = None
        self.nx = None
        self.ny = None
        self.nsrf = None
        self.dx = None
        self.dy = None
        #self.steps = 5
        
        self.tempRange = []
        self.concentrationRange = []
        self.fliqRange = []
        self.uVelRange = []
        self.vVelRange = []
        self.velMagRange = []
        self.pressureRange = []

    def addNode(self, node):
        """
        Adds the Node node to the mesh’s nodes if it is not already there.
        """
        node.x = round(node.x,8)
        node.y = round(node.y,8)
        if node not in self.nodes:   
            found = False
            for otherNode in self.nodes:
                if  node.x < otherNode.x and not found:
                    index = self.nodes.index(otherNode) 
                    self.nodes.insert(index, node)
                    node.index = index+1
                    found = True
                elif  node.x == otherNode.x and not found:   
                    if node.y < otherNode.y:
                        index = self.nodes.index(otherNode) 
                        self.nodes.insert(index, node)
                        node.index = index+1
                        found = True
                elif found:
                    otherNode.index+=1
            if not found:
                self.nodes.append(node)
                node.index = len(self.nodes)
            return node

        i = self.nodes.index(node)
        return self.nodes[i]
            
    def subNode(self, node):
        """
        Subtracts the Node node from the mesh’s nodes, if it is there.
        """
        if node not in self.nodes:
            return
        index = node.index - 1
        for node in self.nodes[index:]:
            node.index -= 1
        self.nodes.pop(index)
        
            
    def addElement(self, element):
        """
        Adds the Element element to the mesh’s elements if it is not already there.
        """
        if element not in self.elements:
            
            for node in element.nodes:
                self.addNode(node)
            
            found = False
            for otherElem in self.elements:
                if element.nodes[0].index < otherElem.nodes[0].index and not found:
                    index = self.elements.index(otherElem)
                    self.elements.insert(index, element)
                    element.index = index+1
                    found = True
                elif found:
                    otherElem.index+=1
                    
            if not found:
                self.elements.append(element)
                element.index = len(self.elements)
            
            
                
            
    def subElement(self, element):
        """
        Subtracts the Element element from the mesh’s elements, if it is there.
        """
        if element not in self.elements:
            return
        index = element.index - 1
        for elem in self.elements[index:]:
            elem.index -= 1
        self.elements.pop(index)
       
            
    def addBoundary(self, boundary):
        """
        Adds the Boundary boundary to the mesh’s boundaries if it is not already there.
        """
        if boundary not in self.boundaries:
            found = False
            for otherBound in self.boundaries:
                if boundary.element.index < otherBound.element.index and not found:
                    index = self.boundaries.index(otherBound)
                    self.boundaries.insert(index,boundary)
                    boundary.index = index+1
                    found = True
                elif found:
                    otherBound.index+=1
            if not found:    
                self.boundaries.append(boundary)
                boundary.index = len(self.boundaries)
            
    def subBoundary(self, boundary):
        """
        Subtracts the Boundary boundary from the mesh’s boundaries, if it is there.
        """
        if boundary not in self.boundaries:
            return
        index = boundary.index - 1
        for bound in self.boundaries[index:]:
            bound.index -= 1
        self.boundaries.pop(index)
            
    def getXY(self):
        """
        Returns the xy-coordinates of the mesh’s nodal points as a 2D numpy array.
        """
        x = np.zeros(len(self.nodes))
        y = np.zeros(len(self.nodes))
        for node in self.nodes:
            x[node.index-1] = node.x
            y[node.index-1] = node.y
  
        return self.split(x),self.split(y)
        
    def getTemperature(self,timeStep):
        """
        Returns the temperature values of the mesh’s nodal points 
        at the given timeStep as a 2D numpy array.
        """
        temp = np.zeros(len(self.nodes))
        for node in self.nodes:
            temp[node.index-1] = node.temperature[timeStep-1]
        return self.split(temp)

    def getUV(self,timeStep):
        """
        Returns the u-Velocity and v-Velocity values of the mesh’s nodal points 
        at the given timeStep as a two 2D numpy arrays.
        """
        u = np.zeros(len(self.nodes))
        v = np.zeros(len(self.nodes))
        for node in self.nodes:
            u[node.index-1] = node.uVelocity[timeStep-1]
            v[node.index-1] = node.vVelocity[timeStep-1]
            
        return  self.split(u), self.split(v)

    def getFLIQ(self,timeStep):
        """
        Returns the liquid fraction values of the mesh’s nodal points 
        at the given timeStep as a 2D numpy array.
        """
        fliq = np.zeros(len(self.nodes))
        for node in self.nodes:
            fliq[node.index-1] = node.liquidFrac[timeStep-1]
        return self.split(fliq)
        
    def getPressure(self,timeStep):
        """
        Returns the pressure values of the mesh’s nodal points 
        at the given timeStep as a 2D numpy array.
        """
        p = np.zeros(len(self.nodes))
        for node in self.nodes:
            p[node.index-1] = node.pressure[timeStep-1]
        return self.split(p)
        
    def split(self, data):
         """
         This method is used to format data in numpy arrays to properly fit the mesh geometry.
         It takes a 1D numpy array and transforms it into a 2D array with the appropriate dimensions.
         """
         length = len(data)
         for i in range(len(self.nodes)):
             if self.nodes[i].x != self.nodes[i+1].x:
                 break
         wanted_parts = int(len(self.nodes)/(i+1))
         return np.transpose(np.array([ data[i*length // wanted_parts: (i+1)*length // wanted_parts] for i in range(wanted_parts) ] ))
            
    def setParam(self, param, toWhat):
        self.params[param] = toWhat

    def calcEntropy(self):
        if self.dx is not None and self.dy is not None:
           
            for i in range(1,self.params['tstep']+1):
                u,v = self.getUV(i);
            
                gradU = np.gradient(u,[self.dx, self.dy])
                dudy = gradU(2)
                
                gradV = np.gradient(v,[self.dx, self.dy])
                dvdx = gradV(1)
          
                mu = self.params['vsc']*self.params['prl']

                T = self.getTemperature(i)
                entropy = mu*np.divide((np.multiply(dudy,dudy)+np.multiply(dvdx,dvdx)) , np.multiply(T,T))
               
                for j in range(length(entropy)):
                    print(entropy[j])
                    self.nodes[j].entropy[i-1] = entropy[j]
            
            
    def sendToAdda(self):
    
        
        X, Y = self.getXY()

        
        
        self.ny, self.nx = X.shape #shape returns the number of rows and the number of columns in the X array
            
        
        caller = PhasesCaller.PhasesCaller()
        
        #TODO: testAdda is currenetly being used for testing and should be changed to just Adda in the final implementation
        tempAdda = caller.testAdda  # Adda will be called frequently so this makes a good shorthand
        
        Adda = lambda a,b,c,d,e: tempAdda(a,b,c,d,e,len(self.nodes),len(self.elements),self.nx,self.ny,len(self.boundaries))
        
        #Adda = lambda a,b,c,d,e: tempAdda(a,b,c,d,e,len(self.nodes),len(self.elements))
        #export params
        Adda(1, 1, 1, 1, self.params['se'])
  
        
        Adda(2, 1, 1, 1, self.params['ao'])
        Adda(3, 1, 1, 1, self.params['tk'])
        Adda(4, 1, 1, 1, self.params['tvk'])
        Adda(5, 1, 1, 1, self.params['tstp'])
        Adda(6, 1, 1, 1, self.params['df'])
        Adda(7, 1, 1, 1, self.params['lr'])
        Adda(8, 1, 1, 1, self.params['wr'])
        Adda(9, 1, 1, 1, self.params['tol'])
        Adda(10, 1, 1, 1, self.params['ur'])
        Adda(11, 1, 1, 1, self.params['pcs'])
        Adda(12, 1, 1, 1, self.params['pcl'])
        Adda(13, 1, 1, 1, self.params['tmlt'])
        Adda(14, 1, 1, 1, self.params['tsol'])
        Adda(15, 1, 1, 1, self.params['bt'])
        Adda(16, 1, 1, 1, self.params['bs'])
        Adda(17, 1, 1, 1, self.params['tmin'])
        Adda(18, 1, 1, 1, self.params['tmax'])
        Adda(19, 1, 1, 1, self.params['vsc'])
        Adda(20, 1, 1, 1, self.params['prs'])
        Adda(21, 1, 1, 1, self.params['prl'])
        Adda(22, 1, 1, 1, self.params['pds0'])
        Adda(23, 1, 1, 1, self.params['pdl0'])
        Adda(24, 1, 1, 1, self.params['pks'])
        Adda(25, 1, 1, 1, self.params['pkl'])
        Adda(26, 1, 1, 1, self.params['phs'])
        Adda(27, 1, 1, 1, self.params['phl'])
        Adda(28, 1, 1, 1, self.params['pl'])
        Adda(29, 1, 1, 1, self.params['prl'])
        
        
        #export mesh data
        
        for i in range(len(self.elements)):
            for j in range(len(self.elements[i].nodes)):
                Adda(31, i+1, j+1, self.elements[i].nodes[j].index, 1) #passes the index of each node to PHASES one at a time in the format of a 2d array
            
                     
        for node in self.nodes:
            
            #clears any data that may have been stored from previous simulations
            node.clearData()

            #this sends x and y coordinates of each node to the phases simulation
            Adda(32,node.index,1,1,node.x) # passes x of each node one at a time
            Adda(33,node.index,1,1,node.y) # passes  y of each node one at a time
                
        #export boundary conditions
        for boundary in self.boundaries:
           
            Adda(41,boundary.index,1,boundary.nodes[0].index, 1) #passes first node of boundary to Adda
            Adda(41,boundary.index,2,boundary.nodes[1].index, 1) #passes second node of boundary to Adda
        
            Adda(42,boundary.index, boundary.element.index, 1, 1)
     
            conditions = deepcopy(boundary.conditions)
     
          
            
            #next three lines get called twice per boundary for each of the 4 ics (total of 8 times per boundary - 24 Adda calls per boundary?????)
            Adda(43, 1, 1, boundary.index*2-1, conditions['concentration'][0][0])
            Adda(43, 1, 2, boundary.index*2-1, conditions['concentration'][0][1])
            Adda(43, 1, 3, boundary.index*2-1, conditions['concentration'][0][2])
            
            Adda(43, 1, 1, boundary.index*2, conditions['concentration'][1][0])
            Adda(43, 1, 2, boundary.index*2, conditions['concentration'][1][1])
            Adda(43, 1, 3, boundary.index*2, conditions['concentration'][1][2])
            
            temp = conditions['temperature'][0]
            temp[2] = (temp[2]-temp[1]*self.params['tmin'])/float((self.params['tmax']-self.params['tmin']))
            Adda(43, 2, 1, boundary.index*2-1, conditions['temperature'][0][0])
            Adda(43, 2, 2, boundary.index*2-1, conditions['temperature'][0][1])
            Adda(43, 2, 3, boundary.index*2-1, conditions['temperature'][0][2])
            
            temp = conditions['temperature'][1]
            temp[2] = (temp[2]-temp[1]*self.params['tmin'])/float((self.params['tmax']-self.params['tmin']))
            Adda(43, 2, 1, boundary.index*2, conditions['temperature'][1][0])
            Adda(43, 2, 2, boundary.index*2, conditions['temperature'][1][1])
            Adda(43, 2, 3, boundary.index*2, conditions['temperature'][1][2])
            
            u = conditions['u-velocity'][0]
            u[2] = u[2]/float(self.params['ur'])
            Adda(43, 3, 1, boundary.index*2-1, conditions['u-velocity'][0][0])
            Adda(43, 3, 2, boundary.index*2-1, conditions['u-velocity'][0][1])
            Adda(43, 3, 3, boundary.index*2-1, conditions['u-velocity'][0][2])
            
            u = conditions['u-velocity'][1]
            u[2] = u[2]/float(self.params['ur'])
            Adda(43, 3, 1, boundary.index*2, conditions['u-velocity'][1][0])
            Adda(43, 3, 2, boundary.index*2, conditions['u-velocity'][1][1])
            Adda(43, 3, 3, boundary.index*2, conditions['u-velocity'][1][2])
            
            v = conditions['v-velocity'][0]
            v[2] = v[2]/float(self.params['ur'])
            Adda(43, 4, 1, boundary.index*2-1, conditions['v-velocity'][0][0])
            Adda(43, 4, 2, boundary.index*2-1, conditions['v-velocity'][0][1])
            Adda(43, 4, 3, boundary.index*2-1, conditions['v-velocity'][0][2])
            
            v = conditions['v-velocity'][1]
            v[2] = v[2]/float(self.params['ur'])
            Adda(43, 4, 1, boundary.index*2, conditions['v-velocity'][1][0])
            Adda(43, 4, 2, boundary.index*2, conditions['v-velocity'][1][1])
            Adda(43, 4, 3, boundary.index*2, conditions['v-velocity'][1][2])
            
            
        #export initial conditions
        for elem in self.elements:
            Adda(51,elem.index,1,1,elem.initialConditions['concentration'])
            Adda(52,elem.index,1,1,elem.initialConditions['temperature'])
            Adda(53,elem.index,1,1,elem.initialConditions['u-velocity'])
            Adda(54,elem.index,1,1,elem.initialConditions['v-velocity'])
            
        
        #execute   
        Adda(60, 1, 1, 1, 1);

    
        self.tempRange = []
        self.concentrationRange = []
        self.fliqRange = []
        self.uVelRange = []
        self.vVelRange = []
        self.velMagRange = []
        self.pressureRange = []
    
        for i in range(int(self.params['tstp'])):
            Adda(0, i+1, 1, 1, 1)
            
        #######get output
        
      
            for node in self.nodes:
                dimTemp = Adda(100, 1, node.index, 1, 1)    
                temp = dimTemp*(self.params['tmax']-self.params['tmin']) + self.params['tmin']   
                
                self.setRange(self.tempRange, temp)    
                node.temperature.append(temp)
    
                concentration = Adda(100, 2, node.index, 1, 1)
                self.setRange(self.concentrationRange, concentration)  
                node.concentration.append(concentration)    

                uVel = Adda(100, 4, node.index, 1, 1)
                uVel = uVel*self.params['ur']
                self.setRange(self.uVelRange, uVel) 
                node.uVelocity.append(uVel)
    
                vVel = Adda(100, 5, node.index, 1, 1)
                vVel = vVel*self.params['ur']
                self.setRange(self.vVelRange, vVel) 
                node.vVelocity.append(vVel)
                
                self.setRange(self.velMagRange, sqrt(vVel**2+uVel**2) )
    
                wVel = Adda(100, 5, node.index, 1, 1)
                node.wVelocity.append(wVel*self.params['ur'])
                
                pressure = wVel * self.params['prl']*self.params['ur']*self.params['ur']
                node.pressure.append(pressure)
                self.setRange(self.pressureRange, pressure)
                
                fliq = Adda(100, 9, node.index, 1, 1)
                node.liquidFrac.append(fliq )
                self.setRange(self.fliqRange, fliq)
            
            
               
        del caller 
        self.calcEntropy()
   
    def setRange(self, rangeToSet, val):
        if rangeToSet:
            rangeToSet[0] = min(rangeToSet[0], val)
            rangeToSet[1] = max(rangeToSet[1], val)
        else: #range is empty and needs to be initialized
            rangeToSet.append(val)
            rangeToSet.append(val)
         

class Node(object):
    """
    An object that represents a nodal point within the mesh.
    
        x           - x-coordinate
        y           - y-coordinate
        index       - an integer label that is required by the phases library to run a simulation
        isBoundary  - a bool stating whether or not this object is on the boundary of the mesh
        
        Output data:
            temperature,concentration,uVelocity,vVelocity,wVelocity,pressure,liquidfrac
            
            these are all lists which contain the output data from a simulation. The nth index gives
            the value during the n+1 timestep. eg. temperature[1] is the calculated temperature during the 2nd timestep at this node
    """
    
    def __init__(self,x,y):
        self.x = x
        self.y = y
        self.index = None
        self.isBoundary = False
        
        self.temperature = []
        self.concentration = []
        self.uVelocity = []
        self.vVelocity = []
        self.wVelocity = []
        self.pressure = []
        self.liquidFrac = []
        self.entropy = []

    def clearData(self):
        """
        Resets and clears the output data of this node
        """
        self.temperature = []
        self.concentation = []
        self.uVelocity =[]
        self.vVelocity =[]
        self.wVelocity=[]
        self.pressure =[]
        self.liquidFrac=[]
        self.entropy = []

    def __str__(self):
        return "Node object at (%s,%s) with index of %s." % (self.x,self.y,self.index)
    
    def __eq__(self, otherNode):
        if round(self.x,8)==round(otherNode.x,8) and round(self.y,8) == round(otherNode.y,8):
            return True
        else:
            return False
        
        
    
    
class Element(object):
    """
    A class representing an element in a mesh.
    
        nodes            - a list of 4 nodes that are the corners of this element in counter clockwise order
        index            - an integer label that is required by the phases library for the simulation
        isBoundary       - a bool stating whether or not this object is on the boundary of the mesh
        initalConditions - a dictionary with keys 'concentration', 'temperature', 'u-velocity','v-velocity'
                           that are associated with floats describing the inital condtions
    """
    
    def __init__(self,nodes):
        assert (len(nodes)==4), 'Error: Mesh element initialized with %n nodes. Element requires exactly 4 nodes.' %len(nodes)
        self.nodes = nodes
        self.index = None
        self.isBoundary = False
        self.boundaries = []
        self.orientation = ' '
        
        self.initialConditions = {
                                  'concentration': None,
                                  'temperature': None,
                                  'u-velocity': None,
                                  'v-velocity': None
                                      }                                  
    def __str__(self):
        info = "Element object with index of %s,nodes at (%s,%s),(%s,%s),(%s,%s),(%s,%s)." %(self.index, self.nodes[0].x,self.nodes[0].y,self.nodes[1].x,self.nodes[1].y,self.nodes[2].x,self.nodes[2].y,self.nodes[3].x,self.nodes[3].y)
        if self.boundaries:
            for boundary in self.boundaries:
                info = info + "Boundary %s." %boundary.index
        return info
                                  
   
        
class Boundary(object):
    """
    A class representing a boundary surface in a mesh.
    
        nodes            - a list of 2 nodes that are the ends of this boundary surface
        elements         - the element that contains this boundary surface
        index            - an integer label that is required by the phases library for the simulation
        conditions       - a dictionary with keys 'concentration', 'temperature', 'u-velocity','v-velocity'
                           that are associated with a list of length 2, as each boundary has two values for each type of condition
    """
    
    def __init__(self,nodes):
        assert (len(nodes)==2), 'Error: Mesh boundary initialized with %n nodes. Boundary requires exactly 2 nodes.' %len(nodes)
        self.nodes = nodes
        self.element = None
        self.index = None
        self.orientation = ' '
        
        self.conditions =   {
                            'concentration' : [],
                            'temperature': [],
                            'u-velocity': [],
                            'v-velocity': [],
                            }
                            
                    
    def __str__(self):
        return "Boundary object with nodes at (%s,%s),(%s,%s), temperature of %s, index of %s, element index of %s" %(self.nodes[0].x,self.nodes[0].y,self.nodes[1].x,self.nodes[1].y,self.conditions['temperature'],self.index,self.element.index)
      
