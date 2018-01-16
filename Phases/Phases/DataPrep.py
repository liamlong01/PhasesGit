"""
This file contains functions that used to read text files
then sort their data into digestable chuncks which
are then used by Mesh.Mesh.loadMeshFromFile to 
initialize a mesh.
"""



import numpy as np 

#split-list takes a 1D array and transforms it into a 2D numpy array .
#alist is the array to be transformed
#wanted_parts 
def split_list(alist, wanted_parts=1):
    length = len(alist)
    return np.array([ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
                                            for i in range(wanted_parts) ] )


#in functions with argument name 'file_directory' remember to preceed the arguments with 'r' i.e outputF(r'C:\Users\eelh73\Desktop\DATA\output2.dat') 
'''
anodes:number of nodes
alem:number of elements
nboni:nodes in boundary i i.e: nbon1 is nodes in boundary 1 
elbon: number of elements in boundary 
elem: number of elements elements with their respective nodes in order 
x_cord,y_cord: x and y coordinates for nodes 
bnodes:number of boundary nodes 
ele_surf:number of elements in surface 
'''
#Reads output files and stores solutions in new arrays, Returns solution arrays. Arguments: file directory and number of nodes 
#New arrays are meant to be used for generating visualization tools in Python 
#uv is a 2d array with x-y velocities components
def outputF(file_directory,anodes,i): 
    
    solution_file=open(file_directory,'r') 
    wholefsol =solution_file.readlines()
    wholefsol = wholefsol[anodes*i:anodes*(i+1)] 
    
    wholefsol =np.array([line.split(',') for line in wholefsol])
    wholefsol = wholefsol.astype(np.float)
    
    nodes = wholefsol[:,0]
    temp  = wholefsol[:,1]
    
    concen= wholefsol[:,2] 
    rho   = wholefsol[:,3]
    
    uvel  = wholefsol[:,4]
    vvel  = wholefsol[:,5]
    wvel  = wholefsol[:,6]
    uv    = wholefsol[:,4:6]
    
    fsol  = wholefsol[:,7]
    fliq  = wholefsol[:,8]
    
        
    solution_file.close()
    
    #print(nodes,'\n',uv) # temp,'\n', concen,'\n',rho,'\n', uvel,'\n', vvel, '\n', wvel,'\n', fsol,'\n', fliq)
    return nodes,temp,concen,rho,uvel,vvel,wvel,fsol,fliq,uv
    
    
#Reads mesh files and stores data in new arrays. Returns mesh arrays and important data such as # of nodes in a mesh and a 2d array of x-y coordinates. 
#First function to be executed. 
#Arguments: file directory. 
def meshF(file_directory):
    
    meshf=open(file_directory,'r')
    wholefm=meshf.readlines()
    
    line1=[int(i) for i in wholefm[0].split(',')]
     
    anodes = line1[0]
    aelem  = line1[1] 
    nbon1  = line1[2] 
    nbon2  = line1[3] 
    nbon3  = line1[4] 
    nbon4  = line1[5]
    elbon  = line1[6]
    
    elem=[]    #placeholder  
    #elem=np.zeros((aelem,4))
    
    
    for line in wholefm[1:aelem+1]:
        nodes = line.split(',')
        for i in range(len(nodes)):
            nodes[i] = int(nodes[i])
        elem.append(nodes)
            
    #i=0
    #while i < aelem: 
        
     #   j=0
      #  new_array2=wholefm[i+1].split(',')
         
       # for digit in new_array2:
            
        #    elem[i][j]=int(digit) 
         #   j+=1
            
       # i+=1  
       
    
    #
    xy = np.array([line3.split(',') for line3 in wholefm[(aelem+1):(aelem+anodes+1)]])
    #convert data to floats
    xy = xy.astype(np.float)
    
    
    
    meshf.close() 
   
    return anodes,aelem,nbon1,nbon2,nbon3,nbon4,elbon,elem,xy


     
    
 #Reads boundary condition file and stores data in new arrays. Return these arrays.  
 #arguments: file direcoty and number of elements in the boundary.
 #A(dG/dn)+B(G)=C where G is the scalar quantity 
def boundaryF(file_directory,elbon):
    boundf=open(file_directory,'r') 
    wholefb=boundf.readlines()  
    i=0
    i1=0
    bnodes_surf=[]
    ele_surf=[] 
    concen_b=[]
    temp_b=[]
    uv_b=[]
    vv_b=[]


    for elemt in wholefb[0:elbon]:
      
        new_array4=elemt.split(',') 
        for i in range(len(new_array4)):
            new_array4[i] = int(new_array4[i])
        bnodes_surf.append(new_array4)
            
        
        
    for elemt2 in wholefb[elbon:(elbon*2)]:
        ele_surf.append(int(elemt2.split(',')[0]))
        
        
    for set1 in wholefb[(elbon*2):(elbon*4)]:
        concen_b.append(set1.split(',')) 
    for set2 in wholefb[(elbon*4):(elbon*6)]: 
        temp_b.append(set2.split(',')) 
    for set3 in wholefb[(elbon*6):(elbon*8)]:
        uv_b.append(set3.split(','))
    for set4 in wholefb[(elbon*8):]:
        vv_b.append(set4.split(','))  
    boundf.close() 
    #print(concen_b,'\n',temp_b,'\n',uv_b,'\n',vv_b)
    return bnodes_surf,ele_surf,concen_b,temp_b,uv_b,vv_b
    

#Read intial condition files and saves data in new arrays. 
#Arguments: files directory and # of elements 
# group 1: Concentration,2: Temperature,3:U velocity,4: V vecolity 
def initialconF(file_directory,aelem): 
    icf=open(file_directory,'r') 
    wholefic=icf.readlines()
    ic_con=np.zeros((aelem,5))
    ic_temp=np.zeros((aelem,5))
    ic_uv=np.zeros((aelem,5))
    ic_vv=np.zeros((aelem,5)) 
    i1,i2,i3,i4=0,0,0,0
    for con in wholefic[:aelem]: 
        new_array5=con.split(',') 
        j1=0
        for digit1 in new_array5:
            ic_con[i1][j1]=float(digit1) 
            j1+=1
        i1+=1
    for temp in wholefic[aelem:(aelem*2)]: 
        new_array6=temp.split(',') 
        j2=0
        for digit2 in new_array6:
            ic_temp[i2][j2]=float(digit2) 
            j2+=1
        i2+=1
    for uv in wholefic[(aelem*2):(aelem*3)]: 
        new_array7=uv.split(',') 
        j3=0 
        for digit3 in new_array7: 
            ic_uv[i3][j3]=float(digit3)
            j3+=1
        i3+=1
    for vv in wholefic[(aelem*3):]:
        new_array8=vv.split(',') 
        j4=0 
        for digit4 in new_array8: 
            ic_vv[i4][j4]=float(digit4) 
            j4+=1
        i4+=1
    icf.close() 
    #print(ic_con,'\n',ic_temp,'\n',ic_uv,'\n', ic_vv) 
    return ic_con,ic_temp,ic_uv,ic_vv
    
    
#Reads sample files (parameter file) and stores data in a dictionary named 'param'. To access a parameter in the dictionary such as number of time steps ('tstp') do: param['tstp']. 
def sampleF(file_directory): 
    samplef=open(file_directory,'r')
    wholefsam=samplef.readlines() 
    param={} 
    param['se']=float(wholefsam[0]) #  there are 3 main solution blocks in cntrl, C(concentration equation), T (temperature) and U/V/P (velocities/ pressures),
# se indicates which block is solved, numbered 1,2,3 or 4 (all together), Options: 
#se = 1(C), 2(T), 3(UVP), 12(CT), 13(CUVP), 23(TUVP), 4(CTUVP)
    param['ao']=float(wholefsam[1]) # an addition capability of droplet flows was added later, ao = 2 means that the droplet flow solver is added, use ao = 1 (without droplets)
    param['tk']=float(wholefsam[2]) #number of times cycling through the temperature equation loop until a convergence tolerance is reached 
    param['tvk']=float(wholefsam[3]) # same as tk, but instead the number of time through the temperature/velocity loop 
    param['tstp']=float(wholefsam[4]) # number of time steps
    param['df']=float(wholefsam[5]) #size of timestep
    param['lr']=float(wholefsam[6]) #domain height 
    param['wr']=float(wholefsam[7])#domain width 
    param['tol']=float(wholefsam[8])#residual tolerance 
    param['ur']=float(wholefsam[9])#boundary velocity reference 
    param['pcs']=float(wholefsam[10]) #concentration of solidus 
    param['pcl']=float(wholefsam[11]) #concentration of liquidus     
    param['tmlt']=float(wholefsam[12])# melting point 
    param['tsol']=float(wholefsam[13]) # temperature of solidus 
    param['bt']=float(wholefsam[14]) #thermal expansivity 
    param['bs']=float(wholefsam[15]) #sollutal expansivity 
    param['tmin']=float(wholefsam[16]) # minimum temperature 
    param['tmax']=float(wholefsam[17]) #maximum temperature 
    param['vsc']=float(wholefsam[18]) #kinematic viscosity 
    param['prs']=float(wholefsam[19]) #solid density 
    param['prl']=float(wholefsam[20]) #fluid density 
    param['pds0']=float(wholefsam[21]) #solid diffusivity 
    param['pdl0']=float(wholefsam[22]) #fluid diffusivity 
    param['pks']=float(wholefsam[23]) #solid conductivity 
    param['pkl']=float(wholefsam[24]) #fluid conductivity 
    param['phs']=float(wholefsam[25]) #solid specific heat
    param['phl']=float(wholefsam[26])  #fluid specific heat
    param['pl']=float(wholefsam[27]) #latent heat of fusion
    samplef.close() 
    
    if param['ur'] == 0:
        param['ur'] = param['pkl']/(param['prl']*param['phl']*param['lr'])
    param['df'] =param['df']/(param['lr']/param['ur'])
    
    
    return param
    
    
    #return se,ao,tk,tvk,tstp,df,lr,wr,tol,ur,pcs,pcl,tmlt,tsol,bt,bs,tmin,tmax,vcs,prs,prl,pds0,pdl0,pks,pkl,phs,phl,pl 
# prepares data for plotting 
def dataprep(mesh_directory,output_directory, i): 
    anodes,aelemn,bon1,nbon2,nbon3,nbon4,elbon,elem,xy=meshF(mesh_directory) 
    
    i=i-1
    nodes,temp,concen,rho,uvel,vvel,wvel,fsol,fliq,uv=outputF(output_directory,anodes,i)
    
    for i in range(anodes): 
        if xy[i][0] != xy[i+1][0]: 
            break 
            
    dimen=int(anodes/(i+1))

    
    X = split_list(xy[:,0],dimen)
    Y = split_list(xy[:,1],dimen)

    U = split_list(uvel,dimen)
    V = split_list(vvel,dimen)
    
    C = split_list(concen,dimen)
    RHO = split_list(rho,dimen)
    
    FSOL = split_list(fsol,dimen)
    FLIQ = split_list(fliq,dimen)
    
    T = split_list(temp,dimen)
    

    return X,Y,T,U,V,C,RHO,FSOL,FLIQ
