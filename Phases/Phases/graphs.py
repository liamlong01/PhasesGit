import matplotlib.pyplot as plt 
from matplotlib import cm
import numpy as np 
from mpl_toolkits.axes_grid1 import make_axes_locatable

"""
 protoype graphing functions
 When implemented in the final application
 they should be used within an mplCanvas.PlotType class
"""   
    

def streamlines(X,Y,U,V, fig, sub):
    np.seterr(divide='ignore',invalid='ignore')
    magn=np.sqrt((U**2)+(V**2))
    strm = sub.streamplot(X,Y,-V,-U, color = magn,cmap=cm.autumn)
    
    divider=make_axes_locatable(sub)
    clbr=fig.colorbar(mappable = strm.lines)
    clbr.set_label(r"Magnitude Reference ($m/s$)")
    
    
    

#plots temperature contours 
def contourtemp(X,Y,T, fig, sub): 
    sub.set_title(' Temperature Contours')
    #plt.title('Temperature Contours in an 11X11-Element Mesh Simulation', fontsize=12.5) 
    
    cont=sub.contourf(X,Y,T,  cmap = cm.hot)
    divider=make_axes_locatable(sub)
    #cax=divider.append_axes('top','3%',pad='3%')
    clbr=fig.colorbar( mappable = cont)
    clbr.set_label(r'Temperature Reference $^\circ K$') 
   

def pressure(X,Y,P,fig,sub):
    sub.set_title(' Pressure Contours')
    #plt.title('Temperature Contours in an 11X11-Element Mesh Simulation', fontsize=12.5) 
    
    cont=sub.contourf(X,Y,P, cmap = cm.bone)
    divider=make_axes_locatable(sub)
    #cax=divider.append_axes('top','3%',pad='3%')
    clbr=fig.colorbar( mappable = cont)
    clbr.set_label(r'Pressure Reference') 
    
    
#plots isotherm lines 
def isothermp(X,Y,T,fig):

    fig.set_title('Isotherms in a Multiphase Flow Simulation ($^\circ K$)')
    fig.set_ylabel('Hot Wall', fontsize=12)
    cont=fig.contour(X,Y,T) 
    fig.clabel(cont,inline=1) 
 
    
    
    
#plots liquid fraction contours 
def  liquidfcont(X, Y, FLIQ,fig,sub):

    
  
    sub.set_title('Liquid Fraction Contours in a Multiphase Flow Simulation')
    sub.set_ylabel('Hot Wall', fontsize=12)
    V = [i/8.0 for i in range(9)]
  
    cont=sub.contourf(X,Y,FLIQ,V) 
    

    divider=make_axes_locatable(sub)
    #cax=divider.append_axes('right','3%',pad='3%')
    clbr=fig.colorbar(mappable=cont)
    clbr.set_label(r'Scalar Reference')
    
    return fig
 
def saveFigs(mesh_directory, output_directory, steps):
    
    for i in range(steps):
        X,Y,T,U,V,C,RHO,FSOL,FLIQ=dataprep(mesh_directory,output_directory, i)
        
        liquidfcont(X,Y, FLIQ)
        
        plt.savefig('fooLFC' + str(i) + '.png')
        
        isothermp(X,Y,T)
        
        plt.savefig('fooIso' + str(i) + '.png')
        
        contourtemp(X,Y,T)
        plt.savefig('fooTcont' + str(i) + '.png')
        
        vectorf(X,Y,U,V)
        plt.savefig('fooVect' + str(i) + '.png')
    
    
    
 







