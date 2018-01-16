from ctypes import *
import platform


class PhasesCaller():

    """
    A class to open and interface with phases library.
    Current implementation only works with .dll files, which are only used on windows.
    Slight adjustments can probably be made to handle other types like .so files
    """
    def __init__(self):
        #TODO: its possible that
        #self.libHandle = windll.kernel32.LoadLibraryA('Resources/phases.dll')
        #self.dll = WinDLL(None, handle = self.libHandle)
      
        #this is where other platform specific code would go 
        #right now there is only windows support for both 32 and 64-bit
        if platform.system() == 'Windows':
            if platform.architecture()[0] == '64bit':
                self.dll = WinDLL('Phases/Resources/phases64.dll')
            elif platform.architecture()[0] =='32bit':
                self.dll = WinDLL('Phases/Resources/phases32.dll')
        

        self.fileData = self.dll['Pyfiledata']
        self.fileData.restype = None
        self.fileData.argtypes = [c_int, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p]
    
        self.Adda = self.dll['PyAdda']
        self.Adda.restype = c_double
        self.Adda.argtypes = [c_int,c_int,c_int,c_int,c_double,
                              c_int,c_int,c_int,c_int,c_int]
        
    def phasesMain(self, files,steps):

        cfiles=files.copy()
    
        for key in cfiles.keys():
            assert (cfiles[key] !=''), '%s File not specified' %(key)
            cfiles[key] = cfiles[key].encode('utf-8')
        try:
            self.fileData(c_int(steps), cfiles['Mesh'],cfiles['Boundary Conditions'],cfiles['Initial Conditions'],cfiles['Output'],cfiles['Project'])
        except:
            pass
        
    def testAdda(self, a,b,c,d,e,f,g,h,i,j):
        if a < 30:
            print(a,b,c,d,e,f,g,h,i,j)
        return self.Adda(a,b,c,d,e,f,g,h,i,j)
    
    
        
        

   # def __del__(self):
    #    del self.dll
     #   windll.kernel32.FreeLibrary(self.libHandle)
        
    #def PhasesAdda(self,w1,w2,w3,w4,w5,nnp,nel,nx,ny,nsrf):

    #    self.Adda(self,w1,w2,w3,w4,w5,nnp,nel,nx,ny,nsrf)

    


