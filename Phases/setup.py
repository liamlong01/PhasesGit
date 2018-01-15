# -*- coding: utf-8 -*-

"""
A script that can be run to build the python project as a standalone executable.
"""
import sys
import matplotlib
import numpy
import mpl_toolkits

from cx_Freeze import setup, Executable

path_platforms = ( r"C:\Users\ljl432\Downloads\WinPython-64bit-3.4.4.2Qt5\python-3.4.4.amd64\Lib\site-packages\PyQt5\plugins\platforms" )
resources = (r"C:\Users\ljl432\Documents\SpyderPhases\Resources")
numpyCore = (r"C:\Users\ljl432\Downloads\WinPython-64bit-3.4.4.2Qt5\python-3.4.4.amd64\Lib\site-packages\numpy\core")


build_exe_options = {
                        "packages"  :   ["os","matplotlib","mpl_toolkits","numpy","platform","ctypes","math","copy"],
                        "excludes"  :   ["tkinter","matplotlib.backends_tkagg"],
                        "includes"  :   ['sys', 'PyQt5', 'PyQt5.QtCore', 'PyQt5.QtGui', 'PyQt5.QtWidgets', 'os',"matplotlib.backends.backend_qt5agg"],
                        "include_files" : [path_platforms,resources,numpyCore]}
base = None 
#if sys.platform == "win32":
 #  base = "Win32GUI"
    
    
setup( name = "PHASES",
       version = "0.1",
       description = "PHASES with GUI capabilities",
       options = {"build_exe": build_exe_options},
       executables = [Executable("QtPHASES.py",base=base)] )
