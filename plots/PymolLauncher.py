#import math

#import pylab
#import matplotlib
from AnnoteFinder import * 
import os 
from pymol import cmd

"""
Inherents from AnnoteFinder
"""
class PymolLauncher(AnnoteFinder):
    def set_native(self, native):
        self.native = native 
    
    def set_pdb_dir(self, pdb_dir):
        self.pdb_dir = pdb_dir
        
    def select_best(self):
        self.drawSpecificAnnote(self.data[0][2])
    
    """
    This is a derived class from AnnoteFinder! 
	
    Overwritten draw annote, extended to start pymol and show predicted/native structures
    
    Edit this file! Implement a functionallity that shows the superimposed native structure and the prediction in pymol
    you can use PyMol commands you already know in python style ( cmd.your_command(arguments) ) to achieve this! 
    
    Make sure that your structures are differently colored, in cartoon representation and zoomed in!
    
    """
    def drawAnnote(self, axis, x, y, annote):
        
        # EDIT: remove last markers
        for markers in self.drawnAnnotations.itervalues():
            for m in markers: # tuple
                #m.set_visible(not m.get_visible())
                m.remove()
        self.axis.figure.canvas.draw()
        self.drawnAnnotations = {}
        
        if (x,y) in self.drawnAnnotations:
            markers = self.drawnAnnotations[(x,y)] 
            for m in markers:
                m.set_visible(not m.get_visible())
            self.axis.figure.canvas.draw()
        
        else:
            """
            Mark data point and show data
            """
            t = axis.text(x,y, "(%3.2f, %3.2f)"%(x,y), )
            m = axis.scatter([x],[y], marker='d', c='r', zorder=100)
            self.drawnAnnotations[(x,y)] =(t,m)
            self.axis.figure.canvas.draw()

        """Your code here!"""
        cmd.delete('all')
        # load pdb files
        pdb_ext = '.pdb'
        cmd.load(self.native)
        cmd.load(os.path.join(self.pdb_dir, annote + pdb_ext))
        # label
        cmd.pseudoatom('foo')
        cmd.hide('all')
        #cmd.label('foo', '"'+annote+'"')
        # color
        native_name = os.path.splitext(os.path.basename(self.native))[0]        
        cmd.color('green', native_name)        
        cmd.set('label_color', 'green', native_name)
        cmd.color('blue', annote)        
        cmd.set('label_color', 'blue', annote)
        cmd.set('label_color', 'blue', 'foo')
        # background
        cmd.bg_color("white")
        cmd.set("depth_cue", 0)
        cmd.set("ray_trace_fog", 0)
        # view
        cmd.show(representation='cartoon')
        # alignment
        cmd.super(annote, native_name)
        # zoom
        cmd.zoom()        
        