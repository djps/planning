"""
InitialPosition.py
Copyright (C) 2013 David Sinden

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""

# Author: David Sinden <david.sinden@icr.ac.uk>
# Copyright (c) 2013, David Sinden / Institute of Cancer Research
# License: GPL.

"""
Information about the code.
"""

################################################################################

import numpy as np

from tvtk.api import tvtk
import vtk

import ICRcolors

from traits.api import HasTraits, Range, Instance, \
                    on_trait_change
from traitsui.api import View, Item, HGroup, Group
from tvtk.pyface.scene_editor import SceneEditor
from mayavi.tools.mlab_scene_model import \
                    MlabSceneModel
                    
from mayavi.core.ui.mayavi_scene import MayaviScene

from numpy import linspace, pi, cos, sin

################################################################################

def position(alpha1, beta1, gamma1):
    target = [270.0, 280.0, 180.0]
    N = 250
    R0 = np.zeros( (3,3) )
    alpha0=0.0
    beta0=0.0
    gamma0=0.0
    c1 = np.cos(alpha0)
    c2 = np.cos(beta0)
    c3 = np.cos(gamma0)
    s1 = np.sin(alpha0)
    s2 = np.sin(beta0)
    s3 = np.sin(gamma0)
    R0 = [[c2, -c3*s2, s2*s3 ], [ c1*s2, c1*c2*c3 - s1*s3,  -c3*s1 - c1*c2*s3 ], [ s1*s2, c1*s3 + c2*c3*s1,  c1*c3 - c2*s1*s3 ] ]
    
    centre      = np.zeros( (N+1,3) )
    new_centre0 = np.zeros( (N+1,3) )  
    new_centre1 = np.zeros( (N+1,3) )  
    
    f = open('../dat/Transducer0.dat', 'r')
    line=f.readline()
    words=line.split()
    # counter is set to 1, rather than zero to define centre point of transducer
    i=1
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    for line in f:
      words=line.split(',')
      # read coordinates and scale to mm
      centre[i,0] = np.float( words[0] )*10.0 
      centre[i,1] = np.float( words[1] )*10.0 
      centre[i,2] = np.float( words[2] )*10.0 
      # apply rotations to centre to form new_centre, by rotated transducer face
      new_centre0[i,:] = np.dot(R0, centre[i,:] )
      i=i+1
    f.close()
    
    centre[0,:] = [0.0,0.0, -np.linalg.norm( centre[1,:]-target ) ]
    new_centre0[0,:] = [0.0,0.0, -np.linalg.norm( new_centre0[1,:]-target ) ]
    
    R1 = np.zeros( (3,3) )
    c1 = np.cos(alpha1)
    c2 = np.cos(beta1)
    c3 = np.cos(gamma1)
    s1 = np.sin(alpha1)
    s2 = np.sin(beta1)
    s3 = np.sin(gamma1)
    R1 = [[c2, -c3*s2, s2*s3 ], [ c1*s2, c1*c2*c3 - s1*s3,  -c3*s1 - c1*c2*s3 ], [ s1*s2, c1*s3 + c2*c3*s1,  c1*c3 - c2*s1*s3 ] ]
    for i in np.arange(0,N+1):
	new_centre1[i,:] = np.dot(R1, new_centre0[i,:] ) + target
	
    return new_centre1

outer_ring = np.zeros( (32,1) )
outer_ring[0] = 8
for i in np.arange(1,31):
  outer_ring[i] = outer_ring[i-1] + 3.0*( 1.0 - (-1.0)**(i) ) + 5.0*( 1.0 + (-1.0)**(i) )
  
################################################################################

class Visualization(HasTraits):
  
    alpha1 = Range(low=0, high=round(2.0*np.pi, 2), value=round(np.pi,2) )
    beta1  = Range(low=0, high=round(2.0*np.pi, 2), value=round(np.pi,2) )
    gamma1 = Range(low=0, high=round(2.0*np.pi, 2), value=round(np.pi,2) )
    
    scene       = Instance(MlabSceneModel, ())
    scene1      = Instance(MlabSceneModel, ())
    scene2      = Instance(MlabSceneModel, ())
    
    editor=SceneEditor(scene_class=MayaviScene)

    def __init__(self):
        
        # Do not forget to call the parent's __init__
	HasTraits.__init__(self)
	
	target = [270.0, 280.0, 180.0]

	t0 = self.scene.mlab.pipeline.open('../dat/test_02_small.stl')
	t1 = self.scene.mlab.pipeline.surface(t0, opacity=1, color=ICRcolors.ICRblue)
	t2 = self.scene.mlab.pipeline.extract_edges(t1)
	self.plot = self.scene.mlab.pipeline.surface(t2, opacity=0.1, color=ICRcolors.ICRblue)
	
	s0 = self.scene.mlab.pipeline.open('../dat/ribs_02_small.stl')
	s1 = self.scene.mlab.pipeline.surface(s0, opacity=1, color=ICRcolors.ICRred)
	s2 = self.scene.mlab.pipeline.extract_edges(s1)
	self.plot = self.scene.mlab.pipeline.surface(s2, opacity=0.1, color=ICRcolors.ICRbrightred)
	
	self.plot = self.scene.mlab.points3d( target[0], target[1], target[2], mode='sphere', scale_factor=10, color=ICRcolors.cyan, figure=self.scene.mayavi_scene) 
        
	L=250
	centre  = np.zeros( (L+1,3) )
	centre0 = np.zeros( (L+1,3) )  
	centre1 = np.zeros( (L+1,3) )  
	
	f = open('../dat/Transducer0.dat', 'r')
	line=f.readline()
	words=line.split()
	# counter is set to 1, rather than zero to define centre point of transducer
	i=1
	f.readline()
	f.readline()
	f.readline()
	f.readline()
	for line in f:
	  words=line.split(',')
	  # read coordinates and scale to mm
	  centre[i,0] = np.float( words[0] )*10.0 
	  centre[i,1] = np.float( words[1] )*10.0 
	  centre[i,2] = np.float( words[2] )*10.0 
	  i=i+1
	f.close()

	centre[0,:] = [0.0, 0.0, -np.linalg.norm( centre[1,:]-target ) ]

        m=10
	centre_alpha = np.zeros( (m+1,3) )
	centre_beta = np.zeros( (m+1,3) )
	centre_gamma = np.zeros( (m+1,3) )
	for j in np.linspace(0, m, num=m+1, endpoint=True):
	  i = int(j)
	  alpha_dummy = j*np.pi/(4.0*m)
	  beta_dummy = round(np.pi/20.0,2) 
	  gamma_dummy = round(np.pi/20.0,2) 
	  c1_dummy = np.cos(alpha_dummy)
	  c2_dummy = np.cos(beta_dummy)
	  c3_dummy = np.cos(gamma_dummy)
	  s1_dummy = np.sin(alpha_dummy)
	  s2_dummy = np.sin(beta_dummy)
	  s3_dummy = np.sin(gamma_dummy)
	  R1_dummy = [[c2_dummy, -c3_dummy*s2_dummy, s2_dummy*s3_dummy ], [ c1_dummy*s2_dummy, c1_dummy*c2_dummy*c3_dummy - s1_dummy*s3_dummy, -c3_dummy*s1_dummy - c1_dummy*c2_dummy*s3_dummy ], [ s1_dummy*s2_dummy, c1_dummy*s3_dummy + c2_dummy*c3_dummy*s1_dummy,  c1_dummy*c3_dummy - c2_dummy*s1_dummy*s3_dummy ] ]
	  centre_alpha[i,:] = np.dot(R1_dummy, centre[0,:] ) + target  
	for j in np.linspace(0, m, num=m+1, endpoint=True):
	  i = int(j)
	  alpha_dummy = round(np.pi/2.0,2) 
	  beta_dummy = j*np.pi/(4.0*m)
	  gamma_dummy = round(np.pi/2.0,2) 
	  c1_dummy = np.cos(alpha_dummy)
	  c2_dummy = np.cos(beta_dummy)
	  c3_dummy = np.cos(gamma_dummy)
	  s1_dummy = np.sin(alpha_dummy)
	  s2_dummy = np.sin(beta_dummy)
	  s3_dummy = np.sin(gamma_dummy)
	  R1_dummy = [[c2_dummy, -c3_dummy*s2_dummy, s2_dummy*s3_dummy ], [ c1_dummy*s2_dummy, c1_dummy*c2_dummy*c3_dummy - s1_dummy*s3_dummy, -c3_dummy*s1_dummy - c1_dummy*c2_dummy*s3_dummy ], [ s1_dummy*s2_dummy, c1_dummy*s3_dummy + c2_dummy*c3_dummy*s1_dummy,  c1_dummy*c3_dummy - c2_dummy*s1_dummy*s3_dummy ] ]
	  centre_beta[i,:] = np.dot(R1_dummy, centre[0,:] ) + target
	for j in np.linspace(0, m, num=m+1, endpoint=True):
	  i = int(j)
	  alpha_dummy = round(np.pi/2.0,2) 
	  beta_dummy = round(np.pi/2.0,2) 
	  gamma_dummy = j*np.pi/(4.0*m)
	  c1_dummy = np.cos(alpha_dummy)
	  c2_dummy = np.cos(beta_dummy)
	  c3_dummy = np.cos(gamma_dummy)
	  s1_dummy = np.sin(alpha_dummy)
	  s2_dummy = np.sin(beta_dummy)
	  s3_dummy = np.sin(gamma_dummy)
	  R1_dummy = [[c2_dummy, -c3_dummy*s2_dummy, s2_dummy*s3_dummy ], [ c1_dummy*s2_dummy, c1_dummy*c2_dummy*c3_dummy - s1_dummy*s3_dummy, -c3_dummy*s1_dummy - c1_dummy*c2_dummy*s3_dummy ], [ s1_dummy*s2_dummy, c1_dummy*s3_dummy + c2_dummy*c3_dummy*s1_dummy,  c1_dummy*c3_dummy - c2_dummy*s1_dummy*s3_dummy ] ]
	  centre_gamma[i,:] = np.dot(R1_dummy, centre[0,:] ) + target
	  
	self.plot = self.scene.mlab.plot3d( centre_alpha[:,0], centre_alpha[:,1], centre_alpha[:,2], opacity=1, tube_radius=2, tube_sides=8, color=ICRcolors.ICRred, figure=self.scene.mayavi_scene)
	
	self.plot = self.scene.mlab.plot3d( centre_beta[:,0], centre_beta[:,1], centre_beta[:,2], opacity=1, tube_radius=1, tube_sides=8, color=ICRcolors.ICRblue, figure=self.scene.mayavi_scene)
	
	self.plot = self.scene.mlab.plot3d( centre_gamma[:,0], centre_gamma[:,1], centre_gamma[:,2], opacity=1, tube_radius=1, tube_sides=8, color=ICRcolors.black, figure=self.scene.mayavi_scene)
	
	centres = position(self.alpha1, self.beta1, self.gamma1)
	x0=centres[:,0]
	y0=centres[:,1]
	z0=centres[:,2] 

	self.plot1 = self.scene.mlab.plot3d( [x0[0],target[0]], [y0[0],target[1]], [z0[0],target[2]], opacity=1, tube_radius=1, tube_sides=8, color=ICRcolors.sienna, figure=self.scene.mayavi_scene)
	
	self.plot2 = self.scene.mlab.points3d(x0, y0, z0, mode='sphere', scale_factor=3, color=ICRcolors.green, figure=self.scene.mayavi_scene )
	
	oa = self.scene.mlab.orientation_axes() 
	oa.marker.set_viewport(0,0,0.4,0.4) 
	
    @on_trait_change('alpha1,beta1,gamma1, scene.activated')
    
    def update_plot1(self):
      
      	# set axes
	xlab = self.scene.mlab.xlabel("x") 
	ylab = self.scene.mlab.ylabel("y")
	zlab = self.scene.mlab.zlabel("z")
	axes = self.scene.mlab.axes()
      
        new_centres = position(self.alpha1, self.beta1, self.gamma1)
        target = [270.0, 280.0, 180.0]
        x1=new_centres[:,0]
        y1=new_centres[:,1]
        z1=new_centres[:,2]
        self.plot1.mlab_source.set(x=[x1[0],target[0]], y=[y1[0],target[1]], z=[z1[0],target[2]])
        new_centres1 = position(self.alpha1, self.beta1, self.gamma1)
        x2=new_centres1[:,0]
        y2=new_centres1[:,1]
        z2=new_centres1[:,2]
        self.plot2.mlab_source.set(x=x2, y=y2, z=z2)
        
        print self.alpha1, self.beta1, self.gamma1

    # the layout of the dialog created
    view = View( Item('scene', editor=SceneEditor(scene_class=MayaviScene),height=350, width=300, show_label=True),
                HGroup('_', 'alpha1', 'beta1', 'gamma1',), resizable=True, )
                
################################################################################

visualization = Visualization()

visualization.configure_traits()