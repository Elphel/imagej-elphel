#!/usr/bin/env python3

'''
/**
 * @file escher_pattern_class.py
 * @brief escher pattern generator
 * @par <b>License</b>:
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
'''

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import numpy as np

from matplotlib.patches import Arc
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
from matplotlib.patches import Polygon
from matplotlib.patches import Wedge
from matplotlib.patches import PathPatch

import matplotlib as mpl

from matplotlib import collections
from matplotlib.path import Path

import copy
import math

import sys

'''
requirements:
  
  - matplotlib, version 3.0.3
  - BUG?: matplotlib, version 2.0.2 draws a little padding around which is not desired

install:

  - sudo pip3 install matplotlib
  - sudo apt install python3-tk
  
'''

class Escher_Pattern:

  # units are mms but coordinates are in pt to compare with php script
  MM2PT = 72.0/25.4
  opts_k = {'facecolor':'black','edgecolor':'r','linewidth':0}
  opts_w = {'facecolor':'white','edgecolor':'r','linewidth':0}
  opts_b = {'facecolor':'blue', 'edgecolor':'r','linewidth':0}
  opts_r = {'facecolor':'red', 'edgecolor':'r','linewidth':0}

  # init
  def __init__(self,
               width    = 270,  # width in units
               height   = 210,  # height in units
               lpm      = 50,   # lines per meter
               escher   = 2.0,  # curvature coefficient
               rotate   = 5,    # degrees
               units    = 'mm',
               ):

    self.width    = int(width*self.MM2PT)
    self.height   = int(height*self.MM2PT)
    self.lpm      = lpm
    self.escher   = escher
    self.angle   = rotate
    self.units = units

    self.basename = 'escher-pattern'
    self.basename += '-ESCHER'+str(self.escher)
    self.basename += '-LPM'+str(self.lpm)
    self.basename += '-ROT'+str(self.angle)
    self.basename += '-PAGE_WIDTH'+str(width)
    self.basename += '-PAGE_HEIGHT'+str(height)

    self.pdf_name = self.basename+".pdf"

    plt.autoscale(tight=True)
    plt.axis('off')
    plt.margins(0.0)
    #plt.rcParams["figure.figsize"] = [60,120]

    self.fig, self.ax = plt.subplots()
    self.fig.set_size_inches(self.width/72, self.height/72)

    self.ax.get_xaxis().set_visible(False)
    self.ax.get_yaxis().set_visible(False)
    self.ax.set_aspect('equal')

    self.ax.spines['top'].set_visible(False)
    self.ax.spines['right'].set_visible(False)
    self.ax.spines['left'].set_visible(False)
    self.ax.spines['bottom'].set_visible(False)

    self.ax.set_ylim([0,self.height])
    self.ax.set_xlim([0,self.width])

    self.ax.set_facecolor('white')

    self.rotation = mpl.transforms.Affine2D().rotate_deg(-self.angle) + self.ax.transData

    #print("init done")


  # generate and place a patch
  def generate_cell(self,x,y,tpts,template,halfAngle,r,ba):

    for k,v in enumerate(template):

      # even-even black cell
      vcp = copy.copy(v)
      if   (type(v)==Wedge):
        vcp.set_center((tpts[k][0]+x,tpts[k][1]+y))
      elif (type(v)==Rectangle):
        vcp.set_xy((tpts[k][0]+x,tpts[k][1]+y))

      self.ax.add_patch(vcp)
      vcp.set_transform(self.rotation)

      # do clipping
      if   (type(v)==Wedge):

        sin0 = math.sin(math.radians(halfAngle))
        cos0 = math.cos(math.radians(halfAngle))

        if ba[k]==0:
          w0 = r - r*cos0
          h0 = 2*r*sin0
          x0 = tpts[k][0] + r - w0
          y0 = tpts[k][1] - h0/2
        if ba[k]==180:
          w0 = r - r*cos0
          h0 = 2*r*sin0
          x0 = tpts[k][0] - r
          y0 = tpts[k][1] - h0/2
        if ba[k]==270:
          w0 = 2*r*sin0
          h0 = r - r*cos0
          x0 = tpts[k][0] - w0/2
          y0 = tpts[k][1] - r
        if ba[k]==90:
          w0 = 2*r*sin0
          h0 = r - r*cos0
          x0 = tpts[k][0] - w0/2
          y0 = tpts[k][1] + r - h0


        # correction, so the shadow line from Rectange is not seen
        if ba[k]==0:
          x0 -= w0
          w0 *= 2
        if ba[k]==180:
          w0 *= 2
        if ba[k]==270:
          h0 *= 2
        if ba[k]==90:
          y0 -= h0
          h0 *= 2

        cp = Rectangle((x0,y0),w0,h0,**self.opts_w, lw=0)
        vcp2 = copy.copy(cp)

        xy = vcp2.get_xy()
        vcp2.set_xy((xy[0]+x,xy[1]+y))
        path = vcp2.get_path()
        transform = vcp2.get_transform()
        path = transform.transform_path(path)
        vcp2 = PathPatch(path, fc='none', ec='none', lw=0)
        self.ax.add_patch(vcp2)
        # rotate here
        vcp2.set_transform(self.rotation)

        vcp.set_clip_path(vcp2)


  # generate the whole pattern
  def generate(self):

    side = 500/self.lpm*self.MM2PT

    # escher pattern
    if (self.escher>0):

      # no rounding
      Size = side
      qSize = Size/4
      hSize = Size/2
      a = self.escher*(math.sqrt(2)-1.0)
      r = (a*a+1)/(2*a)*qSize
      r2 = r*r
      h = math.sqrt(r2-qSize*qSize)
      dc = 2*qSize-h
      halfAngle = math.degrees(math.atan(qSize/h))
      center = dc+r

      ba = [
        None,
        0,
        180,
        270,
        90,
        180,
        0,
        90,
        270
      ]

      tpts = [
          [center-hSize,   center-hSize],
          [center-Size+dc, center-qSize],
          [center-dc,      center+qSize],
          [center-qSize,   center+Size-dc],
          [center+qSize,   center+dc],
          [center+Size-dc, center+qSize],
          [center+dc,      center-qSize],
          [center+qSize,   center-Size+dc],
          [center-qSize,   center-dc]
        ]

      template = [
        Rectangle( tpts[0],  Size, Size, **self.opts_k),
        Wedge(     tpts[1],  r, ba[1]-halfAngle, ba[1]+halfAngle, **self.opts_w),
        Wedge(     tpts[2],  r, ba[2]-halfAngle, ba[2]+halfAngle, **self.opts_k),
        Wedge(     tpts[3],  r, ba[3]-halfAngle, ba[3]+halfAngle, **self.opts_w),
        Wedge(     tpts[4],  r, ba[4]-halfAngle, ba[4]+halfAngle, **self.opts_k),
        Wedge(     tpts[5],  r, ba[5]-halfAngle, ba[5]+halfAngle, **self.opts_w),
        Wedge(     tpts[6],  r, ba[6]-halfAngle, ba[6]+halfAngle, **self.opts_k),
        Wedge(     tpts[7],  r, ba[7]-halfAngle, ba[7]+halfAngle, **self.opts_w),
        Wedge(     tpts[8],  r, ba[8]-halfAngle, ba[8]+halfAngle, **self.opts_k)
      ]

    # checker board
    else:

      halfAngle = 0
      r = 0

      ba = [
        None
      ]

      tpts = [
        [0,0]
      ]

      template = [
        Rectangle(tpts[0],side,side, **self.opts_k)
      ]

    # calc how much more is needed for the rotation
    abs_angle_rad = math.radians(abs(self.angle))
    extra_w = self.height*math.tan(abs_angle_rad) # + side
    extra_h = self.width *math.tan(abs_angle_rad) # + side
    # would like to pass (0,0)
    extra_w = int(extra_w/(2*side)+1)*2*side
    extra_h = int(extra_h/(2*side)+1)*2*side

    a = np.arange(-extra_w, self.width +extra_w, 2*side)
    b = np.arange(-extra_h, self.height+extra_h, 2*side)

    for x, y in [(x,y) for x in a for y in b]:
      self.generate_cell(x     ,y     ,tpts,template,halfAngle,r,ba)
      self.generate_cell(x+side,y+side,tpts,template,halfAngle,r,ba)


  # test 
  def test(self):

    x = [0, 100,   0]
    y = [0,   0, 100]
    self.ax.fill(x, y)

    circ = Circle((0,0),0.01,color="red")
    #circ.set_transform(t_end)
    self.ax.add_patch(circ)
    circ = Circle((100,0),1,color="red")
    #circ.set_transform(t_end)
    self.ax.add_patch(circ)
    circ = Circle((0,100),1,color="red")
    #circ.set_transform(t_end)
    self.ax.add_patch(circ)
    circ = Circle((50,50),1,color="red")
    #circ.set_transform(t_end)
    self.ax.add_patch(circ)

  # save function
  def save(self):

    pp = PdfPages(self.pdf_name)
    self.fig.tight_layout(pad=0)

    #plt.show()

    self.fig.savefig(pp,format='pdf',bbox_inches='tight',pad_inches=0)
    #self.fig.savefig(pp,format='pdf',pad_inches=0)
    pp.close()



if __name__ == "__main__":

  #ep = Escher_Pattern("test.pdf", escher=2.0, lpm=50, rotate=10)
  #http://192.168.0.137/escher/escher_pattern.php?PAGE_WIDTH=1524&PAGE_HEIGHT=3048&LPM=2.705449885575893&ROTATE=14.036243467
  #ep = Escher_Pattern(width= 1524, height= 3048, escher=0, lpm=2.705449885575893, rotate=10)
  ep = Escher_Pattern(width= 1524, height= 3048, escher=2, lpm=2.705449885575893, rotate=14.036243467)
  #ep = Escher_Pattern(width= 160, height= 320, escher=2, lpm=50, rotate=13)
  ep.generate()
  ep.save()
