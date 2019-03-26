#!/usr/bin/env python3

'''
/**
 * @file escher_pattern.py
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

import matplotlib as mpl

from matplotlib import collections
from matplotlib.path import Path

import copy
import math

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

  # init
  def __init__(self,
               filename = "pattern.pdf",
               width    = 270,  # width in units
               height   = 210,  # height in units
               lpm      = 50,   # lines per meter
               escher   = 2.0,  # curvature coefficient
               rotate   = 5,    # degrees
               units    = 'mm',
               ):
    self.filename = filename
    self.width    = int(width*self.MM2PT)
    self.height   = int(height*self.MM2PT)
    self.lpm      = lpm
    self.escher   = escher
    self.angle   = rotate
    self.units = units

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

    print("init done")


  # rotate clock-wise
  def rotate(self,deg):

    t_start = self.ax.transData
    t = mpl.transforms.Affine2D().rotate_deg(-deg)
    t_end = t_start + t
    for x in self.ax.patches + self.ax.collections:
      x.set_transform(t_end)


  def generate(self):

    side = 500/self.lpm*self.MM2PT

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
        Rectangle( tpts[0],  Size, Size,                    facecolor="k",linewidth=0,edgecolor="r"),
        Wedge(     tpts[1],  r,   0-halfAngle,  0+halfAngle,facecolor="w",linewidth=0,edgecolor="r"),
        Wedge(     tpts[2],  r, 180-halfAngle,180+halfAngle,facecolor="k",linewidth=0,edgecolor="r"),
        Wedge(     tpts[3],  r, 270-halfAngle,270+halfAngle,facecolor="w",linewidth=0,edgecolor="r"),
        Wedge(     tpts[4],  r,  90-halfAngle, 90+halfAngle,facecolor="k",linewidth=0,edgecolor="r"),
        Wedge(     tpts[5],  r, 180-halfAngle,180+halfAngle,facecolor="w",linewidth=0,edgecolor="r"),
        Wedge(     tpts[6],  r,   0-halfAngle,  0+halfAngle,facecolor="k",linewidth=0,edgecolor="r"),
        Wedge(     tpts[7],  r,  90-halfAngle, 90+halfAngle,facecolor="w",linewidth=0,edgecolor="r"),
        Wedge(     tpts[8],  r, 270-halfAngle,270+halfAngle,facecolor="k",linewidth=0,edgecolor="r")
      ]

    else:

      tpts = [
        [0,0]
      ]

      template = [
        Rectangle(tpts[0],side,side,facecolor="k",linewidth=0,edgecolor="r")
      ]

    # calc how much more is needed for the rotation
    abs_angle_rad = math.radians(abs(self.angle))
    extra_w = self.height*math.tan(abs_angle_rad) + side
    extra_h = self.width *math.tan(abs_angle_rad) + side

    a = np.arange(-extra_w, self.width +extra_w, 2*side)
    b = np.arange(-extra_h, self.height+extra_h, 2*side)

    for x, y in [(x,y) for x in a for y in b]:

      for k,v in enumerate(template):

        # even-even black cell
        vcp = copy.copy(v)
        if   (type(v)==Wedge):
          vcp.set_center((tpts[k][0]+x,tpts[k][1]+y))
        elif (type(v)==Rectangle):
          vcp.set_xy((tpts[k][0]+x,tpts[k][1]+y))
        self.ax.add_patch(vcp)

        # odd-odd black cell
        vcp = copy.copy(v)
        if   (type(v)==Wedge):
          vcp.set_center((tpts[k][0]+x+side,tpts[k][1]+y+side))
        elif (type(v)==Rectangle):
          vcp.set_xy((tpts[k][0]+x+side,tpts[k][1]+y+side))
        self.ax.add_patch(vcp)

    # now rotate
    self.rotate(self.angle)


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

    pp = PdfPages(self.filename)
    self.fig.tight_layout(pad=0)

    #plt.show()

    self.fig.savefig(pp,format='pdf',bbox_inches='tight',pad_inches=0)
    #self.fig.savefig(pp,format='pdf',pad_inches=0)
    pp.close()



if __name__ == "__main__":

  #ep = Escher_Pattern("test.pdf", escher=2.0, lpm=50, rotate=10)
  #http://192.168.0.137/escher/escher_pattern.php?PAGE_WIDTH=1524&PAGE_HEIGHT=3048&LPM=2.705449885575893&ROTATE=14.036243467
  ep = Escher_Pattern("test.pdf", width= 1524, height= 3048, escher=2.0, lpm=2.705449885575893, rotate=14.036243467)
  ep.generate()
  ep.save()










