import slitlessutils as su
import numpy as np
from scipy.special import expit
from astropy.modeling.models import Sersic2D


import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

Vmax = 250.
Rscl = 5.
rscl = 10.

inc = 45.
pa = 67.


# plt.style.use('/Users/rryan/.matplotlib/idl.mplstyle')

def V(R): return np.full_like(R, Vmax)
# V=lambda R: Vmax*(r/10)
def V(R): return Vmax/(1+np.exp(-R/Rscl))
# V=lambda R: 2*(expit(R/Rscl)-0.5)*Vmax


# R=np.arange(0,50,1.)
# plt.plot(R,V(R))
# plt.show()
#

# y,x=np.mgrid[0:101:0.1,0:101:0.1]
y, x = np.mgrid[-50:51:0.1, -50:51:0.1]

# x0=50.
# y0=50.#
#
# y-=x0
# x-=y0

r = np.maximum(np.hypot(x, y), 0.1)
sint = y/r
cost = x/r
sinp = np.sin(pa*np.pi/180)
cosp = np.cos(pa*np.pi/180)

sina = sinp*cost+cosp*sint


# p=np.arcsin(y/r)
# sinp=np.sin(np.arcsin(y/r)+pa*np.pi/180.)


sini = np.sin(inc*np.pi/180.)


width = 80
aspect = np.cos(inc*np.pi/180.)
height = width*aspect

# v=np.zeros_like(x)


v = V(r)*sina*sini

vmin = np.amin(v)
vmax = np.amax(v)
N = 10
f = 0.99

dv = f*(vmax-vmin)/(N-1)/2.

vlevels = np.arange(f*vmin, f*vmax+dv, dv)
print(len(vlevels))
# vlevels=np.logspace(0,np.log10(vmax-10),num=10)
# vlevels=np.geomspace(100,0.95*vmax,num=5)


fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 5), sharex=True, sharey=True)
a2.set_xlabel('x (pix)')
a1.set_xlabel('x (pix)')
a1.set_ylabel('y (pix)')

extent = [-50, 50, -50, 50]
img1 = a1.imshow(v, cmap=mpl.colormaps['RdBu_r'], extent=extent)
cnt = a1.contour(v, vlevels, colors='k', linestyles='dotted', extent=extent)

ell1 = Ellipse((0, 0), width=width, height=height, angle=90-pa,
               edgecolor='k', facecolor='none', clip_on=True)
ell2 = Ellipse((0, 0), width=width, height=height, angle=90-pa,
               edgecolor='k', facecolor='none', clip_on=True)

a1.add_patch(ell1)

img1.set_clip_path(ell1)


for c in cnt.collections:
    c.set_clip_path(ell1)


sers = Sersic2D(amplitude=1, r_eff=5, n=1, x_0=0, y_0=0, ellip=1-aspect, theta=np.radians(90-pa))


img2 = a2.imshow(sers(x, y), extent=extent,
                 norm=mpl.colors.LogNorm(vmin=1e-5, vmax=5),
                 cmap=mpl.colormaps['gray_r'])
a2.add_patch(ell2)
img2.set_clip_path(ell2)


plt.show()


# plt.plot(R,V)
# plt.show()
