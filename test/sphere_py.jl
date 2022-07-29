# wrapping python code into Julia
using PyCall

#  @pyinclude("sphere.py") # it works for running the script

# magic incantation translated to Julia
np = pyimport("numpy")
@pyimport math as mathh; cos = mathh.cos; sin = mathh.sin; pi = mathh.pi


## define a Python function usable in Julia
py"""
import numpy as np
from math import cos, sin, pi, exp, atan
import matplotlib.pyplot as plt
import matplotlib as mpl

## versions:
# numpy==1.20.1
# matplotlib==3.3.4

#    d_sphere_P(r:np.array,θ:np.array,φ:np.array,x0:float,y0:float,z0:float,μ0:float)
#    In the spherical geometry, e.g. droplet, this function computes the distance from any
#    point M:(r,θ,φ) in the sample to the surface in direction to some point P:(x0,y0,z0).
#    Note: M is given in spherical coordinates and P in Cartesian coordinates.
#    Oz is along the photon beam. Center of the droplet is in the center of coordinates.
#    θ is taken between Oz and the project of M onto the plane xOz, φ between M and Oy, and r = √(x^2+y^2+z^2)
#    μ0 is the radius of the sphere
def d_sphere_P(r, θ, φ, x0, y0, z0, μ0):
    dp = np.zeros((len(r), len(θ), len(φ)))
    for i, R in enumerate(r):
        for j, Ө in enumerate(θ):
            for k, Φ in enumerate(φ):
                ## M coordinates from spherical to cartesian system
                xm = R*cos(Ө)*sin(Φ)
                ym = R*sin(Ө)*sin(Φ)
                zm = R*cos(Φ)

                ## compute all cos/sin
                cosω = (x0-xm) / ((x0-xm)**2 + (y0-ym)**2)**.5
                sinω = (y0-ym) / ((x0-xm)**2 + (y0-ym)**2)**.5
                    
                cosβ = ((x0-xm)**2 + (y0-ym)**2)**.5 / ((x0-xm)**2 + (y0-ym)**2 + (z0-zm)**2)**.5
                sinβ = (z0-zm) / ((x0-xm)**2 + (y0-ym)**2 + (z0-zm)**2)**.5
                    
                cosα = sin(Φ) * cosβ * (
                    cos(Ө) * cosω + sin(Ө) * sinω) + cos(Φ) * sinβ

                dp[i,j,k] =  (μ0**2 - R**2 * (1-cosα**2))**.5 - R * cosα
    return dp


#    sphere_gain_H(r:np.array,θ:np.array,φ:np.array,x0:float,y0:float,z0:float,μ0:float,λe:float)
#    Compute the volume integrales (exact for piecewise linear gain)
#    H_{n,j,k} = ∭ e_i(r)e_j(θ)e_k(φ) e^{-d_P(M)/λe} r^2drdθdφ
#    The arrays r, θ and φ are the discretization subdivisions
#    P:(x0,y0,z0) is the point in Cartesian coordinates used for the computation of the distance d_P
#    μ0 is the radius of the sphere
#    λe is the attenuation length in the Beer-Lambert model
#    Note: does not include the profile of the beam light, but it can be add easily if known
def sphere_gain_H(r, θ, φ, x0, y0, z0, μ0, λe):
    Ari = 0.5 * r**2 * np.diff(r, append=(2*r[-1] - r[-2]))
    Aθj = 0.5 * np.diff(θ, append=(2*θ[-1] - θ[-2]))
    Aφk = 0.5 * np.diff(φ, append=(2*φ[-1] - φ[-2]))
    
    dp = d_sphere_P(r, θ , φ, x0, y0, z0, μ0)
    H_rθφ = np.zeros((len(Ari), len(Aθj), len(Aφk)))
    H_r = np.zeros((len(Ari), len(Aθj), len(Aφk)))
    for i in range(len(Ari)):
        for j in range(len(Aθj)):
            for k in range(len(Aφk)):
                H_rθφ[i,j,k] = exp(-1*dp[i,j,k]/λe)
                H_r[i,j,k] = Ari[i] * Aθj[j] * Aφk[k] * H_rθφ[i,j,k]
    
    return H_r, H_rθφ, Ari, Aθj, Aφk
"""

## alias
d_sphere_P_jl = py"""d_sphere_P"""


## some variables
μ0 = 20.0 # radius
λe = 2.0e-3 # μ0; # EAL
L = 3 * μ0

## discretisation
N = 10
K = 256
J = 64 

r = np.linspace(μ0-10*λe, μ0, 2) 
θ = np.linspace(0, 2*pi, J)
φ = np.linspace(0, pi, J)

# magic angle: atan(sqrt(2.0),1.0)
## near the analyzer
x0_near = py"""21.0*2.0**.5"""
y0_near = 0.0
z0_near = 21.0

## far away from the analyzer
x0_far = py"""200.0*2.0**.5"""
y0_far = 0.0
z0_far = 200.0

## compute the distance
DD = d_sphere_P_jl(r, θ , φ, x0_near, y0_near, z0_near, μ0)

py"""
## transform the coordinate system: spherical->cartesian
x = np.zeros((len(r), len(θ), len(φ)))
y = np.zeros((len(r), len(θ), len(φ)))
z = np.zeros((len(r), len(θ), len(φ)))

for i, R in enumerate(r):
    for j, Ө in enumerate(θ):
        for k, Φ in enumerate(φ):
            x[i,j,k] = R*cos(Ө)*sin(Φ)
            y[i,j,k] = R*sin(Ө)*sin(Φ)
            z[i,j,k] = R*cos(Φ)

## create axes and plot
fig = plt.figure(figsize=(7,7))
ax = plt.axes(projection='3d')

ax.scatter(x0_near, y0_near, z0_near, c='red', s=20)
ax.scatter(x, y, z, c=DD)
ax.set_xlim(-30,30)
ax.set_ylim(-30,30)
ax.set_zlim(-30,30)

#ax.view_init(90, 180)
fig.show()
"""
