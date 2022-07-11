# wrapping python code into Julia
using PyCall

#  @pyinclude("sphere.py") # it works for running the script

# magic incantation translated to Julia
np = pyimport("numpy")
pd = pyimport("pandas")
@pyimport math as mathh; cos = mathh.cos; sin = mathh.sin; pi = mathh.pi
@pyimport alive_progress; alive_bar = alive_progress.alive_bar
go = pyimport("plotly.graph_objects")


## define a Python function usable in Julia
py"""
import numpy as np
import pandas as pd
from math import cos, sin, pi
from alive_progress import alive_bar

## versions:
# numpy==1.20.1
# pandas==1.2.4
# plotly==5.9.0

def d_sphere_P(r, θ, φ, x0, y0, z0, μ0):
    dp = np.zeros((len(r), len(θ), len(φ)))
    dp_dict = {'r':[], 'θ':[], 'φ':[], 'dp':[]}
    
    with alive_bar(len(r)*len(θ)*len(φ), force_tty=True) as bar:
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

                    dp_dict['r'].append(R)
                    dp_dict['θ'].append(Ө)
                    dp_dict['φ'].append(Φ)
                    dp_dict['dp'].append(dp[i,j,k])
                    bar()
    return dp, pd.DataFrame(dp_dict)
"""

## alias
d_sphere_P_jl = py"""d_sphere_P"""


## some variables
μ0 = 20.0 # radius
λe = 2.0e-3 # μ0; # EAL
L = 3 * μ0

N = 10
K = 256
J = 70
r = np.linspace(0.0, μ0, J)
θ = np.linspace(0.0, 2*pi, J)
φ = np.linspace(0.0, pi, J)
y = 0.0


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
DD, DDdf = d_sphere_P_jl(r, θ , φ, x0_near, y0_near, z0_near, μ0)

## transform the coordinate system: spherical->cartesian
DDdf["x"] = DDdf["r"] * np.cos(DDdf["θ"]) * np.sin(DDdf["φ"])
DDdf["y"] = DDdf["r"] * np.sin(DDdf["θ"]) * np.sin(DDdf["φ"])
DDdf["z"] = DDdf["r"] * np.cos(DDdf["φ"])


## plotting
# dict_sphere = Dict("size"=>10,"color"=>DDdf["dp"],"colorscale"=>"Viridis","colorbar"=>Dict("thickness"=>20,"orientation"=>"h"),"opacity"=>1) # orientation is not implemented in the version istalled on my computer
dict_sphere = Dict("size"=>10,"color"=>DDdf["dp"],"colorscale"=>"Viridis","colorbar"=>Dict("thickness"=>20),"opacity"=>1)

fig = go.Figure(go.Scatter3d(x=DDdf["x"], y=DDdf["y"], z=DDdf["z"],
mode="markers",
name="Sphere",
marker=dict_sphere
))

fig.add_trace(
    go.Scatter3d(x=[x0_near],
                 y=[y0_near],
                 z=[z0_near],
                 name="Detector",
                 mode="markers",
                 marker=Dict(
                            "size"=>12,
                            "color"=>"red"          
    )))

fig.update_layout(width=700, height=700, scene_camera_eye=Dict("x"=>1.5, "y"=>-1.5, "z"=>1), scene_aspectmode="data")

fig.show()
