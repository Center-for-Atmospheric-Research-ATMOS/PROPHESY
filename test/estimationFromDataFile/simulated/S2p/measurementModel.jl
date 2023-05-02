# non-geomtrical parameters # there will be a units problem, I'm sure of it
α_all     = zeros(Cdouble,NdataS2p);  # [μm^{-2}] alignment parameter
T_all     = zeros(Cdouble,NdataS2p);  # [dimensionless] transmission function
F_all     = zeros(Cdouble,NdataS2p);  # [# s^{-1}] photon flux photon per seconds
σ_all     = zeros(Cdouble,NdataS2p);  # [Mbarn] total cross section
t_all     = zeros(Cdouble,NdataS2p);  # [s]counting time in seconds
ρB_al     = zeros(Cdouble,NdataS2p);  # [M] bulk molar concentration (it's only one value, but I get it for each dataset)
# D_dark    = 0.0                       # dark current coefficient is a number of electron per unit time [# s^{-1}] 

# for each spectrum, get the alignment parameter, the transmission factor, the photon flux, the total cross section, the integration time and the bulk concentration
for i in 1:NdataS2p
    local symData = Symbol(string("hν_",round(Int64,EssS2p[2*(i-1)+1,1])));
    α_all[i] = dictAllGeomS2p[symbolDictS2p[symData]].α[1]
    T_all[i] = dictAllDataS2p[symData].T[1]
    F_all[i] = dictAllDataS2p[symData].F[1]
    σ_all[i] = dictAllDataS2p[symData].σ_tot[1]
    t_all[i] = dictAllDataS2p[symData].Δt[1]
    ρB_al[i] = dictAllGeomS2p[symbolDictS2p[symData]].ρ[1]
    if FLAG_0004
        ρB_al[i] = maximum(dictAllGeomS2p[symbolDictS2p[symData]].ρ)
    end
end

# geometry factor
x0        = dictAllGeomS2p[symbolDictS2p[Symbol(string("hν_",round(Int64,EssS2p[1,1])))]].x0[1];
y0        = dictAllGeomS2p[symbolDictS2p[Symbol(string("hν_",round(Int64,EssS2p[1,1])))]].y0[1];
z0        = dictAllGeomS2p[symbolDictS2p[Symbol(string("hν_",round(Int64,EssS2p[1,1])))]].z0[1];
μ0        = dictAllGeomS2p[symbolDictS2p[Symbol(string("hν_",round(Int64,EssS2p[1,1])))]].radius[1]; # radius of the cylinder in μm
if BETTER_MODEL # if true: extend the model outside of the sharp edge boundary to include the low density liquid at the edge of the sharp volume
    δr = 0.002 # [μm] 
else
    δr = 0.0
end
Nθ        = 256;            # number of discretization points in the polar angle dimension
Ny        = 257;            # number of discretization points in the cylinder axis dimension
Nr        = 101;             # number of discretization points in the radial dimension

L         = 50.0            # height of the irradiated sample (the vertical extent of the beam is more like 20μm instead of 100μm)
θ0        = atan(x0,z0)     # angular direction of the analyzer
max_depth = 0.02            # [μm] maximum depth of the profile # k0*λe0 # k0 = 10; # λe0 = 0.002;

r = collect(range(μ0-max_depth,μ0+δr,length=Nr));
θ = collect(range(θ0-π/2.0,θ0+π/2.0,Nθ));
y = collect(range(-L/2.0,L/2.0,length=Ny));

Hgeom = zeros(Cdouble,NdataS2p,Nr);
for i in 1:NdataS2p
    Hgeom[i,:],_,_,_,_ = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,1.0e-3λ[i]);
end
