# non-geomtrical parameters # there will be a units problem, I'm sure of it
α_all     = zeros(Cdouble,NdataC1s);  # [μm^{-2}] alignment parameter
T_all     = zeros(Cdouble,NdataC1s);  # [dimensionless] transmission function
F_all     = zeros(Cdouble,NdataC1s);  # [# s^{-1}] photon flux photon per seconds
σ_all     = zeros(Cdouble,NdataC1s);  # [Mbarn] total cross section
t_all     = zeros(Cdouble,NdataC1s);  # [s]counting time in seconds
ρB_al     = zeros(Cdouble,NdataC1s);  # [M] bulk molar concentration (it's only one value, but I get it for each dataset)
# D_dark    = 0.0                       # dark current coefficient is a number of electron per unit time [# s^{-1}] 

# for each spectrum, get the alignment parameter, the transmission factor, the photon flux, the total cross section, the integration time and the bulk concentration
for i in 1:NdataC1s
    local symData = Symbol(string("hν_",round(Int64,EssC1s[3*(i-1)+1,1])));
    α_all[i] = dictAllGeomC1s[symbolDictC1s[symData]].α[1]
    T_all[i] = dictAllDataC1s[symData].T[1]
    F_all[i] = dictAllDataC1s[symData].F[1]
    σ_all[i] = dictAllDataC1s[symData].σ_tot[1]
    t_all[i] = dictAllDataC1s[symData].Δt[1]
    ρB_al[i] = dictAllGeomC1s[symbolDictC1s[symData]].ρ[1]
    if FLAG_0004
        ρB_al[i] = maximum(dictAllGeomC1s[symbolDictC1s[symData]].ρ)
    end
end

# geometry factor
x0        = dictAllGeomC1s[symbolDictC1s[:hν_958]].x0[1];
y0        = dictAllGeomC1s[symbolDictC1s[:hν_958]].y0[1];
z0        = dictAllGeomC1s[symbolDictC1s[:hν_958]].z0[1];
μ0        = dictAllGeomC1s[symbolDictC1s[:hν_958]].radius[1]; # radius of the cylinder in μm
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

Hgeom = zeros(Cdouble,NdataC1s,Nr);
for i in 1:NdataC1s
    Hgeom[i,:],_,_,_,_ = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,1.0e-3λ[i]);
end
