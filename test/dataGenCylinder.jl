## acquisition setup
ħν = 900.0;
μKe = ħν-285.0;
α = 1.0
T = 1.0
Fν = 1.0;
Nke = 200;
Ke = collect(range(μKe-2.0,μKe+2.0,length=Nke));
dKe = Ke[2]-Ke[1]
σν0 = 0.6;
σν = σν0*((0.7/sqrt(2π*0.2^2))*exp.(-0.5*(Ke.-(μKe-0.5)).^2/0.2^2) .+ (0.3/sqrt(2π*0.5^2))*exp.(-0.5*(Ke.-(μKe+0.5)).^2/0.5^2));
λe0 = 0.002;

wsAcq = XPSacq(ħν,μKe,α,T,Fν,Ke,σν,λe0);

## geometry setup
k0 = 10; # 5
Nr = 201; # 101  # 51;
Nr_lowres = 101 # 21; # 51 # 101 # 201 # 401
Nθ = 256;
Ny = 256;
μ0 = 20.0 # 100.0;
L = 200.0;
x0 = sqrt(2.0)*100.0
y0 = 0.0;
z0 = 100.0
r = collect(range(μ0-k0*λe0,μ0,length=Nr))
r_lowres = collect(range(μ0-k0*λe0,μ0,length=Nr_lowres));
θ0 = atan(x0,z0)
θ = collect(range(θ0-π/2.0,θ0+π/2.0,Nθ));
y = collect(range(-L/2.0,L/2.0,length=Ny));

wsGeom = cylinderGeom(x0,y0,z0,μ0,r,θ,y);



# simulate some data (one point in the kinetic energy spectrum for four different concentration profiles)
# ρA_1 = logistic.(1000.0reverse(μ0.-r).-2.0,0.0,1.0,2.0);
# ρA_1 = logistic.(1000.0reverse(μ0.-r).-2.0,0.0,1.0,2.0) .+ 2.0exp.(-(1000.0reverse(μ0.-r).-1.0).^2. /(2.0*0.25^2));
# ρA_1 = logistic.(1000.0reverse(μ0.-r).-2.0,0.0,1.0,2.0) .+ exp.(-(1000.0reverse(μ0.-r).-1.5).^2. /(2.0*0.5^2));
ρA_1 = exp.(-(1000.0reverse(μ0.-r).-2.5).^2. /(2.0*0.5^2));


# measurement operator (only the geometical term since the other comes as a multiplicative scalar estimated from the data)
Ndata = 5 # 5 # 10 # 20 # 50 # 6 # 25
H_better = zeros(Cdouble,Ndata,Nr);
H_lowres = zeros(Cdouble,Ndata,Nr_lowres);
# λbetter0  = 1.0e-3*[1.0; 1.5; 2.0; 2.5; 3.0]; # these are some eal values that would nice to be able to access... but that would not be sufficient to make the uncertainty small enough
λbetter0 = 1.0e-3collect(range(1.3,2.5,Ndata));
# λbetter0 = 1.0e-3collect(range(0.5,5.0,Ndata));
Nλ = 21;
ΓH = Array{Cdouble}(undef,Nr_lowres,Nr_lowres,Ndata);
μH = Array{Cdouble}(undef,Nr_lowres,Ndata);
# λbetter0 = 1.0e-3collect(range(0.6,4.3,Ndata)); #NOTE: if one wants to reconstruct a structure that spans the depth from z1 to z2, then one should use penetration depth in the same range (smaller does not help at all, and bigger helps a bit, but the effect is not dramatic)
# λbetter0 = 1.0e-3collect(range(0.1,7.5,Ndata));


for i in 1:Ndata
    H_better[i,:],_,_,_,_ = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,λbetter0[i]);
    H_lowres[i,:],_,_,_,_ = cylinder_gain_H(r_lowres,θ,y,x0,y0,z0,μ0,0.9λbetter0[i]); # 0.999999999
    # local λrange = collect(range(0.95λbetter0[i],1.05λbetter0[i],length=Nλ))
    # local λrange = collect(range(0.98λbetter0[i],1.02λbetter0[i],length=Nλ))
    local λrange = collect(range(0.99λbetter0[i],1.01λbetter0[i],length=Nλ))
    local Pλ = (1.0/((λrange[2]-λrange[1])*Nλ))*ones(Cdouble,Nλ)
    ΓH[:,:,i], μH[:,i]   = cov_H_cylinder(r_lowres, θ,y,x0,y0,z0,μ0,λrange,Pλ)
end

H_better = reverse(H_better,dims=2);
H_lowres = reverse(H_lowres,dims=2);
ΓH = reverse(ΓH,dims=(1,2));
μH = reverse(μH,dims=1);

# deeper than some distance, the signal is not likely to be disantangled
d0 = 15.0e-3 # NOTE: this value should depend on the penetration depth
N0 = findfirst(r.-μ0.>=-d0) 
N0_lowres = findfirst(r_lowres.-μ0.>=-d0) 
figure(); plot(r.-μ0,H_better'); plot(r.-μ0,ρA_1)
# plot(r[end-(29+50+70):end].-μ0,ρA_1[end-(29+50+70):end])
plot(r[N0+1:end].-μ0,ρA_1[N0+1:end])


# generate some data (data point and covariance)
Nnoise = 20 #0;
σnoise = 0.001*ones(Cdouble,Nnoise);
σnoise = 0.05*ones(Cdouble,Nnoise);

y_data = zeros(Cdouble,Nnoise,Ndata);
ΓI = zeros(Cdouble,Ndata,Ndata,Nnoise);
ΓIsqrt = zeros(Cdouble,Ndata,Ndata,Nnoise);
# ΓIinv = zeros(Cdouble,Ndata,Ndata,Nnoise);
for i in 1:Nnoise
    ΓI[:,:,i] = σnoise[i]^2*diagm(ones(Cdouble,Ndata)); # iid data noise
    ΓIsqrt[:,:,i] = sqrt(ΓI[:,:,i]);
    # ΓIinv[:,:,i] = inv(ΓI[:,:,i]);
    y_data[i,:] = H_better*ρA_1 + ΓIsqrt[:,:,i]*randn(Cdouble,Ndata);
end
y_data[y_data.<0.0] = -y_data[y_data.<0.0];