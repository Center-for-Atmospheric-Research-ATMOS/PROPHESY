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
Nθ = 256;
Ny = 256;
μ0 = 20.0 # 100.0;
L = 200.0;
x0 = sqrt(2.0)*100.0
y0 = 0.0;
z0 = 100.0
r = collect(range(μ0-k0*λe0,μ0,length=Nr))
θ0 = atan(x0,z0)
θ = collect(range(θ0-π/2.0,θ0+π/2.0,Nθ));
y = collect(range(-L/2.0,L/2.0,length=Ny));

wsGeom = cylinderGeom(x0,y0,z0,μ0,r,θ,y);



# simulate some data (one point in the kinetic energy spectrum for four different concentration profiles)
# ρA_1 = logistic.(1000.0reverse(μ0.-r).-2.0,0.0,1.0,2.0);
# ρA_1 = logistic.(1000.0reverse(μ0.-r).-2.0,0.0,1.0,2.0) .+ 2.0exp.(-(1000.0reverse(μ0.-r).-1.0).^2. /(2.0*0.25^2));
ρA_1 = logistic.(1000.0reverse(μ0.-r).-2.0,0.0,1.0,2.0) .+ exp.(-(1000.0reverse(μ0.-r).-1.5).^2. /(2.0*0.5^2));
# ρA_1 = exp.(-(1000.0reverse(μ0.-r).-2.5).^2. /(2.0*0.5^2));


# measurement operator (only the geometical term since the other comes as a multiplicative scalar estimated from the data)
Ndata = 6 # 25
H_better = zeros(Cdouble,Ndata,Nr);
# λbetter0  = 1.0e-3*[1.0; 1.5; 2.0; 2.5; 3.0]; # these are some eal values that would nice to be able to access... but that would not be sufficient to make the uncertainty small enough
λbetter0 = 1.0e-3collect(range(1.3,2.5,Ndata));

for i in 1:Ndata
    H_better[i,:],_,_,_,_ = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,λbetter0[i]);
end

H_better = reverse(H_better,dims=2); #


# deeper than some distance, the signal is not likely to be disantangled
d0 = 15.0e-3 # NOTE: this value should depend on the penetration depth
N0 = findfirst(r.-μ0.>=-d0) 
figure(); plot(r.-μ0,H_better'); plot(r.-μ0,ρA_1)
# plot(r[end-(29+50+70):end].-μ0,ρA_1[end-(29+50+70):end])
plot(r[N0+1:end].-μ0,ρA_1[N0+1:end])


# generate some data (data point and covariance)
Nnoise = 10;
σnoise = 0.1*ones(Cdouble,Nnoise);

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