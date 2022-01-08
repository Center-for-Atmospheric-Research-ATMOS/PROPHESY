"""baseline correction from Baek 2015 (it's an application of IRLS)
Baek SJ, Park A, Ahn YJ, Choo J. Baseline correction using asymmetrically reweighted penalized least squares smoothing. Analyst. 2015;140(1):250-7.
y: raw spectra (as in still with the baseline)
λ: regularization parameter (so far, for my cases, 1.0e5 was looking good)
κ: the asymmetry parameter p is recommended to set between 0.001 and 0.1
D: the second order difference matrix (or any orther regularization matrix)
Nmax: max number of iterations"""
function baseline_removal(y::Array{Cdouble,1},λ::Cdouble,κ::Cdouble,D::Array{Cdouble,2},Nmax::Int64=50)
   H = λ*D'*D
   Ny = length(y)
   w = ones(Cdouble,Ny)
   W = Matrix{Cdouble}(I,Ny,Ny)
   z = zeros(Cdouble,Ny)
   for t in 1:Nmax
      for i in 1:Ny
         W[i,i] = w[i]
      end
      z = inv(W+H)*W*y;
      d = y-z;
      dn = d[d.<0.0]
      m = mean(dn)
      s = std(dn)
      wt = 1.0./(1.0 .+ exp.(2.0*(d.-(2.0s-m))/s))
      if ((norm(w-wt)/norm(w))<κ)
         break;
      end
      w[:] = wt
   end
   z
end

"""X: samples from the probability distribution to be fitted with M (=length(τm)) modes
τm:  the inital probability of a point of the spectra to belong to one of the different peaks
μm:  initial mean of the peaks (the modes)
σm:  initial standard deviation of the peaks (width of the modes)
Nt:  number of iterations for the EM algorithm (se Dempster 77)
Dempster AP, Laird NM, Rubin DB. Maximum likelihood from incomplete data via the EM algorithm. Journal of the Royal Statistical Society: Series B (Methodological). 1977 Sep;39(1):1-22"""
function EM_peaks(X::Array{Cdouble,1},τm::Array{Cdouble,1},μm::Array{Cdouble,1},σm::Array{Cdouble,1},Nt::Int64=10)
   n = length(X)
   M = length(τm)
   T = zeros(Cdouble,M,n);
   F = zeros(Cdouble,M,n);
   τt = copy(τm);
   μt = copy(μm);
   σt = copy(σm);
   for t in 1:Nt
      # compute T
      for i in 1:n
         for m in 1:M
            F[m,i] = (1.0/σt[m])*exp(-0.5*((X[i]-μt[m])/σt[m])^2)
         end
         T[:,i] = τt.*F[:,i]/sum(τt.*F[:,i])
      end
      # update τ
      τt[:] = (1.0/n)*sum(T,dims=2);
      # update mean and covariance
      for m in 1:M
         μt[m] = sum(T[m,:].*X)/sum(T[m,:])
         σt[m] = sqrt.(sum(T[m,:].*((X.-μt[m]).^2))/sum(T[m,:]))
      end
   end
   τt,μt,σt
end

"""Bes: binding energy (one for each point of the spectra)
PES: the spectra
τm:  the inital probability of a point of the spectra to belong to one of the different peaks
μm:  initial mean of the peaks (the modes)
σm:  initial standard deviation of the peaks (width of the modes)
Ns:  number of sample per point in the spectra"""
function EM_peaks(Bes::Array{Cdouble,1},PES::Array{Cdouble,1},τm::Array{Cdouble,1},μm::Array{Cdouble,1},σm::Array{Cdouble,1},Ns::Int64=500)
   Npes = length(PES);
   Xsample = zeros(Cdouble,Npes*Ns);
   Rsample = rand(Npes*Ns);
   cumDist = cumsum(PES)/sum(PES);
   # the spectra probably have negative values and values bigger than 1, so they must be taken care of
   idx_neg = findlast(cumDist.<0.0)
   if size(idx_neg)!=0
      cumDist[1:idx_neg].=0.0;
   end
   idx_pos = findfirst(cumDist.>1.0)
   if size(idx_neg)!=0
      cumDist[idx_pos:end].=1.0;
   end

   # draw samples
   for i in 1:length(Xsample)
      idx = findfirst(cumDist.>Rsample[i])
      Xsample[i] = Bes[idx]
   end
   EM_peaks(Xsample,τm,μm,σm,Ns)
end


"""
I_nbl:   photo-electric signal (corrected for the baseline)
Kes:     kinetic energy discretization points (regularly sub division)
σ_I:     noise level in the measurement (standard deviation)
Nbfgs:   number of iterations in the bounded BFGS loop
Nsearch: maximum number of iteration for the line search
"""
function cross_section_spread_function(I_nbl::Array{Cdouble,1},Kes::Array{Cdouble,1},σ_I::Cdouble;Nbfgs::Int64=1000,Nsearch::Int64=10)
   Ny = length(I_nbl)
   if (length(Kes)!=Ny) throw("not the right amount of kinetic energies") end
   dKe = Kes[2] - Kes[1];
   R = dKe*sum(I_nbl);
   σ_R = dKe*σ_I*sqrt(Ny);
   D_2nd = diagm(Ny-2,Ny,1 => 2ones(Cdouble,Ny-2), 0 => -ones(Cdouble,Ny-2) ,2 => -ones(Cdouble,Ny-2));

   # augmented measurement
   Y = [I_nbl; 1.0; zeros(Cdouble,Ny-2)];

   # augmented measurement operator (augmented => a priori in the operator)
   Rreg = [R*Matrix{Cdouble}(I,Ny,Ny); dKe*ones(Cdouble,Ny)'; D_2nd];

   # covariance matrix of the augmented measurements
   ΓI = zeros(Cdouble,Ny+1+(Ny-2),Ny+1+(Ny-2));
   ΓI[1:Ny,1:Ny] = σ_I^2*Matrix{Cdouble}(I,Ny,Ny);                # noise level in the measurment
   ΓI[Ny+1,Ny+1] = (0.5/(Kes[1]-Kes[end]))^2;                     # accuracy of the numerical integration
   ΓI[Ny+2:end,Ny+2:end] = 0.001^2*Matrix{Cdouble}(I,Ny-2,Ny-2);  # variance of the second order difference of the CS spread function
   ΓI_inv = inv(ΓI);
   W_inv  = inv(Rreg'*ΓI_inv*Rreg);

   # optimization problem
   function F(x::Array{Cdouble,1})
       0.5*(Y-Rreg*x)'*ΓI_inv*(Y-Rreg*x)
   end

   function Fgrad(x::Array{Cdouble,1})
       -Rreg'*ΓI_inv*(Y-Rreg*x)
   end

   X0 = zeros(Cdouble,Ny);
   p0 = zeros(Cdouble,Ny);            # first descent direction
   alpha_min = -4.0                   # smallest value of the length of the step
   alpha_max = 4.0                    # smallest value of the length of the step 2000000.0
   mu = 0.4                           # <0.5 parameter for the line search algorithm
   # H0 = Matrix{Cdouble}(I,Ny,Ny)/R1;  # initial inverse Hessian matrix
   H0 = W_inv
   # define the constraints
   lx = zeros(Cdouble,Ny);            # lower bounds
   ux = Inf*ones(Cdouble,Ny);         # upper bounds
   Xend,Hend,Xpath,Nlast = BFGSB(X0,H0,Nbfgs,alpha_min,alpha_max,mu,lx,ux,F,Fgrad,Nsearch);
   Xend,Hend,Xpath,Nlast,R,σ_R
end
