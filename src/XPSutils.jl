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
