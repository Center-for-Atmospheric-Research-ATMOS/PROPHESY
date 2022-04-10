# Sampling the a posteriori to estimate the uncertainty
#  - P(ρ|I,H): the probability density of the state ρ conditionally to the data I and the model H
#  - P(ρ|I):   the marginalization of P(ρ|I,H) over the measurement model space (using the small perturbation assumption)


"""
    Likelihood of the measurement knowing the state and the model

    likelihood_H(I::Array{Cdouble,1},x::Array{Cdouble,1},H::Array{Cdouble,2},ΓIinv::Array{Cdouble,2},detΓI::Cdouble)

    I: array of data
    x: state of the system
    H: measurement operator
    ΓIinv: inverse of the measurement covariance matrix
    detΓI: determinant of the measurement covariance matrix
"""
function likelihood_H(I::Array{Cdouble,1},x::Array{Cdouble,1},H::Array{Cdouble,2},ΓIinv::Array{Cdouble,2},detΓI::Cdouble) # could become Poisson or something else
    (1.0/sqrt(2π*detΓI))*exp(-0.5*(I-H*x)'*ΓIinv*(I-H*x))
end

function corrCovariance(w::Array{Cdouble,1};cor_len::Cdouble=5.0)
    Nr = length(w);

    Γprior = zeros(Cdouble,Nr,Nr);

    for i in 1:Nr
        Γprior[i,i] = w[i]^2;
        for j in i+1:Nr
            Γprior[i,j] = Γprior[i,i]*exp(-(i-j)^2/(0.5*cor_len^2))
            Γprior[j,i] = Γprior[i,j]
        end
    end

    Γprior
end

function smoothnessCovariance(w::Array{Cdouble,1};cor_len::Cdouble=5.0)
    Nr = length(w);

    Γprior = zeros(Cdouble,Nr,Nr);
    Dprior = D2nd(Nr+2)[:,2:end-1];

    for i in 1:Nr
        Γprior[i,i] = w[i]^2;
        for j in i+1:Nr
            Γprior[i,j] = Γprior[i,i]*exp(-(i-j)^2/(0.5*cor_len^2))
            Γprior[j,i] = Γprior[i,j]
        end
    end

    Bprior = Dprior'*inv(Γprior)*Dprior;
    Dsqrt = real(sqrt(inv(Bprior)));
    Dsqrt, Γprior, Dprior
end

# the transition mechanism generates states that samples the proposal density, which should be fairly similar to the density we hope to sample.
# In our case, the a priori distribution in the a posteriori is fairly easy to sample and so I choose to build the transmission mechanism upon it.
function transmissionMechanismSmooth(x_curr::Array{Cdouble,1},Γsqrt::Array{Cdouble,2})
    x_curr + Γsqrt*randn(Cdouble,length(x_curr))
end

function transmissionMechanismBoundaries(x_curr::Array{Cdouble,1},σB::Array{Cdouble,1})
    x_prop = x_curr;
    if (rand(Cdouble)<0.5)
        x_prop[1] = x_prop[1] + σB[1]*randn(Cdouble);
    else
        x_prop[end] = x_prop[end] + σB[2]*randn(Cdouble)
    end
    x_prop
end

"""
    the transmission mechanism samples the proposal density

    The density is taken as a truncated Gaussian distribution centered on the current state x_curr
    and whose covariance matrix is Γsqrt^2. The truncation comes in to play when negative values 
    appear in the generated state
"""
function transmissionMechanism(x_curr::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},σB::Array{Cdouble,1};psmooth::Cdouble=0.99)
    if (rand(Cdouble)<psmooth)
        x_prop = transmissionMechanismSmooth(x_curr,Γsqrt)
    else
        x_prop = transmissionMechanismBoundaries(x_curr,σB)
    end
    x_prop[x_prop.<0.0] .= 0.0
    x_prop
end


#TODO: optimization SA. A communication kernel that seems like a good idea is one that would start changing only the values closest to the surface and gradually changes deeper and deeper values leaving out the surface


"""
    accepting mechanism based on Metropolis-Hasting algorithm (for symmetric communication mechanisms)
"""
function acceptSample(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},p::Cdouble,y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1})
    # if the posterior probability is larger for the proposed state ρ_prop than the current state ρ_cur, then accept the state, otherwise, reject it with probability p

    # likelihood
    r_cp = 0.5*(y-H*ρ_cur)'*ΓIinv*(y-H*ρ_cur)-0.5*(y-H*ρ_prop)'*ΓIinv*(y-H*ρ_prop);
    # smoothness
    r_cp = r_cp + 0.5ρ_cur'Dprior*ρ_cur - 0.5ρ_prop'Dprior*ρ_prop;
    # known values 
    r_cp = r_cp + 0.5*((ρ_cur[1]  -ρB[1])^2)/(σB[1]^2) - 0.5*((ρ_prop[1]  -ρB[1])^2)/(σB[1]^2);
    r_cp = r_cp + 0.5*((ρ_cur[end]-ρB[2])^2)/(σB[2]^2) - 0.5*((ρ_prop[end]-ρB[2])^2)/(σB[2]^2);
    
    ρ_new = ρ_cur;
    if (r_cp>=0.0)
        # unconditionally accept the new state
        ρ_new = ρ_prop
    else
        # accept the state with probability p
        if (rand()<=p)
            ρ_new = ρ_prop
        end
    end
    ρ_new # maybe it could return the computed values
end


function acceptSampleModelMargin(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},p::Cdouble,y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},ΓH::Array{Cdouble,3},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1})
    # if the posterior probability is larger for the proposed state ρ_prop than the current state ρ_cur, then accept the state, otherwise, reject it with probability p

    # likelihood
    r_cp = 0.5*(y-H*ρ_cur)'*ΓIinv*(y-H*ρ_cur)-0.5*(y-H*ρ_prop)'*ΓIinv*(y-H*ρ_prop);
    # smoothness
    r_cp = r_cp + 0.5ρ_cur'Dprior*ρ_cur - 0.5ρ_prop'Dprior*ρ_prop;
    # known values 
    r_cp = r_cp + 0.5*((ρ_cur[1]  -ρB[1])^2)/(σB[1]^2) - 0.5*((ρ_prop[1]  -ρB[1])^2)/(σB[1]^2);
    r_cp = r_cp + 0.5*((ρ_cur[end]-ρB[2])^2)/(σB[2]^2) - 0.5*((ρ_prop[end]-ρB[2])^2)/(σB[2]^2);
    # marginalization of the model
    for i in 1:length(y)
        r_cp = r_cp + 0.5ΓIinv[i,i]*(ρ_cur'*ΓH[i,:,:]*ρ_cur - ρ_prop'*ΓH[i,:,:]*ρ_prop)
    end

    ρ_new = ρ_cur;
    if (r_cp>=0.0)
        # unconditionally accept the new state
        ρ_new = ρ_prop
    else
        # accept the state with probability p
        if (rand()<=p)
            ρ_new = ρ_prop
        end
    end
    ρ_new # maybe it could return the computed values
end


"""
    Sampling the posterior distribution 
    P(ρ|Y,H) ∝ P(Y|ρ,H)P(ρ)
    retunr an array of samples ρ
"""
function samplePosterior(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},p0::Cdouble,y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1};Ns::Int64=10000,psmooth::Cdouble=0.99)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    ρ_all[1,:] = ρ_start;
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanism(ρ_all[i,:],Γsqrt,σB;psmooth=psmooth)
        
        # accept or reject the sample
        p = p0*(Ns-i)/(Ns-1.0) 
        ρ_all[i+1,:] = acceptSample(ρ_all[i,:],ρ_all[i+1,:],p,y,ΓIinv,H,Dprior,ρB,σB)
    end
    ρ_all
end

function samplePosterior(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},p0::Array{Cdouble,1},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1};Ns::Int64=10000,psmooth::Cdouble=0.99)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    ρ_all[1,:] = ρ_start;
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanism(ρ_all[i,:],Γsqrt,σB;psmooth=psmooth)
        
        # accept or reject the sample
        ρ_all[i+1,:] = acceptSample(ρ_all[i,:],ρ_all[i+1,:],p0[i],y,ΓIinv,H,Dprior,ρB,σB)
    end
    ρ_all
end


function samplePosteriorModelMargin(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},p0::Array{Cdouble,1},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},ΓH::Array{Cdouble,3},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1};Ns::Int64=10000,psmooth::Cdouble=0.99)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    ρ_all[1,:] = ρ_start;
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanism(ρ_all[i,:],Γsqrt,σB;psmooth=psmooth)
        
        # accept or reject the sample
        ρ_all[i+1,:] = acceptSampleModelMargin(ρ_all[i,:],ρ_all[i+1,:],p0[i],y,ΓIinv,H,ΓH,Dprior,ρB,σB)
    end
    ρ_all
end