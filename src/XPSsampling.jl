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
function acceptSample(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},p::Cdouble,y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2})
    # if the posterior probability is larger for the proposed state ρ_prop than the current state ρ_cur, then accept the state, otherwise, reject it with probability p

    # likelihood
    r_cp = 0.5*(y-H*ρ_cur)'*ΓIinv*(y-H*ρ_cur)-0.5*(y-H*ρ_prop)'*ΓIinv*(y-H*ρ_prop);
    # smoothness
    r_cp = r_cp + 0.5*(yd-D*ρ_cur)'*Γdinv*(yd-D*ρ_cur)-0.5*(yd-D*ρ_prop)'*Γdinv*(yd-D*ρ_prop); 
    
    ρ_new = ρ_cur;
    if (r_cp>=0.0)
        # unconditionally accept the new state
        ρ_new = ρ_prop
    else
        # accept the state with probability p
        # if (rand()<=p)
        if (log(rand())<=r_cp)
            ρ_new = ρ_prop
        else
            r_cp = 0.0
        end
    end
    ρ_new,r_cp # maybe it could return the computed values
end

function acceptSampleMargin(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},p::Cdouble,y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2},ΓHΓyinv::Array{Cdouble,2})
    # if the posterior probability is larger for the proposed state ρ_prop than the current state ρ_cur, then accept the state, otherwise, reject it with probability p

    # likelihood
    r_cp = 0.5*(y-H*ρ_cur)'*ΓIinv*(y-H*ρ_cur)-0.5*(y-H*ρ_prop)'*ΓIinv*(y-H*ρ_prop);
    # smoothness
    r_cp = r_cp + 0.5*(yd-D*ρ_cur)'*Γdinv*(yd-D*ρ_cur)-0.5*(yd-D*ρ_prop)'*Γdinv*(yd-D*ρ_prop); 
    # marginalization of errors
    r_cp = r_cp + 0.5*ρ_cur'*ΓHΓyinv*ρ_cur-0.5*ρ_prop'*ΓHΓyinv*ρ_prop; 
    
    ρ_new = ρ_cur;
    if (r_cp>=0.0)
        # unconditionally accept the new state
        ρ_new = ρ_prop
    else
        # accept the state with probability p
        # if (rand()<=p)
        if (log(rand())<=r_cp)
            ρ_new = ρ_prop
        else
            r_cp = 0.0
        end
    end
    ρ_new,r_cp # maybe it could return the computed values
end

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
        # if (rand()<=p)
        if (log(rand())<=r_cp)
            ρ_new = ρ_prop
        else
            r_cp = 0.0
        end
    end
    ρ_new,r_cp # maybe it could return the computed values
end

function acceptSampleSA(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},Tacc::Cdouble,y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1})
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
        if (rand()<=exp(r_cp/Tacc))
            ρ_new = ρ_prop
        end
    end
    ρ_new,r_cp # maybe it could return the computed values
end


function acceptSampleEntropy(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},p::Cdouble,y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},ρ_prior::Array{Cdouble,1},wE::Cdouble)
    # if the posterior probability is larger for the proposed state ρ_prop than the current state ρ_cur, then accept the state, otherwise, reject it with probability p

    # likelihood
    r_cp = 0.5*(y-H*ρ_cur)'*ΓIinv*(y-H*ρ_cur)-0.5*(y-H*ρ_prop)'*ΓIinv*(y-H*ρ_prop);
    # entropy 
    r_cp = r_cp + wE*(sum(ρ_cur-ρ_prop) + sum(ρ_prop.*log.(ρ_prop./ρ_prior) - ρ_cur.*log.(ρ_cur./ρ_prior)))
    
    ρ_new = ρ_cur;
    if (r_cp>=0.0)
        # unconditionally accept the new state
        ρ_new = ρ_prop
    else
        # accept the state with probability p
        # if (rand()<=p)
        if (log(rand())<=r_cp)
            ρ_new = ρ_prop
        end
    end
    ρ_new,r_cp # maybe it could return the computed values
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
        # if (rand()<=p)
        if (log(rand())<=r_cp)
            ρ_new = ρ_prop
        else
            r_cp = 0.0
        end
    end
    ρ_new,r_cp # maybe it could return the computed values
end


"""
    Sampling the posterior distribution 
    P(ρ|Y,H) ∝ P(Y|ρ,H)P(ρ)
    retunr an array of samples ρ
"""
function samplePosterior(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},p0::Array{Cdouble,1},y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2};Ns::Int64=10000)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    ρ_all[1,:] = ρ_start;
    deltaU = zeros(Cdouble,Ns);
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanismSmooth(ρ_all[i,:],Γsqrt)
        ρ_all[i+1,ρ_all[i+1,:].<0.0] .= 0.0  
        
        # accept or reject the sample
        ρ_all[i+1,:],deltaU[i] = acceptSample(ρ_all[i,:],ρ_all[i+1,:],p0[i],y,yd,ΓIinv,Γdinv,H,D)
    end
    ρ_all,deltaU
end

function samplePosteriorMargin(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},p0::Array{Cdouble,1},y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2},ΓHΓyinv::Array{Cdouble,2};Ns::Int64=10000)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    ρ_all[1,:] = ρ_start;
    deltaU = zeros(Cdouble,Ns);
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanismSmooth(ρ_all[i,:],Γsqrt)
        ρ_all[i+1,ρ_all[i+1,:].<0.0] .= 0.0  
        
        # accept or reject the sample
        ρ_all[i+1,:],deltaU[i] = acceptSampleMargin(ρ_all[i,:],ρ_all[i+1,:],p0[i],y,yd,ΓIinv,Γdinv,H,D,ΓHΓyinv)
    end
    ρ_all,deltaU
end

function samplePosterior(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},p0::Cdouble,y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1};Ns::Int64=10000,psmooth::Cdouble=0.99)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    ρ_all[1,:] = ρ_start;
    deltaU = zeros(Cdouble,Ns);
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanism(ρ_all[i,:],Γsqrt,σB;psmooth=psmooth)
        
        # accept or reject the sample
        p = p0*(Ns-i)/(Ns-1.0) 
        ρ_all[i+1,:],deltaU[i] = acceptSample(ρ_all[i,:],ρ_all[i+1,:],p,y,ΓIinv,H,Dprior,ρB,σB)
    end
    ρ_all,deltaU
end


function samplePosterior(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},p0::Array{Cdouble,1},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1};Ns::Int64=10000,psmooth::Cdouble=0.99)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    ρ_all[1,:] = ρ_start;
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanism(ρ_all[i,:],Γsqrt,σB;psmooth=psmooth)
        
        # accept or reject the sample
        ρ_all[i+1,:],_ = acceptSample(ρ_all[i,:],ρ_all[i+1,:],p0[i],y,ΓIinv,H,Dprior,ρB,σB)
    end
    ρ_all
end


function samplePosteriorEntropy(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},p0::Array{Cdouble,1},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},ρ_prior::Array{Cdouble,1},wE::Cdouble;Ns::Int64=10000)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    ρ_all[1,:] = ρ_start;
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanismSmooth(ρ_all[i,:],Γsqrt)
        ρ_all[i+1,ρ_all[i+1,:].<=0.0] .= 1.0e-6;
        
        # accept or reject the sample
        ρ_all[i+1,:],_ = acceptSampleEntropy(ρ_all[i,:],ρ_all[i+1,:],p0[i],y,ΓIinv,H,ρ_prior,wE)
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
        ρ_all[i+1,:],_ = acceptSampleModelMargin(ρ_all[i,:],ρ_all[i+1,:],p0[i],y,ΓIinv,H,ΓH,Dprior,ρB,σB)
    end
    ρ_all
end



##
## why not SA
##


function samplePosteriorDeltaU(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},p0::Array{Cdouble,1},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1};Ns::Int64=10000,psmooth::Cdouble=0.99)
    # all samples
    ρ_all  = zeros(Cdouble,Ns+1,length(ρ_start))
    deltaU = zeros(Cdouble,Ns);
    ρ_all[1,:] = ρ_start;
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanism(ρ_all[i,:],Γsqrt,σB;psmooth=psmooth)
        
        # accept or reject the sample
        ρ_all[i+1,:],deltaU[i] = acceptSample(ρ_all[i,:],ρ_all[i+1,:],p0[i],y,ΓIinv,H,Dprior,ρB,σB)
    end
    ρ_all,deltaU
end

function samplePosteriorEntropyDeltaU(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},p0::Array{Cdouble,1},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},ρ_prior::Array{Cdouble,1},wE::Cdouble;Ns::Int64=10000)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    deltaU = zeros(Cdouble,Ns);
    ρ_all[1,:] = ρ_start;
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanismSmooth(ρ_all[i,:],Γsqrt)
        ρ_all[i+1,ρ_all[i+1,:].<=0.0] .= 1.0e-6;
        
        # accept or reject the sample
        ρ_all[i+1,:],deltaU[i] = acceptSampleEntropy(ρ_all[i,:],ρ_all[i+1,:],p0[i],y,ΓIinv,H,ρ_prior,wE)
    end
    ρ_all,deltaU
end

function samplePosteriorModelMarginDeltaU(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},p0::Array{Cdouble,1},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},ΓH::Array{Cdouble,3},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1};Ns::Int64=10000,psmooth::Cdouble=0.99)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    deltaU = zeros(Cdouble,Ns);
    ρ_all[1,:] = ρ_start;
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanism(ρ_all[i,:],Γsqrt,σB;psmooth=psmooth)
        
        # accept or reject the sample
        ρ_all[i+1,:],deltaU[i] = acceptSampleModelMargin(ρ_all[i,:],ρ_all[i+1,:],p0[i],y,ΓIinv,H,ΓH,Dprior,ρB,σB)
    end
    ρ_all,deltaU
end


function dichotomyTemperature(τ::Cdouble,deltaU::Array{Cdouble,1},Niter::Int64=100)
    Tmin = -minimum(deltaU)/log(τ);    # minimum value for the set
    Tmax = -mean(deltaU)/log(τ);       # maximum value for the set: due to concavity of exp
    f_exp = (x->mean(exp.(-deltaU./x))-τ);
    fval = zeros(Cdouble,Niter+1);
    for i in 1:Niter
        Tmean = 0.5*(Tmin+Tmax);
        fval[i] = f_exp(Tmean);
        if (fval[i]<=0)
            Tmin = Tmean;
        else
            Tmax = Tmean;
        end
    end
    fval[end] = f_exp(0.5*(Tmin+Tmax))
    0.5*(Tmin+Tmax),fval
end


function posteriorTemperature(τ::Cdouble,ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},p0::Array{Cdouble,1},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1};Ns::Int64=10000,psmooth::Cdouble=0.99)
    # generate samples 
    ρ_all,deltaU = samplePosteriorDeltaU(ρ_start,Γsqrt,p0,y,ΓIinv,H,Dprior,ρB,σB;Ns=Ns,psmooth=psmooth)
    # compute temperature (dichotomy) (keep only positive variations of energy)
    dichotomyTemperature(τ,deltaU[deltaU.>0.0],100)
end

function posteriorTemperatureEntropy(τ::Cdouble,ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},p0::Array{Cdouble,1},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},ρ_prior::Array{Cdouble,1},wE::Cdouble;Ns::Int64=10000)
    # generate samples 
    ρ_all,deltaU = samplePosteriorEntropyDeltaU(ρ_start,Γsqrt,p0,y,ΓIinv,H,ρ_prior,wE;Ns=Ns)
    # compute temperature (dichotomy)
    dichotomyTemperature(τ,deltaU[deltaU.>0.0],100)
end

function posteriorTemperatureModelMargin(τ::Cdouble,ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},p0::Array{Cdouble,1},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},ΓH::Array{Cdouble,3},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1};Ns::Int64=10000,psmooth::Cdouble=0.99)
    # generate samples 
    ρ_all,deltaU = samplePosteriorModelMarginDeltaU(ρ_start,Γsqrt,p0,y,ΓIinv,H,ΓH,Dprior,ρB,σB;Ns=Ns,psmooth=psmooth)
    # compute temperature (dichotomy)
    dichotomyTemperature(τ,deltaU[deltaU.>0.0],100)
end


function posteriorSA(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},Tacc::Array{Cdouble,1},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1};Ns::Int64=10000,psmooth::Cdouble=0.99)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    r_cp = zeros(Cdouble,Ns)
    ρ_all[1,:] = ρ_start;
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanism(ρ_all[i,:],Γsqrt,σB;psmooth=psmooth)
        
        # accept or reject the sample
        ρ_all[i+1,:],r_cp[i] = acceptSampleSA(ρ_all[i,:],ρ_all[i+1,:],Tacc[i],y,ΓIinv,H,Dprior,ρB,σB)
    end
    ρ_all,r_cp
end