# with boundary conditions
function acceptSample(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1})
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
        # accept the state with probability e^{r_cp}
        if (log(rand())<=r_cp)
            ρ_new = ρ_prop
        else
            r_cp = 0.0
        end
    end
    ρ_new,r_cp # maybe it could return the computed values
end





# with boundary conditions
function acceptSampleModelMargin(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},ΓH::Array{Cdouble,3},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1})
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
        # accept the state with probability e^{r_cp}
        if (log(rand())<=r_cp)
            ρ_new = ρ_prop
        else
            r_cp = 0.0
        end
    end
    ρ_new,r_cp # maybe it could return the computed values
end







function samplePosterior(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1};Ns::Int64=10000,psmooth::Cdouble=0.99)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    ρ_all[1,:] = ρ_start;
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanism(ρ_all[i,:],Γsqrt,σB;psmooth=psmooth)
        
        # accept or reject the sample
        ρ_all[i+1,:],_ = acceptSample(ρ_all[i,:],ρ_all[i+1,:],y,ΓIinv,H,Dprior,ρB,σB)
    end
    ρ_all
end




function samplePosteriorModelMargin(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},ΓH::Array{Cdouble,3},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1};Ns::Int64=10000,psmooth::Cdouble=0.99)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    ρ_all[1,:] = ρ_start;
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanism(ρ_all[i,:],Γsqrt,σB;psmooth=psmooth)
        
        # accept or reject the sample
        ρ_all[i+1,:],_ = acceptSampleModelMargin(ρ_all[i,:],ρ_all[i+1,:],y,ΓIinv,H,ΓH,Dprior,ρB,σB)
    end
    ρ_all
end
