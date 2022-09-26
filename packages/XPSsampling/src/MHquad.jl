function acceptSample(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2})
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
        # accept the state with probability e^{r_cp}
        if (log(rand())<=r_cp)
            ρ_new = ρ_prop
        else
            r_cp = 0.0
        end
    end
    ρ_new,r_cp # maybe it could return the computed values
end

function samplePosterior(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2};Ns::Int64=10000)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    ρ_all[1,:] = ρ_start;
    deltaU = zeros(Cdouble,Ns);
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanismSmooth(ρ_all[i,:],Γsqrt)
        ρ_all[i+1,ρ_all[i+1,:].<0.0] .= 0.0  
        
        # accept or reject the sample
        ρ_all[i+1,:],deltaU[i] = acceptSample(ρ_all[i,:],ρ_all[i+1,:],y,yd,ΓIinv,Γdinv,H,D)
    end
    ρ_all,deltaU
end

function samplePosteriorMeanAndCov(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2};Ns::Int64=10000,Nburn::Int64=1000)
    # all samples
    ρ_curr = copy(ρ_start);
    ρ_prop = copy(ρ_start);
    μρ_cum = zeros(Cdouble,length(ρ_start));
    Γρ_cum = zeros(Cdouble,length(ρ_start),length(ρ_start));
    deltaU = zeros(Cdouble,Ns);
    for i in 1:Ns
        # global ρ_prop
        # global deltaU
        # draw a new sample from a distribution not to far from the actual one
        ρ_prop = transmissionMechanismSmooth(ρ_curr,Γsqrt)
        ρ_prop[ρ_prop.<0.0] .= 0.0  
        
        # accept or reject the sample
        ρ_prop,deltaU[i] = acceptSample(ρ_curr,ρ_prop,y,yd,ΓIinv,Γdinv,H,D)
        if (i>=Nburn)
            μρ_cum = μρ_cum + ρ_prop; # global 
            Γρ_cum = Γρ_cum + ρ_prop*ρ_prop'; #global 
        end
        ρ_curr = ρ_prop;
    end
    μρ_cum = μρ_cum/(Ns-Nburn+1);
    Γρ_cum = (1.0/(Ns-Nburn))*(Γρ_cum - (Ns-Nburn+1)*μρ_cum*μρ_cum');
    μρ_cum,Γρ_cum,deltaU
end


##
## adding marginalization of model error
##
function acceptSampleMargin(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2},ΓHΓyinv::Array{Cdouble,2})
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
        # accept the state with probability e^{r_cp}
        if (log(rand())<=r_cp)
            ρ_new = ρ_prop
        else
            r_cp = 0.0
        end
    end
    ρ_new,r_cp # maybe it could return the computed values
end

function samplePosteriorMargin(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2},ΓHΓyinv::Array{Cdouble,2};Ns::Int64=10000)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    ρ_all[1,:] = ρ_start;
    deltaU = zeros(Cdouble,Ns);
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanismSmooth(ρ_all[i,:],Γsqrt)
        ρ_all[i+1,ρ_all[i+1,:].<0.0] .= 0.0  
        
        # accept or reject the sample
        ρ_all[i+1,:],deltaU[i] = acceptSampleMargin(ρ_all[i,:],ρ_all[i+1,:],y,yd,ΓIinv,Γdinv,H,D,ΓHΓyinv)
    end
    ρ_all,deltaU
end

