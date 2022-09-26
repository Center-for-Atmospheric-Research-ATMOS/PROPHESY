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