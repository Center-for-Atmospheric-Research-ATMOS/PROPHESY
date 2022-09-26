function model_un(H_std::Array{Cdouble,2},σ_noise::Array{Cdouble,1})
    Nke,M = size(H_std);
    W_inv = zeros(Cdouble,M);
    for p in 1:Nke
        W_inv[:] = W_inv[:] + (1.0/σ_noise[p]^2)*H_std[p,:].^2;
    end
    W_inv
end

function low_rank_reconstruction_un_2(rankK::Int64,H::Array{Cdouble,2},H_std::Array{Cdouble,2},σ_noise::Array{Cdouble,1},Γ_inv::Array{Cdouble,2},D::Array{Cdouble,2},Γ_D_inv::Array{Cdouble,2},B::Array{Cdouble,2},Γ_B_inv::Array{Cdouble,2},y::Array{Cdouble,1},ρ_val::Array{Cdouble,1}) #
    # singular value decomposition
    F = svd(H, full=true);

    # low rank basis
    Vimg = F.Vt[1:rankK,:]'; # the image of the operator (in the right or state space)

    # low rank operators
    H_tilde = H*Vimg;
    D_tilde = D*Vimg;
    B_tilde = B*Vimg;

    #
    if false
        Nke,M = size(H);
        W_inv = zeros(Cdouble,M,M);
        for p in 1:Nke
            W_inv[:,:] = W_inv[:,:] + (1.0/σ_noise[p]^2)*diagm(H_std[p,:].^2);
        end
    else
        W_inv = diagm(model_un(H_std,σ_noise));
    end

    # return
    inv(H_tilde'*Γ_inv*H_tilde + Vimg'*W_inv*Vimg + D_tilde'*Γ_D_inv*D_tilde + B_tilde'*Γ_B_inv*B_tilde)*((H_tilde'*Γ_inv)*y + B_tilde'*Γ_B_inv*ρ_val), F, W_inv
end

function low_rank_reconstruction_un_2(rankK::Int64,F::SVD{Float64,Float64,Array{Float64,2}},H::Array{Cdouble,2},W_inv::Array{Cdouble,2},Γ_inv::Array{Cdouble,2},D::Array{Cdouble,2},Γ_D_inv::Array{Cdouble,2},B::Array{Cdouble,2},Γ_B_inv::Array{Cdouble,2},y::Array{Cdouble,1},ρ_val::Array{Cdouble,1}) #
    # low rank basis
    Vimg = F.Vt[1:rankK,:]'; # the image of the operator (in the right or state space)

    # low rank operators
    H_tilde = H*Vimg;
    D_tilde = D*Vimg;
    B_tilde = B*Vimg;

    # return
    inv(H_tilde'*Γ_inv*H_tilde + Vimg'*W_inv*Vimg + D_tilde'*Γ_D_inv*D_tilde + B_tilde'*Γ_B_inv*B_tilde)*((H_tilde'*Γ_inv)*y + B_tilde'*Γ_B_inv*ρ_val)
end


# need the result of the previous reconstruction (rankK-1)
function one_iteration_nso_un_2(rankK::Int64,z::Array{Cdouble,1},F::SVD{Float64,Float64,Array{Float64,2}},H::Array{Cdouble,2},W_inv::Array{Cdouble,2},Γ_inv::Array{Cdouble,2},D::Array{Cdouble,2},Γ_D_inv::Array{Cdouble,2},B::Array{Cdouble,2},Γ_B_inv::Array{Cdouble,2},y::Array{Cdouble,1},ρ_val::Array{Cdouble,1})
    # low rank basis
    z_bar = z[1:rankK-1]; # this line should be useless
    Vimg = F.Vt[1:rankK-1,:]'; # the image of the operator (in the right or state space)
    V = F.Vt[rankK:end,:]';

    # low rank operators
    H_tilde = H*Vimg;
    D_tilde = D*Vimg;
    B_tilde = B*Vimg;

    H_bar = H*V;
    D_bar = D*V;
    B_bar = B*V;

    y_bar = y - H_tilde*z_bar;
    y_D = -D_tilde*z_bar;
    ρ_bar = ρ_val - B_tilde*z_bar;
    y_H = -Vimg*z_bar;

    # return
    inv(H_bar'*Γ_inv*H_bar + V'W_inv*V + D_bar'*Γ_D_inv*D_bar + B_bar'*Γ_B_inv*B_bar)*((H_bar'*Γ_inv)*y_bar + V'W_inv*y_H + D_bar'*Γ_D_inv*y_D + B_bar'*Γ_B_inv*ρ_bar)
end

# here it is assumed that the reconstruction in the Img space is known
function null_nso_un_2(rankK::Int64,z::Array{Cdouble,1},F::SVD{Float64,Float64,Array{Float64,2}},W_inv::Array{Cdouble,2},D::Array{Cdouble,2},Γ_D_inv::Array{Cdouble,2},B::Array{Cdouble,2},Γ_B_inv::Array{Cdouble,2},ρ_val::Array{Cdouble,1})
    # low rank basis
    z_bar = z[1:rankK]; # this line should be useless
    Vimg = F.Vt[1:rankK,:]';
    V = F.Vt[rankK+1:end,:]'; # the null space of the operator (in the right or state space)

    # low rank operators
    D_tilde = D*Vimg;
    B_tilde = B*Vimg;

    D_bar = D*V;
    B_bar = B*V;

    y_D = -D_tilde*z_bar;
    ρ_bar = ρ_val - B_tilde*z_bar;
    y_H = -Vimg*z_bar;

    # return
    inv(D_bar'*Γ_D_inv*D_bar + V'W_inv*V + B_bar'*Γ_B_inv*B_bar)*(D_bar'*Γ_D_inv*y_D + V'W_inv*y_H + B_bar'*Γ_B_inv*ρ_bar)
end


function iterative_nso_un_2(rankK::Int64,H::Array{Cdouble,2},H_std::Array{Cdouble,2},σ_noise::Array{Cdouble,1},Γ_inv::Array{Cdouble,2},D::Array{Cdouble,2},Γ_D_inv::Array{Cdouble,2},B::Array{Cdouble,2},Γ_B_inv::Array{Cdouble,2},y::Array{Cdouble,1},ρ_val::Array{Cdouble,1})
    # init
    z,F,W_inv = low_rank_reconstruction_un_2(1,H,H_std,σ_noise,Γ_inv,D,Γ_D_inv,B,Γ_B_inv,y,ρ_val);
    z_low_rank = zeros(Cdouble,rankK);
    z_low_rank[1] = z[1];
    # iteration
    for i in 2:rankK
       # assume the z[1:i-1] known
       zz = one_iteration_nso_un_2(i,z_low_rank[1:i-1],F,H,W_inv,Γ_inv,D,Γ_D_inv,B,Γ_B_inv,y,ρ_val);
       z_low_rank[i] = zz[1];
    end
    # the rest of the reconstruction should rely on the a priori only, not the measurement model
    x_nso = null_nso_un_2(rankK,z_low_rank,F,W_inv,D,Γ_D_inv,B,Γ_B_inv,ρ_val);

    # return the reconstruction, the SVD coefficients, the SVD and the diagonal matrix for the model uncertainty
    F.V*[z_low_rank;x_nso], [z_low_rank; x_nso], F, W_inv
end

function iterative_nso_un_2(rankK::Int64,F::SVD{Float64,Float64,Array{Float64,2}},H::Array{Cdouble,2},W_inv::Array{Cdouble,2},Γ_inv::Array{Cdouble,2},D::Array{Cdouble,2},Γ_D_inv::Array{Cdouble,2},B::Array{Cdouble,2},Γ_B_inv::Array{Cdouble,2},y::Array{Cdouble,1},ρ_val::Array{Cdouble,1})
    # init
    z = low_rank_reconstruction_un_2(1,F,H,W_inv,Γ_inv,D,Γ_D_inv,B,Γ_B_inv,y,ρ_val);
    z_low_rank = zeros(Cdouble,rankK);
    z_low_rank[1] = z[1];
    # iteration
    for i in 2:rankK
       # assume the z[1:i-1] known
       zz = one_iteration_nso_un_2(i,z_low_rank[1:i-1],F,H,W_inv,Γ_inv,D,Γ_D_inv,B,Γ_B_inv,y,ρ_val);
       z_low_rank[i] = zz[1];
    end
    # the rest of the reconstruction should rely on the a priori only, not the measurement model
    x_nso = null_nso_un_2(rankK,z_low_rank,F,W_inv,D,Γ_D_inv,B,Γ_B_inv,ρ_val);

    # return the reconstruction
    F.V*[z_low_rank;x_nso]
end


function data_sample_nso_un(Ns::Int64,rankK::Int64,H::Array{Cdouble,2},H_std::Array{Cdouble,2},σ_noise::Array{Cdouble,1},Γ_inv::Array{Cdouble,2},D::Array{Cdouble,2},Γ_D_inv::Array{Cdouble,2},B::Array{Cdouble,2},Γ_B_inv::Array{Cdouble,2},y::Array{Cdouble,1},ρ_val::Array{Cdouble,1})
    Nke,M = size(H);
    # draw sample
    ρ_samples = Array{Cdouble}(undef,Ns,M);
    ρ_samples[1,:], _, F, W_inv = iterative_nso_un_2(rankK,H,H_std,σ_noise,Γ_inv,D,Γ_D_inv,B,Γ_B_inv,y,ρ_val)
    for n in 2:Ns
        ρ_samples[n,:] = iterative_nso_un_2(rankK,F,H,W_inv,Γ_inv,D,Γ_D_inv,B,Γ_B_inv,y+σ_noise.*randn(Nke),ρ_val)
    end
    dropdims(mean(ρ_samples,dims=1),dims=1), dropdims(std(ρ_samples,dims=1),dims=1), ρ_samples, F, W_inv
end


function shuffle_data(Nke::Int64,Npeak::Int64; Nmin::Int64=4)
   # return the idex array to re-arrange the data

   Ns = floor(Int64,Nke/Npeak);
   N_sample_peak = rand(Nmin:Ns,Npeak);
   idx = zeros(Int64,sum(N_sample_peak));
   # draw the samples for each peak
   ii = 1;
   for i in 1:Npeak
      # println(length(collect(ii:ii+N_sample_peak[i]-1)))
      idx[ii:ii+N_sample_peak[i]-1] = rand((i-1)*Ns+1:i*Ns,N_sample_peak[i])
      ii = ii + N_sample_peak[i]
   end
   idx
end

function shuffle_data_sample_nso_un(Ns::Int64,Npeak::Int64,rankK::Int64,H::Array{Cdouble,2},H_std::Array{Cdouble,2},σ_noise::Array{Cdouble,1},Γ_inv::Array{Cdouble,2},D::Array{Cdouble,2},Γ_D_inv::Array{Cdouble,2},B::Array{Cdouble,2},Γ_B_inv::Array{Cdouble,2},y::Array{Cdouble,1},ρ_val::Array{Cdouble,1};Nmin::Int64=4)
    Nke,M = size(H);
    # draw sample
    ρ_samples = Array{Cdouble}(undef,Ns,M);
    ρ_samples[1,:], _, F, W_inv = iterative_nso_un_2(rankK,H,H_std,σ_noise,Γ_inv,D,Γ_D_inv,B,Γ_B_inv,y,ρ_val)
    for n in 2:Ns
        idx = shuffle_data(Nke,Npeak; Nmin=Nmin)
        ρ_samples[n,:], _, _, _ = iterative_nso_un_2(rankK,H[idx,:],H_std[idx,:],σ_noise[idx],Γ_inv[idx,idx],D,Γ_D_inv,B,Γ_B_inv,y[idx],ρ_val)
        # ρ_samples[n,:] = iterative_nso_un_2(rankK,F,H,W_inv,Γ_inv,D,Γ_D_inv,B,Γ_B_inv,y+σ_noise.*randn(Nke),ρ_val)
    end
    dropdims(mean(ρ_samples,dims=1),dims=1), dropdims(std(ρ_samples,dims=1),dims=1), ρ_samples, F, W_inv
end
