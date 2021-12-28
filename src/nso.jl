function low_rank_reconstruction(rankK::Int64,H::Array{Cdouble,2},Γ_inv::Array{Cdouble,2},D::Array{Cdouble,2},Γ_D_inv::Array{Cdouble,2},y::Array{Cdouble,1})
    # singular value decomposition
    F = svd(H, full=true);

    # low rank basis
    Vimg = F.Vt[1:rankK,:]'; # the image of the operator (in the right or state space)

    # low rank operators
    H_tilde = H*Vimg;
    D_tilde = D*Vimg;

    # return
    inv(H_tilde'*Γ_inv*H_tilde + D_tilde'*Γ_D_inv*D_tilde)*(H_tilde'*Γ_inv)*y, F
end



function low_rank_reconstruction(rankK::Int64,H::Array{Cdouble,2},Γ_inv::Array{Cdouble,2},D::Array{Cdouble,2},Γ_D_inv::Array{Cdouble,2},B::Array{Cdouble,2},Γ_B_inv::Array{Cdouble,2},y::Array{Cdouble,1},ρ_val::Array{Cdouble,1}) #
    # singular value decomposition
    F = svd(H, full=true);

    # low rank basis
    Vimg = F.Vt[1:rankK,:]'; # the image of the operator (in the right or state space)

    # low rank operators
    H_tilde = H*Vimg;
    D_tilde = D*Vimg;
    B_tilde = B*Vimg;

    # return
    inv(H_tilde'*Γ_inv*H_tilde + D_tilde'*Γ_D_inv*D_tilde + B_tilde'*Γ_B_inv*B_tilde)*((H_tilde'*Γ_inv)*y + B_tilde'*Γ_B_inv*ρ_val), F
end


function rest_reconstruction(rankK::Int64,rankZ::Int64,H::Array{Cdouble,2},Γ_inv::Array{Cdouble,2},D::Array{Cdouble,2},Γ_D_inv::Array{Cdouble,2},B::Array{Cdouble,2},Γ_B_inv::Array{Cdouble,2},y::Array{Cdouble,1},ρ_val::Array{Cdouble,1}) #
    # first low rank approximation (rankZ) and singular value decomposition
    z,F = low_rank_reconstruction(rankZ,H,Γ_inv,D,Γ_D_inv,B,Γ_B_inv,y,ρ_val) # keep only the rankK first entries of the rankZ approximation

    # low rank basis
    z_bar = z[1:rankK];
    Vimg = F.Vt[1:rankK,:]';
    V = F.Vt[rankK+1:end,:]'; # the image of the operator (in the right or state space)

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

    # return
    inv(H_bar'*Γ_inv*H_bar + D_bar'*Γ_D_inv*D_bar + B_bar'*Γ_B_inv*B_bar)*((H_bar'*Γ_inv)*y_bar + D_bar'*Γ_D_inv*y_D + B_bar'*Γ_B_inv*ρ_bar),z, F
end

# need the result of the previous reconstruction (rankK-1)
function one_iteration_nso(rankK::Int64,z::Array{Cdouble,1},F::SVD{Float64,Float64,Array{Float64,2}},H::Array{Cdouble,2},Γ_inv::Array{Cdouble,2},D::Array{Cdouble,2},Γ_D_inv::Array{Cdouble,2},B::Array{Cdouble,2},Γ_B_inv::Array{Cdouble,2},y::Array{Cdouble,1},ρ_val::Array{Cdouble,1})
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

    # return
    inv(H_bar'*Γ_inv*H_bar + D_bar'*Γ_D_inv*D_bar + B_bar'*Γ_B_inv*B_bar)*((H_bar'*Γ_inv)*y_bar + D_bar'*Γ_D_inv*y_D + B_bar'*Γ_B_inv*ρ_bar)
end

# here it is assumed that the reconstruction in the Img space is known
function null_nso(rankK::Int64,z::Array{Cdouble,1},F::SVD{Float64,Float64,Array{Float64,2}},D::Array{Cdouble,2},Γ_D_inv::Array{Cdouble,2},B::Array{Cdouble,2},Γ_B_inv::Array{Cdouble,2},ρ_val::Array{Cdouble,1})
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

    # return
    inv(D_bar'*Γ_D_inv*D_bar + B_bar'*Γ_B_inv*B_bar)*(D_bar'*Γ_D_inv*y_D + B_bar'*Γ_B_inv*ρ_bar)
end


function iterative_nso(rankK::Int64,H::Array{Cdouble,2},Γ_inv::Array{Cdouble,2},D::Array{Cdouble,2},Γ_D_inv::Array{Cdouble,2},B::Array{Cdouble,2},Γ_B_inv::Array{Cdouble,2},y::Array{Cdouble,1},ρ_val::Array{Cdouble,1})
    # init
    z,F = low_rank_reconstruction(1,H,Γ_inv,D,Γ_D_inv,B,Γ_B_inv,y,ρ_val);
    z_low_rank = zeros(Cdouble,rankK);
    z_low_rank[1] = z[1];
    # iteration
    for i in 2:rankK
       # assume the z[1:i-1] known
       zz = one_iteration_nso(i,z_low_rank[1:i-1],F,H,Γ_inv,D,Γ_D_inv,B,Γ_B_inv,y,ρ_val);
       z_low_rank[i] = zz[1];
    end
    # the rest of the reconstruction should rely on the a priori only, not the measurement model
    x_nso = null_nso(rankK,z_low_rank,F,D,Γ_D_inv,B,Γ_B_inv,ρ_val);

    # return the reconstruction, the SVD coefficients and the SVD
    F.V*[z_low_rank;x_nso], [z_low_rank; x_nso], F
end

##
## introducing operator uncertainty
##
function model_un_cost_img()
end

function model_un_mat_img(σ_noise::Array{Cdouble,1},H_std::Array{Cdouble,2},V::Union{Array{Cdouble,2},Adjoint{Float64,Array{Float64,2}}},rankK::Int64)
    XX = zeros(Cdouble,rankK,rankK);
    Nke = length(σ_noise);
    M = size(H_std,2);

    for p in 1:Nke
        for j in 1:M
            XX[:,:] = XX[:,:] + (H_std[p,j]^2/σ_noise[p]^2)*V[j,1:rankK]* V[j,1:rankK]'
        end
    end
    XX
end


function low_rank_reconstruction_un(rankK::Int64,H::Array{Cdouble,2},H_std::Array{Cdouble,2},σ_noise::Array{Cdouble,1},Γ_inv::Array{Cdouble,2},D::Array{Cdouble,2},Γ_D_inv::Array{Cdouble,2},B::Array{Cdouble,2},Γ_B_inv::Array{Cdouble,2},y::Array{Cdouble,1},ρ_val::Array{Cdouble,1}) #
    # singular value decomposition
    F = svd(H, full=true);

    # low rank basis
    Vimg = F.Vt[1:rankK,:]'; # the image of the operator (in the right or state space)

    # low rank operators
    H_tilde = H*Vimg;
    D_tilde = D*Vimg;
    B_tilde = B*Vimg;

    # model uncertainty
    K_tilde = model_un_mat_img(σ_noise,H_std,F.V,rankK)

    # return
    inv(H_tilde'*Γ_inv*H_tilde + D_tilde'*Γ_D_inv*D_tilde + B_tilde'*Γ_B_inv*B_tilde + K_tilde)*((H_tilde'*Γ_inv)*y + B_tilde'*Γ_B_inv*ρ_val), F
end

#TODO and still to be done
function one_iteration_nso_un(rankK::Int64,z::Array{Cdouble,1},F::SVD{Float64,Float64,Array{Float64,2}},H::Array{Cdouble,2},H_std::Array{Cdouble,2},σ_noise::Array{Cdouble,1},Γ_inv::Array{Cdouble,2},D::Array{Cdouble,2},Γ_D_inv::Array{Cdouble,2},B::Array{Cdouble,2},Γ_B_inv::Array{Cdouble,2},y::Array{Cdouble,1},ρ_val::Array{Cdouble,1})
    # low rank basis
    z_bar = z[1:rankK-1]; # this line should be useless
    Vimg = F.Vt[1:rankK-1,:]'; # the image of the operator (in the right or state space)
    V = F.Vt[rankK:end,:]';

    # low rank operators
    H_tilde = H*Vimg;
    D_tilde = D*Vimg;
    B_tilde = B*Vimg;

    # model uncertainty
    K_tilde = model_un_mat_img(σ_noise,H_std,F.V,rankK)

    H_bar = H*V;
    D_bar = D*V;
    B_bar = B*V;

    y_bar = y - H_tilde*z_bar;
    y_D = -D_tilde*z_bar;
    ρ_bar = ρ_val - B_tilde*z_bar;

    # return
    inv(H_bar'*Γ_inv*H_bar + D_bar'*Γ_D_inv*D_bar + B_bar'*Γ_B_inv*B_bar)*((H_bar'*Γ_inv)*y_bar + D_bar'*Γ_D_inv*y_D + B_bar'*Γ_B_inv*ρ_bar)
end



function nso_left_alone_iter(k::Int64,x::Array{Cdouble,1},F::SVD{Float64,Float64,Array{Float64,2}},H::Array{Cdouble,2},H_std::Array{Cdouble,2},σ_noise::Array{Cdouble,1},Γ_inv::Array{Cdouble,2},D::Array{Cdouble,2},Γ_D_inv::Array{Cdouble,2},B::Array{Cdouble,2},Γ_B_inv::Array{Cdouble,2},y::Array{Cdouble,1},ρ_val::Array{Cdouble,1})
    M = length(x)
    Nke = size(H,1);
    if ((k>M) || (k<1))
        throw(": index error")
    end
    # should check the dimension of all matrices before starting anything... TODO

    # for the sake of clarity, create a few variables
    Vk = F.V[:,k];
    V_bar = [F.V[:,1:k-1] F.V[:,k+1:end]];
    x_bar = [x[1:k-1]; x[k+1:end]];
    Hk = H*Vk;
    Dk = D*Vk;
    Bk = B*Vk;

    # compute the right hand side terms
    y_bar = y - H*V_bar*x_bar;
    ρ_bar = ρ_val - B*V_bar*x_bar;
    y_D_bar = -D*V_bar*x_bar;
    ρ_j_k = F.V[:,1:k-1]*x[1:k-1] + F.V[:,k+1:end]*x[k+1:end];
    y_ρ_k = 0.0
    for p in 1:Nke
        for j in 1:M
            y_ρ_k = y_ρ_k + (H_std[p,j]^2/σ_noise[p]^2)*V[j,k]*ρ_j_k[j]
        end
    end

    # compute the left hand side terms
    h_ρ_k = 0.0
    for p in 1:Nke
        for j in 1:M
            h_ρ_k = h_ρ_k + (H_std[p,j]^2/σ_noise[p]^2)*(V[j,k]^2)
        end
    end

    # return the left alone component
    (Hk'*Γ_inv*y_bar - y_ρ_k + Dk'*Γ_D_inv*y_D_bar + Bk'*Γ_B_inv*ρ_bar)/(Hk'*Γ_inv*Hk + h_ρ_k + Dk'*Γ_D_inv*Dk + Bk'*Γ_B_inv*Bk)
end
