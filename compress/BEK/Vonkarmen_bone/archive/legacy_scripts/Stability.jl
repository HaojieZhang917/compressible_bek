module CRC_STA
    export  COF, Spatial_mode_BEK
    export  assemble_time_mat, assemble_mat, assemble_adjmat
    export  boudary_condition, eig_full, cheb_points
    using LinearAlgebra
    struct COF
            Ta :: Matrix{ComplexF64}
            A :: Matrix{ComplexF64}
            B :: Matrix{ComplexF64}
            C :: Matrix{ComplexF64}
            dC :: Matrix{ComplexF64}
            D1 :: Matrix{ComplexF64}
            Vxx :: Matrix{ComplexF64}
            Vyy :: Matrix{ComplexF64}
            Vzz :: Matrix{ComplexF64}
            dVzz :: Matrix{ComplexF64}
            d2Vzz :: Matrix{ComplexF64}
            Vxy :: Matrix{ComplexF64}
            Vxz :: Matrix{ComplexF64}
            dVxz :: Matrix{ComplexF64}
            Vyz :: Matrix{ComplexF64}
            dVyz :: Matrix{ComplexF64}
    end
    function Spatial_mode_BEK(F,G,H,T,sigma,N_cheb,D,D2,R)
        dT = D * T
        size = N_cheb + 1
        eye = I(N_cheb+1)
        Zero = zeros(N_cheb+1,N_cheb+1)
        # Five disturbance variables: [u, v, w, T, p]
        Ta = zeros(5*size,5*size)
        A = zeros(5*size,5*size)
        B = zeros(5*size,5*size)
        C = zeros(5*size,5*size)
        D1 = zeros(5*size,5*size)
        Vxy = zeros(5*size,5*size)
        Vxz = zeros(5*size,5*size)
        Vyz = zeros(5*size,5*size)
        Vxx = zeros(5*size,5*size)
        Vyy = zeros(5*size,5*size)
        Vzz = zeros(5*size,5*size)
        dVxz = zeros(5*size,5*size)
        dVyz = zeros(5*size,5*size)
        # Boussinesq model consistent with Bone.py:
        #     F'' = F^2 + H*F' - (G - 1)^2 - (T - 1)
        #     G'' = 2*F*G + H*G' - 2*F
        # Bone.py uses T = T_dim/T_inf.  In this local stability scaling,
        # non-derivative body-force/curvature terms enter the disturbance
        # momentum equations with the same 1/R factor as the Coriolis terms.
        # For an ideal gas beta*T_inf = 1, so temperature-disturbance
        # centrifugal buoyancy contributes -(1/R)*theta in the radial equation.
        dF = D * F
        dG = D * G
        # Ta: time-derivative coefficients for [u, v, w, T, p]
        Ta_11 = eye; Ta_12 = Ta_13 = Ta_14 = Ta_15 = Zero
        Ta_22 = eye; Ta_21 = Ta_23 = Ta_24 = Ta_25 = Zero
        Ta_33 = eye; Ta_31 = Ta_32 = Ta_34 = Ta_35 = Zero
        Ta_44 = eye; Ta_41 = Ta_42 = Ta_43 = Ta_45 = Zero
        Ta_51 = Ta_52 = Ta_53 = Ta_54 = Ta_55 = Zero   # p has no time derivative
        # A: alpha coefficients for [u, v, w, T, p]
       A_11 = F .* eye; A_12 = A_13 = A_14 = Zero; A_15 = eye    # u: F*u + p_x
       A_22 = F .* eye; A_21 = A_23 = A_24 = A_25 = Zero         # v: F*v
       A_33 = F .* eye; A_31 = A_32 = A_34 = A_35 = Zero         # w: F*w
        A_44 = F .* eye; A_41 = A_42 = A_43 = A_45 = Zero               # T: F*T
        A_51 = eye; A_52 = A_53 = A_54 = A_55 = Zero                     # p: u_x (
        # B: beta coefficients for [u, v, w, T, p]
        B_11 = (1/R) * G .* eye; B_12 = B_13 = B_14 = B_15 = Zero       # u: G*u/R
        B_22 = (1/R) * G .* eye; B_21 = B_23 = Zero; B_24 = Zero; B_25 = (1/R) .* eye  # v: G*v/R + p_y/R
        B_33 = (1/R) * G .* eye; B_31 = B_32 = B_34 = B_35 = Zero       # w: G*w/R
        B_44 = (1/R) * G .* eye; B_41 = B_43 = B_45 = Zero; B_42 = Zero # T: G*T
        B_51 = B_53 = B_54 = B_55 = Zero; B_52 = (1/R) .* eye            # p: v_y 

        # C: wall-normal advection coefficients for [u, v, w, T, p]
        C_11 = (1/R) * H .* eye; dC_11 = D * diag(C_11) .* eye
        C_12 = C_13 = C_14 = C_15 = Zero; dC_12 = dC_13 = dC_14 = dC_15 = Zero
        C_22 = (1/R) * H .* eye; dC_22 = D * diag(C_22) .* eye
        C_21 = C_23 = C_24 = C_25 = Zero; dC_21 = dC_23 = dC_24 = dC_25 = Zero
        C_33 = (1/R) * H .* eye; dC_33 = D * diag(C_33) .* eye
        C_31 = C_32 = Zero; C_34 = Zero; C_35 = eye;  # w: pressure z-gradient
        dC_31 = dC_32 = Zero; dC_34 = dC_35 = Zero
        C_44 = (1/R) * H .* eye; dC_44 = D * diag(C_44) .* eye
        C_41 = C_42 = C_43 = C_45 = Zero; dC_41 = dC_42 = dC_43 = dC_45 = Zero
        # Continuity: w_z contribution
        C_53 = eye; dC_53 = D * diag(C_53) .* eye
        C_51 = C_52 = C_54 = C_55 = Zero; dC_51 = dC_52 = dC_54 = dC_55 = Zero

        # D1: base-flow shear, curvature, and buoyancy coupling for [u, v, w, T, p]
        D_11 = (1/R) * F .* eye
        D_12 = -(1/R) * 2 .* (G .+ 1) .* eye
        D_13 = dF .* eye
        D_14 = -(1/R) .* eye
        D_15 = Zero
        D_21 = (1/R) * 2 .* (G .+ 1) .* eye
        D_22 = (1/R) * F .* eye
        D_23 = dG .* eye
        D_24 = Zero
        D_25 = Zero
        D_31 = D_32 = Zero; D_34 = Zero  # current base-flow model has no axial thermal-inertial term
        D_33 = (1/R) * (D * H) .* eye; D_35 = Zero
        D_41 = D_42 = Zero; D_44 = Zero; D_45 = Zero
        D_43 = dT .* eye                       # T: T_base' * w
        D_51 = 1/R .* eye; D_52 = D_53 = D_54 = D_55 = Zero  # p: u_r

        # Vxx: alpha-squared diffusion terms for [u, v, w, T, p]
        Vxx_11 = -(1/R) * eye; Vxx_12 = Vxx_13 = Vxx_14 = Vxx_15 = Zero
        Vxx_22 = -(1/R) * eye; Vxx_21 = Vxx_23 = Vxx_24 = Vxx_25 = Zero
        Vxx_33 = -(1/R) * eye; Vxx_31 = Vxx_32 = Vxx_34 = Vxx_35 = Zero
        Vxx_44 = -(1/(sigma*R)) * eye; Vxx_41 = Vxx_42 = Vxx_43 = Vxx_45 = Zero
        Vxx_51 = Vxx_52 = Vxx_53 = Vxx_54 = Vxx_55 = Zero

        # Vyy: beta-squared diffusion terms
        Vyy_11 = -(1/R^3) * eye; Vyy_12 = Vyy_13 = Vyy_14 = Vyy_15 = Zero
        Vyy_22 = -(1/R^3) * eye; Vyy_21 = Vyy_23 = Vyy_24 = Vyy_25 = Zero
        Vyy_33 = -(1/R^3) * eye; Vyy_31 = Vyy_32 = Vyy_34 = Vyy_35 = Zero
        Vyy_44 = -(1/(sigma*R^3)) * eye; Vyy_41 = Vyy_42 = Vyy_43 = Vyy_45 = Zero
        Vyy_51 = Vyy_52 = Vyy_53 = Vyy_54 = Vyy_55 = Zero

        # Vzz: wall-normal diffusion terms
        Vzz_11 = -(1/R) * eye; dVzz_11 = D*diag(Vzz_11).*eye; d2Vzz_11 = D2*diag(Vzz_11).*eye
        Vzz_12 = Vzz_13 = Vzz_14 = Vzz_15 = Zero; dVzz_12 = dVzz_13 = dVzz_14 = dVzz_15 = Zero; d2Vzz_12 = d2Vzz_13 = d2Vzz_14 = d2Vzz_15 = Zero
        Vzz_22 = -(1/R) * eye; dVzz_22 = D*diag(Vzz_22).*eye; d2Vzz_22 = D2*diag(Vzz_22).*eye
        Vzz_21 = Vzz_23 = Vzz_24 = Vzz_25 = Zero; dVzz_21 = dVzz_23 = dVzz_24 = dVzz_25 = Zero; d2Vzz_21 = d2Vzz_23 = d2Vzz_24 = d2Vzz_25 = Zero
        Vzz_33 = -(1/R) * eye; dVzz_33 = D*diag(Vzz_33).*eye; d2Vzz_33 = D2*diag(Vzz_33).*eye
        Vzz_31 = Vzz_32 = Vzz_34 = Vzz_35 = Zero; dVzz_31 = dVzz_32 = dVzz_34 = dVzz_35 = Zero; d2Vzz_31 = d2Vzz_32 = d2Vzz_34 = d2Vzz_35 = Zero
        Vzz_44 = -(1/(sigma*R)) * eye; dVzz_44 = D*diag(Vzz_44).*eye; d2Vzz_44 = D2*diag(Vzz_44).*eye
        Vzz_41 = Vzz_42 = Vzz_43 = Vzz_45 = Zero; dVzz_41 = dVzz_42 = dVzz_43 = dVzz_45 = Zero; d2Vzz_41 = d2Vzz_42 = d2Vzz_43 = d2Vzz_45 = Zero
        Vzz_51 = Vzz_52 = Vzz_53 = Vzz_54 = Vzz_55 = Zero; dVzz_51 = dVzz_52 = dVzz_53 = dVzz_54 = dVzz_55 = Zero; d2Vzz_51 = d2Vzz_52 = d2Vzz_53 = d2Vzz_54 = d2Vzz_55 = Zero

        Vxy = zeros(5*size,5*size)
        dVxy = zeros(5*size,5*size)
        Vxz = zeros(5*size,5*size)
        Vyz = zeros(5*size,5*size)

        # Assemble 5-by-5 block matrices.
        Ta = [Ta_11 Ta_12 Ta_13 Ta_14 Ta_15; Ta_21 Ta_22 Ta_23 Ta_24 Ta_25; Ta_31 Ta_32 Ta_33 Ta_34 Ta_35; Ta_41 Ta_42 Ta_43 Ta_44 Ta_45; Ta_51 Ta_52 Ta_53 Ta_54 Ta_55]
        A  = [A_11  A_12  A_13  A_14  A_15;  A_21  A_22  A_23  A_24  A_25;  A_31  A_32  A_33  A_34  A_35;  A_41  A_42  A_43  A_44  A_45;  A_51  A_52  A_53  A_54  A_55]
        B  = [B_11  B_12  B_13  B_14  B_15;  B_21  B_22  B_23  B_24  B_25;  B_31  B_32  B_33  B_34  B_35;  B_41  B_42  B_43  B_44  B_45;  B_51  B_52  B_53  B_54  B_55]
        C  = [C_11  C_12  C_13  C_14  C_15;  C_21  C_22  C_23  C_24  C_25;  C_31  C_32  C_33  C_34  C_35;  C_41  C_42  C_43  C_44  C_45;  C_51  C_52  C_53  C_54  C_55]
        dC = [dC_11 dC_12 dC_13 dC_14 dC_15; dC_21 dC_22 dC_23 dC_24 dC_25; dC_31 dC_32 dC_33 dC_34 dC_35; dC_41 dC_42 dC_43 dC_44 dC_45; dC_51 dC_52 dC_53 dC_54 dC_55]
        D1 = [D_11  D_12  D_13  D_14  D_15;  D_21  D_22  D_23  D_24  D_25;  D_31  D_32  D_33  D_34  D_35;  D_41  D_42  D_43  D_44  D_45;  D_51  D_52  D_53  D_54  D_55]
        Vxx= [Vxx_11 Vxx_12 Vxx_13 Vxx_14 Vxx_15; Vxx_21 Vxx_22 Vxx_23 Vxx_24 Vxx_25; Vxx_31 Vxx_32 Vxx_33 Vxx_34 Vxx_35; Vxx_41 Vxx_42 Vxx_43 Vxx_44 Vxx_45; Vxx_51 Vxx_52 Vxx_53 Vxx_54 Vxx_55]
        Vyy= [Vyy_11 Vyy_12 Vyy_13 Vyy_14 Vyy_15; Vyy_21 Vyy_22 Vyy_23 Vyy_24 Vyy_25; Vyy_31 Vyy_32 Vyy_33 Vyy_34 Vyy_35; Vyy_41 Vyy_42 Vyy_43 Vyy_44 Vyy_45; Vyy_51 Vyy_52 Vyy_53 Vyy_54 Vyy_55]
        Vzz= [Vzz_11 Vzz_12 Vzz_13 Vzz_14 Vzz_15; Vzz_21 Vzz_22 Vzz_23 Vzz_24 Vzz_25; Vzz_31 Vzz_32 Vzz_33 Vzz_34 Vzz_35; Vzz_41 Vzz_42 Vzz_43 Vzz_44 Vzz_45; Vzz_51 Vzz_52 Vzz_53 Vzz_54 Vzz_55]
        dVzz=[dVzz_11 dVzz_12 dVzz_13 dVzz_14 dVzz_15; dVzz_21 dVzz_22 dVzz_23 dVzz_24 dVzz_25; dVzz_31 dVzz_32 dVzz_33 dVzz_34 dVzz_35; dVzz_41 dVzz_42 dVzz_43 dVzz_44 dVzz_45; dVzz_51 dVzz_52 dVzz_53 dVzz_54 dVzz_55]
        d2Vzz=[d2Vzz_11 d2Vzz_12 d2Vzz_13 d2Vzz_14 d2Vzz_15; d2Vzz_21 d2Vzz_22 d2Vzz_23 d2Vzz_24 d2Vzz_25; d2Vzz_31 d2Vzz_32 d2Vzz_33 d2Vzz_34 d2Vzz_35; d2Vzz_41 d2Vzz_42 d2Vzz_43 d2Vzz_44 d2Vzz_45; d2Vzz_51 d2Vzz_52 d2Vzz_53 d2Vzz_54 d2Vzz_55]
        return COF(Ta,A,B,C,dC,D1,Vxx,Vyy,Vzz,dVzz,d2Vzz,Vxy,Vxz,dVxz,Vyz,dVyz)
    end
    function assemble_time_mat(cof,D,D2,be,alpha,R,N_cheb)
        H0 = cof.D1  + im * R * be * cof.B - be^2 * R^2 * cof.Vyy + (cof.C .+ im * be * R * cof.Vyz) * kron(I(5), D)  + (cof.Vzz) * kron(I(5),D2) + 
        alpha * (im * cof.A - be * R * cof.Vxy + im *  cof.Vxz * kron(I(5),D)) + alpha^2 * (-cof.Vxx)
        H1 = im * cof.Ta
        H0 = H0[setdiff(1:end , (1,N_cheb+1,N_cheb+2,2N_cheb+2,2N_cheb+3,3N_cheb+3,3N_cheb+4,4N_cheb+4)),setdiff(1:end , (1,N_cheb+1,N_cheb+2,2N_cheb+2,2N_cheb+3,3N_cheb+3,3N_cheb+4,4N_cheb+4))]
        H1 = H1[setdiff(1:end , (1,N_cheb+1,N_cheb+2,2N_cheb+2,2N_cheb+3,3N_cheb+3,3N_cheb+4,4N_cheb+4)),setdiff(1:end , (1,N_cheb+1,N_cheb+2,2N_cheb+2,2N_cheb+3,3N_cheb+3,3N_cheb+4,4N_cheb+4))]
        return H0, H1
    end
    function assemble_mat(cof,D,D2,be,omega,R)
        L0 = cof.D1  + im * R * be * cof.B - im * omega * cof.Ta - be^2 * R^2 * cof.Vyy + (cof.C .+ im * be * R * cof.Vyz) * kron(I(5), D)  + (cof.Vzz) * kron(I(5),D2) 
        L1 = im * cof.A - be * R * cof.Vxy + im *  cof.Vxz * kron(I(5),D)
        L2 = -cof.Vxx 
        return L0,L1,L2
    end
    function assemble_adjmat(cof,D,D2,be,omega,R)
        A0_raw = transpose(cof.D1) + (im * R * be * transpose(cof.B)) - (im * omega * transpose(cof.Ta)) - (be^2 * R^2 * transpose(cof.Vyy)) - transpose(cof.dC) - (im *be*transpose(cof.dVyz)) + transpose(cof.d2Vzz) - (transpose(cof.C) + im * be * R * transpose(cof.Vyz) - 2 * transpose(cof.dVzz)) * kron(I(5),D) + transpose(cof.Vzz) * kron(I(5),D2)
        A1_raw = (im * transpose(cof.A)) - (be * R * transpose(cof.Vxy)) - (im * transpose(cof.dVxz)) - (im * transpose(cof.Vxz)) * kron(I(5),D) 
        A2_raw = -transpose(cof.Vxx)
        return A0_raw,A1_raw,A2_raw
    end
    function boudary_condition(L0,L1,L2,N_cheb)
        L0 = L0[setdiff(1:end , (1,N_cheb+1,N_cheb+2,2N_cheb+2,2N_cheb+3,3N_cheb+3,3N_cheb+4,4N_cheb+4,5N_cheb+5)),setdiff(1:end , (1,N_cheb+1,N_cheb+2,2N_cheb+2,2N_cheb+3,3N_cheb+3,3N_cheb+4,4N_cheb+4,5N_cheb+5))]
        L1 = L1[setdiff(1:end , (1,N_cheb+1,N_cheb+2,2N_cheb+2,2N_cheb+3,3N_cheb+3,3N_cheb+4,4N_cheb+4,5N_cheb+5)),setdiff(1:end , (1,N_cheb+1,N_cheb+2,2N_cheb+2,2N_cheb+3,3N_cheb+3,3N_cheb+4,4N_cheb+4,5N_cheb+5))]
        L2 = L2[setdiff(1:end , (1,N_cheb+1,N_cheb+2,2N_cheb+2,2N_cheb+3,3N_cheb+3,3N_cheb+4,4N_cheb+4,5N_cheb+5)),setdiff(1:end , (1,N_cheb+1,N_cheb+2,2N_cheb+2,2N_cheb+3,3N_cheb+3,3N_cheb+4,4N_cheb+4,5N_cheb+5))]
        return L0,L1,L2
    end
    function eig_full(eigvec,N_cheb,num)
        N = N_cheb + 1
        reduced = eigvec[:,num]
        removed = [1, N, N+1, 2N, 2N+1, 3N, 3N+1, 4N, 5N]
        keep = setdiff(collect(1:5N), removed)
        length(reduced) == length(keep) || throw(DimensionMismatch(
            "eigenvector length is inconsistent with boudary_condition",
        ))
        full = zeros(ComplexF64, 5N)
        full[keep] .= reduced
        u = full[1:N]
        v = full[N+1:2N]
        w = full[2N+1:3N]
        T = full[3N+1:4N]
        p = full[4N+1:5N]
        return (u,v,w,T,p)
    end
    function cheb_points(N)
        theta = range(0, stop=pi, length=N+1)
        x = -cos.(theta)
        cc = [2; ones(N-1); 2] .* (-1.0).^(0:N)
        X = repeat(x, 1, N+1)
        dX = X .- X'
        Id = diagm(0 => ones(N+1))
        DM = (cc * (1.0 ./ cc)') ./ (dX .+ Id)
        DM = DM .- diagm(vec(sum(DM, dims=2)))
        a = 2.0; b = 0.6; cmap = 0.5
        for i in 1:N+1
            xi = x[i]
            map_der = b + 3*(1-b)*xi^2 - 2*cmap*(1-b)*xi
            denom = 1 - b*xi - (1-b)*(xi^3 + cmap*(1-xi^2))
            DM[i, :] .*= denom^2 / (2*a*map_der)
        end
        for i in 1:N+1
            xi = x[i]
            denom = 1 - b*xi - (1-b)*(xi^3 + cmap*(1-xi^2))
            x[i] = a * (1 + b*xi + (1-b)*(xi^3 + cmap*(1-xi^2))) / denom
            x[i] = min(x[i], 20.0)
        end
        D2M = DM^2
        return DM, D2M, x
    end
    function cheb_points(N)
        theta = range(0, stop=pi, length=N+1)
        x = -cos.(theta)
        cc = [2; ones(N-1); 2] .* (-1.0).^(0:N)
        X = repeat(x, 1, N+1)
        dX = X .- X'
        Id = diagm(0 => ones(N+1))
        DM = (cc * (1.0 ./ cc)') ./ (dX .+ Id)
        DM = DM .- diagm(vec(sum(DM, dims=2)))
        a = 2.0; b = 0.6; cmap = 0.5
        for i in 1:N+1
            xi = x[i]
            map_der = b + 3*(1-b)*xi^2 - 2*cmap*(1-b)*xi
            denom = 1 - b*xi - (1-b)*(xi^3 + cmap*(1-xi^2))
            DM[i, :] .*= denom^2 / (2*a*map_der)
        end
        for i in 1:N+1
            xi = x[i]
            denom = 1 - b*xi - (1-b)*(xi^3 + cmap*(1-xi^2))
            x[i] = a * (1 + b*xi + (1-b)*(xi^3 + cmap*(1-xi^2))) / denom
            x[i] = min(x[i], 20.0)
        end
        D2M = DM^2
        return DM, D2M, x
    end
end
