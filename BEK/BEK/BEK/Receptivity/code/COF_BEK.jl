    using LinearAlgebra
    using NonlinearEigenproblems
    import .KEB_ADJ
    export KEB_Rec_coff
    function KEB_Rec_coff(β,R,N)
        θ = range(0,length=N+1,stop=pi)
        x = reshape(-cos.(θ), N+1, 1)
        c = [2; ones(N-1, 1) ; 2] .* (-1) .^ (0:N)
        X = repeat(x, 1, N+1);
        dX = X - X';
        D = (c * (1 ./ c)') ./ (dX .+ I(N+1));
        D = D - diagm(vec(sum(D, dims=2))); 
        A0_A,A1_A,A2_A,A3_A,A4_A=KEB_ADJ.KEB_ALST(baseflow,N,ω,β,R,Ro,Co)
        nep = PEP([A0_A,A1_A,A2_A,A3_A,A4_A]); 
        sc=10
        nep1 = shift_and_scale(nep,scale=sc);
        mult_scale = norm(nep1.A[end]);
        nep2 = PEP(nep1.A ./ mult_scale);
        λ1,v1 = iar(nep2,σ=0.05,neigs=20,maxit=500);
        λ_1 = sc*λ1
        min_imag_abs = Inf
        min_index = 0
        for i in 1:20
            eigval = λ_1[i]
            if -0.1<imag(eigval)
                curr_imag_abs = abs(imag(eigval))
                if curr_imag_abs < min_imag_abs
                    min_imag_abs = curr_imag_abs
                    min_index = i
                end
            end
        end
        if min_index==0
            for i in 1:20
                eigval = λ_1[i]
                if -0.1<imag(eigval)<0.1
                    curr_imag_abs = abs(imag(eigval))
                    if curr_imag_abs < min_imag_abs
                        min_imag_abs = curr_imag_abs
                        min_index = i
                    end
                end
            end
        end
        v1=v1[:,min_index]
        A0,A1,A2,A3,A4,dU,dV=KEB_SpatialMode.KEB_LST_ALL(baseflow,N,ω,β,R,Ro,Co)
        nep = PEP([A0,A1,A2,A3,A4]); 
        sc=10
        nep1 = shift_and_scale(nep,scale=sc);
        mult_scale = norm(nep1.A[end]);
        nep2 = PEP(nep1.A ./ mult_scale);
        λ1,v2 = iar(nep2,σ=0.05,neigs=20,maxit=500);
        λ_2 = sc*λ1
        min_imag_abs = Inf
        min_index = 0 
        for i in 1:20
            eigval = λ_2[i]
            if -0.1<imag(eigval) < 0
                curr_imag_abs = abs(imag(eigval))
                if curr_imag_abs < min_imag_abs
                    min_imag_abs = curr_imag_abs
                    min_index = i
                end
            end
        end  
        if min_index==0
            for i in 1:20
                eigval = λ_2[i]
                if -0.2<imag(eigval)<0.2
                    curr_imag_abs = abs(imag(eigval))
                    if curr_imag_abs < min_imag_abs
                        min_imag_abs = curr_imag_abs
                        min_index = i
                    end
                end
            end
        end
        α=λ_2[min_index,1]
        v2=v2[:,min_index]
        L=A1+2*α*A2+3*α^2*A3+4*α^3*A4
        Q=v1'L*v2
        insert!(v1,1,0)
        insert!(v1,2,0)
        insert!(v1,N,0)
        insert!(v1,N+1,0)
        insert!(v1,N+2,0)
        insert!(v1,2N+2,0)
        insert!(v2,1,0)
        insert!(v2,2,0)
        insert!(v2,N,0)
        insert!(v2,N+1,0)
        insert!(v2,N+2,0)
        insert!(v2,2N+2,0)
        q_real=real(v1[1:N+1,1])
        q_imag=imag(v1[1:N+1,1])
        dhi=D*q_imag
        dhr=D*q_real
        η=v1[N+2:2N+2,1]
        ηr=real(η)
        ηi=imag(η)
        A11=α*I(N+1);A12=(1/R)*I(N+1);A13=β*I(N+1);A14=zeros(N+1,N+1)
        A21=A12;A22=-1*A11;A23=A14;A24=-1*A13;
        A31=A24;A32=A34=A14;A33=A11;
        A41=A43=A14;A42=A31;A44=A33
        A=[A11 A12 A13 A14;A21 A22 A23 A24;A31 A32 A33 A34;A41 A42 A43 A44;]
        B1=-1*dhi
        B2=-1*dhr
        B3=ηr
        B4=ηi
        B=[B1;B2;B3;B4]
        t=A\B
        fr=t[1:N+1,1]
        fi=t[N+2:2N+2,1]
        gr=t[2N+3:3N+3,1]
        gi=t[3N+4:4N+4,1]
        g=gr+im*gi
        f=fr+im*fi
        g=g[1:70]
        f=f[1:70]
        df=D[1:70,1:70]*f
        dg=D[1:70,1:70]*g
        x_0=200
        u_wall=-dU[1,1]*exp(-im*x_0*α-2*α^2)
        v_wall=-dV[1,1]*exp(-im*x_0*α-2*α^2)
        w_wall=0
        BC=1/R*(u_wall*df[1,1]+v_wall*dg[1,1])
        Cr=-im*BC/Q
        Cr=abs(Cr)
        return Cr
    end 
