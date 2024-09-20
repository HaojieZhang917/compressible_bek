"""
    This module is used to caculate the adjoint linear stability equation(LST) for rotating disk problems.
    It return to the adjoint operator A0 A1 A2 A3 A4
"""
module KEB_ADJ
    using LinearAlgebra
    using BSplineKit
    using SparseArrays
    export KEB_ALST,KEB_Parameter_Matrix
    function KEB_ALST(baseflow,N,ω,β,R,Ro,Co)
        θ = range(0,length=N+1,stop=pi)
        x = reshape(-cos.(θ), N+1, 1)
        c = [2; ones(N-1, 1) ; 2] .* (-1) .^ (0:N)
        X = repeat(x, 1, N+1);
        dX = X - X';
        D = (c * (1 ./ c)') ./ (dX .+ I(N+1));
        D = D - diagm(vec(sum(D, dims=2))); 
        #interpolation
        f = open( baseflow, "r" )
        n = countlines( f )
        seekstart( f )
        data = zeros(n,6)
        for i = 1:n
            z,w,u,v,du,dv = split( readline( f ), " " ) 
            data[i,1] = parse(Float64,z)
            data[i,2] = parse(Float64,w)
            data[i,3] = parse(Float64,u)
            data[i,4] = parse(Float64,v)  
            data[i,5] = parse(Float64,du)
            data[i,6] = parse(Float64,dv)
        end

        close( f )

        t=data[:,1]
        w=-1*data[:,2]
        u=-1*data[:,3]
        v=-1*data[:,4]
        du=-1*data[:,5]
        dv=-1*data[:,6]
        itpw=itp = interpolate(t, w , BSplineOrder(4))
        itpu=itp = interpolate(t, u , BSplineOrder(4))
        itpv=itp = interpolate(t, v , BSplineOrder(4))
        itpdu=itp = interpolate(t, du , BSplineOrder(4))
        itpdv=itp = interpolate(t, dv , BSplineOrder(4))
        for i=1:N+1
            D[i,:]=D[i,:].*((2*x[i]^3-x[i]^2+3*x[i]-4)^2/(20*(6*x[i]^2-2*x[i]+3)))
        end
        for i=1:N+1
            x[i]=(4*x[i]^3-2*x[i]^2+6*x[i]+12)/(-2*x[i]^3+x[i]^2-3*x[i]+4)
            if x[i]>50
                x[i]=50
            end
        end
        D2=D^2;
        U=zeros(N+1,1)
        V=zeros(N+1,1)
        W=zeros(N+1,1)
        dU=zeros(N+1,1)
        dV=zeros(N+1,1)
        dW=zeros(N+1,1)
        p=zeros(N+1,1)
        for i=1:N+1
            U[i,1]=itpu(x[i])
            V[i,1]=itpv(x[i])
            W[i,1]=itpw(x[i])
            dU[i,1]=itpdu(x[i])
            dV[i,1]=itpdv(x[i])
        end
        dW=-2*U
        ddW=D*dW
        ddU=D*dU;
        L0_1=im*(β^4) .*I(N+1)-R*(β^3)*V.*I(N+1)+R*ω*β^2 .*I(N+1)-im*Ro*ddU .*I(N+1) +(2*R*β*dV .*I(N+1)-im*Ro*β^2 *W .*I(N+1)+im*Ro*ddW .*I(N+1)-2*im*Ro*dU .*I(N+1))*D+(R*β*V .*I(N+1)-R*ω .*I(N+1)-2*im*β^2 .*I(N+1)+2*im*Ro*dW .*I(N+1)-im*Ro*U .*I(N+1))*D2 +im*Ro*W.*D^3+im.*D^4 +im*ddU .*I(N+1) 
        L0_2=im*β*R*dU .*I(N+1)-2*Ro*dV .*I(N+1)-(2*Ro*V .*I(N+1)+Co .*I(N+1))*D
        L0_3=-(2*Ro*V .*I(N+1)+Co .*I(N+1))*D
        L0_4=-im*β^2 .*I(N+1)+R*β*V .*I(N+1)-R*ω .*I(N+1)-im*Ro*U.*I(N+1)+im*Ro*dW .*I(N+1)+im*Ro*W.*D+im*D2
        L1_1=(1/R)*β^2 .*I(N+1)-R*β^2 *U .*I(N+1)+im*β*V .*I(N+1)-im*ω .*I(N+1)+(2*R*dU .*I(N+1)-Ro/R *W .*I(N+1))*D+(R*U .*I(N+1)-(1/R) .*I(N+1))*D2
        L1_2=-im*R*dV .*I(N+1)
        L1_3=zeros(N+1,N+1)
        L1_4=R*U .*I(N+1)
        L2_1=2*im*β^2 .*I(N+1)-R*β*V .*I(N+1)+R*ω.*I(N+1)+im*U .*I(N+1)-im*Ro*W.*D-2*im.*D2
        L2_2=L2_3=L1_3
        L2_4=-im*I(N+1)
        L3_1=(1/R).*I(N+1)-R*U.*I(N+1)
        L3_2=L3_3=L3_4=L1_3
        L4_1=im*I(N+1)
        L4_2=L4_3=L4_4=L1_3
        A0=[L0_1 L0_2 ;L0_3 L0_4 ;]
        A1=[L1_1 L1_2 ;L1_3 L1_4 ;]
        A2=[L2_1 L2_2 ;L2_3 L2_4 ;]
        A3=[L3_1 L3_2 ;L3_3 L3_4 ;]
        A4=[L4_1 L4_2 ;L4_3 L4_4 ;]
        A0=A0[setdiff(1:end , (1,2,N,N+1,N+2,2N+2)),setdiff(1:end , (1,2,N,N+1,N+2,2N+2))]
        A1=A1[setdiff(1:end , (1,2,N,N+1,N+2,2N+2)),setdiff(1:end , (1,2,N,N+1,N+2,2N+2))]
        A2=A2[setdiff(1:end , (1,2,N,N+1,N+2,2N+2)),setdiff(1:end , (1,2,N,N+1,N+2,2N+2))]
        A3=A3[setdiff(1:end , (1,2,N,N+1,N+2,2N+2)),setdiff(1:end , (1,2,N,N+1,N+2,2N+2))]
        A4=A4[setdiff(1:end , (1,2,N,N+1,N+2,2N+2)),setdiff(1:end , (1,2,N,N+1,N+2,2N+2))]
        A0=sparse(A0);
        A1=sparse(A1);
        A2=sparse(A2);
        A3=sparse(A3);
        A4=sparse(A4);
        return A0,A1,A2,A3,A4
    end
    function KEB_Parameter_Matrix(baseflow,N,ω,α,β,R,Ro,Co,z_po)
        θ = range(0,length=N+1,stop=pi)
        x = reshape(-cos.(θ), N+1, 1)
        c = [2; ones(N-1, 1) ; 2] .* (-1) .^ (0:N)
        X = repeat(x, 1, N+1);
        dX = X - X';
        D = (c * (1 ./ c)') ./ (dX .+ I(N+1));
        D = D - diagm(vec(sum(D, dims=2))); 
        #interpolation
        f = open( baseflow, "r" )
        n = countlines( f )
        seekstart( f )
        data = zeros(n,6)
        for i = 1:n
            z,w,u,v,du,dv = split( readline( f ), " " ) 
            data[i,1] = parse(Float64,z)
            data[i,2] = parse(Float64,w)
            data[i,3] = parse(Float64,u)
            data[i,4] = parse(Float64,v)  
            data[i,5] = parse(Float64,du)
            data[i,6] = parse(Float64,dv)
        end

        close( f )

        t=data[:,1]
        w=-1*data[:,2]
        u=-1*data[:,3]
        v=-1*data[:,4]
        du=-1*data[:,5]
        dv=-1*data[:,6]
        itpw=itp = interpolate(t, w , BSplineOrder(4))
        itpu=itp = interpolate(t, u , BSplineOrder(4))
        itpv=itp = interpolate(t, v , BSplineOrder(4))
        itpdu=itp = interpolate(t, du , BSplineOrder(4))
        itpdv=itp = interpolate(t, dv , BSplineOrder(4))
        for i=1:N+1
            D[i,:]=D[i,:].*((2*x[i]^3-x[i]^2+3*x[i]-4)^2/(20*(6*x[i]^2-2*x[i]+3)))
        end
        for i=1:N+1
            x[i]=(4*x[i]^3-2*x[i]^2+6*x[i]+12)/(-2*x[i]^3+x[i]^2-3*x[i]+4)
            if x[i]>50
                x[i]=50
            end
        end
        D2=D^2;
        U=zeros(N+1,1)
        V=zeros(N+1,1)
        W=zeros(N+1,1)
        dU=zeros(N+1,1)
        dV=zeros(N+1,1)
        dW=zeros(N+1,1)
        p=zeros(N+1,1)
        for i=1:N+1
            U[i,1]=itpu(x[i])
            V[i,1]=itpv(x[i])
            W[i,1]=itpw(x[i])
            dU[i,1]=itpdu(x[i])
            dV[i,1]=itpdv(x[i])
        end
        dW=-2*U
        ddU=D*dU;
        ddV=D*dV;
        if z_po==1
            U=U[N+1,1]
            V=V[N+1,1]
            W=W[N+1,1]
            dU=dU[N+1,1]
            dV=dV[N+1,1]
            dW=dW[N+1,1]
            ddU=ddU[N+1,1]
            ddV=ddV[N+1,1]
        end
        if z_po==0
            U=U[1,1]
            V=V[1,1]
            W=W[1,1]
            dU=dU[1,1]
            dV=dV[1,1]
            dW=dW[1,1]
            ddU=ddU[1,1]
            ddV=ddV[1,1]
        end
        α_bar=α-im*(1/R)
        λ=sqrt(α^2 +β^2)
        λ_bar=sqrt(α_bar^2 +β^2)
        A1=im*λ^2 *λ_bar^2 -R*λ_bar^2 *(α*U+β*V)-R*(α_bar*ddU+β*ddV)+im*Ro*λ_bar^2 *dW 
        A2=2*Ro*dV 
        A3=-im*R*(α*dV-β*dU)
        A4=-im*λ^2+R*(α*U+β*V)-im*Ro*U
        B1=im*Ro*λ_bar^2 *W 
        B2=2*Ro*V+Co
        B3=2*Ro*V+Co
        B4=-im*Ro*W 
        C1=-im*(λ^2 +λ_bar^2) +R*(α*U+β*V)-im*Ro*dW-im*Ro*U
        C2=0
        C3=C2
        C4=im
        F1=-im*Ro*W 
        F2=F3=F4=C2
        E1=im
        E2=E3=E4=C2
        dF=[-im*dW*Ro 0;0 0]
        A=[A1 A2;A3 A4]
        B=[B1 B2;B3 B4]
        C=[C1 C2;C3 C4]
        F=[F1 F2;F3 F4]
        E=[E1 E2;E3 E4]

        return A,B,C,F,E
    end
end