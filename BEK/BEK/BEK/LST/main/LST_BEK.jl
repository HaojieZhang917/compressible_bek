"""doc:
    2024.3.4
    ———————————————————————————————————————————————————————————————————————————————————————————————————
    | module_TimeMode is the function to caculate the time eigenvalues problem of the BEK model.      |
    |    The module include three subfunctions.                                                       |
    |    KEB_LST_all:consider the Coriolis force and curvature items which is the Practical problems. |
    |    KEB_LST_OS:neglect the Coriolis force and curvature.                                         |
    |    KEB_LST_Rayleigh:solve the inviscous Rayleigh equation.                                      |
    |input:                                                                                           |
    |    baseflow:velocity profiles of the BEK problem include u v w                                  |
    |    N:Number of discrete points                                                                  |
    |    α:radial wavenumber                                                                          |
    |    β:Circumferential wavenumber                                                                 |
    |    R:indimensional radial                                                                       |
    |    Ro:Rossby Number                                                                             |
    |output:                                                                                          |
    |    Aϕ=ωBϕ                                                                                       |
    |- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
    | module_SpatialMode is the function to caculate the spatial eigenvalues problem of BEK model     |
    |   The module include one subfunctions.                                                          |
    |   KEB_LST_ALL:consider the Coriolis force and curvature items which is the Practical problems.  |
    |   KEB_LST_OS_SPA:neglect the Coriolis force and curvature items
    |input:                                                                                           |
    |    baseflow:velocity profiles of the BEK problem include u v w                                  |
    |    N:Number of discrete points                                                                  |
    |    ω:frequecy of the disturbence                                                                |
    |    β:Circumferential wavenumber                                                                 |
    |    R:indimensional radial                                                                       |
    |    Ro:Rossby Number                                                                             |
    |output:                                                                                          |
    |   (A0 +A1*alpha +A2*alpha^2 +A3*alpha^3 +A4*alpha^4)ϕ=0                                                             |
    ————————————————————————————————————————————————————————————————————————————————————————————————————
    The equation in function is the incompressible N-S equation and simplified by vorticity methods
    Discrete method is chebyshev method 
    2024.3.14
    Add the function KEB_LST_OS_SPA
    """
module KEB_TimeMode

    using LinearAlgebra
    using BSplineKit

    export KEB_LST_all,KEB_LST_OS,KEB_LST_Rayleigh,rayleigh_quotient_iteration
    function  KEB_LST_all(baseflow,N,α,β,R,Ro)

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
        w=data[:,2]
        u=data[:,3]
        v=data[:,4]
        du=data[:,5]
        dv=data[:,6]
        itpw=itp = interpolate(t, w , BSplineOrder(4))
        itpu=itp = interpolate(t, u , BSplineOrder(4))
        itpv=itp = interpolate(t, v , BSplineOrder(4))
        itpdu=itp = interpolate(t, du , BSplineOrder(4))
        itpdv=itp = interpolate(t, dv , BSplineOrder(4))
        #interpolation
        # for i=1:N+1
        #     x[i,1]=10* x[i,1]+ 10
        # end
        # D=0.1* D
        for i=1:N+1
            D[i,:]=D[i,:].*((2*x[i]^3-x[i]^2+3*x[i]-4)^2/(20*(6*x[i]^2-2*x[i]+3)))
        end
        for i=1:N+1
            x[i]=(4*x[i]^3-2*x[i]^2+6*x[i]+12)/(-2*x[i]^3+x[i]^2-3*x[i]+4)
            if x[i]>20
                x[i]=20
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
        α_bar=α-im/R
        λ=sqrt(α^2+β^2)
        λ_bar=sqrt(α*α_bar+β^2)
        A11=im*((D2-(λ^2)*I(N+1))*(D2-(λ_bar^2)*I(N+1)))+R*((α*U+β*V)).*(D2-(λ_bar^2)*I(N+1))-R*(α_bar*ddU+β*ddV).*I(N+1)-im*Ro*W.*D*(D^2-(λ_bar^2)*I(N+1))-im*Ro*dW.*(D2-(λ_bar^2)*I(N+1))-im*Ro*U.*D2
        A12=(2*Ro*V.+ Co).*D+2*Ro*dV.*I(N+1)
        A21=(2*Ro*V.+ Co).*D-im*R*(α*dV-β*dU).*I(N+1)
        A22=im*(D2-(λ^2)*I(N+1))+R*(α*U+β*V).*I(N+1)-im*Ro* W.*D-im*Ro* U.*I(N+1)
        B11=R*(D2-(λ_bar^2)*I(N+1))
        B12=zeros(N+1,N+1)
        B21=B12
        B22=R*I(N+1)
        B22=Matrix(B22)
        #boundary condition
        A=[A11 A12;A21 A22]
        B=[B11 B12;B21 B22]
        A=A[setdiff(1:end , (1,2,N,N+1,N+2,2N+2)),setdiff(1:end , (1,2,N,N+1,N+2,2N+2))]
        B=B[setdiff(1:end , (1,2,N,N+1,N+2,2N+2)),setdiff(1:end , (1,2,N,N+1,N+2,2N+2))]
        return A,B
    end
    function KEB_LST_OS(baseflow,N,α,β_bar,R)
        θ = range(0,length=N+1,stop=pi)
        x = reshape(-cos.(θ), N+1, 1)
        c = [2; ones(N-1, 1) ; 2] .* (-1) .^ (0:N)
        X = repeat(x, 1, N+1);
        dX = X - X';
        D = (c * (1 ./ c)') ./ (dX .+ I(N+1));   # off-diagonal entries
        D = D - diagm(vec(sum(D, dims=2)));      # diagonal entries
        #interpolation
        f = open( baseflow, "r" )
        n = countlines( f )
        seekstart( f )
        data = zeros(n,6)
        for i = 1:n
            z,w,u,v,du,dv = split( readline( f ), " " )  # 读取每一行数据并用split函数将数据“剥离”开来
            data[i,1] = parse(Float64,z)
            data[i,2] = parse(Float64,w)
            data[i,3] = parse(Float64,u)
            data[i,4] = parse(Float64,v)  # 将字符串转化为数字
            data[i,5] = parse(Float64,du)
            data[i,6] = parse(Float64,dv)
        end

        close( f )

        t=data[:,1]
        w=data[:,2]
        u=data[:,3]
        v=data[:,4]
        du=data[:,5]
        dv=data[:,6]
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
            if x[i]>20
                x[i]=20
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
        λ=sqrt(α^2+β_bar^2)
        A=(D2-(λ^2)*I(N+1))^2-im*R*(α*U+β_bar*V).*(D2-(λ^2)*I(N+1))+im*R*(α*ddU+β_bar*ddV).*I(N+1)
        B=-1*im*R*α*(D2-(λ^2)*I(N+1))
        A=A[setdiff(1:end , (1,2,N,N+1,N+2,2N+2)),setdiff(1:end , (1,2,N,N+1,N+2,2N+2))]
        B=B[setdiff(1:end , (1,2,N,N+1,N+2,2N+2)),setdiff(1:end , (1,2,N,N+1,N+2,2N+2))]
        return A,B,x
    end
    function KEB_LST_Rayleigh(baseflow,N,α,β)



        θ = range(0,length=N+1,stop=pi)
        x = reshape(-cos.(θ), N+1, 1)
        c = [2; ones(N-1, 1) ; 2] .* (-1) .^ (0:N)
        X = repeat(x, 1, N+1);
        dX = X - X';
        D = (c * (1 ./ c)') ./ (dX .+ I(N+1));   # off-diagonal entries
        D = D - diagm(vec(sum(D, dims=2)));      # diagonal entries
        #interpolation
        f = open( baseflow, "r" )
        n = countlines( f )
        seekstart( f )
        data = zeros(n,6)
        for i = 1:n
            z,w,u,v,du,dv = split( readline( f ), " " )  # 读取每一行数据并用split函数将数据“剥离”开来
            data[i,1] = parse(Float64,z)
            data[i,2] = parse(Float64,w)
            data[i,3] = parse(Float64,u)
            data[i,4] = parse(Float64,v)  # 将字符串转化为数字
            data[i,5] = parse(Float64,du)
            data[i,6] = parse(Float64,dv)
        end

        close( f )

        t=data[:,1]
        w=data[:,2]
        u=data[:,3]
        v=data[:,4]
        du=data[:,5]
        dv=data[:,6]
        itpw=itp = interpolate(t, w , BSplineOrder(4))
        itpu=itp = interpolate(t, u , BSplineOrder(4))
        itpv=itp = interpolate(t, v , BSplineOrder(4))
        itpdu=itp = interpolate(t, du , BSplineOrder(4))
        itpdv=itp = interpolate(t, dv , BSplineOrder(4))
        psi=zeros(N+1,1)
        for i=1:N+1
            D[i,:]=D[i,:].*((2*x[i]^3-x[i]^2+3*x[i]-4)^2/(20*(6*x[i]^2-2*x[i]+3)))
        end
        for i=1:N+1
            x[i]=(4*x[i]^3-2*x[i]^2+6*x[i]+12)/(-2*x[i]^3+x[i]^2-3*x[i]+4)
            if x[i]>10
                x[i]=10
            end
        end
        D2=D^2;
        D4=D^4;
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
        for i=1:N
            x[i,1]=0.2*(x[i,1]-5)
        end
        dW=-2*U
        ddU=D*dU;
        ddV=D*dV;
        λ=sqrt(α^2+β^2)
        A=(α*U+β*V).*(D2-(λ^2)*I(N+1))-(α*ddU+β*ddV).*I(N+1)
        B=D2-(λ^2)*I(N+1)
        A=A[setdiff(1:end , (1,2,N,N+1,N+2,2N+2)),setdiff(1:end , (1,2,N,N+1,N+2,2N+2))]
        B=B[setdiff(1:end , (1,2,N,N+1,N+2,2N+2)),setdiff(1:end , (1,2,N,N+1,N+2,2N+2))]
        return A,B,x
    
    end
    end
module KEB_SpatialMode

    using LinearAlgebra
    using BSplineKit
    using SparseArrays
    export KEB_LST_ALL
    function KEB_LST_ALL(baseflow,N,ω,β,R,Ro,Co)
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
        w=data[:,2]
        u=data[:,3]
        v=data[:,4]
        du=data[:,5]
        dv=data[:,6]
        itpw=itp = interpolate(t, w , BSplineOrder(4))
        itpu=itp = interpolate(t, u , BSplineOrder(4))
        itpv=itp = interpolate(t, v , BSplineOrder(4))
        itpdu=itp = interpolate(t, du , BSplineOrder(4))
        itpdv=itp = interpolate(t, dv , BSplineOrder(4))
        #interpolation
        # for i=1:N+1
        #     x[i,1]=10* x[i,1]+ 10
        # end
        # D=0.1* D
        for i=1:N+1
            D[i,:]=D[i,:].*((2*x[i]^3-x[i]^2+3*x[i]-4)^2/(20*(6*x[i]^2-2*x[i]+3)))
        end
        for i=1:N+1
            x[i]=(4*x[i]^3-2*x[i]^2+6*x[i]+12)/(-2*x[i]^3+x[i]^2-3*x[i]+4)
            if x[i]>20
                x[i]=20
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
        A=D2-(β^2)*I(N+1)
        L0_1= im*A^2 + R*β*V.*A - R*ω.*A - R*β*ddV.*I(N+1) + im*ddU.*I(N+1) - Ro*im*W.*A*D - Ro*im*dW.*A - Ro*im*U.*D2
        L0_2=(2Ro*V.+Co).*D + 2*Ro*dV.*I(N+1)
        L0_3=(2Ro*V.+Co).*D+ im*R*β*dU.*I(N+1)
        L0_4= im*A + R*β*V.*I(N+1) - R*ω*I(N+1) - im*Ro*W.*D - im*Ro*U.*I(N+1)
        L1_1= -(1/R).*A + R*U.*A + im*β*V.*I(N+1) - im*ω*I(N+1) - R*ddU.*I(N+1) + Ro*(1/R)*W.*D + Ro*(1/R)*dW.*I(N+1)
        L1_2= zeros(N+1,N+1)
        L1_3= -1*im*R*dV.*I(N+1)
        L1_4= R*U.*I(N+1)
        L2_1= -2*im*A + ω*R*I(N+1) - β*R*V.*I(N+1) + im*U.*I(N+1) + Ro*im*W.*D + Ro*im*dW.*I(N+1)
        L2_2= L1_2
        L2_3= L1_2
        L2_4= -im*I(N+1)
        L3_1= (1/R)*I(N+1) - R*U.*I(N+1)
        L3_2= L3_3=L3_4=L1_2
        L4_1= im*I(N+1)
        L4_2= L4_3=L4_4=L1_2
        L2_4= Matrix{ComplexF64}(L2_4)
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
    function KEB_LST_OS_SPA(baseflow,N,ω,β,R)
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
        w=data[:,2]
        u=data[:,3]
        v=data[:,4]
        du=data[:,5]
        dv=data[:,6]
        itpw=itp = interpolate(t, w , BSplineOrder(4))
        itpu=itp = interpolate(t, u , BSplineOrder(4))
        itpv=itp = interpolate(t, v , BSplineOrder(4))
        itpdu=itp = interpolate(t, du , BSplineOrder(4))
        itpdv=itp = interpolate(t, dv , BSplineOrder(4))
        #interpolation
        # for i=1:N+1
        #     x[i,1]=10* x[i,1]+ 10
        # end
        # D=0.1* D
        for i=1:N+1
            D[i,:]=D[i,:].*((2*x[i]^3-x[i]^2+3*x[i]-4)^2/(20*(6*x[i]^2-2*x[i]+3)))
        end
        for i=1:N+1
            x[i]=(4*x[i]^3-2*x[i]^2+6*x[i]+12)/(-2*x[i]^3+x[i]^2-3*x[i]+4)
            if x[i]>20
                x[i]=20
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
        A=D2-(β^2)*I(N+1)
        A0=im*A^2+R*β*V.*A-R*ω.*A-R*β*ddV.*I(N+1)
        A1=R*U.*A-R*ddU.*I(N+1)
        A2=-2*im*A-R*β*V.*I(N+1)+R*ω.*I(N+1)
        A3=-R*U.*I(N+1)
        A4=im*I(N+1)
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
    end
module iteration
    using LinearAlgebra
    function rayleigh_quotient_iteration(A, B, sigma; q0=rand(size(A, 1), 1))

    flg = true
    while flg
        sigma0 = sigma[1]+ 0.0e0im
        q = (A - sigma*B) \ (B*q0)
        q0 = q/maximum(abs.(q))
        sigma = ((q0'*(A*q0))/(q0'*(B*q0)))[1]
        if abs(sigma-sigma0)<=eps(1.0f0)
            flg = false
        end

    end

      return sigma, q0
    end
    end
module io
    using DelimitedFiles
    using LinearAlgebra 
    function inp()
        print("initial,step,time")
        input=readline(stdin)
        inputstring=split(input, ",");
        initial=parse(Float64, inputstring[1]);
        step=parse(Float64, inputstring[2]);
        times=parse(Int64, inputstring[3]);
        Ro=zeros(times+1,1);
        for i=0:times
            Ro[i+1,1]=initial+step*i
        end
        return Ro
    end
    end