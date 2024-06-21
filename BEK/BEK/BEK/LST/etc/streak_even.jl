module op

    import SpecialFunctions
    import QuadGK
    import SparseArrays
    import LinearAlgebra
    import FFTW
    # using DifferentialEquations
    using SparseArrays
    using DelimitedFiles
    using SpecialFunctions
    using QuadGK
    using SparseArrays
    using Printf

    module constvar

        import FFTW
        const y0=0.0;
        const Lz=2*pi;

        # Set parameters
        const eps=0.0;
        const nz=32;
        const z=[i for i=0:nz-1]./nz*Lz.-pi;
        const beta=2*pi/Lz;

        # const Ma=0.0;
        const lambda_B=1.0;
        #const lambda_u=sin.(z).*As.+lambda_B;

        const one=ones(nz);
        const N=div(nz,2)-1;

        #get lambda_uz
        #maybe wrong

        function baseflow(Amp)
            As=complex(zeros(nz));
            As=As_distribute(As,"cos_simple1")

            As=As.*Amp
            #@show As
            lambda_u=FFTW.ifft(As)*nz.+lambda_B;
            spt_lambda_u=im*beta.*[0; [i for i=1:div(nz,2)-1]; 0; [i for i=-div(nz,2)+1:-1]].*FFTW.fft(lambda_u)/nz;
            lambda_uz=FFTW.ifft(spt_lambda_u).*nz;
            return (lambda_u, lambda_uz)
        end

        function As_distribute(As, flag)

            if flag=="cos_simple1"
                As[2]=-0.5e0;
                As[end]=conj(As[2]);
            elseif flag=="cos_simple4"
                As[5]=-0.5e0;
                As[end-3]=conj(As[5]);
            end
            return As;
        end

    end

    using .constvar

    eta(y0, omega, alpha, lambda_u_)=
        eta=(im*alpha.*lambda_u_).^(1.0/3.0).*(constvar.one*y0.-constvar.one*omega./(alpha.*lambda_u_));

    gamma(alpha, n, Ma)=sqrt(alpha^2*(1-Ma^2)
    +constvar.beta^2*(n+constvar.eps)^2);
    g(eta0, kappa, omega, alpha)=
    constvar.one*1.5.+eta0./(2.0*airyai.(eta0)).*(eta0.*kappa.+airyaiprime.(eta0));
    R(g, lambda_u_, lambda_uz_)=
        FFTW.fft(lambda_uz_./lambda_u_.*g)/constvar.nz;
    Q(eta0, kappa, alpha, lambda_u_)=
    FFTW.fft((im*alpha*lambda_u_).^(5//3).*airyaiprime.(eta0)./kappa)/constvar.nz;
    df1(f1, f2, x1, x2)=(f2-f1)/(x2-x1);
    df2(f1, f2, f3, x1, x2, x3)=(f3-f2)/((x3-x2)*(x3-x1))-(f2-f1)/((x2-x1)*(x3-x1));

    get((a, b))=a;

    function initial()
        FFTW.set_num_threads(2)
    end

    function getidx(n, j)
        idx=n-j;
        if idx>=0
            idx=idx+1;
        else
            idx=constvar.nz+1+idx;
        end
        return idx;
    end

    function M1(alpha)
        v=alpha^2*ones(constvar.N+1).+constvar.beta^2*([i for i=0:constvar.N].+ones(constvar.N+1)*constvar.eps).^2
        i_idx=[i for i=1:constvar.N+1];
        return sparse(i_idx, i_idx, v);
    end

    function M2(R, Q, omega, alpha, Ma)
        M2=complex(zeros(constvar.N+1, constvar.N+1));

        for n=0:constvar.N
            M2[n+1, 1]=Q[n+1]/alpha^2*gamma(alpha, 0, Ma);
            for j=1:constvar.N
                if n-j>=-constvar.N && n-j<=constvar.N
                    idx1=getidx(n, j);
                    idx2=getidx(n, -j);
                    M2[n+1, j+1]=
                    (R[idx1]-R[idx2])*im*(constvar.eps+j)*constvar.beta+(Q[idx1]+Q[idx2])/alpha^2*gamma(alpha, j, Ma);
                end
            end
        end

        return M2
    end

    function detM(omega, alpha, lambda_u_, lambda_uz_, Ma)
        eta0=eta(0.0, omega, alpha, lambda_u_);
        kappa=interToInf(eta0, omega, alpha, lambda_u_);
        varg=g(eta0, kappa, omega, alpha);
        varR=R(varg, lambda_u_, lambda_uz_);
        varQ=Q(eta0, kappa, alpha, lambda_u_);
        matM1=M1(alpha);
        matM2=M2(varR, varQ, omega, alpha, Ma);
        M=matM1.+matM2;
        LinearAlgebra.det(M);
        M[div(constvar.nz,2), div(constvar.nz,2)];
    end

    function funlocal(omega, alpha, beta, lambda_u_, lambda_uz_, Ma)
        eta0=eta(0.0, omega, alpha, lambda_u_);
        kappa=interToInf(eta0, omega, alpha, lambda_u_);
        varg=g(eta0, kappa, omega, alpha);
        varR=R(varg, lambda_u_, lambda_uz_);
        varQ=Q(eta0, kappa, alpha, lambda_u_);
        # matM1=M1(alpha);
        # matM2=M2(varR, varQ, omega, alpha, Ma);
        # M=matM1.+matM2;
        # LinearAlgebra.det(M);
        f=alpha^2*(alpha^2+beta^2)+varQ[1]*sqrt(alpha^2*(1-Ma^2)+beta^2);
        # M[div(constvar.nz,2), div(constvar.nz,2)];
    end

    function funlocal_matrix(omega, alpha, lambda_u_, lambda_uz_, Ma)
        # eta0=eta(0.0, omega, alpha, lambda_u_);
        # kappa=interToInf(eta0, omega, alpha, lambda_u_);
        # varg=g(eta0, kappa, omega, alpha);
        # varR=R(varg, lambda_u_, lambda_uz_);
        # varQ=Q(eta0, kappa, alpha, lambda_u_);
        # matM1=M1(alpha);
        # matM2=M2(varR, varQ, omega, alpha, Ma);
        # M=matM1.+matM2;
        M=constructM(omega, alpha, lambda_u_, lambda_uz_, Ma);
        # @show LinearAlgebra.diag(M)
        # f=LinearAlgebra.det(M);
        # @show f
        F=LinearAlgebra.lu(M);
        t=LinearAlgebra.diag(F.U);
        # tt,t=findmin(abs.(p));
        # @show p[t];
        # @show sign(real(f))*abs(real(p[t]))+sign(imag(f))*abs(imag(p[t]))*im
        p=1.0;
        for i=1:15
            p=p*t[i];
        end
        # @show p
        #@show p, abs(LinearAlgebra.det(M))

        return p


        # f=alpha^2*(alpha^2+beta^2)+varQ[1]*sqrt(alpha^2*(1-Ma^2)+beta^2);
        # M[div(constvar.nz,2), div(constvar.nz,2)];
    end
    function gauss_elimination1(A, b)
        (row, col)=size(A)
        n=col;
        for k=1:n-1
            for i=k+1:n
                factor=A[i, k]/A[k, k]
                for j=k+1:n
                    A[i, j]=A[i, j]-factor*A[k, j]
                end
                b[i]=b[i]-factor*b[k]
            end
        end
        return (A, b)
    end

    function gauss_elimination2(A, b)
        (row, col)=size(A)
        n=col;
            x=zeros(n, 1).+zeros(n, 1).*im
            x[n]=1.0
            for i=n-1:-1:1
                sum=b[i]
                for j=i+1: n
                    sum=sum-A[i, j]*x[j]
                end
                x[i]=sum/A[i, i]
            end
            return x[16:-1:1]
    end



    function writeeigen_matrix(omega, alpha, lambda_u_, lambda_uz_, Ma)
        # eta0=eta(0.0, omega, alpha, lambda_u_);
        # kappa=interToInf(eta0, omega, alpha, lambda_u_);
        # varg=g(eta0, kappa, omega, alpha);
        # varR=R(varg, lambda_u_, lambda_uz_);
        # varQ=Q(eta0, kappa, alpha, lambda_u_);
        # matM1=M1(alpha);
        # matM2=M2(varR, varQ, omega, alpha, Ma);
        # M=matM1.+matM2;
        M=constructM(omega, alpha, lambda_u_, lambda_uz_, Ma);
        @show size(M)
        # @show LinearAlgebra.diag(M)
        # f=LinearAlgebra.det(M);
        # @show f
        # F=LinearAlgebra.lu(M);
        # t=LinearAlgebra.diag(F.U);
        # @show t
        # p=1.0;
        # for i=1:15
        #     p=p*t[i];
        # end

        (A1, b1)=gauss_elimination1(M[16:-1:1, 16:-1:1], zeros(16, 1))
        @show A1, b1
        x=gauss_elimination2(A1, b1)
        # @show p
        #@show p, abs(LinearAlgebra.det(M))
        @show x
        y=[x[1:16];0.0;x[16:-1:2]]
        return y
    end


    function constructM(omega, alpha, lambda_u_, lambda_uz_, Ma)
        eta0=eta(0.0, omega, alpha, lambda_u_);
        kappa=interToInf(eta0, omega, alpha, lambda_u_);
        varg=g(eta0, kappa, omega, alpha);
        varR=R(varg, lambda_u_, lambda_uz_);
        varQ=Q(eta0, kappa, alpha, lambda_u_);
        matM1=M1(alpha);
        matM2=M2(varR, varQ, omega, alpha, Ma);
        M=matM1.+matM2;
        # @show size(M)
        return M;

#        LinearAlgebra.det(M);
       # return [M[div(constvar.nz,2), div(constvar.nz,2)]];
    end

    function interToInf(eta1::Array, omega, alpha, lambda_u_)

        dh=2.0;
        y2=40.0;
        eta2=eta(y2, omega, alpha, lambda_u_);

        a=(eta2-eta1)./(y2-constvar.y0);
        b=eta1.-a.*constvar.y0;
        x2=y2; etah=a.*x2.+b;
        int1=quadgk.(airyai, eta1, etah);int1=get.(int1);
        err=ones(constvar.nz);
        while maximum(abs.(err))>1.0e-8
            x2=x2*dh; etah=a*x2+b;
            int2=get.(quadgk.(airyai, eta1, etah));
            err=int2.-int1;
            int1=int2;
        end
        return int1;
    end

    ##################################################
    # julia code implementing Beyn's contour-integral
    # method for nonlinear eigenproblems
    #
    # homer reid 10/2017
    ##################################################

    ##################################################
    # Implement Beyn's method to compute all eigenvalues
    # of a DxD matrix lying within an elliptical contour
    # centered at z0, with horizontal/vertical radius Rx/Ry,
    # using N rectangular-rule quadrature points, assuming
    # the contour contains no more than L eigenvalues.
    #
    # The matrix is only ever referenced abstractly through
    # the user-supplied UserBeynFunction, which has prototype
    #
    #  UserBeynFunction(z, VHat, MInvVHat)
    #
    # and which should set MInvVHat = M(z) \ VHat and
    # return MInvVHat.
    #
    # BeynSolve returns a tuple (Lambda,V)
    # where Lambda[k] = kth eigenvalue,
    #       V[:,k]    = kth eigenvector
    ##################################################
    function BeynSolve(z0, Rx, Ry, UserBeynFunction, D; L=10, N=50)

      MInvVHat   = im*zeros(D,L);
      A0         = im*zeros(D,L);
      A1         = im*zeros(D,L);
      VHat       = randn(D,L) + im*randn(D,L);

      #################################################################
      # evaluate contour integrals
      #################################################################
      DeltaTheta = 2.0*pi/N;
      for n=0:N-1
        Theta      = n*DeltaTheta;
        CT         = cos(Theta);
        ST         = sin(Theta);
        z1         = Rx*CT + im*Ry*ST;
        dz         = (im*Rx*ST + Ry*CT)/N;
        MInvVHat   = UserBeynFunction(z0+z1, VHat, MInvVHat);
        A0        += dz*MInvVHat;
        A1        += z1*dz*MInvVHat;
      end

      #################################################################
      # linear-algebra postprocessing to extract eigenvalues/vectors
      #################################################################
      (V,Sigma,W)=LinearAlgebra.svd(A0);
      # @show size(A0)
      # @show size(V)
      # @show size(W)

      SigmaThreshold=1.0e-8;
      logical=abs.(Sigma) .> SigmaThreshold
      K=findlast(logical)
      if K!=nothing
          @printf("Found %i nonzero singular values.\n",K);
      else
          K=0;
          @printf("Found no nonzero singular values.\n");
      end

      if (K==0)
        Lambda=V=im*zeros(0,0);
        return (Lambda,V)
      end

      if (K==L)
        @printf("** Warning: K=L=%i in BeynMethod (repeat with higher L)",K);
      end

      Vk = V[:,1:K];
      Sk = Sigma[1:K];
      Wk = W[:,1:K];
      B  = Vk' * A1 * Wk * LinearAlgebra.diagm( 0=>1.0 ./ Sk );

      #(Lambda,V)=LinearAlgebra.eig(B)
      Lambda=LinearAlgebra.eigvals(B);
      order=sortperm(Lambda, by=abs);

      V=LinearAlgebra.eigvecs(B);
      # @show Lambda[order]
      # @show order
      # @show size(V[:, 1]), V[:, 1]
      # @show size(V[:, 2]), V[:, 2]


      Lambda += z0*ones(size(Lambda));
      Lambda=Lambda[order];
      V=Vk*V[:, order];


      return (Lambda, V);

    end

    ##################################################
    # alternative entry point to BeynSolve for circular contour
    # of radius R centered at z0
    ##################################################
    function BeynSolve(z0, R, UserBeynFunction, D; L=10, N=50)
      return BeynSolve(z0, R, R, UserBeynFunction, D; L=L, N=N)
    end

    function callinterface(z, VHat, MInvVHat, omega, lambda_u_, lambda_uz_, Ma)

        M=constructM(omega, z, lambda_u_, lambda_uz_, Ma);
        MInvVHat=M\VHat;
        MInvVHat;

    end

    function mull(omega, alpha, lambda_u_, lambda_uz_, Ma; R=1.0e0)

        L=10;
        N=100;
        D=constvar.N+1;
        # R=0.01*abs(alpha)
        R=abs(imag(alpha))
        @show R
        (lambda, V)=BeynSolve(alpha, R, R, (z,V,MIV)->callinterface(z,V,MIV, omega, lambda_u_, lambda_uz_, Ma), D; L=L, N=N)

#       (lambda, V)=BeynSolve(alpha, 0.4*abs(real(alpha)), 2.0*abs(imag(alpha)), (z,V,MIV)->callinterface(z,V,MIV, omega, lambda_u_, lambda_uz_, Ma), D; L=L, N=N)
       while true
           if length(lambda)>0
               break;
           else
               R=R*2.0;
               (lambda, V)=BeynSolve(alpha, R, R, (z,V,MIV)->callinterface(z,V,MIV, omega, lambda_u_, lambda_uz_, Ma), D; L=L, N=N)
           end
       end
       @show V[1:3]
       # (lambda, V)=BeynSolve(alpha, R, R, (z,V,MIV)->callinterface(z,V,MIV, omega, lambda_u_, lambda_uz_, Ma), D; L=L, N=N)
       err=Array{Float64}(undef, length(lambda))
           for i=1:length(lambda)
               M=constructM(omega, lambda[i], lambda_u_, lambda_uz_, Ma)
               # D=LinearAlgebra.diag(M);
               # @show D
               # err[i]=maxval(LinearAlgebra.norm(M*V[:, 1]);)
               err[i]=LinearAlgebra.norm(M*V[:, i]);

           end
           @show err
           @show lambda
        #return (lambda[1], V[:, 1])
        return (lambda[1], err[1], V[:, 1])
        #@show V
    end

    # function mull(omega, alpha, lambda_u_, lambda_uz_, Ma)
    #
    #     R=
    #
    # end

    function mull_old(omega, alpha, beta, lambda_u_, lambda_uz_, Ma)

        d_alf=1e-4;

        alpha1=alpha+d_alf*im;
        alpha2=alpha;
        alpha3=alpha-d_alf*im;

        f1=funlocal(omega, alpha1, beta, lambda_u_, lambda_uz_, Ma);
        f2=funlocal(omega, alpha2, beta, lambda_u_, lambda_uz_, Ma);
        f3=funlocal(omega, alpha3, beta, lambda_u_, lambda_uz_, Ma);

        error=1.0;
        while error>1.0e-8
            d12=df1(f1, f2, alpha1, alpha2);
            d13=df1(f1, f3, alpha1, alpha3);
            d23=df1(f2, f3, alpha2, alpha3);
            w=d12+d13-d23;
            delta=sqrt(w^2-4.0*f1*df2(f1, f2, f3, alpha1, alpha2, alpha3));
            p1=w+delta;
            p2=w-delta;
            if abs(p1)>=abs(p2)
             alpha=alpha1-2.0*f1/p1;
            else
             alpha=alpha1-2.0*f1/p2;
            end;
            f3=f2; f2=f1;
            alpha3=alpha2; alpha2=alpha1;
            alpha1=alpha;
            f1=funlocal(omega, alpha1, beta, lambda_u_, lambda_uz_, Ma);
            error=max(abs(alpha1-alpha2), abs(alpha2-alpha3), abs(alpha1-alpha3));
#            println(omega, alpha);
        end
        # @show alpha
        alf=alpha;
    end

    function mull_matrix_old(omega, alpha, lambda_u_, lambda_uz_, Ma)

        d_alf=1e-8;

        alpha1=alpha+d_alf*im;
        alpha2=alpha;
        alpha3=alpha-d_alf*im;

        f1=funlocal_matrix(omega, alpha1, lambda_u_, lambda_uz_, Ma);
        f2=funlocal_matrix(omega, alpha2, lambda_u_, lambda_uz_, Ma);
        f3=funlocal_matrix(omega, alpha3, lambda_u_, lambda_uz_, Ma);

        error=1.0;
        while error>1.0e-10
            d12=df1(f1, f2, alpha1, alpha2);
            d13=df1(f1, f3, alpha1, alpha3);
            d23=df1(f2, f3, alpha2, alpha3);
            w=d12+d13-d23;
            delta=sqrt(w^2-4.0*f1*df2(f1, f2, f3, alpha1, alpha2, alpha3));
            p1=w+delta;
            p2=w-delta;
            if abs(p1)>=abs(p2)
             alpha=alpha1-2.0*f1/p1;
            else
             alpha=alpha1-2.0*f1/p2;
            end;
            f3=f2; f2=f1;
            alpha3=alpha2; alpha2=alpha1;
            alpha1=alpha;
            f1=funlocal_matrix(omega, alpha1, lambda_u_, lambda_uz_, Ma);
            @show f1
            error=max(abs(alpha1-alpha2), abs(alpha2-alpha3), abs(alpha1-alpha3));
#            println(omega, alpha);
        end
        # @show alpha
        alf, err=alpha, error;
    end

    function get_w(omega, alpha, lambda_u_, y, p)

        n=div(constvar.nz,2)-1;
        nz=[i for i in -n:n];
        pz=im.*(nz.+constvar.eps).*constvar.beta.*p;
        pz=[pz[n+1:end]; 0.0e0; pz[n:-1:1]];
        pz=FFTW.ifft(pz).*constvar.nz;
        tspan=(y[1], y[end]);
        dy=y[2]-y[1];

        lv=1.0/dy^2;
        rv=1.0/dy^2;
        ny,=size(y);
        w1=zeros(complex(Float64), ny, constvar.nz)

        for nspan=1:constvar.nz
            cv=-2.0/dy^2*ones(ny)+im*(omega.-alpha.*lambda_u_[nspan].*y);
            cv[1]=1.0e0; cv[end]=1.0e0;
            mlv=lv.*ones(ny-1);
            mrv=rv.*ones(ny-1);
            MA=SparseArrays.spdiagm(-1=>mlv, 1=>mrv, 0=>cv);
            b=ones(ny).*pz[nspan];
            b[1]=0.0e0; b[end]=-1.e0/(im*alpha*lambda_u_[nspan]*y[end])*pz[nspan];
            w1[:, nspan]=MA\b;
        end
        # for nspan in 1:1 #constvar.nz
        #     @show pz[nspan], lambda_u_[nspan]
        #     bvp2 = TwoPointBVProblem((du, u, p, t)->pendulum!(du, u, p, t, omega, alpha, lambda_u_[nspan], pz[nspan]),
        #     (residual, u, p, t)->bc2!(residual, u, p, t, alpha, lambda_u_[nspan], y, pz[nspan]),
        #     [0.0e0,0.0e0], tspan);
        #     sol2 = solve(bvp2, RK4(), dt=0.01);
        #
        # end


        return w1
    end

    function plot_inner(p)
        @show p
        z=FFTW.ifft(p)
        @show z
        return z
    end

    function plot_outter(p, alpha, M)
        @show p, size(p)
        pp=zeros(32, 200).+zeros(32, 200)*im
        p01=zeros(2, 200).+zeros(2, 200)*im

        dy=0.1
        for j=1:200
            for i=0:15
                k=sqrt(alpha^2*(1-M^2)+1^2*(constvar.eps+i)^2)
                @show i, k
                pp[i+1, j]=p[i+1]*exp(-k*(j-1)*dy)
            end
            pp[17, j]=0.0e0
            for i=2:16
                k=sqrt(alpha^2*(1-M^2)+1^2*(constvar.eps+i-1)^2)
                pp[34-i, j]=-p[i]*exp(-k*(j-1)*dy)
            end

            if j==200
            @show j, pp[1:32, j]
            end
            p01[1:2, j]=pp[1:2, j]
            pp[1:32, j]=FFTW.ifft(pp[1:32, j])
        end
        return pp,p01

    end
    #
    # function pendulum!(du, u, p, t, omega, alpha, lambda_u, pz)
    #     w1=u[1];
    #     u0=u[2];
    #     du[1]=u0;
    #     du[2]=pz-im*(omega-alpha*lambda_u*t)*w1;
    # end
    #
    # function bc2!(residual, u, p, t, alpha, lambda_u, y, pz)
    #     residual[1]=u[1][1]+0.0e0;
    #     residual[2]=u[end][1]+1.0/(im*alpha*lambda_u*y[end])*pz;
    # end



end

function backwardsearch(omega_start, alpha, omega_end, domega, As, Ma)
    omega_group=[];
    alpha_group=[];
    lambda_u_=copy(As);
    lambda_uz_=copy(As);
    (lambda_u_, lambda_uz_)=op.constvar.baseflow(As);
    for omega in (omega_start:domega:omega_end)
        alpha,=op.mull(omega, alpha, lambda_u_, lambda_uz_, Ma);
        println((omega, alpha));
        omega_group=[omega_group; omega];
        alpha_group=[alpha_group; alpha];
    end
    return (omega_group, alpha_group);
end

function forwardsearch(omega_start, alpha, omega_end, domega, As, Ma)
    lambda_u_=copy(As);
    lambda_uz_=copy(As);
    (lambda_u_, lambda_uz_)=op.constvar.baseflow(As);

    omega_group=[];
    alpha_group=[];
    for omega in (omega_start:domega:omega_end)
        alpha,=op.mull(omega, alpha, lambda_u_, lambda_uz_, Ma);
        println((omega, alpha));
        omega_group=[omega;omega_group];
        alpha_group=[alpha;alpha_group];
    end
    return (omega_group, alpha_group);
end

function search(omega, alpha, As, Ma)
    lambda_u_=copy(As);
    lambda_uz_=copy(As);
    (lambda_u_, lambda_uz_)=op.constvar.baseflow(As)
    for alphalist=3.0:-0.1:0.0
        alpha=alphalist
        alpha=op.mull(omega, alpha, lambda_u_, lambda_uz_, Ma)
        println((omega, alpha))
        exit()
    end
end
function search_max_omega(omega, alpha, lambda_u_, lambda_uz_, Ma, domega)

    alpha1,=op.mull(omega, alpha, lambda_u_, lambda_uz_, Ma);
    alpha2,=op.mull(omega+domega, alpha, lambda_u_, lambda_uz_, Ma);
    alpha3,=op.mull(omega-domega, alpha, lambda_u_, lambda_uz_, Ma);


    while imag(alpha1-alpha2)*imag(alpha1-alpha3)<0.0
        if imag(alpha1)<imag(alpha2) && imag(alpha3)<imag(alpha1)
            @show omega, domega
            omega=omega-domega;
            alpha2=alpha1;
            alpha1=alpha3;
            alpha3,=op.mull(omega-domega, alpha, lambda_u_, lambda_uz_, Ma);
            @show (omega, alpha1)
            @show (omega+domega, alpha2)
            @show (omega-domega, alpha3)
        elseif imag(alpha1)>imag(alpha2) && imag(alpha3)>imag(alpha1)
            omega=omega+domega;
            alpha3=alpha1;
            alpha1=alpha2;
            alpha2,=op.mull(omega+domega, alpha, lambda_u_, lambda_uz_, Ma);
            @show (omega, alpha1)
            @show (omega+domega, alpha2)
            @show (omega-domega, alpha3)
        end
    end

    return (omega, alpha1)

end



function backsearch_with_Ma(omega, alpha, As, Mastart, Maend, dMa)

    lambda_u_=copy(As);
    lambda_uz_=copy(As);
    (lambda_u_, lambda_uz_)=op.constvar.baseflow(As);

    Ma_group=[];
    omega_group=[];
    alpha_group=[];
    Maseries=[Ma for Ma in 0.0:0.1:0.6];
    Maseries=[Maseries; [Ma for Ma in 0.65:0.05:1.0]];
    Maseries=[Maseries; [Ma for Ma in 1.02:0.02:3.2]];
    Maseries=[Maseries; [Ma for Ma in 1.22:0.02:3.2]];
    Maseries=[Maseries; [Ma for Ma in 3.22:0.01:4.5]];
    # Maseries=[];
    Maseries=[Maseries; [Ma for Ma in 3.450:0.01:4.5]];
    Maseries=[Maseries; [Ma for Ma in 4.505:0.005:6.0]];
    io=open("Ma_info.txt", "a");
    for Ma in Maseries
        #alpha,=op.mull(omega, alpha, lambda_u_, lambda_uz_, Ma);
        (omega, alpha)=search_max_omega(omega, alpha, lambda_u_, lambda_uz_, Ma, 0.001*omega)

        println((Ma, omega, alpha));
        omega_group=[omega_group; omega]
        Ma_group=[Ma_group; Ma];
        alpha_group=[alpha_group; alpha];
        writedlm(io, [Ma real(alpha) imag.(alpha) real.(omega)]);
        flush(io)

    end
    close(io)
    return (Ma_group, omega_group, alpha_group);


end

function betasearch_old(omega, alpha, beta, As, Ma)
    (lambda_u_, lambda_uz_)=op.constvar.baseflow(As)
    #(alpha, p)=op.mull(omega, alpha, lambda_u_, lambda_uz_, Ma)
    alphawatch=[];
    alphaold=alpha;
    for beta in 3.0:0.5:15.5
        alpha=op.mull_old(omega, alpha, beta, lambda_u_, lambda_uz_, Ma)
        @show (beta,alpha)
        alphawatch=[alphawatch;alpha];
    end
    alpha=alphaold;
    for beta in 2.5:-0.5:0.5
        alpha=op.mull_old(omega, alpha, beta, lambda_u_, lambda_uz_, Ma)
        @show (beta,alpha)
        alphawatch=[alpha;alphawatch];
    end
    return alphawatch
end
function masearch_old(omega, alpha, beta, As, Ma)
    (lambda_u_, lambda_uz_)=op.constvar.baseflow(As)
    #(alpha, p)=op.mull(omega, alpha, lambda_u_, lambda_uz_, Ma)
    alphawatch=[];
    for Ma in 0.0e0:0.01e0:1.5e0
    alpha=op.mull_old(omega, alpha, beta, lambda_u_, lambda_uz_, Ma)
    @show (Ma,alpha)
    alphawatch=[alphawatch;alpha];
    end
    return alphawatch
end

import .others
import DelimitedFiles

function search_by_input(omega, alpha, beta, As, Ma)

    io=open("output.plt","a")
    (lambda_u_, lambda_uz_)=op.constvar.baseflow(As)
    # @show lambda_u_
    # @show omega, alpha, beta, Ma

    alpha=op.mull_old(omega, alpha, beta, lambda_u_, lambda_uz_, Ma)
    @show omega, alpha, beta, Ma
    DelimitedFiles.writedlm(io, [real(omega) imag(omega) real(alpha) imag(alpha) real(beta) imag(beta) Ma]);
    flush(io);
    inputstring=["0"; "0"; "0.001"];
    inputno=[0, 0, 0.00];
    while true
        inputstring=split(others.input("k, ki, value="), ",");
        inputno[1]=parse(Int64, inputstring[1]);
        inputno[2]=parse(Int64, inputstring[2]);
        inputno[3]=parse(Float64, inputstring[3]);

        k=inputno[1]; ki=inputno[2]; delta=inputno[3];
        if k==3
            for i=1:ki
                beta=beta+delta;
                alpha=op.mull_old(omega, alpha, beta, lambda_u_, lambda_uz_, Ma)
                @show omega, beta, Ma
                @show alpha
                DelimitedFiles.writedlm(io, [real(omega) imag(omega) real(alpha) imag(alpha) real(beta) imag(beta) Ma]);
                flush(io);
            end
        elseif k==4
            for i=1:ki
                omega=omega+delta;
                alpha=op.mull_old(omega, alpha, beta, lambda_u_, lambda_uz_, Ma)
                @show omega, beta, Ma
                @show alpha
                DelimitedFiles.writedlm(io, [real(omega) imag(omega) real(alpha) imag(alpha) real(beta) imag(beta) Ma]);
                flush(io);
            end
        elseif k==5
            for i=1:ki
                Ma=Ma+delta;
                alpha=op.mull_old(omega, alpha, beta, lambda_u_, lambda_uz_, Ma)
                @show omega, beta, Ma
                @show alpha
                DelimitedFiles.writedlm(io, [real(omega) imag(omega) real(alpha) imag(alpha) real(beta) imag(beta) Ma]);
                flush(io);
            end
        elseif k==-1
            break;
        else
            println("Please input a correct value!")
            continue;
        end
    end

    close(io)


end

function search_by_auto_new(omega, alpha, beta, As_end, Ma, filename)

    fn="omega"*string(omega)*"_"*filename
    @show fn
    io=open(fn,"a")
    As=0.0e0
    (lambda_u_, lambda_uz_)=op.constvar.baseflow(As)
    # @show lambda_u_
    # @show omega, alpha, beta, Ma

    alpha=op.mull_old(omega, alpha, beta, lambda_u_, lambda_uz_, Ma)
    @show omega, alpha, beta, Ma
    alpha_old=alpha;
    err=0.0e0;
    DelimitedFiles.writedlm(io, [As real(omega) imag(omega) real(alpha) imag(alpha) real(beta) imag(beta) Ma err]);
    flush(io);

    while omega<10.0e0

        while As<=As_end
            As=As+0.02e0;
            (lambda_u_, lambda_uz_)=op.constvar.baseflow(As);
            alpha,err,=op.mull_matrix_old(omega, alpha, lambda_u_, lambda_uz_, Ma)
            @show omega, Ma, As
            @show alpha
            if abs(As-0.3e0)<=1e-5
                DelimitedFiles.writedlm(io, [As real(omega) imag(omega) real(alpha) imag(alpha) real(beta) imag(beta) Ma err]);
                flush(io);
            elseif abs(As-0.6e0)<=1e-5
                DelimitedFiles.writedlm(io, [As real(omega) imag(omega) real(alpha) imag(alpha) real(beta) imag(beta) Ma err]);
                flush(io);
            end
        end
        omega=omega+0.05e0;
        alpha=alpha_old;
        As=0.0e0;
        (lambda_u_, lambda_uz_)=op.constvar.baseflow(As);
        # @show lambda_u_
        # @show omega, alpha, beta, Ma

        alpha=op.mull_old(omega, alpha, beta, lambda_u_, lambda_uz_, Ma)
        @show omega, alpha, beta, Ma
        err=0.0e0;
        alpha_old=alpha;

        DelimitedFiles.writedlm(io, [As real(omega) imag(omega) real(alpha) imag(alpha) real(beta) imag(beta) Ma err]);
        flush(io);
    end



end

function search_by_input_new(omega, alpha, As, Ma, filename)

    fn="omega"*string(omega)*"_"*filename
    @show fn
    io=open(fn,"a")
    (lambda_u_, lambda_uz_)=op.constvar.baseflow(As)
    # @show lambda_u_
    # @show omega, alpha, beta, Ma

    # alpha,err,=op.mull(omega, alpha, lambda_u_, lambda_uz_, Ma)
    alpha, err=op.mull_matrix_old(omega, alpha, lambda_u_, lambda_uz_, Ma)
    @show omega, Ma, As
    @show alpha
    DelimitedFiles.writedlm(io, [As real(omega) imag(omega) real(alpha) imag(alpha) real(beta) imag(beta) Ma err]);
    flush(io);
    inputstring=["0"; "0"; "0.001"];
    inputno=[0, 0, 0.00];
    while true
        inputstring=split(others.input("k, ki, value="), ",");
        inputno[1]=parse(Int64, inputstring[1]);
        inputno[2]=parse(Int64, inputstring[2]);
        inputno[3]=parse(Float64, inputstring[3]);

        k=inputno[1]; ki=inputno[2]; delta=inputno[3];
        if k==1
            for i=1:ki
                As=As+delta;
                (lambda_u_, lambda_uz_)=op.constvar.baseflow(As);
                alpha,err,=op.mull_matrix_old(omega, alpha, lambda_u_, lambda_uz_, Ma)
                @show omega, Ma, As
                @show alpha
                DelimitedFiles.writedlm(io, [As real(omega) imag(omega) real(alpha) imag(alpha) real(beta) imag(beta) Ma err]);
                flush(io);
            end
        elseif k==4
            for i=1:ki

                omega=omega+delta;
                alpha,err,=op.mull_matrix_old(omega, alpha, lambda_u_, lambda_uz_, Ma)

                @show omega, Ma, As
                @show alpha
                DelimitedFiles.writedlm(io, [As real(omega) imag(omega) real(alpha) imag(alpha) real(beta) imag(beta) Ma err]);
                flush(io);
            end
        elseif k==41
            for i=1:ki

                omega=omega+sign(delta)*max(0.5*abs(imag(alpha)),abs(delta));
                alpha,err,=op.mull(omega, alpha, lambda_u_, lambda_uz_, Ma)
                @show sign(delta)*max(0.5*abs(imag(alpha)),abs(delta))
                @show omega, Ma, As
                @show alpha
                DelimitedFiles.writedlm(io, [As real(omega) imag(omega) real(alpha) imag(alpha) real(beta) imag(beta) Ma err]);
                flush(io);
            end
        elseif k==5
            for i=1:ki
                Ma=Ma+delta;
                alpha,err,=op.mull(omega, alpha, lambda_u_, lambda_uz_, Ma)
                @show omega, Ma, As
                @show alpha
                DelimitedFiles.writedlm(io, [As real(omega) imag(omega) real(alpha) imag(alpha) real(beta) imag(beta) Ma err]);
                flush(io);
            end
        elseif k==-1
            break;
        elseif k==0
            f1=op.writeeigen_matrix(omega, alpha, lambda_u_, lambda_uz_, Ma);
            @show f1
            p=op.plot_inner(f1)
            io_inner=open("function_inner.plt","w")
            for i=1:32
                DelimitedFiles.writedlm(io_inner, [i real(p[i]) imag(p[i]) abs(p[i])]);
            end
            DelimitedFiles.writedlm(io_inner, [33 real(p[1]) imag(p[1]) abs(p[1])]);
            close(io_inner)
            @show alpha, Ma
            pp,p01=op.plot_outter(f1, alpha, Ma)
            io_outter=open("function_outter.plt","w")
            for j=1:200
                for i=1:32
                    DelimitedFiles.writedlm(io_outter, [i j-1 real(pp[i,j]) imag(pp[i, j]) abs(pp[i, j])]);
                end
                DelimitedFiles.writedlm(io_outter, [33 j-1 real(pp[1, j]) imag(pp[1, j]) abs(pp[1, j])]);
            end
            close(io_outter)
            io_01=open("function_outter_01.plt","w")
            for j=1:200
                for i=1:2
                    DelimitedFiles.writedlm(io_01, [i j-1 real(p01[i,j]) imag(p01[i, j]) abs(p01[i, j])]);
                end
            end
            close(io_01)
        else
            println("Please input a correct value!")
            continue;
        end
    end

    close(io)


end

module others

@doc """
          input(prompt::String="")::String

      Read a string from STDIN. The trailing newline is stripped.

      The prompt string, if given, is printed to standard output without a
      trailing newline before reading input.
      """ ->
      function input(prompt::String="")
          print(prompt)
          return readline(stdin)
      end

end



import .op
#import PyCall
using PyCall
@pyimport matplotlib
matplotlib.use("TkAgg")
@pyimport matplotlib.pyplot as plt

As=0.0e0;
Ma=1.5e0;
Ny=200;
ymax=10;
# omega=7.25;
# omega=4.0;
omega=4.0e0
alpha=0.9710066271585998	-0.1528922495786666im;
# alpha=1.22735-0.0256407im;

omega_end=20.0;
omega_begin=0.01;
domega=0.2;
op.initial();
beta=0.0e0;
@show omega, As, Ma, alpha

# alphawatch=betasearch_old(omega, alpha, beta, As, Ma)
# alphawatch=masearch_old(omega, alpha, beta, As, Ma)
 # @show alphawatch
 # alpha=op.mull_old(omega, alpha, 3.0,lambda_u_, lambda_uz_, Ma)
 # @show alpha

 As=0.3e0
 omega=1.32e0
 # alpha=0.2940294737880188	-0.04086253198398886im
 #0.7087936128382296	-0.05625760414966105
 alpha=0.416696728449895	-9.37516999932633e-5im;
 beta=1.0
 Ma=1.5
 #
 # alpha=0.21192346336157836	-0.017584860959858364im
 # omega=1.0e0
# (lambda_u_, lambda_uz_)=op.constvar.baseflow(As)
# alpha, V=op.mull(omega, alpha, lambda_u_, lambda_uz_, Ma)
# @show alpha
# @show V
 # search_by_input(omega, alpha, beta, As, Ma)
 search_by_input_new(omega, alpha, As, Ma, "Ma"*string(Ma)*"_As"*string(As)*"_beta"*string(beta)*"_omega_vs_alpha_even.plt")
 # search_by_input_new(omega, alpha, As, Ma, "Ma"*string(Ma)*"_omega"*string(omega)*"_beta"*string(beta)*"_As_vs_alpha_even.plt")
 # search_by_auto_new(omega, alpha, beta, As, Ma, "Ma"*string(Ma)*"_omega"*string(omega)*"_beta"*string(beta)*"_As_vs_alpha_even_auto.plt")


#  # alpha=1.2356455778148874 - 0.25607803494449904im;
#  alpha=0.9869306649201584	-0.15296330283107484im;
# @show omega, Ma, alpha
 # (alpha, p)=op.mull(omega, alpha, lambda_u_, lambda_uz_, Ma)
 # @show alpha


# (Maplt1, omegaplt1, alphaplt1)=backsearch_with_Ma(omega, alpha, As, Ma, 3.0, 0.05)

       # omegaplt=[omegaplt2;omegaplt1];
# alphaplt=[alphaplt2;alphaplt1];
# plt.plot(Maplt1, real(alphaplt1))
# plt.savefig("Ma_vs_max_alpha_r.png")
# plt.close()
# plt.plot(Maplt1, imag(alphaplt1))
# plt.savefig("Ma_vs_max_alpha_i.png")
# plt.close()
# plt.plot(Maplt1, real(omegaplt1))
# plt.savefig("Ma_vs_max_omega_r.png")
# plt.close()
#
# io=open("Ma_info.txt", "a");
# for i=1:length(Maplt1)
#     writedlm(io, [Maplt1[i] real(alphaplt1[i]) imag.(alphaplt1[i]) real.(omegaplt1[i])])
# end;

#y=[i for i in 0:Ny]/Ny*ymax
#op.get_w(omega, alpha, lambda_u_, y, p)
# w=get_w()
# aaa=[2,2,2];
#search(omega, alpha, As, Ma)





# @time (omegaplt1, alphaplt1)=backwardsearch(omega, alpha, omega_end, domega, As, Ma);
# @time (omegaplt2, alphaplt2)=forwardsearch(omega-domega, alpha, omega_begin, -domega, As, Ma);
# omegaplt=[omegaplt2;omegaplt1];
# alphaplt=[alphaplt2;alphaplt1];
# plt.plot(omegaplt, real(alphaplt))
# plt.savefig("omega_vs_alpha_r.png")
# plt.close()
# plt.plot(omegaplt, imag(alphaplt))
# plt.savefig("omega_vs_alpha_i.png")
# plt.close()
