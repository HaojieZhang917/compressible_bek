module FDqScheme

    using Polynomials, GenericLinearAlgebra, SparseArrays, LinearAlgebra

    struct LagrangePolynomial
        order::Int64
        p_0
    end

    function construct_lagrange_polynomial(x)
        p_0 = Array{Polynomial, 1}(undef, length(x))
        for i in 1:length(x)
            p_0[i] = Polynomial([-x[i], 1])
        end
        lagrange_polynomial=LagrangePolynomial(length(x)-1, p_0)
        return lagrange_polynomial
    end

    function polyval(poly::LagrangePolynomial, x)
        return prod([poly.p_0[i](x) for i in 1:length(poly.p_0)])
    end

    #返回x0到x1中的极值
    function find_extreme(poly::LagrangePolynomial, x0, x1)

        xrange=range(x0, x1, length=11)
        idx = 0
        while (xrange[end]-xrange[begin])>eps(1.0e0)*20.0e0
            y=[polyval(poly, x) for x in xrange]
            value, idx=findmax(abs.(y))
            xrange=range(xrange[max(idx-1, 1)], xrange[min(idx+1, 11)], length=11)
        end
        return xrange[idx]
    end


    function lagrange_pi(poly::LagrangePolynomial, x, xi)
        # lagrange_pi=LagrangePolynomial(poly.order, poly.p_0)
        lagrange_pi = LagrangePolynomial(poly.order, [poly.p_0[i] for i in 1:poly.order+1])
        lagrange_pi.p_0[findall(x->x==xi, x)].=Polynomial([1.0e0])
        return lagrange_pi
    end

    function derivative_lagrange_pi(poly::LagrangePolynomial, x)

        lagrange_pi=LagrangePolynomial(poly.order, [poly.p_0[i] for i in 1:poly.order+1])
        derivative_lagrange_pi = 0.0e0
        for i in 1:poly.order+1
            derivative_langrange_pi_local = 1.0e0
            for j in 1:poly.order+1
                if j!=i
                    derivative_langrange_pi_local = derivative_langrange_pi_local*poly.p_0[j](x)
                    # @show j, poly.p_0[j](x)
                else
                    derivative_langrange_pi_local = derivative_langrange_pi_local*derivative(poly.p_0[j])(x)
                    # @show j, derivative(poly.p_0[j])(x)
                end
            end
            # @show derivative_langrange_pi_local
            derivative_lagrange_pi=derivative_lagrange_pi + derivative_langrange_pi_local
        end
        return derivative_lagrange_pi
    end

    function LagrangeScheme(x, x0)
        pi0=construct_lagrange_polynomial(x)
        iloc = findall(x->x==x0, x)[1]
        LagrangeScheme=[]
        for i in 1:length(x)
            push!(LagrangeScheme, 
                derivative_lagrange_pi(lagrange_pi(pi0, x, x[i]), x[iloc])/polyval(lagrange_pi(pi0, x, x[i]), x[i]))
        end
        return LagrangeScheme
    end

    #construct q order interpolation of degree q-1
    function construct_si(N, q)
        N = N - 1
        si=[]
        if isodd(q)
            for i in 1:Int32((q-1)/2)
                push!(si, 0)
            end
            for i in 0:N-q
                push!(si, i)
            end
            for i in 1:Int32((q+1)/2)
                push!(si, N-q)
            end
        else
            for i in 1:Int32(q/2)
                push!(si, 0)
            end
            for i in 0:N-q
                push!(si, i)
            end
            for i in 1:Int32(q/2)
                push!(si, N-q)
            end
        end
        return reshape(Int64.(si), N+1, :)
    end

    # return a matrix for a distributed nodes in a section
    function construct_schemes(x::Vector, q)
        n = length(x)
        si = construct_si(n, q) .+ 1
        scheme_matrix = zeros(n, n)
        for i in 1:n
            scheme_matrix[i, si[i]:si[i]+q]=LagrangeScheme(x[si[i]:si[i]+q], x[i])
        end
        return (scheme_matrix)
    end

    #construct_polynomial factor pi_y
    function construct_polynomial(y, N, q, si)
        pi_y = Array{LagrangePolynomial}(undef, N-1)
        for i in 1:N-1
            pi_y[i] = construct_lagrange_polynomial(y[si[i]+1:si[i]+q+1])
        end
        return pi_y
    end

    # search_extrema_xi
    function search_extrema_xi(y, pi_y, N)
        x=Array{Float64}(undef, N+1)
        x[1] = -1#*Float64(N)
        x[end] = 1#*Float64(N)
        for i in 1:N-1
            x[i+1] = find_extreme(pi_y[i], y[i], y[i+1])
        end
        return x
    end

    #get_f
    function get_f(y, N, q, si, x, flg=true)
        #construct_polynomial factor pi_y
            pi_y = construct_polynomial(y, N, q-1, si)
        # search extrema xi
        if flg
            x_extrema = search_extrema_xi(y, pi_y, N)
        else
            x_extrema = x
        end
        #return abs(pi(x)) at xi
        abs_pi_x = []

        push!(abs_pi_x, abs(polyval(pi_y[1], x_extrema[1])))
        for i in 1:N-1
            push!(abs_pi_x, abs(polyval(pi_y[i], x_extrema[i+1])))
        end
        push!(abs_pi_x, abs(polyval(pi_y[N-1], x_extrema[N+1])))

        f= []
        for i in 1:N
            push!(f, abs_pi_x[i+1]-abs_pi_x[i])
        end
        return Float64.(f), Float64.(x_extrema), pi_y
    end

    function construct_J(y, x, f, N, q, si)
        epsilon = sqrt(eps(1.0e0))
        ymin = 1.0e-6
        J = zeros(N, N)

        for i in 1:N
            h = zeros(N)
            if abs(y[i])>ymin
                h[i] = epsilon*y[i]
            else
                h[i] = ymin*sign(y[i]+eps(1.0e0))
            end
            y_new=y+h
            f_new, = get_f(y_new, N, q, si, x, false)
            J[:, i] = (f_new-f)/h[i]
        end
        return sparse(J)
    end


    function FD_q_Scheme(n, q_target)

        N = n - 1

        # guess y_i using q=N, equal space
        y = Array{Float64}(undef, N)
        for i in 1:N
            y[i] = (2.0e0/Float64(2N)-1.0e0+2.0e0/Float64(N)*Float64(i-1))
        end

        x=zeros(N+1)

        for q in 2:1:q_target
            si = construct_si(N, q-1)

            f, x, pi_y= get_f(y, N, q, si, x)

            maxdy = 1.0e0
            while maxdy > 2.2e-6
                J = construct_J(y, x, f, N, q, si)
                dy = - J\f

                if(y+dy!=sort(vec(y+dy))) 
                    @show J
                    @show J[1, 1], J[end, end]
                    sleep()
                end
                y = y+ dy
                f, x= get_f(y, N, q, si, x)
                maxdy = maximum(abs.(dy))
            end
        end

        return x, construct_schemes(collect(x), q_target)
    end

end