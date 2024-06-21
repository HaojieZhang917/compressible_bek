function rayleigh_quotient_iteration(A, B, sigma; q0=rand(size(A, 1), 1))

    flg = true
    i=1
    while flg
        i=i+1
        sigma0 = real(sigma[1]) + abs(imag(sigma[1]))im + 0.0e0im
        q = (A - sigma*B) \ (B*q0)
        q0 = q/maximum(abs.(q))
        sigma = ((q0'*(A*q0))/(q0'*(B*q0)))[1]
        if abs(sigma-sigma0)<=eps(1.0f0)
            flg = false
        end
        if i==20
            flg=false
        end
    end

      return sigma, q0
end