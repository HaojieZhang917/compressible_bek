module BlackburnStability

using LinearAlgebra
using NonlinearEigenproblems

export blackburn_spatial_matrices, blackburn_temporal_matrices
export apply_homogeneous_boundaries, spatial_spectrum, temporal_spectrum
export reconstruct_reduced_mode
export eigsol_blackburn
export tracked_spatial_mode, find_neutral_R, track_neutral_curve
export neutral_curve_matrix

const C = ComplexF64

complex_eye(n) = Matrix{C}(I, n, n)
complex_zero(n) = zeros(C, n, n)
diag_matrix(values) = Matrix(Diagonal(C.(values)))

function check_inputs(F, GL, H, T, D, D2)
    n = length(F)
    length(GL) == n || throw(DimensionMismatch("GL has the wrong length"))
    length(H) == n || throw(DimensionMismatch("H has the wrong length"))
    length(T) == n || throw(DimensionMismatch("T has the wrong length"))
    size(D) == (n, n) || throw(DimensionMismatch("D has the wrong size"))
    size(D2) == (n, n) || throw(DimensionMismatch("D2 has the wrong size"))
    return n
end

"""
    blackburn_spatial_matrices(F, GL, H, T, Pr, D, D2, R, beta, omega;
                              frame=:disk)

Build `(L0,L1,L2)` for

    (L0 + alpha*L1 + alpha^2*L2) q = 0,

where `q=[u,v,w,theta,p]`. `GL` is the azimuthal deficit of the shared
steady generalized-Boussinesq base flow, with
`GL(0)=0` and `GL(infinity)=1`. The inertial azimuthal profile is `Q=1-GL`.

For `frame=:disk`, `omega` is the disk-fixed frequency and the internal
inertial frequency is `omega + beta` in the present local scaling. For
`frame=:inertial`, `omega` is used directly. Blackburn et al.'s canonical
model weights the complete material acceleration by `chi=2-T`; consequently
the three momentum frequency blocks are `-im*omega_inertial*diag(chi)`.
"""
function blackburn_spatial_matrices(
    F, GL, H, T, Pr, D, D2, R, beta, omega; frame::Symbol=:disk,
)
    n = check_inputs(F, GL, H, T, D, D2)
    R > 0 || throw(ArgumentError("R must be positive"))
    Pr > 0 || throw(ArgumentError("Pr must be positive"))
    frame in (:disk, :inertial) || throw(ArgumentError("frame must be :disk or :inertial"))

    eye = complex_eye(n)
    zero = complex_zero(n)
    Dc = C.(D)
    D2c = C.(D2)

    F = vec(Float64.(F))
    GL = vec(Float64.(GL))
    H = vec(Float64.(H))
    T = vec(Float64.(T))
    chi = 2.0 .- T
    Q = 1.0 .- GL
    Fp = D * F
    Qp = D * Q
    Hp = D * H
    Tp = D * T

    Ar = F.^2 .+ H .* Fp .- Q.^2
    Atheta = 2.0 .* F .* Q .+ H .* Qp
    # Keeping Q inertial while shifting frequency separately makes the
    # chi-weighted material acceleration transform consistently between the
    # inertial and disk-fixed frames.
    omega_i = frame == :disk ? omega + beta : omega

    # Malik's local parallel system retains terms through O(1/R). Therefore
    # both the cylindrical vector-Laplacian corrections and the base axial
    # acceleration H*H'/R^2 are omitted from the disturbance operator.
    lap0 = D2c - beta^2 .* eye
    momentum_transport = (
        -im * omega_i .* diag_matrix(chi)
        + im * beta .* diag_matrix(chi .* Q)
        + diag_matrix(chi .* H ./ R) * Dc
        - lap0 ./ R
    )
    radial_transport = momentum_transport + diag_matrix(chi .* F ./ R)
    axial_transport = momentum_transport + diag_matrix(chi .* Hp ./ R)
    thermal_transport = (
        -im * omega_i .* eye
        + im * beta .* diag_matrix(Q)
        + diag_matrix(H ./ R) * Dc
        - lap0 ./ (Pr * R)
    )

    L0 = [
        radial_transport                 diag_matrix(-2.0 .* chi .* Q ./ R)  diag_matrix(chi .* Fp)  diag_matrix(-Ar ./ R)                 zero;
        diag_matrix(2.0 .* chi .* Q ./ R) radial_transport                  diag_matrix(chi .* Qp)  diag_matrix(-Atheta ./ R)             im * beta .* eye;
        zero                             zero                              axial_transport          zero                                    Dc;
        zero                             zero                              diag_matrix(Tp)           thermal_transport                       zero;
        eye ./ R                         im * beta .* eye                  Dc                        zero                                    zero
    ]

    L1 = [
        im .* diag_matrix(chi .* F)  zero                         zero                         zero  im .* eye;
        zero                         im .* diag_matrix(chi .* F)  zero                         zero  zero;
        zero                         zero                         im .* diag_matrix(chi .* F)  zero  zero;
        zero                         zero                         zero                         im .* diag_matrix(F)  zero;
        im .* eye                    zero                         zero                         zero  zero
    ]

    L2 = [
        eye ./ R  zero      zero      zero              zero;
        zero      eye ./ R  zero      zero              zero;
        zero      zero      eye ./ R  zero              zero;
        zero      zero      zero      eye ./ (Pr * R)  zero;
        zero      zero      zero      zero              zero
    ]

    metadata = (
        chi=chi, Q=Q, Fp=Fp, Qp=Qp, Hp=Hp, Tp=Tp,
        Ar=Ar, Atheta=Atheta, omega_inertial=omega_i,
        acceleration_model=:blackburn2021,
    )
    return L0, L1, L2, metadata
end

"""
    blackburn_temporal_matrices(F, GL, H, T, Pr, D, D2, R, alpha, beta;
                               frame=:disk)

Build `(A,B)` for `A*q = omega*B*q`. The returned eigenvalue is in the frame
selected by `frame`. The velocity blocks of `B` are `im*diag(chi)`, the
temperature block is `im*I`, and the pressure block is zero.
"""
function blackburn_temporal_matrices(
    F, GL, H, T, Pr, D, D2, R, alpha, beta; frame::Symbol=:disk,
)
    # First build the operator with zero frequency in the requested frame.
    L0, L1, L2, metadata = blackburn_spatial_matrices(
        F, GL, H, T, Pr, D, D2, R, beta, 0.0; frame=frame,
    )
    A = L0 + alpha .* L1 + alpha^2 .* L2
    n = length(F)
    eye = complex_eye(n)
    zero = complex_zero(n)
    chi_mass = diag_matrix(metadata.chi)
    M = [
        chi_mass zero zero zero zero;
        zero chi_mass zero zero zero;
        zero zero chi_mass zero zero;
        zero zero zero eye zero;
        zero zero zero zero zero
    ]
    B = im .* M
    return A, B, metadata
end

function homogeneous_boundary_indices(n; far_pressure::Bool=true)
    indices = [1, n, n + 1, 2n, 2n + 1, 3n, 3n + 1, 4n]
    far_pressure && push!(indices, 5n)
    return indices
end

"""
    apply_homogeneous_boundaries(matrices, n; far_pressure=true)

Eliminate wall and far-field degrees of freedom for `u,v,w,theta`. By default,
also impose the far-field decay condition `p(infinity)=0`; wall pressure is
retained and is determined by momentum and continuity. Set `far_pressure=false`
only for pressure-boundary sensitivity checks.
"""
function apply_homogeneous_boundaries(
    matrices::Tuple, n::Int; far_pressure::Bool=true,
)
    removed = homogeneous_boundary_indices(n; far_pressure=far_pressure)
    keep = setdiff(collect(1:5n), removed)
    reduced = map(matrix -> matrix[keep, keep], matrices)
    return reduced..., keep
end

function reconstruct_reduced_mode(q, keep, n)
    full = zeros(C, 5n)
    full[keep] .= q
    return full
end

function polynomial_residual(L0, L1, L2, alpha, q)
    numerator = norm((L0 + alpha .* L1 + alpha^2 .* L2) * q)
    scale = (norm(L0) + abs(alpha) * norm(L1) + abs2(alpha) * norm(L2)) * norm(q)
    return numerator / max(scale, eps(Float64))
end

"""Solve the dense companion generalized eigenproblem for spatial alpha."""
function spatial_spectrum(L0, L1, L2; keep=nothing, n_full=nothing)
    m = size(L0, 1)
    size(L1) == (m, m) || throw(DimensionMismatch("L1 has the wrong size"))
    size(L2) == (m, m) || throw(DimensionMismatch("L2 has the wrong size"))
    eye = complex_eye(m)
    zero = complex_zero(m)
    companion_A = [zero eye; -L0 -L1]
    companion_B = [eye zero; zero L2]
    decomposition = eigen(companion_A, companion_B)

    modes = NamedTuple[]
    for j in eachindex(decomposition.values)
        alpha = decomposition.values[j]
        isfinite(real(alpha)) && isfinite(imag(alpha)) || continue
        q = decomposition.vectors[1:m, j]
        norm(q) > 1.0e-12 || continue
        residual = polynomial_residual(L0, L1, L2, alpha, q)
        companion_error = norm(decomposition.vectors[m+1:2m, j] - alpha .* q) /
                          max(norm(alpha .* q), eps(Float64))
        thermal_fraction = NaN
        if keep !== nothing && n_full !== nothing
            full = reconstruct_reduced_mode(q, keep, n_full)
            velocity_norm = norm(full[1:3n_full])
            thermal_norm = norm(full[3n_full+1:4n_full])
            thermal_fraction = thermal_norm / max(velocity_norm + thermal_norm, eps(Float64))
        end
        push!(modes, (
            alpha=alpha, residual=residual, companion_error=companion_error,
            thermal_fraction=thermal_fraction, vector=q,
        ))
    end
    return modes
end

"""
    eigsol_blackburn(F, GL, H, T, R, omega, beta, N_cheb, D, D2, c, num;
                     Pr=0.72, frame=:disk, tol=1e-10, maxit=100,
                     initial_vector=nothing, full_eigenvectors=false,
                     return_info=false)

Target the `num` spatial eigenvalues nearest the shift `c` for the Blackburn
canonical generalized-Boussinesq operator. The positional arguments
intentionally match the existing notebook `eigsol` interface. `GL` is the
azimuthal deficit of the shared steady base flow, with `GL(0)=0` and
`GL(infinity)=1`.

The default eigenvectors contain only the unconstrained collocation degrees of
freedom, which makes the return value directly compatible with the current
neutral-curve tracker. Set `full_eigenvectors=true` to reinsert the nine zero
boundary values for `u,v,w,theta` and `p(infinity)`.

For continuation, use `num=1` and update `c` with the previous eigenvalue. If
several eigenvalues are requested and IAR only converges a subset, the function
returns the subset whose original polynomial residual is converged and reports
`info.partial=true` when `return_info=true`.
"""
function eigsol_blackburn(
    F, GL, H, T, R, omega, beta, N_cheb, D, D2, c, num;
    Pr::Real=0.72,
    frame::Symbol=:disk,
    tol::Real=1.0e-10,
    maxit::Integer=500,
    initial_vector=nothing,
    full_eigenvectors::Bool=false,
    return_info::Bool=false,
)
    n = check_inputs(F, GL, H, T, D, D2)
    n == N_cheb + 1 || throw(DimensionMismatch(
        "N_cheb must be one less than the number of collocation points",
    ))
    num >= 1 || throw(ArgumentError("num must be positive"))

    L0, L1, L2, metadata = blackburn_spatial_matrices(
        F, GL, H, T, Pr, D, D2, R, beta, omega; frame=frame,
    )
    L0b, L1b, L2b, keep = apply_homogeneous_boundaries(
        (L0, L1, L2), n,
    )
    problem = PEP([L0b, L1b, L2b])
    shift = ComplexF64(c)
    # NonlinearEigenproblems names its shift keyword with Greek sigma (U+03C3).
    # Constructing the keyword dynamically keeps this source ASCII.
    iar_options = Dict{Symbol, Any}(
        Symbol(Char(0x03c3)) => shift,
        :neigs => num,
        :maxit => maxit,
        :tol => tol,
    )
    if initial_vector !== nothing
        length(initial_vector) == size(L0b, 1) || throw(DimensionMismatch(
            "initial_vector must use the reduced boundary-condition layout",
        ))
        iar_options[:v] = ComplexF64.(initial_vector)
    end
    partial = false
    try
        values, vectors = iar(problem; iar_options...)
    catch exception
        if exception isa NonlinearEigenproblems.NEPCore.NoConvergenceException
            values = getfield(exception, 1)
            vectors = getfield(exception, 2)
            partial = true
        else
            rethrow()
        end
    end

    values = ComplexF64.(vec(values))
    order = sortperm(eachindex(values); by=i -> abs(values[i] - shift))
    values = values[order]
    vectors = ComplexF64.(vectors[:, order])
    residuals = [
        polynomial_residual(L0b, L1b, L2b, values[j], vectors[:, j])
        for j in eachindex(values)
    ]
    if partial
        converged = findall(residual -> residual <= max(100tol, 1.0e-9), residuals)
        isempty(converged) && throw(ErrorException(
            "IAR did not return a polynomial eigenpair with a converged residual",
        ))
        values = values[converged]
        vectors = vectors[:, converged]
        residuals = residuals[converged]
    end

    if full_eigenvectors
        full = zeros(C, 5n, length(values))
        for j in eachindex(values)
            full[:, j] .= reconstruct_reduced_mode(vectors[:, j], keep, n)
        end
        vectors = full
    end

    if return_info
        info = (
            residuals=residuals,
            keep=keep,
            metadata=metadata,
            reduced_size=size(L0b, 1),
            partial=partial,
        )
        return values, vectors, info
    end
    return values, vectors
end

function mode_overlap(left, right)
    denominator = norm(left) * norm(right)
    denominator > eps(Float64) || return 0.0
    return abs(dot(left, right)) / denominator
end

"""
    tracked_spatial_mode(F, GL, H, T, R, omega, beta, N_cheb, D, D2,
                         alpha_seed; vector_seed=nothing, kwargs...)

Solve a small shifted eigenproblem and select the mode that is continuous with
the supplied eigenpair. A candidate must pass the polynomial-residual,
eigenvalue-jump, and eigenvector-overlap checks.
"""
function tracked_spatial_mode(
    F, GL, H, T, R, omega, beta, N_cheb, D, D2, alpha_seed;
    vector_seed=nothing,
    num_candidates::Integer=2,
    min_overlap::Real=0.60,
    max_alpha_jump::Real=0.05,
    residual_tol::Real=1.0e-8,
    eig_tol::Real=1.0e-10,
    maxit::Integer=800,
    frame::Symbol=:disk,
)
    num_candidates >= 1 || throw(ArgumentError("num_candidates must be positive"))
    values, vectors, info = eigsol_blackburn(
        F, GL, H, T, R, omega, beta, N_cheb, D, D2,
        alpha_seed, num_candidates;
        frame=frame,
        tol=eig_tol,
        maxit=maxit,
        initial_vector=vector_seed,
        return_info=true,
    )
    isempty(values) && error("the shifted eigensolver returned no eigenvalues")

    distances = abs.(values .- alpha_seed)
    overlaps = vector_seed === nothing ? fill(NaN, length(values)) : [
        mode_overlap(vector_seed, vectors[:, j]) for j in eachindex(values)
    ]
    valid = findall(info.residuals .<= residual_tol)
    isempty(valid) && error(
        "no candidate satisfies residual_tol=$residual_tol; " *
        "minimum residual=$(minimum(info.residuals))",
    )

    if vector_seed === nothing
        index = valid[argmin(distances[valid])]
    else
        continuous = filter(
            j -> overlaps[j] >= min_overlap && distances[j] <= max_alpha_jump,
            valid,
        )
        isempty(continuous) && error(
            "mode-continuity check failed: best overlap=$(maximum(overlaps[valid])), " *
            "smallest alpha jump=$(minimum(distances[valid]))",
        )
        index = continuous[argmax(overlaps[continuous])]
    end

    distance = distances[index]
    distance <= max_alpha_jump || error(
        "eigenvalue jump $distance exceeds max_alpha_jump=$max_alpha_jump",
    )
    vector = copy(vectors[:, index])
    vector ./= norm(vector)
    overlap = vector_seed === nothing ? NaN : overlaps[index]
    if vector_seed !== nothing
        phase = dot(vector_seed, vector)
        abs(phase) > eps(Float64) && (vector .*= conj(phase) / abs(phase))
    end
    return (
        alpha=values[index], vector=vector,
        residual=info.residuals[index], overlap=overlap,
        alpha_jump=distance, candidate_index=index,
        candidate_values=values, candidate_vectors=vectors,
        candidate_residuals=info.residuals,
        candidate_overlaps=overlaps,
    )
end

function with_parameters(mode, R, beta; scan_iterations=0, refine_iterations=0)
    return merge(mode, (
        R=Float64(R), beta=Float64(beta),
        scan_iterations=scan_iterations,
        refine_iterations=refine_iterations,
    ))
end

function try_tracked_mode(args...; kwargs...)
    try
        return tracked_spatial_mode(args...; kwargs...)
    catch exception
        exception isa ErrorException || rethrow()
        return nothing
    end
end

function refine_neutral_bracket(
    F, GL, H, T, omega, beta, N_cheb, D, D2, left, right;
    neutral_tol, R_tol, max_refine, mode_options,
)
    left.R <= right.R || ((left, right) = (right, left))
    imag(left.alpha) * imag(right.alpha) <= 0 || error(
        "neutral refinement requires a sign-changing bracket",
    )

    for iteration in 1:max_refine
        y_left = imag(left.alpha)
        y_right = imag(right.alpha)
        abs(y_left) <= neutral_tol && return merge(left, (refine_iterations=iteration - 1,))
        abs(y_right) <= neutral_tol && return merge(right, (refine_iterations=iteration - 1,))

        width = right.R - left.R
        width <= R_tol && return merge(
            abs(y_left) <= abs(y_right) ? left : right,
            (refine_iterations=iteration - 1,),
        )
        R_trial = (left.R * y_right - right.R * y_left) / (y_right - y_left)
        margin = 0.1 * width
        R_trial = clamp(R_trial, left.R + margin, right.R - margin)
        seed = abs(R_trial - left.R) <= abs(right.R - R_trial) ? left : right
        trial_mode = try_tracked_mode(
            F, GL, H, T, R_trial, omega, beta, N_cheb, D, D2, seed.alpha;
            vector_seed=seed.vector, mode_options...,
        )
        if trial_mode === nothing
            R_trial = 0.5 * (left.R + right.R)
            seed = abs(R_trial - left.R) <= abs(right.R - R_trial) ? left : right
            trial_mode = tracked_spatial_mode(
                F, GL, H, T, R_trial, omega, beta, N_cheb, D, D2, seed.alpha;
                vector_seed=seed.vector, mode_options...,
            )
        end
        trial = with_parameters(
            trial_mode, R_trial, beta;
            scan_iterations=max(left.scan_iterations, right.scan_iterations),
            refine_iterations=iteration,
        )
        abs(imag(trial.alpha)) <= neutral_tol && return trial
        if imag(left.alpha) * imag(trial.alpha) <= 0
            right = trial
        else
            left = trial
        end
    end
    return abs(imag(left.alpha)) <= abs(imag(right.alpha)) ? left : right
end

"""
    find_neutral_R(F, GL, H, T, omega, beta, N_cheb, D, D2;
                   R_guess, alpha_seed, vector_seed=nothing, kwargs...)

Track one spatial mode while scanning in `R`, bracket `imag(alpha)=0`, and
refine the neutral point. Failed or discontinuous eigensolver steps are rejected.
"""
function find_neutral_R(
    F, GL, H, T, omega, beta, N_cheb, D, D2;
    R_guess::Real,
    alpha_seed,
    vector_seed=nothing,
    R_step::Real=1.0,
    preferred_direction::Integer=0,
    max_scan_steps::Integer=80,
    max_refine::Integer=30,
    neutral_tol::Real=1.0e-7,
    R_tol::Real=1.0e-4,
    max_R_deviation::Real=Inf,
    num_candidates::Integer=2,
    min_overlap::Real=0.60,
    max_alpha_jump::Real=0.05,
    residual_tol::Real=1.0e-8,
    eig_tol::Real=1.0e-10,
    maxit::Integer=800,
    frame::Symbol=:disk,
)
    R_guess > 0 || throw(ArgumentError("R_guess must be positive"))
    R_step > 0 || throw(ArgumentError("R_step must be positive"))
    preferred_direction in (-1, 0, 1) || throw(ArgumentError(
        "preferred_direction must be -1, 0, or 1",
    ))
    mode_options = (
        num_candidates=num_candidates,
        min_overlap=min_overlap,
        max_alpha_jump=max_alpha_jump,
        residual_tol=residual_tol,
        eig_tol=eig_tol,
        maxit=maxit,
        frame=frame,
    )
    center_mode = tracked_spatial_mode(
        F, GL, H, T, R_guess, omega, beta, N_cheb, D, D2, alpha_seed;
        vector_seed=vector_seed, mode_options...,
    )
    center = with_parameters(center_mode, R_guess, beta)
    abs(imag(center.alpha)) <= neutral_tol && return center

    neighbors = NamedTuple[]
    roots = NamedTuple[]
    for direction in (-1, 1)
        R_next = R_guess + direction * R_step
        R_next > 0 || continue
        next_mode = try_tracked_mode(
            F, GL, H, T, R_next, omega, beta, N_cheb, D, D2, center.alpha;
            vector_seed=center.vector, mode_options...,
        )
        next_mode === nothing && continue
        next = with_parameters(next_mode, R_next, beta; scan_iterations=1)
        if imag(center.alpha) * imag(next.alpha) <= 0
            root = refine_neutral_bracket(
                F, GL, H, T, omega, beta, N_cheb, D, D2, center, next;
                neutral_tol=neutral_tol, R_tol=R_tol, max_refine=max_refine,
                mode_options=mode_options,
            )
            abs(root.R - R_guess) <= R_step && return root
            push!(roots, root)
            continue
        end
        push!(neighbors, merge(next, (direction=direction,)))
    end
    isempty(neighbors) && error("both initial R-continuation steps were rejected")

    sort!(neighbors; by=neighbor -> (
        preferred_direction == 0 || neighbor.direction == preferred_direction ? 0 : 1,
        abs(imag(neighbor.alpha)),
    ))
    for neighbor in neighbors
        previous = center
        current = neighbor
        direction = neighbor.direction
        for iteration in 2:max_scan_steps
            R_next = current.R + direction * R_step
            R_next > 0 || break
            next_mode = try_tracked_mode(
                F, GL, H, T, R_next, omega, beta, N_cheb, D, D2, current.alpha;
                vector_seed=current.vector, mode_options...,
            )
            next_mode === nothing && break
            next = with_parameters(next_mode, R_next, beta; scan_iterations=iteration)
            if imag(current.alpha) * imag(next.alpha) <= 0
                root = refine_neutral_bracket(
                    F, GL, H, T, omega, beta, N_cheb, D, D2, current, next;
                    neutral_tol=neutral_tol, R_tol=R_tol, max_refine=max_refine,
                    mode_options=mode_options,
                )
                push!(roots, root)
                break
            end
            previous = current
            current = next
            if iteration >= 4 &&
               abs(imag(current.alpha)) > abs(imag(previous.alpha)) &&
               abs(imag(previous.alpha)) > abs(imag(center.alpha))
                break
            end
        end
    end
    isempty(roots) && error(
        "failed to bracket a neutral point at beta=$beta from R_guess=$R_guess",
    )
    root = roots[argmin([abs(candidate.R - R_guess) for candidate in roots])]
    deviation = abs(root.R - R_guess)
    deviation <= max_R_deviation || error(
        "nearest neutral root is $deviation away from R_guess=$R_guess, " *
        "exceeding max_R_deviation=$max_R_deviation",
    )
    return root
end

"""
    track_neutral_curve(F, GL, H, T, omega, N_cheb, D, D2; kwargs...)

Continue neutral points with a fixed beta step. The next beta is attempted only
after the current neutral point has converged. The result contains `points` and
an optional `failure` record.
"""
function track_neutral_curve(
    F, GL, H, T, omega, N_cheb, D, D2;
    beta_start::Real,
    beta_end::Real,
    beta_step::Real=-8.0e-4,
    R_start::Real,
    alpha_start,
    vector_start=nothing,
    initial_R_delta::Real=0.0,
    initial_alpha_delta=0.0 + 0.0im,
    verbose::Bool=true,
    on_point::Function=(point, index) -> nothing,
    kwargs...,
)
    beta_step != 0 || throw(ArgumentError("beta_step must be nonzero"))
    direction = sign(beta_step)
    direction * (beta_end - beta_start) >= 0 || throw(ArgumentError(
        "beta_step points away from beta_end",
    ))

    points = NamedTuple[]
    failure = nothing
    beta = Float64(beta_start)
    R_guess = Float64(R_start)
    alpha_seed = ComplexF64(alpha_start)
    vector_seed = vector_start
    preferred_direction = 0
    reached_end(value) = direction > 0 ? value > beta_end + 10eps(Float64) :
                                        value < beta_end - 10eps(Float64)

    while !reached_end(beta)
        point = try
            find_neutral_R(
                F, GL, H, T, omega, beta, N_cheb, D, D2;
                R_guess=R_guess,
                alpha_seed=alpha_seed,
                vector_seed=vector_seed,
                preferred_direction=preferred_direction,
                kwargs...,
            )
        catch exception
            failure = (
                beta=beta, R_guess=R_guess,
                message=sprint(showerror, exception),
            )
            verbose && println(
                "Neutral continuation stopped at beta=$beta: $(failure.message)",
            )
            break
        end
        push!(points, point)
        on_point(point, length(points))
        verbose && println(
            "beta=$(round(beta; digits=7)), " *
            "R=$(round(point.R; digits=6)), " *
            "alpha=$(point.alpha), residual=$(point.residual), " *
            "overlap=$(point.overlap)",
        )

        if length(points) >= 2
            current = points[end]
            previous = points[end - 1]
            preferred_direction = Int(sign(current.R - previous.R))
            R_guess = current.R + (current.R - previous.R)
            alpha_seed = current.alpha + (current.alpha - previous.alpha)
        else
            R_guess = point.R + initial_R_delta
            alpha_seed = point.alpha + initial_alpha_delta
            preferred_direction = Int(sign(initial_R_delta))
        end
        vector_seed = point.vector
        beta += beta_step
    end
    return (points=points, failure=failure)
end

function neutral_curve_matrix(result, omega)
    data = zeros(Float64, length(result.points), 7)
    for (row, point) in enumerate(result.points)
        data[row, :] .= (
            omega, point.R, point.beta, real(point.alpha), imag(point.alpha),
            point.residual, point.overlap,
        )
    end
    return data
end

"""Solve the dense descriptor generalized eigenproblem for temporal omega."""
function temporal_spectrum(A, B; keep=nothing, n_full=nothing)
    decomposition = eigen(A, B)
    modes = NamedTuple[]
    for j in eachindex(decomposition.values)
        omega = decomposition.values[j]
        isfinite(real(omega)) && isfinite(imag(omega)) || continue
        q = decomposition.vectors[:, j]
        denominator = (norm(A) + abs(omega) * norm(B)) * norm(q)
        residual = norm((A - omega .* B) * q) / max(denominator, eps(Float64))
        thermal_fraction = NaN
        if keep !== nothing && n_full !== nothing
            full = reconstruct_reduced_mode(q, keep, n_full)
            velocity_norm = norm(full[1:3n_full])
            thermal_norm = norm(full[3n_full+1:4n_full])
            thermal_fraction = thermal_norm / max(velocity_norm + thermal_norm, eps(Float64))
        end
        push!(modes, (
            omega=omega, residual=residual,
            thermal_fraction=thermal_fraction, vector=q,
        ))
    end
    return modes
end

end # module BlackburnStability
