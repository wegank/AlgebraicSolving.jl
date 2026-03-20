algebraic_solving_lock = ReentrantLock()

function _real_roots(f::Vector{QQMPolyRingElem})::Vector{Vector{Vector{QQFieldElem}}}
    rs = nothing
    lock(algebraic_solving_lock) do
        rs = real_solutions(AlgebraicSolving.Ideal(f), interval=true)
    end
    return rs
end

function _from_univariate(R::QQMPolyRing, f::QQPolyRingElem)::QQMPolyRingElem
    @assert length(gens(R)) == 1
    M = MPolyBuildCtx(R)
    for i in 0:length(f)
        push_term!(M, coeff(f, i), [i])
    end
    finish(M)
end

function _sample_points(f::Vector{QQPolyRingElem})::Vector{QQFieldElem}
    R, _ = polynomial_ring(QQ, [:x])
    fₓ = map(p -> _from_univariate(R, p), f)
    factors = unique([p for fₓ′ in fₓ for (p, _) in factor(fₓ′)])
    # We map each factor to its roots, and each root to its factor
    roots_by_factor = Dict{QQMPolyRingElem,Vector{Vector{Vector{QQFieldElem}}}}()
    factors_by_root = Dict{Vector{Vector{QQFieldElem}},QQMPolyRingElem}()
    for factor in factors
        roots = _real_roots([factor])
        roots_by_factor[factor] = roots
        for r in roots
            factors_by_root[r] = factor
        end
    end
    # We order intervals by their left endpoint
    roots = sort(collect(keys(factors_by_root)), by=x -> x[1][1])
    # If intervals are not ordered by their right endpoint,
    # we merge the offending factors
    while (!issorted(roots, by=x -> x[1][2]))
        i = findfirst(i -> roots[i][1][2] > roots[i+1][1][2], 1:length(roots)-1)
        bad_factors = [factors_by_root[roots[i]], factors_by_root[roots[i+1]]]
        # Remove old factors and roots
        for bad_factor in bad_factors
            for r in roots_by_factor[bad_factor]
                delete!(factors_by_root, r)
            end
            delete!(roots_by_factor, bad_factor)
        end
        # Replace with merged factor and its roots
        merged_factor = prod(bad_factors)
        merged_roots = _real_roots([merged_factor])
        roots_by_factor[merged_factor] = merged_roots
        for r in merged_roots
            factors_by_root[r] = merged_factor
        end
        roots = sort(collect(keys(factors_by_root)), by=r -> r[1][1])
    end
    # Now the intervals are properly ordered, we can sample points
    points = Vector{QQFieldElem}()
    if length(roots) == 0
        push!(points, QQ(0))
        return points
    end
    push!(points, floor(roots[1][1][1]) - 1)
    for i in 1:length(roots)-1
        push!(points, (roots[i][1][2] + roots[i+1][1][1]) // 2)
    end
    push!(points, ceil(roots[end][1][2]) + 1)
    return points
end

function _sample_points_0(f::Vector{QQMPolyRingElem})::Vector{Vector{QQFieldElem}}
    @assert all(map(is_constant, f))
    return [Vector{QQFieldElem}()]
end

function _sample_points_1(f::Vector{QQMPolyRingElem})::Vector{Vector{QQFieldElem}}
    @assert all(map(is_univariate, f))
    R, _ = polynomial_ring(QQ, :x)
    return [[p] for p in _sample_points(map(p -> to_univariate(R, p), f))]
end

function _sample_points_desc(n::Int)::String
    if n == 0
        return "Constant"
    elseif n == 1
        return "Univariate"
    elseif n == 2
        return "Bivariate"
    elseif n == 3
        return "Trivariate"
    else
        return "Multivariate"
    end
end

function _sample_points_2(
    f::Vector{QQMPolyRingElem},
    xs::Vector{QQMPolyRingElem};
    nr_thrds::Int=1,
    show_progress::Bool=false,
    desc::String="$(_sample_points_desc(length(xs))) sample points"
)::Vector{Vector{QQFieldElem}}
    @assert length(xs) >= 2
    x₁ = xs[1:end-1]
    x₂ = xs[end]
    factors = unique([p for fₓ in f for (p, _) in factor(fₓ)])
    v = Vector{QQMPolyRingElem}()
    for i in eachindex(factors)
        if !isempty(intersect(x₁, vars(factors[i])))
            push!(v, leading_coefficient(factors[i], length(xs)))
        end
        if x₂ in vars(factors[i])
            push!(v, Interpolation.discriminant(factors[i], length(xs); nr_thrds=nr_thrds))
        end
    end
    for i in 1:length(factors)-1
        for j in i+1:length(factors)
            if x₂ in vars(factors[i]) || x₂ in vars(factors[j])
                push!(v, Interpolation.resultant(factors[i], factors[j], length(xs); nr_thrds=nr_thrds))
            end
        end
    end
    p₁ = length(xs) == 2 ? _sample_points_1(v) : _sample_points_2(v, x₁; show_progress=show_progress)
    prog = Progress.ProgressBar(total=length(p₁); desc=desc, enabled=show_progress)
    Progress.update!(prog, 0)
    function _points_chunk(p::Vector{Vector{QQFieldElem}})::Vector{Vector{QQFieldElem}}
        res_chunk = Vector{Vector{QQFieldElem}}()
        for i in eachindex(p)
            p₂ = _sample_points_1(map(p′ -> evaluate(p′, x₁, p[i]), f))
            for j in eachindex(p₂)
                push!(res_chunk, [p[i]; p₂[j][1]])
            end
            Progress.next!(prog)
        end
        res_chunk
    end
    chunk_size = ceil(Int, length(p₁) / nr_thrds)
    chunks = [p₁[i:min(i + chunk_size - 1, end)] for i in 1:chunk_size:length(p₁)]
    tasks = [Threads.@spawn _points_chunk(chunk) for chunk in chunks]
    res = sort(vcat(fetch.(tasks)...))
    Progress.finish!(prog)
    res
end

function _unique_by_sign(p::Vector{Vector{QQFieldElem}}, f::Vector{QQMPolyRingElem})::Vector{Vector{QQFieldElem}}
    signs = Set{Vector{QQFieldElem}}()
    p′ = Vector{Vector{QQFieldElem}}()
    for pᵢ in p
        s = [sign(evaluate(q, pᵢ)) for q in f]
        if !(s in signs)
            push!(signs, s)
            push!(p′, pᵢ)
        end
    end
    p′
end

function _sample_points(
    f::Vector{QQMPolyRingElem};
    distinct_signs::Bool=false,
    nr_thrds::Int=1,
    show_progress::Bool=false
)::Vector{Vector{QQFieldElem}}
    if any(is_zero, f)
        error("Cannot sample points for polynomials containing zero polynomial")
    end
    p = Vector{Vector{QQFieldElem}}([])
    if length(gens(parent(f[1]))) == 0
        p = _sample_points_0(f)
    elseif length(gens(parent(f[1]))) == 1
        p = _sample_points_1(f)
    else
        p = _sample_points_2(f, gens(parent(f[1])); nr_thrds=nr_thrds, show_progress=show_progress)
    end
    if distinct_signs
        p = _unique_by_sign(p, f)
    end
    return p
end
