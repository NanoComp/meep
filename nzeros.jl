using Gridap, Gridap.Geometry, Gridap.Fields, GridapGmsh
using LinearAlgebra

using ChainRulesCore, Zygote
import ChainRulesCore: rrule
using DSP
using NLopt
using CairoMakie, GridapMakie
using KrylovKit

using JLD

n = 3
filtertype = "stop"
profile = "ellip"
pass = 0.25
stop = 25
band_l = 0.98
band_r = 1.02

#freqs = [0.983147947076 + 0.004516989107im, 0.999545615922 + 0.022551312414im, 1.016712594700 + 0.004671199008im]
function calc_poles(n,pass,stop,band_l, band_r,profile,filtertype)
    fs = max(2*band_r, 10)
    if filtertype == "bandpass"
        responsetype = Bandpass(band_l, band_r, fs=fs)
    elseif filtertype == "stop"
        responsetype = Bandstop(band_l, band_r, fs=fs)
    end
    if profile == "butter"
        designmethod = Butterworth(n)
    elseif profile == "cheby"
        designmethod = Chebyshev1(n, pass)
    elseif profile == "icheby"
        designmethod = Chebyshev2(n, stop)
    elseif profile == "ellip"
        designmethod = Elliptic(1n, pass, stop)
    end

    fil=analogfilter(responsetype, designmethod)
    p = sort(-fil.p*0.5im*fs, by=real)
    return p[n+1:end]
end

freqs = calc_poles(n, pass, stop, band_l, band_r, profile, filtertype)

Γ = maximum(imag.(freqs))
C_freqs = 1 .+ [-10, -8, 8, 10]*Γ

numfreqs = length(freqs)
num_C_freqs = length(C_freqs)

function calc_Sbar(freq, poles, σ)#poles and σ are row vec
    N = length(poles)
    sPoles = 1im*poles;
    D = [ones(1,N); σ]
    M = (D'*D)./(sPoles' .+ sPoles);
    K = conj(D/M); #% This is K without C.

    #%% Calculate B=\bar{S}, its derivative and its pole residues.
    Nw = length(freq); #% number of calculation frequencies
    P = 2#number of ports
    B = zeros(ComplexF64, P,P,Nw);

    #if ~isempty(wPoles_abs), sPoles = 1i*wPoles_abs; end % Pole perturbation due to absorption.

    for p = 1:P
        for q = 1:P
            rBpq = D[p,:].*K[q,:];
            for i1 = 1:Nw
                wi = freq[i1];
                B[p,q,i1] = (p==q) + sum(rBpq./(1im*wi .- sPoles));
                #dB(p,q,i1) = 1i*sum(rBpq./(1i*wi-sPoles).^2);
            end
        end
    end
    return B
end

σs = [-(-1)^i for i=1:numfreqs]
Sbar = calc_Sbar(C_freqs, freqs',  σs')
Sbar_inv = [inv(Sbar[:,:,i]) for i=1:num_C_freqs]

λs = 1 ./ freqs
ks = 2*π * freqs
C_ks = 2*π * C_freqs
L = 2+2.15*0.40362   # Width of the numerical cell (excluding PML) (μm)
H = 0.2   # Height of the numerical cell (excluding PML) (μm)
dpml = 0.5   # Thickness of the PML (μm)

n_Si = sqrt(11) # Silicon refractive index
n_air = 1    # Air refractive index
μ = 1        # Magnetic permeability
R = 1e-20
LH=(L, H)

model = GmshDiscreteModel("/Users/mochen/desktop/research/Filter/2d_short.msh")
order = 1
reffe = ReferenceFE(lagrangian, Float64, order)
V = TestFESpace(model, reffe, vector_type = Vector{ComplexF64},dirichlet_tags=["DirichletEdge", "Dnode"])
U = V   # mathematically equivalent to TrialFESpace(V,0)
uh = FEFunction(V,rand(ComplexF64,num_free_dofs(V)))
all_phys_params = [(; k=ks[i], σ = -(-1)^i, n_Si, n_air, μ, R, dpml, LH) for i=1:numfreqs]
all_phys_params_C = [(; k=C_ks[i],σ = 1, n_Si, n_air, μ, R, dpml, LH) for i=1:num_C_freqs]

degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

Γ_ls = BoundaryTriangulation(model; tags = ["L_Src"]) # Left Source line
dΓ_ls = Measure(Γ_ls, degree)
Γ_rs = BoundaryTriangulation(model; tags = ["R_Src"]) # Right Source line
dΓ_rs = Measure(Γ_rs, degree)

Γ_s = BoundaryTriangulation(model; tags = ["L_Src", "R_Src"]) # Right Source line
dΓ_s = Measure(Γ_s, degree)

Ω_d = Triangulation(model, tags="Design")
dΩ_d = Measure(Ω_d, degree)

Ω_l = Triangulation(model, tags="Left")
dΩ_l = Measure(Ω_l, degree)
Ω_r = Triangulation(model, tags="Right")
dΩ_r = Measure(Ω_r, degree)


p_reffe = ReferenceFE(lagrangian, Float64, 0)
Q = TestFESpace(Ω_d, p_reffe, vector_type = Vector{Float64})
P = Q
np = num_free_dofs(P)

pf_reffe = ReferenceFE(lagrangian, Float64, 1)
Qf = TestFESpace(Ω_d, pf_reffe, vector_type = Vector{Float64})
Pf = Qf

fem_params = (; V, U, Q, P, Qf, Pf, np, Ω, dΩ, dΩ_d, dΩ_l, dΩ_r, dΓ_ls, dΓ_rs, dΓ_s)

labels = get_face_labeling(model)
dimension = 2
tags = get_face_tag(labels,dimension)

const l_tag = get_tag_from_name(labels,"Left")
const r_tag = get_tag_from_name(labels,"Right")

function subtract_outgoing(u,u_ls, u_rs,tag)
  if tag == l_tag
    return u - u_ls
  elseif tag == r_tag
    return u - u_rs
  else
    return u
  end
end

function s_PML(x; phys_params)
    σ = -3 / 4 * log(phys_params.R) / phys_params.dpml / phys_params.n_air
    u = abs.(Tuple(x)).-phys_params.LH./2
    return @. ifelse(u > 0,  1 + (1im * σ / phys_params.k) * (u / phys_params.dpml)^2, $(1.0+0im))
end

function ds_PML(x; phys_params)
    σ = -3 / 4 * log(phys_params.R) / phys_params.dpml / phys_params.n_air
    u = abs.(Tuple(x)).-phys_params.LH./2
    ds = @. ifelse(u > 0, (2im * σ / phys_params.k) * (1 / phys_params.dpml)^2 * u, $(0.0+0im))
    return ds.*sign.(Tuple(x))
end

struct Λ{PT} <: Function
    phys_params::PT
end

function (Λf::Λ)(x)
    s_x,s_y = s_PML(x; Λf.phys_params)
    return VectorValue(1/s_x, 1/s_y)
end

Fields.∇(Λf::Λ) = x -> TensorValue{2, 2, ComplexF64}(-(Λf(x)[1])^2 * ds_PML(x; Λf.phys_params)[1], 0, 0, -(Λf(x)[2])^2 * ds_PML(x; Λf.phys_params)[2])

r = 0.002               # Filter radius
β = 32.0                    # β∈[1,∞], threshold sharpness
η = 0.5                     # η∈[0,1], threshold center

a_f(r, u, v) = r^2 * (∇(v) ⋅ ∇(u))
function Filter(p0; r, fem_params)
    ph = FEFunction(fem_params.P, p0)
    op = AffineFEOperator(fem_params.Pf, fem_params.Qf) do u, v
        ∫(a_f(r, u, v))fem_params.dΩ_d + ∫(v * u)fem_params.dΩ_d, ∫(v * ph)fem_params.dΩ_d
      end
    pfh = solve(op)
    return get_free_dof_values(pfh)
end

function Threshold(pfh; β, η)
    return ((tanh(β * η) + tanh(β * (pfh - η))) / (tanh(β * η) + tanh(β * (1.0 - η))))
end

ξd(p, n_air, n_Si)= (n_air + (n_Si - n_air) * p)^2 - n_air^2 # in the design region
a_base(u, v; phys_params) =  ((∇ .* (Λ(phys_params) * v)) ⊙ (Λ(phys_params) .* ∇(u))) - (phys_params.n_air^2) * (phys_params.k^2 * phys_params.μ * (v * u))
a_design(u, v, pth; phys_params) = - ((p -> ξd(p, phys_params.n_air, phys_params.n_Si)) ∘ pth) * (phys_params.k^2 * phys_params.μ * (v * u))
a_design0(u, v, pth; phys_params) = ((p -> ξd(p, phys_params.n_air, phys_params.n_Si)) ∘ pth) * (v * u)

function MatrixA(pth; phys_params, fem_params)
    A_mat = assemble_matrix(fem_params.U, fem_params.V) do u, v
        ∫(a_base(u, v; phys_params))fem_params.dΩ + ∫(a_design(u, v, pth; phys_params))fem_params.dΩ_d
    end
    return lu(A_mat)
end

function MatrixB(pth; phys_params, fem_params)
    B_mat = assemble_matrix(fem_params.U, fem_params.V) do u, v
        ∫(v * u)fem_params.dΩ + ∫(a_design0(u, v, pth; phys_params))fem_params.dΩ_d
    end
    return B_mat
end

p0 = zeros(fem_params.np)  # Here we make p=0 everywhere just for illustration purpose
pf_vec = Filter(p0;r, fem_params)
pfh = FEFunction(fem_params.Pf, pf_vec)
pth = (pf -> Threshold(pf; β, η)) ∘ pfh

u_vecs_ls = Vector{Vector{ComplexF64}}()
u_vecs_rs = Vector{Vector{ComplexF64}}()
b_ls = assemble_vector(v->(∫(v)fem_params.dΓ_ls), fem_params.V)
b_rs = assemble_vector(v->(∫(v)fem_params.dΓ_rs), fem_params.V)
for i=1:numfreqs
    A_mat_zero_i = MatrixA(pth; phys_params=all_phys_params[i], fem_params)
    u_vec_ls_i = A_mat_zero_i \ b_ls
    u_vec_rs_i = A_mat_zero_i \ (all_phys_params[i].σ * b_rs)
    push!(u_vecs_ls, u_vec_ls_i)
    push!(u_vecs_rs, u_vec_rs_i)
end
incidences = [(; u_ls = u_vecs_ls[i], u_rs = u_vecs_rs[i]) for i=1:numfreqs]

function szeros(pf_vec; freqi, r, β, η, fem_params)
    pfh = FEFunction(fem_params.Pf, pf_vec)
    pth = (pf -> Threshold(pf; β, η)) ∘ pfh

    phys_params_i = (; k=2*π*freq, σ = 1, n_Si, n_air, μ, R, dpml, LH)
    A_mat_zero_i = MatrixA(pth; phys_params=phys_params_i, fem_params)
    u_vec_ls_i = A_mat_zero_i \ b_ls#left incidence
    u_vec_rs_i = A_mat_zero_i \ (phys_params_i.σ * b_rs)#right incidence

    A_mat = MatrixA(pth; phys_params=phys_params_i, fem_params)
    u_vec = A_mat \ (b_ls + phys_params_i.σ * b_rs)
    Sm2 = dot(b_rs, u_vec-u_vec_rs_i)
    Sm1 = dot(b_ls, u_vec-u_vec_ls_i)
    return Sm1, Sm2
end






function C21(pf_vec; Ci, β, η, fem_params)
    pfh = FEFunction(fem_params.Pf, pf_vec)
    pth = (pf -> Threshold(pf; β, η)) ∘ pfh
    pzero = FEFunction(fem_params.Pf, zeros(fem_params.np))
    A0 = MatrixA(pzero; phys_params=all_phys_params_C[Ci], fem_params)
    u_vec0 = A0 \ b_ls

    A_mat_C = MatrixA(pth; phys_params=all_phys_params_C[Ci], fem_params)
    u_vec_C = A_mat_C \ b_ls

    Sm2 = dot(b_rs, u_vec_C)
    Sp1 = dot(b_ls, u_vec0)
    Sm1 = dot(b_ls, u_vec_C-u_vec0)
    C21 = Sbar_inv[Ci][1,1] * (Sm1/Sp1) + Sbar_inv[Ci][1,2] * (Sm2/Sp1)#C11
    #C21 = Sbar_inv[Ci][2,1] * (Sm1/Sp1) + Sbar_inv[Ci][2,2] * (Sm2/Sp1)
    return abs2(C21)
end


function MatrixOf_l(pth; fem_params)
    return assemble_matrix(fem_params.U, fem_params.V) do u, v
        ∫(u ⋅ v)fem_params.dΩ_l
    end
end

function MatrixOf_r(pth; fem_params)
    return assemble_matrix(fem_params.U, fem_params.V) do u, v
        ∫(u ⋅ v)fem_params.dΩ_r
    end
end

NO_FIELDS = ZeroTangent()
Dptdpf(pf, β, η) = β * (1.0 - tanh(β * (pf - η))^2) / (tanh(β * η) + tanh(β * (1.0 - η)))
Dξdpf(pf, n_air, n_Si, β, η)= 2 * (n_Si - n_air) * (n_air + (n_Si - n_air) * Threshold(pf; β, η)) * Dptdpf(pf, β, η)
DAdpf(u, v, pfh; phys_params, β, η) = -  (phys_params.k^2 * phys_params.μ * (v * ((p -> Dξdpf(p, phys_params.n_air, phys_params.n_Si, β, η)) ∘ pfh) * u))


function gf_pf(pf_vec; β, η, phys_params, fem_params, incidence)
    pfh = FEFunction(fem_params.Pf, pf_vec)
    pth = (pf -> Threshold(pf; β, η)) ∘ pfh
    A_mat = MatrixA(pth; phys_params, fem_params)
    b_vec = assemble_vector(v->(∫(v)fem_params.dΓ_ls + phys_params.σ * ∫(v)fem_params.dΓ_rs), fem_params.V)
    u_vec = A_mat \ b_vec
    O_mat_l = MatrixOf_l(pth;fem_params)
    O_mat_r = MatrixOf_r(pth;fem_params)
    real((u_vec - incidence.u_ls)' * O_mat_l * (u_vec - incidence.u_ls) + (u_vec - incidence.u_rs)' * O_mat_r * (u_vec - incidence.u_rs))
end

function rrule(::typeof(gf_pf), pf_vec; β, η, phys_params, fem_params, incidence)
    function U_pullback(dgdg)
      NO_FIELDS, dgdg * Dgfdpf(pf_vec; β, η, phys_params, fem_params, incidence)
    end
    gf_pf(pf_vec; β, η, phys_params, fem_params, incidence), U_pullback
end

function Dgfdpf(pf_vec; β, η, phys_params, fem_params, incidence)
    pfh = FEFunction(fem_params.Pf, pf_vec)
    pth = (pf -> Threshold(pf; β, η)) ∘ pfh
    A_mat = MatrixA(pth; phys_params, fem_params)
    b_vec = assemble_vector(v->(∫(v)fem_params.dΓ_ls + phys_params.σ * ∫(v)fem_params.dΓ_rs), fem_params.V)
    u_vec = A_mat \ b_vec
    O_mat_l = MatrixOf_l(pth;fem_params)
    O_mat_r = MatrixOf_r(pth;fem_params)

    uh = FEFunction(fem_params.U, u_vec)
    w_vec =  A_mat' \ (O_mat_l * (u_vec - incidence.u_ls) + O_mat_r * (u_vec - incidence.u_rs))
    wconjh = FEFunction(fem_params.U, conj(w_vec))
    uconjh = FEFunction(fem_params.U, conj(u_vec))

    l_temp(dp) = ∫(real(-2*DAdpf(uh, wconjh, pfh; phys_params, β, η)) * dp)fem_params.dΩ_d
    dgfdpf = assemble_vector(l_temp, fem_params.Pf)
    return dgfdpf
end

function rrule(::typeof(C21), pf_vec; Ci, β, η, fem_params)
    function U_pullback(dgdg)
      NO_FIELDS, dgdg * DC21(pf_vec; Ci, β, η, fem_params)
    end
    C21(pf_vec; Ci, β, η, fem_params), U_pullback
end

function DC21(pf_vec; Ci, β, η, fem_params)
    pfh = FEFunction(fem_params.Pf, pf_vec)
    pth = (pf -> Threshold(pf; β, η)) ∘ pfh
    pzero = FEFunction(fem_params.Pf, zeros(fem_params.np))
    A0 = MatrixA(pzero; phys_params=all_phys_params_C[Ci], fem_params)
    u_vec0 = A0 \ b_ls
    A_mat_C = MatrixA(pth; phys_params=all_phys_params_C[Ci], fem_params)
    u_vec_C = A_mat_C \ b_ls
    Sm2 = dot(b_rs, u_vec_C)
    Sp1 = dot(b_ls, u_vec0)
    Sm1 = dot(b_ls, u_vec_C-u_vec0)
    C21 = Sbar_inv[Ci][1,1] * (Sm1/Sp1) + Sbar_inv[Ci][1,2] * (Sm2/Sp1)

    uh = FEFunction(fem_params.U, u_vec_C)
    rhs = C21 .* conj(Sbar_inv[Ci][1,1] * (b_ls/Sp1) .+ Sbar_inv[Ci][1,2] * (b_rs/Sp1))
    w_vec =  A_mat_C' \ rhs
    wconjh = FEFunction(fem_params.U, conj(w_vec))
    uconjh = FEFunction(fem_params.U, conj(u_vec_C))

    l_temp(dp) = ∫(real(-2*DAdpf(uh, wconjh, pfh; phys_params=all_phys_params_C[Ci], β, η)) * dp)fem_params.dΩ_d
    dc21dpf = assemble_vector(l_temp, fem_params.Pf)
    return dc21dpf
end

function pf_p0(p0; r, fem_params)#filter
    pf_vec = Filter(p0; r, fem_params)
    pf_vec
end

function rrule(::typeof(pf_p0), p0; r, fem_params)#rrule of filter
  function pf_pullback(dgdpf)
    NO_FIELDS, Dgdp(dgdpf; r, fem_params)
  end
  pf_p0(p0; r, fem_params), pf_pullback
end

function Dgdp(dgdpf; r, fem_params)#adjoint of filter
    Af = assemble_matrix(fem_params.Pf, fem_params.Qf) do u, v
        ∫(a_f(r, u, v))fem_params.dΩ_d + ∫(v * u)fem_params.dΩ_d
    end
    wvec = Af' \ dgdpf
    wh = FEFunction(fem_params.Pf, wvec)
    l_temp(dp) = ∫(wh * dp)fem_params.dΩ_d
    return assemble_vector(l_temp, fem_params.P)
end

function gf_p(p0::Vector; r, β, η, phys_params, fem_params, incidence)
    pf_vec = pf_p0(p0; r, fem_params)
    gf_pf(pf_vec; β, η, phys_params, fem_params, incidence)
end


function sum_Of(p0::Vector, grad::Vector; r, β, η, fem_params)
    if length(grad) > 0
        dgdp, = Zygote.gradient(p -> gf_p(p; r, β, η, phys_params=all_phys_params[1], fem_params, incidence=incidences[1]), p0)
        for i=2:numfreqs
            dgdpi, = Zygote.gradient(p -> gf_p(p; r, β, η, phys_params=all_phys_params[i], fem_params, incidence=incidences[i]), p0)
            dgdp = dgdp .+ dgdpi
        end
        grad[:] = dgdp
    end
    gvalues = zeros(numfreqs)
    for i = 1:numfreqs
        gvalues[i] = gf_p(p0::Vector; r, β, η, phys_params=all_phys_params[i], fem_params, incidence=incidences[i])
    end
    push!(history, gvalues)
    return sum(gvalues) - 1e-8
end

function C21_0(p0::Vector; Ci, r, β, η, fem_params)
    pf_vec = pf_p0(p0; r, fem_params)
    C21(pf_vec; Ci, β, η, fem_params)
end

function constraint_C(p0::Vector, grad::Vector; Ci, r, β, η, fem_params)
    if length(grad) > 0
        dcdp, = Zygote.gradient(p -> C21_0(p; Ci, r, β, η, fem_params), p0)
        grad[:] = dcdp
    end
    cvalue = C21_0(p0::Vector; Ci, r, β, η, fem_params)
    push!(C_history[Ci], cvalue)
    return cvalue - 0.02
end

function constraint_allC(p0::Vector, grad::Vector; r, β, η, fem_params)
    dcdp, = Zygote.gradient(p -> C21_0(p; Ci=1, r, β, η, fem_params), p0)
    for Ci = 2:num_C_freqs
        dcdp_i, = Zygote.gradient(p -> C21_0(p; Ci, r, β, η, fem_params), p0)
        dcdp = dcdp .+ dcdp_i
    end
    grad[:] = dcdp
    cvalue = C21_0(p0::Vector; Ci=1, r, β, η, fem_params)
    push!(C_history[1], cvalue)
    for Ci = 2:num_C_freqs
        cvalue_i = C21_0(p0::Vector; Ci, r, β, η, fem_params)
        push!(C_history[Ci], cvalue_i)
        cvalue += cvalue_i
    end
    return cvalue #- 0.05
end


ub = ones(fem_params.np)
lb = zeros(fem_params.np)

history = Vector{Vector{Float64}}()
C_history = [Vector{Float64}() for i=1:num_C_freqs]

function gf_p_optimize(p_init; r, β, η, TOL = 1e-12, MAX_ITER = 1500, fem_params)
    ##################### Optimize #################
    opt = Opt(:LD_MMA, fem_params.np)
    opt.lower_bounds = lb
    opt.upper_bounds = ub
    opt.ftol_rel = TOL
    opt.maxeval = MAX_ITER
    opt.min_objective = (x, grad) -> constraint_allC(x, grad; r, β, η, fem_params)
    #opt.min_objective = (x, grad) -> sum_Of(x, grad; r, β, η,  fem_params)
    #for ci = 1:num_C_freqs
    inequality_constraint!(opt, (x, grad) -> sum_Of(x, grad; r, β, η, fem_params), 1e-8)
    #end
    opt.params["inner_maxeval"] = 10

    (g_opt, p_opt, ret) = optimize(opt, p_init)
    @show g_opt
    @show numevals = opt.numevals # the number of function evaluations
    return g_opt, p_opt
end



psaved = load("/Users/mochen/desktop/research/Filter/2d.jld", "arr")
d = load("p_opt_init.jld")
p_opt = d["p_stop_opt"]
#p_opt = fill(0.4, fem_params.np)   # Initial guess
#p_opt[1] = 0
#p_opt[2:fem_params.np+1]=psaved
β_list = [4.0, 8.0, 16.0, 32.0]
