# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Julia 1.10.2
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# # KP方程式のN-soliton解
#
# See https://arxiv.org/abs/1208.2904

# %%
using Combinatorics
using ForwardDiff
using Plots
default(fmt=:png, colorbar=false)
color = :gist_earth
using SymPy

nan2zero(x) = isnan(x) ? zero(x) : x

_eta(u, v, x, y, t) = (u - v) * x + (u^2 - v^2) * y + (u^3 - v^3) * t

function _soliton_exp_factor(u, v, c, i, x, y, t)
    n = length(i)
    C = prod(c[i[k]] for k in 1:n; init=1.0)
    E = sum(_eta(u[i[k]], v[i[k]], x, y, t) for k in 1:n; init=0.0)
    C * exp(E)
end

function _soliton_tau_term(u, v, c, i, x, y, t)
    n = length(i)
    A = prod((u[i[k]] - u[i[l]]) * (v[i[k]] - v[i[l]]) /
        ((u[i[k]] - v[i[l]]) * (v[i[k]] - u[i[l]])) for k in 1:n for l in k+1:n; init=1.0)
    A * _soliton_exp_factor(u, v, c, i, x, y, t)
end

function make_soliton_tau(u, v, c)
    N = length(c)
    tau(x, y, t) = sum(_soliton_tau_term(u, v, c, i, x, y, t) for n in 0:N for i in combinations(1:N, n))
    tau
end

make_soliton_tau(u, v) = make_soliton_tau(u, v, ones(length(u)))
make_soliton_tau_AB(k, P, c) = make_soliton_tau((P+k)/2, (P-k)/2, c)
make_soliton_tau_AB(k, P) = make_soliton_tau_AB(k, P, ones(length(k)))

first_derivative(f, x) = ForwardDiff.derivative(f, x)
second_derivative(f, x) = ForwardDiff.derivative(x -> ForwardDiff.derivative(f, x), x)

function tau2sol(tau)
    tau_x(x, y, t) = first_derivative(x -> tau(x, y, t), x)
    tau_xx(x, y, t) = second_derivative(x -> tau(x, y, t), x)
    sol(x, y, t) = 2(tau_xx(x, y, t) / tau(x, y, t) - (tau_x(x, y, t) / tau(x, y, t))^2)
    sol
end

function _make_soliton_solution(u, v, c)
    _tau =  make_soliton_tau(u, v, c)
    @syms x y t
    τ = _tau(x, y, t)
    τ_x = sympy.diff(τ, x)
    τ_xx = sympy.diff(τ_x, x)
    tau = eval(Meta.parse("(x, y, t) -> $τ"))
    tau_x = eval(Meta.parse("(x, y, t) -> $τ_x"))
    tau_xx = eval(Meta.parse("(x, y, t) -> $τ_xx"))
    sol(x, y, t) = 2(tau_xx(x, y, t) / tau(x, y, t) - (tau_x(x, y, t) / tau(x, y, t))^2)
    (; sol, tau, tau_x, tau_xx)
end

make_soliton_solution(u, v, c) = _make_soliton_solution(u, v, c).sol
make_soliton_solution(u, v) = make_soliton_solution(u, v, ones(length(u)))
make_soliton_solution_AB(k, P, c) = make_soliton_solution((P+k)/2, (P-k)/2, c)
make_soliton_solution_AB(k, P) = make_soliton_solution((P+k)/2, (P-k)/2, ones(length(k)))

function show_tau_and_gif_sol(sol, u, v, xs, ys, ts; fn="tmp.gif", fps=20, disp=false)
    @syms x y t
    tau = make_soliton_tau(u, v)
    if disp
        print("tau(x, y, t) = ")
        display(tau(x, y, t))
    else
        @show tau(x, y, t)
    end
    flush(stdout)
    
    solmax = maximum(nan2zero(sol(x, y, t)) for x in xs, y in ys, t in ts)
    anim = @animate for t in ts
        surface(xs, ys, (x, y) -> sol(x, y, t); camera=(15, 80), colorbar=false, color)
        plot!(xguide="x", yguide="y", zguide="u")
        plot!(size=(600, 600))
        plot!(zlim=(-0.02solmax, 1.05solmax))
        plot!(margin=-10Plots.mm)
    end

    gif(anim, fn; fps)    
end

show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="tmp.gif", fps=20, disp=false) =
    show_tau_and_gif_sol(sol, (P+k)/2, (P-k)/2, xs, ys, ts; fn, fps, disp)

# %%
u, v = [1.0], [-0.8]

@syms x y t
tau = make_soliton_tau(u, v)
2SymPy.diff(log(tau(x, y, t)), x, x).simplify() |> display
sol = make_soliton_solution(u, v)
sol(x, y, t).simplify() |> display

sol = make_soliton_solution(u, v)
xs = range(-25, 25, 251)
ys = range(-25, 25, 251)
ts = range(-25, 25, 101)
show_tau_and_gif_sol(sol, u, v, xs, ys, ts; fn="1-soliton.gif", disp=true)

# %%
u, v = [1.0, 1.5], [-0.8, -1.0]
sol = make_soliton_solution(u, v)
xs = range(-25, 25, 251)
ys = range(-25, 25, 251)
ts = range(-25, 25, 101)
show_tau_and_gif_sol(sol, u, v, xs, ys, ts; fn="2-soliton.gif", disp=true)

# %%
k, P = [0.5, 0.5], [0.66, -0.66]
sol = make_soliton_solution_AB(k, P)
xs = range(-50, 50, 251)
ys = range(-50, 50, 251)
ts = range(-200, 200, 101)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="AB_FIG1.gif", disp=true)

# %%
k, P = [0.5, 1.0], [0.75, 0.25]
sol = make_soliton_solution_AB(k, P)
xs = range(-50, 50, 251)
ys = range(-75, 75, 251)
ts = range(-200, 200, 101)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="AB_FIG2.gif", disp=true)

# %%
k, P = [0.5, 0.5], [-0.26, 0.75]
sol = make_soliton_solution_AB(k, P)
xs = range(-50, 50, 251)
ys = range(-75, 75, 251)
ts = range(-200, 200, 101)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="AB_FIG3.gif", disp=true)

# %%
k, P = [0.5, 0.5], [-0.5+1e-10, 0.5+1e-9]
sol = make_soliton_solution_AB(k, P)
xs = range(-50, 50, 251)
ys = range(-75, 75, 251)
ts = range(-200, 200, 101)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="AB_FIG4.gif", disp=true)

# %%
k, P = [1.0, 0.5], [0.5-1e-7, 0.0]
sol = make_soliton_solution_AB(k, P)
xs = range(-100, 100, 251)
ys = range(-150, 150, 251)
ts = range(-400, 400, 101)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="AB_FIG5.gif", disp=true)

# %%
k, P = [1.0, 2.0, 3.0], [-0.333, -0.667, -1.667]
sol = make_soliton_solution_AB(k, P)
xs = range(-25, 25, 251)
ys = range(-25, 50, 251)
ts = range(-25, 25, 201)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="AB_FIG6.gif")

# %%
k, P = [1.0, 2.0, 3.0], [-0.333, -0.667, -1.66]
sol = make_soliton_solution_AB(k, P)
xs = range(-25, 25, 251)
ys = range(-25, 50, 251)
ts = range(-25, 25, 201)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="3-soliton_a.gif")

# %%
k, P = [1.0, 2.0, 3.0], [-0.333, -0.667, -1.5]
sol = make_soliton_solution_AB(k, P)
xs = range(-25, 25, 251)
ys = range(-25, 50, 251)
ts = range(-25, 25, 201)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="3-soliton_b.gif")

# %%
k, P = [0.5, 1.0, 2.3], [0.75, 0.25, -1.0]
sol = make_soliton_solution_AB(k, P)
xs = range(-25, 25, 251)
ys = range(-25, 50, 251)
ts = range(-25, 25, 201)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="3-soliton_c.gif")

# %%
k, P = [0.5, 1.0, 1.5], [0.75, 0.25, -0.25]
sol = make_soliton_solution_AB(k, P)
xs = range(-25, 25, 251)
ys = range(-25, 50, 251)
ts = range(-25, 50, 201)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="3-soliton_d.gif")

# %%
k, P = [0.5, 1.0, 2.25], [0.75, 0.25, -1.0]
sol = make_soliton_solution_AB(k, P)
xs = range(-25, 25, 251)
ys = range(-25, 50, 251)
ts = range(-25, 25, 201)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="3-soliton_e.gif")

# %%
