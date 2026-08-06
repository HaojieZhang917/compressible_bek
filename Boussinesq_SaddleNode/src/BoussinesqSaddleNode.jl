module BoussinesqSaddleNode

using CSV
using Dates
using DelimitedFiles
using JSON3
using LinearAlgebra
using NPZ
using Printf
using Statistics

include("Utils.jl")
include("VonKarman.jl")
include("RotorStator.jl")

export project_root, parse_cli, read_csv_rows, write_csv_rows, write_json,
       chebyshev_operators, barycentric_interpolate, linear_interpolate,
       locate_quadratic_extrema, bracketed_roots, simpson_uniform,
       VKProfile, VKBranch, vk_state, solve_vk_isothermal, solve_vk_fixed_h,
       trace_vk_h_branch, vk_turning_points, vk_roots_at_tw,
       vk_similarity_eigenvalues, vk_newton_diagnostics, vk_shoot,
       RotorConfig, RotorProfile, RotorBranch, solve_rotor_isothermal,
       continue_rotor_pressure, continue_rotor_temperature,
       solve_rotor_fixed_temperature, rotor_turning_points,
       rotor_roots_at_temperature, rotor_state, rotor_validation, pressure_grid,
       clean_name, quote_tecplot, numeric_value, format_number

end
