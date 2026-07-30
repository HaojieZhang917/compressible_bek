#!/usr/bin/env python3
"""
Shooting method for incompressible von Kármán rotating disk base flow
with axial temperature variation under the Boussinesq approximation.
Outputs profiles compatible with Bone.py format.

Usage:
    python ShootBone.py          # uses TARGET_TW below
    python ShootBone.py 0.75     # specify Tw from command line
"""

import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ═══════════════════════════════════════
#  User settings
# ═══════════════════════════════════════
TARGET_TW = 1.0      # wall temperature ratio Tw=T0/Tinf
Pr        = 0.72      # Prandtl number
ZMAX      = 20.0      # domain height
N_OUT     = 2000      # output resolution

if len(sys.argv) > 1:
    TARGET_TW = float(sys.argv[1])

# ═══════════════════════════════════════
#  ODE system
# ═══════════════════════════════════════
def von_karman_ode(z, y):
    H, Fp, F, Gp, G, Tp, T = y
    return [
        -2.0 * F,                              # H'
        F**2 + H*Fp - (G - 1)**2 - T + 1,      # F'' = F²+HF'-(G-1)²-(T-1)
        Fp,                                     # F'
        2*F*G + H*Gp - 2*F,                    # G''
        Gp,                                     # G'
        Pr * H * Tp,                            # T''
        Tp,                                     # T'
    ]

# ═══════════════════════════════════════
#  Shooting: error at infinity
# ═══════════════════════════════════════
def shooting_error(guess, Tw):
    """guess = [F'(0), G'(0), T'(0)]  →  error = [F(∞), G(∞)-1, T(∞)-1]"""
    y0 = [0.0, guess[0], 0.0, guess[1], 0.0, guess[2], Tw]
    sol = solve_ivp(von_karman_ode, (0, ZMAX), y0,
                    method='RK45', rtol=1e-10, atol=1e-12, max_step=0.05)
    return [sol.y[2, -1], sol.y[4, -1] - 1.0, sol.y[6, -1] - 1.0]

# ═══════════════════════════════════════
#  Full integration with converged guess
# ═══════════════════════════════════════
def integrate_full(guess, Tw):
    """Returns (z, y) from full integration"""
    y0 = [0.0, guess[0], 0.0, guess[1], 0.0, guess[2], Tw]
    z_eval = np.linspace(0, ZMAX, N_OUT)
    sol = solve_ivp(von_karman_ode, (0, ZMAX), y0,
                    method='RK45', rtol=1e-10, atol=1e-12,
                    max_step=0.05, t_eval=z_eval)
    return sol.t, sol.y

# ═══════════════════════════════════════
#  Continuation solver
# ═══════════════════════════════════════
def solve_shooting(Tw_target, verbose=True):
    # Get accurate starting guess from solve_bvp at Tw=1.0
    from scipy.integrate import solve_bvp as sbvp
    z_bvp = np.linspace(0, ZMAX, 500)
    g0 = np.zeros((7, 500))
    g0[0] = -0.8845*(1-np.exp(-0.8*z_bvp))
    g0[1] =  0.5102*(1-0.8*z_bvp)*np.exp(-0.8*z_bvp)
    g0[2] =  0.5102*z_bvp*np.exp(-0.8*z_bvp)
    g0[3] =  0.6159*np.exp(-0.8*z_bvp)
    g0[4] =  1.0 - np.exp(-0.8*z_bvp)
    g0[5] =  0.0;  g0[6] = 1.0
    
    def bvp_bc(ya, yb, Tw):
        return [ya[0], ya[2], ya[4], ya[6]-Tw, yb[2], yb[4]-1, yb[6]-1]
    
    sol0 = sbvp(von_karman_ode, lambda ya, yb: bvp_bc(ya, yb, 1.0),
                z_bvp, g0, tol=1e-8, max_nodes=200000)
    
    # Refine with shooting at Tw=1.0
    guess_bvp = np.array([sol0.sol(0)[1], sol0.sol(0)[3], sol0.sol(0)[5]])
    guess = fsolve(lambda g: shooting_error(g, 1.0), guess_bvp,
                   xtol=1e-12, maxfev=100)
    
    if verbose:
        y0 = [0.0, guess[0], 0.0, guess[1], 0.0, guess[2], 1.0]
        s = solve_ivp(von_karman_ode, (0, ZMAX), y0,
                      method='RK45', rtol=1e-10, atol=1e-12, max_step=0.05)
        print(f"[1/2] Tw=1.0  H(∞)={s.y[0,-1]:.6f}  "
              f"F'(0)={guess[0]:.6f}  G'(0)={guess[1]:.6f}  ✅")
    
    if abs(Tw_target - 1.0) < 1e-10:
        return finalize(guess, 1.0, verbose)
    
    # Continuation: use solve_bvp solution as initial guess at each step
    cg = sol0.sol(z_bvp)
    direction = 1.0 if Tw_target > 1.0 else -1.0
    dTw = 0.002 if Tw_target > 1.0 else 0.002
    
    if verbose:
        wt = "热壁" if Tw_target > 1.0 else "冷壁"
        print(f"[2/2] 延拓 Tw → {Tw_target:.4f} (dTw={dTw}, {wt}) ...",
              end=" ", flush=True)
    
    for tw in np.arange(1.0 + direction*dTw, Tw_target + direction*dTw*0.5,
                        direction*dTw):
        # solve_bvp step as initial guess provider
        r = sbvp(von_karman_ode, lambda ya, yb: bvp_bc(ya, yb, tw),
                 z_bvp, cg, tol=1e-6, max_nodes=200000)
        if not r.success:
            if verbose: print(f"❌\n    solve_bvp failed at {tw:.4f}")
            return None, None, None, None, None
        cg = r.sol(z_bvp)
        
        # Shooting refinement with bvp-provided guess
        g_bvp = np.array([r.sol(0)[1], r.sol(0)[3], r.sol(0)[5]])
        try:
            gs = fsolve(lambda g: shooting_error(g, tw), g_bvp,
                        xtol=1e-12, maxfev=100)
        except Exception:
            gs = g_bvp  # fallback to bvp guess
        
        # Verify
        y0 = [0.0, gs[0], 0.0, gs[1], 0.0, gs[2], tw]
        s = solve_ivp(von_karman_ode, (0, ZMAX), y0,
                      method='RK45', rtol=1e-10, atol=1e-12, max_step=0.05)
        Hinf_s = s.y[0, -1]
        
        # Compare with solve_bvp
        Hinf_b = r.sol(ZMAX)[0]
        if abs(Hinf_s - Hinf_b) > 0.05:
            if verbose:
                print(f"❌\n    伪解: shoot H∞={Hinf_s:.4f} vs bvp H∞={Hinf_b:.4f}")
            return None, None, None, None, None
        
        guess = gs.copy()
    
    if verbose:
        y0 = [0.0, guess[0], 0.0, guess[1], 0.0, guess[2], tw]
        s = solve_ivp(von_karman_ode, (0, ZMAX), y0,
                      method='RK45', rtol=1e-10, atol=1e-12, max_step=0.05)
        print(f"✅")
        print(f"    H(∞)={s.y[0,-1]:.6f}  F'(0)={guess[0]:.6f}  "
              f"G'(0)={guess[1]:.6f}  T'(0)={guess[2]:.6f}")
    
    return finalize(guess, tw, verbose)

def finalize(guess, Tw, verbose):
    y0 = [0.0, guess[0], 0.0, guess[1], 0.0, guess[2], Tw]
    s = solve_ivp(von_karman_ode, (0, ZMAX), y0,
                  method='RK45', rtol=1e-10, atol=1e-12, max_step=0.05)
    Hinf = s.y[0, -1]
    Fp0  = guess[0]
    Gp0  = guess[1]
    Tp0  = guess[2]
    if verbose:
        print(f"✅")
        print(f"    H(∞)={Hinf:.6f}  F'(0)={Fp0:.6f}  "
              f"G'(0)={Gp0:.6f}  T'(0)={Tp0:.6f}")
    return guess, Hinf, Fp0, Gp0, Tp0

# ═══════════════════════════════════════
#  Save output
# ═══════════════════════════════════════
def save_output(Tw, guess):
    """Save full profile to .dat file"""
    z, y = integrate_full(guess, Tw)
    
    data = np.column_stack([
        z,
        y[0],   # H
        y[2],   # F
        y[4],   # G
        y[6],   # T
        y[1],   # dF
        y[3],   # dG
        y[5],   # dT
    ])
    
    tag = f"Tw{Tw:.4f}".replace('.', 'p')
    fname = f"baseflow_shoot_{tag}.dat"
    header = (f"# Shooting method, Tw={Tw:.4f}, Pr={Pr}\n"
              f"# columns: z H F G T dF dG dT")
    np.savetxt(fname, data, header=header, fmt="%.8e")
    print(f"\n📁 {fname}  ({N_OUT} pts)")
    return z, y

def plot_profile(Tw, guess):
    """Generate and save 4-panel profile plot"""
    z, y = integrate_full(guess, Tw)
    H, F, G, T = y[0], y[2], y[4], y[6]

    z_plot = z[z <= 10]
    idx = len(z_plot)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(f'von Karman Base Flow (Shooting)  Tw = {Tw:.4f}', fontsize=13)

    ax = axes[0, 0]
    ax.plot(z_plot, F[:idx], 'b-', lw=1.5)
    ax.set_xlabel('z'); ax.set_ylabel('F')
    ax.set_title('Radial velocity F(z)')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='gray', ls='--', lw=0.5)

    ax = axes[0, 1]
    ax.plot(z_plot, G[:idx], 'r-', lw=1.5)
    ax.set_xlabel('z'); ax.set_ylabel('G')
    ax.set_title('Azimuthal velocity G(z)')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=1, color='gray', ls='--', lw=0.5, label='G(inf)=1')
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.plot(z_plot, H[:idx], 'g-', lw=1.5)
    ax.set_xlabel('z'); ax.set_ylabel('H')
    ax.set_title('Axial velocity H(z)')
    ax.grid(True, alpha=0.3)
    Hinf_v = H[-1]
    ax.axhline(y=Hinf_v, color='green', ls=':', lw=0.8, label=f'H(inf)={Hinf_v:.4f}')
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    ax.plot(z_plot, T[:idx], 'm-', lw=1.5)
    ax.set_xlabel('z'); ax.set_ylabel('T')
    ax.set_title('Temperature T(z)')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=1.0, color='gray', ls='--', lw=0.5, label='T(inf)=1')
    ax.axhline(y=Tw, color='red', ls=':', lw=0.8, label=f'Tw={Tw:.4f}')
    ax.legend(fontsize=8)

    plt.tight_layout()
    tag = f"Tw{Tw:.4f}".replace('.', 'p')
    png = f"profile_shoot_{tag}.png"
    plt.savefig(png, dpi=150)
    plt.close()
    print(f"📊 {png}")

# ═══════════════════════════════════════
#  Main
# ═══════════════════════════════════════
if __name__ == "__main__":
    print("=" * 55)
    print("  Shooting method: von Kármán base flow")
    print(f"  Tw = {TARGET_TW:.4f},  Pr = {Pr},  zmax = {ZMAX}")
    print("=" * 55)
    
    result = solve_shooting(TARGET_TW)
    
    if result[0] is not None:
        guess, Hinf, Fp0, Gp0, Tp0 = result
        save_output(TARGET_TW, guess)
        plot_profile(TARGET_TW, guess)
        print(f"\n  壁面: F'(0)={Fp0:.6f}  G'(0)={Gp0:.6f}  T'(0)={Tp0:.6f}")
        print(f"  远场: H(∞) ={Hinf:.6f}")
    else:
        print(f"\n❌ 求解失败 (Tw={TARGET_TW:.4f})")
        sys.exit(1)
