# ================================
# Load libraries
# ================================
import sys as sys
import scikits.bvp_solver as scbs
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
font = {'family': 'serif', 'weight': 'medium', 'size': 40}
mpl.rc('font', **font)
mpl.rcParams['axes.linewidth'] = 5.0
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ==================================
# Parse user input
# ==================================
if len(sys.argv) == 1 or (len(sys.argv) == 2 and 'help' in sys.argv[2]):
    exit()
else:
    for arg in sys.argv:
     if arg[:4] == '--re':
        reynolds = float(arg.split('=')[1])
     elif arg[:4] == '--ka':
        kappa = float(arg.split('=')[1])
        if not abs(kappa) <= 1.0:
           exit()
     elif arg[:5] == '--ver':
         verbose = int(arg.split('=')[1])
         if verbose != 0 and verbose != 1 and verbose != 2:
             exit()
     elif arg[:4] == '--rk':
         rk = int(arg.split('=')[1])
         if rk != 2 and rk != 4 and rk != 6:
             exit()
     elif arg[:5] == '--tol':
        tol = float(arg.split('=')[1])
     elif arg[:5] == '--num':
        disks = int(arg.split('=')[1])
        if disks != 1 and disks != 2:
            exit()

# ================================
# Define the ODEs
# ================================
def oneDiskODE(z, y):
    # Y0 = H, Y1 = F', Y2 = F, Y3 = G', Y4 = G
    dydz = np.array([-
                     2.0*
                     y[2], kappa * 
                     (y[2] *
                      y[2] +
                         y[0] *
                         y[1] -
                         (y[4] - 
                          1.0)) -
                     (2.0 -
                         kappa - 
                         kappa**2) *
                     (y[4] -
                        1.0), y[1], kappa *
                     (2.0 * 
                         y[2] *
                         y[4] +
                         y[0] *
                         y[3]) +
                     (2.0 -
                          kappa -
                          kappa**2) *
                     y[2], y[3]])
    return dydz    

def oneDiskBC(y0, yinf):
    res0 = np.array([y0[0],
                     y0[2],
                      y0[4]])
    
    resinf = np.array([yinf[2],
                       yinf[4] - 1.0])
    return (res0, resinf)


def twoDiskODE(z, y, k):
# Y0 = H, Y1 = F', Y2 = F, Y3 = G', Y4 = G
    dydz = np.array([-
                         2.0 *
                         np.sqrt(reynolds) *
                         y[2], reynolds *
                         (y[0] *
                             y[1] /
                             np.sqrt(reynolds) +
                          y[2] *
                          y[2] -
                          (y[4] +
                           1.0)**2 +
                          k), y[1], reynolds *
                         (y[0] *
                          y[3] /
                          np.sqrt(reynolds) +
                           2.0 *
                          y[2] *
                          (y[4] +
                            1)), y[3]])
    return dydz

def twoDiskBC(y0, yinf, k):
    res0 = np.array([y0[0],
                     y0[2],
                     y0[4]])
            
    resinf = np.array([yinf[0],
                       yinf[2],
                       yinf[4] + 1.0])
    
    return (res0, resinf)

     # ============================
# define the problem
# ============================

if disks == 1:
    problem = scbs.ProblemDefinition(num_ODE=5,
                                     num_parameters=0,
                                     num_left_boundary_conditions=3,
                                     boundary_points=(0.0, np.sqrt(reynolds)),
                                     function=oneDiskODE,
                                     boundary_conditions=oneDiskBC)
    
    if kappa == -1:
         solguess = [1.0, 0.0, 0.0, 0.0, 1.0]
    elif kappa == 1:
         solguess = [1.2, 0.0, 0.0, 0.0, 1.0]
    else:
         solguess = [1.0, 0.0, 0.0, 0.0, 1.0]
    solguess = np.array(solguess)
    
solution = scbs.solve(
    problem,
    solution_guess=solguess,
    initial_mesh=np.linspace(
        0.0,
        np.sqrt(reynolds),
        2000),
    method=int(rk),
    tolerance=float(tol),
    max_subintervals=1000000,
    trace=int(verbose)
)

vecZ = np.linspace(0.0, np.sqrt(reynolds), 2000)


if disks == 2:
    solution = None
    solguess = np.array([0.0, 0.0, 0.0, 0.0, 0.3])

    problem = scbs.ProblemDefinition(num_ODE=5,
                                     num_parameters=1,
                                     num_left_boundary_conditions=3,
                                     boundary_points=(0.0, 1.0),
                                     function=twoDiskODE,
                                     boundary_conditions=twoDiskBC)

    solution = scbs.solve(problem,
                          solution_guess=solguess,
                          parameter_guess=0.09,
                          initial_mesh=np.linspace(0.0, 1.0, 2000),
                          method=int(rk),
                          tolerance=float(tol),
                          max_subintervals=1000000,
                          trace=int(verbose))

    vecZ = np.linspace(0.0, 1.0, 2000)

# ================================
# Plot the solution
# ================================

if disks == 2:
    print (solution.parameters)

vecH, vecFp, vecF, vecGp, vecG = solution(vecZ)

fig = plt.figure(figsize=(30, 20), dpi=300)
ax = plt.subplot(111)
plt.plot(vecH, vecZ, '-k', lw=4, label='Axial')
plt.plot(vecG, vecZ, '-b', lw=4, label='Azimutal')
plt.plot(vecF, vecZ, '-r', lw=4, label='Radial')
plt.legend(loc='best')
plt.savefig('meanflow.png', bbox_inches='tight')
plt.close()


# ================================
# Export the solution
# ================================
fid = open('mean.dat', 'w')
for idx in xrange(vecZ.size):
    fid.write('%+.12f %+.12f %+.12f %+.12f %+.12f %+.12f %+.12f\n' %
             (vecF[idx], vecG[idx], vecH[idx], vecFp[idx], vecGp[idx], -
              2.0 *
              np.sqrt(reynolds) *
                 vecF[idx], vecZ[idx]))
fid.close()
