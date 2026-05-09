"""
Initial condition generators for the 2-D WENO Euler solver.

Modified from Manuel A. Diaz's WENO code (github.com/wme7/WENO).

Two classic test problems are provided:

1.  Double Mach Reflection (DMR) — Woodward & Colella 1984.
    A Mach-10 shock at 30° incidence; the standard benchmark for
    high-resolution shock-capturing schemes.

2.  2-D Riemann problems — Kurganov & Tadmor 2002 (19 configurations).
    Each configuration divides the unit square into four constant-state
    quadrants.  The interaction produces complex wave patterns including
    shocks, rarefactions, contact discontinuities, and vortex sheets.
"""

import numpy as np


# ---------------------------------------------------------------------------
# Double Mach Reflection
# ---------------------------------------------------------------------------

def ic_double_mach_reflection(nx, ny, gamma=1.4):
    """
    Initial condition for the Double Mach Reflection test problem.

    The domain is [0, 4] × [0, 1].  A Mach-10 planar shock is inclined at
    30° to the x-axis and initially positioned at x = 1/6.

    Parameters
    ----------
    nx, ny : int — *physical* cells (before ghost-cell padding).
    gamma  : float — ratio of specific heats (default 1.4).

    Returns
    -------
    r0, u0, v0, p0 : ndarray, shape (ny, nx) — primitive variables.
    preshock       : ndarray, shape (4,) — conserved state ahead of shock.
    postshock      : ndarray, shape (4,) — conserved state behind shock.
    shock_speed    : float — shock propagation speed in x.
    """
    print("Initial condition: Double Mach Reflection (Mach 10, 30° wedge)")

    # State ahead of the shock (undisturbed gas)
    r1, u1, v1, p1 = 1.4, 0.0, 0.0, 1.0
    E1 = p1 / (gamma - 1.0) + 0.5*r1*(u1**2 + v1**2)
    preshock = np.array([r1, r1*u1, r1*v1, E1])

    # State behind the shock via Rankine–Hugoniot relations (Ms = 10, θ = 30°)
    Ms  = 10.0
    c1  = np.sqrt(gamma*p1/r1)
    Tau = (gamma + 1.0) / (gamma - 1.0)

    p2  = p1 * (2.0*gamma*Ms**2 - (gamma - 1.0)) / (gamma + 1.0)
    r2  = r1 * (Tau*(p2/p1) + 1.0) / (Tau + (p2/p1))
    vn2 = Ms * (1.0 - ((gamma-1.0)*Ms**2 + 2.0) / ((gamma+1.0)*Ms**2)) * c1
    u2  =  vn2 * np.cos(np.pi/6.0)
    v2  = -vn2 * np.sin(np.pi/6.0)
    E2  = p2 / (gamma - 1.0) + 0.5*r2*(u2**2 + v2**2)
    postshock   = np.array([r2, r2*u2, r2*v2, E2])
    shock_speed = Ms / np.cos(np.pi/6.0)

    # Build physical grid
    dx = 4.0/nx;  dy = 1.0/ny
    xc = np.linspace(dx/2.0, 4.0 - dx/2.0, nx)
    yc = np.linspace(dy/2.0, 1.0 - dy/2.0, ny)
    x, y = np.meshgrid(xc, yc)

    # Initial field: undisturbed everywhere except where shock has passed
    r0 = r1 * np.ones((ny, nx))
    u0 = u1 * np.ones((ny, nx))
    v0 = v1 * np.ones((ny, nx))
    p0 = p1 * np.ones((ny, nx))

    # Region behind the initial shock front x < x0 + y*tan(30°)
    x0   = 1.0/6.0
    mask = x < x0 + y*np.tan(np.pi/6.0)
    r0[mask] = r2;  u0[mask] = u2;  v0[mask] = v2;  p0[mask] = p2

    return r0, u0, v0, p0, preshock, postshock, shock_speed


# ---------------------------------------------------------------------------
# 2-D Riemann problems
# ---------------------------------------------------------------------------

# Each configuration is (p, r, u, v) for regions [reg1, reg2, reg3, reg4]
# following the quadrant layout:
#
#   1.0 +-----------+-----------+
#       |  reg 2    |  reg 1    |
#   0.5 +-----------+-----------+
#       |  reg 3    |  reg 4    |
#   0.0 +-----------+-----------+
#       0.0        0.5         1.0
#
_RIEMANN_CONFIGS = {
    1:  {'p': [1.0,   0.4,    0.0439, 0.15  ],
         'r': [1.0,   0.5197, 0.1072, 0.2579],
         'u': [0.0,  -0.7259,-0.7259, 0.0   ],
         'v': [0.0,   0.0,   -1.4045,-1.4045]},
    2:  {'p': [1.0,   0.4,    1.0,    0.4   ],
         'r': [1.0,   0.5197, 1.0,    0.5197],
         'u': [0.0,  -0.7259,-0.7259, 0.0   ],
         'v': [0.0,   0.0,   -0.7259,-0.7259]},
    3:  {'p': [1.5,   0.3,    0.029,  0.3   ],
         'r': [1.5,   0.5323, 0.138,  0.5323],
         'u': [0.0,   1.206,  1.206,  0.0   ],
         'v': [0.0,   0.0,    1.206,  1.206 ]},
    4:  {'p': [1.1,   0.35,   1.1,    0.35  ],
         'r': [1.1,   0.5065, 1.1,    0.5065],
         'u': [0.0,   0.8939, 0.8939, 0.0   ],
         'v': [0.0,   0.0,    0.8939, 0.8939]},
    5:  {'p': [1.0,   1.0,    1.0,    1.0   ],
         'r': [1.0,   2.0,    1.0,    3.0   ],
         'u': [-0.75,-0.75,   0.75,   0.75  ],
         'v': [-0.5,  0.5,    0.5,   -0.5   ]},
    6:  {'p': [1.0,   1.0,    1.0,    1.0   ],
         'r': [1.0,   2.0,    1.0,    3.0   ],
         'u': [0.75,  0.75,  -0.75,  -0.75  ],
         'v': [-0.5,  0.5,    0.5,   -0.5   ]},
    7:  {'p': [1.0,   0.4,    0.4,    0.4   ],
         'r': [1.0,   0.5197, 0.8,    0.5197],
         'u': [0.1,  -0.6259, 0.1,    0.1   ],
         'v': [0.1,   0.1,    0.1,   -0.6259]},
    8:  {'p': [0.4,   1.0,    1.0,    1.0   ],
         'r': [0.5197, 1.0,   0.8,    1.0   ],
         'u': [0.1,  -0.6259, 0.1,    0.1   ],
         'v': [0.1,   0.1,    0.1,   -0.6259]},
    9:  {'p': [1.0,   1.0,    0.4,    0.4   ],
         'r': [1.0,   2.0,    1.039,  0.5197],
         'u': [0.0,   0.0,    0.0,    0.0   ],
         'v': [0.3,  -0.3,   -0.8133,-0.4259]},
    10: {'p': [1.0,   1.0,    0.3333, 0.3333],
         'r': [1.0,   0.5,    0.2281, 0.4562],
         'u': [0.0,   0.0,    0.0,    0.0   ],
         'v': [0.4297, 0.6076,-0.6076,-0.4297]},
    11: {'p': [1.0,   0.4,    0.4,    0.4   ],
         'r': [1.0,   0.5313, 0.8,    0.5313],
         'u': [0.1,   0.8276, 0.1,    0.1   ],
         'v': [0.0,   0.0,    0.0,    0.7276]},
    12: {'p': [0.4,   1.0,    1.0,    1.0   ],
         'r': [0.5313, 1.0,   0.8,    1.0   ],
         'u': [0.0,   0.7276, 0.0,    0.0   ],
         'v': [0.0,   0.0,    0.0,    0.7276]},
    13: {'p': [1.0,   1.0,    0.4,    0.4   ],
         'r': [1.0,   2.0,    1.0625, 0.5313],
         'u': [0.0,   0.0,    0.0,    0.0   ],
         'v': [-0.3,  0.3,    0.8145, 0.4276]},
    14: {'p': [8.0,   8.0,    2.6667, 2.6667],
         'r': [2.0,   1.0,    0.4736, 0.9474],
         'u': [0.0,   0.0,    0.0,    0.0   ],
         'v': [-0.5606,-1.2172, 1.2172, 1.1606]},
    15: {'p': [1.0,   0.4,    0.4,    0.4   ],
         'r': [1.0,   0.5197, 0.8,    0.5313],
         'u': [0.1,  -0.6259, 0.1,    0.1   ],
         'v': [-0.3, -0.3,   -0.3,    0.4276]},
    16: {'p': [0.4,   1.0,    1.0,    1.0   ],
         'r': [0.5313, 1.0222, 0.8,   1.0   ],
         'u': [0.1,  -0.6179, 0.1,    0.1   ],
         'v': [0.1,   0.1,    0.1,    0.8276]},
    17: {'p': [1.0,   1.0,    0.4,    0.4   ],
         'r': [1.0,   2.0,    1.0625, 0.5197],
         'u': [0.0,   0.0,    0.0,    0.0   ],
         'v': [-0.4, -0.3,    0.2145,-1.1259]},
    18: {'p': [1.0,   1.0,    0.4,    0.4   ],
         'r': [1.0,   2.0,    1.0625, 0.5197],
         'u': [0.0,   0.0,    0.0,    0.0   ],
         'v': [1.0,  -0.3,    0.2145, 0.2741]},
    19: {'p': [1.0,   1.0,    0.4,    0.4   ],
         'r': [1.0,   2.0,    1.0625, 0.5197],
         'u': [0.0,   0.0,    0.0,    0.0   ],
         'v': [0.3,  -0.3,    0.2145,-0.4259]},
    # Special 1-D Riemann tests mapped to 2-D
    'Sod_x': {'p': [0.1,   1.0,    1.0,    0.1  ],
              'r': [0.125, 1.0,    1.0,    0.125],
              'u': [0.0,   0.0,    0.0,    0.0  ],
              'v': [0.0,   0.0,    0.0,    0.0  ]},
    'Sod_y': {'p': [1.0,   1.0,    0.1,    0.1  ],
              'r': [1.0,   1.0,    0.125,  0.125],
              'u': [0.0,   0.0,    0.0,    0.0  ],
              'v': [0.0,   0.0,    0.0,    0.0  ]},
}


def ic_riemann2d(x, y, config_id, gamma=1.4):
    """
    Initial condition for a 2-D Riemann problem (Kurganov & Tadmor 2002).

    Parameters
    ----------
    x, y      : ndarray, shape (ny, nx) — cell-centre coordinates.
    config_id : int (1–19) or str ('Sod_x', 'Sod_y').
    gamma     : float — ratio of specific heats.

    Returns
    -------
    r0, u0, v0, p0 : ndarray, shape (ny, nx).
    """
    if config_id not in _RIEMANN_CONFIGS:
        raise ValueError(f"Unknown Riemann configuration '{config_id}'. "
                         f"Available: 1–19, 'Sod_x', 'Sod_y'.")

    cfg = _RIEMANN_CONFIGS[config_id]
    p_v = cfg['p'];  r_v = cfg['r'];  u_v = cfg['u'];  v_v = cfg['v']

    print(f"Initial condition: 2-D Riemann Configuration {config_id}")
    print(f"{'':8s} {'reg1':>8s} {'reg2':>8s} {'reg3':>8s} {'reg4':>8s}")
    for lbl, vals in zip(('density', 'x-vel', 'y-vel', 'pressure'),
                         (r_v, u_v, v_v, p_v)):
        print(f"  {lbl:8s}: {vals[0]:8.4f} {vals[1]:8.4f} "
              f"{vals[2]:8.4f} {vals[3]:8.4f}")

    reg1 = (x >= 0.5) & (y >= 0.5)
    reg2 = (x <  0.5) & (y >= 0.5)
    reg3 = (x <  0.5) & (y <  0.5)
    reg4 = (x >= 0.5) & (y <  0.5)

    r0 = r_v[0]*reg1 + r_v[1]*reg2 + r_v[2]*reg3 + r_v[3]*reg4
    u0 = u_v[0]*reg1 + u_v[1]*reg2 + u_v[2]*reg3 + u_v[3]*reg4
    v0 = v_v[0]*reg1 + v_v[1]*reg2 + v_v[2]*reg3 + v_v[3]*reg4
    p0 = p_v[0]*reg1 + p_v[1]*reg2 + p_v[2]*reg3 + p_v[3]*reg4

    return r0, u0, v0, p0
