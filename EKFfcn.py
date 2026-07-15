"""
EKF Vehicle State Estimator — Python implementation for Simulink "Python Code" block.

Converted from:
  EKFfcn.m               — main EKF filter
  DiscreteStateEstimator.m — state equations & Jacobian (LittleF / BigF)

State vector  X = [Vy, Vx, r, Psi, Ax, Rbias, Abias]   (7 x 1)
Measurement   z = [VYmeas, VXmeas, rmeas, Axmeas]        (4 x 1)

MODERNIZATION:
- Uses SymPy to symbolically compute process Jacobian (BigF) and measurement Jacobian (Hk).
- Eliminates manual derivative code (`_dmagic_dalpha`) entirely.
- Uses `sp.lambdify` with the NumPy backend for zero-overhead execution during timesteps.
"""

import numpy as np
import sympy as sp

class EKFfcn:
    # ---------------------------------------------------------------------------
    # Vehicle / Pacejka tyre constants  (match MATLAB values exactly)
    # ---------------------------------------------------------------------------
    _a  = 0.3036
    _b  = 0.2484
    _B  = 10.0
    _C  = 1.6
    _D  = 73.1580
    _E  = 1.0
    _Iz = 1.0
    _Mv = 15.7
    _Ts = 0.001


    # ---------------------------------------------------------------------------
    # SymPy Symbolic Generation (Executed ONCE at module import/initialization)
    # ---------------------------------------------------------------------------
    # Define symbols for states, inputs, and constants
    Vy_s, Vx_s, r_s, Psi_s, Ax_s, Rbias_s, Abias_s, delta_s = sp.symbols(
        'Vy_s Vx_s r_s Psi_s Ax_s Rbias_s Abias_s delta_s'
    )
    a_s, b_s, B_s, C_s, D_s, E_s, Iz_s, Mv_s, Ts_s = sp.symbols(
        'a_s b_s B_s C_s D_s E_s Iz_s Mv_s Ts_s'
    )

    # Symbolic Pacejka Magic Formula helper
    def _magic_sym(alpha, B, C, D, E):
        Ba = B * alpha
        return D * sp.sin(C * sp.atan(Ba - E * (Ba - sp.atan(Ba))))

    # Slip angles (including the Rbias fix)
    AlphaF_s = delta_s - sp.atan((-Vy_s + a_s * (-r_s + Rbias_s)) / Vx_s)
    AlphaR_s = - sp.atan((-Vy_s - b_s * (-r_s + Rbias_s)) / Vx_s)

    Ff_s = _magic_sym(AlphaF_s, B_s, C_s, D_s, E_s)
    Fr_s = _magic_sym(AlphaR_s, B_s, C_s, D_s, E_s)

    # Non-linear state transition equations f(x)
    Vy_next  = Vy_s  + ((Fr_s + Ff_s * sp.cos(delta_s)) / Mv_s - Vx_s * r_s) * Ts_s
    Vx_next  = Vx_s  + ((Ax_s + Abias_s) - Vy_s * r_s) * Ts_s
    r_next   = r_s   + ((Ff_s * sp.cos(delta_s) * a_s - Fr_s * b_s) * Ts_s) / Iz_s
    Psi_next = Psi_s + r_s * Ts_s
    Ax_next  = Ax_s  + Abias_s
    Rbias_next = Rbias_s
    Abias_next = Abias_s

    f_symbolic = sp.Matrix([Vy_next, Vx_next, r_next, Psi_next, Ax_next, Rbias_next, Abias_next])
    X_symbolic = sp.Matrix([Vy_s, Vx_s, r_s, Psi_s, Ax_s, Rbias_s, Abias_s])

    # Compute State Jacobian (Big F) symbolically
    Fk_symbolic = f_symbolic.jacobian(X_symbolic)

    # Lambdify State Jacobian into high-performance NumPy code
    _jacobian_Fk_lam = staticmethod(sp.lambdify(
        (Vy_s, Vx_s, r_s, Psi_s, Ax_s, Rbias_s, Abias_s, delta_s, a_s, b_s, B_s, C_s, D_s, E_s, Iz_s, Mv_s, Ts_s),
        Fk_symbolic,
        modules='numpy'
    ))

    # Measurement equations h(x)
    h_symbolic = sp.Matrix([
        Vx_s * sp.sin(Psi_s) + Vy_s * sp.cos(Psi_s),
        Vx_s * sp.cos(Psi_s) - Vy_s * sp.sin(Psi_s),
        r_s,
        Ax_s
    ])

    # Compute Measurement Jacobian (Hk) symbolically
    Hk_symbolic = h_symbolic.jacobian(X_symbolic)

    # Lambdify Measurement Jacobian into high-performance NumPy code
    _jacobian_Hk_lam = staticmethod(sp.lambdify(
        (Vy_s, Vx_s, r_s, Psi_s, Ax_s, Rbias_s, Abias_s),
        Hk_symbolic,
        modules='numpy'
    ))


    # ---------------------------------------------------------------------------
    # Helper: Pacejka magic formula (Numeric version for _little_f)
    # ---------------------------------------------------------------------------
    def _magic(self, alpha, B, C, D, E):
        """D * sin(C * atan(B*alpha - E*(B*alpha - atan(B*alpha))))"""
        Ba  = B * alpha
        return D * np.sin(C * np.arctan(Ba - E * (Ba - np.arctan(Ba))))


    # ---------------------------------------------------------------------------
    # LittleF — nonlinear state-transition  f(x_k)
    # ---------------------------------------------------------------------------
    def _little_f(self, Vypri, Vxpri, rpri, Psipri, Axpri, Rbiaspri, Abiaspri, delta,
                a=_a, b=_b, B=_B, C=_C, D=_D, E=_E, Iz=_Iz, Mv=_Mv, Ts=_Ts):
        """
        Returns fk = [Vy, Vx, r, Psi, Ax, Rbias, Abias]  (7-element 1-D array)
        """
        AlphaF = delta - np.arctan((-Vypri + a * (-rpri + Rbiaspri)) / Vxpri)
        AlphaR =       - np.arctan((-Vypri - b * (-rpri + Rbiaspri)) / Vxpri)

        Ff = self._magic(AlphaF, B, C, D, E)
        Fr = self._magic(AlphaR, B, C, D, E)

        Vy    = Vypri  + ((Fr + Ff * np.cos(delta)) / Mv - Vxpri * rpri) * Ts
        Vx    = Vxpri  + ((Axpri + Abiaspri) - Vypri * rpri) * Ts
        r     = rpri   + ((Ff * np.cos(delta) * a - Fr * b) * Ts) / Iz   
        Psi   = Psipri + rpri * Ts                                         
        Ax    = Axpri  + Abiaspri
        Rbias = Rbiaspri
        Abias = Abiaspri

        return np.array([Vy, Vx, r, Psi, Ax, Rbias, Abias], dtype=float)


    # ---------------------------------------------------------------------------
    # BigF — Jacobian of f w.r.t. state  (7 x 7) via SymPy Compiled Function
    # ---------------------------------------------------------------------------
    def _big_f(self, Vypri, Vxpri, rpri, Psipri, Axpri, Rbiaspri, Abiaspri, delta,
            a=_a, b=_b, B=_B, C=_C, D=_D, E=_E, Iz=_Iz, Mv=_Mv, Ts=_Ts):
        """
        Returns Fk  (7 x 7 numpy array) — Jacobian of _little_f w.r.t. X_prior.
        """
        Fk = self._jacobian_Fk_lam(
            Vypri, Vxpri, rpri, Psipri, Axpri, Rbiaspri, Abiaspri, delta,
            a, b, B, C, D, E, Iz, Mv, Ts
        )
        return Fk


    # ---------------------------------------------------------------------------
    # Main EKF step — mirrors EKFfcn.m exactly
    # ---------------------------------------------------------------------------
    def ekf_step(self, VYmeas, VXmeas, rmeas, Axmeas,
                Xcov, Ycov,
                Ppri_flat, Xpri_arr,
                startflag, SA, Ts):
        """
        Extended Kalman Filter step.
        """
        delta = float(SA)

        rvar  = 0.01
        Axvar = 0.001

        Rk = np.array([
            [Ycov, 0.0,  0.0,   0.0  ],
            [0.0,  Xcov, 0.0,   0.0  ],
            [0.0,  0.0,  rvar,  0.0  ],
            [0.0,  0.0,  0.0,   Axvar],
        ], dtype=float)

        Qk = np.diag([0.02, 0.02, 0.01, 0.01, 0.0001, 0.0001, 0.0001])

        zk = np.array([VYmeas, VXmeas, rmeas, Axmeas], dtype=float)

        Xpri = np.asarray(Xpri_arr, dtype=float).ravel()
        Ppri = np.asarray(Ppri_flat, dtype=float).reshape(7, 7)

        if startflag:
            Ppredict = np.diag([0.1, 0.1, 0.1, 0.1, 0.001, 0.001, 0.001])
            Xpri     = np.array([0.5, -1.9, -0.49, 0.0, 0.34, -0.022, -0.2], dtype=float)
        else:
            Ppredict = Ppri.copy()

        Vypri    = Xpri[0]
        Vxpri    = Xpri[1]
        rpri     = Xpri[2]
        Psipri   = Xpri[3]
        Axpri    = Xpri[4]
        Rbiaspri = Xpri[5]
        Abiaspri = Xpri[6]

        # ---- predict ----
        fk = self._little_f(Vypri, Vxpri, rpri, Psipri, Axpri, Rbiaspri, Abiaspri, delta, Ts=Ts)
        Fk = self._big_f(   Vypri, Vxpri, rpri, Psipri, Axpri, Rbiaspri, Abiaspri, delta, Ts=Ts)

        Ppost = Fk @ Ppredict @ Fk.T + Qk

        # ---- measurement model  h(fk) ----
        Vy_p  = fk[0]
        Vx_p  = fk[1]
        r_p   = fk[2]
        Psi_p = fk[3]
        Ax_p  = fk[4]

        hk = np.array([
            Vx_p * np.sin(Psi_p) + Vy_p * np.cos(Psi_p),
            Vx_p * np.cos(Psi_p) - Vy_p * np.sin(Psi_p),
            r_p,
            Ax_p,
        ], dtype=float)

        # ---- Jacobian of h via SymPy Compiled Function ----
        Hk = self._jacobian_Hk_lam(*fk)

        # ---- update ----
        res      = zk - hk
        Sk       = Hk @ Ppost @ Hk.T + Rk
        Kk       = Ppost @ Hk.T @ np.linalg.inv(Sk)
        Xpri_out = fk + Kk @ res
        Ppri_out = (np.eye(7) - Kk @ Hk) @ Ppost

        return (
            Ppri_out.ravel(),       
            Xpri_out,
            res,
            fk,
            float(Xpri_out[5]),
            float(Xpri_out[6]),
        )


    # ---------------------------------------------------------------------------
    # Simulink "Python Code" block entry point
    # ---------------------------------------------------------------------------
    def simulink_step(self, u):
        u = list(u)

        VYmeas    = u[0]
        VXmeas    = u[1]
        rmeas     = u[2]
        Axmeas    = u[3]
        Xcov      = u[4]
        Ycov      = u[5]

        Ppri_col  = np.array(u[6:55], dtype=float)
        Ppri_flat = Ppri_col.reshape(7, 7, order='F').ravel(order='C')  

        Xpri_arr  = u[55:62]
        startflag = u[62]
        SA        = u[63]

        Ppri_out_rm, Xpri_out, res, fk, Rb, Ab = self.ekf_step(
            VYmeas, VXmeas, rmeas, Axmeas,
            Xcov, Ycov,
            Ppri_flat, Xpri_arr,
            startflag, SA,
        )

        Ppri_out_cm = Ppri_out_rm.reshape(7, 7, order='C').ravel(order='F')

        y = (list(Ppri_out_cm)   
        + list(Xpri_out)      
        + list(res)           
        + list(fk)            
        + [Rb, Ab])           
        return y