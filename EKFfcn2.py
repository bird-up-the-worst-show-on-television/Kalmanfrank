import numpy as np

class EKFfcn2:
    """
    EKF Vehicle State Estimator — Python implementation for Simulink "Python Code" block.

    Converted from:
    EKFfcn.m               — main EKF filter
    DiscreteStateEstimator.m — state equations & Jacobian (LittleF / BigF)

    State vector  X = [Vy, Vx, r, Psi, Ax, Rbias, Abias]   (7 x 1)
    Measurement   z = [VYmeas, VXmeas, rmeas, Axmeas]        (4 x 1)
                      [Inertial Frame velocity, Inertial frame velocity, angular velocity, longitudinal acceleration]
    FIXES vs original translation
    ------------------------------
    1. Ppri is converted from Simulink column-major to row-major on input,
    and back to column-major on output, preventing state cross-coupling.
    2. r and Psi state equations use rpri only (not rpri+Rbiaspri), preventing
    Rbias from accumulating into the integrated states each timestep.
    3. Corresponding Jacobian rows 2 and 3 corrected (no spurious +1 on Rbias col,
    dPsi/dRbias set to 0).
    """

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

    def __init__(self):
        return

    # ---------------------------------------------------------------------------
    # Helper: Pacejka magic formula
    # ---------------------------------------------------------------------------
    def _magic(self, alpha, B, C, D, E):
        """D * sin(C * atan(B*alpha - E*(B*alpha - atan(B*alpha))))"""
        Ba  = B * alpha
        return D * np.sin(C * np.arctan(Ba - E * (Ba - np.arctan(Ba))))


    def _dmagic_dalpha(self, alpha, B, C, D, E):
        """Derivative of magic formula w.r.t. alpha."""
        Ba      = B * alpha
        atn     = np.arctan(Ba)
        inner   = Ba - E * (Ba - atn)
        d_inner = 1.0 - E * (1.0 - 1.0 / (1.0 + Ba**2))
        return D * C * np.cos(C * np.arctan(inner)) / (1.0 + inner**2) * d_inner * B


    # ---------------------------------------------------------------------------
    # LittleF — nonlinear state-transition  f(x_k)
    #
    # FIX: r and Psi integrate from rpri only.
    #      Rbiaspri enters only through the tyre slip angles (AlphaF, AlphaR),
    #      which is where the original MATLAB symbolic model uses it.
    # ---------------------------------------------------------------------------
    def _little_f(self, Vypri, Vxpri, rpri, Psipri, Axpri, Rbiaspri, Abiaspri, delta,
                a=_a, b=_b, B=_B, C=_C, D=_D, E=_E, Iz=_Iz, Mv=_Mv, Ts=_Ts):
        """
        Returns fk = [Vy, Vx, r, Psi, Ax, Rbias, Abias]  (7-element 1-D array)
        """
        AlphaF = delta - np.arctan((-Vypri + a * (rpri + Rbiaspri)) / Vxpri)
        AlphaR =       - np.arctan((-Vypri - b * (rpri + Rbiaspri)) / Vxpri)

        Ff = self._magic(AlphaF, B, C, D, E)
        Fr = self._magic(AlphaR, B, C, D, E)

        Vy    = Vypri  + ((Fr + Ff * np.cos(delta)) / Mv - Vxpri * rpri) * Ts
        Vx    = Vxpri  + ((Axpri + Abiaspri) - Vypri * rpri) * Ts
        r     = rpri   + ((Ff * np.cos(delta) * a - Fr * b) * Ts) / Iz   # FIX: rpri only
        Psi   = Psipri + rpri * Ts                                         # FIX: rpri only
        Ax    = Axpri  + Abiaspri
        Rbias = Rbiaspri
        Abias = Abiaspri

        return np.array([Vy, Vx, r, Psi, Ax, Rbias, Abias], dtype=float)


    # ---------------------------------------------------------------------------
    # BigF — Jacobian of f w.r.t. state  (7 x 7)
    #
    # FIX: Row 2 col 5: no spurious +1 (Rbias does not directly integrate into r).
    #      Row 3 col 5: dPsi/dRbias = 0 (Rbias no longer in Psi equation).
    # ---------------------------------------------------------------------------
    def _big_f(self, Vypri, Vxpri, rpri, Psipri, Axpri, Rbiaspri, Abiaspri, delta,
            a=_a, b=_b, B=_B, C=_C, D=_D, E=_E, Iz=_Iz, Mv=_Mv, Ts=_Ts):
        """
        Returns Fk  (7 x 7 numpy array) — Jacobian of _little_f w.r.t. X_prior.
        State order: [Vy(0), Vx(1), r(2), Psi(3), Ax(4), Rbias(5), Abias(6)]
        """
        r_eff = rpri + Rbiaspri

        num_F = -Vypri + a * r_eff
        num_R = -Vypri - b * r_eff

        denom_F = Vxpri**2 + num_F**2
        denom_R = Vxpri**2 + num_R**2

        # Partials of AlphaF
        dAlphaF_dVy = Vxpri  / denom_F
        dAlphaF_dVx = -num_F / denom_F
        dAlphaF_dr  = -a * Vxpri / denom_F
        dAlphaF_dRb = dAlphaF_dr          # same dependence on r_eff

        # Partials of AlphaR
        dAlphaR_dVy = -Vxpri / denom_R
        dAlphaR_dVx =  num_R / denom_R
        dAlphaR_dr  =  b * Vxpri / denom_R
        dAlphaR_dRb = dAlphaR_dr

        # dFf/dAlphaF,  dFr/dAlphaR
        dFf_dAF = self._dmagic_dalpha(delta - np.arctan(num_F / Vxpri), B, C, D, E)
        dFr_dAR = self._dmagic_dalpha(       -np.arctan(num_R / Vxpri), B, C, D, E)

        cosd = np.cos(delta)

        # Chain-rule partials of Ff and Fr w.r.t. each state
        dFf_dVy = dFf_dAF * dAlphaF_dVy
        dFf_dVx = dFf_dAF * dAlphaF_dVx
        dFf_dr  = dFf_dAF * dAlphaF_dr
        dFf_dRb = dFf_dAF * dAlphaF_dRb

        dFr_dVy = dFr_dAR * dAlphaR_dVy
        dFr_dVx = dFr_dAR * dAlphaR_dVx
        dFr_dr  = dFr_dAR * dAlphaR_dr
        dFr_dRb = dFr_dAR * dAlphaR_dRb

        Fk = np.zeros((7, 7), dtype=float)

        # Row 0: dVy/d(...)
        # Vy = Vypri + ((Fr + Ff*cos(d))/Mv - Vx*r) * Ts
        Fk[0, 0] = 1.0 + Ts * (dFr_dVy + dFf_dVy * cosd) / Mv
        Fk[0, 1] = Ts  * ((dFr_dVx + dFf_dVx * cosd) / Mv - rpri)
        Fk[0, 2] = Ts  * ((dFr_dr  + dFf_dr  * cosd) / Mv - Vxpri)
        Fk[0, 3] = 0.0
        Fk[0, 4] = 0.0
        Fk[0, 5] = Ts  * (dFr_dRb + dFf_dRb * cosd) / Mv
        Fk[0, 6] = 0.0

        # Row 1: dVx/d(...)
        # Vx = Vxpri + ((Ax + Abias) - Vy*r) * Ts
        Fk[1, 0] = -Ts * rpri
        Fk[1, 1] = 1.0
        Fk[1, 2] = -Ts * Vypri
        Fk[1, 3] = 0.0
        Fk[1, 4] = Ts
        Fk[1, 5] = 0.0
        Fk[1, 6] = Ts

        # Row 2: dr/d(...)
        # r = rpri + (Ff*cos(d)*a - Fr*b) * Ts / Iz
        # FIX: no +1 on Rbias column — Rbias enters only via dFf/dRb and dFr/dRb
        Fk[2, 0] = Ts * (dFf_dVy * cosd * a - dFr_dVy * b) / Iz
        Fk[2, 1] = Ts * (dFf_dVx * cosd * a - dFr_dVx * b) / Iz
        Fk[2, 2] = 1.0 + Ts * (dFf_dr  * cosd * a - dFr_dr  * b) / Iz
        Fk[2, 3] = 0.0
        Fk[2, 4] = 0.0
        Fk[2, 5] = Ts * (dFf_dRb * cosd * a - dFr_dRb * b) / Iz   # FIX: no +1
        Fk[2, 6] = 0.0

        # Row 3: dPsi/d(...)
        # Psi = Psipri + rpri * Ts  (Rbias no longer here)
        Fk[3, 0] = 0.0
        Fk[3, 1] = 0.0
        Fk[3, 2] = Ts
        Fk[3, 3] = 1.0
        Fk[3, 4] = 0.0
        Fk[3, 5] = 0.0    # FIX: was Ts, now 0 since Rbias removed from Psi
        Fk[3, 6] = 0.0

        # Row 4: dAx/d(...)
        Fk[4, 4] = 1.0
        Fk[4, 6] = 1.0

        # Row 5: dRbias/d(...)
        Fk[5, 5] = 1.0

        # Row 6: dAbias/d(...)
        Fk[6, 6] = 1.0
        return Fk


    # ---------------------------------------------------------------------------
    # Main EKF step — mirrors EKFfcn.m exactly
    # ---------------------------------------------------------------------------
    def ekf_step(self, VYmeas, VXmeas, rmeas, Axmeas,
                Xcov, Ycov,
                Ppri_flat, Xpri_arr,
                startflag, SA):
        """
        Extended Kalman Filter step.

        Parameters
        ----------
        VYmeas, VXmeas, rmeas, Axmeas : float
        Xcov, Ycov                    : float   measurement noise variances
        Ppri_flat : array-like, 49 floats       prior covariance (7x7, ROW-major)
        Xpri_arr  : array-like,  7 floats       prior state
        startflag : int  — 0 on first call (triggers initialisation)
        SA        : float — steering angle (radians)

        Returns
        -------
        Ppri_out : np.ndarray (49,)   updated covariance, ROW-major flattened
        Xpri_out : np.ndarray (7,)    updated state
        res      : np.ndarray (4,)    innovation residual
        fk       : np.ndarray (7,)    predicted state (before update)
        Rbiaspri : float
        Abiaspri : float
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

        if int(startflag) == 0:
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
        fk = self._little_f(Vypri, Vxpri, rpri, Psipri, Axpri, Rbiaspri, Abiaspri, delta)
        Fk = self._big_f(   Vypri, Vxpri, rpri, Psipri, Axpri, Rbiaspri, Abiaspri, delta)

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

        # ---- Jacobian of h ----
        # hk[0] = Vx*sin(Psi) + Vy*cos(Psi)  →  dh0/dVy=cos, dh0/dVx=sin, dh0/dPsi=Vx*cos-Vy*sin
        # hk[1] = Vx*cos(Psi) - Vy*sin(Psi)  →  dh1/dVy=-sin, dh1/dVx=cos, dh1/dPsi=-Vx*sin-Vy*cos
        Hk = np.array([
            [ np.cos(Psi_p),  np.sin(Psi_p), 0., Vx_p*np.cos(Psi_p) - Vy_p*np.sin(Psi_p), 0., 0., 0.],
            [-np.sin(Psi_p),  np.cos(Psi_p), 0., -Vx_p*np.sin(Psi_p) - Vy_p*np.cos(Psi_p), 0., 0., 0.],
            [ 0.,             0.,            1.,  0.,                                          0., 0., 0.],
            [ 0.,             0.,            0.,  0.,                                          1., 0., 0.],
        ], dtype=float)

        # ---- update ----
        res      = zk - hk
        Sk       = Hk @ Ppost @ Hk.T + Rk
        Kk       = Ppost @ Hk.T @ np.linalg.inv(Sk)
        Xpri_out = fk + Kk @ res
        Ppri_out = (np.eye(7) - Kk @ Hk) @ Ppost

        return (
            Ppri_out.ravel(),       # row-major; simulink_step converts to col-major
            Xpri_out,
            res,
            fk,
            float(Xpri_out[5]),
            float(Xpri_out[6]),
        )


    # ---------------------------------------------------------------------------
    # Simulink "Python Code" block entry point
    #
    # FIX: Ppri arrives from Simulink Unit Delay in column-major order.
    #      Convert to row-major before use, and back to column-major on output.
    #
    # Input vector u layout (64 elements total):
    #   u[0]      VYmeas
    #   u[1]      VXmeas
    #   u[2]      rmeas
    #   u[3]      Axmeas
    #   u[4]      Xcov
    #   u[5]      Ycov
    #   u[6:55]   Ppri  (7x7, COLUMN-major from Simulink)
    #   u[55:62]  Xpri  (7 states)
    #   u[62]     startflag
    #   u[63]     SA
    #
    # Output vector y layout (69 elements total):
    #   y[0:49]   Ppri_out  (7x7, COLUMN-major for Simulink Unit Delay)
    #   y[49:56]  Xpri_out  (7 states)
    #   y[56:60]  res       (4 residuals)
    #   y[60:67]  fk        (7 predicted states)
    #   y[67]     Rbiaspri
    #   y[68]     Abiaspri
    # ---------------------------------------------------------------------------
    def simulink_step(self, u):
        u = list(u)

        VYmeas    = u[0]
        VXmeas    = u[1]
        rmeas     = u[2]
        Axmeas    = u[3]
        Xcov      = u[4]
        Ycov      = u[5]

        # FIX: Simulink stores matrices column-major; reshape accordingly then
        #      convert to row-major for internal numpy operations
        Ppri_col  = np.array(u[6:55], dtype=float)
        Ppri_flat = Ppri_col.reshape(7, 7, order='F').ravel(order='C')  # col→row major

        Xpri_arr  = u[55:62]
        startflag = u[62]
        SA        = u[63]

        Ppri_out_rm, Xpri_out, res, fk, Rb, Ab = self.ekf_step(
            VYmeas, VXmeas, rmeas, Axmeas,
            Xcov, Ycov,
            Ppri_flat, Xpri_arr,
            startflag, SA,
        )

        # FIX: convert Ppri back to column-major for Simulink Unit Delay
        Ppri_out_cm = Ppri_out_rm.reshape(7, 7, order='C').ravel(order='F')

        y = (list(Ppri_out_cm)   # 49  col-major
        + list(Xpri_out)      #  7
        + list(res)           #  4
        + list(fk)            #  7
        + [Rb, Ab])           #  2
        return y
