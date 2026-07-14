import utm
import pandas as pd
import numpy as np
from EKFfcn import EKFfcn
import matplotlib.pyplot as plt
import math

GPS_FILE = "/home/lpschexn/Kalmanfrank/res/Lot24CRCGPS.csv"
INS_FILE = "/home/lpschexn/Kalmanfrank/res/Lot24CRCINS.csv"

NUM_STATES = 7
NUM_MEASUREMENTS = 4


# Moving-mean low pass filter window size (samples).
# Because INS is sampled 4x faster than GPS, the INS window is scaled
# so both filters cover the same amount of wall-clock time.
GPS_WINDOW = 5
INS_WINDOW = 100

gps_data = pd.read_csv(GPS_FILE)
ins_data = pd.read_csv(INS_FILE)

# Pre-compute low-pass-filtered measurement columns.
# min_periods=1 avoids NaNs at the start of the run (partial windows are used
# until the window fills). Rolling mean is causal by default: at index i it
# averages samples [i-window+1 ... i], so no future data leaks in.
# gps_vy_lpf = gps_data["field.twist.twist.linear.y"].rolling(window=GPS_WINDOW, min_periods=1).mean()
# gps_vx_lpf = gps_data["field.twist.twist.linear.x"].rolling(window=GPS_WINDOW, min_periods=1).mean()
ins_wz_lpf = ins_data["field.angular_velocity.z"].rolling(window=INS_WINDOW, min_periods=1).mean()
ins_ax_lpf = ins_data["field.linear_acceleration.x"].rolling(window=INS_WINDOW, min_periods=1).mean()

ins_gps_pub_ratio = math.floor(ins_wz_lpf.shape[0] / gps_data.shape[0])

start_flag = True

ekf = EKFfcn()

P_pri = np.diag(np.ones(NUM_STATES))
P_pri = P_pri.ravel(order="C")
X_pri = [0.5, 1.9, -0.49, 0.0, 0.34, -0.022, -0.2]
residual = np.zeros(NUM_MEASUREMENTS)
predicted_state = np.zeros(NUM_STATES)
R_bias_prior = 0.0
A_bias_prior = 0.0

states = pd.DataFrame()
residuals = pd.DataFrame()

steering_angle = 0.135

for i in range(0, gps_data.shape[0]):
    gps = gps_data.iloc[i]
    ins_idx = i * ins_gps_pub_ratio
    # ins = ins_data.iloc[ins_idx]

    P_pri, X_pri, residual, predicted_state, R_bias_prior, A_bias_prior = ekf.ekf_step(
        gps["field.twist.twist.linear.y"],gps["field.twist.twist.linear.x"],
        ins_wz_lpf.iloc[ins_idx], ins_ax_lpf.iloc[ins_idx],
        gps["field.twist.covariance0"], gps["field.twist.covariance0.1"],
        P_pri, X_pri,
        start_flag, steering_angle)

    if start_flag:
        start_flag = False

    new_states = pd.DataFrame({"Vy":    X_pri[0], "Vx":    X_pri[1],
                               "r":     X_pri[2], "Psi":   X_pri[3],
                               "Ax":    X_pri[4], "Rbias": X_pri[5],
                               "Abias": X_pri[6]}, index=[i])
    residual_df = pd.DataFrame({"Vy": residual[0], "Vx": residual[1],
                                "r":  residual[2], "Ax": residual[3]}, index=[i])
    states = pd.concat([states, new_states])
    residuals = pd.concat([residuals, residual_df])

fig, axes = plt.subplots(nrows=2, ncols=4)
plt.suptitle("State Estimates")
axes = axes.ravel()
for i in range(0, NUM_STATES):
    axes[i].plot(states.iloc[:, i])
    axes[i].set_xlabel(states.columns[i])

fig, axes_res = plt.subplots(nrows=1, ncols=4)
plt.suptitle("Residuals")
axes_res = axes_res.ravel()
for i in range(0, NUM_MEASUREMENTS):
    axes_res[i].plot(residuals.iloc[:, i])
    axes_res[i].set_xlabel(residuals.columns[i])

fig, axes_in = plt.subplots(nrows=1, ncols=4)
plt.suptitle("Inputs (raw vs. LPF)")
axes_in[0].plot(gps_data["field.twist.twist.linear.y"], alpha=0.4, label="raw")
axes_in[0].plot(gps["field.twist.twist.linear.y"], label="LPF")
axes_in[0].set_xlabel("GNSS linear vel y"); axes_in[0].legend()
axes_in[1].plot(gps_data["field.twist.twist.linear.x"], alpha=0.4, label="raw")
axes_in[1].plot(gps["field.twist.twist.linear.x"], label="LPF")
axes_in[1].set_xlabel("GNSS linear vel x"); axes_in[1].legend()
axes_in[2].plot(ins_data["field.angular_velocity.z"], alpha=0.4, label="raw")
axes_in[2].plot(ins_wz_lpf, label="LPF")
axes_in[2].set_xlabel("INS angular vel z"); axes_in[2].legend()
axes_in[3].plot(ins_data["field.linear_acceleration.x"], alpha=0.4, label="raw")
axes_in[3].plot(ins_ax_lpf, label="LPF")
axes_in[3].set_xlabel("INS linear acceleration x"); axes_in[3].legend()
plt.show()