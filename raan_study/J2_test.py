import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def RAAN_from_rv_vec(r, v, tol=1e-10):
    """
    Compute Right Ascension of the Ascending Node (RAAN) from
    inertial position and velocity vectors.

    Parameters
    ----------
    r : array_like, shape (3,) or (N,3)
        Position vector(s) [m]
    v : array_like, shape (3,) or (N,3)
        Velocity vector(s) [m/s]
    tol : float
        Tolerance for detecting equatorial orbits

    Returns
    -------
    RAAN : ndarray
        RAAN angle(s) in radians
        (undefined -> np.nan)
    """

    # Ensure 2D arrays
    r = np.atleast_2d(r)
    v = np.atleast_2d(v)

    # Angular momentum vector: h = r × v
    h = np.cross(r, v)

    # Node vector: n = k × h
    k = np.array([0.0, 0.0, 1.0])
    n = np.cross(k, h)

    n_norm = np.linalg.norm(n, axis=1)

    # Initialize RAAN with NaN (undefined cases)
    RAAN = np.full(r.shape[0], np.nan)

    # Valid where node vector exists (non-equatorial orbit)
    valid = n_norm > tol

    if np.any(valid):
        RAAN[valid] = np.arctan2(n[valid, 1], n[valid, 0])
        RAAN[valid] = np.mod(RAAN[valid], 2 * np.pi)

    # Return scalar if single vector input
    return RAAN[0] if RAAN.size == 1 else RAAN


df = pd.read_csv("../downloads/NbodyWalkerDelta/data/WalkerDelta.csv")

print(df.head())


def deltaRAAN(raan1, raan2):

    return ((raan2 - raan1 + 180) % 360) - 180


def plotRAAN(satname, df):
    dfs = df[df["name"] == satname].copy()
    dfs = dfs.sort_values("time_s")

    # Sample the data
    # sample_interval = 200
    # dfs_sampled = dfs.iloc[::sample_interval].copy()

    rs = dfs[["x", "y", "z"]].to_numpy()
    vs = dfs[["vx", "vy", "vz"]].to_numpy()

    RAAN = RAAN_from_rv_vec(rs, vs)

    RAAN_unwrapped = np.unwrap(RAAN)

    slope = np.gradient(RAAN_unwrapped, dfs["time_s"].values)

    dfsample = pd.DataFrame(
        {"raan": np.degrees(RAAN_unwrapped), "time_s": dfs["time_s"].values}
    )
    """

    rs = dfs_sampled[["x", "y", "z"]].to_numpy()
    vs = dfs_sampled[["vx", "vy", "vz"]].to_numpy()
    dfs_sampled["RAAN"] = RAAN_from_rv_vec(rs, vs)

    # Unwrap to handle 0/360 degree discontinuities
    raan_unwrapped = np.unwrap(dfs_sampled["RAAN"].values)
    raan_deg_unwrapped = np.degrees(raan_unwrapped)

    # Calculate ABSOLUTE change (taking absolute value of the drift)
    raan_change = raan_deg_unwrapped - raan_deg_unwrapped[0]
    raan_abs_change = np.abs(raan_change)

    # Net change over entire period
    net_change_deg = raan_change[-1]  # Signed change
    abs_change_deg = np.abs(net_change_deg)  # Absolute change
    time_days = (dfs_sampled['time_s'].values[-1] - dfs_sampled['time_s'].values[0]) / 86400

    print(f"{satname}: Initial RAAN = {raan_deg_unwrapped[0]:.2f}°, "
          f"Total change = {net_change_deg:.2f}° (|{abs_change_deg:.2f}°|) over {time_days:.2f} days")

    """

    dfsample["time"] = pd.to_datetime(dfsample["time_s"], unit="s")
    dfsample = dfsample.set_index("time")

    dfsample = dfsample["raan"].resample("1D").first()
    # dfsample.plot()

    plt.plot(dfs["time_s"].values, np.degrees(RAAN_unwrapped), label="RAAN")

    plt.ylabel("RAAN [deg]")
    # plt.plot(dfs['time_s'].values,np.degrees(dfs['RAAN'].values))
    # Plot ABSOLUTE change over time
    # time_days_array = (dfs_sampled['time_s'].values - dfs_sampled['time_s'].values[0]) / 86400
    # plt.plot(time_days_array, raan_abs_change, label=satname, marker='o', markersize=3)
    # plt.plot(time_days_array,raan_deg_unwrapped,label=satname,marker="o",markersize=3)
    # plt.plot(dfs["time_s"], dfs["RAAN"])
    plt.xlabel("Time (days)")
    # plt.ylabel('|RAAN Change| (degrees)')
    plt.title("Absolute RAAN Precession over Time")
    plt.legend()
    plt.grid(True)

    return


"""
def plotRAAN(satname, df):
    dfs = df[df["name"] == satname].copy()
    dfs = dfs.sort_values("time_s")

    rs = dfs[["x", "y", "z"]].to_numpy()
    vs = dfs[["vx", "vy", "vz"]].to_numpy()

    # Diagnostics
    r_mag = np.linalg.norm(rs, axis=1)
    orbital_period = (
        2 * np.pi * np.sqrt((r_mag.mean() / 1000) ** 3 / 398600.4418)
    )  # seconds
    time_steps = np.diff(dfs["time_s"].values)

    print(f"\n=== Diagnostics for {satname} ===")
    print(f"Data points: {len(dfs)}")
    print(f"Time span: {(dfs['time_s'].max() - dfs['time_s'].min()) / 86400:.1f} days")
    print(
        f"Avg time step: {time_steps.mean():.1f} s ({time_steps.mean() / 60:.1f} min)"
    )
    print(f"Min/Max time step: {time_steps.min():.1f} / {time_steps.max():.1f} s")
    print(f"Orbital period: ~{orbital_period / 60:.1f} minutes")
    print(f"Altitude: {r_mag.mean() / 1000 - 6371:.1f} km")
    print(f"r_mag mean: {r_mag.mean()} m = {r_mag.mean() / 1000} km")
    print(f"r_mag values (first 5): {r_mag[:5]}")

    # Calculate inclination to understand J2 behavior
    h = np.cross(rs, vs)
    h_mag = np.linalg.norm(h, axis=1)
    inc = np.arccos(h[:, 2] / h_mag)
    print(f"Inclination: {np.degrees(inc.mean()):.2f}°")

    RAAN = RAAN_from_rv_vec(rs, vs)

    # Check for NaNs
    nan_count = np.isnan(RAAN).sum()
    if nan_count > 0:
        print(f"WARNING: {nan_count} NaN values in RAAN!")

    RAAN_unwrapped = np.unwrap(RAAN)

    # Plot both raw and daily-sampled
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

    # Raw data (possibly subsampled for visibility)
    sample_rate = max(1, len(dfs) // 5000)
    time_days = (dfs["time_s"].values - dfs["time_s"].values[0]) / 86400
    ax1.plot(
        time_days[::sample_rate],
        np.degrees(RAAN_unwrapped[::sample_rate]),
        ".",
        markersize=1,
        alpha=0.5,
        label="Raw data",
    )
    ax1.set_ylabel("RAAN [deg]")
    ax1.set_xlabel("Time (days)")
    ax1.set_title(f"{satname} - Raw RAAN")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Daily averaged
    dfsample = pd.DataFrame(
        {"raan": np.degrees(RAAN_unwrapped), "time_s": dfs["time_s"].values}
    )
    dfsample["time"] = pd.to_datetime(dfsample["time_s"], unit="s")
    dfsample = dfsample.set_index("time")
    dfsample_daily = dfsample["raan"].resample("1D").first()

    dfsample_daily.plot(ax=ax2)
    ax2.set_ylabel("RAAN [deg]")
    ax2.set_xlabel("Time")
    ax2.set_title(f"{satname} - Daily Sampled RAAN")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    return fig

"""

# plotRAAN("service_1", df)
plotRAAN("sat_10", df)
plotRAAN("sat_40", df)
plotRAAN("sat_100", df)

plt.legend()
plt.show()
