from itertools import combinations

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation


def RAAN_from_rv_vec(r, v, tol=1e-10):
    """
    Compute Right Ascension of the Ascending Node (RAAN)
    from inertial position and velocity vectors.
    """

    r = np.atleast_2d(r)
    v = np.atleast_2d(v)

    # Angular momentum vector
    h = np.cross(r, v)

    # Node vector
    k = np.array([0.0, 0.0, 1.0])
    n = np.cross(k, h)

    n_norm = np.linalg.norm(n, axis=1)

    RAAN = np.full(r.shape[0], np.nan)

    valid = n_norm > tol

    if np.any(valid):
        RAAN[valid] = np.arctan2(n[valid, 1], n[valid, 0])
        RAAN[valid] = np.mod(RAAN[valid], 2 * np.pi)

    return RAAN[0] if RAAN.size == 1 else RAAN


def animate_raan_distance_network(df, satellite_names, interval=500):

    df = df[df["name"].isin(satellite_names)].copy()
    times = sorted(df["time_s"].unique())

    fig, ax = plt.subplots(figsize=(8, 8))
    frame_state = {"frame": 0}

    def on_key(event):
        if event.key == "right":
            frame_state["frame"] = (frame_state["frame"] + 1) % len(times)
            update(frame_state["frame"])
            plt.draw()

        elif event.key == "left":
            frame_state["frame"] = (frame_state["frame"] - 1) % len(times)
            update(frame_state["frame"])
            plt.draw()

    def angular_diff_rad(a, b):
        """Smallest angular difference"""
        diff = abs(a - b) % (2 * np.pi)
        return min(diff, 2 * np.pi - diff)

    def update(frame):
        ax.clear()

        t = times[frame]
        current = df[df["time_s"] == t]

        rs = current[["x", "y", "z"]].to_numpy()
        vs = current[["vx", "vy", "vz"]].to_numpy()

        RAAN = RAAN_from_rv_vec(rs, vs)
        raan_dict = dict(zip(current["name"], RAAN))

        G = nx.Graph()
        G.add_nodes_from(satellite_names)

        for sat1, sat2 in combinations(satellite_names, 2):
            if sat1 in raan_dict and sat2 in raan_dict:
                diff_rad = angular_diff_rad(
                    raan_dict[sat1],
                    raan_dict[sat2]
                )
                diff_deg = np.degrees(diff_rad)

                # Invert weight so larger RAAN → longer edge
                layout_weight = 1.0 / (diff_deg + 1e-6)

                G.add_edge(
                    sat1,
                    sat2,
                    weight=layout_weight,  # layout physics
                    label=diff_deg         # displayed value
                )

        # Layout
        pos = nx.spring_layout(G, weight="weight", seed=42)

        # Edge colors (blue if connected to service_1)
        edge_colors = [
            "blue" if (u == "service_1" or v == "service_1") else "gray"
            for u, v in G.edges()
        ]

        # Draw graph
        nx.draw(
            G,
            pos,
            with_labels=True,
            node_size=800,
            node_color="lightblue",
            edge_color=edge_colors,
            ax=ax,
        )

        # Edge labels (use label, not weight)
        edge_labels = nx.get_edge_attributes(G, "label")
        edge_labels = {k: f"{v:.1f}°" for k, v in edge_labels.items()}

        nx.draw_networkx_edge_labels(
            G,
            pos,
            edge_labels=edge_labels,
            ax=ax
        )

        ax.set_title(f"RAAN Distance Network (t = {t})")
        ax.set_axis_off()

    ani = FuncAnimation(
        fig,
        update,
        frames=len(times),
        interval=1e9  # effectively paused
    )
    ani.event_source.stop()

    fig.canvas.mpl_connect("key_press_event", on_key)

    update(0)
    plt.show()

    return ani


# ------------------------
# Load data
# ------------------------

df = pd.read_csv("../downloads/NbodyWalkerDelta/data/WalkerDelta.csv")

sat_list = ["service_1", "sat_1", "sat_36", "sat_100", "sat_200"]

animate_raan_distance_network(df, sat_list, interval=300)
