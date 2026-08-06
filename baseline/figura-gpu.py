#!/usr/bin/env python3
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
SUP = "#fcfcfb"
TINTA, TINTA2 = "#1a1a19", "#5c5b55"

nps = [1, 2, 4, 6, 8, 10, 12]
v0 = [146.49, 87.53, 64.26, 70.97, 69.06, None, 69.70]
v5 = [122.02, 69.82, 40.50, 39.50, 38.84, None, 38.77]
gpu = [35.65, 26.28, 24.43, 26.71, 25.08, 30.87, 37.79]

fig, ax1 = plt.subplots(1, 1, figsize=(8, 4.6), dpi=200)
fig.patch.set_facecolor(SUP)
ax1.set_facecolor(SUP)

ax1.plot([n for i, n in enumerate(nps) if v0[i]], [t for t in v0 if t], color="#2a78d6", lw=2, marker="o", ms=7, label="V0 (CPU Base)")
ax1.plot([n for i, n in enumerate(nps) if v5[i]], [t for t in v5 if t], color="#eb6834", lw=2, marker="o", ms=7, label="V5 (CPU Máx)")
ax1.plot(nps, gpu, color="#2ea043", lw=2, marker="o", ms=7, label="GPU (Fase 3)")

ax1.axvline(6, color="#d8d7d0", lw=1, ls="--", zorder=1)
ax1.text(6.15, 140, "6 núcleos físicos", fontsize=8, color=TINTA2, va="top")

ax1.set_xlabel("ranks MPI (-np)", color=TINTA)
ax1.set_ylabel("tempo de parede (s)", color=TINTA)
ax1.set_title("Otimização CPU vs GPU", color=TINTA, fontsize=12, loc="left", pad=16)
ax1.set_xticks(nps)
ax1.set_ylim(0, 150)
ax1.legend(frameon=False, fontsize=10, labelcolor=TINTA)

ax1.grid(axis="y", color="#ebeae3", lw=0.8)
ax1.set_axisbelow(True)
for lado in ("top", "right"):
    ax1.spines[lado].set_visible(False)
for lado in ("left", "bottom"):
    ax1.spines[lado].set_color("#d8d7d0")
ax1.tick_params(colors=TINTA2, labelsize=9)

fig.tight_layout()
saida = os.path.join(AQUI, "gpu-v0-v5.png")
fig.savefig(saida, facecolor=fig.get_facecolor(), bbox_inches="tight")
print("gravado:", saida)
