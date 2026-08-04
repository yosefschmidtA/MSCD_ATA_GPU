#!/usr/bin/env python3
"""Curva de escalabilidade V0 x V2.

Le baseline/escala.csv (versao,np,rep,wall,factors) e gera
baseline/escalabilidade.png.

A linha e' o MINIMO das repeticoes -- estimador padrao em benchmark, porque
ruido so' atrasa, nunca acelera. Os pontos claros sao as medicoes individuais,
para que a dispersao fique visivel em vez de escondida.
"""
import os, csv
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
V0C, V2C = "#2a78d6", "#eb6834"
TINTA, TINTA2 = "#1a1a19", "#5c5b55"

todas = defaultdict(list)
with open(os.path.join(AQUI, "escala.csv")) as f:
    for r in csv.DictReader(f):
        todas[(r["versao"], int(r["np"]))].append(float(r["wall"]))

nps = sorted({k[1] for k in todas})
minimo = {v: [min(todas[(v, n)]) for n in nps] for v in ("V0", "V2")}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4), dpi=200)
fig.patch.set_facecolor("#fcfcfb")

# --- painel A: tempo absoluto ------------------------------------------
for v, cor, rot in (("V0", V0C, "V0 base"), ("V2", V2C, "V2 otimizado")):
    for n in nps:                                  # medicoes individuais
        for w in todas[(v, n)]:
            ax1.plot([n], [w], marker="o", ms=4, mfc=cor, mec="none",
                     alpha=0.28, ls="none", zorder=2)
    ax1.plot(nps, minimo[v], color=cor, lw=2, marker="o", ms=7,
             mec="#fcfcfb", mew=1.4, label=rot, zorder=3)

ax1.axvline(6, color="#d8d7d0", lw=1, ls="--", zorder=1)
ax1.text(6.15, ax1.get_ylim()[1] * 0.985, "6 núcleos físicos",
         fontsize=8, color=TINTA2, va="top")
ax1.set_xlabel("ranks MPI (-np)", color=TINTA)
ax1.set_ylabel("tempo de parede (s)", color=TINTA)
ax1.set_title("Tempo total", color=TINTA, fontsize=12, loc="left", pad=10)
ax1.text(0.0, 1.005, "linha = mínimo · pontos claros = medições individuais",
         transform=ax1.transAxes, color=TINTA2, fontsize=8.5, va="bottom")
ax1.set_xticks(nps)
ax1.legend(frameon=False, fontsize=9, labelcolor=TINTA)

# --- painel B: ganho do V2 ---------------------------------------------
ganho = [minimo["V0"][i] / minimo["V2"][i] for i in range(len(nps))]
ax2.axhline(1.0, color="#d8d7d0", lw=1, zorder=1)
ax2.bar([str(n) for n in nps], ganho, color=V2C, width=0.55, zorder=2)
for i, g in enumerate(ganho):
    ax2.text(i, g + 0.02, f"{g:.2f}×", ha="center", fontsize=9, color=TINTA2)
ax2.set_ylim(0, max(ganho) * 1.22)
ax2.set_xlabel("ranks MPI (-np)", color=TINTA)
ax2.set_ylabel("V0 / V2", color=TINTA)
ax2.set_title("Ganho do V2 sobre a base", color=TINTA, fontsize=12,
              loc="left", pad=10)
ax2.text(0.0, 1.005, "cresce com o np: quanto menor o laço paralelo, "
         "mais pesa o preparo serial",
         transform=ax2.transAxes, color=TINTA2, fontsize=8.5, va="bottom")

for ax, eixo in ((ax1, "y"), (ax2, "y")):
    ax.set_facecolor("#fcfcfb")
    ax.grid(axis=eixo, color="#ebeae3", lw=0.8)
    ax.set_axisbelow(True)
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    for lado in ("left", "bottom"):
        ax.spines[lado].set_color("#d8d7d0")
    ax.tick_params(colors=TINTA2, labelsize=9)

fig.tight_layout()
saida = os.path.join(AQUI, "escalabilidade.png")
fig.savefig(saida, facecolor=fig.get_facecolor())
for i, n in enumerate(nps):
    print(f"np={n:2d}  V0={minimo['V0'][i]:6.2f}  V2={minimo['V2'][i]:6.2f}"
          f"  ganho={ganho[i]:.2f}x")
print("gravado:", saida)
