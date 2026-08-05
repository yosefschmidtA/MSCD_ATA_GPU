#!/usr/bin/env python3
"""Figura V0 x V5: a curva nao muda, o tempo cai -- e o ganho vive em np=1.

Uso: python3 baseline/figura-v5.py
Le baseline/escala.csv (linhas versao=V0 e versao=V5) e gera baseline/v0-v5.png.

DUAS JANELAS, E ISSO NAO E' DETALHE. O V0 so' foi medido em 2026-08-04; o V5 em
2026-08-05. A deriva entre janelas neste projeto ja chegou a 5% com o MESMO
binario -- medida com o braco de controle da campanha de binding: 132,45 s (V4,
janela da madrugada) contra 125,69 s (V5 com o OpenMP neutralizado pelo binding,
janela da tarde). Entao o ganho desenhado aqui e' um LIMITE SUPERIOR, e o rodape
da figura diz isso. Remedir o V0 na janela do V5 e' o que fecharia a conta.

Paleta: azul #2a78d6 (V0) e laranja #eb6834 (V5), os dois primeiros slots
categoricos ja validados neste projeto. Reconferido com o validate_palette.js do
skill dataviz, modo light sobre #fcfcfb, --pairs all:
  CVD       dE 24,7 (protanopia) · 32,7 (tritanopia)   PASS
  normal    dE 33,6                                     PASS
  faixa de luminancia, piso de croma, contraste         PASS
Nao troque a cor sem rodar o validador.
"""
import os, csv
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
COR = {"V0": "#2a78d6", "V5": "#eb6834"}
ROTULO = {"V0": "V0 base", "V5": "V5"}
TINTA, TINTA2 = "#1a1a19", "#5c5b55"
SUP = "#fcfcfb"

todas = defaultdict(list)
janelas = defaultdict(set)
with open(os.path.join(AQUI, "escala.csv")) as f:
    for r in csv.DictReader(f):
        ver = r["versao"]
        if ver not in COR:
            continue
        todas[(ver, int(r["np"]))].append(float(r["wall"]))
        janelas[ver].add(r["janela"])

versoes = [v for v in ("V0", "V5") if v in janelas]
if len(versoes) < 2:
    raise SystemExit("faltam linhas V0 ou V5 no escala.csv -- rode "
                     "baseline/campanha-v5-escala.sh")

# So' os np medidos nas DUAS versoes: um ponto solto de um lado so' viraria
# um degrau falso na linha.
nps = sorted({n for (v, n) in todas if v == "V0"} &
             {n for (v, n) in todas if v == "V5"})
minimo = {v: [min(todas[(v, n)]) for n in nps] for v in versoes}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6), dpi=200)
fig.patch.set_facecolor(SUP)

# --- painel A: tempo absoluto ------------------------------------------
for v in versoes:
    cor = COR[v]
    for n in nps:                                  # medicoes individuais
        for w in todas[(v, n)]:
            ax1.plot([n], [w], marker="o", ms=4, mfc=cor, mec="none",
                     alpha=0.28, ls="none", zorder=2)
    ax1.plot(nps, minimo[v], color=cor, lw=2, marker="o", ms=7,
             mec=SUP, mew=1.4, label=ROTULO[v], zorder=3)
    ax1.annotate(ROTULO[v], xy=(nps[-1], minimo[v][-1]),
                 xytext=(6, 0), textcoords="offset points",
                 color=cor, fontsize=9.5, fontweight="bold", va="center")

ax1.axvline(6, color="#d8d7d0", lw=1, ls="--", zorder=1)
# No topo, nao embaixo: com as duas linhas juntas em ~39 s o rodape do eixo esta
# ocupado, e o rotulo caia em cima da linha do V5.
ax1.text(6.15, ax1.get_ylim()[1] * 0.985, "6 núcleos físicos",
         fontsize=8, color=TINTA2, va="top")
ax1.set_xlabel("ranks MPI (-np)", color=TINTA)
ax1.set_ylabel("tempo de parede (s)", color=TINTA)
ax1.set_title("Tempo total", color=TINTA, fontsize=12, loc="left", pad=16)
ax1.text(0.0, 1.005, "linha = mínimo · pontos claros = medições individuais",
         transform=ax1.transAxes, color=TINTA2, fontsize=8.5, va="bottom")
ax1.set_xticks(nps)
ax1.set_xlim(min(nps) - 0.4, max(nps) + 1.6)
ax1.legend(frameon=False, fontsize=9, labelcolor=TINTA)

# --- painel B: ganho do V5 sobre o V0 ----------------------------------
ganho = [minimo["V0"][i] / minimo["V5"][i] for i in range(len(nps))]
ax2.axhline(1.0, color="#d8d7d0", lw=1, zorder=1)
ax2.bar(range(len(nps)), ganho, color=COR["V5"], width=0.62, zorder=2)
for i, g in enumerate(ganho):
    ax2.text(i, g + 0.03, f"{g:.2f}×", ha="center", fontsize=9, color=TINTA2)
ax2.set_xticks(range(len(nps)))
ax2.set_xticklabels([str(n) for n in nps])
ax2.set_xlabel("ranks MPI (-np)", color=TINTA)
ax2.set_ylabel("V0 / V5", color=TINTA)
ax2.set_title("Ganho sobre a base", color=TINTA, fontsize=12, loc="left",
              pad=16)
# O ganho aqui e' o ACUMULADO das quatro otimizacoes (V1/V2 no symtrivert, V4 no
# csum, V5 no pathcut), nao o efeito isolado do V5 -- que vale ~4% e so' em np=1.
ax2.text(0.0, 1.005, "acumulado das quatro otimizações · cresce com o np: "
         "quanto menor o laço, mais pesa o preparo serial",
         transform=ax2.transAxes, color=TINTA2, fontsize=8.5, va="bottom")
ax2.set_ylim(0, max(ganho) * 1.22)

for ax in (ax1, ax2):
    ax.set_facecolor(SUP)
    ax.grid(axis="y", color="#ebeae3", lw=0.8)
    ax.set_axisbelow(True)
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    for lado in ("left", "bottom"):
        ax.spines[lado].set_color("#d8d7d0")
    ax.tick_params(colors=TINTA2, labelsize=9)

jan0 = ", ".join(sorted(janelas["V0"]))
jan5 = ", ".join(sorted(janelas["V5"]))
fig.text(0.5, -0.005,
         f"V0 medido em {jan0}; V5 em {jan5}. O controle na janela do V5 "
         f"(mesmo binário, OpenMP neutralizado) deu 125,69 s em np=1 contra "
         f"132,45 s do V4 na janela anterior — ~5% de deriva. "
         f"O ganho acima é um limite superior.",
         ha="center", fontsize=8, color=TINTA2)

fig.tight_layout()
saida = os.path.join(AQUI, "v0-v5.png")
fig.savefig(saida, facecolor=fig.get_facecolor(), bbox_inches="tight")

print("np   " + "".join(f"{ROTULO[v]:>9}" for v in versoes) + "    ganho")
for i, n in enumerate(nps):
    print(f"{n:<5}" + "".join(f"{minimo[v][i]:9.2f}" for v in versoes)
          + f"{ganho[i]:9.2f}×")
print("\ngravado:", saida)
