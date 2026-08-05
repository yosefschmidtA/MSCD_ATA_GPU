#!/usr/bin/env python3
"""Curva de escalabilidade V0 x V2 x V4.

Le baseline/escala.csv (janela,versao,np,rep,wall,factors) e gera
baseline/escalabilidade.png.

A linha e' o MINIMO das repeticoes -- estimador padrao em benchmark, porque
ruido so' atrasa, nunca acelera. Os pontos claros sao as medicoes individuais,
para que a dispersao fique visivel em vez de escondida.

DUAS JANELAS DE MEDICAO. O V0 so' foi medido em 2026-08-04; V2 foi medido nas
duas; V4 so' em 2026-08-05. O script usa a janela mais recente onde a versao
existe, e imprime a deriva do V2 entre as duas janelas -- que e' o que diz se
comparar o V4 de hoje com o V0 de ontem e' legitimo. Se a deriva for grande,
a comparacao com V0 nao vale e o texto tem de dizer isso.

Paleta validada com o validate_palette.js do skill dataviz, modo light sobre
a superficie #fcfcfb, criterio --pairs all (com 3 linhas simultaneas qualquer
par pode ser confundido, nao so' os adjacentes):
  azul #2a78d6 / laranja #eb6834 / carmim #c2255c
  pior par: carmim-laranja, dE 16,0 (deuteranopia) e 17,1 (visao normal).
Verde e roxo foram testados e REPROVARAM: verde-laranja da dE 6,1 em
protanopia, e roxo-azul da 5,5 em deuteranopia (o roxo #7048e8 reprova ate
para visao normal, dE 12,9). Nao troque a cor sem rodar o validador.
"""
import os, csv
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
COR = {"V0": "#2a78d6", "V2": "#eb6834", "V4": "#c2255c"}
ROTULO = {"V0": "V0 base", "V2": "V2", "V4": "V4"}
TINTA, TINTA2 = "#1a1a19", "#5c5b55"

todas = defaultdict(list)
janelas = defaultdict(set)
with open(os.path.join(AQUI, "escala.csv")) as f:
    for r in csv.DictReader(f):
        jan, ver, np_ = r["janela"], r["versao"], int(r["np"])
        todas[(jan, ver, np_)].append(float(r["wall"]))
        janelas[ver].add(jan)

nps = sorted({k[2] for k in todas})
recente = {v: max(janelas[v]) for v in janelas}          # janela usada por versao
versoes = [v for v in ("V0", "V2", "V4") if v in janelas]

minimo = {v: [min(todas[(recente[v], v, n)]) for n in nps] for v in versoes}

# --- deriva do V2 entre as duas janelas --------------------------------
deriva = None
if len(janelas.get("V2", ())) > 1:
    velha, nova = sorted(janelas["V2"])[0], sorted(janelas["V2"])[-1]
    pares = [(n, min(todas[(velha, "V2", n)]), min(todas[(nova, "V2", n)]))
             for n in nps if (velha, "V2", n) in todas and (nova, "V2", n) in todas]
    if pares:
        piora = [100.0 * (b - a) / a for _, a, b in pares]
        deriva = (velha, nova, pares, sum(piora) / len(piora),
                  max(piora, key=abs))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6), dpi=200)
fig.patch.set_facecolor("#fcfcfb")

# --- painel A: tempo absoluto ------------------------------------------
for v in versoes:
    cor = COR[v]
    for n in nps:                                  # medicoes individuais
        for w in todas[(recente[v], v, n)]:
            ax1.plot([n], [w], marker="o", ms=4, mfc=cor, mec="none",
                     alpha=0.28, ls="none", zorder=2)
    ax1.plot(nps, minimo[v], color=cor, lw=2, marker="o", ms=7,
             mec="#fcfcfb", mew=1.4, label=f"{ROTULO[v]} ({recente[v][5:]})",
             zorder=3)
    # rotulo direto no fim da linha: identidade nao fica so' na cor
    ax1.annotate(ROTULO[v], xy=(nps[-1], minimo[v][-1]),
                 xytext=(6, 0), textcoords="offset points",
                 color=cor, fontsize=9.5, fontweight="bold", va="center")

ax1.axvline(6, color="#d8d7d0", lw=1, ls="--", zorder=1)
# rotulo embaixo: em cima ele colide com a caixa de legenda
ax1.text(6.15, ax1.get_ylim()[0] + 1.5, "6 núcleos físicos",
         fontsize=8, color=TINTA2, va="bottom")
ax1.set_xlabel("ranks MPI (-np)", color=TINTA)
ax1.set_ylabel("tempo de parede (s)", color=TINTA)
ax1.set_title("Tempo total", color=TINTA, fontsize=12, loc="left", pad=16)
ax1.text(0.0, 1.005, "linha = mínimo · pontos claros = medições individuais",
         transform=ax1.transAxes, color=TINTA2, fontsize=8.5, va="bottom")
ax1.set_xticks(nps)
ax1.set_xlim(min(nps) - 0.4, max(nps) + 1.6)
ax1.legend(frameon=False, fontsize=9, labelcolor=TINTA)

# --- painel B: ganho sobre o V0 ----------------------------------------
alvos = [v for v in versoes if v != "V0"]
larg = 0.8 / max(len(alvos), 1)
ax2.axhline(1.0, color="#d8d7d0", lw=1, zorder=1)
for k, v in enumerate(alvos):
    ganho = [minimo["V0"][i] / minimo[v][i] for i in range(len(nps))]
    desl = (k - (len(alvos) - 1) / 2.0) * (larg + 0.03)
    pos = [i + desl for i in range(len(nps))]
    ax2.bar(pos, ganho, color=COR[v], width=larg, zorder=2, label=ROTULO[v])
    for i, g in enumerate(ganho):
        ax2.text(pos[i], g + 0.03, f"{g:.2f}×", ha="center", fontsize=8,
                 color=TINTA2)
ax2.set_xticks(range(len(nps)))
ax2.set_xticklabels([str(n) for n in nps])
ax2.set_xlabel("ranks MPI (-np)", color=TINTA)
ax2.set_ylabel("V0 / versão", color=TINTA)
ax2.set_title("Ganho sobre a base", color=TINTA, fontsize=12, loc="left",
              pad=10)
ax2.text(0.0, 1.005, "cresce com o np: quanto menor o laço paralelo, "
         "mais pesa o preparo serial",
         transform=ax2.transAxes, color=TINTA2, fontsize=8.5, va="bottom")
ax2.legend(frameon=False, fontsize=9, labelcolor=TINTA, loc="upper left")
ax2.set_ylim(0, None)

for ax in (ax1, ax2):
    ax.set_facecolor("#fcfcfb")
    ax.grid(axis="y", color="#ebeae3", lw=0.8)
    ax.set_axisbelow(True)
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    for lado in ("left", "bottom"):
        ax.spines[lado].set_color("#d8d7d0")
    ax.tick_params(colors=TINTA2, labelsize=9)

if deriva:
    velha, nova, _, media, pior = deriva
    fig.text(0.5, -0.005,
             f"V0 medido em {velha}; V2 e V4 em {nova}. "
             f"Deriva do V2 entre as janelas: {media:+.1f}% em média, "
             f"{pior:+.1f}% no pior np.",
             ha="center", fontsize=8, color=TINTA2)

fig.tight_layout()
saida = os.path.join(AQUI, "escalabilidade.png")
fig.savefig(saida, facecolor=fig.get_facecolor(), bbox_inches="tight")

cab = "np   " + "".join(f"{ROTULO[v]:>9}" for v in versoes)
print(cab)
for i, n in enumerate(nps):
    print(f"{n:<5}" + "".join(f"{minimo[v][i]:9.2f}" for v in versoes))
if deriva:
    velha, nova, pares, media, pior = deriva
    print(f"\nderiva do V2, {velha} -> {nova}: "
          f"media {media:+.1f}%, pior {pior:+.1f}%")
    for n, a, b in pares:
        print(f"  np={n:<3} {a:7.2f} -> {b:7.2f}  {100*(b-a)/a:+6.1f}%")
print("\ngravado:", saida)
