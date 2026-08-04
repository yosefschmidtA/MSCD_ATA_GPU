#!/usr/bin/env python3
"""Figura antes/depois: a curva nao muda, o tempo cai.

Uso: python3 baseline/figura.py
Gera baseline/antes-depois.png a partir de baseline/saida.txt (V0) e
baseline/v1/saida.txt (V1).
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
V0, V1, V2 = "#2a78d6", "#eb6834", "#1baf7a"   # slots categoricos 1-3, validados
TINTA, TINTA2 = "#1a1a19", "#5c5b55"

def ler(caminho):
    """Devolve [(phi, chical)] por bloco de theta."""
    blocos, atual = [], []
    for linha in open(caminho):
        campos = linha.split()
        if len(campos) == 5:
            try:
                phi, _inten, _fundo, chical, _exp = map(float, campos)
            except ValueError:
                continue
            atual.append((phi, chical))
        elif atual:
            blocos.append(atual); atual = []
    if atual:
        blocos.append(atual)
    return blocos

b0 = ler(os.path.join(AQUI, "saida.txt"))
b1 = ler(os.path.join(AQUI, "v1", "saida.txt"))
b2 = ler(os.path.join(AQUI, "v2", "saida.txt"))

todos0 = [v for bl in b0 for _, v in bl]
todos1 = [v for bl in b1 for _, v in bl]
todos2 = [v for bl in b2 for _, v in bl]
assert len(todos0) == len(todos1) == len(todos2)
residuo = max(max(abs(a - b) for a, b in zip(todos0, t))
              for t in (todos1, todos2))
npontos = len(todos0)

fases  = ["A dedup\n$natoms^3$", "B ordena\nbaldes", "C indexa\n$natoms^3$",  # noqa
          "symdblvert", "precutable", "sendjobs"]
t0     = [6.75, 1.02, 13.09, 0.97, 11.20, 0.88]
t1     = [7.43, 1.07,  0.19, 0.90, 10.89, 0.98]
t2     = [1.38, 0.05,  0.16, 0.92,  9.67, 0.94]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4), dpi=200)
fig.patch.set_facecolor("#fcfcfb")

# --- painel A: a curva -------------------------------------------------
bloco = 0
phi0 = [p for p, _ in b0[bloco]]
chi0 = [v for _, v in b0[bloco]]
phi1 = [p for p, _ in b1[bloco]][1::4]
chi1 = [v for _, v in b1[bloco]][1::4]
phi2 = [p for p, _ in b2[bloco]][3::4]
chi2 = [v for _, v in b2[bloco]][3::4]

ax1.set_facecolor("#fcfcfb")
ax1.axhline(0, color="#d8d7d0", lw=1, zorder=1)
ax1.plot(phi0, chi0, color=V0, lw=2, label="V0 base", zorder=2)
ax1.plot(phi1, chi1, ls="none", marker="o", ms=5.5, mfc=V1, mec="#fcfcfb",
         mew=1.2, label="V1", zorder=3)
ax1.plot(phi2, chi2, ls="none", marker="s", ms=5.5, mfc=V2, mec="#fcfcfb",
         mew=1.2, label="V2", zorder=4)
ax1.set_xlabel("$\\phi$ (graus)", color=TINTA)
ax1.set_ylabel("$\\chi$ calculado", color=TINTA)
ax1.set_title("A curva não muda", color=TINTA, fontsize=12, loc="left", pad=10)
resid = "0 exato" if residuo == 0 else f"{residuo:.1e}"
ax1.text(0.0, 1.005, f"$\\theta$ = 18$^\\circ$ · resíduo máximo "
         f"{resid} nos {npontos} pontos",
         transform=ax1.transAxes, color=TINTA2, fontsize=8.5, va="bottom")
ax1.legend(frameon=False, fontsize=9, labelcolor=TINTA)
ax1.grid(axis="y", color="#ebeae3", lw=0.8)
ax1.set_axisbelow(True)
for lado in ("top", "right"):
    ax1.spines[lado].set_visible(False)
for lado in ("left", "bottom"):
    ax1.spines[lado].set_color("#d8d7d0")
ax1.tick_params(colors=TINTA2, labelsize=9)

# --- painel B: o tempo -------------------------------------------------
y = range(len(fases))
alt = 0.26
ax2.set_facecolor("#fcfcfb")
ax2.barh([i - alt - 0.015 for i in y], t0, height=alt, color=V0,
         label="V0 base", zorder=2)
ax2.barh([i for i in y], t1, height=alt, color=V1, label="V1", zorder=2)
ax2.barh([i + alt + 0.015 for i in y], t2, height=alt, color=V2,
         label="V2", zorder=2)
for i, (a, b, c) in enumerate(zip(t0, t1, t2)):
    ax2.text(a + 0.25, i - alt - 0.015, f"{a:.2f}", va="center",
             fontsize=8, color=TINTA2)
    ax2.text(b + 0.25, i, f"{b:.2f}", va="center", fontsize=8, color=TINTA2)
    ax2.text(c + 0.25, i + alt + 0.015, f"{c:.2f}", va="center",
             fontsize=8, color=TINTA2)
ax2.set_yticks(list(y)); ax2.set_yticklabels(fases, fontsize=8.5, color=TINTA)
ax2.invert_yaxis()
ax2.set_xlabel("segundos (rank 0)", color=TINTA)
ax2.set_xlim(0, 15.6)
ax2.set_title("O tempo cai", color=TINTA, fontsize=12, loc="left", pad=10)
ax2.text(0.0, 1.005, "preparo serial 33,9 → 21,5 → 13,0 s · "
         "total 70 → 55,6 → 48,4 s",
         transform=ax2.transAxes, color=TINTA2, fontsize=8.5, va="bottom")
ax2.legend(frameon=False, fontsize=9, labelcolor=TINTA, loc="lower right")
ax2.grid(axis="x", color="#ebeae3", lw=0.8)
ax2.set_axisbelow(True)
for lado in ("top", "right", "left"):
    ax2.spines[lado].set_visible(False)
ax2.spines["bottom"].set_color("#d8d7d0")
ax2.tick_params(colors=TINTA2, labelsize=9)

fig.tight_layout()
saida = os.path.join(AQUI, "antes-depois.png")
fig.savefig(saida, facecolor=fig.get_facecolor())
print(f"residuo maximo em {npontos} pontos: {residuo}")
print(f"gravado: {saida}")
