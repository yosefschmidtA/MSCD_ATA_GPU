# Mapa equação ↔ código

As equações são as do `MANUAL_MSCD_CEA.pdf` (numeração original do manual, com a
página onde aparecem). Os arquivos são deste diretório.

**Para que serve:** o código chama as coisas de `cxa`, `algam`, `tevenelem`,
`xb`. O manual chama de Γ, F, t_l, W_c. Quem lê um não reconhece o outro. Esta
ponte não está em lugar nenhum dos dois — nem no grafo — e some se não for
escrita.

**Atenção ao arquivo certo:** o `makefile` linka **`mscdrunb_not_reanalize.o`**,
não `mscdrunb.o`. As linhas citadas aqui para `symtrivert`/`symdblvert` são as do
arquivo que de fato compila.

---

## Índice rápido

| Eq. | pág. | o que é | onde está |
|---|---|---|---|
| (1) | 3 | I ∝ \|φ₀ + Σφ_sj\|² | `mscdrund.cpp:149` |
| (2) | 4 | propagador livre G_LL'(ρ) | implícito — nunca montado (ver eq. 6) |
| (3) | 4, 9 | série exata de espalhamento múltiplo | substituída pela eq. (12) |
| (4) | 9 | intensidade de fotoemissão | `mscdrund.cpp:128-151` |
| (5) | 10 | elemento de matriz dipolar m_lf,c | `mscdrune.cpp:20` `matrixelement()` |
| (6) | 10 | propagador separável de Rehr–Albers | base de tudo; `radim` = nº de λ |
| (7)(8) | 10 | Γ^L_λ(ρ) e Γ̃^L_λ(ρ) | `msfuncs.cpp:108` `Hankel::makecurve` |
| (9) | 10 | normalização N_lμ | dentro de `Hankel::makecurve` |
| (10) | 10 | C_l^(ν)(z) = dᵛ/dzᵛ C_l(z) | `msfuncs.cpp:108` |
| (11) | 10 | rotação dos harmônicos esféricos | `rotamat.cpp:175` |
| (12) | 10 | série R-A separável | `mscdrund.cpp:59-126` |
| (13) | 11 | **matriz de amplitude F_λλ'(ρ,ρ')** | `mscdrunc.cpp:134` `alltrievent()` |
| (14) | 13 | R^j_mm'(α,β,γ) = e^(-imα) r^j_mm'(β) e^(-im'γ) | `rotamat.cpp:175` |
| (15) | 13 | fórmula de Wigner para r^j_mm'(β) | `rotamat.cpp:175` |
| (16) | 13 | simetrias de r^j_mm' | `rotamat.cpp:175` |
| (17) | 13 | r¹(β) explícita (l=1) | `rotamat.cpp:175` |
| (18) | 14 | F com matriz de rotação composta | `mscdrunc.cpp:134` |
| (19) | 14 | matriz de rotação composta R^l_μμ'(Ω_ρρ') | `mscdrunc.cpp:712` `allrotation()` |
| (20) | 14 | γ^l_μν e γ̃^l_μν | `msfuncs.cpp:108` |
| (21)(22) | 14, 15 | derivação dos ângulos compostos | — (álgebra, não vira código) |
| (23) | 15 | **cos β = cosθcosθ' + senθsenθ'cos(φ'−φ)** | `mscdrunb_not_reanalize.cpp:397`; refeita em `mscdrund.cpp:600` `maketripar()` |
| (24) | 15 | **ângulos de Euler (α, β, γ) da rotação composta** | `mscdrunc.cpp:643` `onerotation()`, chamada por `allrotation()` (`:712`) |
| (25) | 15 | sign(1.0, cos β) | dentro de `onerotation()` |
| (26) | 16 | **F com absorção inelástica e vibração** | `mscdrunc.cpp:183-186` |
| (27)(28) | 16 | **fórmula TPP-2 do livre caminho médio** | `meanpath.cpp:143` `finvpath()` |
| (29) | 16 | k = 0,512331·√E | `mscdrund.cpp:648` `kcoeff` |
| (30) | 16 | λ = k·Eᵐ (Wagner–Davis–Riggs) | alternativa em `meanpath.cpp` |
| (31) | 18 | **W_c = exp[−k²(1−cosβ)σ_c²]** | `mscdrunc.cpp:185` |
| (32) | 18 | MSRD σ_c² (Debye correlacionado) | `vibrate.cpp:149` `fvibmsrd()` |
| (33) | 18 | q_D = (6π²N/V)^(1/3) | `vibrate.cpp` |
| (34)(35) | 18 | massa efetiva de superfície | `vibrate.cpp` |
| (36)(37) | 19 | **correção de potencial interno (refração)** | `mscdrund.cpp:646` `kinside()` + `:654` `thetainside()` |
| (38) | 19 | **I_abc = (2I_a + I_b + I_c + I_d + I_e)/6** | `mscdrund.cpp:231-245` |
| (39) | 29 | fator de confiabilidade de intensidade R_I | `pdintena.cpp:408` |
| (41) | 33 | **fator-R: R = Σ(χ_c−χ_e)²/(χ_c²+χ_e²)** | `pdintena.cpp:408` `reliability()` |
| (42)(43) | 34 | forma matricial Γ̃·F·F·Γ | `mscdrund.cpp:59-126` |
| (44)(45) | 36 | somatório iterativo (frente / trás) | `mscdrund.cpp:59` (é a versão "para trás") |
| (46)(47) | 36 | **fatoração que justifica o `symtrivert`** | `mscdrunb_not_reanalize.cpp:348` |

---

## O eixo do programa: eq. (46)–(47)

Tudo no MSCD gira em torno de uma observação só, na página 36 do manual. A
matriz de amplitude de espalhamento se fatora:

> **(46)**  F_λλ'(ρ,ρ') = exp(−iμα) · f_λλ'(ρ,ρ',β) · exp(−iμ'γ)
>
> **(47)**  f_λλ'(ρ,ρ',β) = exp(−a'/2λ − k²(1−cosβ)σ_c²) · [exp(iρ')/ρ'] ·
> Σ_l t_l γ^l_μα(ρ) d^l_μμ'(β) γ̃^l_μ'ν'(ρ')

O caro é `f`, e `f` **não depende da orientação absoluta** do par de vetores —
só de |ρ|, |ρ'| e do ângulo β entre eles. O barato é o par de fases
exp(−iμα)·exp(−iμ'γ), que carrega toda a orientação.

Daí sai a arquitetura inteira do programa, e o motivo de existir uma fase serial
gorda antes do cálculo:

| passo | equação | código | quando roda |
|---|---|---|---|
| enumerar os N³ trios e reduzir aos distintos | (47) | `symtrivert()` | **1×**, serial |
| refazer a assinatura (a', 1/a', 1/a, cos β) | (23) | `maketripar()` | 1× (e a cada passo de ajuste) |
| calcular α e γ de cada trio **real** | (24) | `allrotation()` | 1× (e a cada passo de ajuste) |
| calcular `f` para cada trio **distinto** | (47) | `alltrievent()` | 1× por \|k\| |
| aplicar exp(−iμα−iμ'γ) e somar | (46), (12) | `summation()` | **3895×**, paralelo |

A separação entre as duas linhas do meio é justamente a eq. (46): `allrotation`
percorre os N³ trios reais porque α e γ dependem da orientação absoluta;
`alltrievent` percorre só os 494 297 distintos porque `f` não depende.

Com `natoms=205` são 205³ = 8 615 125 eventos de espalhamento possíveis, que
colapsam em **494 297 geometrias distintas** — fator 17 de economia. O preço é
os ~9 s de `symtrivert` e os ~190 MB de tabela que precisam ser difundidos para
cada rank.

O manual, mesma página: *"we further pre-calculate γ̃(ρ), γ(ρ) and r(β) for a
preselected series of ρ and β values, then use interpolation"*. É o que fazem as
classes com `makecurve()` — `Hankel`, `Rotamat`, `Expix` — todas construídas em
`precutable()`.

---

## A tabela `tevenpar`: a chave de tudo

A assinatura que o `symtrivert` deduplica está declarada num comentário em
`mscdrunb_not_reanalize.cpp:339`:

```
tevenpar   0        1            2             3       4
          lengtha  inverse_lena  inverse_lengb cosbeta atom_kind
          5        6       7     8     9
          eledim   memadd  ia    ib   ic
```

Lendo contra a Fig. 2 do manual (pág. 11), o átomo `a` emite, `b` espalha, `c`
recebe:

| campo | manual | significado |
|---|---|---|
| `tevenpar[j*10+0]` | a' | distância a→b, em Å |
| `tevenpar[j*10+1]` | 1/a' | recíproco, pré-calculado |
| `tevenpar[j*10+2]` | 1/a | recíproco da distância b→c |
| `tevenpar[j*10+3]` | cos β | eq. (23) |
| `tevenpar[j*10+4]` | — | espécie química do espalhador (escolhe t_l) |
| `tevenpar[j*10+5]` | — | `eledim`: dimensão da matriz R-A (1, 3, 6, 10 ou 15) |
| `tevenpar[j*10+6]` | — | deslocamento dentro de `tevenelem` |

Note que `ρ = k·a` e `ρ' = k·a'` do manual são **adimensionais**; a tabela guarda
as distâncias em Å e a multiplicação por `k` acontece em `alltrievent`
(`ka=akin*xa`). Por isso a tabela serve para qualquer energia.

---

## Eq. (26) linha a linha

O trecho mais denso do código é `mscdrunc.cpp:183-186`, dentro de
`alltrievent()`. É a eq. (26)/(47) inteira em quatro linhas:

```cpp
xc = meanpath->finvpath(akin);                              // 1/λ        eq.(27)
xd = vibrate->fvibmsrd(xa, aweight[akind-1]);               // σ_c²       eq.(32)
xb = (float)exp(-0.5*xa*xc - akin*akin*xd*(1.0-cosbeta))/ka;
//        └── −a'/(2λ) ──┘   └── −k²(1−cosβ)σ_c² ──┘        eq.(26), eq.(31)
//                                                    └ 1/ρ'
cxa = xb * expix->fexpix(ka/radian);                        // ×exp(iρ')
```

- `xa` é a' (distância interatômica em Å), `ka = k·a' = ρ'`
- `-0.5*xa*xc` é o decaimento inelástico exp(−a'/2λ) — o fator ½ é porque λ é
  definido para a **intensidade** e aqui se propaga **amplitude**
- o segundo termo é o Debye–Waller correlacionado W_c da eq. (31)
- `expix->fexpix` recebe graus (daí `/radian`): é a tabela interpolada de e^(ix)

E o par de fases da eq. (46), aplicado só na hora do uso, está em
`mscdrund.cpp:112-113` dentro de `summation()`:

```cpp
xc = -p*xb - q*xa;              // p=μ, q=μ', xb=tgamma[id]=γ, xa=talpha[id]=α
algam[t] = expix->fexpix(xc);   // exp(−iμγ − iμ'α)
```

`talpha[]` e `tgamma[]` são dois vetores de `natoms³` floats (34 MB cada) que
guardam os ângulos de Euler da eq. (24) para cada trio **real** — enquanto
`tevenelem` guarda `f` só para os trios **distintos**.

---

## `msorder` e `raorder` não são a mesma coisa

Confusão fácil, e o manual só separa direito na Tabela 4 (pág. 29):

- **`msorder`** = ordem de espalhamento múltiplo = quantos átomos há no caminho.
  É o `for (m=msorder; m>=2; --m)` de `mscdrund.cpp:59`, a série da eq. (12).
  `msorder=8` no `Cov0.txt`.
- **`raorder`** = ordem da expansão separável de Rehr–Albers = **tamanho da
  matriz** F. Vira `radim` em `mscdruna.cpp:472-476`:

  | `raorder` | 0 | 1 | 2 | 3 | 4 |
  |---|---|---|---|---|---|
  | `radim` | 1 | 3 | 6 | 10 | 15 |

  São os tamanhos (1×1), (3×3), (6×6), (10×10), (15×15) da Tabela 1 do manual.
  `raorder=4` no `Cov0.txt` — o máximo.

Os índices λ = (μ,ν) de cada posição da matriz estão nas Tabelas 5–8 do manual
(pág. 35), e aparecem no código como a tabela `lamda[]` construída à mão no
início de `summation()` (`mscdrund.cpp:30-45`) e de `alltrievent()`
(`mscdrunc.cpp:149-164`) — `lamda[j]` é μ e `lamda[32+j]` é ν.

## O `pathcut` e a esparsidade adaptativa

Manual, pág. 29: o maior F₀₀(ρ,ρ') vira referência; todo elemento menor que
`pathcut` vezes esse máximo é declarado desprezível, *"a scattering matrix is
automatically reduced as appropriate to a lower-order R-A event with smaller
matrix size"*.

No código isso é `precutable()` (`mscdrunc.cpp:285`), e o resultado vai para dois
vetores:

- **`tevencut[m][ia][ib]`** — o caminho de ordem `m` por (ia,ib) sobrevive? Se 0,
  `summation()` pula (`mscdrund.cpp:79`).
- **`tevendim[ia][ib][ic]`** — a dimensão da matriz R-A para *cada ordem* `m`,
  empacotada em **4 bits por ordem** dentro de um único `int`
  (`mscdrunc.cpp:432`: `tevendim[id] += (15<<((m-2)*4))`). Na leitura,
  `mscdrund.cpp:83-86` desempacota com `evedim >>= (m-2)*4; evedim &= 15;`.

A referência contra a qual o `pathcut` compara é o `pemeven` de
`mscdrunc.cpp:385-386`: o maior `real(tevenelem[memadd+1])` sobre todos os trios.
Vale entender o que são esses seis primeiros elementos, porque o `tevenelem`
**muda de significado** entre as duas passadas:

1. **Passada de pré-corte.** `precutable` aloca `tevenelem` com `ntrieven*6`
   (`mscdrunc.cpp:339,347`) e `alltrievent(1,kmin)` grava, por trio distinto,
   `[memadd+0]` = o elemento (0,0) da matriz, e `[memadd+1..+5]` = os máximos
   acumulados de \|elemento\| sobre os blocos 1, 3, 6, 10 e 15
   (`mscdrunc.cpp:271-277`). São sondas de magnitude, não a matriz.
2. **Dimensionamento.** `mscdrunc.cpp:462-489` percorre todos os trios reais e
   todas as ordens `m`, e para cada trio *distinto* guarda em
   `tevenpar[k*10+5]` a **maior** dimensão R-A que algum uso dele exigiu. Então
   `ntrielem = Σ eledim²` e `tevenpar[j*10+6]` vira o deslocamento acumulado.
3. **Passada real.** `tevenelem` é realocado com esse `ntrielem`
   (`mscdrunc.cpp:617-624`) e passa a guardar as matrizes de verdade.

Ou seja: o `tevenelem` final é um **vetor irregular** — cada trio distinto ocupa
exatamente `eledim²` complexos, com `eledim` ∈ {1, 3, 6, 10, 15}.

**Consequência para a GPU:** o tamanho da matriz varia por trio *e* por ordem de
espalhamento, e o armazenamento é irregular por construção. Não é um batch
homogêneo; é um batch com cinco formas misturadas e deslocamentos que só se
conhecem depois do pré-corte. Agrupar por `eledim` antes de despachar
provavelmente vale mais que qualquer otimização aritmética.

---

## O laço de 5 pontos: eq. (38)

`accepang = 1.5°` no `Cov0.txt` é o semi-ângulo de aceitação do analisador. O
manual (pág. 19) manda calcular a intensidade em cinco direções da abertura e
combinar com peso duplo no centro:

> **(38)**  I_abc = (2 I_a + I_b + I_c + I_d + I_e) / 6

O código (`mscdrund.cpp:231-245`):

```cpp
for (k=0; k<5; ++k)
{ thetainside(k, ...);                        // k=0 -> centro; k=1..4 -> borda
  error = summation(...);
  if ((k==0) && (accepang>=1.0e-3))
  { suminten += suminten; bakinten += bakinten; }   // dobra o centro
  else if (accepang<1.0e-3) break;                   // sem abertura: 1 ponto só
}
if (k>0) { suminten /= float(k+1.0); bakinten /= float(k+1.0); }  // k=5 -> /6
```

Em `thetainside` (`mscdrund.cpp:661-664`), `accepnum=0` dá (0,0) e
`accepnum=1..4` dá (accepang, 0°/90°/180°/270°) — os quatro pontos da borda a 90°
um do outro. O manual escreve os azimutes como (θ,0) (θ,θ/2) (θ,θ) (θ,3θ/2), com
2θ valendo a volta completa; é a mesma coisa.

**É por isso que cada ponto da varredura custa 5 resoluções de espalhamento
múltiplo**, e não uma. São 779 × 5 = 3895 no total.

---

## Da amplitude à observável

Fim do `summation()` (`mscdrund.cpp:128-151`) — é a eq. (4) e a eq. (1):

```cpp
cxa += asum[id+j]*csum[j];              // onda espalhada  (Σ φ_sj)
cxb = onemidetec(akin,ia,alf,am,...);   // onda direta     (φ₀)
cxc = matrixelement(ali,alf,am,akin);   // m_lf,c e^(iδ)   eq.(5)
dsum += cxa*cxc;  esum += cxb*cxc;
*suminten += emiter*norm(dsum+esum);    // I   = |φ₀+Σφ_sj|²   eq.(1)
*bakinten += emiter*norm(esum);         // I₀  = |φ₀|²
```

`dsum` é o total, `esum` é só a onda direta. A modulação sai logo depois, em
`mscdrund.cpp:246-249`:

```cpp
netinten = suminten/bakinten - 1.0f;    // χ = I/I₀ − 1
```

com corte em ±10. E o `factors = 0.6724 0.8647` impresso no fim são os dois
fatores-R da eq. (41), calculados em `Pdintensity::reliability`
(`pdintena.cpp:408`).

O `t_l = sin(δ_l)·exp(iδ_l)` da pág. 4 do manual é `Phaseshift::fsinexp`
(`phase.cpp:294`), alimentado pelos arquivos `psl9.txt` / `psAg111-slab.txt`; a
parte radial R_Ekin,lf da eq. (5) vem de `rml9.txt` via a classe `Radialmatrix`.

---

## O que o manual **não** cobre

- **A extensão ATA** (`Mscdrun::ATAevenelem`, `mscdrunc.cpp:45`) não está neste
  manual: é adição posterior, da UNICAMP, com t efetivo
  t = (1−w)·t₁ + w·t₂ para ligas aleatórias. Referência: Soares, de Siervo,
  Landers e Kleiman, *Surface Science* **497** (2002) 205–213 (`MSCD_ATA.pdf`).
  Entra dentro de `alltrievent` (`mscdrunc.cpp:222`) e `alldblevent`
  (`mscdrunc.cpp:779`), ou seja, na eq. (47), trocando o `t_l`.
- **O paralelismo.** O manual só diz que existe (Cray T3E, COMPS). Como o
  trabalho é dividido, e o que custa o quê, está medido no `README.md`.
- **`eq. (40)`** não existe na numeração do manual — pula de (39) para (41).
