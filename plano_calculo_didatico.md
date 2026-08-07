# Plano de Voo: Cálculo de Intensidade (Rehr-Albers na Unha)

Este é o nosso roteiro para provar, na ponta do lápis, como o código do MSCD usa o formalismo de Rehr-Albers para calcular a difração de fotoelétrons, garantindo que nossos números batam exatamente com o gabarito do computador.

## Objetivo
Calcular a intensidade de emissão para **qualquer ângulo escolhido** (por exemplo $\theta = 30^\circ, \phi = 90^\circ$ ou $\theta = 72^\circ, \phi = 120^\circ$) em um aglomerado bidimensional perfeito de **4 átomos**, usando *Single Scattering*.

## O "Hack" (A Sacada Matemática)
Para fugirmos da multiplicação insana de matrizes $16 \times 16$, usaremos um **Átomo Didático**.
Se o espalhamento for restrito à onda `s` (momento angular $l=0$), o formalismo de Rehr-Albers colapsa. Todas as matrizes perdem suas dimensões e viram matrizes $1 \times 1$ (ou seja, números escalares comuns). A física contínua válida, mas a álgebra cabe numa folha A4.

---

## Os 5 Passos do Experimento

### Passo 1: Preparar a Física do Átomo (Fase `l=0`) [CONCLUÍDO]
- **Ação:** Criamos o arquivo de phase shift (`didatico_co.ph`) para o Cobalto, mantendo apenas a coluna de momento angular $l=0$ não-nula, e zerando todas as outras ($l \ge 1$).
- **Consequência:** A T-Matrix do espalhamento virará um único número complexo $t_0 = \sin(\delta_0)e^{i\delta_0}$.

### Passo 2: Construir o Cluster Quadrado (4 Átomos) [CONCLUÍDO]
- **Ação:** Criamos o input (`4atomos_didatico.in`) definindo exatamente 4 átomos no plano XY: $(0,0,0)$, $(4.086,0,0)$, $(0,4.086,0)$ e $(4.086,4.086,0)$. O segredo para o MSCD não criar uma esfera foi setar `radius = 0.0` e usar `latoms(0, 1, 0, 1)`.
- **Configuração no MSCD:** 
  - `msorder = 1` e `raorder = 1` (Apenas Espalhamento Simples - Single Scattering).
  - Emissor fixo na origem $(0,0,0)$.
  - Modo de varredura `scanmode = 223` (Energy fixed, Angle scan do arquivo experimental).
  - **Idealização:** Zeramos `tdebye`, `tsample` e `vinner` para remover fatores térmicos (Debye-Waller) e de refração, criando um "vácuo ideal" perfeito para o cálculo de papel.

### Passo 3: O Gabarito do MSCD [CONCLUÍDO]
- **Ação:** A simulação foi executada com sucesso.
- **Resultado:** O MSCD cuspiu as intensidades para uma varredura completa de ângulos no arquivo `saida_didatico_co.txt`.
- **Objetivo Atual:** Escolher qualquer linha (par de $\theta$ e $\phi$) dessa saída para se tornar a nossa "verdade absoluta" que teremos que alcançar no papel.

### Passo 4: A Matemática no Papel [PRÓXIMO PASSO]
Com base na energia extraída do arquivo de fase ($k = 13.63$):
1. Calcular o vetor do detector $\hat{k}$ para o ângulo escolhido.
2. Calcular a fase viajante ($kr$) para os 3 caminhos de espalhamento possíveis.
3. Extrair a T-Matrix $t_0$ do arquivo `didatico_co.ph` para o momento $k=13.63$.
4. Calcular a defasagem geométrica projetada no detector ($-\vec{k} \cdot \vec{r}$).
5. Somar as 4 amplitudes (Onda Direta + Espalhada 1 + Espalhada 2 + Espalhada 3).
6. Elevar a norma complexa ao quadrado para achar a Intensidade $I_{papel}$.

### Passo 5: O Check de Sanidade
- Comparamos $I_{papel}$ com a intensidade lida na linha correspondente de `saida_didatico_co.txt`. 
- Se os números baterem na décima casa decimal, a física e o formalismo de Rehr-Albers (na sua versão escalar mais simples) estarão 100% dominados!
