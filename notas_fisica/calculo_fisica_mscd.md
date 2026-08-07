# Cálculo Didático da Intensidade (XPD) para 3 Átomos

Para evitar a loucura das funções de Hankel e harmônicos esféricos complexos do MSCD, vamos calcular o caso mais simples possível no papel (Aproximação de Onda Plana).

Temos 3 átomos:
- **Átomo 0 (Emissor):** Na origem $(0,0,0)$.
- **Átomo 1 (Espalhador A):** Na posição $\vec{r}_1$.
- **Átomo 2 (Espalhador B):** Na posição $\vec{r}_2$.

O detector está muito longe, na direção definida pelo vetor unitário $\hat{k}$.

---

## 1. A Onda Direta ($A_0$)
A onda que sai do emissor e vai direto para o detector não sofre nenhum desvio. Podemos representá-classificar a amplitude dela simplesmente como a onda viajando na direção $\hat{k}$:

$$ A_0 = e^{i \vec{k} \cdot \vec{r}_0} $$

Como $\vec{r}_0 = 0$ (ele está na origem), temos:
$$ A_0 = 1 $$

*(Nota: Na vida real, multiplicamos pelo elemento de matriz inicial do fóton, mas vamos ignorar isso para focar na interferência).*

---

## 2. Espalhamento Simples no Átomo 1 ($A_1$)
O elétron viaja até o Átomo 1, bate nele e vai para o detector. 
- Ele viaja uma distância $r_1$, acumulando uma fase $e^{i k r_1}$.
- Ao bater, a força com que o átomo rebate a onda na direção $\hat{k}$ é dada pelo Fator de Espalhamento Complexo $f(\theta_1)$, onde $\theta_1$ é o ângulo de desvio.
- Do átomo 1 até o detector, ele acumula uma fase geométrica $-i \vec{k} \cdot \vec{r}_1$.

A amplitude da onda que bate no Átomo 1 e vai pro detector é:

$$ A_1 = \frac{e^{i k r_1}}{r_1} \cdot f(\theta_1) \cdot e^{-i \vec{k} \cdot \vec{r}_1} $$

E a mesma coisa acontece para o **Átomo 2**, de forma independente:

$$ A_2 = \frac{e^{i k r_2}}{r_2} \cdot f(\theta_2) \cdot e^{-i \vec{k} \cdot \vec{r}_2} $$

---

## 3. A Intensidade Total
Se desligarmos o espalhamento duplo (Átomo 1 rebatendo pro 2), a física quântica diz que as ondas se somam antes de serem detectadas (Princípio da Superposição). 

A Amplitude Total ($A_{total}$) chegando no detector é:

$$ A_{total} = A_0 + A_1 + A_2 $$

Substituindo tudo:

$$ A_{total} = 1 + \left( \frac{e^{i k r_1}}{r_1} f(\theta_1) e^{-i \vec{k} \cdot \vec{r}_1} \right) + \left( \frac{e^{i k r_2}}{r_2} f(\theta_2) e^{-i \vec{k} \cdot \vec{r}_2} \right) $$

Para achar o que você realmente enxerga na tela do computador (a **Intensidade** $I$), você tira o Módulo ao Quadrado da Amplitude:

$$ I(\hat{k}) = | A_{total} |^2 = A_{total} \cdot A_{total}^* $$

Quando você faz essa multiplicação, nascem os termos cruzados (como $A_0 A_1^*$). São esses termos cruzados que dependem da diferença de caminho $(\Delta r = k r_1 - \vec{k} \cdot \vec{r}_1)$ e geram os picos construtivos e vales destrutivos no gráfico de difração!

---

## O Desafio Oculto
Para colocar números nisso na caneta:
1. Você escolheria o vetor do detector $\hat{k} = (\sin\theta\cos\phi, \sin\theta\sin\phi, \cos\theta)$.
2. Calcularia os produtos escalares $\vec{k} \cdot \vec{r}_1$.
3. Teria que ter o valor complexo de $f(\theta_1)$ pronto (algo como $0.5 + 0.2i$).
4. Faria a soma complexa e tiraria o módulo.

Mesmo super simplificado, fazer para 1 único ângulo $\theta, \phi$ ocuparia uma folha de caderno inteira! Agora imagine o código repetindo isso milhões de vezes.
