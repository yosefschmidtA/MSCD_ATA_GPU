# MSCDAUTOMATOR: O que é e como funciona

O **MSCDAUTOMATOR** (localizado em `/home/yosef/MSCD_AUTOMATOR`) é um pipeline em Python conteinerizado que automatiza a geração dos arquivos de Phase Shifts (`ps*.txt`) e Matrizes Radiais (`rm*.txt`) exigidos pelo MSCD, orquestrando os antigos códigos em Fortran do Van Hove.

## O Fluxo de Execução

A automação é disparada pelo script `run_half.sh` (ou `run_all.sh`), que compila os binários necessários (`phsh0`, `poconv`, `randmscd_parallel`) e depois executa a seguinte cadeia de scripts Python dentro da pasta `arquivos/`:

1. **`gerador.py`**:
   - Lê o arquivo `input_cluster.txt` (que contém um cabeçalho expandido com as configurações atômicas e a energia do orbital).
   - Extrai o parâmetro de rede e converte a geometria para coordenadas que o Fortran entende.
   - Gera o arquivo `cluster.i` (geometria) e concatena os arquivos da base de orbitais (`atelem.*`) para formar o `atomic.i`.

2. **`muff.py`**:
   - **Cálculo Bulk**: Roda o executável `phsh1` em modo Bulk (`opção 0`) para calcular o *Muffin-tin zero* (MTZ).
   - **Cálculo Slab**: Com o MTZ calculado, roda o `phsh1` em modo Slab (`opção 1`) para gerar o arquivo contendo os potenciais Muffin-tin de todos os elementos (`mufftin.d`).
   - Corta o `mufftin.d` em arquivos individuais (`mufftin1.d`, `mufftin2.d`, etc.) para cada elemento.
   - Roda o programa em C++ **`poconv`**, que converte esses potenciais no formato intermediário `ps*.txt`.
   - Gera dinamicamente os arquivos de entrada (`psrmin*.txt`) e roda o programa **`psrm.x`**, que de fato resolve a equação de Schrödinger radial, produzindo os **Phase Shifts finais** (`ps*.1.txt`) e as **Matrizes Radiais** (`rm*.txt`).

3. **`leitoF.py` / `criador_final.py`** e ferramentas de plotagem (`teo.py`, `exp.py`):
   - Organizam as saídas, movem arquivos para os lugares corretos e processam os resultados das simulações numéricas para exibir gráficos no frontend Flask (`app.py`).

## Quando usar?
Este pipeline é extremamente útil quando precisarmos simular um novo material (com um novo elemento ou nova rede cristalina) e não tivermos os phase shifts prontos. Ele abstrai toda a complexidade de rodar o `phsh1` e o `psrm.x` manualmente.
