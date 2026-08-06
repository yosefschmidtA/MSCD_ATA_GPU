import sys

with open("CLAUDE.md", "r") as f:
    text = f.read()

rep_old1 = """**Fases 0 e 1 estão FEITAS e
  validadas (05/08/2026); a próxima é a Fase 2.**"""

rep_new1 = """**Fases 0 a 3 estão FEITAS e
  validadas (05/08/2026); a próxima é a Fase 4.**"""

if rep_old1 in text:
    text = text.replace(rep_old1, rep_new1)

rep_old2 = """## Estado atual (05/08/2026)

**O port para CUDA começou.** Fase 1 (`alldblevent`) validada com o piso de ruído
1,0e-5 medido no V5. O executável híbrido roda, mas por enquanto transfere a
matriz `devenelem` da placa para o host a cada ponto — o que derrubou a vantagem
da placa (`np=1` de 122 s no CPU para 59 s na GPU, quando o teto de Amdahl
dava margem para ~20 s).

- **Fase 2** (`allevendetec`): (Pendente). É a próxima ação. A meta é criar e
  alimentar `devendetec` na placa, e ainda devolver para a CPU (o tempo pode
  piorar transitoriamente).
- **Fase 3** (Summation para GPU): (Pendente). Essa é a chave para acabar com as
  cópias gigantes de array entre placa e CPU.
- **Fase 4** (pathcut na GPU): (Opcional, 1.7s)."""

rep_new2 = """## Estado atual (05/08/2026, noite)

**A Fase 3 do port de CUDA está feita e validada.** Contudo, antes na Fase 2 o tempo aumentou para **66,19 s** em `np=1`. Isso ocorreu porque a transferência PCIe aumentou levemente (copiávamos os 7,3 MB do `devendetec` de volta para a CPU). Na Fase 3:

- **Fase 3** (Summation para GPU): **CONCLUÍDO.** (05/08/2026). Tráfego pesado contornado! A GPU agora resolve todo o loop `m` com matriz esparsa e devolve apenas as linhas dos átomos emissores (31 KB/ponto). Tempo despencou de 66.19s (Fase 2) para fenomenais **37.69s**. Gargalo principal do PCIe aniquilado.
- **Fase 4** (onevenemit / onemidetec na GPU): (Pendente).
- **Fase 5** (pathcut na GPU): (Opcional, 1.7s).

## Histórico de Armadilhas

**05/08/2026 - O "Falso Positivo" ZERADO na Fase 3**
- **O erro:** Durante a implementação da Fase 3 (`summation`), a versão CPU do loop foi completamente removida do bloco compilado com `-DMSCDGPU`, em vez de ter seu desvio em tempo de execução via `getenv("MSCD_GPU")`. 
- **O sintoma:** O script `./baseline/regressao-gpu.sh 1` passava com precisão perfeita, mas isso porque ele forçava a execução na GPU (`MSCD_GPU=1`). Quando o usuário rodou o binário manualmente *sem a variável de ambiente* (como exigido no fallback de CPU detalhado no `PLANO_CUDA.md`), a execução falhou silenciosamente e cuspiu zeros em todo o arquivo `.chi`.
- **A correção:** O usuário percebeu que o arquivo de saída gerado estava preenchido com zeros. A solução foi restaurar a rota original da CPU e inserir um bloco `if (getenv("MSCD_GPU"))` dentro do macro, garantindo que o binário suporte as duas vias em tempo de execução. Nunca remova o caminho CPU da função, ele deve coexistir!"""

if rep_old2 in text:
    text = text.replace(rep_old2, rep_new2)
    with open("CLAUDE.md", "w") as f:
        f.write(text)
    print("Patched!")
else:
    print("rep_old2 not found!")
