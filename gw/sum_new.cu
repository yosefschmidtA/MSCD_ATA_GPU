__constant__ int c_lamda[64];

/* ---------------- Fase 3 refeita (08/10/2026): summation ----------------
   O kernel antigo usava UMA thread por par (ia,ib) e cada thread varria os
   natoms valores de ic em serie. Com 367 atomos sobravam ~1500 pares por
   passo de m -- uma dezena de blocos numa placa que segura ~36 mil threads --
   e 97% dos ic eram descartados dentro do laco. Medido: 345 s de 425 s.

   O desenho novo tem tres pecas, e NENHUMA muda a ordem de uma soma:

   1. Thread por (par, j). csum[j] so' depende da linha j de algam e de
      tevenelem, entao as radim (15) linhas sao independentes. Cada thread
      percorre os trios na MESMA ordem de ic e soma k na MESMA ordem que a
      CPU: csum[j] += (algam*bsum)*tevenelem. Bit a bit igual.

   2. Lista de trios compactada uma vez no setup. evedim, eegdim, mevadd,
      talpha e tgamma nao mudam na corrida (energia fixa, allrotation roda
      uma vez), entao o filtro "ic==ib ou evedim<1" sai do laco.

   3. algam em forma fechada. A recorrencia da CPU (mscdrund.cpp:160) so'
      aplica conj, que e' exato, a uma entrada vizinha cujo (p,q) e' o
      oposto. Pela tabela lamda, toda entrada com p<0, ou p==0 e q<0, e'
      conj da entrada (-p,-q), e essa e' sempre calculada direto por fexpix.
      Entao cada linha sai sozinha, com a mesma expressao xc=-p*xb-q*xa que
      a CPU avaliou -- mesmos bits.

   E a copia bsum<-asum (16 MB por passo de m) sumiu: o passo le asum, que
   so' muda depois, quando k_sum_scatter espalha o resultado. */

typedef struct { int ic, evedim, eegdim, mevadd; float xa, xb; } Gtrio;

static Gcplx *d_tevenelem = NULL;
static Gcplx *d_asum = NULL;
static Gcplx *d_sumtmp = NULL;   /* radim por par sobrevivente */
static short2 *d_pairs = NULL;
static int *d_toff = NULL;        /* inicio dos trios de cada par, +1 */
static Gtrio *d_trios = NULL;
static Gcplx *d_emit = NULL;      /* linhas dos emissores, contiguas */
static int *d_emitidx = NULL;
static Gcplx *h_emit = NULL;      /* pinned */
static int h_emitidx[1250];
static int g_nemit = 0;
static int g_ntrielem = 0;
static int h_pair_offset[16] = {0};
static int h_pair_count[16] = {0};
static int g_maxpairs = 0;
static long g_ntrios = 0;

extern "C" int mscdgpu_setup_summation(
    const int *tevencut, const int *tevendim, const int *tevenadd,
    const float *tevenpar, const float *talpha, const float *tgamma,
    int ntrieven, int ntrielem, const float *patom, int msorder)
{
    int n = K.natoms;
    g_ntrielem = ntrielem;
    if (K.radim > 16) { snprintf(g_err,sizeof(g_err),"radim %d > 16",K.radim); return 1; }
    if (msorder > 15) { snprintf(g_err,sizeof(g_err),"msorder %d > 15",msorder); return 1; }

    CK(cudaMalloc((void**)&d_tevenelem, (size_t)ntrielem * sizeof(Gcplx)));
    CK(cudaMalloc((void**)&d_asum, (size_t)n * n * K.radim * sizeof(Gcplx)));

    int lamda[64];
    for (int j=0; j<32; ++j) {
        int k, m;
        if ((j==0)||(j==3)||(j==10)) k=0;
        else if (j<3) k=3-j*2;
        else if (j<6) k=18-j*4;
        else if (j<8) k=13-j*2;
        else if (j<10) k=51-j*6;
        else if (j<13) k=46-j*4;
        else if (j<15) k=108-j*8;
        else k=0;
        if (j==10) m=2;
        else if ((j==3)||(j==6)||(j==7)||(j==11)||(j==12)) m=1;
        else m=0;
        lamda[j] = k; lamda[32+j] = m;
    }
    CK(cudaMemcpyToSymbol(c_lamda, lamda, 64 * sizeof(int)));

    /* Pares sobreviventes por m (mesma ordem da CPU) e, para cada par, os
       trios com evedim>=1 na ordem crescente de ic. */
    size_t cap = (size_t)(msorder + 1) * n * n;
    short2 *hp = (short2*)malloc(cap * sizeof(short2));
    int *ho = (int*)malloc((cap + 1) * sizeof(int));
    size_t tcap = 1 << 20, nt = 0;
    Gtrio *ht = (Gtrio*)malloc(tcap * sizeof(Gtrio));
    if (!hp || !ho || !ht) { snprintf(g_err,sizeof(g_err),"malloc setup_summation"); return 1; }

    int tot = 0;
    g_maxpairs = 0;
    for (int m = 0; m < 16; ++m) h_pair_count[m] = h_pair_offset[m] = 0;
    for (int m = 2; m <= msorder; ++m) {
        h_pair_offset[m] = tot;
        int count = 0;
        for (int ia = 0; ia < n; ++ia) {
            if (m == 2 && patom[ia * 12 + 7] == 0.0f) continue;
            for (int ib = 0; ib < n; ++ib) {
                if (ib == ia) continue;
                if (tevencut[(m-1)*n*n + ia*n + ib] == 0) continue;
                short2 p; p.x = ia; p.y = ib;
                hp[tot + count] = p;
                ho[tot + count] = (int)nt;
                for (int ic = 0; ic < n; ++ic) {
                    size_t id = (size_t)ia*n*n + (size_t)ib*n + ic;
                    int evedim = tevendim[id];
                    if (m > 8) evedim >>= 24;          /* sizeof(int)==4 */
                    else evedim >>= (m - 2) * 4;
                    evedim &= 15;
                    if (ic == ib || evedim < 1) continue;
                    int k = tevenadd[id];
                    Gtrio t;
                    t.ic = ic; t.evedim = evedim;
                    t.eegdim = (int)tevenpar[k * 10 + 5];
                    t.mevadd = (int)tevenpar[k * 10 + 6];
                    t.xa = talpha ? talpha[id] : 0.0f;
                    t.xb = tgamma ? tgamma[id] : 0.0f;
                    if (evedim > 1 && (t.eegdim >= 16 || !talpha || !tgamma)) {
                        snprintf(g_err,sizeof(g_err),"erro 901 no trio %d %d %d",ia,ib,ic);
                        return 901;
                    }
                    if (nt == tcap) {
                        tcap *= 2;
                        ht = (Gtrio*)realloc(ht, tcap * sizeof(Gtrio));
                        if (!ht) { snprintf(g_err,sizeof(g_err),"realloc trios"); return 1; }
                    }
                    ht[nt++] = t;
                }
                count++;
            }
        }
        h_pair_count[m] = count;
        if (count > g_maxpairs) g_maxpairs = count;
        tot += count;
    }
    ho[tot] = (int)nt;
    g_ntrios = (long)nt;

    if (tot > 0) {
        CK(cudaMalloc((void**)&d_pairs, (size_t)tot * sizeof(short2)));
        CK(cudaMemcpy(d_pairs, hp, (size_t)tot * sizeof(short2), cudaMemcpyHostToDevice));
        CK(cudaMalloc((void**)&d_toff, (size_t)(tot + 1) * sizeof(int)));
        CK(cudaMemcpy(d_toff, ho, (size_t)(tot + 1) * sizeof(int), cudaMemcpyHostToDevice));
        CK(cudaMalloc((void**)&d_sumtmp, (size_t)g_maxpairs * K.radim * sizeof(Gcplx)));
    }
    if (nt > 0) {
        CK(cudaMalloc((void**)&d_trios, nt * sizeof(Gtrio)));
        CK(cudaMemcpy(d_trios, ht, nt * sizeof(Gtrio), cudaMemcpyHostToDevice));
    }
    free(hp); free(ho); free(ht);

    g_nemit = 0;
    for (int ia = 0; ia < n; ++ia)
        if (patom[ia * 12 + 7] != 0.0f) h_emitidx[g_nemit++] = ia;
    if (g_nemit > 0) {
        size_t sz = (size_t)g_nemit * n * K.radim * sizeof(Gcplx);
        CK(cudaMalloc((void**)&d_emit, sz));
        CK(cudaMallocHost((void**)&h_emit, sz));
        CK(cudaMalloc((void**)&d_emitidx, g_nemit * sizeof(int)));
        CK(cudaMemcpy(d_emitidx, h_emitidx, g_nemit * sizeof(int), cudaMemcpyHostToDevice));
    }
    fprintf(stderr, "GPU summation: %d pares, %ld trios, %d emissores, %.1f MB de trios\n",
        tot, g_ntrios, g_nemit, nt * sizeof(Gtrio) / 1048576.0);
    return 0;
}

__global__ static void k_init_asum(Gcplx *asum, const Gcplx *devendetec, int natoms, int radim, int msorder)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int ic = blockIdx.y * blockDim.y + threadIdx.y;
    int ib = blockIdx.z * blockDim.z + threadIdx.z;
    if (j >= radim || ic >= natoms || ib >= natoms) return;

    int id = ib * natoms * radim + ic * radim + j;
    if (msorder > 0 && ic != ib) {
        asum[id] = devendetec[id];
    } else {
        asum[id].re = 0.0f;
        asum[id].im = 0.0f;
    }
}

/* blockDim = (16, PPB): threadIdx.x e' j, threadIdx.y escolhe o par. */
__global__ static void k_sum_step(
    int count, const short2 *pairs, const int *toff, const Gtrio *trios,
    const Gcplx *asum, const Gcplx *devendetec, const Gcplx *tevenelem,
    const Gcplx *cexpix, int natoms, int radim, int exndata, int exmdata,
    Gcplx *out)
{
    int j = threadIdx.x;
    int pidx = blockIdx.x * blockDim.y + threadIdx.y;
    if (pidx >= count || j >= radim) return;

    int ia = pairs[pidx].x;
    int ib = pairs[pidx].y;
    Gcplx c = devendetec[ia * natoms * radim + ib * radim + j];
    int p = c_lamda[j];
    int t1 = toff[pidx + 1];

    for (int t = toff[pidx]; t < t1; ++t) {
        Gtrio T = trios[t];
        if (j >= T.evedim) continue;
        int megadd = ib * natoms * radim + T.ic * radim;
        if (T.evedim == 1) {
            c = cadd(c, cmul(asum[megadd], tevenelem[T.mevadd]));
        } else {
            float xa = T.xa, xb = T.xb;
            const Gcplx *te = tevenelem + T.mevadd + j * T.eegdim;
            for (int k = 0; k < T.evedim; ++k) {
                int q = c_lamda[k];
                Gcplx al;
                if (p == 0 && q == 0) {
                    al.re = 1.0f; al.im = 0.0f;
                } else if (p < 0 || (p == 0 && q < 0)) {
                    int p2 = -p, q2 = -q;
                    float xc = -p2 * xb - q2 * xa;
                    Gcplx e = d_fexpix(cexpix, exndata, exmdata, xc);
                    al.re = e.re; al.im = -e.im;
                } else {
                    float xc = -p * xb - q * xa;
                    al = d_fexpix(cexpix, exndata, exmdata, xc);
                }
                c = cadd(c, cmul(cmul(al, asum[megadd + k]), te[k]));
            }
        }
    }
    out[pidx * radim + j] = c;
}

__global__ static void k_sum_scatter(int count, const short2 *pairs,
    const Gcplx *tmp, Gcplx *asum, int natoms, int radim)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= count * radim) return;
    int pidx = i / radim, j = i - pidx * radim;
    asum[pairs[pidx].x * natoms * radim + pairs[pidx].y * radim + j] = tmp[i];
}

__global__ static void k_gather_emit(int nemit, const int *emitidx,
    const Gcplx *asum, Gcplx *out, int rowlen)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nemit * rowlen) return;
    int e = i / rowlen, r = i - e * rowlen;
    out[i] = asum[(size_t)emitidx[e] * rowlen + r];
}

#define SUM_PPB 8
static float last_akin = -1.0f;
extern "C" int mscdgpu_summation(float akin, const Gcplx *tevenelem, Gcplx *asum_host, const float *patom)
{
    if (!g_ready) { snprintf(g_err,sizeof(g_err),"setup nao chamado"); return 1; }

    /* So' reenvia tevenelem quando a energia muda (scanmode=223, k fixo). */
    if (akin != last_akin) {
        CK(cudaMemcpy(d_tevenelem, tevenelem, (size_t)g_ntrielem * sizeof(Gcplx), cudaMemcpyHostToDevice));
        last_akin = akin;
    }

    int n = K.natoms, rd = K.radim;
    dim3 threads_init(16, 8, 8);
    dim3 blocks_init((rd + 15)/16, (n + 7)/8, (n + 7)/8);
    k_init_asum<<<blocks_init, threads_init>>>(d_asum, D.devendetec, n, rd, K.msorder);
    CK(cudaGetLastError());

    for (int m = K.msorder; m >= 2; --m) {
        int count = h_pair_count[m];
        if (count <= 0) continue;
        int off = h_pair_offset[m];
        dim3 tb(16, SUM_PPB);
        k_sum_step<<<(count + SUM_PPB - 1) / SUM_PPB, tb>>>(
            count, d_pairs + off, d_toff + off, d_trios,
            d_asum, D.devendetec, d_tevenelem, D.cexpix, n, rd,
            K.exndata, K.exmdata, d_sumtmp);
        CK(cudaGetLastError());
        int tot = count * rd;
        k_sum_scatter<<<(tot + 255) / 256, 256>>>(count, d_pairs + off, d_sumtmp, d_asum, n, rd);
        CK(cudaGetLastError());
    }

    if (g_nemit > 0) {
        int rowlen = n * rd, tot = g_nemit * rowlen;
        k_gather_emit<<<(tot + 255) / 256, 256>>>(g_nemit, d_emitidx, d_asum, d_emit, rowlen);
        CK(cudaGetLastError());
        CK(cudaMemcpy(h_emit, d_emit, (size_t)tot * sizeof(Gcplx), cudaMemcpyDeviceToHost));
        for (int e = 0; e < g_nemit; ++e)
            memcpy(asum_host + (size_t)h_emitidx[e] * rowlen, h_emit + (size_t)e * rowlen,
                (size_t)rowlen * sizeof(Gcplx));
    }
    return 0;
}

extern "C" void mscdgpu_teardown(void)
{ if (!g_ready) return;
  cudaFree(D.patom); cudaFree(D.devenpar); cudaFree(D.rotmata);
  cudaFree(D.rotmatc); cudaFree(D.thermat); cudaFree(D.aweight);
  cudaFree(D.hankmat_a); cudaFree(D.hankarg_b); cudaFree(D.phasec);
  cudaFree(D.cexpix); cudaFree(D.pairgeo); cudaFree(D.pairkind);
  cudaFree(D.alnum);
  cudaFree((void*)D.devenadd); cudaFree(D.devendetec);

  cudaFree(d_tevenelem); cudaFree(d_asum); cudaFree(d_sumtmp);
  cudaFree(d_pairs); cudaFree(d_toff); cudaFree(d_trios);
  cudaFree(d_emit); cudaFree(d_emitidx); cudaFreeHost(h_emit);
  d_tevenelem=d_asum=d_sumtmp=d_emit=NULL; d_pairs=NULL; d_toff=NULL;
  d_trios=NULL; d_emitidx=NULL; h_emit=NULL;

  last_akin = -1.0f;
  memset(&D,0,sizeof(D)); g_ready=0;
}
