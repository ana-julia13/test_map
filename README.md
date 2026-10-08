# Mapas dinâmicos (N corpos + J2/J4 + FFT)

O código foi dividido em duas partes:

| Arquivo | O que é | Precisa mexer? |
|---|---|---|
| `mapa_base.f` | integrador RA15, forças, J2/J4, conversões, FFT, saídas | **não** |
| `dimensoes.inc` | número máximo de corpos (`NMAX=20`) | só se tiver mais de 20 corpos |
| `<caso>.f` (`naiad.f`, `victor.f`, ...) | grade, passo, planeta, massas e condições iniciais | **sim, é aqui que se trabalha** |
| `originais/` | os 4 programas antigos, só para consulta | não |

## Como compilar e rodar

```sh
gfortran -o mapa_naiad naiad.f
./mapa_naiad
```

Cada arquivo de caso termina com `include 'mapa_base.f'`, então ele já puxa
o código base sozinho (os três arquivos têm que estar na mesma pasta).
Use `-O2` para rodar mais rápido. Ou use o Makefile: `make run CASO=naiad`.

Casos que já existem (convertidos dos programas antigos):

| Caso | Veio de | Corpos | Grade |
|---|---|---|---|
| `naiad` | `mapa_naiad.f` | 8 (Netuno) | a × e de Naiad (corpo 1) |
| `map_i` | `map_i.f` | 8 (Netuno) | a × **I** de Naiad (corpo 1) |
| `victor` | `mapa_victor.f` | 9 | a × e de Perdita (corpo 9) |
| `mapa` | `mapa.f` | 7, em posição/velocidade | a × e do corpo 2 |

## Como criar um caso novo

Copie um arquivo parecido (`cp naiad.f meucaso.f`), edite
e compile com `gfortran -o mapa_meucaso meucaso.f`. O arquivo de condições tem só duas subrotinas:

**`CONFIG`**: os parâmetros da rodada:

```fortran
      N=8                 ! numero de corpos
      NALVO=1             ! corpo analisado (emax, amax, imax, FFT)
      NLIN=262144         ! numero de saidas (potencia de 2)
      TTF=0.2943D0/10.0d0 ! intervalo entre saidas (dias)
      X0=48290.0d0        ! semi-eixo inicial da grade (km)
      XF=48300.0d0        ! semi-eixo final (km)
      NX=101              ! numero de pontos em a
      Y0=0.0d0            ! 2a variavel: inicial
      YF=0.005d0          !              final
      NY=101              !              numero de pontos
      LL=12               ! precisao do RA15
```

**`INICIAIS(SEMI,YV,...)`**: é chamada a cada ponto da grade com `SEMI` (km)
e `YV` (a 2ª variável). Ela define o planeta, as massas e as condições
iniciais. Quem decide onde `YV` entra é você: `e(1)=YV` para um mapa a × e,
`di(1)=YV*conv` para um mapa a × I etc.

Rotinas de ajuda (já estão no código base):

| Rotina | Faz |
|---|---|
| `GRAVIT(dmplaneta,requat,G)` | G nas unidades do programa (raios equatoriais, massas do planeta, dias) |
| `MASSA(i,dmi,dm0,G,dm,ami,dmimi)` | massa do corpo i (em massas do planeta), já calcula `ami` e `dmimi` |
| `ELEM2XYZ(N,a,e,di,w,om,am,ami,xpla)` | elementos → cartesianas, todos os corpos |
| `ELEM2XYZ1(i,...)` | elementos → cartesianas, só o corpo i |
| `POSVEL(i,x,y,z,vx,vy,vz,xpla)` | dá o corpo i em posição/velocidade |
| `XYZ2ELEM1(i,xpla,ami,a,e,di,w,om,am)` | cartesianas → elementos do corpo i |

No fim, `INICIAIS` tem que ter preenchido `xpla` de todos os corpos e
`a`, `e`, `di` do corpo `NALVO`.

## Arquivos de saída

São os mesmos nomes de antes (`emax.dat`, `emax.grd`, `amax.*`, `imax.*`,
`*-t.dat`, `Amedio.*`, `Emedio.*`, `Imedio.*`, `EcosVmedio.*`,
`dn<lim><var>.grd`, `e1e2dn<lim><var>.dat`). Os cabeçalhos DSAA dos `.grd`
(número de pontos e intervalos) agora são gerados a partir do `CONFIG`, e
não precisam mais ser trocados à mão.

## O que mudou em relação aos programas antigos

Comparei com os originais na mesma grade pequena: `naiad` e `victor` dão
**exatamente** os mesmos `emax`, `imax`, valores médios e números espectrais.
O que foi corrigido de propósito:

1. **`t` não era inicializado**: compilado sem `-finit-local-zero` o original
   pode travar. Agora `t`, `TTT`, `TM`, `XL` começam em zero.
2. **`amax` nunca mudava**: comparava `a` (em raios do planeta) com `semi`
   (em km), então `amax.dat` saía sempre igual a `semi`. Agora começa de
   `a(NALVO)` e fica em raios do planeta, como o `Amedio`.
3. **Cabeçalhos `.grd` errados** (ex.: o victor dizia `51 71` e
   `198010.514 198209.3984`). Agora são automáticos.
4. **`e1e2dn10i.dat` era aberto duas vezes** (inclinação e e·cos ϖ), o que
   dá erro no gfortran. O de e·cos ϖ virou `e1e2dn10evarpi.dat`. Também os
   unidades 4490/44901 eram reusadas, e `e1e2dn1i.dat`/`e1e2dn3i.dat`
   ficavam vazios. Agora cada arquivo tem a sua unidade.
5. **e·cos ϖ**: o vetor da 4ª FFT nunca era preenchido (FFT de zeros). Agora
   é `e*cos(w+om)` do corpo `NALVO`, como dizia a linha comentada.
6. **`AmedioEmedioEVmedioImedio.DAT`**: uma linha por ponto da grade (antes
   eram 4 linhas, com valores lixo).
7. **naiad**: `ami(7)`, `dmimi(7)` usavam `dm(6)` e `ami(8)`, `dmimi(8)`
   usavam `dm(7)`. Corrigido (no `map_i.f` já estava certo).
8. **map_i**: a inclinação do loop (`dii`) nunca chegava às condições
   iniciais, então o mapa não variava em y. Agora `di(1)=YV*conv` (YV em
   graus) e `e(1)=0`, que era o valor que o original usava na prática.
9. **mapa.f**: `fitzpaorb` usava `pi`, `pi2`, `pi05`, `pi15` de commons que
   nunca eram preenchidos (valiam 0). Agora usa `COMMON/IN`.

## Coisas para conferir nos dados (não mudei)

- **mapa.f**: o original (e o novo) param logo no início com
  `exc > =1, parar`: a órbita do corpo 2 sai hiperbólica com esses números.
  A posição/velocidade do corpo 1 parece estar em UA e UA/dia, e não em km.
  O último número tem `6.043612890766993-4` sem o `d`. A massa do corpo 1
  era 225.8d15 kg no programa principal e 3.75d22 kg na `FORCE`; agora vale
  225.8d15 para os dois. O corpo 7 tem massa maior que a do planeta.
- **naiad / map_i**: `a(7)` (Proteus) e `a(8)` (Triton) estão em **metros**
  divididos por `requat` em km, então ficam 1000× mais longe.
- **gaussj** usa `REAL big,dum,pivinv` (precisão simples) ao calcular os
  momenta iniciais, por isso as condições iniciais têm só ~7 dígitos.
  Mantive assim para reproduzir os resultados antigos.
- O filtro da FFT (só conta picos com período > 5 dias) está em
  `PARAMETER (PERMIN=5.0d0)` na rotina `FFT` do código base.
