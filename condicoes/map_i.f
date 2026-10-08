c=======================================================================
c     CONDICOES: MAP_I  (sistema de Netuno, 8 corpos)
c     Grade: semi-eixo x INCLINACAO de Naiad (corpo 1)
c     (vem de map_i.f)
c
c     ATENCAO: no map_i.f original a inclinacao do loop (dii) nunca
c     chegava ao ENTRE, entao a grade nao variava nada na direcao y
c     (e a excentricidade usada era a variavel exc nao inicializada).
c     Aqui a inclinacao de Naiad e' de fato YV (em GRAUS) e e(1)=0,
c     que era o valor efetivo no original. Confira o intervalo Y0..YF!
c=======================================================================

c-----------------------------------------------------------------------
c     PARAMETROS DA GRADE E DA INTEGRACAO
c       N      : numero de corpos (sem contar o planeta)
c       NALVO  : corpo cujos a, e, I sao analisados (max e FFT)
c       NLIN   : numero de saidas por integracao (potencia de 2)
c       TTF    : intervalo entre saidas (dias)
c       X0,XF,NX : semi-eixo inicial, final (km) e numero de pontos
c       Y0,YF,NY : 2a variavel da grade (aqui inclinacao, em graus)
c       LL     : precisao do RA15 (tolerancia 10**(-LL))
c-----------------------------------------------------------------------
      SUBROUTINE CONFIG(N,NALVO,NLIN,TTF,X0,XF,NX,Y0,YF,NY,LL)
      IMPLICIT REAL *8(a-h,o-z)

      N=8
      NALVO=1

c     2**18
      NLIN=262144
c     periodo medio de Naiad / 10
      TTF=0.2943D0/10.0d0

      X0=48290.0d0
      XF=48300.0d0
      NX=101

      Y0=0.0d0
      YF=0.005d0
      NY=101

      LL=12
      return
      end

c-----------------------------------------------------------------------
c     MASSAS E CONDICOES INICIAIS
c     Chamada para cada ponto da grade, com:
c       SEMI = semi-eixo do ponto (km)
c       YV   = 2a variavel do ponto (aqui inclinacao, em graus)
c     Tem que devolver dm0, G, dj2, dj4, as massas (via MASSA),
c     a,e,di,w,om,am (pelo menos do corpo NALVO) e xpla.
c-----------------------------------------------------------------------
      SUBROUTINE INICIAIS(SEMI,YV,N,dm0,G,dj2,dj4,dm,ami,dmimi,
     |                    a,e,di,w,om,am,xpla)
      IMPLICIT REAL *8(a-h,o-z)
      include 'dimensoes.inc'
      DIMENSION dm(NMAX),ami(NMAX),dmimi(NMAX),a(NMAX),e(NMAX),
     |di(NMAX),w(NMAX),om(NMAX),am(NMAX),xpla(NMAX,6)

      pi=4.0d0*datan(1.0d0)
      conv=pi/180.0d0

c     Planeta: massa (kg) e raio equatorial (km)
      dm0=1.0d0
      dmplaneta=102.4092d24
      requat=24764.0d0
      CALL GRAVIT(dmplaneta,requat,G)

c     J2, J4  (French et al. 2024)
      dj2= 3536.3d-6
      dj4=-36.0d-6

c Naiad  (corpo da grade)
      CALL MASSA(1,1.2d17/1.024092d26,dm0,G,dm,ami,dmimi)
      a(1)=SEMI/requat
      e(1)=0.0d0
      di(1)=YV*conv
      w(1)=3.451702741265048d+02*conv
      om(1)=2.257981437775343d+01*conv
      am(1)=1.377006934371379d+00*conv

c Thalassa
      CALL MASSA(2,3.54d17/1.024092d26,dm0,G,dm,ami,dmimi)
      a(2)=5.013968111887394d+04/requat
      e(2)=1.484152758507258d-03
      di(2)=5.894669689234114d-01*conv
      om(2)= 1.181264843636133d+01*conv
      w(2)=3.534651544354857d+02*conv
      am(2)=3.574559972800677d+02*conv

c Despina
      CALL MASSA(3,2.1d18/1.024092d26,dm0,G,dm,ami,dmimi)
      a(3)=5.258820669225104d+04/requat
      e(3)=1.550192775154805d-03
      di(3)=4.602061682258348d-01*conv
      w(3)=2.142808120389758d+02*conv
      om(3)=3.552217165239591d+02*conv
      am(3)=3.584333608541766d+02*conv

c Galatea
      CALL MASSA(4,1.94d18/1.024092d26,dm0,G,dm,ami,dmimi)
      a(4)=6.200528852643136d+04/requat
      e(4)=1.046394401605568d-03
      di(4)=5.085670339398293d-01*conv
      w(4)=2.250808351391018d+02*conv
      om(4)=3.586238854863935d+02*conv
      am(4)=3.586475239710762d+02*conv

c Larissa
      CALL MASSA(5,4.2d18/1.024092d26,dm0,G,dm,ami,dmimi)
      a(5)=7.359220611324030d+04/requat
      e(5)=7.583741858377249d-04
      di(5)=7.582271536412154d-01*conv
      w(5)=1.822601917448936d+02*conv
      om(5)=3.552768352859704d+02 *conv
      am(5)=1.297912784496201d+02*conv

c Hippocamp
      CALL MASSA(6,2.097d16/1.024092d26,dm0,G,dm,ami,dmimi)
      a(6)=1.053143448007733d+05/requat
      e(6)=5.804868704543895d-04
      di(6)=7.795096858029827d-01*conv
      w(6)=3.059981393949778d+02*conv
      om(6)=3.593479472191129d+02*conv
      am(6)=6.37751942368533d+01*conv

c Proteus
      CALL MASSA(7,4.4d19/1.024092d26,dm0,G,dm,ami,dmimi)
      a(7)=1.176756915329071d+08/requat
      e(7)=3.491769246372516d-04
      di(7)=1.028618541299257d+00*conv
      w(7)=1.232604547938490d+02*conv
      om(7)=3.533625940682254d+02*conv
      am(7)=1.075988743969018d+02*conv

c Triton
      CALL MASSA(8,2.14d22/1.024092d26,dm0,G,dm,ami,dmimi)
      a(8)=3.547658968420361d+08/requat
      e(8)=2.520758828858335d-05
      w(8)=2.831380537336281d+02 *conv
      om(8)=1.724421099342844d+02*conv
      am(8)=3.449417938724240d+01*conv
      di(8)=1.568252139565797d+02*conv

c     elementos -> cartesianas (xpla) para todos os corpos
      CALL ELEM2XYZ(N,a,e,di,w,om,am,ami,xpla)

      return
      end
