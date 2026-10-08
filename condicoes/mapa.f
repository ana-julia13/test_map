c=======================================================================
c     CONDICOES: MAPA  (7 corpos dados em posicao/velocidade)
c     Grade: semi-eixo x excentricidade do corpo 2
c     (vem de mapa.f)
c
c     Aqui os corpos sao dados em coordenadas cartesianas (km, km/dia?).
c     O corpo da grade (2) e' convertido para elementos, troca-se
c     a e e pelos valores da grade, e volta-se para cartesianas.
c=======================================================================

c-----------------------------------------------------------------------
c     PARAMETROS DA GRADE E DA INTEGRACAO  (ver condicoes/naiad.f)
c-----------------------------------------------------------------------
      SUBROUTINE CONFIG(N,NALVO,NLIN,TTF,X0,XF,NX,Y0,YF,NY,LL)
      IMPLICIT REAL *8(a-h,o-z)

      N=7
      NALVO=2

c     2**20
      NLIN=1048576
c     periodo medio de Bianca / 6
      TTF=0.942422d0/6.0d0

      X0=48000.0d0
      XF=48400.0d0
      NX=301

      Y0=0.0d0
      YF=0.00025d0
      NY=101

      LL=12
      return
      end

c-----------------------------------------------------------------------
c     MASSAS E CONDICOES INICIAIS  (ver condicoes/naiad.f)
c       SEMI = semi-eixo do ponto (km),  YV = excentricidade
c-----------------------------------------------------------------------
      SUBROUTINE INICIAIS(SEMI,YV,N,dm0,G,dj2,dj4,dm,ami,dmimi,
     |                    a,e,di,w,om,am,xpla)
      IMPLICIT REAL *8(a-h,o-z)
      include 'dimensoes.inc'
      DIMENSION dm(NMAX),ami(NMAX),dmimi(NMAX),a(NMAX),e(NMAX),
     |di(NMAX),w(NMAX),om(NMAX),am(NMAX),xpla(NMAX,6)

c     Planeta: massa (kg) e raio equatorial (km)
      dm0=1.0d0
      dmplaneta=1.02409d26
      requat=24764.0d0
      CALL GRAVIT(dmplaneta,requat,G)

c     J2, J4  (French et al. 2024)
      dj2= 0.0035363d0
      dj4=-0.000036d0

c     POSVEL(i, x,y,z, vx,vy,vz, xpla)   (tudo dividido por requat)

c NAIAD
c     ATENCAO: estes numeros parecem estar em UA e UA/dia, nao em km
c     (3.18d-4 UA = 47600 km, a orbita de Naiad), mas sao divididos
c     por requat como se fossem km. E o ultimo valor esta escrito
c     "6.043612890766993-4" (sem o d, ou seja 6.04-4 = 2.04).
c     Mantido igual ao original - confira!
c     A massa no programa principal era 225.8d15 kg e na FORCE era
c     3.75d22 kg; aqui vale a mesma para tudo (225.8d15).
      CALL MASSA(1,225.8d15/dmplaneta,dm0,G,dm,ami,dmimi)
      CALL POSVEL(1,
     | 3.182199056200447d-4/requat,
     | 5.144655193209696d-5/requat,
     | -6.773700274705878d-6/requat,
     | -1.080962242181789d-3/requat,
     | 6.767104688172076d-3/requat,
     | 6.043612890766993-4/requat, xpla)

c Corpo 2  (corpo da grade)
      CALL MASSA(2,1.5d22/dmplaneta,dm0,G,dm,ami,dmimi)
      CALL POSVEL(2,
     | 1.887941912157757d5/requat,
     | -0.5841476965655777d5/requat,
     | 0.06214488374540088d3/requat,
     | 3.532381786778981d5/requat,
     | 11.45283973897300d5/requat,
     | 0.2317890517421091d3/requat, xpla)

c Desdemona
      CALL MASSA(3,10.805d22/dmplaneta,dm0,G,dm,ami,dmimi)
      CALL POSVEL(3,
     | 1.315147962185523d5/requat,
     | 1.979561102139357d5/requat,
     | -0.02025080133460870d3/requat,
     | -9.079252699042778d5/requat,
     | 6.089688696211173d5/requat,
     | 0.07746824354835020d3/requat, xpla)

c Juliet
      CALL MASSA(4,61.76d22/dmplaneta,dm0,G,dm,ami,dmimi)
      CALL POSVEL(4,
     | 2.627135680595500d5/requat,
     | -1.334148467583099d5/requat,
     | 2.872359096420754d3/requat,
     | 4.440804163548058d5/requat,
     | 8.743554362614603/requat,
     | -16.09720050079744d3/requat, xpla)

c Portia
      CALL MASSA(5,230.9d22/dmplaneta,dm0,G,dm,ami,dmimi)
      CALL POSVEL(5,
     | -2.552969495926177d5/requat,
     | -4.607326039145249d5/requat,
     | 2.371330898317865d3/requat,
     | 6.414327329812178d5/requat,
     | -3.558051547767367d5/requat,
     | 3.137987742929639d3/requat, xpla)

c Rosalinda
      CALL MASSA(6,109.572d22/dmplaneta,dm0,G,dm,ami,dmimi)
      CALL POSVEL(6,
     | -3.526019345195424d5/requat,
     | -1.364575008702126d5/requat,
     | 0.09870384249549534d3/requat,
     | 3.131938716837554d5/requat,
     | -8.062054103890521d5/requat,
     | 0.3889053098793468d3/requat, xpla)

c S/2025 U1
      CALL MASSA(7,13455.3d22/dmplaneta,dm0,G,dm,ami,dmimi)
      CALL POSVEL(7,
     | -11.87700282145936d5/requat,
     | 1.222192648444817d5/requat,
     | -8.050225814608742d3/requat,
     | -0.4194237112800376d5/requat,
     | -4.907872640333754d5/requat,
     | 1.024400940643259d3/requat, xpla)

c     Corpo da grade: cartesianas -> elementos, troca a e e pelos
c     valores da grade e volta para cartesianas
      CALL XYZ2ELEM1(2,xpla,ami,a,e,di,w,om,am)
      a(2)=SEMI/requat
      e(2)=YV
      CALL ELEM2XYZ1(2,a,e,di,w,om,am,ami,xpla)

      return
      end
