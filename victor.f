c=======================================================================
c     CONDICOES: VICTOR  (9 corpos; corpo da grade = 9, Perdita)
c     Grade: semi-eixo x excentricidade de Perdita
c     (vem de mapa_victor.f)
c=======================================================================

c-----------------------------------------------------------------------
c     PARAMETROS DA GRADE E DA INTEGRACAO  (ver condicoes/naiad.f)
c-----------------------------------------------------------------------
      SUBROUTINE CONFIG(N,NALVO,NLIN,TTF,X0,XF,NX,Y0,YF,NY,LL)
      IMPLICIT REAL *8(a-h,o-z)

      N=9
      NALVO=9

c     2**16
      NLIN=65536
c     periodo medio de Bianca / 6
      TTF=0.435007d0/6.0d0

      X0=76416.0d0
      XF=76417.0d0
      NX=101

      Y0=0.0d0
      YF=0.01d0
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

      pi=4.0d0*datan(1.0d0)
      conv=pi/180.0d0

c     Planeta: massa (kg) e raio equatorial (km)
      dm0=1.0d0
      dmplaneta=1898.125d24
      requat=71492.0d0
      CALL GRAVIT(dmplaneta,requat,G)

c     J2, J4 jup
      dj2=14696.5d-6
      dj4=-586.6d-6

c: lambda=am+varpi, logo am=lambda-varpi

c Bianca
      CALL MASSA(1,7.35d-10,dm0,G,dm,ami,dmimi)
      a(1)=59165.0d0/requat
      e(1)=0.00027d0
      di(1)=0.1811d0*conv
      w(1)=(274.65569d0-207.44872d0)*conv
      om(1)=207.44872d0*conv
      am(1)=(137.70117d0-274.65569d0)*conv

c Cressida
      CALL MASSA(2,21.18d-10,dm0,G,dm,ami,dmimi)
      a(2)=61767.0d0/requat
      e(2)= 0.00056d0
      di(2)=0.05023d0*conv
      w(2)=(274.65569d0-207.44872d0)*conv
      om(2)=207.44872d0*conv
      am(2)=(137.70117d0-274.65569d0)*conv

c Desdemona
      CALL MASSA(3,14.25d-10,dm0,G,dm,ami,dmimi)
      a(3)=62659.0d0/requat
      e(3)=0.00073d0
      di(3)=0.04673d0*conv
      w(3)=(157.93454d0-122.42866d0)*conv
      om(3)=122.42866d0*conv
      am(3)=(352.89892-157.93454d0)*conv

c Juliet
      CALL MASSA(4,44.59d-10,dm0,G,dm,ami,dmimi)
      a(4)=64358.0d0/requat
      e(4)= 0.00122d0
      di(4)=0.03518d0*conv
      w(4)=(89.97098d0-144.28899d0)*conv
      om(4)=144.28899d0*conv
      am(4)=(76.17466d0-89.97098d0)*conv

c Portia
      CALL MASSA(5,134.44d-10,dm0,G,dm,ami,dmimi)
      a(5)= 66097.0d0/requat
      e(5)= 0.00051d0
      di(5)=0.01643d0*conv
      w(5)=(6.91932d0-31.75335d0)*conv
      om(5)= 31.75335d0*conv
      am(5)=(50.43991d0-6.91932d0) *conv

c Rosalinda
      CALL MASSA(6,20.26d-10,dm0,G,dm,ami,dmimi)
      a(6)= 69927.0d0/requat
      e(6)=0.00090d0
      di(6)=0.04922d0*conv
      w(6)=(239.72554d0-249.40149d0)*conv
      om(6)=249.40149d0*conv
      am(6)=(139.85182d0-239.72554d0) *conv

c Cupid
c Cuk 2022     dm(7)=(28.46d-10/125.0d0)
c French et al. 2015 via Charalambous
      CALL MASSA(7,0.295297d16/dmplaneta,dm0,G,dm,ami,dmimi)
      a(7)= 74392.338d0/requat
      e(7)=0.00047d0
      di(7)=0.07028d0*conv
      w(7)=(320.20439d0-59.71798d0)*conv
      om(7)=59.71798d0*conv
      am(7)=(109.06840d0-320.20439d0) *conv

c Belinda
      CALL MASSA(8,28.46d-10,dm0,G,dm,ami,dmimi)
      a(8)= 75255.0d0/requat
      e(8)= 0.00079d0
      di(8)=0.00172d0*conv
      w(8)=(38.78082d0-232.28659d0)*conv
      om(8)=232.28659d0*conv
      am(8)=(324.66202d0-38.78082d0)*conv

c Perdita  (corpo da grade)
c Cuk 2022       dm(9)=(28.46d-10/27.0d0)
c French et al. 2015 via Charalambous
      CALL MASSA(9,0.985470d16/dmplaneta,dm0,G,dm,ami,dmimi)
      a(9)=SEMI/requat
      e(9)=YV
      di(9)=0.05889d0*conv
      w(9)= (338.81201d0-248.66870d0)*conv
      om(9)=248.66870d0*conv
      am(9)=(67.48953d0-338.81201d0) *conv

c     elementos -> cartesianas (xpla) para todos os corpos
      CALL ELEM2XYZ(N,a,e,di,w,om,am,ami,xpla)

      return
      end

c=======================================================================
c     Puxa o codigo base (nao mexer). Para compilar basta:
c        gfortran -o mapa_victor victor.f
c=======================================================================
      include 'mapa_base.f'
