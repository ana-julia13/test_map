c=======================================================================
c     MAPA  -  CODIGO BASE (NAO PRECISA MEXER AQUI)
c
c     Integra N corpos em torno de um planeta achatado (J2, J4) com o
c     RA15, para cada ponto de uma grade (semi-eixo x variavel Y) do
c     corpo NALVO, e calcula:
c       - emax, amax, imax do corpo NALVO
c       - valores medios e numeros espectrais (FFT) de a, e, I, e*cos(varpi)
c
c     Tudo que muda de um caso para outro (numero de corpos, massas,
c     condicoes iniciais, grade, passo, J2/J4 ...) fica num arquivo
c     separado em  condicoes/<caso>.f , que deve conter:
c        SUBROUTINE CONFIG(...)    -> parametros da grade/integracao
c        SUBROUTINE INICIAIS(...)  -> massas e condicoes iniciais
c     Ver README.md.
c=======================================================================
      PROGRAM MAPA
      IMPLICIT REAL *8(a-h,o-z)
      include 'dimensoes.inc'
      DIMENSION X(6*NMAX),V(6*NMAX)
      dimension a(NMAX),e(NMAX),di(NMAX),w(NMAX),om(NMAX),am(NMAX),
     |ami(NMAX),xpla(NMAX,6),bx(NMAX,1),by(NMAX,1),bz(NMAX,1)
      REAL*8, ALLOCATABLE :: SERIE(:,:)
      INTEGER NESPEC(7)
      dimension xmedio(4)
      COMMON/MASSAS/dm0,G,dj2,dj4
      COMMON/IN/ikepler,inicio,pi,pi05,pi15,pi2
      COMMON/CONST/conv
      COMMON/NUMERO/N,neq
      COMMON/CORPOS/dm(NMAX),dmimi(NMAX)

c**** PARAMETROS DO CASO (arquivo de condicoes)
      CALL CONFIG(N,NALVO,NLIN,TTF,X0,XF,NX,Y0,YF,NY,LL)

      if(N.lt.1 .or. N.gt.NMAX) then
         write(*,*)'ERRO: N=',N,' fora de 1..NMAX=',NMAX,
     |   ' (aumente NMAX em base/dimensoes.inc)'
         stop
      endif
      if(NALVO.lt.1 .or. NALVO.gt.N) then
         write(*,*)'ERRO: NALVO=',NALVO,' fora de 1..N=',N
         stop
      endif
      if(NLIN.lt.4 .or. iand(NLIN,NLIN-1).ne.0) then
         write(*,*)'ERRO: NLIN=',NLIN,' tem que ser potencia de 2'
         stop
      endif
      if(NX.lt.1 .or. NY.lt.1) then
         write(*,*)'ERRO: NX e NY tem que ser >= 1'
         stop
      endif

      neq=6*N
      tf=ttf
      ALLOCATE(SERIE(NLIN,4))

      dx=0.0d0
      dy=0.0d0
      if(NX.gt.1) dx=(XF-X0)/dble(NX-1)
      if(NY.gt.1) dy=(YF-Y0)/dble(NY-1)

      CALL ARQQ(NX,NY,X0,XF,Y0,YF)

      TM=0.0d0
      XL=0.0d0

c**** LOOPS DA GRADE: Y (fora) x SEMI-EIXO (dentro)
      yv=Y0
      do 900 iy=1,NY
      xv=X0
      do 800 ix=1,NX

c     constantes numericas
      ikepler=1
      inicio=1
      kk=1
      pi=4.0d0*datan(1.0d0)
      pi2=pi+pi
      pi05=0.5d0*pi
      pi15=1.5d0*pi
      conv=pi/180.0d0

c**** MASSAS E CONDICOES INICIAIS (arquivo de condicoes)
      CALL INICIAIS(xv,yv,N,dm0,G,dj2,dj4,dm,ami,dmimi,
     |              a,e,di,w,om,am,xpla)

c     comeca com emax=excentricidade inicial, amax=..., imax=...
      emax=e(NALVO)
      amax=a(NALVO)
      dimax=di(NALVO)
      ttt=0.0d0

c     COORDENADAS a serem utilizadas na subroutine FORCE:
      k=1
      do 2 i=1,N
      do 1 j=4,6
      x(k)=xpla(i,j)
      k=k+1
1     continue
      k=k+3
2     continue

c     Velocidades relativas -> momenta canonicos (gaussj)
      do 4 i=1,N
      bx(i,1)=xpla(i,1)
      by(i,1)=xpla(i,2)
      bz(i,1)=xpla(i,3)
4     continue
      CALL gaussj(N,bx,dmimi)
      CALL gaussj(N,by,dmimi)
      CALL gaussj(N,bz,dmimi)

      k=4
      do 5 l=1,N
      x(k)=bx(l,1)
      x(k+1)=by(l,1)
      x(k+2)=bz(l,1)
      k=k+6
5     continue

      CALL XPLAX(N,X,xpla)
      CALL XYZORB(N,xpla,a,e,di,w,om,am,ami,lei,kk)
      CALL GUARDA(SERIE,NLIN,1,NALVO,a,e,di,w,om)

c**** INTEGRACAO NUMERICA
      t=0.0d0
      do 700 jjj=2,NLIN
      CALL RA15(TM,X,V,TF,LL,XL,neq,1)
      t=t+tf

      CALL XPLAX(N,X,xpla)
      CALL XYZORB(N,xpla,a,e,di,w,om,am,ami,lei,kk)

      IF(e(NALVO).GT.EMAX) then
      EMAX=e(NALVO)
      TTT=T
      ENDIF
      IF(a(NALVO).GT.AMAX) then
      AMAX=a(NALVO)
      TTT=T
      ENDIF
      IF(di(NALVO).GT.DIMAX) then
      DIMAX=di(NALVO)
      TTT=T
      ENDIF

      CALL GUARDA(SERIE,NLIN,jjj,NALVO,a,e,di,w,om)
700   continue

c**** SAIDAS: maximos
      write(923,*)ttt,ttt/365.25d0,xv,yv,emax
      write(9234,*)ttt,ttt/365.25d0,xv,yv,amax
      write(92345,*)ttt,ttt/365.25d0,xv,yv,dimax/conv

      write(123,*)xv,yv,emax
      write(1234,*)xv,yv,amax
      write(12345,*)xv,yv,dimax/conv

      write(1239,*)emax
      write(12349,*)amax
      write(23459,*)dimax/conv

c**** SAIDAS: FFT de a (1), e (2), I (3), e*cos(varpi) (4)
      do 600 iv=1,4
      CALL FFT(SERIE(1,iv),NLIN,TTF,xm,NESPEC)
      xmedio(iv)=xm
      if(iv.ge.3) xm=xm/conv
      write(3000+iv,*)xv,yv,xm
      write(3010+iv,*)xm
      do 550 it=1,7
      write(1000+10*iv+it,*)NESPEC(it)
      write(2000+10*iv+it,*)yv,xv,NESPEC(it)
550   continue
600   continue
      write(3020,*)xmedio(1),xmedio(2),xmedio(4),xmedio(3)/conv,xv,yv

      xv=xv+dx
800   continue
      yv=yv+dy
900   continue

      DEALLOCATE(SERIE)
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Converte o vetor X (coordenadas e momenta) em xpla:
c     xpla(i,1..3)=velocidades, xpla(i,4..6)=coordenadas
      SUBROUTINE XPLAX(N,X,xpla)
      IMPLICIT REAL *8(a-h,o-z)
      include 'dimensoes.inc'
      DIMENSION X(6*NMAX),xpla(NMAX,6)
      COMMON/CORPOS/dm(NMAX),dmimi(NMAX)
      j=1
      do 25 idef=1,N
      do 24 jdef=1,3
      xpla(idef,jdef)=x(j+3)/dm(idef)
      xpla(idef,jdef+3)=x(j)
      j=j+1
24    continue
      j=j+3
25    continue
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Guarda os elementos do corpo NALVO na linha jjj das series da FFT
      SUBROUTINE GUARDA(SERIE,NLIN,jjj,NALVO,a,e,di,w,om)
      IMPLICIT REAL *8(a-h,o-z)
      include 'dimensoes.inc'
      DIMENSION SERIE(NLIN,4),a(NMAX),e(NMAX),di(NMAX),w(NMAX),om(NMAX)
      SERIE(jjj,1)=a(NALVO)
      SERIE(jjj,2)=e(NALVO)
      SERIE(jjj,3)=di(NALVO)
      SERIE(jjj,4)=e(NALVO)*dcos(w(NALVO)+om(NALVO))
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Abre os arquivos de saida e escreve os cabecalhos DSAA (Surfer)
      SUBROUTINE ARQQ(NX,NY,X0,XF,Y0,YF)
      IMPLICIT REAL *8(a-h,o-z)
      CHARACTER*8 SUF(4),LAB(7)
      CHARACTER*12 VMED(4)
      DATA SUF/'a','','i','evarpi'/
      DATA LAB/'01','05','1','3','5','10','15'/
      DATA VMED/'Amedio','Emedio','Imedio','EcosVmedio'/

      OPEN(123,FILE='emax.dat',STATUS='UNKNOWN')
      OPEN(1234,FILE='amax.dat',STATUS='UNKNOWN')
      OPEN(12345,FILE='imax.dat',STATUS='UNKNOWN')

      OPEN(1239,FILE='emax.grd',STATUS='UNKNOWN')
      OPEN(12349,FILE='amax.grd',STATUS='UNKNOWN')
      OPEN(23459,FILE='imax.grd',STATUS='UNKNOWN')

      OPEN(923,FILE='emax-t.dat',STATUS='UNKNOWN')
      OPEN(9234,FILE='amax-t.dat',STATUS='UNKNOWN')
      OPEN(92345,FILE='imax-t.dat',STATUS='UNKNOWN')

      CALL DSAA(1239,NX,NY,X0,XF,Y0,YF)
      CALL DSAA(12349,NX,NY,X0,XF,Y0,YF)
      CALL DSAA(23459,NX,NY,X0,XF,Y0,YF)

c     valores medios (termo de frequencia zero da FFT)
      do 10 iv=1,4
      OPEN(3000+iv,FILE=TRIM(VMED(iv))//'.DAT',STATUS='UNKNOWN')
      OPEN(3010+iv,FILE=TRIM(VMED(iv))//'.grd',STATUS='UNKNOWN')
      CALL DSAA(3010+iv,NX,NY,X0,XF,Y0,YF)
10    continue
      OPEN(3020,FILE='AmedioEmedioEVmedioImedio.DAT',STATUS='UNKNOWN')

c     numeros espectrais: dn<lim><var>.grd e e1e2dn<lim><var>.dat
c     var: a (semi-eixo), '' (excentricidade), i (inclinacao),
c          evarpi (e*cos(varpi))
      do 30 iv=1,4
      do 20 it=1,7
      OPEN(1000+10*iv+it,
     |FILE='dn'//TRIM(LAB(it))//TRIM(SUF(iv))//'.grd',STATUS='UNKNOWN')
      OPEN(2000+10*iv+it,
     |FILE='e1e2dn'//TRIM(LAB(it))//TRIM(SUF(iv))//'.dat',
     |STATUS='UNKNOWN')
      CALL DSAA(1000+10*iv+it,NX,NY,X0,XF,Y0,YF)
20    continue
30    continue

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE DSAA(iu,NX,NY,X0,XF,Y0,YF)
      IMPLICIT REAL *8(a-h,o-z)
      write(iu,'(A)')'DSAA'
      write(iu,*)NX,NY
      write(iu,*)X0,XF
      write(iu,*)Y0,YF
      write(iu,*)'1 100'
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     FFT de uma serie de NN pontos (NN potencia de 2) espacados de DELTAT.
c     Saida: XMEDIO = valor medio (frequencia zero)
c            NESPEC(1..7) = numeros espectrais (picos com periodo > PERMIN
c            e amplitude > 0.1, 0.5, 1, 3, 5, 10, 15 % da maior amplitude)
      SUBROUTINE FFT(SERIE,NN,DELTAT,XMEDIO,NESPEC)
      IMPLICIT REAL *8(a-h,o-z)
      INTEGER NN,NESPEC(7)
      DIMENSION SERIE(NN),FATOR(7)
      REAL*8, ALLOCATABLE :: data(:),gg(:),FRE(:),AMP(:),FFRE(:),AAMP(:)
      PARAMETER (PERMIN=5.0d0)
      DATA FATOR/0.1d0,0.5d0,1.0d0,3.0d0,5.0d0,10.0d0,15.0d0/

      ALLOCATE(data(2*NN),gg(NN+2),FRE(NN+2),AMP(NN+2),FFRE(NN+2),
     |AAMP(NN+2))
      data=0.0d0
      gg=0.0d0
      FRE=0.0d0
      AMP=0.0d0
      FFRE=0.0d0
      AAMP=0.0d0

      i=1
      do 10 n=0,NN-1
      data(i)=SERIE(n+1)
      data(i+1)=0.0d0
      if(n.eq.0) gg(n+1)=0.0d0
      if(n.ge.1. and. n.le.NN/2-1) gg(n+1)=n/(NN*deltat)
      if(n.eq.NN/2) gg(n+1)=1.0d0/(2.0d0*deltat)
      i=i+2
10    continue

      isign=1
      CALL FOUR1(data,NN,isign)

      i=1
      do 20 n=0,NN-1
      tf=2.0d0*((data(i)**2+data(i+1)**2)**0.5d0)/NN
      if(n.eq.0) XMEDIO=tf/2.0d0
      if(n.ge.1.and.n.le.NN/2) then
      FRE(n+1)=gg(n+1)
      AMP(n+1)=tf
      endif
      i=i+2
20    continue
      ii=NN+1

cccccccccc PROCURA PELOS PICOS DE FREQUENCIA
      J=1
      do jj=2,ii
      a1=AMP(jj-1)
      a =AMP(jj)
      a2=AMP(jj+1)
      if(a.gt.a1. and. a2.lt.a) then
      AAMP(J)=AMP(jj)
      FFRE(J)=FRE(jj)
      J=J+1
      endif
      enddo
      NPICOS=J-1

cccccccccc PROCURA PELA amplitude MAXIMA (so periodos > PERMIN)
      ampmax=-1.
      do JJ=1,NPICOS
      peri=1.0d0/FFRE(JJ)
      ampl=AAMP(JJ)
      if(peri.gt.PERMIN. and. ampl.GT.ampmax) ampmax=ampl
      enddo
      amplim=ampmax/100.0d0

cccccccccc NUMEROS ESPECTRAIS
      do k=1,7
      NESPEC(k)=0
      enddo
      do J=1,NPICOS
      ampl=AAMP(J)
      peri=1.0d0/FFRE(J)
      do k=1,7
      if(peri.gt.PERMIN.and.ampl.gt.(FATOR(k)*amplim))
     |NESPEC(k)=NESPEC(k)+1
      enddo
      enddo

      DEALLOCATE(data,gg,FRE,AMP,FFRE,AAMP)
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     EQUACOES DE MOVIMENTO CANONICAS PARA N CORPOS
c     EM VARIAVEIS CANONICAS RELATIVAS AO PLANETA
c     x(6k-5..6k-3) = coordenadas do corpo k
c     x(6k-2..6k)   = momenta do corpo k
      SUBROUTINE FORCE (TM,X,V,F1)
      IMPLICIT REAL *8(a-h,o-z)
      include 'dimensoes.inc'
      DIMENSION X(6*NMAX),V(6*NMAX),F1(6*NMAX)
      DIMENSION DDJ2(3),DDJ4(3)
      INTEGER II
      COMMON/MASSAS/dm0,G,dj2,dj4
      COMMON/NUMERO/N,neq
      COMMON/CORPOS/dm(NMAX),dmimi(NMAX)

      l=1
      i=1

1     continue
      if(l+5.gt.neq. and .i.gt.N) go to 2
      F1(l)=dmimi(i)*x(l+3)
      F1(l+1)=dmimi(i)*x(l+4)
      F1(l+2)=dmimi(i)*x(l+5)

      j=1
      smx=0.0d0
      smy=0.0d0
      smz=0.0d0
      do while(j+5.le.neq)
      if (j.eq.l) go to 21
      smx=smx+x(j+3)
      smy=smy+x(j+4)
      smz=smz+x(j+5)
21    continue
      j=j+6
      enddo

      F1(l)=F1(l)+smx/dm0
      F1(l+1)=F1(l+1)+smy/dm0
      F1(l+2)=F1(l+2)+smz/dm0

      ri3=(x(l)**2+x(l+1)**2+x(l+2)**2)**(1.5d0)

      F1(l+3)=-G*dm(i)*dm0*x(l)/ri3
      F1(l+4)=-G*dm(i)*dm0*x(l+1)/ri3
      F1(l+5)=-G*dm(i)*dm0*x(l+2)/ri3

      j=1
      do while(j+5.le.neq)
      if(j.eq.l) go to 12
      m=(j-1)/6+1

      dij3=( (x(l)-x(j))**2 + (x(l+1)-x(j+1))**2
     | + (x(l+2)-x(j+2))**2 )**(1.5d0)

      F1(l+3)=F1(l+3)-G*dm(i)*dm(m)*(x(l)-x(j))/dij3
      F1(l+4)=F1(l+4)-G*dm(i)*dm(m)*(x(l+1)-x(j+1))/dij3
      F1(l+5)=F1(l+5)-G*dm(i)*dm(m)*(x(l+2)-x(j+2))/dij3

12    continue
      j=j+6
      enddo

      II=i
      CALL ACHAT(X,DDJ2,DDJ4,II)

      F1(l+3)=F1(l+3)-DDJ2(1)-DDJ4(1)
      F1(l+4)=F1(l+4)-DDJ2(2)-DDJ4(2)
      F1(l+5)=F1(l+5)-DDJ2(3)-DDJ4(3)

      l=l+6
      i=i+1
      go to 1
2     continue

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Forcas devidas ao achatamento (J2, J4) sobre o corpo II
      SUBROUTINE ACHAT(X,DDJ2,DDJ4,II)
      IMPLICIT REAL *8(a-h,o-z)
      include 'dimensoes.inc'
      DIMENSION X(6*NMAX),DDJ2(3),DDJ4(3)
      INTEGER II
      COMMON/MASSAS/dm0,G,dj2,dj4
      COMMON/CORPOS/dm(NMAX),dmimi(NMAX)

      x1=x(6*II-5)
      y1=x(6*II-4)
      z1=x(6*II-3)
      r1=dsqrt(x1**2.0d0+y1**2.0d0+z1**2.0d0)
      dMm=dm(II)

c     Derivadas que entram no achatamento

c     coordenada x
      r17=r1**7.0d0
      DPX1=-5.0d0*x1/r17
      r15=r1**5.0d0
      DPX2=-3.0d0*x1/r15
      r11=r1**11.0d0
      DPX3=-9.0d0*x1/r11
      r19=r1**9.0d0
      DPX4=-7.0d0*x1/r19

c     coordenada y
      DPY1=-5.0d0*y1/r17
      DPY2=-3.0d0*y1/r15
      DPY3=-9.0d0*y1/r11
      DPY4=-7.0d0*y1/r19

c     coordenada z
      DPZ1=-5.0d0*z1/r17
      DPZ2=-3.0d0*z1/r15
      DPZ3=-9.0d0*z1/r11
      DPZ4=-7.0d0*z1/r19

c     Termos que aparecem do achatamento
      cj2=G*dm0*dMm*dj2/2.0d0
      cj4=G*dm0*dMm*dj4/8.0d0

      dj2x1=3.0d0*z1*z1*DPX1-DPX2
      DDJ2(1)=cj2*dj2x1

      dj2y1=3.0d0*z1*z1*DPY1-DPY2
      DDJ2(2)=cj2*dj2y1

      dj2z1=3.0d0*( z1*z1*DPZ1+2.0d0*z1/(r1**5.0d0) )-DPZ2
      DDJ2(3)=cj2*dj2z1

      dj4x1=35.0d0*(z1**4.0d0)*DPX3-
     -      30.0D0*(z1**2.0d0)*DPX4+
     +      3.0d0*DPX1
      DDJ4(1)=cj4*dj4x1

      dj4y1=35.0d0*(z1**4.0d0)*DPY3-
     -      30.0D0*(z1**2.0d0)*DPY4+
     +      3.0d0*DPY1
      DDJ4(2)=cj4*dj4y1

      dj4z1=35.0d0*( (z1**4.0d0)*DPZ3+4.0d0*(z1**3.0d0)/(r1**9.0d0) )-
     ~      30.0D0*( (z1**2.0d0)*DPZ4+2.0d0*z1/(r1**7.0d0) )+
     ~      3.0d0*DPZ1
      DDJ4(3)=cj4*dj4z1

      return
      end

c=======================================================================
c     ROTINAS DE AJUDA PARA USAR NO ARQUIVO DE CONDICOES
c=======================================================================

c     G nas unidades do programa: distancia em raios equatoriais do
c     planeta (requat, em km), massa em massas do planeta (dmplaneta,
c     em kg), tempo em dias.
      SUBROUTINE GRAVIT(dmplaneta,requat,G)
      IMPLICIT REAL *8(a-h,o-z)
      pi=4.0d0*datan(1.0d0)
      grav=4.0d0*pi*pi
c     GRAV acima em UA**3/(Massa Solar*Ano**2)
      ua=1.49597870691d11
      dmsol=1988500.0d24
      dmsatt=dmsol/dmplaneta
      esc=ua/(requat*1000.0d0)
      grav=grav*esc*esc*esc
      G=grav/(dmsatt*365.25d0*365.25d0)
      return
      end

c     Define a massa do corpo i (em massas do planeta) e calcula
c     ami(i)=G*(dm0+dm(i)) e dmimi(i)=(dm0+dm(i))/(dm0*dm(i))
      SUBROUTINE MASSA(i,dmi,dm0,G,dm,ami,dmimi)
      IMPLICIT REAL *8(a-h,o-z)
      include 'dimensoes.inc'
      DIMENSION dm(NMAX),ami(NMAX),dmimi(NMAX)
      dm(i)=dmi
      ami(i)=G*(dm0+dm(i))
      dmimi(i)=(dm0+dm(i))/(dm0*dm(i))
      return
      end

c     Elementos orbitais (a,e,di,w,om,am) -> xpla para os corpos 1..N
      SUBROUTINE ELEM2XYZ(N,a,e,di,w,om,am,ami,xpla)
      IMPLICIT REAL *8(a-h,o-z)
      include 'dimensoes.inc'
      DIMENSION a(NMAX),e(NMAX),di(NMAX),w(NMAX),om(NMAX),am(NMAX),
     |ami(NMAX),xpla(NMAX,6),f(NMAX)
      CALL ORBXYZ(N,a,e,di,w,om,am,ami,xpla,f)
      return
      end

c     Elementos orbitais -> xpla so para o corpo i
      SUBROUTINE ELEM2XYZ1(i,a,e,di,w,om,am,ami,xpla)
      IMPLICIT REAL *8(a-h,o-z)
      include 'dimensoes.inc'
      DIMENSION a(NMAX),e(NMAX),di(NMAX),w(NMAX),om(NMAX),am(NMAX),
     |ami(NMAX),xpla(NMAX,6),f(NMAX)
      CALL ORBXYZ1(i,a,e,di,w,om,am,ami,xpla,f)
      return
      end

c     Posicao (x,y,z) e velocidade (vx,vy,vz) do corpo i -> xpla
      SUBROUTINE POSVEL(i,px,py,pz,vx,vy,vz,xpla)
      IMPLICIT REAL *8(a-h,o-z)
      include 'dimensoes.inc'
      DIMENSION xpla(NMAX,6)
      xpla(i,1)=vx
      xpla(i,2)=vy
      xpla(i,3)=vz
      xpla(i,4)=px
      xpla(i,5)=py
      xpla(i,6)=pz
      return
      end

c     xpla do corpo i -> elementos orbitais do corpo i (Fitzpatrick)
      SUBROUTINE XYZ2ELEM1(i,xpla,ami,a,e,di,w,om,am)
      IMPLICIT REAL *8(a-h,o-z)
      include 'dimensoes.inc'
      DIMENSION a(NMAX),e(NMAX),di(NMAX),w(NMAX),om(NMAX),am(NMAX),
     |ami(NMAX),xpla(NMAX,6),xk(6)
      do 1 k=1,6
      xk(k)=xpla(i,k)
1     continue
      CALL FITZPAORB(xk,a(i),e(i),di(i),w(i),alat,am(i),u,om(i),ami(i))
      return
      end

c=======================================================================
c     ROTINAS ORIGINAIS (conversoes, Kepler, gaussj, FFT, RA15)
c=======================================================================
      SUBROUTINE four1(data,nn,isign)
      INTEGER isign,nn
c      REAL data(2*nn)
	REAL*8 data(2*nn)
      INTEGER i,istep,j,m,mmax,n
c      REAL*8 tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp,tempi,tempr
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2

2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep

      goto 2
      endif
      return
      END


cccccccccccccccccccccccccccccccccccccccccccc     

      SUBROUTINE gaussj(N,b,dmimi)
      implicit real *8(a-h,o-z)
      include 'dimensoes.inc'
      dimension dmimi(NMAX),a(NMAX,NMAX),b(NMAX,1)
      INTEGER m,n,i,icol,irow,j,k,l,ll,indxc(NMAX),
     +indxr(NMAX),ipiv(NMAX)
      REAL big,dum,pivinv


C          LINEAR EQUATION SOLUTION BY GAUSS-JORDAN ELIMINATION
C          *****************NUMERICAL RECIPES******************
      
C     a(1:n,1:n) e' uma matriz de entrada de dimensao np X np
C     b(1:n,1:m) e' uma matriz de entrada de dimensao np X mp

c     Saidas:  a(1:n,1:n) e' substituida por sua inversa
c              b(1:n,1:m) e' substituida pela solucao 

      i=0
      j=0
	k=1
	do 2 i=1,n
	do 1 j=1,n

	if(i.eq.j) then
	a(i,j)=dmimi(k)
	k=k+1
	else
	a(i,j)=1.0d0
	endif

1     continue
2     continue               
      i=0
      j=0
      k=0       
cccccccccccc
      m=1
CCCCCCCCCCCC

      do j=1,n
          ipiv(j)=0
      enddo
      
      
      do i=1,n
      
          big=0.
          do j=1,n
              if(ipiv(j).ne.1) then
                  do k=1,n
                      if (ipiv(k).eq.0) then
                          if (abs(a(j,k)).ge.big) then
                              big=abs(a(j,k))
                              irow=j
                              icol=k
                          endif
                      else if (ipiv(k).gt.1) then
c                          write(*,*) 'singular matriz in gaussj'
c                          pause
                      endif
                  enddo
               endif
           enddo
     
      ipiv(icol)=ipiv(icol)+1
      if(irow.ne.icol) then
          do l=1,n
              dum=a(irow,l)
              a(irow,l)=a(icol,l)
              a(icol,l)=dum
          enddo
          do l=1,m
              dum=b(irow,l) 
              b(irow,l)=b(icol,l)
              b(icol,l)=dum
          enddo
      endif
      indxr(i)=irow
      indxc(j)=icol
      if  (a(icol,icol).eq.0) then
c                              write(*,*)'singular matriz in gaussj II'
c                              pause
	endif
      pivinv=1./a(icol,icol)
      a(icol,icol)=1.
      do l=1,n 
          a(icol,l)=a(icol,l)*pivinv
      enddo
      do l=1,m
          b(icol,l)=b(icol,l)*pivinv
      enddo
      do ll=1,n
          if(ll.ne.icol) then
              dum=a(ll,icol)
              a(ll,icol)=0.
              do l=1,n
                  a(ll,l)=a(ll,l)-a(icol,l)*dum
              enddo
              do l=1,m
                  b(ll,l)=b(ll,l)-b(icol,l)*dum
              enddo
          endif
      enddo
      
      enddo
c      do l=n,1,-1
c	write(*,*) 'caca final'
c	write(*,*) n, l, indxr(l), indxc(l)
c          if (indxr(l).ne.indxc(l)) then
c              do k=1,n
c                  dum=a(k,indxr(l))   
c                  a(k,indxr(l))=a(k,indxc(l))
c                  a(k,indxc(l))=dum   
c              enddo
c          endif
c      enddo
      return
      END    

      subroutine orbxyz(N,a,e,dincl,w,om,am,ami,x,f)
      implicit real *8(a-h,o-z)
      include 'dimensoes.inc'
      COMMON/IN/ikepler,inicio,pi,pi05,pi15,pi2
      dimension a(NMAX),e(NMAX),w(NMAX),om(NMAX),am(NMAX),ami(NMAX),
     |f(NMAX),x(NMAX,6),u(NMAX),dincl(NMAX)
 
c: rotina que passa de elementos orbitais para cartesianos
c: preciso tambem das derivadas no tempo das cartesianas(velocidades)
c: referencia: Fitzpatrick (em velocidades parece que nao esta ok)

c: entrada : elementos orbitais , saida : conversao para cartesianos
c: u=anom.excent. (sera achado via rotina Kepler)
c: am=anom. media
c: om=longitude do nodo (indefinido se I=0,mas tomo zero para que as formulas
c: do caso espacial fiquem validas no caso plano)
c: w=argumento do pericentro (se I=0 este w representara longitude do pericen.)
c: alat=latitude=w+f (representara wtil+f no caso I=0) (wtil=wpi)

c: obtencao de anom. excent. (u) e verdadeira (f) via Kepler

      call kepler(f,u,am,e,N)

      do 50 i=1,N
      inicio=0
      alat=w(i)+f(i)
      r=a(i)*(1.d0-e(i)*dcos(u(i)))
      rp=dsqrt(ami(i)*a(i))*e(i)*dsin(u(i))/r
      rfp=dsqrt(ami(i)*a(i)*(1.d0-e(i)*e(i)))/r
c: velocidades (x(1),x(2),x(3) no caso tridimensional
c: Nesta rotina se I = zero ou e= 0 nao ha nenhum problema, embora
c: para entrar com angulos eu entre com Wtil (no lugar de w)
c: e em lugar de om=nodo=indefinido tomo=0
      cfw=dcos(alat)
      sfw=dsin(alat)
      coom=dcos(om(i))
      siom=dsin(om(i))
      coi=dcos(dincl(i))
      somcoi=siom*coi
      comcoi=coom*coi
      sini=dsin(dincl(i))
      cosi=dcos(dincl(i))
c: no caso plano basta tomar I=dincl=0 e definir om=0.d0
c: tambem w daqui deve ser interpretado como sendo o wpi=long. pericentro
c: tambem alat deve ser pensado como alat=wpi+f

c: velocidades

      x(i,1)=rp*cfw*coom-rfp*sfw*coom-rp*sfw*somcoi-rfp*cfw*somcoi
      
      x(i,2)=rp*cfw*siom-rfp*sfw*siom+rp*sfw*comcoi+rfp*cfw*comcoi
      
      x(i,3)=rp*sfw*sini+rfp*cfw*sini

c: coordenadas

      x(i,4)=r*(cfw*coom-sfw*cosi*siom)    
     
      x(i,5)=r*(cfw*siom+sfw*cosi*coom)
     
      x(i,6)=r*sini*sfw

50    continue      

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
C SUBROTINA PARA RESOLVER EQCAO DE KEPLER.  ADO AM ACHO U DEPOIS F
C FIZ TODOS OS TESTES : ESTA RODANDO REDONDO , U1 ,F (anomalias)
      SUBROUTINE KEPLER (f,u1,AM,E0,N)
      IMPLICIT REAL *8(A-H,O-Z) 
      include 'dimensoes.inc'
      COMMON/IN/ikepler,inicio,pi,pi05,pi15,pi2                       
      dimension f(NMAX),fa(NMAX),am(NMAX),e0(NMAX),u1(NMAX)

C      READ(*,*) AM,E0
C: metodo de Newton Raphson  
c:ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c: este bloco eu usava p/ garantir grande precisao na solucao da eq.Kepler
c: apenas p/ converter no inicio os elementos orbitais em x,y,z (qdo tomo kepl=1)
c: voce pode deixar sempre kepl=1, e garante sempre grande precisao, mas fica
c: mais lento.Fica a seu criterio tirar, mas tem que dar o LITERA e PREC)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      DO 80 j=1,N
      
      litera=15
      prec=1.d-10
      if(ikepler.eq.1) then
                           litera=35
                           prec=1.d-12
                           end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                           
      U0=AM(j)+E0(j)*DSIN(AM(j))
      DO 10 I=1,litera
      U1(j)=U0-((U0-E0(j)*DSIN(U0)-
     /AM(j)))/(1.D0-E0(j)*DCOS(U0))
      TEST1=DABS(U1(j)-U0)
      DSU1=DSIN(U1(j))
      TEST2=DABS(U1(j)-E0(j)*DSU1-AM(j))
      IF(TEST1.LE.prec.AND.TEST2.LE.prec) GO TO 11
      U0=U1(j)
 10   CONTINUE
c**      WRITE(*,*)' NAO CONVERGENCIA EM KEPLER'
c**      write(*,*)'test2    excent', test2, e0,u1
 11   CONTINUE
      RCF=DCOS(U1(j))-E0(j)
c: rcf=0 ou quase se f=anom. verdad.= proximo de pi/2,pi*3/2.Se rcf<=1.d-11
c: suporemos f =pi/2 ou pi*3/2, mas se isto ocorrer ja no inicio (primeira)
c: chamada de Kepler, posso (opcionalmente) parar a integracao entrando com
c: nova condicao inicial para o Planeta (esperando que nao ocorra o mesmo)
c: (se porem o parametro INICIO e'=0,apenas considero f= pi/2 ou pi*3/2.Esta 
c: aproximacao so' tem problema se estamos usando Kepler para o Planeta,
c: pois no caso do satelite, esta transformacao e' apenas um dado de saida.
      if(dabs(rcf).gt.1.d-11) go to 20
      if(inicio.eq.1)then
      open (10,file='keple',status='unknown')
c      write(10,*)'melhor mudar AM ou E0 inic do Planeta,para que f seja
c     +melhor determinado'
                         stop 
                         end if
      if(dsu1.gt.0.d0) f(j)=pi05
      if(dsu1.lt.0.d0) f(j)=pi15
      return
  20  continue

      RSF=DSQRT(1.D0-E0(j)*E0(j))*DSU1
      FA(j)=DATAN(RSF/RCF)
      F(j)=FA(j)
      
      IF(f(j).ge.0.D0.AND.rcf.LT.0.D0) F(j)=FA(j)+PI
      IF(f(j).lt.0.D0.AND.rcf.lt.0.D0) F(j)=FA(j)+PI
      if(f(j).lt.0.d0.and.rcf.gt.0.d0) F(j)=FA(j)+PI2
      
80    continue                                        

c: apos orbxyz chamar Kepler, ele redefine inicio=0
      RETURN
      END 
      
c: SUBROTINA XYZORB(....) 

c: NELSON, NELSONNNNNNNNNNNNNNNNNNN (13/01/98)

c: Pelo que me lembro, o programa deve funcinar p/ I=0, mas eu quase nao testei pois
c: o nosso problema era espacial. Mas se seu I=0 (sempre plano) e'facilimo mexer
c: neste aqui p/ que funcione sem galho algum. QQr coisa podemos falar depois. 
c: No primeiro caso qdo I=quase zero, zero, e nao zero o programa sempre avisa
c: onde esta a indeterminacao. Idem se excent=zero ou quase. O aviso vem com
c: lei=1 (se lei=zero, esta tudo ok). Sugiro que sempre teste. 
c: Aqui ainda existem alguns PAUSE, mas voce pode tira-los (num processo REMOTO)
c: pois, o aviso de que deu problema voce sabe com o parametro LEI (que voce
c: sempre  deve imprimir , eu sugiro).
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      subroutine xyzorb(x,a,e,dincl,w,om,am,ami,lei,kk)
c     /wtif,f,rcosf,rsinf)      
      subroutine xyzorb(N,x,a,e,dincl,w,om,am,ami,lei,kk)
c      ,wtif,alat,f,rcosf,ww)
      implicit real *8(a-h,o-z)
      include 'dimensoes.inc'
      dimension x(NMAX,6),am(NMAX),e(NMAX),a(NMAX),dincl(NMAX),w(NMAX),
     |om(NMAX),ami(NMAX),f(NMAX),u(NMAX)
      COMMON/IN/ikepler,inicio,pi,pi05,pi15,pi2   

c: referencias : Fitzpatrick e Brouwer  
c: entrada: r,v (x(1)....x(6)) , saida : a,e,incl,w,om,am)
c: se kk=1 : converto (acho) tudo. Se kk=0, entao so acho Excent. e Inclin. 
c: ami=G*(Mj+msat), msat=0. ou ami=Amilu=G*(Mj+dml),etc.  depende do corpo    
c: extrv=produto escalar r.v
c: lei=1 se excent. ou inclinacao sao nulos ou quase nulos e lei entao
c: pode ser usado para controlar a impressao de dados de saida.Nlei=0 se
c: nenhum problema de singularidade.

      do 70 i=1,N

      lei=0
      extrv=x(i,1)*x(i,4)+x(i,2)*x(i,5)+x(i,3)*x(i,6)
      r2=x(i,4)*x(i,4)+x(i,5)*x(i,5)+x(i,6)*x(i,6)
      v2=x(i,1)*x(i,1)+x(i,2)*x(i,2)+x(i,3)*x(i,3)
  20  r=dsqrt(r2)
      v=dsqrt(v2)
      amido=ami(i)+ami(i)
      a(i)=ami(i)/(amido/r-v2)
c: achar angulo teta entre r e v (Fitzpatrick p.71)
c: teta serve tambem p/ definir quadrante de u=anom. excentrica (ver p.71)
c: O processo de Fitzpatrick, para achar u e excent. esta' bom (abaixo),mas,
c: prefiro usar o processo de Brouwer para achar excen/ e u (Brouwer p.48)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c: esta parte (acha u , e ) esta ok,bom mas prefiro Brouwer
c:       costeta=extrv/(r*v)
c:       sinteta2=1.d0-costeta*costeta
c:       efitz=dsqrt(1.d0-r*v2*(2.d0-r*v2/ami)*sinteta2/ami)
c:       afitz=r/(2.d0-r*v2/ami)
c:       cosu=(a-r)/(a*e)
c:       write(*,*)'a,r,e,cosu',a,r,e,cosu
c:       pause
c:       ufitz=dacos(cosu)
c:       if(costeta.lt.0.d0) ufitz=-ufitz
c:       u=ufitz
c:       e=efitz
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c:              abaixo , Brouwer :   (parece ser melhor)

      esinu=extrv/dsqrt(ami(i)*a(i))
      ecosu=r*v2/ami(i)-1.d0
      eb=dsqrt(esinu*esinu+ecosu*ecosu)
      e(i)=eb
c      open(85,file='xyzprob',status='unknown')        
      if(kk.eq.0) go to 30
      if(e(i).le.1.d-11. or. ecosu.eq.0.d0) then  
                                        w(i)=9.d10
                                        f(i)=9.d10
                                        am(i)=9.d10
                                        lei=1           
c                                        write(85,*)'e=0 ou quase'
c                                        write(*,*)'e=0 ou quase'
c                                        pause
                                        go to 30
                                        end if
      ub=datan(esinu/ecosu)
c: em principio,parece ser esta (abaixo) a melhor forma de definir o quadrante

c*****************************************************************************
C*******MODIFIQUEI A DEFINICAO DO QUADRANTE, FAZENDO DE 0-360 GRAUS***********
C*****************************************************************************
      
      if(ub.ge.0.d0 . and . ecosu.lt.0.d0) ub=ub+pi
      if(ub.lt.0.d0 . and . ecosu.lt.0.d0) ub=ub+pi
      if(ub.lt.0.d0 . and . ecosu.gt.0.d0) ub=ub+pi2          
c:    write(*,*)'ub,ufitz,eb,efitz',ub,ufitz,eb,efitz
      u(i)=ub
      rsinf=dsqrt(1.d0-e(i)*e(i))*dsin(u(i))
      rcosf=dcos(u(i))-e(i)
      if(dabs(rcosf).gt.1.d-11) go to 31
      if(rsinf.gt.0.d0) f(i)=pi05
      if(rsinf.lt.0.d0) f(i)=pi15
      go to 29
  31  continue
      f(i)=datan(rsinf/rcosf)
      if(f(i).ge.0.d0 . and . rcosf.lt.0.d0) f(i)=f(i)+pi
      if(f(i).lt.0.d0 . and . rcosf.lt.0.d0) f(i)=f(i)+pi
      if(f(i).lt.0.d0 . and . rcosf.gt.0.d0) f(i)=f(i)+pi2      
C      sinf=dsin(f)
C      cosf=dcos(f)
C      ff=datan(sinf/cosf)
   29 continue 



      am(i)=u(i)-e(i)*dsin(u(i))
      


  30  continue
c: determinacao da inclinacao
      h=dsqrt(ami(i)*a(i)*(1.d0-e(i)*e(i)))
      cosi=(x(i,4)*x(i,2)-x(i,5)*x(i,1))/h
      if(kk.eq.0) then
                      dincl(i)=cosi
                      return
      end if   

C     Daqui para baixo, modifiquei praticamente toda a subrotina por
c     causa do caso PLANO.      
      difcosi=cosi-1.d0
ccccccccccccccccccc
      IF (difcosi.gt.0) THEN
c     Se o cosi>>>1, entao ha' problemas serios. Testes mostraram
c     que mesmo nos casos corretos, APESAR DE I=0, difcosi pode atingir valores 
C     da ordem de 10d-16, talvez por "imprecisao" do fortran (??). No entanto, aqui, 
c     apesar de absurdo, aproximo para i=0 QUANDO:      
                          IF((dabs(difcosi)*1.0e+16).lt.9 ) then                      
                                       dincl(i)=0.0d0
                                       om(i)=0.0d0
                                       lei=5
                          ELSE                   
C     se difcosi>9.d-16, ai considero cosi>>1:                          
c                          WRITE(*,*)'problemas serios com cosi=',cosi
c                          WRITE(85,*)'problemas serios com cosi=',cosi                          
c                          pause
                          end if
      
c      write(85,*)'cosi>1',cosi,difcosi
      else                                              


                          DINCL(i)=DACOS(COSI)      



c     Para o caso plano:
C     Do mesmo modo: Se cosi<1: alguns testes mostram que qdo cosi<1, mas bem
c     proximo de 1, a diferenca  difcosi=cosi-1=-1.0e-16,-2.0e-16,... . 
c     Assim, quando:
C                  IF( (dabs(difcosi)*1.0e+16).gt.1 ) then       
c  aproximo para i=0 nestes casos tambem.      
c: pequenos erros (ao fazer I=0) nao afetam nada, pois este erro nao
c: sera carregado adiante, e' apenas saida de dado.
C                      lei=6
C                      dincl(i)=0.0d0
C                      om(i)=0.d0
C                  END IF

      end if
cccccccccccccc
      if(e(i).le.5.d-11) then
      w(i)=0.0D0
      lei=7
      return 
      end if
cccccccccccccc
      if(dincl(i).eq.0.0d0) then
c     como i=0, w agora e' a longitude do periastro
      
      if(DABS(X(i,4)).gt.1.D-11) go to 34 
c      WRITE(*,*)'PROBLEMAS COM X(i,4)=',X(i,4)
c      WRITE(85,*)'PROBLEMAS COM X(i,4)=',X(i,4)                               
      lei=8          
      w(i)=88
      return
  34  continue                                   

                      wtif=datan(x(i,5)/x(i,4))

                    if(wtif.ge.0.d0. and . x(i,4).lt.0.d0) wtif=wtif+pi
                    if(wtif.lt.0.d0. and . x(i,4).lt.0.d0) wtif=wtif+pi
                    if(wtif.lt.0.d0. and . x(i,4).gt.0.d0) wtif=wtif+pi2
                      
                      ww=wtif-f(i)               
c     A reducao abaixo e' para os casos onde wtif avanca em relacao a f (ou 
c     vice-versa.                      

C 1)    Reduzindo w aos "senos e cossenos" abaixo, posso evitar problemas nos 
c     casos onde |wtif|<|f| e -360<w<-270 (por exemplo, w=-340 graus).
c     Por exemplo, apos a 1a. volta de wtif: se w=wtif-f=10-350=-340 graus.
c     Este w=-340 e', na verdade, igual a +20, pois o wtif ja' deu uma volta.
c     Para resolver isto, basta fazer: 
                      sinw=dsin(ww)
                      cosw=dcos(ww)
                      w(i)=datan(sinw/cosw)                       
c     Ou seja, sinw=sin(-340)=sin(+20) e cosw=cos(-340)=cos(+20), e portanto
c     o w=arctan(sin20/cos20)=+20>0, como queriamos.

C  2)   Outro problema que esta reducao conserta: se wtif=352 e f=355, w=-2. Se
c     f agora avanca 1 volta �frente de wtif (p.ex, wtif vai p/ 354 e f vai
c     para +4, teremos: w=354-4=350. Com a reducao: sin(350)=-sin(10)<0,
c     cos(350)=cos(10)>0, ou seja, w=arctan(-sin10/cos10)<0=-10 graus, 
c     como queriamos.  

c  3)   Se w<0, mas agora, 0<|w|<90, nao ha' problema, o angulo continua a 
c     sendo w<0, com essa reducao.
                      

c: aqui, em que I=0 ,este w de saida representa  w+om,ie, w e'long. pericentro
c: tambem om nao esta definido (nao existe,mas tomarei zero)
      endif


c      if(x(4).gt.0.0d0.and.x(5).gt.0.0d0) iqwtif=1
c      if(x(4).lt.0.0d0.and.x(5).gt.0.0d0) iqwtif=2
c      if(x(4).lt.0.0d0.and.x(5).lt.0.0d0) iqwtif=3      
c      if(x(4).gt.0.0d0.and.x(5).lt.0.0d0) iqwtif=4
      
c      if(rcosf.gt.0.0d0.and.rsinf.gt.0.0d0) iqf=11
c      if(rcosf.lt.0.0d0.and.rsinf.gt.0.0d0) iqf=22
c      if(rcosf.lt.0.0d0.and.rsinf.lt.0.0d0) iqf=33      
c      if(rcosf.gt.0.0d0.and.rsinf.lt.0.0d0) iqf=44
          
cccccccccccccccccccccccc

c:Aqui Inclin.nao e'zero nem quase,logo,pelo menos omega=om posso achar.
      sinomi=(x(i,5)*x(i,3)-x(i,6)*x(i,2))                                   
      cosomi=-(x(i,6)*x(i,1)-x(i,4)*x(i,3))
      if(dabs(cosomi).gt.1.d-11) go to 45
      if(sinomi.gt.0.d0) om(i)=pi05
      if(sinomi.lt.0.d0) om(i)=pi15
      go to 41
  45  continue
      om(i)=datan(sinomi/cosomi)
      if(om(i).ge.0.d0. and . cosomi.lt.0.d0) om(i)=om(i)+pi
      if(om(i).lt.0.d0. and . cosomi.lt.0.d0) om(i)=om(i)+pi
      if(om(i).lt.0.d0. and . cosomi.gt.0.d0) om(i)=om(i)+pi2      
  41  continue    
      if(eb.le.5.d-11) return 
C: Aqui excent. nao e' nula, posso achar w tambem 

      coslat=x(i,4)*dcos(om(i))+x(i,5)*dsin(om(i))
      sinlat=(x(i,5)*dcos(om(i))-x(i,4)*dsin(om(i)))*dcos(dincl(i))+
     |x(i,6)*dsin(dincl(i))
      if(dabs(coslat).gt.5.d-11) go to 46
      if(sinlat.gt.0.d0) alat=pi05
      if(sinlat.lt.0.d0) alat=pi15
      go to 47
  46  continue
      alat=datan(sinlat/coslat)
      if(alat.ge.0.d0. and . coslat.lt.0.d0) alat=alat+pi
      if(alat.lt.0.d0. and . coslat.lt.0.d0) alat=alat+pi
      if(alat.lt.0.d0. and . coslat.gt.0.d0) alat=alat+pi2      
   47 continue    
      w(i)=alat-f(i)


70    continue      
c: a formula do sinlat vem do McCuskey (posso deduzir facilmente das eqcoes
c: x,y,z  que relacionam elementos orbitais, idem coslat)
      return
      end






      subroutine orbxyz1(i,a,e,dincl,w,om,am,ami,x,f)
      implicit real *8(a-h,o-z)
      include 'dimensoes.inc'
      COMMON/IN/ikepler,inicio,pi,pi05,pi15,pi2
      dimension a(NMAX),e(NMAX),w(NMAX),om(NMAX),am(NMAX),ami(NMAX),
     |f(NMAX),x(NMAX,6),u(NMAX),dincl(NMAX)

 

c: rotina que passa de elementos orbitais para cartesianos
c: preciso tambem das derivadas no tempo das cartesianas(velocidades)
c: referencia: Fitzpatrick (em velocidades parece que nao esta ok)

c: entrada : elementos orbitais , saida : conversao para cartesianos
c: u=anom.excent. (sera achado via rotina Kepler)
c: am=anom. media
c: om=longitude do nodo (indefinido se I=0,mas tomo zero para que as formulas
c: do caso espacial fiquem validas no caso plano)
c: w=argumento do pericentro (se I=0 este w representara longitude do pericen.)
c: alat=latitude=w+f (representara wtil+f no caso I=0) (wtil=wpi)

c: obtencao de anom. excent. (u) e verdadeira (f) via Kepler
c     (so o corpo i)
      call kepler1(i,f,u,am,e)

c      do 50 i=1,N
      inicio=0
      alat=w(i)+f(i)
      r=a(i)*(1.d0-e(i)*dcos(u(i)))
      rp=dsqrt(ami(i)*a(i))*e(i)*dsin(u(i))/r
      rfp=dsqrt(ami(i)*a(i)*(1.d0-e(i)*e(i)))/r
c: velocidades (x(1),x(2),x(3) no caso tridimensional
c: Nesta rotina se I = zero ou e= 0 nao ha nenhum problema, embora
c: para entrar com angulos eu entre com Wtil (no lugar de w)
c: e em lugar de om=nodo=indefinido tomo=0
      cfw=dcos(alat)
      sfw=dsin(alat)
      coom=dcos(om(i))
      siom=dsin(om(i))
      coi=dcos(dincl(i))
      somcoi=siom*coi
      comcoi=coom*coi
      sini=dsin(dincl(i))
      cosi=dcos(dincl(i))
c: no caso plano basta tomar I=dincl=0 e definir om=0.d0
c: tambem w daqui deve ser interpretado como sendo o wpi=long. pericentro
c: tambem alat deve ser pensado como alat=wpi+f

c: velocidades

      x(i,1)=rp*cfw*coom-rfp*sfw*coom-rp*sfw*somcoi-rfp*cfw*somcoi
      
      x(i,2)=rp*cfw*siom-rfp*sfw*siom+rp*sfw*comcoi+rfp*cfw*comcoi
      
      x(i,3)=rp*sfw*sini+rfp*cfw*sini

c: coordenadas

      x(i,4)=r*(cfw*coom-sfw*cosi*siom)    
     
      x(i,5)=r*(cfw*siom+sfw*cosi*coom)
     
      x(i,6)=r*sini*sfw

    

      return
      end

      
      SUBROUTINE KEPLER1 (j,f,u1,AM,E0)
      IMPLICIT REAL *8(A-H,O-Z) 
      include 'dimensoes.inc'
      COMMON/IN/ikepler,inicio,pi,pi05,pi15,pi2                       
      dimension f(NMAX),fa(NMAX),am(NMAX),e0(NMAX),u1(NMAX)

C      READ(*,*) AM,E0
C: metodo de Newton Raphson  
c:ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c: este bloco eu usava p/ garantir grande precisao na solucao da eq.Kepler
c: apenas p/ converter no inicio os elementos orbitais em x,y,z (qdo tomo kepl=1)
c: voce pode deixar sempre kepl=1, e garante sempre grande precisao, mas fica
c: mais lento.Fica a seu criterio tirar, mas tem que dar o LITERA e PREC)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c     (so o corpo j)
      
      litera=15
      prec=1.d-10
      if(ikepler.eq.1) then
                           litera=35
                           prec=1.d-12
                           end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                           
      U0=AM(j)+E0(j)*DSIN(AM(j))
      DO 10 I=1,litera
      U1(j)=U0-((U0-E0(j)*DSIN(U0)-
     /AM(j)))/(1.D0-E0(j)*DCOS(U0))
      TEST1=DABS(U1(j)-U0)
      DSU1=DSIN(U1(j))
      TEST2=DABS(U1(j)-E0(j)*DSU1-AM(j))
      IF(TEST1.LE.prec.AND.TEST2.LE.prec) GO TO 11
      U0=U1(j)
 10   CONTINUE
c**      WRITE(*,*)' NAO CONVERGENCIA EM KEPLER'
c**      write(*,*)'test2    excent', test2, e0,u1
 11   CONTINUE
      RCF=DCOS(U1(j))-E0(j)
c: rcf=0 ou quase se f=anom. verdad.= proximo de pi/2,pi*3/2.Se rcf<=1.d-11
c: suporemos f =pi/2 ou pi*3/2, mas se isto ocorrer ja no inicio (primeira)
c: chamada de Kepler, posso (opcionalmente) parar a integracao entrando com
c: nova condicao inicial para o Planeta (esperando que nao ocorra o mesmo)
c: (se porem o parametro INICIO e'=0,apenas considero f= pi/2 ou pi*3/2.Esta 
c: aproximacao so' tem problema se estamos usando Kepler para o Planeta,
c: pois no caso do satelite, esta transformacao e' apenas um dado de saida.
      if(dabs(rcf).gt.1.d-11) go to 20
      if(inicio.eq.1)then
      open (10,file='keple',status='unknown')
c      write(10,*)'melhor mudar AM ou E0 inic do Planeta,para que f seja
c     +melhor determinado'
                         stop 
                         end if
      if(dsu1.gt.0.d0) f(j)=pi05
      if(dsu1.lt.0.d0) f(j)=pi15
      return
  20  continue

      RSF=DSQRT(1.D0-E0(j)*E0(j))*DSU1
      FA(j)=DATAN(RSF/RCF)
      F(j)=FA(j)
      
      IF(f(j).ge.0.D0.AND.rcf.LT.0.D0) F(j)=FA(j)+PI
      IF(f(j).lt.0.D0.AND.rcf.lt.0.D0) F(j)=FA(j)+PI
      if(f(j).lt.0.d0.and.rcf.gt.0.d0) F(j)=FA(j)+PI2
      
                                      

c: apos orbxyz chamar Kepler, ele redefine inicio=0
      RETURN
      END 
      subroutine fitzpaorb(x,a,e,dincl,w,alat,am,u,om,ami) ! nao uso, lei nem kk
      

c: x1,x2,x3=velocidades (do corpo em questao)  x4,x5,x6= coordenadas
 
      implicit real *8(a-h,o-z)
      dimension x(6)    ! no ifort uso x(1)     

      COMMON/IN/ikepler,inicio,pi,pi05,pi15,pi2
      common/ie/iecentr, ipular       ! 15/11/2014:  iecentr p/ parar, ipular retornar sem  transformar
            
c: referencias : Fitzpatrick e Brouwer  
c: entrada: r,v (x(1)....x(6)) , saida : (a,e,incl,w,om,am)

c: ami=G*(Mj+msat), msat=0. ou ami=Amilu=G*(Mj+dml),etc.  depende do corpo    
c: extrv=produto escalar r.v


      lei=0     
      iecentr=0  ! em 06/01/2002 : pus para controlar caso onde e > =1, para nao travar.
      ipular=0   ! em 14/nov/2014
      
      extrv=x(1)*x(4)+x(2)*x(5)+x(3)*x(6)
      r2=x(4)*x(4)+x(5)*x(5)+x(6)*x(6)
      v2=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
  20  r=dsqrt(r2)
      v=dsqrt(v2)
      amido=ami+ami
     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c: 07/01/2002: testar se excentricidade e' maior que 1. Se sim, retorno, sem calcular w,om, am.
c:            Aqui, se energia e' maior ou igual a zero, a orbita nao e' mais eliptica
c:            Neste caso so calculo a inclinacao de forma separada, aqui 
      energia=0.5d0*v2-ami/r   

      if(energia.ge.0.d0) then
                          iecentr=555   
                          e=1.d0  ! Eventualmente pode ser maior que 1, mas 1 ja esta ok.
                          lei=1
                          write(*,*)'exc > =1,parar,energia=', energia
      write(85,*) x(1),x(2),x(3),x(4),x(5),x(6)       
                          stop
                          end if
cc** daqui ate o proximo cc** cc, posso ignorar, so ativo se quero saber a inclinacao deste caso                              
                              hx=x(5)*x(3)-x(2)*x(6) ! componentes x,y,z do momento angular h
                              hy=-x(4)*x(3)+x(1)*x(6)! Danby p.201 ou Fitzpatrick
                              hz=x(4)*x(2)-x(1)*x(5)
                              hmodulo=dsqrt(hx*hx+hy*hy+hz*hz)
                              cosi=hz/hmodulo  ! cosseno da inclinacao (sempre entre 0,180) 
                     abcosi=dabs(cosi)
                     dif=dabs(abcosi-1.d0)
                     if(dabs(cosi).gt.1.d0 .or. dif.le.1.d-13) then
                              ipular=666
                              lei=1
                              write(85,*) 'Problema cosi', cosi

                              return ! volto em OUTPUT, sem calcular elementos orbitais
                                                                 end if
                              dincl=dacos(cosi) ! nunca ha problema de quadrante para dincl                 

c:                              w=7.d9  ! nao me interessam w,om, am neste caso (e >=0).
c:                              om=7.d9 ! dou valores absurdos so para identificar que e> = 0
c:                              am=7.d9 ! Note que dincl e excentr., sao verdadeiros  
cc** cc  ate aqui posso ignorar, se nao me interessa a inclinacao deste caso particular                              
                          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      
c:      a=ami/(amido/r-v2) ! do Brouwer
c: achar angulo teta entre r e v (Fitzpatrick p.71)
c: teta serve tambem p/ definir quadrante de u=anom. excentrica (ver p.71)
c: O processo de Fitzpatrick, para achar u e excent. esta' bom (abaixo),mas,
c: Na rotina inicial XYZORB uso o processo de Brouwer para achar excen/ e u (Brouwer p.48)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c: esta parte (acha excent e  u=anomalia excentrica (ufitz=E do Brouwer)
c: conforme Fitzpatrick,p.71: facil ver: teta=ang. entre vetor r e vetor v, entao teta so pode estar
c:  em [0,pi], pois v=vetor velocidade é sempre tangente aa orbita. Ainda:
 
c: 1-qdo costeta <0 : teta estaria no quadrante 2 ou 3, mas devido acima, teta automaticamente esta no quadrante 2
c: 2-qdo costeta >0 ; teta estaria no quadrante 1 ou 4, mas devido acima, teta automaticamente esta no quadrante 1
c: Na verdade isso e'apenas para saber o quadrante de teta (por ora isso nao sera preciso ) Importante eh o abaixo:

       costeta=extrv/(r*v)
       sinteta2=1.d0-costeta*costeta
       efitz=dsqrt(1.d0-r*v2*(2.d0-r*v2/ami)*sinteta2/ami)
       afitz=r/(2.d0-r*v2/ami)
       cosu=(afitz-r)/(afitz*efitz)
       ufitz0=dacos(cosu)
       
       ecosf=afitz*(1.d0-efitz*efitz)/r-1.d0     ! da eq 44 do Brouwer p.21

c: 31-0502018: Na minha opiniao, se existisse uma elipse osculadora, ecosf deveria variar de (-e,+e). Mostro
c:             que o Anthe nao permite isso, ou seja nao existe nenhum r, com a,e,fixados tal que f= varia de 0,2PI
c:             Assim novamente mostro que este f so pode librar, assim como u e também am. Por esta razao, o w varia com
c:             uma perturbacao de curto periodo igual a am, ao mesmo tempo que am libra em torno de zero.

c: Preciso agora definir o quadrante de ufitz. Para isso: sabemos trigonometria que:

c: se cosu >0 , necessariamente ufitz estara no quadrante 1 ou 4. 
c: se cosu <0, entao necessariamente  ufitz estara no quadrante 2 ou 3

c: Agora ajunto a informacao : veja abaixo (3):

c: 3-Se costeta <0 entao r =modulo=distancia do sat ao foco, esta diminuindo, pois: dr/dt={modulo da velocidade do sat}*costeta
c: esta equacao dr/dt eu deduzi rigorosamente, ver verso pag.71 do Fizpatrick. 
c: entao: se r=modulo da distancia ao foco esta decaindo, entao o sat esta no quadrante 3 ou 4 (so fazer o desenho da orbita eliptica).
c: Analogamente se costeta >0 entao sat esta se afastando, logo sat estara no quadrante 1 ou 2

c: Resumindo tudo: 
c: costeta < 0 : sat em [3,4]   - sat aproximando do foco 
c: costeta > 0 : sat em [1,2]   - sat afastando do foco
c: cosu < 0 : sat em [2,3]
c: cosu > 0 : sat em [1,4]  
c: fazendo a combinacao de todos os casos resultam 4  e uFitz é possivel de ser achado de forma única, facil, so fazer

      if(costeta .lt. 0.d0 .and. cosu.ge. 0.d0) uFitz=-uFitz0   ! aqui uFitz vai p 4-quad, pois uFitz0 esta Quad 1
      if(costeta.lt. 0.d0 .and. cosu .le. 0.d0) uFitz=pi2-uFitz0  ! uFitz vai p/3-quadrante, pois uFitz0=esta Quad 2
      if(costeta .gt. 0.d0 . and. cosu .ge.0.d0) uFitz=uFitz0 ! aqui nao precisa nada, uFitz0 esta Quad 1
      if(costeta .gt. 0.d0 . and. cosu .le.0.d0) uFitz=uFitz0 ! aqui nao precisa fazer nada, uFtz0 esta Quad 2


       u=ufitz  ! anomalia excentrica 
       e=efitz
       a=afitz
c: achar anomalia média 
       am=u-e*dsin(u)

      rsinf=dsqrt(1.d0-e*e)*dsin(u)
      rcosf=dcos(u)-e
      if(dabs(rcosf).gt.1.d-13) go to 31
      if(rsinf.gt.0.d0) f=pi05
      if(rsinf.lt.0.d0) f=pi15
      go to 29
  31  continue
      f=datan(rsinf/rcosf)
      if(f.ge.0.d0 . and . rcosf.lt.0.d0) f=f+pi
      if(f.lt.0.d0 . and . rcosf.lt.0.d0) f=f+pi
   29 continue 
      amantigo=u-e*dsin(u)  ! esta eh a anomalia que sempre usei, ate 25-05-2018 (agora uso uFitz) 

  40  continue
c:aqui Inclin.nao e'zero nem quase,logo, omega=om posso achar.
      sinomi=(x(5)*x(3)-x(6)*x(2))                                   
      cosomi=-(x(6)*x(1)-x(4)*x(3))
      if(dabs(cosomi).gt.1.d-12) go to 45
      if(sinomi.gt.0.d0) om=pi05
      if(sinomi.lt.0.d0) om=pi15
      go to 41
  45  continue
      om=datan(sinomi/cosomi)
      if(om.ge.0.d0. and . cosomi.lt.0.d0) om=om+pi
      if(om.lt.0.d0. and . cosomi.lt.0.d0) om=om+pi
  41  continue    
       
c: Aqui excent. nao e' nula, posso achar w tambem 

      coslat=x(4)*dcos(om)+x(5)*dsin(om)
      sinlat=(x(5)*dcos(om)-x(4)*dsin(om))*dcos(dincl)+x(6)*dsin(dincl)
      if(dabs(coslat).gt.1.d-13) go to 46
      if(sinlat.gt.0.d0) alat=pi05
      if(sinlat.lt.0.d0) alat=pi15
      go to 47
  46  continue
      alat=datan(sinlat/coslat)
      if(alat.ge.0.d0. and . coslat.lt.0.d0) alat=alat+pi
      if(alat.lt.0.d0. and . coslat.lt.0.d0) alat=alat+pi
   47 continue    
      w=alat-f
c: a formula do sinlat vem do McCuskey (posso deduzir facilmente das eqcoes
c: x,y,z  que relacionam elementos orbitais, idem coslat)
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      SUBROUTINE RA15(TM,X,V,TF,LL,XL,NV,NCLASS)

C  Integrador RADAU de Everhart - Version de orden 15
C  Integra ecuaciones de segundo orden:
C    y' = F(y,t)    es NCLASS = 1
C    y" = F(y',y,t) es NCLASS = 2  
C    y" = F(y,t)    es NCLASS = -2
C  TF es t(final) - t(initial). Debe ser negativo para 
C    integracion hacia atras
C  NV es el numero de ecuaciones diferenciales simultaneas
C  LL controla el tamanio de paso, monitoreando el error en los
C    terminos en base a la tolerancia SS = 10**(-LL). Un valor
C    tipico para LL esta entre 6 y 12
C  XL es el tamanio de paso constante si LL < 0. Si no, es el
C    tamanio de paso inicial
C  TM, X y V son tiempo, coordenadas y velocidades iniciales en
C    la entrada. En la salida son devueltas como tiempo, 
C    coordenadas y velocidades finales.
      
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimensoes.inc'
      DIMENSION X(6*NMAX),V(6*NMAX),F1(6*NMAX),FJ(6*NMAX),
     >     C(21),D(21),R(21),Y(6*NMAX),Z(6*NMAX),
     >     B(7,6*NMAX),G(7,6*NMAX),E(7,6*NMAX),BD(7,6*NMAX),H(8),
     >     W(7),U(7),NW(8)
      LOGICAL NSF,NPER,NPQ,NCL,NES
      DATA NW/0,0,1,3,6,10,15,21/
      DATA ZERO, HALF, ONE, SR/0.0D0, 0.5D0, 1.0D0, 1.4D0/
      DATA H/         0.D0, .05626256053692215D0, .18024069173689236D0,
     >.35262471711316964D0, .54715362633055538D0, .73421017721541053D0,
     >.88532094683909577D0, .97752061356128750D0/
C  Estos valores de H son los espaciadores de Gauss - Radau,
C  escalados en el intervalo [0,1]

      NPER=.FALSE.
      NSF=.FALSE.
      NCL=NCLASS.EQ.1
      NPQ=NCLASS.LT.2
      NES=LL.LT.0
C  NCLASS =  1   ==>   NCL = .TRUE.    NPQ = .TRUE.
C  NCLASS = -2   ==>   NCL = .FALSE.   NPQ = .TRUE.
C  NCLASS =  2   ==>   NCL = .FALSE.   NPQ = .FALSE.
C  NPER es .TRUE. solo en la ultima secuencia de integracion.
C  NSF es .FALSE. solo en la secuencia inicial
C  NES es .TRUE. solo si LL es negativo. Entonces el tamanio de
C  paso es XL.
      
      IF (TF.LT.ZERO) THEN
       DIR=-ONE
      ELSE
       DIR=ONE
      END IF
      
      XL=DIR*DABS(XL)
      PW=1./9.
      
      DO N=2,8
       WW=N+N*N
       IF (NCL) WW=N
       W(N-1)=ONE/WW
       WW=N
       U(N-1)=ONE/WW
      END DO
      
      DO K=1,NV
       IF (NCL) V(K)=ZERO
       DO L=1,7
	BD(L,K)=ZERO
	B(L,K)=ZERO
       END DO
      END DO
      
      W1=HALF
      IF (NCL) W1=ONE
      C(1)=-H(2)
      D(1)=H(2)
      R(1)=ONE/(H(3)-H(2))
      LA=1
      LC=1
      
      DO K=3,7
       LB=LA
       LA=LC+1
       LC=NW(K+1)
       C(LA)=-H(K)*C(LB)
       C(LC)=C(LA-1)-H(K)
       D(LA)=H(2)*D(LB)
       D(LC)=-C(LC)
       R(LA)=ONE/(H(K+1)-H(2))
       R(LC)=ONE/(H(K+1)-H(K))
       
       IF (K.NE.3) THEN
	DO L=4,K
	 LD=LA+L-3
	 LE=LB+L-4
	 C(LD)=C(LE)-H(K)*C(LE+1)
	 D(LD)=D(LE)+H(L-1)*D(LE+1)
	 R(LD)=ONE/(H(K+1)-H(L-1))
	END DO
       END IF
      
      END DO
      
      SS=10.**(-LL)
C  Las instrucciones anteriores son calculadas solo una vez
C  en una integracion para inicializar las constantes

C  A continuacion se inicializa el tamanio de paso inicial TP
      IF (NES) THEN
       TP=XL
      ELSE IF (XL.NE.ZERO) THEN
       TP=XL
      ELSE
       TP=0.1D0*DIR
      END IF
      
      IF (TP/TF.GT.HALF) TP=HALF*TF
      NCOUNT=0

C  La linea 1000 es el comienzo de la primera secuencia.
C  NS es el numero de secuencia. 
C  NF es el numero de llamados a la subrrutina FORCE. 
C  NI es el numero de iteraciones en cada secuencia
1000  NS=0
      NF=0
      NI=6
      TM=ZERO
      CALL FORCE (TM, X, V, F1)
      NF=NF+1

C  La linea 2000 es el comienzo de cada secuencia despues de
C  la primera
2000  DO K=1,NV
       G(1,K)=B(1,K)+D( 1)*B(2,K)+D(2)*B(3,K)+
     >   D(4)*B(4,K)+D( 7)*B(5,K)+D(11)*B(6,K)+D(16)*B(7,K)
       G(2,K)=             B(2,K)+D(3)*B(3,K)+
     >   D(5)*B(4,K)+D( 8)*B(5,K)+D(12)*B(6,K)+D(17)*B(7,K)
       G(3,K)=                         B(3,K)+
     >   D(6)*B(4,K)+D( 9)*B(5,K)+D(13)*B(6,K)+D(18)*B(7,K)
       G(4,K)=B(4,K)+D(10)*B(5,K)+D(14)*B(6,K)+D(19)*B(7,K)
       G(5,K)=             B(5,K)+D(15)*B(6,K)+D(20)*B(7,K)
       G(6,K)=                          B(6,K)+D(21)*B(7,K)
       G(7,K)=                                       B(7,K)
      END DO
      
      T=TP
      T2=T*T
      IF (NCL) T2=T
      TVAL=DABS(T)

C  Comienzo de las iteraciones para cada secuencia. Realiza 6
C  iteraciones en la primera secuencia (NI=6)
      DO M=1,NI
       DO J=2,8
	JD=J-1
	JDM=J-2
	S=H(J)
	Q=S
	IF (NCL) Q=ONE
	
C  Aqui calcula los predictores de posicion y velocidad para
C  cada subsecuencia        
	DO K=1,NV
	 A=W(3)*B(3,K)+S*(W(4)*B(4,K)+S*(W(5)*B(5,K)+S*(W(6)*
     >     B(6,K)+S*W(7)*B(7,K))))
	 Y(K)=X(K)+Q*(T*V(K)+T2*S*(F1(K)*W1+S*(W(1)*B(1,K)+
     >        S*(W(2)*B(2,K)+S*A))))
	 
	 IF (.NOT.NPQ) THEN
C         Esta parte solo es calculada cuando NCLASS = 2          
	  A=U(3)*B(3,K)+S*(U(4)*B(4,K)+S*(U(5)*B(5,K)+
     >      S*(U(6)*B(6,K)+S*U(7)*B(7,K))))
	  Z(K)=V(K)+S*T*(F1(K)+S*(U(1)*B(1,K)+S*(U(2)*B(2,K)+
     >      S*A)))
	 END IF
	
	END DO
	
C  Calcula la fuerza FJ al final de cada subsecuencia
	CALL FORCE(TM+S*T, Y, Z, FJ)
	NF=NF+1
	
C  Calcula los nuevos valores de G y de B para la fuerza FJ        
	DO K=1,NV
	 TEMP=G(JD,K)
	 GK=(FJ(K)-F1(K))/S
	 
	 GO TO (102,102,103,104,105,106,107,108),J
 102        G(1,K)=GK
	    GO TO 100
 103        G(2,K)=(GK-G(1,K))*R(1)
	    GO TO 100
 104        G(3,K)=((GK-G(1,K))*R(2)-G(2,K))*R(3)
	    GO TO 100
 105        G(4,K)=(((GK-G(1,K))*R(4)-G(2,K))*R(5)-G(3,K))*
     >             R(6)
	    GO TO 100
 106        G(5,K)=((((GK-G(1,K))*R(7)-G(2,K))*R(8)-G(3,K))*
     >             R(9)-G(4,K))*R(10)
	    GO TO 100
 107        G(6,K)=(((((GK-G(1,K))*R(11)-G(2,K))*R(12)-
     >             G(3,K))*R(13)-G(4,K))*R(14)-G(5,K))*R(15)
	    GO TO 100
 108        G(7,K)=((((((GK-G(1,K))*R(16)-G(2,K))*R(17)-
     >             G(3,K))*R(18)-G(4,K))*R(19)-G(5,K))*R(20)-
     >             G(6,K))*R(21)
 100        CONTINUE
	 
	 TEMP=G(JD,K)-TEMP
	 B(JD,K)=B(JD,K)+TEMP
	 
	 GO TO (200,200,203,204,205,206,207,208),J
 203        B(1,K)=B(1,K)+C(1)*TEMP
	    GO TO 200
 204        B(1,K)=B(1,K)+C(2)*TEMP
	    B(2,K)=B(2,K)+C(3)*TEMP
	    GO TO 200
 205        B(1,K)=B(1,K)+C(4)*TEMP
	    B(2,K)=B(2,K)+C(5)*TEMP
	    B(3,K)=B(3,K)+C(6)*TEMP
	    GO TO 200
 206        B(1,K)=B(1,K)+C(7)*TEMP
	    B(2,K)=B(2,K)+C(8)*TEMP
	    B(3,K)=B(3,K)+C(9)*TEMP
	    B(4,K)=B(4,K)+C(10)*TEMP
	    GO TO 200
 207        B(1,K)=B(1,K)+C(11)*TEMP
	    B(2,K)=B(2,K)+C(12)*TEMP
	    B(3,K)=B(3,K)+C(13)*TEMP
	    B(4,K)=B(4,K)+C(14)*TEMP
	    B(5,K)=B(5,K)+C(15)*TEMP
	    GO TO 200
 208        B(1,K)=B(1,K)+C(16)*TEMP
	    B(2,K)=B(2,K)+C(17)*TEMP
	    B(3,K)=B(3,K)+C(18)*TEMP
	    B(4,K)=B(4,K)+C(19)*TEMP
	    B(5,K)=B(5,K)+C(20)*TEMP
	    B(6,K)=B(6,K)+C(21)*TEMP
 200        CONTINUE
	
	END DO
       END DO
       
C  Final de la secuencia. Calculo del control (HV) del tamanio
C  de paso       
       IF (.NOT.NES.OR.M.GE.NI) THEN
	HV=ZERO
	DO K=1,NV
	 HV=DMAX1(HV,DABS(B(7,K)))
	END DO
	TVAL1=TVAL*TVAL*TVAL*TVAL*TVAL*TVAL*TVAL
	HV=HV*W(7)/TVAL1
       END IF
      
      END DO
      
C  Calcula el nuevo tamanio de paso (TP) y lo compara con HV.
C  Si es menor que lo previsto, reinicia la secuencia anterior
C  usando un paso 0.8 veces menor
      IF (.NOT.NSF) THEN
       
       IF (NES) THEN
	TP=XL
       ELSE
	TP=(SS**PW)/(HV**PW)*DIR
	IF (TP/T.LE.ONE) THEN
	 TP=.8D0*TP
	 NCOUNT=NCOUNT+1
	 IF (NCOUNT.GT.10) RETURN
	 GO TO 1000
	END IF
       END IF
       
       NSF=.TRUE.
      END IF

C  Coordenadas y velocidades y tiempo al final de la secuencia
      DO K=1,NV
       X(K)=X(K)+V(K)*T+T2*(F1(K)*W1+B(1,K)*W(1)+B(2,K)*W(2)+
     >      B(3,K)*W(3)+B(4,K)*W(4)+B(5,K)*W(5)+B(6,K)*W(6)+
     >      B(7,K)*W(7))
       IF (.NOT.NCL) THEN
	V(K)=V(K)+T*(F1(K)+B(1,K)*U(1)+B(2,K)*U(2)+B(3,K)*U(3)
     >       +B(4,K)*U(4)+B(5,K)*U(5)+B(6,K)*U(6)+B(7,K)*U(7))
       END IF
      END DO
      
      TM=TM+T
      NS=NS+1
      
C  Si termino retorna. Caso contrario controla el tamanio de
C  la siguiente secuencia y la ajusta para cubrir exactamente
C  el intervalo de total de integracion (TF)
      IF (NPER) RETURN
      CALL FORCE (TM, X, V, F1)
      NF=NF+1
      
      IF (NES) THEN
       TP=XL
      ELSE
       TP=DIR*(SS**PW)/(HV**PW)
       IF (TP/T.GT.SR) TP=T*SR
      END IF
      
      IF (DIR*(TM+TP).GE.DIR*TF-1.D-8) THEN
       TP=TF-TM
       NPER=.TRUE.
      END IF
  
C  Ahora predice los nuevos valores de B para utilizarlos en
C  la proxima secuencia      
      Q=TP/T
      
      DO K=1,NV
       
       IF (NS.NE.1) THEN
	DO J=1,7
	 BD(J,K)=B(J,K)-E(J,K)  
	END DO
       END IF
       
       E(1,K)=     Q*(B(1,K)+ 2.D0*B(2,K)+ 3.D0*B(3,K)+
     >           4.D0*B(4,K)+ 5.D0*B(5,K)+ 6.D0*B(6,K)+ 
     >           7.D0*B(7,K))
       E(2,K)=               Q**2*(B(2,K)+ 3.D0*B(3,K)+
     >           6.D0*B(4,K)+10.D0*B(5,K)+15.D0*B(6,K)+
     >          21.D0*B(7,K))
       E(3,K)=                            Q**3*(B(3,K)+
     >           4.D0*B(4,K)+10.D0*B(5,K)+20.D0*B(6,K)+
     >          35.D0*B(7,K))
       E(4,K)=  Q**4*(B(4,K)+ 5.D0*B(5,K)+15.D0*B(6,K)+
     >          35.D0*B(7,K))
       E(5,K)=               Q**5*(B(5,K)+ 6.D0*B(6,K)+
     >          21.D0*B(7,K))
       E(6,K)=                            Q**6*(B(6,K)+
     >           7.D0*B(7,K))
       E(7,K)=   Q**7*B(7,K)
       
       DO L=1,7
	B(L,K)=E(L,K)+BD(L,K)
       END DO
      
      END DO

C  A partir de la segunda secuencia, solo realiza dos 
C  iteraciones por cada secuencia
      NI=2
      GO TO 2000

      END
