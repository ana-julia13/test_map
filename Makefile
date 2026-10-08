# Uso:
#   make CASO=naiad        -> compila mapa_naiad  (usa condicoes/naiad.f)
#   make CASO=victor       -> compila mapa_victor (usa condicoes/victor.f)
#   make run CASO=naiad    -> compila e roda dentro de rodadas/naiad/
#
# O codigo base (base/mapa_base.f) e o mesmo para todos os casos.

CASO ?= naiad
FC ?= gfortran
FFLAGS ?= -O2 -ffixed-line-length-none -std=legacy -Ibase

EXE = mapa_$(CASO)

$(EXE): base/mapa_base.f base/dimensoes.inc condicoes/$(CASO).f
	$(FC) $(FFLAGS) -o $@ base/mapa_base.f condicoes/$(CASO).f

run: $(EXE)
	mkdir -p rodadas/$(CASO)
	cd rodadas/$(CASO) && ../../$(EXE)

clean:
	rm -f mapa_*

.PHONY: run clean
