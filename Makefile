# Tambem da para compilar na mao:  gfortran -o mapa_naiad naiad.f
#   make CASO=naiad        -> compila mapa_naiad
#   make run CASO=naiad    -> compila e roda dentro de rodadas/naiad/

CASO ?= naiad
FC ?= gfortran
FFLAGS ?= -O2

EXE = mapa_$(CASO)

$(EXE): mapa_base.f dimensoes.inc $(CASO).f
	$(FC) $(FFLAGS) -o $@ $(CASO).f

run: $(EXE)
	mkdir -p rodadas/$(CASO)
	cd rodadas/$(CASO) && ../../$(EXE)

clean:
	rm -f mapa_naiad mapa_victor mapa_map_i mapa_mapa

.PHONY: run clean
