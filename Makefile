CC      ?= gcc
CFLAGS  ?= -O3 -march=native -fopenmp -Wall
LDLIBS  := -lm
BIN     := bin
PROGS   := $(BIN)/search $(BIN)/caseA2 $(BIN)/caseA3 $(BIN)/caseA $(BIN)/audit $(BIN)/frontier_audit
PROTOS  := $(BIN)/routed_join $(BIN)/mmap_bloom

all: $(PROGS)

$(BIN):
	mkdir -p $(BIN)

$(BIN)/audit: src/audit.c | $(BIN)
	$(CC) -O2 -o $@ $< $(LDLIBS)

$(BIN)/%: src/%.c | $(BIN)
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

# k=5 mode must rediscover Lander-Parkin: 27^5+84^5+110^5+133^5 = 144^5
validate: $(BIN)/search $(BIN)/audit
	@echo "== fifth-power validation (expect f=144 solution) =="
	@./$(BIN)/search 5 2 150
	@echo "== candidate enumeration audit (expect all OK) =="
	@./$(BIN)/audit

# Control run inside the known-clear region: expect 124 candidates, 0 found
control: $(BIN)/caseA2
	@echo "== control run 700k-730k (expect candidates=124 found=0) =="
	@./$(BIN)/caseA2 700000 730000

# Bucketed variant must evaluate exactly the same oracle nodes as monolithic
equiv: $(BIN)/caseA3
	@echo "== equivalence: monolithic vs 8 buckets (evaluated counts must match) =="
	@./$(BIN)/caseA3 700000 730000 16 -b 1 2>&1 | grep -E 'j2_|done:'
	@./$(BIN)/caseA3 700000 730000 16 -b 8 2>&1 | grep -E 'j2_|done:'

differential: $(BIN)/caseA2 $(BIN)/caseA3
	@python3 tests/differential.py

prototypes: $(PROTOS)
	@./$(BIN)/routed_join 2000 8 20000
	@./$(BIN)/mmap_bloom 5000 12 20000

frontier-audit: $(BIN)/frontier_audit
	@./$(BIN)/frontier_audit

direct-checks:
	@python3 tools/astra_geometry.py
	@python3 tests/direct_verifiers.py
	@python3 tools/replay_surface_divisor.py results/astra_direct_2026_09_04

clean:
	rm -rf $(BIN)

.PHONY: all validate control equiv differential prototypes frontier-audit direct-checks clean
