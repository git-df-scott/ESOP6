CC      ?= gcc
CFLAGS  ?= -O3 -march=native -fopenmp -Wall
LDLIBS  := -lm
BIN     := bin
PROGS   := $(BIN)/search $(BIN)/caseA2 $(BIN)/caseA $(BIN)/audit

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

clean:
	rm -rf $(BIN)

.PHONY: all validate control clean
