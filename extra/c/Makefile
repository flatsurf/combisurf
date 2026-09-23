# Benchmark of combisurf/src/conjugate_tree.c on word files.
#
#   make bench      build build/bench
#   make clean
#
# Then: build/bench <word file> [rounds] [--with-inverse], with word files
# written by extra/word_files.py (into words/, say).

CC ?= cc

SRC = ../../combisurf/src/conjugate_tree.c
HDR = ../../combisurf/src/conjugate_tree.h
FLAGS = -std=c99 -Wall -Wextra -Werror -O2
BUILD = build

.PHONY: bench clean

bench: $(BUILD)/bench

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/bench: bench_conjugate_tree.c $(SRC) $(HDR) | $(BUILD)
	$(CC) $(FLAGS) -o $@ bench_conjugate_tree.c $(SRC)

clean:
	rm -rf $(BUILD)
