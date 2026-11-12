# Tests of combisurf/src/conjugate_tree.c in C.
#
#   make check      optimized build, run
#   make debug      -O0 -DCT_DEBUG (local checks after each step), run
#   make asan       AddressSanitizer and UBSan with -DCT_DEBUG, run
#   make valgrind   -O0 build, run under valgrind
#   make clean
#
# The seed of the random cases is SEED (default 1): make check SEED=7

CC ?= cc
SEED ?= 1
VALGRIND ?= valgrind

SRC = ../../combisurf/src/conjugate_tree.c
HDR = ../../combisurf/src/conjugate_tree.h
TEST = test_conjugate_tree.c
FLAGS = -std=c99 -Wall -Wextra -Werror
BUILD = build

.PHONY: all check debug asan valgrind clean

all: check

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/test-check: $(TEST) $(SRC) $(HDR) | $(BUILD)
	$(CC) $(FLAGS) -O2 -o $@ $(TEST) $(SRC)

$(BUILD)/test-debug: $(TEST) $(SRC) $(HDR) | $(BUILD)
	$(CC) $(FLAGS) -O0 -g -DCT_DEBUG -o $@ $(TEST) $(SRC)

$(BUILD)/test-asan: $(TEST) $(SRC) $(HDR) | $(BUILD)
	$(CC) $(FLAGS) -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined \
		-fno-sanitize-recover=all -DCT_DEBUG -o $@ $(TEST) $(SRC)

$(BUILD)/test-valgrind: $(TEST) $(SRC) $(HDR) | $(BUILD)
	$(CC) $(FLAGS) -O0 -g -o $@ $(TEST) $(SRC)

check: $(BUILD)/test-check
	./$< $(SEED)

debug: $(BUILD)/test-debug
	./$< $(SEED)

asan: $(BUILD)/test-asan
	./$< $(SEED)

valgrind: $(BUILD)/test-valgrind
	$(VALGRIND) --error-exitcode=1 --leak-check=full --errors-for-leak-kinds=all ./$< $(SEED)

clean:
	rm -rf $(BUILD)
