# Tests in C of combisurf/src/conjugate_tree.c (the conjugate trees) and of
# combisurf/src/isomorphism.c (canonical keys and automorphisms of maps).
#
#   make check      optimized build, run
#   make debug      -O0 -DCT_DEBUG -DISO_DEBUG (local checks), run
#   make asan       AddressSanitizer and UBSan with the debug checks, run
#   make valgrind   -O0 build, run under valgrind
#   make clean
#
# The seed of the random cases is SEED (default 1): make check SEED=7

CC ?= cc
SEED ?= 1
VALGRIND ?= valgrind

CT_SRC = ../../combisurf/src/conjugate_tree.c
CT_HDR = ../../combisurf/src/conjugate_tree.h
CT_TEST = test_conjugate_tree.c
ISO_SRC = ../../combisurf/src/isomorphism.c
ISO_HDR = ../../combisurf/src/isomorphism.h
ISO_TEST = test_isomorphism.c
FLAGS = -std=c99 -Wall -Wextra -Werror
DEBUG_FLAGS = -DCT_DEBUG -DISO_DEBUG
ASAN_FLAGS = -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined -fno-sanitize-recover=all
BUILD = build

.PHONY: all check debug asan valgrind clean

all: check

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/test-check: $(CT_TEST) $(CT_SRC) $(CT_HDR) | $(BUILD)
	$(CC) $(FLAGS) -O2 -o $@ $(CT_TEST) $(CT_SRC)

$(BUILD)/test-debug: $(CT_TEST) $(CT_SRC) $(CT_HDR) | $(BUILD)
	$(CC) $(FLAGS) -O0 -g $(DEBUG_FLAGS) -o $@ $(CT_TEST) $(CT_SRC)

$(BUILD)/test-asan: $(CT_TEST) $(CT_SRC) $(CT_HDR) | $(BUILD)
	$(CC) $(FLAGS) $(ASAN_FLAGS) $(DEBUG_FLAGS) -o $@ $(CT_TEST) $(CT_SRC)

$(BUILD)/test-valgrind: $(CT_TEST) $(CT_SRC) $(CT_HDR) | $(BUILD)
	$(CC) $(FLAGS) -O0 -g -o $@ $(CT_TEST) $(CT_SRC)

$(BUILD)/iso-check: $(ISO_TEST) $(ISO_SRC) $(ISO_HDR) | $(BUILD)
	$(CC) $(FLAGS) -O2 -o $@ $(ISO_TEST) $(ISO_SRC)

$(BUILD)/iso-debug: $(ISO_TEST) $(ISO_SRC) $(ISO_HDR) | $(BUILD)
	$(CC) $(FLAGS) -O0 -g $(DEBUG_FLAGS) -o $@ $(ISO_TEST) $(ISO_SRC)

$(BUILD)/iso-asan: $(ISO_TEST) $(ISO_SRC) $(ISO_HDR) | $(BUILD)
	$(CC) $(FLAGS) $(ASAN_FLAGS) $(DEBUG_FLAGS) -o $@ $(ISO_TEST) $(ISO_SRC)

$(BUILD)/iso-valgrind: $(ISO_TEST) $(ISO_SRC) $(ISO_HDR) | $(BUILD)
	$(CC) $(FLAGS) -O0 -g -o $@ $(ISO_TEST) $(ISO_SRC)

check: $(BUILD)/test-check $(BUILD)/iso-check
	./$(BUILD)/test-check $(SEED)
	./$(BUILD)/iso-check $(SEED)

debug: $(BUILD)/test-debug $(BUILD)/iso-debug
	./$(BUILD)/test-debug $(SEED)
	./$(BUILD)/iso-debug $(SEED)

asan: $(BUILD)/test-asan $(BUILD)/iso-asan
	./$(BUILD)/test-asan $(SEED)
	./$(BUILD)/iso-asan $(SEED)

valgrind: $(BUILD)/test-valgrind $(BUILD)/iso-valgrind
	$(VALGRIND) --error-exitcode=1 --leak-check=full --errors-for-leak-kinds=all ./$(BUILD)/test-valgrind $(SEED)
	$(VALGRIND) --error-exitcode=1 --leak-check=full --errors-for-leak-kinds=all ./$(BUILD)/iso-valgrind $(SEED)

clean:
	rm -rf $(BUILD)
