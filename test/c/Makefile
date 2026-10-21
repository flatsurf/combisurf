# Tests in C of combisurf/src/conjugate_tree.c (the conjugate trees) and of
# combisurf/src/omap.c (canonical keys and automorphisms of maps, which
# uses combisurf/src/perm.c).
#
#   make check      optimized build, run
#   make debug      -O0 -DCT_DEBUG -DOMAP_DEBUG (local checks), run
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
OMAP_SRC = ../../combisurf/src/omap.c
OMAP_HDR = ../../combisurf/src/omap.h
OMAP_TEST = test_omap.c
PERM_SRC = ../../combisurf/src/perm.c
PERM_HDR = ../../combisurf/src/perm.h
FLAGS = -std=c99 -Wall -Wextra -Werror
DEBUG_FLAGS = -DCT_DEBUG -DOMAP_DEBUG
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

$(BUILD)/omap-check: $(OMAP_TEST) $(OMAP_SRC) $(OMAP_HDR) $(PERM_SRC) $(PERM_HDR) | $(BUILD)
	$(CC) $(FLAGS) -O2 -o $@ $(OMAP_TEST) $(OMAP_SRC) $(PERM_SRC)

$(BUILD)/omap-debug: $(OMAP_TEST) $(OMAP_SRC) $(OMAP_HDR) $(PERM_SRC) $(PERM_HDR) | $(BUILD)
	$(CC) $(FLAGS) -O0 -g $(DEBUG_FLAGS) -o $@ $(OMAP_TEST) $(OMAP_SRC) $(PERM_SRC)

$(BUILD)/omap-asan: $(OMAP_TEST) $(OMAP_SRC) $(OMAP_HDR) $(PERM_SRC) $(PERM_HDR) | $(BUILD)
	$(CC) $(FLAGS) $(ASAN_FLAGS) $(DEBUG_FLAGS) -o $@ $(OMAP_TEST) $(OMAP_SRC) $(PERM_SRC)

$(BUILD)/omap-valgrind: $(OMAP_TEST) $(OMAP_SRC) $(OMAP_HDR) $(PERM_SRC) $(PERM_HDR) | $(BUILD)
	$(CC) $(FLAGS) -O0 -g -o $@ $(OMAP_TEST) $(OMAP_SRC) $(PERM_SRC)

check: $(BUILD)/test-check $(BUILD)/omap-check
	./$(BUILD)/test-check $(SEED)
	./$(BUILD)/omap-check $(SEED)

debug: $(BUILD)/test-debug $(BUILD)/omap-debug
	./$(BUILD)/test-debug $(SEED)
	./$(BUILD)/omap-debug $(SEED)

asan: $(BUILD)/test-asan $(BUILD)/omap-asan
	./$(BUILD)/test-asan $(SEED)
	./$(BUILD)/omap-asan $(SEED)

valgrind: $(BUILD)/test-valgrind $(BUILD)/omap-valgrind
	$(VALGRIND) --error-exitcode=1 --leak-check=full --errors-for-leak-kinds=all ./$(BUILD)/test-valgrind $(SEED)
	$(VALGRIND) --error-exitcode=1 --leak-check=full --errors-for-leak-kinds=all ./$(BUILD)/omap-valgrind $(SEED)

clean:
	rm -rf $(BUILD)
