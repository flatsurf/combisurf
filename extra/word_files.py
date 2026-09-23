r"""
Files of words, to test and benchmark conjugate trees from Python and from C.

A word file is a text file whose first line is ``# alphabet <n>``, followed
by free ``#`` comment lines (the generator and its parameters), then one word
per line, its letters written as decimal integers separated by single
spaces. The C benchmark ``extra/c/bench_conjugate_tree.c`` reads the same
files.

Usage::

    python extra/word_files.py <generator> <n> <L> <count> <seed> <path> [<param>] [--free-group]
    python extra/word_files.py bench <path> [<rounds>] [--with-inverse]

The first form writes ``count`` words of length ``L`` over the alphabet
`\{0, \ldots, n - 1\}` to ``path``, with the generator:

- ``random``: random cyclically reduced words of the free group, whose letter
  ``h ^ 1`` is the inverse of ``h`` (:func:`random_cyclically_reduced`);
- ``zipf <s>``: i.i.d. letters, the letter `j` with probability proportional
  to `1 / (j + 1)^s` (:func:`random_words.zipf`);
- ``fibonacci``, ``ruler``: the Fibonacci word over ``0, 1`` and the ruler
  sequence, the same for each of the ``count`` words;
- ``noisy-fibonacci <fraction>``: the Fibonacci word with ``round(fraction *
  L)`` positions replaced by letters drawn uniformly from
  `\{2, \ldots, \min(n, 16) - 1\}`.

With ``--free-group``, each word is generated over the abstract alphabet
`\{0, \ldots, n/2 - 1\}` and sent through :func:`random_words.to_free_group`,
which makes it a cyclically reduced word of the free group (and relabels the
deterministic words differently each time).

The second form builds a fresh ``ConjugateTree(n)`` and adds every word of
the file with ``process`` (or, with ``--with-inverse``, with
``crossing_arcs.tree_add_with_inverse``, which also adds its inverse), and
prints the best time over the rounds, per tree and per letter.
"""
import random
import sys
import time
from array import array

import random_words


def write_words(path, n, words, comment=None):
    r"""
    Write ``words`` over the alphabet of size ``n`` to the file ``path``.

    INPUT:

    - ``path`` -- file name

    - ``n`` -- size of the alphabet; every letter must be in `\{0, \ldots, n - 1\}`

    - ``words`` -- iterable of words (sequences of integers)

    - ``comment`` -- optional string, written as ``#`` lines after the header
    """
    with open(path, "w") as f:
        f.write(f"# alphabet {n}\n")
        if comment:
            for line in comment.splitlines():
                f.write(f"# {line}\n")
        for w in words:
            if any(x < 0 or x >= n for x in w):
                raise ValueError(f"letter outside the alphabet of size {n}")
            f.write(" ".join(map(str, w)))
            f.write("\n")


def read_words(path):
    r"""
    Return ``(n, words)`` read from the word file ``path``: the size of the
    alphabet and the list of words, each an ``array('i')``.

    INPUT:

    - ``path`` -- file name
    """
    with open(path) as f:
        header = f.readline().split()
        if header[:2] != ["#", "alphabet"] or len(header) != 3:
            raise ValueError(f"{path}: the first line must be '# alphabet <n>'")
        n = int(header[2])
        words = [array('i', map(int, line.split()))
                 for line in f if line.strip() and not line.startswith("#")]
    return n, words


def random_cyclically_reduced(n, L, rng):
    r"""
    Return a uniformly random cyclically reduced word of length ``L`` of the
    free group over the letters `\{0, \ldots, n - 1\}`, where ``h ^ 1`` is the
    inverse of ``h``.

    Each letter is drawn uniformly among those that do not cancel with the
    previous one, and the word is drawn again when its last letter cancels
    with its first: the draws, hence the words of a seed, do not depend on
    anything else.

    INPUT:

    - ``n`` -- even size of the alphabet

    - ``L`` -- length of the word

    - ``rng`` -- a ``random.Random``
    """
    while True:
        w = [rng.randrange(n)]
        for _ in range(L - 1):
            h = rng.randrange(n)
            while h == w[-1] ^ 1:
                h = rng.randrange(n)
            w.append(h)
        if w[0] != w[-1] ^ 1:
            return array('i', w)


def generate(generator, n, L, count, rng, param=None, free_group=False):
    r"""
    Return a list of ``count`` words of length ``L`` over `\{0, \ldots, n - 1\}`
    made by ``generator`` (see the module documentation).

    INPUT:

    - ``generator`` -- one of ``"random"``, ``"zipf"``, ``"fibonacci"``,
      ``"ruler"`` and ``"noisy-fibonacci"``

    - ``n``, ``L``, ``count`` -- alphabet size, length and number of words

    - ``rng`` -- a ``random.Random``

    - ``param`` -- the exponent of ``"zipf"``, the fraction of
      ``"noisy-fibonacci"``

    - ``free_group`` -- whether to send each word through
      :func:`random_words.to_free_group`
    """
    k = n // 2 if free_group else n
    if generator == "random":
        if free_group:
            raise ValueError("random words are free group words already")
        return [random_cyclically_reduced(n, L, rng) for _ in range(count)]
    if generator == "zipf":
        if param is None:
            raise ValueError("zipf needs its exponent s")
        words = [random_words.zipf(L, k, float(param), rng) for _ in range(count)]
    elif generator == "fibonacci":
        words = [random_words.fibonacci(L) for _ in range(count)]
    elif generator == "ruler":
        words = [random_words.ruler(L) for _ in range(count)]
    elif generator == "noisy-fibonacci":
        if param is None:
            raise ValueError("noisy-fibonacci needs the fraction of noisy positions")
        m = min(k, 16)
        if m <= 2:
            raise ValueError("noisy-fibonacci needs an alphabet of at least 3 letters")
        num = round(float(param) * L)
        words = []
        for _ in range(count):
            noise = random_words.word_random(num, list(range(2, m)), [1] * (m - 2), rng)
            words.append(random_words.noisify(random_words.fibonacci(L), noise, num, rng))
    else:
        raise ValueError(f"unknown generator {generator!r}")
    if free_group:
        words = [random_words.to_free_group(w, n, rng) for w in words]
    elif any(max(w) >= n for w in words):
        raise ValueError(f"{generator} words of length {L} need more than {n} letters")
    return words


def bench(path, rounds=5, with_inverse=False, min_time=0.01):
    r"""
    Return ``(n, letters, seconds)`` for the words of the file ``path``: the
    size of the alphabet, the number of letters of the file, and the best
    time over ``rounds`` rounds of the build of one fresh tree holding them.

    A round repeats the build until it has lasted ``min_time`` seconds and
    counts the mean.

    INPUT:

    - ``path`` -- a word file

    - ``rounds`` -- number of rounds

    - ``with_inverse`` -- whether to add each word with its inverse, with
      ``crossing_arcs.tree_add_with_inverse``; the alphabet must have even size

    - ``min_time`` -- the least duration of a round, in seconds
    """
    from combisurf.conjugate_tree import ConjugateTree
    from combisurf.crossing_arcs import tree_add_with_inverse

    n, words = read_words(path)
    if with_inverse and n % 2:
        raise ValueError(f"{path}: an alphabet of odd size {n} has no inverses")
    letters = sum(len(w) for w in words)

    def build():
        T = ConjugateTree(n)
        if with_inverse:
            for w in words:
                tree_add_with_inverse(T, w)
        else:
            process = T.process
            for w in words:
                process(w, False)

    best = float("inf")
    for _ in range(rounds):
        reps = 0
        t0 = time.perf_counter()
        while True:
            build()
            reps += 1
            t = time.perf_counter() - t0
            if t >= min_time:
                break
        best = min(best, t / reps)
    return n, letters, best


def _fmt(t):
    if t < 1e-3:
        return "%.3f us" % (t * 1e6)
    if t < 1:
        return "%.3f ms" % (t * 1e3)
    return "%.3f s" % t


def main(argv):
    flags = {a for a in argv if a.startswith("--")}
    args = [a for a in argv if not a.startswith("--")]
    if args and args[0] == "bench":
        if not 2 <= len(args) <= 3 or flags - {"--with-inverse"}:
            raise SystemExit("usage: word_files.py bench <path> [<rounds>] [--with-inverse]")
        rounds = int(args[2]) if len(args) == 3 else 5
        n, letters, t = bench(args[1], rounds, "--with-inverse" in flags)
        print(f"{args[1]}: n = {n}, {letters} letters, "
              f"{'with inverse' if '--with-inverse' in flags else 'plain'}: "
              f"{_fmt(t)} per tree, {t / letters * 1e9:.1f} ns per letter")
        return
    if not 6 <= len(args) <= 7 or flags - {"--free-group"}:
        raise SystemExit("usage: word_files.py <generator> <n> <L> <count> <seed> <path> "
                         "[<param>] [--free-group]")
    generator, n, L, count, seed, path = args[:6]
    n, L, count, seed = int(n), int(L), int(count), int(seed)
    param = args[6] if len(args) == 7 else None
    free_group = "--free-group" in flags
    rng = random.Random(seed)
    words = generate(generator, n, L, count, rng, param, free_group)
    comment = (f"generator {generator}" + (f" {param}" if param is not None else "")
               + f", n = {n}, L = {L}, count = {count}, seed = {seed}"
               + (", free group" if free_group else ""))
    write_words(path, n, words, comment)


if __name__ == "__main__":
    main(sys.argv[1:])
