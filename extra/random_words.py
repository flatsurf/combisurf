r"""
Words whose conjugate trees have skewed shapes.

The generators produce abstract words over `\{0, \ldots, k - 1\}`, with no
free group meaning: they can be given to a conjugate tree as they are. To
turn one into a cyclically reduced word of the free group, whose letter
``h ^ 1`` is the inverse of the letter ``h``, use :func:`to_free_group`.

Every function returns an ``array('i')``, and every random one takes a
``random.Random`` as ``rng`` so that a seed reproduces its words.
"""
from array import array


def fibonacci(length, a=0, b=1):
    r"""
    Return the prefix of the given length of the Fibonacci word over ``a``
    and ``b``.

    The Fibonacci word is the fixed point of the substitution ``a -> ab``,
    ``b -> a``. It has exactly `n + 1` factors of length `n`, the fewest a
    non-periodic word can have.

    INPUT:

    - ``length`` -- length of the output

    - ``a``, ``b`` -- letters to use

    EXAMPLES::

        sage: fibonacci(13)
        array('i', [0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1])
        sage: fibonacci(13, 4, 6)
        array('i', [4, 6, 4, 4, 6, 4, 6, 4, 4, 6, 4, 4, 6])
    """
    w = array('i', [a, b])
    i = 1
    while len(w) < length:
        w.append(a)
        if w[i] == a:
            w.append(b)
        i += 1
    return w[:length]


def infinibonacci(length):
    r"""
    Return the prefix of the given length of the ruler sequence
    ``0 1 0 2 0 1 0 3 0 1 0 2 0 1 0 4 ...``.

    The letter at position `j`, counted from `1`, is the 2-adic valuation of
    `j`. It is the fixed point of the substitution ``k -> 0 (k + 1)``, over an
    infinite alphabet: a prefix of length `L` has `\lfloor \log_2 L \rfloor + 1`
    letters.

    INPUT:

    - ``length`` -- length of the output

    EXAMPLES::

        sage: infinibonacci(16)
        array('i', [0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4])
        sage: ruler(5)
        array('i', [0, 1, 0, 2, 0])
    """
    w = array('i', [0, 1])
    i = 1
    while len(w) < length:
        w.append(0)
        w.append(w[i] + 1)
        i += 1
    return w[:length]


ruler = infinibonacci


def word_random(length, alphabet, distribution, rng):
    r"""
    Return a word of i.i.d. letters drawn from ``alphabet`` with the
    probabilities ``distribution``.

    INPUT:

    - ``length`` -- length of the output

    - ``alphabet`` -- list of letters

    - ``distribution`` -- list of non-negative weights, one per letter (they
      need not sum to one)

    - ``rng`` -- a ``random.Random``

    EXAMPLES::

        sage: import random
        sage: word_random(10, [0, 3, 7], [0.1, 0.6, 0.3], random.Random(0))
        array('i', [7, 7, 3, 3, 3, 3, 7, 3, 3, 3])
    """
    if len(alphabet) != len(distribution):
        raise ValueError("alphabet and distribution must have the same length")
    return array('i', rng.choices(alphabet, weights=distribution, k=length))


def zipf(length, k, s, rng):
    r"""
    Return a word of i.i.d. letters of `\{0, \ldots, k - 1\}`, the letter
    `j` having probability proportional to `1 / (j + 1)^s`.

    With ``s = 0`` the letters are uniform, which gives the widest conjugate
    trees; as ``s`` grows, the letter ``0`` takes over and the tree becomes
    nearly unary.

    INPUT:

    - ``length`` -- length of the output

    - ``k`` -- size of the alphabet

    - ``s`` -- non-negative real, the exponent

    - ``rng`` -- a ``random.Random``

    EXAMPLES::

        sage: import random
        sage: zipf(12, 4, 1.0, random.Random(0))
        array('i', [2, 2, 0, 0, 1, 0, 2, 0, 0, 1, 3, 1])
        sage: zipf(12, 4, 8.0, random.Random(0))
        array('i', [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    """
    return word_random(length, list(range(k)), [1 / (j + 1) ** s for j in range(k)], rng)


def noisify(w_source, w_noise, k, rng):
    r"""
    Return a copy of ``w_source`` in which ``k`` distinct positions, drawn
    uniformly, get the first ``k`` letters of ``w_noise``, in order.

    INPUT:

    - ``w_source`` -- the word to modify

    - ``w_noise`` -- a word of length at least ``k``

    - ``k`` -- number of positions to change, at most ``len(w_source)``

    - ``rng`` -- a ``random.Random``

    EXAMPLES::

        sage: import random
        sage: rng = random.Random(0)
        sage: w_source = fibonacci(20)
        sage: w_noise = word_random(3, [3, 7, 9], [1, 1, 1], rng)
        sage: w_noise
        array('i', [9, 9, 7])
        sage: noisify(w_source, w_noise, 3, rng)
        array('i', [0, 1, 0, 0, 1, 0, 1, 0, 9, 1, 0, 0, 1, 0, 1, 9, 7, 1, 0, 1])
    """
    if k > len(w_noise):
        raise ValueError("w_noise must have at least k letters")
    w = array('i', w_source)
    for j, pos in enumerate(sorted(rng.sample(range(len(w)), k))):
        w[pos] = w_noise[j]
    return w


def to_free_group(w, n, rng):
    r"""
    Return the image of the abstract word ``w`` under a random injection of
    its letters into the letters `\{0, \ldots, n - 1\}` of the free group,
    that never uses both a letter ``h`` and its inverse ``h ^ 1``.

    The abstract letter ``j`` goes to ``2 * perm[j] + sign[j]``, where
    ``perm`` is a random injection into `\{0, \ldots, n/2 - 1\}` and
    ``sign[j]`` a random bit, both drawn once per call. The image contains no
    pair of inverse letters, so it is cyclically reduced.

    INPUT:

    - ``w`` -- a word over non-negative integers

    - ``n`` -- even size of the free group alphabet

    - ``rng`` -- a ``random.Random``

    EXAMPLES::

        sage: import random
        sage: to_free_group(fibonacci(8), 8, random.Random(0))
        array('i', [6, 3, 6, 6, 3, 6, 3, 6])
        sage: to_free_group(ruler(8), 4, random.Random(0))
        Traceback (most recent call last):
        ...
        ValueError: the word has 4 distinct letters, more than n/2 = 2
    """
    if n % 2:
        raise ValueError("n must be even")
    letters = sorted(set(w))
    if len(letters) > n // 2:
        raise ValueError(f"the word has {len(letters)} distinct letters, more than n/2 = {n // 2}")
    perm = rng.sample(range(n // 2), len(letters))
    image = {j: 2 * p + rng.randrange(2) for j, p in zip(letters, perm)}
    return array('i', [image[x] for x in w])
