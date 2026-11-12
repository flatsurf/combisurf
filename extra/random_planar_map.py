import subprocess

from combisurf.oriented_map import OrientedMap

PLANAR_MAP_BIN = "/home/vincent/programming/flatsurf/PlanarMap/planarmap"

def random_planar_map(*args, as_list=False):
    r"""
    Calls Schaeffer planarmap program with the given arguments.

    map selection:

    -C<2,3,4>: 2,3,4-edge-connected cubic maps
    -Q<1,2,3>: 2,4,6-edge-connected quartic maps
    -Q4      : bipartite quartic maps
    -M<1,2,3>: 1,2,3-connected maps
    -M4      : bipartite bicolor maps
    -B<1,2>  : 2,3-edge-connected bipartite cubic

    Parameters:

    -N<nb>: number of maps generated (default is one)
    -E<nb>: number of edges
    -V<nb>: number of vertices (=edges for -Mk)
    -F<nb>: number of faces
    -R<nb>: number of red faces in Q2 (=vertices for -M2)
    -G<nb>: number of green faces in Q2 (=faces for -M2)
    -I<nb>: approximate size [N+/-I]

    Methods:

    -l extraction by largest component (default)
    -c extraction by core
    -s suppress pic optimisation in extraction
    """
    params = {}
    map_selection = None

    for arg in args:
        if not arg:
            # ignore empty argument
            continue

        if not isinstance(arg, str) and arg.startswith("-") and len(arg) >= 2:
            raise TypeError("each argument must be a string")

        c = arg[1]
        if c in "CQMB":
            # - C<2,3,4>: 2,3,4-edge-connected cubic maps
            # - Q<1,2,3>: 2,4,6-edge-connected quartic maps
            # -Q4      : bipartite quartic maps
            # -M<1,2,3>: 1,2,3-connected maps
            # -M4      : bipartite bicolor maps
            # -B<1,2>  : 2,3-edge-connected bipartite cubic
            if map_selection is not None:
                raise ValueError("multiple map selectors provided")
            map_selection = arg[1:]
        elif c in "NEVFRGI":
            # -N<nb>: number of maps generated (default is one)
            if c in params:
                raise ValueError(f"-{c}<nb> specified more than once")
            value = arg[2:]
            if not value or not value.isdigit():
                raise ValueError("option -N must be followed by an integer")
            params[c] = int(arg[2:])

    # TODO: analyze arguments and check output
    s = subprocess.run((PLANAR_MAP_BIN,) + args + ("-O3", "-p"), capture_output=True)
    if s.returncode != 0:
        raise ValueError(s.stderr)
    text = s.stdout

    ans = []

    vp = None
    for line in text.splitlines():
        if line.startswith(b"Map"):
            if vp is not None:
                ans.append(OrientedMap(vp="".join(vp)))
            vp = []
            s = set()
        else:
            line = line[2:-2].strip()
            while line and line[-1] in b" ]":
                line = line[:-1]
            i = line.find(b",")
            vert = [int(x.strip()) for x in line[i+3:].split(b",")]
            s.update(vert)
            vert = [f"{x - 1}" if x > 0 else f"~{-x - 1}" for x in vert]
            vp.append("(" + ",".join(map(str,vert)) + ")")
    if vp is not None:
        ans.append(OrientedMap(vp="".join(vp)))

    if not as_list and len(ans) == 1:
        return ans[0]
    else:
        return ans
