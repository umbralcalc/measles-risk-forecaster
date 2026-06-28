#!/usr/bin/env python3
"""Build the UTLA adjacency graph for the CAR spatial prior (gating check #4).

No geometry library is needed: ONS boundary products are topologically
consistent — where two areas share a border, both polygons trace it through the
*same* vertices. So two areas are rook-adjacent iff they share at least one
boundary segment (an ordered pair of consecutive vertices, normalised for
direction). This is exact for these boundaries, not an approximation.

Inputs : dat/utla_boundaries.geojson  (from fetch_geography.sh)
Outputs: dat/utla_nodes.csv           (code,name — one row per area)
         dat/utla_adjacency.csv       (code_a,code_b — one undirected edge per row)

Islands (no shared land border, e.g. Isle of Wight, Isles of Scilly) come out as
isolated nodes; that is geographically correct and is reported here so the CAR
prior can apply a nearest-neighbour fallback rather than silently leaving them
unpooled.
"""

import csv
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
GEOJSON = os.path.join(HERE, "utla_boundaries.geojson")
NODES_OUT = os.path.join(HERE, "utla_nodes.csv")
ADJ_OUT = os.path.join(HERE, "utla_adjacency.csv")

# Round vertices to ~0.1 m so identical shared vertices match despite float
# representation, without collapsing genuinely distinct points.
PRECISION = 6


def rings(geometry):
    """Yield every linear ring (exterior and holes) of a (Multi)Polygon."""
    gtype = geometry["type"]
    coords = geometry["coordinates"]
    if gtype == "Polygon":
        polys = [coords]
    elif gtype == "MultiPolygon":
        polys = coords
    else:
        return
    for poly in polys:
        for ring in poly:
            yield ring


def segments(ring):
    """Yield normalised (direction-independent) segments of a ring."""
    pts = [(round(x, PRECISION), round(y, PRECISION)) for x, y, *_ in ring]
    for a, b in zip(pts, pts[1:]):
        if a == b:
            continue
        yield (a, b) if a <= b else (b, a)


def main():
    with open(GEOJSON) as f:
        gj = json.load(f)

    nodes = {}  # code -> name
    seg_owners = defaultdict(set)  # segment -> set of codes touching it
    for feat in gj["features"]:
        props = feat["properties"]
        code, name = props["CTYUA24CD"], props["CTYUA24NM"]
        nodes[code] = name
        geom = feat.get("geometry")
        if geom is None:
            continue
        for ring in rings(geom):
            for seg in segments(ring):
                seg_owners[seg].add(code)

    # Two areas are adjacent iff they co-own at least one boundary segment.
    adj = defaultdict(set)
    for owners in seg_owners.values():
        if len(owners) < 2:
            continue
        for a in owners:
            for b in owners:
                if a != b:
                    adj[a].add(b)

    # --- write outputs ---
    with open(NODES_OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["code", "name"])
        for code in sorted(nodes):
            w.writerow([code, nodes[code]])

    edges = set()
    for a, nbrs in adj.items():
        for b in nbrs:
            edges.add((a, b) if a < b else (b, a))
    with open(ADJ_OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["code_a", "code_b"])
        for a, b in sorted(edges):
            w.writerow([a, b])

    # --- stats & sanity ---
    degrees = {c: len(adj[c]) for c in nodes}
    isolated = sorted(c for c in nodes if degrees[c] == 0)
    nz = [d for d in degrees.values() if d > 0]
    print(f"nodes: {len(nodes)}   edges: {len(edges)}")
    if nz:
        print(
            f"degree (non-isolated): mean={sum(nz)/len(nz):.2f} "
            f"min={min(nz)} max={max(nz)}"
        )
    print(f"isolated (islands): {len(isolated)} -> "
          f"{[nodes[c] for c in isolated]}")

    # Connected components (BFS).
    seen, comps = set(), []
    for start in nodes:
        if start in seen:
            continue
        stack, comp = [start], []
        seen.add(start)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        comps.append(comp)
    comps.sort(key=len, reverse=True)
    print(f"connected components: {len(comps)} "
          f"(largest={len(comps[0])}, singletons={sum(1 for c in comps if len(c)==1)})")

    # Known-pair sanity: print neighbours of a few well-understood UTLAs.
    name_to_code = {v: k for k, v in nodes.items()}
    for probe in ["Westminster", "Birmingham", "Cornwall"]:
        code = name_to_code.get(probe)
        if code:
            nbrs = sorted(nodes[c] for c in adj[code])
            print(f"  {probe} neighbours ({len(nbrs)}): {nbrs}")


if __name__ == "__main__":
    sys.exit(main())
