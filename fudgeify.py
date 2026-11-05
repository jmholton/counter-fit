#!/usr/bin/env python3
"""
fudgeify.py  - James Holton  Oct-26

Read the fullgeo output of molprobify and generate phenix geometry
correction restraints, accumulating sigma-tightening state in fudge_table.txt.

Each call is one ground_factory cycle: the worst maxadd outliers of each type
(or bonds=N, angles=N, torsions=N, planes=N, crowbars=N for one category)
are added or tightened (sigma/1.5, down to the floors reweight_opts uses), and
every other entry is relaxed (sigma*1.01); NONBOND and BOND entries are
dropped once relaxed to sigma 5.  A new BOND, ANGLE or TORSION entry changes
phenix's own restraint, so it starts at natural sigma/1.5.  Changes to
backbone bonds only take effect with phenix's CDL off (ground_factory
cdl=False): CDL rewrites them every macro-cycle.

Contacts get crowbars (NONBOND entries): a bond that only pushes, from dmin
out to dmin+2 (distance_ideal dmin+1, slack 1).  The worst probe clashes
(<prefix>_clashes.txt) and the worst scored contacts (fullgeo NONBOND
lines) are chosen separately, crowbars=N of each per cycle (default maxadd).  A hydrogen is replaced by its parent
heavy atom, since crowbarring a hydrogen only bends it, and dmin is the
parents' contact distance plus the overlap (plus 0.1 per hydrogen replaced,
for scored contacts).  Parents that phenix does not treat as a contact
(bonded or angle neighbours) get no crowbar.  A new crowbar starts at
sigma 0.4.  A crowbar chosen again once it is at the sigma floor is
lengthened instead: dmin grows by 0.1 per cycle, as in gf7.

phenix puts a crowbar in its bond graph, so the contacts between each
crowbar atom and the other atom's bonded neighbours (1-3) are switched off,
and the next shell (1-4) is scaled down.  Those contacts are still scored;
each cycle they are listed under their crowbar in <prefix>_crowbars.txt.

Planes out of plane get a second planarity restraint (PLANE entries), on the
heavy atoms of phenix's plane (default.geo), measured in the model
(<prefix>.pdb).  phenix refuses a plane it already has, so a heavy-atom copy
of one of its own planes drops its best-fitting atom; a 4-atom copy (a
peptide's CA, C, O, N) is restrained instead by an improper dihedral on the
same 4 atoms (a TORSION entry), which is the one way 4 atoms leave a plane.
Like crowbars, they add a restraint, so they start at phenix's sigma.

Table format (fudge_table.txt):
  NONBOND  C1 R1 A1 F1  C2 R2 A2 F2  dmin sigma
  BOND     C1 R1 A1 F1  C2 R2 A2 F2  value sigma
  ANGLE    C1 R1 A1 F1  C2 R2 A2 F2  C3 R3 A3 F3  value sigma
  TORSION  C1 R1 A1 F1  C2 R2 A2 F2  C3 R3 A3 F3  C4 R4 A4 F4  value sigma period
  PLANE    sigma  C1 R1 A1 F1  C2 R2 A2 F2  ...
A NONBOND with a symmetry mate ends in sym=OP, the operator applied to the
second atom.

Atom fields: C=chain, R=resi, A=atomname, F=altid  ("_" = blank)
"""

import sys
import os
import re
import math
import shutil
import argparse
from fractions import Fraction
from functools import lru_cache

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fullgeo", help="fullgeo file from molprobify")
    p.add_argument("--table",       default="fudge_table.txt")
    p.add_argument("--maxadd",      type=int,   default=1)
    p.add_argument("--crowbars",    type=int,   default=-1,
                   help="clashes and contacts crowbarred per cycle, each (default maxadd)")
    p.add_argument("--bonds",       type=int,   default=-1, help="bonds per cycle (default maxadd)")
    p.add_argument("--angles",      type=int,   default=-1, help="angles per cycle (default maxadd)")
    p.add_argument("--torsions",    type=int,   default=-1,
                   help="torsions, CB deviations included, per cycle (default maxadd)")
    p.add_argument("--planes",      type=int,   default=-1, help="planes per cycle (default maxadd)")
    p.add_argument("--sigma_fudge", type=float, default=1.0)
    p.add_argument("--outfile",     default="phenix_fudge.eff")
    p.add_argument("--geo",         default="default.geo")
    p.add_argument("--cbetadev",    default="")
    p.add_argument("--overridefile",default="")
    p.add_argument("--clashes",     default="", help="probe clashes (default <prefix>_clashes.txt)")
    p.add_argument("--watch",       default="", help="crowbar side effects (default <prefix>_crowbars.txt)")
    p.add_argument("--pdb",         default="", help="model, for plane deviations (default <prefix>.pdb)")

    # Accept "key=val" positional arguments in tcsh style
    argv_clean = []
    kv_args = {}
    for a in sys.argv[1:]:
        if "=" in a and not a.startswith("-"):
            k, v = a.split("=", 1)
            kv_args[k] = v
        else:
            argv_clean.append(a)

    args = p.parse_args(argv_clean)
    for k, v in kv_args.items():
        if hasattr(args, k):
            attr = getattr(args, k)
            if isinstance(attr, int):
                setattr(args, k, int(v))
            elif isinstance(attr, float):
                setattr(args, k, float(v))
            else:
                setattr(args, k, v)
        elif k == "fullgeo":
            args.fullgeo = v
    return args

# ---------------------------------------------------------------------------
# Atom field helpers
# ---------------------------------------------------------------------------

def blank(s):
    s = s.strip()
    return "_" if s == "" else s

def parse_atom(spec_word_list):
    """Given split words from one atom in the fullgeo pipe section:
       ATOM ALTID RESNAME CHAIN RESI
    Return (chain, resi, atomname, altid) normalized, resi as int string.
    """
    if len(spec_word_list) < 5:
        return None
    atom = spec_word_list[0].strip()
    altid = blank(spec_word_list[1])
    chain = blank(spec_word_list[3])
    resi  = spec_word_list[4].strip()
    return (chain, resi, atom, altid)

def parse_atom_specs(pipe_section):
    """Split the pipe-section into per-atom tuples.
    Atoms are separated by ' - '.
    Returns list of (C, R, A, F) tuples.
    """
    parts = re.split(r'\s+-\s+', pipe_section.strip())
    atoms = []
    for p in parts:
        words = p.split()
        a = parse_atom(words)
        if a:
            atoms.append(a)
    return atoms

def is_heavy(atomname):
    return not atomname.startswith("H")

def load_dihedrals(geo_path):
    """Atoms of every dihedral phenix already restrains, from a phenix.refine
    .geo file, as tuples of (C, R, A, F).  phenix refuses to *add these."""
    dihedrals = set()
    if not geo_path or not os.path.exists(geo_path):
        return dihedrals
    with open(geo_path) as f:
        lines = iter(f)
        for line in lines:
            if not line.startswith("dihedral pdb="):
                continue
            ids = [line] + [next(lines, "") for i in range(3)]
            atoms = []
            for s in ids:
                w = s.split('"')
                if len(w) < 2 or len(w[1]) < 15:
                    break
                pid = w[1]
                atoms.append((blank(pid[9]), pid[10:15].strip(), pid[0:4].strip(), blank(pid[4])))
            if len(atoms) == 4:
                dihedrals.add(tuple(atoms))
    return dihedrals

# ---------------------------------------------------------------------------
# Step 1: Extract proposals from fullgeo
# ---------------------------------------------------------------------------

def entry_key(fields, attrs):
    """Identity of an entry (atoms and symmetry, not sigma/value), as text."""
    t = fields[0]
    if t == "NONBOND" and not attrs.get("sym"):
        k = [t] + min(fields[1:5], fields[5:9]) + max(fields[1:5], fields[5:9])
    elif t in ("NONBOND", "BOND"):
        k = fields[:9]      # type + 2 atoms (4 fields each)
    elif t == "ANGLE":
        k = fields[:13]     # type + 3 atoms
    elif t == "TORSION":
        k = fields[:17]     # type + 4 atoms
    elif t == "PLANE":
        k = [t] + fields[2:]    # type sigma then atoms; key = atoms only
    else:
        k = fields
    if attrs.get("sym"):
        k = k + [attrs["sym"]]
    return " ".join(k)

def sigma_index(fields):
    """0-based index of sigma field in a fields list."""
    t = fields[0]
    if t == "TORSION":
        return len(fields) - 2
    elif t == "PLANE":
        return 1
    else:
        return len(fields) - 1

def extract_proposals(fullgeo_path, sigma_fudge, geo_path="", cbetadev_path="", override_path=""):
    """Return list of (fields, attrs) proposals from fullgeo, worst first."""
    proposals = []

    # Overridden torsions:  TORSION OVERRIDE target sigma | A_resi A_resi A_resi A_resi
    overridden = set()
    if override_path and os.path.exists(override_path):
        with open(override_path) as f:
            for line in f:
                if line.startswith("TORSION") and "OVERRIDE" in line and "|" in line:
                    overridden.add(tuple(line.split("|", 1)[1].split()))

    with open(fullgeo_path) as f:
        for line in f:
            line = line.rstrip()
            if not line or "|" not in line:
                continue
            head, pipe = line.split("|", 1)
            fields_head = head.split()
            if not fields_head:
                continue
            typ = fields_head[0]
            # contacts are handled by crowbar_proposals
            if typ not in ("BOND", "ANGLE", "TORSION", "PLANE"):
                continue

            try:
                energy = float(fields_head[1])
                delta  = float(fields_head[2])
                obs    = float(fields_head[3])
                ideal  = float(fields_head[4])
                sigma_nat = float(fields_head[5])
            except (IndexError, ValueError):
                continue

            if energy <= sigma_fudge:
                continue

            # contact with a symmetry mate: operator applies to the second atom
            sym = ""
            m = re.search(r"\s+symop=(\S+)", pipe)
            if m:
                sym = m.group(1)
                pipe = pipe[:m.start()]

            atoms = parse_atom_specs(pipe)

            if typ == "BOND":
                if len(atoms) < 2:
                    continue
                # phenix cannot *change a bond between symmetry mates
                if sym:
                    continue
                c1,r1,a1,f1 = atoms[0]
                c2,r2,a2,f2 = atoms[1]
                if not (is_heavy(a1) and is_heavy(a2)):
                    continue
                if a1 == "SG" and a2 == "SG":
                    continue
                sigma_use = max(sigma_nat, 0.001)
                proposals.append((["BOND", c1,r1,a1,f1, c2,r2,a2,f2, f"{ideal:.4f}", f"{sigma_use:.6f}"], {}))

            elif typ == "ANGLE":
                if len(atoms) < 3:
                    continue
                c1,r1,a1,f1 = atoms[0]
                c2,r2,a2,f2 = atoms[1]
                c3,r3,a3,f3 = atoms[2]
                if not all(is_heavy(a) for a in (a1,a2,a3)):
                    continue
                sigma_use = max(sigma_nat, 0.001)
                proposals.append((["ANGLE", c1,r1,a1,f1, c2,r2,a2,f2, c3,r3,a3,f3, f"{ideal:.4f}", f"{sigma_use:.6f}"], {}))

            elif typ == "TORSION":
                if len(atoms) < 4:
                    continue
                c1,r1,a1,f1 = atoms[0]
                c2,r2,a2,f2 = atoms[1]
                c3,r3,a3,f3 = atoms[2]
                c4,r4,a4,f4 = atoms[3]
                if not all(is_heavy(a) for a in (a1,a2,a3,a4)):
                    continue
                # Skip N/C-CA-CB (CBdev)
                if a3 == "CA" and a4 == "CB" and a1 in ("N","C") and a2 in ("N","C"):
                    continue
                # Skip torsions named in the override file (either direction)
                ov = tuple(f"{a}_{r}" for c, r, a, f in atoms[:4])
                if ov in overridden or ov[::-1] in overridden:
                    continue
                sigma_use = max(sigma_nat, 0.001)
                proposals.append((["TORSION", c1,r1,a1,f1, c2,r2,a2,f2, c3,r3,a3,f3, c4,r4,a4,f4,
                                   f"{ideal:.4f}", f"{sigma_use:.6f}", "0"], {}))

    # CBdev pseudo-dihedrals
    if cbetadev_path and os.path.exists(cbetadev_path):
        ideal_ncab = {
            "ala":122.9, "pro":115.1, "val":123.4, "thr":123.4,
            "ile":123.4, "gly":121.6,
        }
        ideal_cnab = {
            "ala":-122.6, "pro":-120.7, "val":-122.0, "thr":-122.0,
            "ile":-122.0, "gly":-121.6,
        }
        with open(cbetadev_path) as f:
            for line in f:
                cols = line.split(":")
                if len(cols) < 5:
                    continue
                try:
                    en = float(cols[1])
                except ValueError:
                    continue
                if en <= sigma_fudge:
                    continue
                alt = blank(cols[2]); res = cols[3].strip().lower()
                c   = blank(cols[4]); r   = cols[5].strip() if len(cols)>5 else "_"
                ncab = ideal_ncab.get(res, 122.8)
                cnab = ideal_cnab.get(res, -122.6)
                for v, desc in [(ncab,"N"), (cnab,"C")]:
                    proposals.append((["TORSION", c,r,desc,alt, c,r,("C" if desc=="N" else "N"),alt,
                                       c,r,"CA",alt, c,r,"CB",alt, f"{v:.4f}", "2.5", "0"], {}))

    return proposals

# ---------------------------------------------------------------------------
# Planes
# ---------------------------------------------------------------------------

def load_planes(geo_path):
    """Every planarity phenix restrains, from a phenix.refine .geo file, as
    lists of ((C, R, A, F), sigma)."""
    planes, cur = [], None
    if not geo_path or not os.path.exists(geo_path):
        return planes
    with open(geo_path) as f:
        for line in f:
            if line.startswith("plane pdb="):
                cur = []
                planes.append(cur)
            elif cur is None or not line.startswith("      pdb="):
                cur = None
                continue
            w = line.split('"')
            pid = w[1]
            cur.append(((blank(pid[9]), pid[10:15].strip(), pid[0:4].strip(), blank(pid[4])),
                        float(w[2].split()[1])))
    return planes

def load_sites(pdb_path):
    """Coordinates of each atom (C, R, A, F) in a PDB file."""
    sites = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom = (blank(line[21]), line[22:27].strip(), line[12:16].strip(), blank(line[16]))
                sites[atom] = tuple(float(line[i:i+8]) for i in (30, 38, 46))
    return sites

def sub(a, b):
    return [x - y for x, y in zip(a, b)]

def dot(a, b):
    return sum(x * y for x, y in zip(a, b))

def cross(a, b):
    return [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]

def smallest_eigenvector(m):
    """Unit eigenvector of the smallest eigenvalue of a symmetric 3x3 matrix,
    by Jacobi rotations."""
    a = [list(row) for row in m]
    v = [[float(i == j) for j in range(3)] for i in range(3)]
    for sweep in range(50):
        if sum(a[i][j]**2 for i in range(3) for j in range(3) if i != j) < 1e-30:
            break
        for p in range(2):
            for q in range(p + 1, 3):
                if a[p][q] == 0:
                    continue
                theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
                t = (1 if theta >= 0 else -1) / (abs(theta) + math.sqrt(theta*theta + 1))
                c = 1 / math.sqrt(t*t + 1)
                s = t * c
                for k in range(3):
                    a[k][p], a[k][q] = c*a[k][p] - s*a[k][q], s*a[k][p] + c*a[k][q]
                for k in range(3):
                    a[p][k], a[q][k] = c*a[p][k] - s*a[q][k], s*a[p][k] + c*a[q][k]
                for k in range(3):
                    v[k][p], v[k][q] = c*v[k][p] - s*v[k][q], s*v[k][p] + c*v[k][q]
    i = min(range(3), key=lambda i: a[i][i])
    return [v[k][i] for k in range(3)]

def plane_deviations(xyz, weights):
    """Distance of each site from the weighted least-squares plane, as phenix
    computes it."""
    w = sum(weights)
    centre = [sum(wi * p[k] for wi, p in zip(weights, xyz)) / w for k in range(3)]
    d = [sub(p, centre) for p in xyz]
    m = [[sum(wi * di[r] * di[c] for wi, di in zip(weights, d)) for c in range(3)] for r in range(3)]
    n = smallest_eigenvector(m)
    return [dot(di, n) for di in d]

def dihedral(p1, p2, p3, p4):
    b1, b2, b3 = sub(p2, p1), sub(p3, p2), sub(p4, p3)
    n1, n2 = cross(b1, b2), cross(b2, b3)
    m = cross(n1, [x / math.sqrt(dot(b2, b2)) for x in b2])
    return math.degrees(math.atan2(dot(m, n2), dot(n1, n2)))

def plane_proposals(fullgeo_path, planes, sites, neighbours, dihedrals, threshold):
    """Restraints for the planes whose atoms are out of plane, worst first.

    An atom is out of plane when both its score in fullgeo (after overrides)
    and its distance from the plane, measured here, exceed threshold: fullgeo
    does not say which of an atom's planes it is out of.  Planes are ranked
    by the fullgeo score, since molprobify scores its own hydrogens, which
    can move a plane's fit away from the one measured here.  The restraint is
    the plane's heavy atoms.  phenix refuses a plane it already has, so a
    copy of one of its own planes loses its best-fitting atom; a 4-atom copy
    (a peptide's CA, C, O, N) has one way out of plane, and the same 4 atoms
    are restrained by an improper dihedral instead: 0 or 180 degrees, sigma
    the plane's sigma over the 4th atom's distance from the dihedral axis.
    """
    scored = {}
    with open(fullgeo_path) as f:
        for line in f:
            if not line.startswith("PLANE") or "|" not in line:
                continue
            head, pipe = line.split("|", 1)
            try:
                energy = float(head.split()[1])
            except (IndexError, ValueError):
                continue
            atoms = parse_atom_specs(pipe)
            if atoms:
                scored[atoms[0]] = max(energy, scored.get(atoms[0], 0.0))

    own = {frozenset(a for a, s in plane) for plane in planes}
    found = []
    for plane in planes:
        atoms = [a for a, s in plane]
        if not all(a in sites for a in atoms):
            continue
        sigma = min(s for a, s in plane)
        dev = plane_deviations([sites[a] for a in atoms], [1 / s**2 for a, s in plane])
        energy = {a: scored.get(a, 0.0) if (d / s)**2 > threshold else 0.0
                  for (a, s), d in zip(plane, dev)}
        heavy = [a for a in atoms if is_heavy(a[2])]
        worst = max((energy[a] for a in heavy), default=0.0)
        if worst <= threshold or len(heavy) < 4:
            continue
        if frozenset(heavy) in own and len(heavy) > 4:
            heavy.remove(min(heavy, key=lambda a: abs(dev[atoms.index(a)])))
        if frozenset(heavy) not in own:
            fields = ["PLANE", f"{sigma:.6f}"] + [x for a in sorted(heavy) for x in a]
        elif len(heavy) == 4:
            # improper: the atom bonded to the other three goes last
            quad = tuple(sorted(heavy, key=lambda a: len(neighbours.get(a, set()) & set(heavy))))
            if quad in dihedrals or quad[::-1] in dihedrals:
                continue
            xyz = [sites[a] for a in quad]
            axis = sub(xyz[2], xyz[1])
            arm = sub(xyz[3], xyz[1])
            h = math.sqrt(dot(cross(axis, arm), cross(axis, arm)) / dot(axis, axis))
            ideal = 0.0 if abs(dihedral(*xyz)) < 90 else 180.0
            fields = (["TORSION"] + [x for a in quad for x in a] +
                      [f"{ideal:.4f}", f"{math.degrees(sigma / h):.6f}", "0"])
        else:
            continue
        found.append((worst, (fields, {"_via": "PLANE"})))
    found.sort(key=lambda e: -e[0])
    return [p for worst, p in found]

# ---------------------------------------------------------------------------
# Crowbars: restraints that push contacts apart
# ---------------------------------------------------------------------------

CROWBAR_SIGMA = 0.4   # sigma of a new crowbar
CROWBAR_STEP  = 0.1   # dmin added each cycle a crowbar at the floor is chosen

@lru_cache(maxsize=None)
def parse_symop(op):
    """A symmetry operator such as "-x+1/2,y,z+1" as a comparable
    (rotation, translation) pair; "" is the identity."""
    rot, trans = [], []
    for comp in (op or "x,y,z").lower().replace(" ", "").split(","):
        row, t = [0, 0, 0], Fraction(0)
        for sign, term in re.findall(r"([+-]?)([^+-]+)", comp):
            s = -1 if sign == "-" else 1
            if term in ("x", "y", "z"):
                row["xyz".index(term)] += s
            else:
                t += s * Fraction(term)
        rot.append(tuple(row))
        trans.append(t)
    return tuple(rot), tuple(trans)

IDENTITY = parse_symop("")

@lru_cache(maxsize=None)
def inverse_symop(op):
    """Inverse of a parse_symop operator (integer rotation, determinant +-1)."""
    (a, b, c), (d, e, f), (g, h, i) = op[0]
    det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)
    adj = ((e*i - f*h, c*h - b*i, b*f - c*e),
           (f*g - d*i, a*i - c*g, c*d - a*f),
           (d*h - e*g, b*g - a*h, a*e - b*d))
    rot = tuple(tuple(x // det for x in row) for row in adj)
    trans = tuple(-sum(rot[r][k] * op[1][k] for k in range(3)) for r in range(3))
    return rot, trans

def load_geometry(fullgeo_path):
    """Bonded neighbours of each atom, and every scored contact, from fullgeo.

    Contacts are keyed (atom1, atom2, operator applied to atom2) and hold
    (energy, delta, obs, ideal, sym text, fullgeo line).
    """
    neighbours, contacts = {}, {}
    with open(fullgeo_path) as f:
        for line in f:
            typ = line.split(None, 1)[0] if line.strip() else ""
            if typ not in ("BOND", "NONBOND") or "|" not in line:
                continue
            head, pipe = line.rstrip().split("|", 1)
            sym = ""
            m = re.search(r"\s+symop=(\S+)", pipe)
            if m:
                sym = m.group(1)
                pipe = pipe[:m.start()]
            atoms = parse_atom_specs(pipe)
            if len(atoms) < 2:
                continue
            x, y = atoms[0], atoms[1]
            if typ == "BOND":
                if not sym:
                    neighbours.setdefault(x, set()).add(y)
                    neighbours.setdefault(y, set()).add(x)
                continue
            try:
                energy, delta, obs, ideal = (float(w) for w in head.split()[1:5])
            except ValueError:
                continue
            contacts[(x, y, parse_symop(sym))] = (energy, delta, obs, ideal, sym, line.rstrip())
    return neighbours, contacts

def find_contact(contacts, x, y, op):
    """The scored contact between x and y (op applied to y) as
    (atom1, atom2, entry) in fullgeo's order, or None."""
    hit = contacts.get((x, y, op))
    if hit:
        return x, y, hit
    hit = contacts.get((y, x, inverse_symop(op)))
    if hit:
        return y, x, hit
    return None

def parent_atoms(neighbours):
    """Parent heavy atom of each hydrogen (the one in the hydrogen's own
    altloc when there is a choice)."""
    parent = {}
    for h, bonded in neighbours.items():
        heavy = sorted(a for a in bonded if is_heavy(a[2]))
        if not is_heavy(h[2]) and heavy:
            same = [a for a in heavy if a[3] == h[3]]
            parent[h] = (same or heavy)[0]
    return parent

def parse_probe_atom(s):
    """(C, R, A, F) from a probe atom id such as " A   3 BVAL HG21"."""
    s = s.ljust(16)
    return (blank(s[1]), s[2:6].strip(), s[12:16].strip(), blank(s[7]))

def crowbar_proposals(neighbours, contacts, clashes_path, threshold):
    """Crowbars for the probe clashes and for the scored contacts, each list
    worst first and tagged, so accumulate() picks from each separately."""
    parent = parent_atoms(neighbours)

    def crowbar(x, y, op, overlap, per_h):
        found = find_contact(contacts, parent.get(x, x), parent.get(y, y), op)
        if not found:
            return None
        (c1,r1,a1,f1), (c2,r2,a2,f2), hit = found
        dmin = hit[3] + overlap + per_h * ((x in parent) + (y in parent))
        attrs = {"sym": hit[4]} if hit[4] else {}
        return ["NONBOND", c1,r1,a1,f1, c2,r2,a2,f2, f"{dmin:.4f}", f"{CROWBAR_SIGMA:.4f}"], attrs

    lists = {"CLASH": [], "NONBOND": []}
    if clashes_path and os.path.exists(clashes_path):
        with open(clashes_path) as f:
            for line in f:
                w = line.split("|")
                if not line.startswith("CLASH") or len(w) < 3:
                    continue
                try:
                    energy, overlap = (float(v) for v in w[0].split()[1:3])
                except ValueError:
                    continue
                p = crowbar(parse_probe_atom(w[1]), parse_probe_atom(w[2]), IDENTITY, overlap, 0.0)
                if p:
                    lists["CLASH"].append((energy, p))
    for (x, y, op), (energy, delta, obs, ideal, sym, line) in contacts.items():
        if energy > threshold:
            p = crowbar(x, y, op, delta, 0.1)
            if p:
                lists["NONBOND"].append((energy, p))

    proposals = []
    for via, found in lists.items():
        found.sort(key=lambda e: -e[0])
        proposals += [(fields, dict(attrs, _via=via)) for energy, (fields, attrs) in found]
    return proposals

def watch_crowbars(entries, neighbours, contacts, watch_path, threshold):
    """Write, under each crowbar, the scored contacts phenix no longer
    restrains because of it: each crowbar atom against the other's bonded
    neighbours (1-3, switched off) and the next shell (1-4, scaled down).
    Contacts a crowbar covers (between crowbar atoms or their hydrogens) and
    contacts already 1-4 through covalent bonds are left out.  Returns how
    many were listed, and those scored above threshold, worst first."""
    def shells(a):
        n1 = neighbours.get(a, set())
        n2 = set().union(*(neighbours.get(b, set()) for b in n1)) - n1 - {a}
        return n1, n2

    def covalent_1_4(x, y):
        n1, n2 = shells(x)
        return y in n1 | n2 | set().union(*(neighbours.get(b, set()) for b in n2))

    parent = parent_atoms(neighbours)
    covered = set()
    for fields, attrs in entries:
        if fields[0] == "NONBOND":
            a, b = tuple(fields[1:5]), tuple(fields[5:9])
            op = parse_symop(attrs.get("sym", ""))
            covered |= {(a, b, op), (b, a, inverse_symop(op))}

    out, listed, close = [], 0, []
    for fields, attrs in entries:
        if fields[0] != "NONBOND":
            continue
        a, b = tuple(fields[1:5]), tuple(fields[5:9])
        op = parse_symop(attrs.get("sym", ""))
        a1, a2 = shells(a)
        b1, b2 = shells(b)
        pairs = ([("1-3", a, y) for y in b1] + [("1-3", x, b) for x in a1] +
                 [("1-4", a, y) for y in b2] + [("1-4", x, b) for x in a2] +
                 [("1-4", x, y) for x in a1 for y in b1])
        hits = []
        for kind, x, y in pairs:
            if (parent.get(x, x), parent.get(y, y), op) in covered:
                continue
            if kind == "1-4" and op == IDENTITY and covalent_1_4(x, y):
                continue
            found = find_contact(contacts, x, y, op)
            if found:
                hits.append((found[2][0], kind, found[2][5]))
        now = find_contact(contacts, a, b, op)
        out.append("CROWBAR " + " ".join(fields[1:] + [f"{k}={v}" for k, v in attrs.items()])
                   + (f"  now {now[2][2]:.3f} (contact distance {now[2][3]:.3f})" if now else ""))
        for energy, kind, line in sorted(hits, reverse=True):
            out.append(f"  {kind}  {line}")
            if energy > threshold:
                close.append((energy, kind, line))
        listed += len(hits)
    with open(watch_path, "w") as f:
        f.write("\n".join(out) + "\n")
    return listed, sorted(close, reverse=True)

# ---------------------------------------------------------------------------
# Step 2: Accumulate into fudge_table.txt
# ---------------------------------------------------------------------------

def load_table(table_path):
    """Return list of (fields, attrs) from existing table."""
    entries = []
    if not os.path.exists(table_path):
        return entries
    with open(table_path) as f:
        for line in f:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            words  = line.split()
            fields = [w for w in words if "=" not in w]
            attrs  = dict(w.split("=", 1) for w in words if "=" in w)
            if len(fields) >= 4:
                entries.append((fields, attrs))
    return entries

# Tightening floors, the same as reweight_opts_runme.com except crowbars (NONBOND),
# which can get as tight as a torsion before they are lengthened
FLOORS = {"NONBOND": 0.01, "BOND": 0.001, "ANGLE": 0.01, "TORSION": 0.01, "PLANE": 0.001}

def tighten(fields, sig):
    """sigma/1.5, not below the floor for its type, and never looser."""
    return max(sig / 1.5, min(sig, FLOORS.get(fields[0], 0.001)))

def accumulate(proposals, table_entries, maxadd, limits={}):
    """Tighten or add the worst proposals of each category: limits[category]
    of them, default maxadd.  The categories are the types, except that the
    crowbar lists are CLASH and NONBOND and planes are PLANE.  Relax the rest
    1%.  A crowbar tightened past its floor is lengthened instead.

    Returns the new entries and how many were tightened, lengthened and added.
    """
    chosen = {}
    count  = {}
    for p in proposals:
        t = p[1].get("_via", p[0][0])     # crowbars: clash or contact list; planes
        limit = limits.get(t, maxadd)
        k = entry_key(*p)
        if count.get(t, 0) < limit and k not in chosen:
            chosen[k] = p
            count[t] = count.get(t, 0) + 1

    new_entries = []
    in_table = set()
    nlengthened = 0
    for fields, attrs in table_entries:
        fields = list(fields)
        k   = entry_key(fields, attrs)
        t   = fields[0]
        si  = sigma_index(fields)
        sig = float(fields[si])
        in_table.add(k)
        if k in chosen:
            if t == "NONBOND" and sig / 1.5 < FLOORS[t]:
                fields[9] = f"{float(fields[9]) + CROWBAR_STEP:.4f}"
                nlengthened += 1
            sig = tighten(fields, sig)
        else:
            # Relax; drop loose nonbond/bond/plane
            sig *= 1.01
            if t in ("NONBOND","BOND","PLANE") and sig >= 5.0:
                continue
        fields[si] = f"{sig:.6f}"
        new_entries.append((fields, attrs))
    ntightened = len(in_table & chosen.keys())

    # Add new proposals not yet in table.  Crowbars and planes add a restraint;
    # the others change an existing one, so start them one step tighter than
    # natural or the change does nothing.
    nadded = 0
    for k, (fields, attrs) in chosen.items():
        if k not in in_table:
            fields = list(fields)
            adds = attrs.get("_via") in ("CLASH", "NONBOND", "PLANE")
            attrs = {a: v for a, v in attrs.items() if not a.startswith("_")}
            if not adds:
                si = sigma_index(fields)
                fields[si] = f"{tighten(fields, float(fields[si])):.6f}"
            new_entries.append((fields, attrs))
            nadded += 1

    return new_entries, ntightened, nlengthened, nadded

def write_table(entries, table_path):
    if os.path.exists(table_path):
        shutil.copy(table_path, table_path + ".bak")
    with open(table_path, "w") as f:
        for fields, attrs in entries:
            words = list(fields) + [f"{k}={v}" for k, v in attrs.items()]
            f.write(" ".join(words) + "\n")

# ---------------------------------------------------------------------------
# Step 3: Convert table to phenix_fudge.eff
# ---------------------------------------------------------------------------

def asel(c, r, a, f):
    """Build phenix atom_selection string, omitting chain/altid when blank."""
    s = f"name {a} and resseq {r}"
    if c != "_":
        s += f" and chain {c}"
    if f != "_":
        s += f" and altid {f}"
    return s

def write_phenix(entries, outfile, dihedrals=frozenset()):
    lines = ["refinement {", "  geometry_restraints.edits {"]

    for fields, attrs in entries:
        t = fields[0]

        if t == "NONBOND":
            c1,r1,a1,f1 = fields[1],fields[2],fields[3],fields[4]
            c2,r2,a2,f2 = fields[5],fields[6],fields[7],fields[8]
            dmin  = float(fields[9])
            sigma = max(float(fields[10]), 0.001)
            lines += [
                "    bond {",
                "      action = *add",
                f'      atom_selection_1 = "{asel(c1,r1,a1,f1)}"',
                f'      atom_selection_2 = "{asel(c2,r2,a2,f2)}"',
            ]
            if attrs.get("sym"):
                lines += [f"      symmetry_operation = {attrs['sym']}"]
            lines += [
                f"      distance_ideal = {dmin+1:.4f}",
                f"      sigma = {sigma:.6f}",
                "      slack = 1",
                "    }",
            ]

        elif t == "BOND":
            c1,r1,a1,f1 = fields[1],fields[2],fields[3],fields[4]
            c2,r2,a2,f2 = fields[5],fields[6],fields[7],fields[8]
            value = float(fields[9]); sigma = max(float(fields[10]), 0.001)
            lines += [
                "    bond {",
                "      action = *change",
                f'      atom_selection_1 = "{asel(c1,r1,a1,f1)}"',
                f'      atom_selection_2 = "{asel(c2,r2,a2,f2)}"',
                f"      distance_ideal = {value:.4f}",
                f"      sigma = {sigma:.6f}",
                "    }",
            ]

        elif t == "ANGLE":
            c1,r1,a1,f1 = fields[1],fields[2],fields[3],fields[4]
            c2,r2,a2,f2 = fields[5],fields[6],fields[7],fields[8]
            c3,r3,a3,f3 = fields[9],fields[10],fields[11],fields[12]
            value = float(fields[13]); sigma = max(float(fields[14]), 0.001)
            lines += [
                "    angle {",
                "      action = *change",
                f'      atom_selection_1 = "{asel(c1,r1,a1,f1)}"',
                f'      atom_selection_2 = "{asel(c2,r2,a2,f2)}"',
                f'      atom_selection_3 = "{asel(c3,r3,a3,f3)}"',
                f"      angle_ideal = {value:.4f}",
                f"      sigma = {sigma:.6f}",
                "    }",
            ]

        elif t == "TORSION":
            c1,r1,a1,f1 = fields[1],fields[2],fields[3],fields[4]
            c2,r2,a2,f2 = fields[5],fields[6],fields[7],fields[8]
            c3,r3,a3,f3 = fields[9],fields[10],fields[11],fields[12]
            c4,r4,a4,f4 = fields[13],fields[14],fields[15],fields[16]
            value  = float(fields[17])
            sigma  = max(float(fields[18]), 0.001)
            # change a dihedral phenix already has (it refuses a duplicate add);
            # periodicity 0 = harmonic (phenix's default of 1 is sinusoidal)
            quad = ((c1,r1,a1,f1), (c2,r2,a2,f2), (c3,r3,a3,f3), (c4,r4,a4,f4))
            action = "change" if (quad in dihedrals or quad[::-1] in dihedrals) else "add"
            lines += [
                "    dihedral {",
                f"      action = *{action}",
                f'      atom_selection_1 = "{asel(c1,r1,a1,f1)}"',
                f'      atom_selection_2 = "{asel(c2,r2,a2,f2)}"',
                f'      atom_selection_3 = "{asel(c3,r3,a3,f3)}"',
                f'      atom_selection_4 = "{asel(c4,r4,a4,f4)}"',
                f"      angle_ideal = {value:.4f}",
                f"      sigma = {sigma:.6f}",
                f"      periodicity = {fields[19]}",
                "    }",
            ]

        elif t == "PLANE":
            sigma = max(float(fields[1]), 0.001)
            atom_fields = fields[2:]
            atom_sels = []
            for i in range(0, len(atom_fields)-3, 4):
                c,r,a,f = atom_fields[i],atom_fields[i+1],atom_fields[i+2],atom_fields[i+3]
                atom_sels.append(asel(c,r,a,f))
            sel_str = "\n        or ".join(f"({s})" for s in atom_sels)
            lines += [
                "    planarity {",
                f'      atom_selection = "{sel_str}"',
                f"      sigma = {sigma:.6f}",
                "    }",
            ]

    lines += ["  }", "}"]
    with open(outfile, "w") as f:
        f.write("\n".join(lines) + "\n")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    if not os.path.exists(args.fullgeo):
        print(f"ERROR: fullgeo file not found: {args.fullgeo}", file=sys.stderr)
        sys.exit(1)
    prefix = re.sub(r"_fullgeo\.txt$", "", args.fullgeo)

    neighbours, contacts = load_geometry(args.fullgeo)
    dihedrals = load_dihedrals(args.geo)
    proposals = extract_proposals(args.fullgeo, args.sigma_fudge,
                                  args.geo, args.cbetadev, args.overridefile)
    proposals += crowbar_proposals(neighbours, contacts,
                                   args.clashes or prefix + "_clashes.txt", args.sigma_fudge)
    pdb = args.pdb or prefix + ".pdb"
    if os.path.exists(pdb):
        proposals += plane_proposals(args.fullgeo, load_planes(args.geo), load_sites(pdb),
                                     neighbours, dihedrals, args.sigma_fudge)
    else:
        print(f"fudgeify: no {pdb}, so no planes", file=sys.stderr)
    limits = {t: n for t, n in (("CLASH", args.crowbars), ("NONBOND", args.crowbars),
                                ("BOND", args.bonds), ("ANGLE", args.angles),
                                ("TORSION", args.torsions), ("PLANE", args.planes)) if n >= 0}
    new_entries, ntightened, nlengthened, nadded = accumulate(proposals, load_table(args.table),
                                                             args.maxadd, limits)

    write_table(new_entries, args.table)
    write_phenix(new_entries, args.outfile, dihedrals)
    listed, close = watch_crowbars(new_entries, neighbours, contacts,
                                   args.watch or prefix + "_crowbars.txt", args.sigma_fudge)

    ncrowbars = sum(1 for fields, attrs in new_entries if fields[0] == "NONBOND")
    print(f"fudgeify: wrote {args.outfile} from {args.fullgeo} ({len(new_entries)} table entries,"
          f" {nadded} added, {ntightened} tightened, {nlengthened} crowbars lengthened)")
    worst = ""
    if close:
        worst = f"; worst {' '.join(close[0][2].split('|', 1)[1].split())} (energy {close[0][0]:g})"
    print(f"fudgeify: {ncrowbars} crowbars switch off or weaken {listed} other scored contacts,"
          f" {len(close)} of them now too close{worst}")

if __name__ == "__main__":
    main()
