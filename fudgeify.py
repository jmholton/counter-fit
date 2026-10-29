#!/usr/bin/env python3
"""
fudgeify.py  - James Holton  Oct-26

Read the fullgeo output of molprobify and generate phenix geometry
correction restraints, accumulating sigma-tightening state in fudge_table.txt.

Table format (fudge_table.txt):
  NONBOND  C1 R1 A1 F1  C2 R2 A2 F2  dmin sigma
  BOND     C1 R1 A1 F1  C2 R2 A2 F2  value sigma
  ANGLE    C1 R1 A1 F1  C2 R2 A2 F2  C3 R3 A3 F3  value sigma
  TORSION  C1 R1 A1 F1  C2 R2 A2 F2  C3 R3 A3 F3  C4 R4 A4 F4  value sigma period
  PLANE    sigma  C1 R1 A1 F1  C2 R2 A2 F2  ...

Atom fields: C=chain, R=resi, A=atomname, F=altid  ("_" = blank)
"""

import sys
import os
import re
import shutil
import argparse

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fullgeo", help="fullgeo file from molprobify")
    p.add_argument("--table",       default="fudge_table.txt")
    p.add_argument("--maxadd",      type=int,   default=1)
    p.add_argument("--sigma_fudge", type=float, default=1.0)
    p.add_argument("--outfile",     default="phenix_fudge.eff")
    p.add_argument("--geo",         default="")
    p.add_argument("--cbetadev",    default="")
    p.add_argument("--overridefile",default="")

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

# ---------------------------------------------------------------------------
# Step 1: Extract proposals from fullgeo
# ---------------------------------------------------------------------------

SIGMA_COLS = {
    "NONBOND": -1,   # last field
    "BOND":    -1,
    "ANGLE":   -1,
    "TORSION": -2,   # second-to-last (period is last)
    "PLANE":    1,   # index 1 (0-based) after split
}

def table_key(fields):
    """Return the identity key for an entry (atom coordinates, not sigma/value)."""
    t = fields[0]
    if t in ("NONBOND", "BOND"):
        return tuple(fields[:9])    # type + 2 atoms (4 fields each)
    elif t == "ANGLE":
        return tuple(fields[:13])   # type + 3 atoms
    elif t == "TORSION":
        return tuple(fields[:17])   # type + 4 atoms
    elif t == "PLANE":
        return tuple(fields[2:])    # type sigma then atoms; key = atoms only
    return tuple(fields)

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
    """Return list of field-lists (proposals) from fullgeo."""
    proposals = []

    # Load override residues to skip torsions
    overridden = set()
    if override_path and os.path.exists(override_path):
        with open(override_path) as f:
            for line in f:
                if "TORSION" in line and "OVERRIDE" in line:
                    overridden.add(line.strip())

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
            if typ not in ("NONBOND", "BOND", "ANGLE", "TORSION", "PLANE"):
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

            atoms = parse_atom_specs(pipe)

            if typ == "NONBOND":
                if len(atoms) < 2:
                    continue
                c1,r1,a1,f1 = atoms[0]
                c2,r2,a2,f2 = atoms[1]
                if not (is_heavy(a1) and is_heavy(a2)):
                    continue
                dmin  = ideal
                sigma = 0.4 / (abs(delta) + 1.1)
                proposals.append(["NONBOND", c1,r1,a1,f1, c2,r2,a2,f2, f"{dmin:.4f}", f"{sigma:.4f}"])

            elif typ == "BOND":
                if len(atoms) < 2:
                    continue
                c1,r1,a1,f1 = atoms[0]
                c2,r2,a2,f2 = atoms[1]
                if not (is_heavy(a1) and is_heavy(a2)):
                    continue
                if a1 == "SG" and a2 == "SG":
                    continue
                sigma_use = max(sigma_nat, 0.001)
                proposals.append(["BOND", c1,r1,a1,f1, c2,r2,a2,f2, f"{ideal:.4f}", f"{sigma_use:.6f}"])

            elif typ == "ANGLE":
                if len(atoms) < 3:
                    continue
                c1,r1,a1,f1 = atoms[0]
                c2,r2,a2,f2 = atoms[1]
                c3,r3,a3,f3 = atoms[2]
                if not all(is_heavy(a) for a in (a1,a2,a3)):
                    continue
                sigma_use = max(sigma_nat, 0.001)
                proposals.append(["ANGLE", c1,r1,a1,f1, c2,r2,a2,f2, c3,r3,a3,f3, f"{ideal:.4f}", f"{sigma_use:.6f}"])

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
                # Check overrides (simple: skip if any override line mentions this residue+chain)
                skip = False
                for ov in overridden:
                    if r1 in ov and c1 in ov:
                        skip = True
                        break
                if skip:
                    continue
                sigma_use = max(sigma_nat, 0.001)
                proposals.append(["TORSION", c1,r1,a1,f1, c2,r2,a2,f2, c3,r3,a3,f3, c4,r4,a4,f4,
                                   f"{ideal:.4f}", f"{sigma_use:.6f}", "0"])

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
                    proposals.append(["TORSION", c,r,desc,alt, c,r,("C" if desc=="N" else "N"),alt,
                                       c,r,"CA",alt, c,r,"CB",alt, f"{v:.4f}", "2.5", "0"])

    # PLANE from .geo file (only first atom per plane given in fullgeo, so skip for now unless geo provided)
    # geo parsing is complex; stub for future implementation
    if geo_path and os.path.exists(geo_path):
        pass  # TODO: parse phenix .geo plane records

    return proposals

# ---------------------------------------------------------------------------
# Step 2: Accumulate into fudge_table.txt
# ---------------------------------------------------------------------------

def load_table(table_path):
    """Return list of field-lists from existing table."""
    entries = []
    if not os.path.exists(table_path):
        return entries
    with open(table_path) as f:
        for line in f:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) >= 4:
                entries.append(fields)
    return entries

def accumulate(proposals, table_entries, maxadd):
    """Escalate existing sigmas, tighten matches, add new entries."""
    # Index proposals by key (up to maxadd per type)
    prop_by_key = {}
    prop_count  = {}
    for p in proposals:
        t = p[0]
        k = table_key(p)
        if prop_count.get(t, 0) < maxadd:
            prop_by_key[k] = p
            prop_count[t] = prop_count.get(t, 0) + 1

    matched = set()
    new_entries = []

    for fields in table_entries:
        if len(fields) < 4:
            continue
        t  = fields[0]
        si = sigma_index(fields)
        # Escalate sigma 1%
        sig = float(fields[si]) * 1.01
        # Drop loose nonbond/bond
        if t in ("NONBOND","BOND") and sig >= 5.0:
            continue
        fields[si] = f"{sig:.6f}"
        k = table_key(fields)
        if k in prop_by_key:
            floor = {"ANGLE": 0.5, "TORSION": 5.0}.get(t, 0.001)
            sig = max(sig / 1.5, floor)
            fields[si] = f"{sig:.6f}"
            matched.add(k)
        new_entries.append(fields)

    # Add new proposals not yet in table
    for k, p in prop_by_key.items():
        if k not in matched:
            new_entries.append(p)

    return new_entries

def write_table(entries, table_path):
    with open(table_path + ".bak", "w") as f:
        pass  # backup will be overwritten by shutil.copy below
    if os.path.exists(table_path):
        shutil.copy(table_path, table_path + ".bak")
    with open(table_path, "w") as f:
        for fields in entries:
            f.write(" ".join(str(x) for x in fields) + "\n")

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

def write_phenix(entries, outfile):
    lines = ["refinement {", "  geometry_restraints.edits {"]

    for fields in entries:
        t = fields[0]

        if t == "NONBOND":
            c1,r1,a1,f1 = fields[1],fields[2],fields[3],fields[4]
            c2,r2,a2,f2 = fields[5],fields[6],fields[7],fields[8]
            dmin  = float(fields[9])
            sigma = float(fields[10])
            lines += [
                "    bond {",
                "      action = *add",
                f'      atom_selection_1 = "{asel(c1,r1,a1,f1)}"',
                f'      atom_selection_2 = "{asel(c2,r2,a2,f2)}"',
                f"      distance_ideal = {dmin+1:.4f}",
                f"      sigma = {sigma:.6f}",
                "      slack = 1",
                "    }",
            ]

        elif t == "BOND":
            c1,r1,a1,f1 = fields[1],fields[2],fields[3],fields[4]
            c2,r2,a2,f2 = fields[5],fields[6],fields[7],fields[8]
            value = float(fields[9]); sigma = float(fields[10])
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
            value = float(fields[13]); sigma = float(fields[14])
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
            sigma  = float(fields[18])
            lines += [
                "    dihedral {",
                "      action = *add",
                f'      atom_selection_1 = "{asel(c1,r1,a1,f1)}"',
                f'      atom_selection_2 = "{asel(c2,r2,a2,f2)}"',
                f'      atom_selection_3 = "{asel(c3,r3,a3,f3)}"',
                f'      atom_selection_4 = "{asel(c4,r4,a4,f4)}"',
                f"      angle_ideal = {value:.4f}",
                f"      sigma = {sigma:.6f}",
                "    }",
            ]

        elif t == "PLANE":
            sigma = float(fields[1])
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

    proposals     = extract_proposals(args.fullgeo, args.sigma_fudge,
                                       args.geo, args.cbetadev, args.overridefile)
    table_entries = load_table(args.table)
    new_entries   = accumulate(proposals, table_entries, args.maxadd)
    write_table(new_entries, args.table)
    write_phenix(new_entries, args.outfile)

    print(f"fudgeify: wrote {args.outfile} from {args.fullgeo} ({len(new_entries)} table entries)")

if __name__ == "__main__":
    main()
