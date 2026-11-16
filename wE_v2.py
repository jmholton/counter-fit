#!/usr/bin/env python3
"""wE_v2.py: Python port of molprobify_v2.com, the version-2 weighted energy (wE) of the Untangle Challenge.

Reproduces claude_numbers/wE_v2/molprobify_v2.com (md5 c461987c, 23 Oct 2026) term by term, including the awk/gnuplot
rounding of intermediate numbers to 6 significant digits, so the printed wE matches the tcsh scorer.  Documentation:
WE_V2.md next to this file.

Backends (backend=auto|mmtbx|cli):
    - In-process mmtbx (default under phenix.python when keep/workdir are not set): builds the model and restraints once
      in memory and calls ramalyze, rotalyze, cbetadev and clashscore directly (~7-10x faster than spawning 5 subprocesses
      and writing/parsing a multi-MB .geo file, while matching molprobify_v2.com to the printed digit).
    - Subprocess CLI (backend=cli, or when keep=1 / workdir=DIR is requested): runs phenix.ramalyze, rotalyze, cbetadev,
      molprobity and geometry_minimization (optionally in parallel with jobs=N), then scores the output files.

Usage (same key=value style as the tcsh scorer):
    wE_v2.py model.pdb [model2.pdb ...] [ligand.cif ...] [nbmode=occ] [nbN=occ] [cbsigma=0.032] ...
        scores one or more PDB models (batch mode amortizes Python/cctbx/table startup across models).
    wE_v2.py model.pdb geo=x.geo rama=r.log rota=o.log cbdev=c.log [coot=molprobity_coot.py]
        scores outputs you already have (no phenix needed; the pdb is read for occupancies only).
Python API:
    import wE_v2
    r = wE_v2.score_model('model.pdb')                        # in-process mmtbx if available, else phenix CLI
    rs = wE_v2.score_models(['m1.pdb', 'm2.pdb'])             # batch scoring in one warm process
    scorer = wE_v2.ModelScorer('template.pdb')                # build topology/restraints once, re-score new coordinates
    r = scorer.score_sites(new_sites_cart)
    r = wE_v2.score_files(pdb, geo, rama, rota, cbdev, coot)  # existing outputs (~5x faster .geo parser)
    wE_v2.score_items(...)                                    # bring your own deviates (see WE_V2.md)
Needs only the Python 3.6+ standard library for score_files / score_items / backend=cli.
"""
import contextlib, math, os, re, sys, json, shutil, subprocess, tempfile

DEFAULTS = dict(
    nbmode='occ',      # occ (default) | bonds | stat | physical | both
    nbN='occ',         # occ: total atom occupancy | vdw120 | vdw100 (occupancy-weighted contacts within 1.2 / 1.0 x vdW)
    nboutw=0,          # 1: weight nonbond OUTLIER items by occupancy too
    nbP0=1.045e-3,     # nbmode=stat tail: P(LJ' >= LJ) = nbP0 exp(-nbk LJ)
    nbk=0.667,
    cbsigma=0.032,     # per-component sigma (A) of the C-beta deviation
    cbdof=3,           # 3: Maxwell (3D distance); 2: Rayleigh
    pnn='exact',       # exact: erf(z/sqrt2)**N; soft: v1's softPnna
    csda=0.6,          # const_shrink_donor_acceptor for phenix.geometry_minimization
    cdl=True,          # conformation-dependent library for phenix.geometry_minimization
    cdl_interp=0,      # 1: interpolate CDL ideals between phi/psi bins; 0: phenix default
    sigomega=4.01398721805631,
    ignore='ARGH',     # comma list of ARGH (H-H nonbonds between ARG), HOHH (anything with a water H), HBON (bonds to H)
    round=1,           # 1: emulate awk/gnuplot/.geo intermediate rounding; 0 (or round_numbers=False): full float64
    backend='auto',    # auto (mmtbx if available and not keep/workdir, else cli) | mmtbx | cli
    jobs=1,            # max concurrent phenix subprocesses when backend=cli (1..5)
    overridefile='',   # TORSION/BOND/... OVERRIDE lines (counter-fit): new ideal and sigma for named restraints
)

CATEGORIES = ['BOND', 'ANGLE', 'TORSION', 'CHIR', 'PLANE', 'OMEGA', 'RAMA', 'ROTA', 'CBDEV', 'NONBOND']
_AWKNUM_RE = re.compile(r'\s*[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?')

# ---------------------------------------------------------------------------------------------------------------- options
def _as_bool(v, name):
    if isinstance(v, bool):
        return v
    if isinstance(v, int) and v in (0, 1):
        return bool(v)
    s = str(v).strip().lower()
    if s in ('1', 'true', 'yes', 'on'):
        return True
    if s in ('0', 'false', 'no', 'off'):
        return False
    raise ValueError('invalid boolean %s=%r' % (name, v))

def _normalize_opts(opts):
    o = dict(DEFAULTS)
    for k, v in (opts or {}).items():
        if k == 'round_numbers':
            k = 'round'
        elif k == 'override' or k.startswith('override'):
            k = 'overridefile'
        if k not in DEFAULTS:
            raise ValueError('unknown option %r' % (k,))
        o[k] = v
    o['nboutw'] = int(_as_bool(o['nboutw'], 'nboutw'))
    o['round'] = int(_as_bool(o['round'], 'round'))
    o['cdl'] = _as_bool(o['cdl'], 'cdl')
    o['cdl_interp'] = int(_as_bool(o['cdl_interp'], 'cdl_interp'))
    o['cbdof'] = int(o['cbdof'])
    o['jobs'] = int(o['jobs'])
    if o['cbdof'] not in (2, 3):
        raise ValueError('invalid cbdof=%r (expected 2 or 3)' % (o['cbdof'],))
    if o['jobs'] < 1:
        raise ValueError('invalid jobs=%r (expected >= 1)' % (o['jobs'],))
    if o['nbmode'] not in ('occ', 'bonds', 'stat', 'physical', 'both'):
        raise ValueError('invalid nbmode=%r (expected occ|bonds|stat|physical|both)' % (o['nbmode'],))
    if o['nbN'] not in ('occ', 'vdw120', 'vdw100'):
        raise ValueError('invalid nbN=%r (expected occ|vdw120|vdw100)' % (o['nbN'],))
    if o['pnn'] not in ('exact', 'soft'):
        raise ValueError('invalid pnn=%r (expected exact|soft)' % (o['pnn'],))
    if o['backend'] not in ('auto', 'mmtbx', 'cli'):
        raise ValueError('invalid backend=%r (expected auto|mmtbx|cli)' % (o['backend'],))
    for k in ('nbP0', 'nbk', 'cbsigma', 'sigomega', 'csda'):
        o[k] = float(o[k])
    o['overridefile'] = str(o['overridefile'] or '')
    return o

# ---------------------------------------------------------------------------------------------------------------- numbers
def awkfmt(x, digits=6):
    """A number as awk prints it (integers as integers, otherwise %.6g), read back as a float."""
    if x != x or x in (float('inf'), float('-inf')):
        return x
    if x == int(x) and abs(x) < 1e16:
        return float(int(x))
    return float('%.*g' % (digits, x))

def g6(x):
    return awkfmt(x, 6)

def g15(x):        # gnuplot's print
    return awkfmt(x, 15)

def _identity(x):
    return x

def awknum(s):
    """awk's string-to-number: direct float when clean, else the longest numeric prefix, else 0."""
    if '_' not in s:
        try:
            v = float(s)
            if v == v and v not in (float('inf'), float('-inf')):
                return v
        except ValueError:
            pass
    m = _AWKNUM_RE.match(s)
    return float(m.group(0)) if m else 0.0

def soft(e):
    """energies above 10 are softened to 10 + ln(e/10), as in v1"""
    return 10 + math.log(e) - math.log(10) if e > 10 else e

try:
    from statistics import NormalDist
    _ND = NormalDist()
except ImportError:          # Python < 3.8
    _ND = None

def ndtri(p):
    """inverse of the standard normal CDF (gnuplot invnorm)"""
    if p <= 0.0:
        return float('-inf')
    if p >= 1.0:
        return float('inf')
    if _ND is not None:
        return _ND.inv_cdf(p)
    lower = p < 0.5
    q = p if lower else 1.0 - p
    x = -math.sqrt(-2.0 * math.log(q))                      # Newton on log CDF, lower tail
    for _ in range(200):
        c = 0.5 * math.erfc(-x / math.sqrt(2.0))
        if c <= 0:
            x += 1.0
            continue
        dx = (math.log(c) - math.log(q)) / (math.exp(-0.5 * x * x) / math.sqrt(2 * math.pi) / c)
        x -= dx
        if abs(dx) < 1e-15 * (1 + abs(x)):
            break
    return x if lower else -x

def inverf(y):
    """inverse error function (gnuplot inverf)"""
    if y >= 1.0:
        return float('inf')
    if y <= -1.0:
        return float('-inf')
    if abs(y) < 1e-5:
        y2 = y * y
        return 0.5 * math.sqrt(math.pi) * y * (1.0 + (math.pi / 12.0) * y2 + (7.0 * math.pi * math.pi / 480.0) * y2 * y2)
    if y >= 0:
        return -ndtri((1.0 - y) / 2.0) / math.sqrt(2.0)
    return ndtri((1.0 + y) / 2.0) / math.sqrt(2.0)

def igamma(a, x):
    """regularized lower incomplete gamma P(a, x) (gnuplot igamma)"""
    if x <= 0:
        return 0.0
    lnpre = -x + a * math.log(x) - math.lgamma(a)
    if x < a + 1:
        ap, s, d = a, 1.0 / a, 1.0 / a
        for _ in range(1000000):
            ap += 1
            d *= x / ap
            s += d
            if abs(d) < abs(s) * 1e-16:
                break
        return s * math.exp(lnpre)
    tiny = 1e-300                                           # continued fraction for Q, modified Lentz
    b = x + 1 - a
    c = 1 / tiny
    d = 1 / b
    h = d
    for i in range(1, 1000000):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = tiny if abs(d) < tiny else d
        c = b + an / c
        c = tiny if abs(c) < tiny else c
        d = 1 / d
        de = d * c
        h *= de
        if abs(de - 1) < 1e-16:
            break
    return 1.0 - math.exp(lnpre) * h

def chisq(x, k):
    """P(chi-square with k dof < x): the weight of an average term"""
    if k < 1:
        return 0.0
    if k < 3200:
        return igamma(k / 2.0, x / 2.0)
    return math.erf(math.sqrt(x) - math.sqrt(k)) / 2 + 0.5

def pnn(z, n):
    """probability that a deviate z among n is not noise"""
    return math.erf(abs(z) / math.sqrt(2)) ** n

def soft_pnna(z, n):
    """v1's approximation of pnn (still used for the probe CLASH term)"""
    a2, a1, a0, mx, my, y0 = -0.0192266, 0.751694, 1.12482, 0.21805, 0.736621, 1.0
    logn = 0.0 if n <= 0 else math.log10(n)
    asig = a2 * logn ** 2 + a1 * logn + a0
    expo = my * logn + mx * z + y0
    v = -(abs(z / asig) ** expo)
    v = max(-1000.0, min(1000.0, v))
    return 1 - 2.0 ** v

_TWO_M16 = 2.0 ** (-1.0 / 6.0)

def lj0(r, r0):
    if r == 0:
        return 1e40
    s = r0 * _TWO_M16 / r
    s6 = s * s * s * s * s * s
    return 4.0 * (s6 * s6 - s6)

def lj(r, r0):
    """the scorer's Lennard-Jones energy in well depths: minimum (-1) at r0 = phenix's vdw distance, 0 at 6 A"""
    return lj0(r, r0) - lj0(6.0, r0)

# -------------------------------------------------------------------------------------------------------------- parsers
def _quoted(line):
    w = line.split('"')
    return w[1] if len(w) > 1 else ''

def _atom(pid):
    """phenix pdb label ' CB  TYR A  14 ' -> (name, altloc, resname, chain, resnum) as the scorer's awk splits it"""
    a, f, t, c, r = pid[0:4], pid[4:5], pid[5:9], pid[9:10], pid[10:15]
    f = '_' if f in (' ', '') else f
    c = '_' if c in (' ', '') else c
    a_s, t_s, r_s = a.strip(), t.strip(), r.strip()
    if a_s and t_s and r_s and ' ' not in a_s and ' ' not in t_s and ' ' not in r_s:
        return (a_s, f, t_s, c, r_s)
    return tuple(('%4s %s %4s %s %5s' % (a, f, t, c, r)).split())

def _atom_has_digit(at):
    for tok in at:
        for ch in tok:
            if '0' <= ch <= '9':
                return True
    return False

def parse_geo(path, round_numbers=True, _fast_nb=False):
    """phenix .geo -> records dict(key, energy, dev, obs, ideal, sigma, atoms); NONBOND energy = lj(model, vdw)"""
    r6 = g6 if round_numbers else _identity
    lines = open(path).read().splitlines()
    recs = []
    i, n = 0, len(lines)
    atom_cache = {}

    def get_atom(pid):
        hit = atom_cache.get(pid)
        if hit is None:
            at = _atom(pid)
            hit = (at, _atom_has_digit(at))
            atom_cache[pid] = hit
        return hit

    def make_atoms(ids):
        out = [get_atom(ids[0])[0]]
        for p in ids[1:]:
            at, hd = get_atom(p)
            if hd:
                out.append(at)
        return out

    def nxt():
        nonlocal i
        if i + 1 < n:
            i += 1
        return lines[i]

    while i < n:
        line = lines[i]
        if 'pdb=' not in line:
            i += 1
            continue
        rec = None
        if 'nonbonded pdb=' in line:
            ids = [_quoted(line), _quoted(nxt())]
            nxt()
            v = nxt().split()
            obs, ideal = awknum(v[0]), awknum(v[1])
            if _fast_nb:
                if obs >= 0.92 * ideal and 0.0 < ideal <= 5.5:
                    e_val = -1.0
                else:
                    e_raw = lj(obs, ideal)
                    e_val = r6(e_raw) if e_raw > 0.0 else e_raw
                rec = dict(key='NONBOND', energy=e_val, dev=ideal - obs, obs=obs, ideal=ideal, sigma=1.0, ids=ids)
            else:
                rec = dict(key='NONBOND', energy=r6(lj(obs, ideal)), dev=r6(ideal - obs), obs=obs, ideal=ideal, sigma=1.0, ids=ids)
            line = lines[i]
        if 'bond pdb=' in line:
            ids = [_quoted(line), _quoted(nxt())]
            nxt()
            v = nxt().split()
            ideal, obs, sigma, e = awknum(v[0]), awknum(v[1]), awknum(v[3]), awknum(v[5])
            rec = dict(key='BOND', energy=r6(e), dev=r6(ideal - obs), obs=obs, ideal=ideal, sigma=r6(sigma), ids=ids)
            line = lines[i]
        if 'angle pdb=' in line:
            ids = [_quoted(line), _quoted(nxt()), _quoted(nxt())]
            nxt()
            v = nxt().split()
            ideal, obs, sigma, e = awknum(v[0]), awknum(v[1]), awknum(v[3]), awknum(v[5])
            rec = dict(key='ANGLE', energy=r6(e), dev=r6(ideal - obs), obs=obs, ideal=ideal, sigma=r6(sigma), ids=ids)
            line = lines[i]
        if 'dihedral pdb=' in line:
            ids = [_quoted(line), _quoted(nxt()), _quoted(nxt()), _quoted(nxt())]
            nxt()
            v = nxt().split()
            ideal, obs, sigma, e = awknum(v[0]), awknum(v[1]), awknum(v[4]), awknum(v[6])
            rec = dict(key='TORSION', energy=r6(e), dev=r6(ideal - obs), obs=obs, ideal=ideal, sigma=r6(sigma), ids=ids)
            line = lines[i]
        if 'chirality pdb=' in line:
            ids = [_quoted(line), _quoted(nxt()), _quoted(nxt()), _quoted(nxt())]
            nxt()
            v = nxt().split()
            ideal, obs, sigma, e = awknum(v[1]), awknum(v[2]), awknum(v[4]), awknum(v[6])
            rec = dict(key='CHIR', energy=r6(e), dev=r6(ideal - obs), obs=obs, ideal=ideal, sigma=r6(sigma), ids=ids)
            line = lines[i]
        if 'plane pdb=' in line:
            v = line.split()
            obs, sigma = awknum(v[-5]), awknum(v[-4])
            while line.split() and 'sigma' not in line:      # one record per plane atom, energy (delta/sigma)^2
                v = line.split()
                if 'plane' not in line:
                    obs, sigma = awknum(v[-3]), awknum(v[-2])
                recs.append(dict(key='PLANE', energy=r6((obs / sigma) ** 2), dev=r6(-obs), obs=obs, ideal=0.0,
                                 sigma=r6(sigma), atoms=[get_atom(_quoted(line))[0]]))
                if i + 1 >= n:
                    break
                line = nxt()
        if rec is not None:
            rec['atoms'] = make_atoms(rec.pop('ids'))
            recs.append(rec)
        i += 1
    return recs

def filter_geo(recs, ignore):
    """the scorer's ignore filters (default ARGH: drop H-H nonbonds between two ARG)"""
    ig = set(x.strip().upper() for x in re.split('[ ,]+', ignore or '') if x.strip())
    argh, hohh, hbon = 'ARGH' in ig, 'HOHH' in ig, 'HBON' in ig
    if not (argh or hohh or hbon):
        return list(recs)
    out = []
    for r in recs:
        a = r['atoms']
        a1 = a[0] if len(a) > 0 else ('', '', '', '', '')
        a2 = a[1] if len(a) > 1 else ('', '', '', '', '')
        if argh and r['key'] == 'NONBOND' and a1[0].startswith('H') and a2[0].startswith('H') \
                and a1[2] == 'ARG' and a2[2] == 'ARG':
            continue
        if hohh and (a1[0].startswith('H') or a2[0].startswith('H')) and (a1[2] == 'HOH' or a2[2] == 'HOH'):
            continue
        if hbon and r['key'] == 'BOND' and (a1[0].startswith('H') or a2[0].startswith('H')):
            continue
        out.append(r)
    return out

def parse_rama(path, round_numbers=True):
    """phenix.ramalyze log -> energies z^2 with p = score%/100 two-sided (lines of chain ' A' only, as the scorer)"""
    out = []
    for line in open(path):
        line = line.rstrip('\r\n')
        f = line.split(':')
        if line.startswith('SUMM') or line.startswith('resid') or not line.startswith(' A') or len(f) < 2:
            continue
        out.append(dict(energy=_prob_energy(1 - awknum(f[1]) / 100, round_numbers), label=line))
    return out

def parse_rota(path, round_numbers=True):
    """phenix.rotalyze log -> energies (residues with occupancy > 0)"""
    out = []
    for line in open(path):
        line = line.rstrip('\r\n')
        f = line.split(':')
        if line.startswith('SUMM') or line.startswith('resid') or len(f) <= 8 or not awknum(f[1]) > 0:
            continue
        out.append(dict(energy=_prob_energy(1 - awknum(f[2]) / 100, round_numbers), label=line))
    return out

def _prob_energy(P, round_numbers=True):
    """P = 1 - probability-OK; energy = 2 inverf(P)^2 = z^2 with p = erfc(z/sqrt2); 99 if P prints as 1"""
    if round_numbers:
        if float('%g' % P) == 1:
            return 99.0
        return g6(g15(2 * inverf(P - 1e-16) ** 2))
    if P >= 1.0:
        return 99.0
    if P <= 0.0:
        return 0.0
    return 2.0 * inverf(P) ** 2

def _cbdev_energy(dev, cbsigma=0.032, cbdof=3, round_numbers=True):
    r6 = g6 if round_numbers else _identity
    r15 = g15 if round_numbers else _identity
    x = r6(dev / cbsigma)
    if x > 37:
        e = x * x
    else:
        Q = math.exp(-x ** 2 / 2.0) if cbdof == 2 else \
            math.erfc(x / math.sqrt(2.0)) + math.sqrt(2.0 / math.pi) * x * math.exp(-x ** 2 / 2.0)
        e = ndtri(Q / 2.0) ** 2
    return r6(r15(e))

def parse_cbetadev(path, cbsigma=0.032, cbdof=3, round_numbers=True):
    """phenix.cbetadev log -> energies: x = dev/cbsigma is a chi(cbdof) deviate, energy = Gaussian z^2 of equal tail"""
    cbdof = int(cbdof)
    if cbdof not in (2, 3):
        raise ValueError('invalid cbdof=%r (expected 2 or 3)' % (cbdof,))
    out = []
    for line in open(path):
        line = line.rstrip('\r\n')
        f = line.split(':')
        if line.startswith('SUMM') or line.startswith('pdb') or len(f) <= 5:
            continue
        out.append(dict(energy=_cbdev_energy(awknum(f[5]), cbsigma, cbdof, round_numbers), label=line))
    return out

def parse_coot_clashes(path, round_numbers=True):
    """molprobity_coot.py (phenix.molprobity output.coot=True) -> probe clashes, energy lj(3 + overlap, 3)"""
    if not path or not os.path.exists(path):
        return []
    r6 = g6 if round_numbers else _identity
    text = '\n'.join(open(path).read().splitlines()[-2:]) + '\n'   # tail -n 2
    out, atoms = [], ''
    for rec in text.split("'"):
        if len(rec) == 16:
            atoms += '|' + rec
        if rec.startswith(',') and len(rec.split()) > 1:
            ov = awknum(rec.replace(',', '').split()[0])
            out.append(dict(energy=r6(lj(3 + ov, 3)), dev=r6(-ov), atoms=atoms))
            atoms = ''
    return out

def _parse_occupancies_lines(lines, filter_xp=True):
    occ, tot = {}, 0.0
    for line in lines:
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue
        if filter_xp and line.split()[-1:] == ['XP']:
            continue
        nm = line[12:16].replace(' ', '')
        al = line[16:17] if line[16:17] not in (' ', '') else '_'
        rn = line[17:20].replace(' ', '')
        ch = line[21:22] if line[21:22] not in (' ', '') else '_'
        rs = int(awknum(line[22:26]))
        o = awknum(line[54:60])
        occ[(nm, al, rn, ch, rs)] = o
        tot += o
    return occ, tot

def read_occupancies(pdb):
    """(name, altloc, resname, chain, resnum) -> occupancy, and the total occupancy of all non-XP ATOM/HETATM records"""
    with open(pdb) as f:
        return _parse_occupancies_lines(f, filter_xp=True)

# -------------------------------------------------------------------------------------------------------------- scoring
def score_items(items, averages, clashes=(), nb_phys=0.0, opts=None):
    """The v2 sum.
    items:    {category: [energy z^2 of every item, zeros and negatives included]}; each item with energy E > 0 adds
              Pnn(sqrt(min(E,400)), N) * soft(E), N = len(items[category]).
    averages: {name: (mean, n)}; each adds chisq(min(max(mean,0),400) * n, n) * soft(max(mean,0)) (+ n for 'clash').
    clashes:  probe clash LJ energies; adds softPnna(sqrt(min(1 + worst*n, 400)), n) * (soft(1 + worst*n) + n).
    nb_phys:  added as is (nbmode physical/both).
    Returns dict(wE, items=[(value, category, k, N, Pnn, energy)], weighted=[(value, name)])."""
    o = _normalize_opts(opts)
    rnd = bool(o['round'])
    r6 = g6 if rnd else _identity
    r15 = g15 if rnd else _identity
    out_items = []
    for c, es in items.items():
        N = len(es)
        for k, E in enumerate(sorted(es, reverse=True), 1):
            if E <= 0:
                continue
            Ec = min(E, 400.0)
            z = r6(math.sqrt(Ec))
            P = pnn(z, N) if o['pnn'] == 'exact' else soft_pnna(z, N)
            if P == 0.0:                      # gnuplot drops an underflowed result ("undefined value"): item not listed
                continue
            P = r15(P)
            if rnd and P < 1e-30:
                P = 0.0
            out_items.append((r6(P * soft(E)), c, k, N, r6(P), E))
    weighted = []
    for name, (mean, n) in averages.items():
        e = 0.0 if (mean is None or mean <= 0) else min(mean, 400.0)
        w = r15(chisq(r6(e * n), n))
        ee = 0.0 if (mean is None or mean <= 0) else soft(mean)
        if name.lower() == 'clash':
            ee += n
        weighted.append((r6(ee * w), 'avg_' + name))
    if clashes:
        n = len(clashes)
        worst = max(clashes)
        raw_e1 = max(0.0, 1.0 + worst * n)
        e1 = g6(math.sqrt(raw_e1) ** 2) if rnd else raw_e1
        w = r15(soft_pnna(r6(math.sqrt(min(e1, 400.0))), n))
        weighted.append((r6((soft(e1) + n) * w), 'clash'))
    wE = r6(sum(v for v, *_ in out_items) + sum(v for v, _ in weighted) + nb_phys)
    out_items.sort(key=lambda x: (x[0], x[4], x[5], x[1], -x[2]), reverse=True)
    return dict(wE=wE, items=out_items, weighted=weighted)

def _mean(vals, r6=g6):
    return (r6(sum(vals) / len(vals)), len(vals)) if vals else (None, 0)

def _score_parsed(full, ramas, rotas, cbs, clashes, get_occ_fn, o, ovr=None):
    rnd = bool(o['round'])
    r6 = g6 if rnd else _identity
    r15 = g15 if rnd else _identity
    filt = filter_geo(full, o['ignore'])
    geo_noH = [r for r in filt if r['key'] == 'NONBOND' or not any(a[0].startswith('H') for a in r['atoms'])]
    if ovr:
        geo_noH = override_records(geo_noH, ovr, rnd)   # the scorer's ${t}_geo.txt; nbmode=occ still reads ${t}_filtgeo.txt
    by = {c: [] for c in ('BOND', 'ANGLE', 'TORSION', 'CHIR', 'PLANE')}
    nonbond = []
    for r in geo_noH:
        k = r['key']
        if k == 'NONBOND':
            nonbond.append(r)
        elif k in by:
            by[k].append(r['energy'])

    # omegas: every CA-x-x-CA torsion of the full .geo, divided by 3 (5) next to one (two) prolines
    pro = set()
    for r in full:
        for a in r['atoms']:
            if len(a) >= 5 and ('PRO' in a[2] or 'PRO' in a[3] or 'PRO' in a[4]):
                pro.add(int(awknum(a[4])))
    sig = math.radians(o['sigomega'])
    sin_sig = math.sin(sig)
    omegas = []
    for r in full:
        if r['key'] == 'TORSION':
            a = r['atoms']
            if len(a) >= 4 and a[0][0] == 'CA' and a[3][0] == 'CA':
                n = int(math.fmod(int(awknum(a[0][4])) - 1, 9999)) + 1
                if ovr and str(n) in ovr['skip']['omega']:
                    continue
                prox = (n - 1 in pro) + (n + 1 in pro)
                om = math.radians(r['obs'])
                omegas.append(r6(((math.sin(om) / sin_sig) ** 2 + (1 + math.cos(om)) ** 10) / (prox * 2 + 1)))

    averages = dict(bond=_mean(by['BOND'], r6), clash=_mean(clashes, r6), angle=_mean(by['ANGLE'], r6),
                    torsion=_mean(by['TORSION'], r6), chiral=_mean(by['CHIR'], r6), plane=_mean(by['PLANE'], r6),
                    omega=_mean(omegas, r6), rota=_mean(rotas, r6), rama=_mean(ramas, r6), CBdev=_mean(cbs, r6))
    if averages['clash'][0] is None:
        averages['clash'] = (None, 0)

    # nonbonds
    pos = [r['energy'] for r in nonbond if r['energy'] > 0]
    nb_phys = (r6(sum(soft(e) for e in pos)), len(pos))
    nbitems = []
    nbinfo = {}
    if o['nbmode'] == 'bonds':
        nbonds = len(by['BOND'])
        nbitems = pos + [0.0] * max(0, nbonds - len(pos))
        if nbonds > 0:
            averages['nonbond'] = (r6(sum(pos) / nbonds), nbonds)
        nbinfo = dict(N=nbonds)
    elif o['nbmode'] == 'occ':
        occ, tot = get_occ_fn()
        w_lj = []
        s12 = s10 = 0.0
        resnum_cache = {}
        for r in filt:
            if r['key'] != 'NONBOND':
                continue
            a, b = r['atoms'][0], r['atoms'][1]
            ra_s, rb_s = a[4], b[4]
            ra_i = resnum_cache.get(ra_s)
            if ra_i is None:
                ra_i = int(awknum(ra_s))
                resnum_cache[ra_s] = ra_i
            rb_i = resnum_cache.get(rb_s)
            if rb_i is None:
                rb_i = int(awknum(rb_s))
                resnum_cache[rb_s] = rb_i
            ka = (a[0], a[1], a[2], a[3], ra_i)
            kb = (b[0], b[1], b[2], b[3], rb_i)
            oa, ob = occ.get(ka, 1.0), occ.get(kb, 1.0)
            w = oa if oa < ob else ob
            ideal = r['ideal']
            ratio = r6(r['obs'] / ideal) if ideal > 0 else 9
            if ratio < 1.2:
                s12 += w
            if ratio < 1.0:
                s10 += w
            if r['energy'] > 0:
                w_lj.append((r6(w), r['energy']))
        N = r6({'vdw120': s12, 'vdw100': s10}.get(o['nbN'], tot))
        Nint = int(N + 0.5)
        if N > 0:
            averages['nonbond'] = (r6(sum(w * e for w, e in w_lj) / N), N)
        posw = [(r6(w * e) if o['nboutw'] else e) for w, e in w_lj]
        nbitems = posw + [0.0] * max(0, Nint - len(posw))
        nbinfo = dict(N=N, total_occupancy=r6(tot))
    elif o['nbmode'] in ('stat', 'both'):
        nbitems = [0.0] * (len(nonbond) - len(pos))
        for e_in in pos:
            E = r6(min(e_in, 400.0))
            nbitems.append(r6(r15(ndtri(o['nbP0'] * math.exp(-o['nbk'] * E) / 2.0) ** 2)))

    items = dict(BOND=by['BOND'], ANGLE=by['ANGLE'], TORSION=by['TORSION'], CHIR=by['CHIR'], PLANE=by['PLANE'],
                 OMEGA=omegas, RAMA=ramas, ROTA=rotas, CBDEV=cbs)
    if nbitems:
        items['NONBOND'] = nbitems
    add = nb_phys[0] if o['nbmode'] in ('physical', 'both') else 0.0
    res = score_items(items, averages, clashes, add, o)

    terms = {}
    for v, c, k, N, P, E in res['items']:
        t = terms.setdefault(c, dict(outliers=0.0, n_gt_001=0, N=N, average=0.0))
        t['outliers'] += v
        t['n_gt_001'] += v > 0.01
    amap = dict(avg_bond='BOND', avg_clash='CLASH', avg_angle='ANGLE', avg_torsion='TORSION', avg_chiral='CHIR',
                avg_plane='PLANE', avg_omega='OMEGA', avg_rota='ROTA', avg_rama='RAMA', avg_cbdev='CBDEV',
                avg_nonbond='NONBOND', clash='CLASH')
    for v, name in res['weighted']:
        terms.setdefault(amap[name.lower()], dict(outliers=None, n_gt_001=None, N=None, average=0.0))['average'] += v
    res.update(terms=terms, nb_phys=nb_phys, nonbond=nbinfo, n_clashes=len(clashes),
               options={k: o[k] for k in DEFAULTS})
    return res

def score_files(pdb, geo, rama, rota, cbdev, coot=None, fullgeo=None, clashes=None, cbdevlog=None, outprefix=None,
                **opts):
    """Score existing phenix outputs exactly as molprobify_v2.com does.
    fullgeo / clashes / cbdevlog / outprefix: also write counter-fit's files (see write_counterfit_files)."""
    o = _normalize_opts(opts)
    rnd = bool(o['round'])
    ovr = read_overrides(o['overridefile'])
    full = parse_geo(geo, rnd, _fast_nb=True)
    ramas = [x['energy'] for x in parse_rama(rama, rnd)
             if not (ovr and _resnum_mod(x['label'], 64) in ovr['skip']['rama'])]
    rotas = [x['energy'] for x in parse_rota(rota, rnd)
             if not (ovr and _resnum_mod(x['label'], 9999) in ovr['skip']['rota'])]
    cbs = [x['energy'] for x in parse_cbetadev(cbdev, o['cbsigma'], o['cbdof'], rnd)]
    probe = [x['energy'] for x in parse_coot_clashes(coot, rnd)]
    res = _score_parsed(full, ramas, rotas, cbs, probe, lambda: read_occupancies(pdb), o, ovr)
    write_counterfit_files(geo, coot, cbdev, o['overridefile'], fullgeo=fullgeo, clashes=clashes, cbdevlog=cbdevlog,
                           outprefix=outprefix)
    return res

# ------------------------------------------------------------------------- counter-fit files: fullgeo, clashes, overrides
# Ports of the awk in counter-fit's code_v2/molprobify_runme.com (md5 457d87ea): <prefix>_fullgeo.txt (one line per .geo
# restraint, H kept, after overrides, " symop=<op>" on symmetry-mate contacts), <prefix>_clashes.txt (probe clashes) and the
# raw phenix.cbetadev log, which counter-fit's fudgeify.py reads; and the override file, which changes the score.
def _awk_out(x):
    """a number as awk's print writes it: integral values as integers, others by OFMT %.6g; strings as they are"""
    if isinstance(x, str):
        return x
    if x == x and x not in (float('inf'), float('-inf')) and x == int(x) and abs(x) < 9.2e18:
        return '%d' % int(x)
    return '%.6g' % x

def _awk_pow(x, n):
    """gawk x^n: an integral exponent by repeated squaring (calc_exp_posint), otherwise pow()"""
    if n != int(n):
        return math.pow(x, n)
    n = int(n)
    if n == 0:
        return 1.0
    neg, n, mult = n < 0, abs(n), 1.0
    while n > 1:
        if n % 2 == 1:
            mult *= x
        x *= x
        n //= 2
    return 1.0 / (mult * x) if neg else mult * x

def _awk_lj(r, r0):
    """the scorer's lj() evaluated as gawk does it, for byte-identical printed energies"""
    def lj0(r, r0):
        if r == 0:
            return 1e40
        s = r0 * math.pow(2.0, -1.0 / 6) / r
        return 4 * (_awk_pow(s, 12) - _awk_pow(s, 6))
    return lj0(r, r0) - lj0(6, r0)

_NUM_PREFIX = re.compile(r'\s*[-+]?(inf(inity)?|nan|(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?)', re.I)

def _sort_k2gr(lines):
    """GNU sort -k2gr in the C locale: field 2 as a general number, largest first (non-numbers last); ties by the whole line,
    ascending, in ASCII byte order (the last-resort comparison does not take the key's r)"""
    def key(line):
        f = line.split()
        m = _NUM_PREFIX.match(f[1]) if len(f) > 1 else None
        v = float(m.group(0)) if m else float('-inf')
        return (0 if m else 1, -v if m else 0.0, line)
    return sorted(lines, key=key)

def _resnum_mod(label, modulo):
    """the scorer's residue tag of a ramalyze/rotalyze line: (substr($0,3,4)-1)%modulo+1, as a string"""
    return str(int(math.fmod(int(awknum(label[2:6])) - 1, modulo)) + 1)

def fullgeo_lines_text(geo_text):
    """phenix .geo text -> the scorer's fullgeo lines (unsorted): key energy delta obs ideal sigma | atom alt resname chain
    resnum - atom ...  (code_v2/molprobify_runme.com lines 352-411, including its symop= suffix)"""
    lines = geo_text.splitlines()
    out = []
    n = len(lines)
    st = dict(i=0)

    def nxt():                      # awk getline: at end of file $0 stays as it was
        if st['i'] + 1 < n:
            st['i'] += 1
        return lines[st['i']]

    def fld(line, k):               # awk $k (1-based; negative k counts from NF)
        f = line.split()
        if k < 0:
            k = len(f) + 1 + k
        return f[k - 1] if 1 <= k <= len(f) else ''

    def nid(pid):
        a, f, t, c, r = pid[0:4], pid[4:5], pid[5:9], pid[9:10], pid[10:15]
        f, c = f.replace(' ', '_'), c.replace(' ', '_')
        return '%4s %s %4s %s %5s' % (a, f, t, c, r)

    key = obs = ideal = sigma = ''
    energy = 0.0
    ids_raw, sym = [], ''
    while st['i'] < n:
        line = lines[st['i']]
        if 'nonbonded pdb=' in line:
            key = 'NONBOND'
            ids_raw = [_quoted(line), _quoted(nxt())]
            nxt()
            line = nxt()
            obs, ideal, sigma = fld(line, 1), fld(line, 2), 1.0
            energy = _awk_lj(awknum(obs), awknum(ideal))
            if ',' in fld(line, 3):
                sym = fld(line, 3)
        if 'bond pdb=' in line:
            key = 'BOND'
            ids_raw = [_quoted(line), _quoted(nxt())]
            nxt()
            line = nxt()
            ideal, obs, sigma, energy = fld(line, 1), fld(line, 2), fld(line, 4), fld(line, 6)
            if ',' in fld(line, -1):
                sym = fld(line, -1)
        if 'angle pdb=' in line:
            key = 'ANGLE'
            ids_raw = [_quoted(line), _quoted(nxt()), _quoted(nxt())]
            nxt()
            line = nxt()
            ideal, obs, sigma, energy = fld(line, 1), fld(line, 2), fld(line, 4), fld(line, 6)
        if 'dihedral pdb=' in line:
            key = 'TORSION'
            ids_raw = [_quoted(line), _quoted(nxt()), _quoted(nxt()), _quoted(nxt())]
            nxt()
            line = nxt()
            ideal, obs, sigma, energy = fld(line, 1), fld(line, 2), fld(line, 5), fld(line, 7)
        if 'chirality pdb=' in line:
            key = 'CHIR'
            ids_raw = [_quoted(line), _quoted(nxt()), _quoted(nxt()), _quoted(nxt())]
            nxt()
            line = nxt()
            ideal, obs, sigma, energy = fld(line, 2), fld(line, 3), fld(line, 5), fld(line, 7)
        if 'plane pdb=' in line:
            key, ideal = 'PLANE', 0.0
            obs, sigma = fld(line, -5), fld(line, -4)
            while line.split() and 'sigma' not in line:
                if 'plane' not in line:
                    obs, sigma = fld(line, -3), fld(line, -2)
                o_, s_ = awknum(obs), awknum(sigma)
                energy = _awk_pow(o_ / s_, 2)
                out.append(' '.join(['PLANE', _awk_out(energy), _awk_out(0.0 - o_), obs, '0', _awk_out(s_), '|',
                                     nid(_quoted(line))]))
                if st['i'] + 1 >= n:
                    break
                line = nxt()
        if ids_raw:
            ns = [nid(p) for p in ids_raw]
            ids = ns[0]
            for x in ns[1:]:
                if any('0' <= ch <= '9' for ch in x):
                    ids += ' - ' + x
            if sym != '' and sym != 'x,y,z':
                ids += ' symop=' + sym
            out.append(' '.join([key, _awk_out(_num(energy)), _awk_out(_num(ideal) - _num(obs)), obs, ideal,
                                 _awk_out(_num(sigma)), '|', ids]))
            ids_raw, sym = [], ''
        st['i'] += 1
    return out

def fullgeo_lines(geo_path):
    """phenix .geo file -> the scorer's fullgeo lines (unsorted)"""
    with open(geo_path) as f:
        return fullgeo_lines_text(f.read())

def read_overrides(path):
    """the scorer's override file -> dict(lines=raw lines, by_full={ent: [(ideal, sigma)]}, by_geo={...}, skip=...);
    None when path is empty or missing (the scorer then ignores it)"""
    if not path or not os.path.exists(path):
        return None
    raw = open(path).read().split('\n')
    if raw and raw[-1] == '':
        raw.pop()
    by_full, by_geo = {}, {}
    skip = dict(rama=set(), rota=set(), omega=set())
    for line in raw:
        f = line.split()
        if len(f) < 2 or f[1] != 'OVERRIDE':
            continue
        after = line.split('|')[1] if '|' in line else ''
        t, s = (f[2] if len(f) > 2 else ''), (f[3] if len(f) > 3 else '')
        by_full.setdefault(f[0] + ' ' + (after[1:] if after.startswith(' ') else after), []).append((t, s))
        by_geo.setdefault(f[0] + ' ' + after.lstrip(' '), []).append((t, s))
    for line in raw:                # "edit other tables": TORSION overrides drop that residue's rama/rota/omega items
        if not line.startswith('TORS'):
            continue
        atoms = line[line.index('|') + 1:] if '|' in line else line
        a = atoms.split()
        w = a[0].split('_') if a else ['']
        resn = w[1] if len(w) > 1 else ''
        names = re.sub('[ 0-9]', '', atoms).replace('_', ' ').split()
        typ = 'rama'
        if any(x not in ('CA', 'N', 'C') for x in names[0:4]) or len(names) < 4:
            typ = 'rota'
        if names and names[0] == 'CA' and names[-1] == 'CA':
            typ = 'omega'
        skip[typ].add(resn)
    return dict(lines=raw, by_full=by_full, by_geo=by_geo, skip=skip)

def _num(x):
    return x if isinstance(x, float) else awknum(x)

def _apply_override(key, energy, obs, ideal, sigma, ovs):
    """the scorer's override loop: each (ideal, sigma) that lowers the energy replaces it; returns (energy, delta, ideal,
    sigma, changed) with energy/delta numbers when changed (inputs: awk fields as strings, or numbers)"""
    changed, delta = False, None
    for t, s in ovs:
        t = ideal if t == '' else t
        s = sigma if s == '' else s
        dev = _num(t) - _num(obs)
        newE = _awk_pow(dev / _num(s), 2) + _num(energy) / 1e6
        if key == 'NONBOND':
            if _num(t) > _num(obs):
                t = ideal
            newE = _awk_lj(_num(obs), _num(t))
        if newE < _num(energy):
            ideal, energy, delta, sigma, changed = t, newE, dev, s, True
    return energy, delta, ideal, sigma, changed

def override_fullgeo(lines, ovr):
    """the scorer's second override pass, on fullgeo lines (override-file lines that are not OVERRIDE pass through)"""
    if not ovr:
        return list(lines)
    out = []
    src = [l for l in ovr['lines'] if not (len(l.split()) > 1 and l.split()[1] == 'OVERRIDE')] + list(lines)
    for line in src:
        f = line.split()
        f += [''] * (6 - len(f))
        atoms = line.split('|')[1] if '|' in line else ''
        if atoms.startswith(' '):
            atoms = atoms[1:]
        w = atoms.split()
        ent = f[0] + ' ' + ' '.join(w[i] + '_' + (w[i + 4] if i + 4 < len(w) else '') for i in range(0, len(w) - 1, 6))
        ovs = ovr['by_full'].get(ent)
        if not ovs:
            out.append(line)
            continue
        energy, delta, ideal, sigma, changed = _apply_override(f[0], f[1], f[3], f[4], f[5], ovs)
        if changed:
            out.append(' '.join([f[0], _awk_out(energy), _awk_out(delta), f[3], ideal, sigma, '|', atoms]))
        else:
            out.append(' '.join([f[0], f[1], f[2], f[3], f[4], f[5], '|', atoms]))
    return out

def override_records(recs, ovr, round_numbers=True):
    """the scorer's first override pass, on ${t}_geo.txt (scored energies): records whose key + name_resnum atoms match an
    override get the lowest energy it allows"""
    r6 = g6 if round_numbers else _identity
    out = []
    for r in recs:
        ent = r['key'] + ' ' + ' '.join(a[0] + '_' + a[4] for a in r['atoms'])
        ovs = ovr['by_geo'].get(ent)
        if not ovs:
            out.append(r)
            continue
        e = r['energy']
        if r['key'] == 'NONBOND':
            e = r6(lj(r['obs'], r['ideal']))
        energy, _, _, _, changed = _apply_override(r['key'], float(e), float(r['obs']), float(r['ideal']),
                                                   float(r['sigma']), ovs)
        if changed:
            r = dict(r, energy=r6(energy))
        out.append(r)
    return out

def clash_lines_text(tail_text):
    """last lines of molprobity_coot.py -> the scorer's clashes lines (unsorted): CLASH ljenergy -overlap |atom16|atom16"""
    out, atoms = [], ''
    for line in tail_text.replace("'", '\n').split('\n'):                   # awk RS="[']" then one record per line
        if len(line) == 16:
            atoms += '|' + line
        if line.startswith(',') and len(line.split()) > 1:
            v = awknum(line.replace(',', '').split()[0])
            out.append(' '.join(['CLASH', _awk_out(_awk_lj(3 + v, 3)), _awk_out(-v), atoms]))
            atoms = ''
    return out

def clash_lines(coot_path):
    """molprobity_coot.py -> the scorer's clashes lines (unsorted): CLASH ljenergy -overlap |atom16|atom16"""
    if not coot_path or not os.path.exists(coot_path):
        return []
    with open(coot_path) as f:
        text = '\n'.join(f.read().splitlines()[-2:]) + '\n'                 # tail -n 2
    return clash_lines_text(text)

def _write_fullgeo_text(geo_text, ovr, out_path):
    lines = _sort_k2gr(fullgeo_lines_text(geo_text))
    if ovr:
        lines = _sort_k2gr(override_fullgeo(lines, ovr))
    with open(out_path, 'w') as f:
        f.write(''.join(l + '\n' for l in lines))

def write_counterfit_files(geo, coot, cbdev, overridefile='', fullgeo=None, clashes=None, cbdevlog=None, outprefix=None):
    """Write counter-fit's scorer files. outprefix=P is shorthand for P_fullgeo.txt, P_clashes.txt and P_cbetadev.log."""
    if outprefix:
        fullgeo = fullgeo or outprefix + '_fullgeo.txt'
        clashes = clashes or outprefix + '_clashes.txt'
        cbdevlog = cbdevlog or outprefix + '_cbetadev.log'
    if fullgeo:
        with open(geo) as gf:
            _write_fullgeo_text(gf.read(), read_overrides(overridefile), fullgeo)
    if clashes:
        with open(clashes, 'w') as f:
            f.write(''.join(l + '\n' for l in _sort_k2gr(clash_lines(coot))))
    if cbdevlog and os.path.abspath(cbdevlog) != os.path.abspath(cbdev):
        shutil.copyfile(cbdev, cbdevlog)

# ------------------------------------------------------------------------------------------------------- mmtbx backend
_MMTBX = None

def _load_mmtbx():
    """Lazy one-time import of cctbx/mmtbx modules when running under phenix.python; returns module bundle or None."""
    global _MMTBX
    if _MMTBX is not None:
        return _MMTBX if _MMTBX is not False else None
    try:
        import iotbx.pdb
        import iotbx.cif
        from cctbx import uctbx
        from scitbx.array_family import flex
        import mmtbx.model
        import mmtbx.command_line.geometry_minimization as gm_cli
        import mmtbx.conformation_dependent_library as cdl_mod
        from libtbx.utils import null_out
        from mmtbx.validation.ramalyze import ramalyze
        from mmtbx.validation.rotalyze import rotalyze
        from mmtbx.validation.cbetadev import cbetadev
        from mmtbx.validation.clashscore import clashscore
        from mmtbx.conformation_dependent_library.cdl_setup import setup_restraints as cdl_setup
        from cctbx.geometry_restraints.linking_class import linking_class
        _MMTBX = dict(
            pdb=iotbx.pdb, cif=iotbx.cif, uctbx=uctbx, flex=flex, model=mmtbx.model,
            gm_cli=gm_cli, gm_master=gm_cli.master_params(), cdl_mod=cdl_mod, null_out=null_out,
            ramalyze=ramalyze, rotalyze=rotalyze, cbetadev=cbetadev, clashscore=clashscore,
            cdl_setup=cdl_setup, origin_ids=linking_class(),
        )
        return _MMTBX
    except ImportError:
        _MMTBX = False
        return None

@contextlib.contextmanager
def _cdl_interp_context(M, enabled):
    if not enabled:
        yield
        return
    cdl_mod = M['cdl_mod']
    orig = cdl_mod.update_restraints
    def _interp_update(*args, **kwds):
        kwds['interpolate'] = True
        return orig(*args, **kwds)
    cdl_mod.update_restraints = _interp_update
    try:
        yield
    finally:
        cdl_mod.update_restraints = orig

def _extract_restraint_records(grm, sites_cart, site_labels, origin_ids, round_numbers=True):
    """Extract BOND, ANGLE, TORSION, CHIR, PLANE and NONBOND records directly from mmtbx restraints_manager."""
    r6 = g6 if round_numbers else _identity
    recs = []
    atom_cache = {}

    def get_atom_lbl(lbl):
        hit = atom_cache.get(lbl)
        if hit is None:
            at = _atom(lbl[5:-1])
            hit = (at, _atom_has_digit(at))
            atom_cache[lbl] = hit
        return hit

    def make_atoms_lbls(lbls):
        out = [get_atom_lbl(lbls[0])[0]]
        for lbl in lbls[1:]:
            at, hd = get_atom_lbl(lbl)
            if hd:
                out.append(at)
        return out

    atom_tuples = [_atom(lbl[5:-1]) for lbl in site_labels]
    atom_has_dig = [_atom_has_digit(at) for at in atom_tuples]

    pair_proxies = grm.pair_proxies(sites_cart=sites_cart)
    nbp = pair_proxies.nonbonded_proxies
    if nbp is not None:
        sorted_table, _ = nbp.get_sorted(
            by_value='delta', sites_cart=sites_cart, site_labels=None,
            max_items=None, include_proxy=False)
        for _, i_seq, j_seq, delta, vdw_distance, _, _ in sorted_table:
            obs = round(delta, 3) if round_numbers else delta
            ideal = round(vdw_distance, 3) if round_numbers else vdw_distance
            if obs >= 0.92 * ideal and 0.0 < ideal <= 5.5:
                e_val = -1.0
            else:
                e_raw = lj(obs, ideal)
                e_val = r6(e_raw) if e_raw > 0.0 else e_raw
            a1 = atom_tuples[i_seq]
            recs.append(dict(key='NONBOND', energy=e_val, dev=ideal - obs, obs=obs, ideal=ideal,
                             sigma=1.0, atoms=[a1, atom_tuples[j_seq]] if atom_has_dig[j_seq] else [a1]))

    bond_oids = [origin_ids.get_origin_id(k) for k in origin_ids.get_bond_origin_id_labels()]
    for pname, key in (('bond_proxies', 'BOND'), ('angle_proxies', 'ANGLE'),
                       ('dihedral_proxies', 'TORSION'), ('chirality_proxies', 'CHIR'),
                       ('planarity_proxies', 'PLANE')):
        proxies = getattr(pair_proxies, pname, getattr(grm, pname, None))
        if proxies is None:
            continue
        oids = bond_oids if pname in ('bond_proxies', 'angle_proxies') else (None,)
        for oid in oids:
            r = proxies.show_sorted(by_value='residual', sites_cart=sites_cart, site_labels=site_labels,
                                    f=_MMTBX['null_out'](), origin_id=oid, return_result=True)
            if not r:
                continue
            if key == 'PLANE':
                for v in r.value_list:
                    obs = round(v.delta, 3) if round_numbers else v.delta
                    sig = float('%.2e' % v.sigma) if round_numbers else v.sigma
                    recs.append(dict(key='PLANE', energy=r6((obs / sig) ** 2), dev=r6(-obs), obs=obs,
                                     ideal=0.0, sigma=r6(sig), atoms=[get_atom_lbl(v.labels[0])[0]]))
            elif key == 'TORSION':
                for v in r.value_list:
                    obs = round(v.model, 2) if round_numbers else v.model
                    ideal = round(v.ideal, 2) if round_numbers else v.ideal
                    sig = float('%.2e' % v.sigma) if round_numbers else v.sigma
                    e = float('%.2e' % v.residual) if round_numbers else v.residual
                    recs.append(dict(key='TORSION', energy=r6(e), dev=ideal - obs, obs=obs, ideal=ideal,
                                     sigma=r6(sig), atoms=make_atoms_lbls(v.labels)))
            else:
                dec = 3 if key == 'BOND' else 2
                for v in r.value_list:
                    obs = round(v.model, dec) if round_numbers else v.model
                    ideal = round(v.ideal, dec) if round_numbers else v.ideal
                    sig = float('%.2e' % v.sigma) if round_numbers else v.sigma
                    e = float('%.2e' % v.residual) if round_numbers else v.residual
                    recs.append(dict(key=key, energy=r6(e), dev=ideal - obs, obs=obs, ideal=ideal,
                                     sigma=r6(sig), atoms=make_atoms_lbls(v.labels)))
    return recs

def _extract_validation_lists(M, raw_ph, cl_ph, o, ovr=None, include_clashes=True, log=sys.stderr):
    """Run ramalyze, rotalyze, cbetadev on raw_ph and (optionally) clashscore on cl_ph in-process."""
    rnd = bool(o['round'])
    r6 = g6 if rnd else _identity
    null = M['null_out']()

    ra_res = M['ramalyze'](pdb_hierarchy=raw_ph, outliers_only=False, out=null)
    ramas = []
    for r in ra_res.results:
        old = r.format_old()
        if not old.startswith(' A'):
            continue
        if ovr and _resnum_mod(old, 64) in ovr['skip']['rama']:
            continue
        p_ok = (awknum(old.split(':')[1]) / 100.0) if rnd else (r.score / 100.0)
        ramas.append(_prob_energy(1.0 - p_ok, rnd))

    ro_res = M['rotalyze'](pdb_hierarchy=raw_ph, outliers_only=False, out=null)
    rotas = []
    for r in ro_res.results:
        old = r.format_old() if (rnd or ovr) else None
        if ovr and _resnum_mod(old, 9999) in ovr['skip']['rota']:
            continue
        if rnd:
            f = old.split(':')
            if len(f) <= 8 or not awknum(f[1]) > 0:
                continue
            p_ok = awknum(f[2]) / 100.0
        else:
            if not r.occupancy > 0:
                continue
            p_ok = r.score / 100.0
        rotas.append(_prob_energy(1.0 - p_ok, rnd))

    cb_res = M['cbetadev'](pdb_hierarchy=raw_ph, outliers_only=False, out=null)
    cbs = [_cbdev_energy(round(r.deviation, 3) if rnd else r.deviation, o['cbsigma'], o['cbdof'], rnd)
           for r in cb_res.results]

    clashes, cl_res = [], None
    if include_clashes:
        try:
            cl_res = M['clashscore'](pdb_hierarchy=cl_ph, fast=False, condensed_probe=False,
                                     keep_hydrogens=False, nuclear=False)
            for r in cl_res.results:
                if r.is_outlier():
                    ov = round(r.overlap, 3) if rnd else r.overlap
                    clashes.append(r6(lj(3.0 + ov, 3.0)))
        except Exception as exc:
            if log:
                log.write('WARNING: in-process clashscore failed (%s); scoring with 0 probe clashes\n' % exc)

    return ramas, rotas, cbs, clashes, cb_res, cl_res

def _write_counterfit_mmtbx(pdb, model, cl_res, cb_res, ovr, fullgeo=None, clashes=None, cbdevlog=None, outprefix=None):
    if outprefix:
        fullgeo = fullgeo or outprefix + '_fullgeo.txt'
        clashes = clashes or outprefix + '_clashes.txt'
        cbdevlog = cbdevlog or outprefix + '_cbetadev.log'
    if fullgeo:
        geo_text = '# Geometry restraints after refinement\n' + model.restraints_as_geo()
        _write_fullgeo_text(geo_text, ovr, fullgeo)
    if clashes:
        probe_data = cl_res.as_coot_data() if cl_res is not None else []
        tail = "data['probe'] = %s\ngui = coot_molprobity_todo_list_gui(data=data)\n" % probe_data
        with open(clashes, 'w') as f:
            f.write(''.join(l + '\n' for l in _sort_k2gr(clash_lines_text(tail))))
    if cbdevlog:
        prefix = os.path.splitext(os.path.basename(pdb))[0]
        with open(cbdevlog, 'w') as f:
            cb_res.show_old_output(out=f, verbose=True, prefix=prefix)

def _parse_cif_objects(M, cifs):
    return [(os.path.abspath(c), M['cif'].reader(file_path=os.path.abspath(c), strict=False).model()) for c in cifs]

def _score_model_mmtbx(pdb, cif_objects=(), o=None, cf=None, log=sys.stderr):
    M = _load_mmtbx()
    if M is None:
        raise RuntimeError('backend=mmtbx requires phenix.python (mmtbx/iotbx not importable)')
    ovr = read_overrides(o['overridefile'])
    lines = open(pdb).read().splitlines()
    if any(l.split()[-1:] == ['XP'] for l in lines):
        lines = [l for l in lines if l.split()[-1:] != ['XP']]
    occ_tot = _parse_occupancies_lines(lines, filter_xp=False)

    pdb_inp = M['pdb'].input(lines=M['flex'].split_lines('\n'.join(lines)), source_info=None)
    raw_ph = pdb_inp.construct_hierarchy()

    cs = pdb_inp.crystal_symmetry()
    if cs is None or cs.unit_cell() is None:
        box = M['uctbx'].non_crystallographic_unit_cell_with_the_sites_in_its_center(
            sites_cart=pdb_inp.atoms().extract_xyz(), buffer_layer=10)
        cs = box.crystal_symmetry()

    p = M['gm_master'].extract()
    p.stop_for_unknowns = False
    p.pdb_interpretation.allow_polymer_cross_special_position = True
    p.pdb_interpretation.const_shrink_donor_acceptor = o['csda']
    p.pdb_interpretation.restraints_library.cdl = o['cdl']

    null = M['null_out']()
    with _cdl_interp_context(M, bool(o['cdl_interp'])):
        model = M['model'].manager(
            model_input=pdb_inp, crystal_symmetry=cs, restraint_objects=list(cif_objects),
            stop_for_unknowns=False, log=null)
        model.process(pdb_interpretation_params=p, make_restraints=True)
        cl_ph = model.get_hierarchy().deep_copy()

        model.setup_riding_h_manager(idealize=True)
        grm = model.get_restraints_manager().geometry
        if o['cdl']:
            cdl_proxies = M['cdl_setup'](grm)
            M['cdl_mod'].update_restraints(model.get_hierarchy(), grm, cdl_proxies=cdl_proxies, log=null, verbose=False)
        model.set_sites_cart_from_hierarchy()

    sites_cart = model.get_sites_cart()
    site_labels = [a.id_str() for a in model.get_hierarchy().atoms()]

    full = _extract_restraint_records(grm, sites_cart, site_labels, M['origin_ids'], bool(o['round']))
    ramas, rotas, cbs, clashes, cb_res, cl_res = _extract_validation_lists(M, raw_ph, cl_ph, o, ovr=ovr, log=log)
    res = _score_parsed(full, ramas, rotas, cbs, clashes, lambda: occ_tot, o, ovr=ovr)
    if cf and any(cf.values()):
        _write_counterfit_mmtbx(pdb, model, cl_res, cb_res, ovr, **cf)
    return res

class ModelScorer(object):
    """Build mmtbx.model + restraints once for a topology, then score new coordinate sets in memory (~1.6-3.7 s/call)."""
    def __init__(self, pdb, cifs=(), log=sys.stderr, **opts):
        M = _load_mmtbx()
        if M is None:
            raise RuntimeError('ModelScorer requires phenix.python (mmtbx/iotbx not importable)')
        self._M = M
        self._pdb = pdb
        self.log = log
        self.opts = _normalize_opts(opts)
        self._ovr = read_overrides(self.opts['overridefile'])
        lines = open(pdb).read().splitlines()
        if any(l.split()[-1:] == ['XP'] for l in lines):
            lines = [l for l in lines if l.split()[-1:] != ['XP']]
        self._occ_tot = _parse_occupancies_lines(lines, filter_xp=False)

        pdb_inp = M['pdb'].input(lines=M['flex'].split_lines('\n'.join(lines)), source_info=None)
        self._raw_ph = pdb_inp.construct_hierarchy()
        cs = pdb_inp.crystal_symmetry()
        if cs is None or cs.unit_cell() is None:
            box = M['uctbx'].non_crystallographic_unit_cell_with_the_sites_in_its_center(
                sites_cart=pdb_inp.atoms().extract_xyz(), buffer_layer=10)
            cs = box.crystal_symmetry()
        cif_objects = _parse_cif_objects(M, cifs)

        p = M['gm_master'].extract()
        p.stop_for_unknowns = False
        p.pdb_interpretation.allow_polymer_cross_special_position = True
        p.pdb_interpretation.const_shrink_donor_acceptor = self.opts['csda']
        p.pdb_interpretation.restraints_library.cdl = self.opts['cdl']

        null = M['null_out']()
        with _cdl_interp_context(M, bool(self.opts['cdl_interp'])):
            self.model = M['model'].manager(
                model_input=pdb_inp, crystal_symmetry=cs, restraint_objects=cif_objects,
                stop_for_unknowns=False, log=null)
            self.model.process(pdb_interpretation_params=p, make_restraints=True)
            self._cl_ph = self.model.get_hierarchy().deep_copy()
            self.model.setup_riding_h_manager(idealize=True)
            self.grm = self.model.get_restraints_manager().geometry
            self._cdl_proxies = M['cdl_setup'](self.grm) if self.opts['cdl'] else None
            if self._cdl_proxies is not None:
                M['cdl_mod'].update_restraints(self.model.get_hierarchy(), self.grm, cdl_proxies=self._cdl_proxies,
                                               log=null, verbose=False)
            self.model.set_sites_cart_from_hierarchy()
        self._site_labels = [a.id_str() for a in self.model.get_hierarchy().atoms()]

    def score_sites(self, sites_cart=None, include_clashes=True, fullgeo=None, clashes=None, cbdevlog=None,
                    outprefix=None, **override_opts):
        """Score current model coordinates (or a new sites_cart flex.vec3_double / list of (x,y,z) in model atom order)."""
        M = self._M
        o = _normalize_opts(dict(self.opts, **override_opts)) if override_opts else self.opts
        ovr = read_overrides(o['overridefile']) if o['overridefile'] != self.opts['overridefile'] else self._ovr
        cf = dict(fullgeo=fullgeo, clashes=clashes, cbdevlog=cbdevlog, outprefix=outprefix)
        null = M['null_out']()
        if sites_cart is not None:
            if not hasattr(sites_cart, 'size'):
                sites_cart = M['flex'].vec3_double(sites_cart)
            self.model.set_sites_cart(sites_cart)
            cl_ph = self.model.get_hierarchy().deep_copy() if include_clashes else None
            self.model.idealize_h_riding()
            if self._cdl_proxies is not None:
                with _cdl_interp_context(M, bool(o['cdl_interp'])):
                    M['cdl_mod'].update_restraints(self.model.get_hierarchy(), self.grm, cdl_proxies=self._cdl_proxies,
                                                   log=null, verbose=False)
            self.model.set_sites_cart_from_hierarchy()
            val_ph = self.model.get_hierarchy()
        else:
            val_ph = self._raw_ph
            cl_ph = self._cl_ph
        sc = self.model.get_sites_cart()
        full = _extract_restraint_records(self.grm, sc, self._site_labels, M['origin_ids'], bool(o['round']))
        ramas, rotas, cbs, probe, cb_res, cl_res = _extract_validation_lists(
            M, val_ph, cl_ph, o, ovr=ovr, include_clashes=include_clashes, log=self.log)
        res = _score_parsed(full, ramas, rotas, cbs, probe, lambda: self._occ_tot, o, ovr=ovr)
        if any(cf.values()):
            _write_counterfit_mmtbx(self._pdb, self.model, cl_res, cb_res, ovr, **cf)
        return res

# -------------------------------------------------------------------------------------------------------------- phenix CLI
_CDL_INTERP_CLI = (
    'import sys, runpy, mmtbx.conformation_dependent_library as c; '
    'u = c.update_restraints; '
    'c.update_restraints = lambda *a, **k: u(*a, **dict(k, interpolate=True)); '
    'sys.argv = sys.argv[1:]; '
    'runpy.run_module("mmtbx.command_line.geometry_minimization", run_name="__main__", alter_sys=True)'
)

def _mkdtemp():
    """Prefer /dev/shm/$USER on Linux when TMPDIR is unset (avoids node-local /tmp on octamus1), as molprobify_v2.com does."""
    if not any(os.environ.get(k) for k in ('TMPDIR', 'TEMP', 'TMP')) and os.path.isdir('/dev/shm'):
        user = os.environ.get('USER') or os.environ.get('LOGNAME') or (str(os.getuid()) if hasattr(os, 'getuid') else 'user')
        base = os.path.join('/dev/shm', user)
        try:
            os.makedirs(base, exist_ok=True)
            return tempfile.mkdtemp(prefix='wE_v2_', dir=base)
        except OSError:
            pass
    return tempfile.mkdtemp(prefix='wE_v2_')

def run_phenix(pdb, workdir, cifs=(), csda=0.6, cdl=True, cdl_interp=0, jobs=1, log=sys.stderr):
    """Run the five phenix programs molprobify_v2.com runs; returns the paths score_files needs."""
    pdb = os.path.abspath(pdb)
    cifs = [os.path.abspath(c) for c in cifs]
    env = dict(os.environ, PHENIX_OVERWRITE_ALL='1')
    for stale in ('ramalyze.log', 'rotalyze.log', 'cbetadev.log', 'molprobity.log', 'molprobity_coot.py', 'geom.log', 'geom.geo'):
        p = os.path.join(workdir, stale)
        if os.path.exists(p):
            os.remove(p)
    # discard amber extra points, as the scorer does
    lines = open(pdb).read().splitlines()
    has_xp = any(l.split()[-1:] == ['XP'] for l in lines)
    if has_xp:
        lines = [l for l in lines if l.split()[-1:] != ['XP']]
    this = os.path.join(workdir, 'this.pdb')
    open(this, 'w').write('\n'.join(lines) + '\n')
    if has_xp:
        in_dir = os.path.join(workdir, '_in')
        os.makedirs(in_dir, exist_ok=True)
        in_pdb = os.path.join(in_dir, os.path.basename(pdb))
        shutil.copyfile(this, in_pdb)
    else:
        in_pdb = pdb

    gm_args = [in_pdb] + cifs + [
        'macro_cycles=0', 'stop_for_unknowns=false', 'allow_polymer_cross_special_position=True',
        'const_shrink_donor_acceptor=%s' % csda, 'cdl=%s' % ('True' if cdl else 'False'),
        'output_file_name_prefix=' + os.path.join(workdir, 'geom'),
    ]
    if cdl_interp:
        gm_cmd = ['phenix.python', '-c', _CDL_INTERP_CLI, 'phenix.geometry_minimization'] + gm_args
    else:
        gm_cmd = ['phenix.geometry_minimization'] + gm_args

    tasks = [
        (['phenix.ramalyze', in_pdb], 'phenix.ramalyze', 'ramalyze.log', True, False),
        (['phenix.rotalyze', in_pdb], 'phenix.rotalyze', 'rotalyze.log', True, False),
        (['phenix.cbetadev', in_pdb], 'phenix.cbetadev', 'cbetadev.log', True, False),
        (['phenix.molprobity', 'flip_symmetric_amino_acids=True'] + cifs +
         ['outliers_only=False', 'output.probe_dots=False', 'output.coot=True', 'this.pdb'],
         'phenix.molprobity', 'molprobity.log', False, True),
        (gm_cmd, 'phenix.geometry_minimization', 'geom.log', True, True),
    ]
    resolved = []
    for cmd, prog, out, req, merge_err in tasks:
        exe = shutil.which(cmd[0])
        if exe is None:
            raise RuntimeError('%s not found on PATH' % cmd[0])
        resolved.append(([exe] + cmd[1:], prog, os.path.join(workdir, out), req, merge_err))

    max_jobs = max(1, min(int(jobs), len(resolved)))
    if max_jobs == 1:
        for full_cmd, prog, out_path, req, merge_err in resolved:
            if log:
                log.write(prog + '\n')
            err_fh = None if merge_err else open(out_path + '.err', 'w')
            try:
                with open(out_path, 'w') as f:
                    rc = subprocess.call(full_cmd, stdout=f, stderr=subprocess.STDOUT if merge_err else err_fh,
                                         cwd=workdir, env=env)
            finally:
                if err_fh is not None:
                    err_fh.close()
            if rc != 0:
                if req:
                    raise RuntimeError('%s failed (exit %d): see %s' % (prog, rc, out_path if merge_err else out_path + '.err'))
                if log:
                    log.write('WARNING: %s exited with status %d (see %s)\n' % (prog, rc, out_path))
    else:
        # Launch heaviest tasks first so critical path (geometry_minimization, molprobity) starts immediately
        order = [4, 3, 1, 2, 0]
        queue = [resolved[i] for i in order]
        running = []
        results = []
        try:
            while queue or running:
                while queue and len(running) < max_jobs:
                    full_cmd, prog, out_path, req, merge_err = queue.pop(0)
                    if log:
                        log.write(prog + '\n')
                    fh = open(out_path, 'w')
                    efh = None if merge_err else open(out_path + '.err', 'w')
                    proc = subprocess.Popen(full_cmd, stdout=fh, stderr=subprocess.STDOUT if merge_err else efh,
                                            cwd=workdir, env=env)
                    running.append((proc, fh, efh, prog, out_path, req, merge_err))
                still = []
                for proc, fh, efh, prog, out_path, req, merge_err in running:
                    rc = proc.poll()
                    if rc is None:
                        still.append((proc, fh, efh, prog, out_path, req, merge_err))
                    else:
                        fh.close()
                        if efh is not None:
                            efh.close()
                        results.append((rc, prog, out_path, req, merge_err))
                if len(still) == len(running) and still:
                    proc, fh, efh, prog, out_path, req, merge_err = still.pop(0)
                    rc = proc.wait()
                    fh.close()
                    if efh is not None:
                        efh.close()
                    results.append((rc, prog, out_path, req, merge_err))
                running = still
        finally:
            for proc, fh, efh, *_ in running:
                proc.terminate()
                fh.close()
                if efh is not None:
                    efh.close()
        for rc, prog, out_path, req, merge_err in results:
            if rc != 0:
                if req:
                    raise RuntimeError('%s failed (exit %d): see %s' % (prog, rc, out_path if merge_err else out_path + '.err'))
                if log:
                    log.write('WARNING: %s exited with status %d (see %s)\n' % (prog, rc, out_path))

    coot = os.path.join(workdir, 'molprobity_coot.py')
    if not os.path.exists(coot) and log:
        log.write('WARNING: phenix.molprobity did not produce molprobity_coot.py; scoring with 0 probe clashes\n')
    geo = os.path.join(workdir, 'geom.geo')
    if not os.path.exists(geo):
        raise RuntimeError('phenix.geometry_minimization failed: see %s' % os.path.join(workdir, 'geom.log'))
    return dict(pdb=this, geo=geo, rama=os.path.join(workdir, 'ramalyze.log'), rota=os.path.join(workdir, 'rotalyze.log'),
                cbdev=os.path.join(workdir, 'cbetadev.log'), coot=coot)

def score_model(pdb, cifs=(), workdir=None, keep=False, log=sys.stderr, fullgeo=None, clashes=None, cbdevlog=None,
                outprefix=None, **opts):
    """Score a PDB model (uses in-process mmtbx when available unless backend=cli or workdir/keep is requested)."""
    o = _normalize_opts(opts)
    cf = dict(fullgeo=fullgeo, clashes=clashes, cbdevlog=cbdevlog, outprefix=outprefix)
    use_mmtbx = (o['backend'] == 'mmtbx') or (o['backend'] == 'auto' and not workdir and not keep and _load_mmtbx() is not None)
    if use_mmtbx:
        M = _load_mmtbx()
        if M is None:
            raise RuntimeError('backend=mmtbx requires phenix.python (mmtbx/iotbx not importable)')
        cif_objs = _parse_cif_objects(M, cifs)
        return _score_model_mmtbx(pdb, cif_objs, o, cf=cf, log=log)
    tmp = workdir or _mkdtemp()
    os.makedirs(tmp, exist_ok=True)
    try:
        f = run_phenix(pdb, tmp, cifs, o['csda'], o['cdl'], o['cdl_interp'], o['jobs'], log=log)
        return score_files(**f, **cf, **o)
    finally:
        if not keep and not workdir:
            shutil.rmtree(tmp, ignore_errors=True)

def score_models(pdbs, cifs=(), log=sys.stderr, **opts):
    """Score a sequence of PDB models in one warm process (parses ligand CIFs and imports mmtbx once)."""
    o = _normalize_opts(opts)
    use_mmtbx = (o['backend'] == 'mmtbx') or (o['backend'] == 'auto' and _load_mmtbx() is not None)
    if use_mmtbx:
        M = _load_mmtbx()
        if M is None:
            raise RuntimeError('backend=mmtbx requires phenix.python (mmtbx/iotbx not importable)')
        cif_objs = _parse_cif_objects(M, cifs)
        return [_score_model_mmtbx(p, cif_objs, o, log=log) for p in pdbs]
    return [score_model(p, cifs=cifs, log=log, **o) for p in pdbs]

# -------------------------------------------------------------------------------------------------------------- output
def report(res, out=sys.stdout):
    out.write('%-8s %12s %8s %8s %10s\n' % ('category', 'outliers', 'n>0.01', 'N', 'average'))
    for c, t in sorted(res['terms'].items()):
        if t['outliers'] is None:
            out.write('%-8s %12s %8s %8s %10.4f\n' % (c, '-', '-', '-', t['average']))
        else:
            out.write('%-8s %12.4f %8d %8s %10.4f\n' % (c, t['outliers'], t['n_gt_001'], '%g' % t['N'], t['average']))
    out.write('%-8s %12.4f %8d %8s %10s  (physical: sum of repulsive LJ, no weight; in wE only if nbmode physical/both)\n'
              % ('NB_PHYS', res['nb_phys'][0], res['nb_phys'][1], 'contacts', '-'))
    fmt = '%.6g' if res.get('options', {}).get('round', 1) else '%.10g'
    out.write('weighted energy (wE): %s\n' % (fmt % res['wE']))

def _json_summary(res, pdb=None):
    d = dict(wE=res['wE'], terms=res['terms'], nb_phys=res['nb_phys'], nonbond=res['nonbond'],
             n_clashes=res['n_clashes'], options=res['options'])
    if pdb is not None:
        d['pdb'] = pdb
    return d

def main(argv):
    pdbs, cifs, opts, files = [], [], {}, {}
    for arg in argv[1:]:
        if '=' in arg:
            k, v = arg.split('=', 1)
            if k in ('pdb', 'pdbfile'):
                pdbs.append(v)
            elif k in ('geo', 'rama', 'rota', 'cbdev', 'coot', 'json', 'items', 'workdir', 'keep',
                       'fullgeo', 'clashes', 'cbdevlog', 'outprefix'):
                files[k] = v
            elif k in DEFAULTS or k == 'round_numbers':
                opts[k] = v
            elif k == 'override' or k.startswith('override'):      # the tcsh takes any override* key
                opts['overridefile'] = v
            else:
                sys.exit('unknown option %s' % k)
        elif arg.endswith(('.pdb', '.ent')):
            pdbs.append(arg)
        elif arg.endswith('.cif'):
            cifs.append(arg)
        elif arg in ('-h', '--help'):
            print(__doc__)
            return
        else:
            sys.exit('unknown argument %s' % arg)
    if not pdbs:
        sys.exit(__doc__)
    req_files = ('geo', 'rama', 'rota', 'cbdev')
    try:
        _normalize_opts(opts)
    except ValueError as e:
        sys.exit(str(e))

    if len(pdbs) > 1:
        if any(k in files for k in req_files + ('workdir', 'items', 'fullgeo', 'clashes', 'cbdevlog', 'outprefix')):
            sys.exit('batch mode (multiple PDB files) cannot be combined with single-model file/workdir/items/outprefix options')
        results = score_models(pdbs, cifs, **opts)
        for p, res in zip(pdbs, results):
            sys.stdout.write('=== %s ===\n' % p)
            report(res)
        if 'json' in files:
            with open(files['json'], 'w') as jf:
                json.dump([_json_summary(res, p) for p, res in zip(pdbs, results)], jf, indent=1)
        return

    pdb = pdbs[0]
    cf = {k: files[k] for k in ('fullgeo', 'clashes', 'cbdevlog', 'outprefix') if k in files}
    if any(k in files for k in req_files):
        missing = [k for k in req_files if k not in files]
        if missing:
            sys.exit('score_files mode missing required file option(s): %s' % ', '.join(missing))
        res = score_files(pdb, files['geo'], files['rama'], files['rota'], files['cbdev'], files.get('coot'), **cf, **opts)
    else:
        res = score_model(pdb, cifs, files.get('workdir'), keep=files.get('keep', '').lower() in ('1', 'true', 'yes', 'on'),
                          **cf, **opts)
    report(res)
    if 'items' in files:
        with open(files['items'], 'w') as f:
            for v, c, k, N, P, E in res['items']:
                f.write('%.6g %s %d %s %.6g %s\n' % (v, c, k, '%g' % N, P, '%g' % E))
    if 'json' in files:
        with open(files['json'], 'w') as jf:
            json.dump(_json_summary(res), jf, indent=1)

if __name__ == '__main__':
    main(sys.argv)
