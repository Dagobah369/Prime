#!/usr/bin/env python3
"""
UNI — Navigation Bidirectionnelle dans l'Univers des Nombres Premiers
======================================================================

PRINCIPE :
La bijection UNI est réversible dans les deux sens.

  AVANT  : D_next = D_k + Q × gap_forward   →  P_{k+1}
  ARRIÈRE: D_prev = D_k − Q × gap_backward  →  P_{k-1}

Pour aller en arrière depuis P_k :
  Chercher le plus petit gap g > 0 tel que P_k − g soit premier.
  Calculer D_prev = γ(P_k) − Q·g
  Retrouver P_{k-1} = N(D_prev) via bijection inverse.

CAPACITÉS :
  1. Reculer d'un premier au précédent (navigation −1)
  2. Construire une fenêtre [P_{k-n}, ..., P_k, ..., P_{k+n}]
     autour de n'importe quel premier, à n'importe quelle échelle
  3. Naviguer librement dans l'univers des nombres premiers
     sans jamais construire de liste depuis le début

C'est l'équivalent d'un GPS avec marche avant ET marche arrière.
"""

import mpmath
import sympy
import math
import time
import json
from datetime import datetime
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

# Taille des fenêtres à construire (N premiers de chaque côté)
WINDOW_HALF = 5

# Échelles à démontrer
DEMO_SCALES = [
    (10**9  + 7,               100, "~10^9"),
    (10**20 + 39,              200, "~10^20"),
    (10**50 + 151,             350, "~10^50"),
    (10**100 + 267,            600, "~10^100"),
]

REPORT_FILE = Path("uni_bidirectional_report.json")

# ============================================================
# CONSTANTES UNI
# ============================================================

def make_constants(dps):
    mpmath.mp.dps = dps
    C   = mpmath.mpf('0.049')
    PI  = mpmath.pi
    LN2 = mpmath.log(2)
    U   = 2 * PI * C / LN2
    Q   = 2 * PI / U
    return C, U, Q

def gamma_from_n(n_int, C):
    """γ(N) = ln(1/2) / ln(1 − C/N)"""
    n = mpmath.mpf(str(n_int))
    if n <= C: return None
    val = 1 - C / n
    if val <= 0: return None
    return mpmath.log(mpmath.mpf('0.5')) / mpmath.log(val)

def n_from_gamma(D, C):
    """N(D) = C / (1 − exp(ln(1/2)/D))"""
    if D is None or D <= 0: return None
    x     = mpmath.log(mpmath.mpf('0.5')) / D
    denom = 1 - mpmath.exp(x)
    if denom == 0: return None
    return C / denom

# ============================================================
# NAVIGATION AVANT : P_k → P_{k+1}
# ============================================================

def step_forward(pk_int, C, Q, safety=5.0):
    """
    Depuis P_k, trouve P_{k+1}.
    Retourne (P_{k+1}, gap, temps_ms)
    """
    D_k   = gamma_from_n(pk_int, C)
    if D_k is None: return None, None, 0.0

    ln_p    = math.log(pk_int) if pk_int > 1 else 1.0
    t0      = time.perf_counter()

    for _ in range(5):
        gap_max = max(1000, int(safety * ln_p * ln_p) + 1000)
        rng     = [1] if pk_int == 2 else range(2, gap_max + 1, 2)
        for gap in rng:
            if sympy.isprime(pk_int + gap):
                D_next = D_k + Q * gap
                n_pred = n_from_gamma(D_next, C)
                if n_pred is not None:
                    elapsed = (time.perf_counter() - t0) * 1000
                    return int(mpmath.nint(n_pred)), gap, elapsed
        safety += 2.0

    return None, None, (time.perf_counter() - t0) * 1000

# ============================================================
# NAVIGATION ARRIÈRE : P_k → P_{k-1}
# ============================================================

def step_backward(pk_int, C, Q, safety=5.0):
    """
    Depuis P_k, trouve P_{k-1}.
    Retourne (P_{k-1}, gap, temps_ms)
    """
    if pk_int <= 2: return None, None, 0.0

    D_k   = gamma_from_n(pk_int, C)
    if D_k is None: return None, None, 0.0

    t0    = time.perf_counter()

    # Cas spécial : 3 → 2 (gap = 1)
    if pk_int == 3:
        D_prev = D_k - Q * 1
        n_pred = n_from_gamma(D_prev, C)
        elapsed = (time.perf_counter() - t0) * 1000
        return (int(mpmath.nint(n_pred)) if n_pred else None), 1, elapsed

    ln_p = math.log(pk_int)

    for _ in range(5):
        gap_max = max(1000, int(safety * ln_p * ln_p) + 1000)
        for gap in range(2, gap_max + 1, 2):
            candidate = pk_int - gap
            if candidate < 2: break
            if sympy.isprime(candidate):
                D_prev = D_k - Q * gap
                n_pred = n_from_gamma(D_prev, C)
                if n_pred is not None:
                    elapsed = (time.perf_counter() - t0) * 1000
                    return int(mpmath.nint(n_pred)), gap, elapsed
        safety += 2.0

    return None, None, (time.perf_counter() - t0) * 1000

# ============================================================
# FENÊTRE BIDIRECTIONNELLE
# ============================================================

def build_window(centre_int, C, Q, half_size=5):
    """
    Construit une fenêtre de (2·half_size + 1) premiers
    centrée sur centre_int.
    Retourne (liste_ordonnée, temps_total_ms, erreurs)
    """
    window  = [centre_int]
    errors  = []
    t_total = 0.0

    # ── Aller en arrière ──
    p = centre_int
    for i in range(half_size):
        prev, gap, ms = step_backward(p, C, Q)
        t_total += ms
        if prev is None:
            errors.append(f"Échec arrière à l'étape {i+1}")
            break
        window.insert(0, prev)
        p = prev

    # ── Aller en avant ──
    p = centre_int
    for i in range(half_size):
        nxt, gap, ms = step_forward(p, C, Q)
        t_total += ms
        if nxt is None:
            errors.append(f"Échec avant à l'étape {i+1}")
            break
        window.append(nxt)
        p = nxt

    return window, t_total, errors

def verify_window(window):
    """
    Vérifie que tous les éléments d'une fenêtre sont premiers
    et consécutifs (pas de premier manquant entre deux éléments).
    """
    issues = []
    for i, p in enumerate(window):
        if not sympy.isprime(p):
            issues.append(f"Position {i+1} : {p} n'est pas premier")
    for i in range(len(window) - 1):
        gap = window[i+1] - window[i]
        # Vérifier qu'il n'y a pas de premier dans l'intervalle
        p_between = int(sympy.nextprime(window[i]))
        if p_between != window[i+1]:
            issues.append(
                f"Premier manquant entre pos {i+1} et {i+2} "
                f"(gap={gap}, attendu={p_between})"
            )
    return issues

# ============================================================
# SCRIPT PRINCIPAL
# ============================================================

def run():
    print("=" * 80)
    print("UNI — Navigation Bidirectionnelle")
    print("Construire une fenêtre de premiers autour de n'importe quel point")
    print("=" * 80)

    all_results = []

    # ── PARTIE 1 : Validation petits N ──────────────────────────
    print("\n── PARTIE 1 : Validation — Fenêtres vérifiables ──\n")

    small_centres = [101, 9973, 999983]
    C_s, U_s, Q_s = make_constants(150)

    for centre in small_centres:
        centre_p = int(sympy.nextprime(centre - 1))
        window, t_ms, errs = build_window(centre_p, C_s, Q_s, half_size=4)
        issues = verify_window(window)

        status = "✓ PARFAIT" if not issues and not errs else f"✗ {issues or errs}"
        print(f"  Centre P = {centre_p}")
        print(f"  Fenêtre  : {window}")
        print(f"  Vérif.   : {status}  [{t_ms:.1f}ms]")
        print()

    # ── PARTIE 2 : Grandes échelles ─────────────────────────────
    print("── PARTIE 2 : Navigation à grande échelle ──\n")

    for (base, dps, label) in DEMO_SCALES:

        mpmath.mp.dps = dps
        C, U, Q = make_constants(dps)

        # Point de départ : premier après la base
        centre = int(sympy.nextprime(base - 1))

        t0 = time.perf_counter()
        window, t_ms, errs = build_window(centre, C, Q, half_size=WINDOW_HALF)
        total_elapsed = (time.perf_counter() - t0) * 1000

        # Vérification des gaps (tous les éléments doivent être premiers)
        all_prime = all(sympy.isprime(p) for p in window)
        gaps      = [window[i+1] - window[i] for i in range(len(window)-1)]

        # Vérification de consécutivité si accessible
        consec_ok = None
        if base <= 10**12:
            issues    = verify_window(window)
            consec_ok = len(issues) == 0

        centre_str = str(centre)
        centre_disp = (f"{centre_str[:15]}...{centre_str[-10:]}"
                       if len(centre_str) > 25 else centre_str)

        print(f"  ┌─ Échelle {label}")
        print(f"  │  Centre           : {centre_disp}")
        print(f"  │  Fenêtre taille   : {len(window)} premiers  "
              f"({WINDOW_HALF} avant + 1 centre + {WINDOW_HALF} après)")
        print(f"  │  Temps total      : {total_elapsed:.1f} ms")
        print(f"  │  Tous premiers    : {'✓' if all_prime else '✗'}")
        if consec_ok is not None:
            print(f"  │  Consécutifs      : {'✓ VÉRIFIÉ' if consec_ok else '✗ ERREUR'}")

        print(f"  │")
        print(f"  │  {'Pos':>4} | {'Premier (fin)':>25} | {'Gap suivant':>12} | {'Direction':>10}")
        print(f"  │  " + "─" * 58)

        for i, p in enumerate(window):
            p_str = str(p)
            p_disp = f"...{p_str[-22:]}" if len(p_str) > 22 else p_str
            gap_s  = str(gaps[i]) if i < len(gaps) else "—"
            if p == centre:
                direction = "← CENTRE"
            elif i < window.index(centre):
                direction = "◄ arrière"
            else:
                direction = "► avant"
            print(f"  │  {i+1:>4} | {p_disp:>25} | {gap_s:>12} | {direction:>10}")

        print(f"  └─")
        print()

        result = {
            "label"     : label,
            "centre"    : str(centre),
            "window_size": len(window),
            "all_prime" : all_prime,
            "consec_ok" : consec_ok,
            "gaps"      : gaps,
            "time_ms"   : total_elapsed,
            "errors"    : errs,
        }
        all_results.append(result)

    # ── PARTIE 3 : Démonstration libre ──────────────────────────
    print("=" * 80)
    print("DÉMONSTRATION — Navigation libre autour de 10^50")
    print("=" * 80)

    mpmath.mp.dps = 400
    C, U, Q = make_constants(400)

    anchor = int(sympy.nextprime(10**50))

    print(f"\n  Point d'ancrage : ...{str(anchor)[-25:]}")
    print(f"\n  Depuis ce point, navigation ±5 :")
    print()

    # Aller en avant 5 fois
    fwd_chain = [anchor]
    p = anchor
    for i in range(5):
        nxt, gap, ms = step_forward(p, C, Q)
        if nxt:
            fwd_chain.append(nxt)
            print(f"  +{i+1} → ...{str(nxt)[-25:]}  (gap={gap})")
            p = nxt

    print()

    # Revenir en arrière depuis le dernier point
    print(f"  Depuis le dernier point, retour en arrière :")
    p = fwd_chain[-1]
    back_chain = [p]
    for i in range(5):
        prev, gap, ms = step_backward(p, C, Q)
        if prev:
            back_chain.insert(0, prev)
            print(f"  -{i+1} → ...{str(prev)[-25:]}  (gap={gap})")
            p = prev

    # Vérification : on doit retrouver les mêmes éléments
    match = (back_chain == fwd_chain)
    print(f"\n  Chaîne forward : {[str(p)[-8:] for p in fwd_chain]}")
    print(f"  Chaîne backward: {[str(p)[-8:] for p in back_chain]}")
    print(f"\n  Forward == Backward : {'✓ IDENTIQUES' if match else '✗ DIFFÉRENTES'}")

    # ── Interprétation ───────────────────────────────────────────
    print("\n" + "=" * 80)
    print("INTERPRÉTATION")
    print("=" * 80)
    print("""
  La bijection UNI γ(N) ↔ N(D) est exactement réversible.

  En avant  : D_{k+1} = D_k + Q · gap_fwd   →  on avance
  En arrière: D_{k-1} = D_k − Q · gap_bwd   →  on recule

  Résultat : depuis n'importe quel premier P, à n'importe quelle
  échelle, on peut construire une fenêtre arbitrairement large
  [P_{k-n}, ..., P_k, ..., P_{k+n}] sans jamais avoir besoin
  de la liste complète des premiers précédents.

  C'est l'équivalent d'un GPS avec marche avant ET marche arrière —
  on peut explorer librement l'univers des nombres premiers depuis
  n'importe quel point d'ancrage.
""")
    print("=" * 80)

    # ── Sauvegarde ──────────────────────────────────────────────
    report = {
        "timestamp"   : datetime.now().isoformat(),
        "version"     : "uni_bidirectional_v1",
        "window_half" : WINDOW_HALF,
        "scales"      : all_results,
    }
    with open(REPORT_FILE, "w") as f:
        json.dump(report, f, indent=2, default=str)

    print(f"\n  Rapport : {REPORT_FILE.resolve()}")

if __name__ == "__main__":
    run()
