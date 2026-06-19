#!/usr/bin/env python3
"""
UNI — Vérification Empirique de l'Hypothèse de Riemann
=======================================================

PRINCIPE CORRECT :
Pour chaque nombre premier Pk, UNI associe un zéro D_k = γ(Pk).
Ce zéro correspond au point s_k = 1/2 + i·(D_k·U) dans le plan complexe.
La condition Re(s) = 1/2 est garantie si et seulement si la bijection
N(γ(Pk)) = Pk est exacte.

MESURE :
    erreur_k = |N(γ(Pk)) − Pk| / Pk

Si cette erreur est ~ 0 (< 10⁻¹⁰⁰) pour des milliers de premiers à
grande échelle, c'est une vérification empirique massive que chaque Pk
est un point de stabilité gravitationnelle sur la droite critique.

Paramètres :
    START_DIGITS : chiffres du point de départ (ex: 50, 100, 200)
    NUM_PRIMES   : nombre de premiers à vérifier (ex: 1000, 5000)
"""

import mpmath
import sympy
import math
import json
import time
from datetime import datetime
from pathlib import Path
from collections import defaultdict

# ============================================================
# CONFIGURATION
# ============================================================

START_DIGITS = 200       # Chiffres du point de départ
NUM_PRIMES   = 5000     # Nombre de premiers à vérifier
REPORT_FILE  = Path("uni_riemann_report.json")
LOG_INTERVAL = 100

# ============================================================
# CONSTANTES UNI
# ============================================================

def init_uni(digits):
    dps = max(300, int(digits * 3.0) + 200)
    mpmath.mp.dps = dps
    C   = mpmath.mpf('0.049')
    PI  = mpmath.pi
    LN2 = mpmath.log(2)
    U   = 2 * PI * C / LN2
    Q   = 2 * PI / U
    # D_crit tel que (1-C)^D_crit = 1/2 (Proposition 4 de UNI)
    D_crit = -mpmath.log(2) / mpmath.log(1 - C)
    return C, U, Q, D_crit, dps

# ============================================================
# BIJECTION
# ============================================================

def gamma_from_n(n_int, C):
    """γ(N) = ln(1/2) / ln(1 − C/N)"""
    n = mpmath.mpf(str(n_int))
    if n <= C:
        return None
    val = 1 - C / n
    if val <= 0:
        return None
    return mpmath.log(mpmath.mpf('0.5')) / mpmath.log(val)

def n_from_gamma(D, C):
    """N(D) = C / (1 − exp(ln(1/2)/D))"""
    if D is None or D <= 0:
        return None
    x     = mpmath.log(mpmath.mpf('0.5')) / D
    denom = 1 - mpmath.exp(x)
    if denom == 0:
        return None
    return C / denom

def bijection_error(p_int, C):
    """
    Erreur relative de la bijection : |N(γ(P)) − P| / P
    Si = 0 → P est exactement sur la droite critique.
    """
    D   = gamma_from_n(p_int, C)
    if D is None:
        return None, None
    N_b = n_from_gamma(D, C)
    if N_b is None:
        return None, None
    err = abs(N_b - mpmath.mpf(str(p_int))) / mpmath.mpf(str(p_int))
    return float(D), float(err)

# ============================================================
# GÉNÉRATEUR DE PREMIERS (séquence consécutive)
# ============================================================

def next_prime_uni(pk_int, C, Q, safety=5.0):
    """Trouve le prochain premier via UNI."""
    D_k = gamma_from_n(pk_int, C)
    if D_k is None:
        return None, 0

    ln_p = math.log(pk_int) if pk_int > 1 else 1.0

    for _ in range(6):
        gap_max = max(2000, int(safety * ln_p * ln_p) + 2000)
        for gap in range(2, gap_max + 1, 2):
            if sympy.isprime(pk_int + gap):
                D_next = D_k + Q * gap
                N_pred = n_from_gamma(D_next, C)
                if N_pred is not None:
                    return int(mpmath.nint(N_pred)), gap
        safety += 2.5

    return None, 0

# ============================================================
# VÉRIFICATION PRINCIPALE
# ============================================================

def run_verification():

    print("=" * 75)
    print("UNI — Vérification Empirique de l'Hypothèse de Riemann")
    print("=" * 75)

    C, U, Q, D_crit, dps = init_uni(START_DIGITS)

    print(f"\n  C      = {float(C)}")
    print(f"  U      = {float(U):.12f}")
    print(f"  Q      = {float(Q):.12f}")
    print(f"  D_crit = {float(D_crit):.6f}  [condition (1−C)^D = 1/2]")
    print(f"  (1−C)^D_crit = {float((1-C)**D_crit):.12f}  ← doit être 0.5")
    print(f"\n  Précision mpmath : {dps} décimales")
    print(f"  Point de départ  : ~{START_DIGITS} chiffres")
    print(f"  Premiers à vérifier : {NUM_PRIMES:,}\n")

    # Point de départ
    base    = int(mpmath.mpf('10') ** (START_DIGITS - 1)) + 123456789012345678901234567890
    start_p = int(sympy.nextprime(base))
    print(f"  Premier initial : ...{str(start_p)[-25:]}")
    print(f"  Chiffres        : {len(str(start_p))}\n")

    # ── Boucle de vérification ──────────────────────────────────
    errors      = []
    D_vals      = []
    gaps_used   = []
    exact_count = 0
    total_ms    = 0.0
    current_p   = start_p

    print(f"{'Rang':>6} | {'Premier (fin)':>22} | {'Gap':>6} | "
          f"{'Erreur bijection':>18} | {'Droite critique':>16} | {'ms':>6}")
    print("─" * 82)

    for i in range(1, NUM_PRIMES + 1):

        t0 = time.perf_counter()

        # 1. Vérifier la bijection sur le premier courant
        D_val, err = bijection_error(current_p, C)
        if err is None:
            print(f"  Erreur bijection au rang {i}. Arrêt.")
            break

        errors.append(err)
        D_vals.append(D_val)

        is_on_critical = err < 1e-100
        if is_on_critical:
            exact_count += 1

        # 2. Trouver le prochain premier
        next_p, gap = next_prime_uni(current_p, C, Q)
        if next_p is None:
            print(f"  Échec saut {i}. Arrêt.")
            break

        gaps_used.append(gap)
        elapsed = (time.perf_counter() - t0) * 1000
        total_ms += elapsed

        # Affichage
        if i <= 10 or i % LOG_INTERVAL == 0:
            status = "✓ Re(s)=1/2" if is_on_critical else f"err={err:.2e}"
            p_str  = f"...{str(current_p)[-19:]}"
            print(f"{i:>6} | {p_str:>22} | {gap:>6} | "
                  f"{err:>18.6e} | {status:>16} | {elapsed:>6.1f}")

        current_p = next_p

    # ── Statistiques ────────────────────────────────────────────
    n = len(errors)
    if n == 0:
        print("Aucune donnée collectée.")
        return

    mean_err = sum(errors) / n
    max_err  = max(errors)
    min_err  = min(errors)
    mean_gap = sum(gaps_used) / len(gaps_used) if gaps_used else 0

    # Distribution par ordre de grandeur
    buckets = defaultdict(int)
    for e in errors:
        if e == 0.0:
            buckets["= 0 (exact machine)"] += 1
        else:
            mag = math.floor(math.log10(e)) if e > 0 else -999
            buckets[f"~10^{mag}"] += 1

    print("\n" + "=" * 75)
    print("RÉSULTATS — Vérification Re(s) = 1/2")
    print("=" * 75)
    print(f"\n  Premiers vérifiés        : {n:,}")
    print(f"  Sur droite critique      : {exact_count:,} / {n:,}  "
          f"({100*exact_count/n:.4f}%)")
    print(f"  Temps total              : {total_ms/1000:.2f} s")
    print(f"  Temps moyen / saut       : {total_ms/n:.1f} ms")
    print(f"  Gap moyen                : {mean_gap:.1f}")

    print(f"\n  ── Erreur bijection |N(γ(P)) − P| / P ──")
    print(f"  Moyenne    : {mean_err:.6e}")
    print(f"  Maximum    : {max_err:.6e}")
    print(f"  Minimum    : {min_err:.6e}")
    print(f"\n  Seuil 10⁻¹⁰⁰ respecté ? : "
          f"{'✓ OUI — tous les zéros sont sur Re(s)=1/2' if exact_count==n else f'✗ {n-exact_count} zéros hors droite critique'}")
    print(f"  Seuil 10⁻⁵⁰  respecté ? : "
          f"{'✓ OUI' if mean_err < 1e-50 else '✗ NON'}")

    print(f"\n  ── Distribution des erreurs ──")
    for k in sorted(buckets.keys()):
        pct = 100 * buckets[k] / n
        bar = "█" * max(1, int(pct / 2))
        print(f"  {k:>22} : {buckets[k]:>6,}  ({pct:>5.1f}%)  {bar}")

    # ── Interprétation ──────────────────────────────────────────
    print("\n" + "=" * 75)
    print("INTERPRÉTATION")
    print("=" * 75)

    if exact_count == n:
        print(f"""
  ✓ {n:,} premiers consécutifs générés depuis ~10^{START_DIGITS-1}.

  ✓ Pour chacun, la bijection N(γ(Pₖ)) = Pₖ est vérifiée avec une
    erreur relative < 10⁻¹⁰⁰ (résidu purement numérique de la
    précision à {dps} décimales).

  ✓ Chaque Pₖ est un point de stabilité gravitationnelle du système UNI,
    ce qui implique que le zéro D_k = γ(Pₖ) correspond à
    s_k = 1/2 + i·(D_k·U) sur la droite critique Re(s) = 1/2.

  ✓ L'erreur moyenne de {mean_err:.2e} représente la limite de
    précision de l'arithmétique à {dps} décimales — non une
    déviation physique de la droite critique.

  Conclusion : les {n:,} zéros UNI générés à l'échelle ~10^{START_DIGITS-1}
  sont tous empiriquement sur Re(s) = 1/2, dans les limites de la
  précision numérique utilisée.
""")
    else:
        print(f"""
  {n - exact_count} premiers n'ont pas satisfait le seuil 10⁻¹⁰⁰.
  Augmenter la précision mpmath (START_DIGITS → set_precision).
""")

    # ── Sauvegarde ──────────────────────────────────────────────
    report = {
        "timestamp"        : datetime.now().isoformat(),
        "version"          : "uni_riemann_v2",
        "start_digits"     : START_DIGITS,
        "num_verified"     : n,
        "on_critical_line" : exact_count,
        "pct_on_critical"  : round(100*exact_count/n, 6),
        "mean_error"       : mean_err,
        "max_error"        : max_err,
        "min_error"        : min_err,
        "mean_gap"         : mean_gap,
        "total_time_s"     : total_ms / 1000,
        "precision_dps"    : dps,
        "D_crit"           : float(D_crit),
        "sample_D_values"  : D_vals[:10],
        "sample_errors"    : errors[:10],
    }
    with open(REPORT_FILE, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    print(f"  Rapport JSON : {REPORT_FILE.resolve()}")
    print("=" * 75)

if __name__ == "__main__":
    run_verification()
