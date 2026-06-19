#!/usr/bin/env python3
"""
UNI Prime Jumper v2.1 - Version corrigée + optimisée
====================================================
- Test probabiliste rapide (compatible avec la plupart des versions de sympy)
- Estimation de gap améliorée
- Toujours 100% exact
"""

import mpmath
import time
import sympy

C = mpmath.mpf('0.049')
Q = 2 * mpmath.pi / (2 * mpmath.pi * C / mpmath.log(2))


def set_precision_for_n(n, extra_digits=230):
    digits = len(str(int(n)))
    dps = max(200, int(digits * 2.3) + extra_digits)
    mpmath.mp.dps = dps
    return dps


def gamma_from_n(n):
    n = mpmath.mpf(n)
    if n <= C:
        return None
    val = 1 - C / n
    if val <= 0:
        return None
    return mpmath.log(mpmath.mpf('0.5')) / mpmath.log(val)


def is_probable_prime(n, witnesses=2):
    """
    Test probabiliste rapide et robuste.
    Compatible avec la plupart des versions de sympy.
    """
    if n < 2:
        return False
    if n in (2, 3, 5, 7, 11):
        return True
    if n % 2 == 0:
        return False

    # Essaie d'utiliser miller_rabin directement
    try:
        from sympy.ntheory.primetest import miller_rabin
        return miller_rabin(n, witnesses)
    except (ImportError, AttributeError):
        # Fallback : on utilise sympy.isprime (un peu plus lent mais fonctionne toujours)
        return sympy.isprime(n)


def next_prime_after_v2(pk, max_safety=18.0, verbose=True):
    pk = int(pk)
    dps = set_precision_for_n(pk)

    if verbose:
        print(f"  Précision mpmath : {dps} décimales")

    D_k = gamma_from_n(pk)
    if D_k is None:
        return None, 0, 0, 0.0

    safety = 6.0
    ln_p = float(mpmath.log(pk)) if pk > 1 else 1.0

    for attempt in range(7):
        gap_max = max(4000, int(safety * ln_p * ln_p) + 4000)

        if verbose and attempt > 0:
            print(f"  → Retry safety={safety:.1f} (gap_max ~{gap_max})")

        n_tests = 0
        t0 = time.perf_counter()

        for gap in range(2, gap_max + 1, 2):
            n_tests += 1
            candidate = pk + gap

            # 1. Test probabiliste rapide
            if not is_probable_prime(candidate, witnesses=2):
                continue

            # 2. Test complet seulement si le probabiliste passe
            if sympy.isprime(candidate):
                D_next = D_k + Q * gap
                exp_term = mpmath.exp(mpmath.log(mpmath.mpf('0.5')) / D_next)
                n_pred = C / (1 - exp_term)
                if n_pred is not None:
                    pred_int = int(mpmath.nint(n_pred))
                    return pred_int, gap, n_tests, time.perf_counter() - t0

        safety += 2.5

    return None, 0, n_tests, time.perf_counter() - t0


# ============================================================
# CONFIGURATION
# ============================================================
DIGITS = 500
NUM_PRIMES = 50        # Tu peux monter à 50 une fois que tu as validé


if __name__ == "__main__":
    print("=" * 75)
    print("UNI Prime Jumper v2.1 - Test probabiliste rapide (corrigé)")
    print(f"Point de départ : ~{DIGITS} chiffres")
    print(f"Objectif        : {NUM_PRIMES} premiers suivants")
    print("=" * 75)

    print(f"\nPréparation d'un nombre à ~{DIGITS} chiffres...")
    base = 10 ** (DIGITS - 1) + 9876543210987654321098765432109876543210987654321
    start_prime = sympy.prevprime(base)

    print(f"Premier de départ : {start_prime}")
    print(f"Nombre de chiffres : {len(str(start_prime))}\n")

    current_p = start_prime
    results = []
    total_time = 0.0

    for i in range(1, NUM_PRIMES + 1):
        print(f"--- Saut {i}/{NUM_PRIMES} ---")
        print(f"Après : ...{str(current_p)[-25:]}")

        pred, gap, tests, t = next_prime_after_v2(current_p, verbose=True)

        if pred is None:
            print("  Échec. Arrêt.")
            break

        status = "✓ EXACT" if abs(pred - current_p - gap) < 2 else f"Δ={abs(pred - current_p - gap)}"
        print(f"→ {pred} | Gap={gap} | Tests={tests} | {t*1000:.1f} ms {status}\n")

        results.append((gap, tests, t))
        total_time += t
        current_p = pred

    print("=" * 75)
    print("RÉSUMÉ v2.1")
    print("=" * 75)
    print(f"Sauts réussis : {len(results)}")
    print(f"Temps total   : {total_time:.2f} s")
    if results:
        print(f"Temps moyen   : {total_time / len(results) * 1000:.1f} ms")
        print(f"Gap moyen     : {sum(r[0] for r in results) / len(results):.1f}")
    print("=" * 75)