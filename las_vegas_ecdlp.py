#!/usr/bin/env python3
import sys
import argparse
import random
import time
import multiprocessing as mp
from itertools import combinations
from math import comb
from queue import Empty
from ecdsa import ellipticcurve, curves

# Parameters for the secp256k1 curve
def decompress_pubkey(pubkey_hex):
    if pubkey_hex[:2] not in ('02', '03') or len(pubkey_hex) != 66:
        raise ValueError(f"Invalid format (hex len={len(pubkey_hex)}), expected 66 characters")
    prefix = int(pubkey_hex[:2], 16)
    x = int(pubkey_hex[2:], 16)
    p = curves.SECP256k1.curve.p()
    alpha = (pow(x, 3, p) + 7) % p
    beta = pow(alpha, (p + 1) // 4, p)
    y = beta if (beta % 2) == (prefix & 1) else p - beta
    curve = curves.SECP256k1.curve
    if not curve.contains_point(x, y):
        raise ValueError("Point is not on the secp256k1 curve")
    return x, y

# Generate monomials x^i * y^j for i + j <= n_prime
def get_monomials(n_prime):
    monos = []
    for i in range(n_prime + 1):
        for j in range(n_prime + 1 - i):
            monos.append((i, j))
    return monos

# Evaluate the monomials at a given point P
def eval_monomials(P, monomials):
    x, y = P.x(), P.y()
    p = curves.SECP256k1.curve.p()
    return [(pow(x, i, p) * pow(y, j, p)) % p for i, j in monomials]

# Compute the left-kernel (nullspace) in GF(p)
def nullspace_mod_p(M, mod):
    nrows = len(M)
    if nrows == 0:
        return []
    ncols = len(M[0])
    # Transpose to work by columns
    A = [[M[r][c] for r in range(nrows)] for c in range(ncols)]
    pivots = []
    row = 0
    for col in range(nrows):
        for r in range(row, ncols):
            if A[r][col] % mod != 0:
                A[row], A[r] = A[r], A[row]
                inv = pow(A[row][col], -1, mod)
                for c in range(col, nrows):
                    A[row][c] = (A[row][c] * inv) % mod
                for rr in range(ncols):
                    if rr != row and A[rr][col] % mod != 0:
                        factor = A[rr][col]
                        for c in range(col, nrows):
                            A[rr][c] = (A[rr][c] - factor * A[row][c]) % mod
                pivots.append(col)
                row += 1
                break
        if row == ncols:
            break
    free_cols = [c for c in range(nrows) if c not in pivots]
    basis = []
    for free in free_cols:
        vec = [0] * nrows
        vec[free] = 1
        for i, piv in enumerate(pivots):
            s = sum(A[i][j] * vec[j] for j in range(piv + 1, nrows)) % mod
            vec[piv] = (-s) % mod
        basis.append(vec)
    return basis

# Search for a sparse vector with exactly l zeros
def find_sparse_vector(basis, l, mod):
    dim = len(basis)
    if dim == 0:
        return None
    length = len(basis[0])
    # Try all combinations of basis vectors up to full dimension
    for r in range(1, dim + 1):
        for combo in combinations(range(dim), r):
            vec = [0] * length
            for idx in combo:
                for i in range(length):
                    vec[i] = (vec[i] + basis[idx][i]) % mod
            if sum(1 for x in vec if x == 0) == l:
                return vec
    return None

# Parallel worker with progress notifications
def worker(n_prime, monomials, Qx, Qy, start_time, counter, total, progress_interval, queue):
    curve = curves.SECP256k1.curve
    G = curves.SECP256k1.generator
    n = curves.SECP256k1.order
    p = curve.p()
    k = l = 3 * n_prime
    Q = ellipticcurve.Point(curve, Qx, Qy)
    iteration = 0
    while True:
        iteration += 1
        # Build matrix M for this attempt
        M, I_list, J_list = [], [], []
        for _ in range(k - 1):
            r = random.randrange(1, n)
            I_list.append(r)
            M.append(eval_monomials(r * G, monomials))
        for _ in range(l + 1):
            r = random.randrange(1, n)
            J_list.append(r)
            M.append(eval_monomials(-(r * Q), monomials))
        # Update counter and possibly report progress
        with counter.get_lock():
            counter.value += 1
            if counter.value % progress_interval == 0:
                queue.put(('progress', counter.value))
        # Solve nullspace and look for sparse vector
        basis = nullspace_mod_p(M, p)
        vec = find_sparse_vector(basis, l, p)
        if vec:
            A = sum(I_list[i] for i in range(k - 1) if vec[i] != 0) % n
            B = sum(J_list[j] for j in range(l + 1) if vec[(k - 1) + j] != 0) % n
            m = (A * pow(B, -1, n)) % n
            elapsed = time.time() - start_time
            queue.put(('result', counter.value, elapsed, m))
            return

# Main function
def main():
    parser = argparse.ArgumentParser(description="Solve ECDLP with multiprocessing and progress reporting")
    parser.add_argument("pubkey", help="Compressed public key in hex (33 bytes)")
    parser.add_argument("--nprime", type=int, default=2, help="Parameter n' for monomial bounds")
    parser.add_argument("--progress-interval", type=int, default=10, help="Iteration interval to report progress")
    args = parser.parse_args()

    print("[*] Decompressing public key...")
    try:
        Qx, Qy = decompress_pubkey(args.pubkey)
        print(f"[+] Q.x = {hex(Qx)}")
        print(f"[+] Q.y = {hex(Qy)}\n")
    except ValueError as e:
        print(f"[!] Pubkey error: {e}")
        sys.exit(1)

    n_prime = args.nprime
    monomials = get_monomials(n_prime)
    cols = len(monomials)
    k = l = 3 * n_prime
    total = comb(k + l, l)
    chance = total / curves.SECP256k1.order
    print(f"Parameters: Nprime={n_prime}. Success likelihood ~1 in {int(1/chance):,}\n")

    manager = mp.Manager()
    counter = mp.Value('L', 0)
    queue = manager.Queue()
    workers = []
    # Spawn one worker per CPU core
    for _ in range(mp.cpu_count()):
        p = mp.Process(
            target=worker,
            args=(n_prime, monomials, Qx, Qy,
                  time.time(), counter, total, args.progress_interval, queue)
        )
        p.start()
        workers.append(p)

    try:
        while True:
            tag, *data = queue.get()
            if tag == 'progress':
                iters = data[0]
                print(f"[...] Progress: {iters} iterations", end='\r', flush=True)
            else:
                iters, elapsed, m = data  # unpack result data
                # Corrected f-string: ensure quotes and content are on the same line
                print(f"\n[✔] Found after {iters} iterations in {elapsed:.1f}s: 0x{m:x}")
                break
    finally:
        for p in workers:
            p.terminate()

if __name__ == '__main__':
    main()