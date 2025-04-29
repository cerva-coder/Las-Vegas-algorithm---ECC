# Las Vegas ECDLP Solver (secp256k1)

This repository provides a Python script to solve the Elliptic Curve Discrete Logarithm Problem (ECDLP) on secp256k1 using a randomized Las Vegas algorithm (https://arxiv.org/pdf/1703.07544) and multiprocessing for performance.

## Overview

- **Algorithm**: Implements a Las Vegas approach, generating random combinations of monomials evaluated on points to build a matrix whose nullspace reveals the discrete logarithm.
- **Parallelization**: Spawns one worker process per CPU core to test independent random trials concurrently.
- **Progress Reporting**: Periodically outputs the number of iterations completed.

## How It Works

1. **Decompress Public Key**: Converts a compressed hex public key (33 bytes) to an (x, y) point on secp256k1.
2. **Monomial Generation**: Builds all monomials \(x^i y^j\) such that \(i + j \le n'\).
3. **Matrix Assembly**: For each trial, constructs a matrix \(M\) by evaluating monomials on random multiples of the generator \(G\) and the target point \(Q\).
4. **Nullspace in GF(p)**: Solves for the left nullspace (kernel) of \(M\) over the prime field.
5. **Sparse Vector Search**: Searches for a nullspace vector with exactly \(l\) zeros; leads to two scalar sums whose ratio gives the discrete log.
6. **Las Vegas Guarantee**: Always correct on success; randomized trials continue until a solution (private key) is found.

## Parameters

- `--nprime` (integer): Denoted \(n'\), bounds the total degree of monomials. Larger values increase the dimension of the search space and probability of success per trial, at the cost of heavier linear algebra.
- `l` and `k`: Internally set to \(3 n'\). The matrix dimensions and sparsity conditions depend on this choice.

## Usage

```bash
# Basic usage
python las_vegas_ecdlp.py <compressed_pubkey_hex> --nprime 2

# Example
python las_vegas_ecdlp.py 02e0a8b039282faf6fe0fd769cfbc4b6b4cf8758ba68220eac420e32b91ddfa673 --nprime 5
