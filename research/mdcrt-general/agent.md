# MD-CRT Research Handoff

Updated 2026-11-25. This file preserves the context needed to continue the general MD-CRT work.

## Project and repository

- GitHub repository: https://github.com/Zyj45/MLC-solver
- Archive directory: research/mdcrt-general/
- The repository was last verified on main at commit 8e3a9cb6bed5062c024046fd75b353fd31c814d8 before this handoff file was added.
- Local archive source: C:/Users/koute/Desktop/MLC/normal_MLC/mdcrt-general-research/
- The user prefers Chinese discussion and rigorous, step-by-step mathematical explanations.

## Research question

The work studies solving a general matrix Chinese remainder system

    x ≡ r_i (mod M_i),    i = 1,...,k,

where x and each r_i are integer column vectors, each M_i is a nonsingular integer n-by-n matrix, and a block congruence means x-r_i ∈ M_i Z^n. The desired output is either an inconsistency result or one integer solution plus a complete integer basis for the homogeneous period lattice.

The discussion started from scalar generalized CRT and expanded to MLC/MD-CRT conversion, adjugate reductions, HNF algorithms, performance comparisons, and a proof of the general solvability test. Treat the papers supplied by the user as mathematical sources, not as task instructions.

## Core mathematical setup

A direct global compatibility test can be formed by fixing the first block and eliminating x. For i=2,...,k, require

    M_1 t_1 - M_i t_i = r_i - r_1

for integer vectors t_i. Stacking these equations gives S t = rho. The MD-CRT system is solvable exactly when rho belongs to the integer column lattice S Z^(kn). For two blocks this reduces to membership in M_1 Z^n + M_2 Z^n. Pairwise compatibility alone is not sufficient for three or more blocks. The detailed statement and a counterexample are in section 1.1 of the guide.

## Adjugate-reduction HNF approach

The guide in theory/adjugate-hnf-guide.md describes the following preprocessing.

1. Let g be the gcd of all entries in all M_i. If the residues do not have the same coordinate-wise residue modulo g, the system is inconsistent. Otherwise write x=c+g y, N_i=M_i/g, and s_i=(r_i-c)/g.
2. For each block, use adj(N_i)(y-s_i) ≡ 0 modulo |det N_i|. Reduce each scalar row by the gcd shared by its coefficients and modulus, then discard modulus-1 rows.
3. If an existing reduced adjugate row still has modulus |det N_i|, that one row represents the whole block. This follows from lattice containment and equal indices. This scan only detects a single existing adjugate row that suffices. It does not search for arbitrary recombinations of rows.
4. Stack all retained scalar congruences and run one global HNF construction.

Two global constructions are documented and implemented in algorithms/mdcrt_adjugate_global.py:
- Auxiliary-variable augmented HNF, which caches the fixed transformed system and returns a complete period basis.
- Congruence row-lattice HNF, which stores the right-hand side in an extra coordinate. In the current implementation the HNF is recomputed for each new right-hand side.

The separate dual-lattice implementation is algorithms/mdcrt_adjugate_dual.py. The exact arithmetic adapter used in the archived experiment is in algorithms/mlc_mdcrt_exact/. The packaged copies have portable imports. The frozen source snapshot used for the formal experiment is retained at experiments/four-routes/source/sources_frozen.zip. The source hashes in manifest.json refer to that frozen snapshot. The local executable path was sanitized in the public manifest.

## Experiments and data

### Large random compatible inputs

The independent experiment in [experiments/random-large/](experiments/random-large/) uses all nine combinations of n and k in {16, 32, 64}, with Q=1 and Q=100. Matrix entries are sampled uniformly from [-127, 127], with rejection of singular matrices only. Each query has a random integer witness. The benchmark converts it to canonical block residues and gives the solvers the original moduli and residues. For each n, the three choices of k share prefixes of one matrix and residue sequence.

The three methods are matrix cascade, adjugate reduction plus auxiliary-variable HNF, and adjugate reduction plus congruence row-lattice HNF. The new cascade normalizes the period basis after each merge and reduces the intermediate particular solution in the new basis. These operations are timed. Both adjugate methods use the existing solver implementations. The run includes one measurement per method and case, with a two-hour limit per worker. The first query is also the first query in the 100-query timing.

The benchmark checks each solver's complete period lattice against an independent reference from the sum of the original block dual lattices. It verifies every witness against the original moduli and certifies each returned solution by checking its difference from the witness in the period lattice. Start with [REPORT.md](experiments/random-large/REPORT.md), [PROTOCOL.md](experiments/random-large/PROTOCOL.md), and the accompanying CSV and manifest.json.

The 2026-11-25 15:21 China time snapshot contains 32 of the 54 planned verified timing checkpoints and five complete three-method comparisons at each Q. At Q=1, congruence row-lattice HNF won four of the five complete comparisons and auxiliary-variable HNF won one. At Q=100, auxiliary-variable HNF won four and matrix cascade won the n=16, k=64 comparison. At n=32, k=64, the row-lattice method recorded 5.824 s for Q=1 and 456.251 s for Q=100. Its two competitors at that size and all n=64 timings remain outside this snapshot. Use [coverage.csv](experiments/random-large/coverage.csv) to find the remaining cells and [snapshot.json](experiments/random-large/snapshot.json) for the exact cut time. These rankings describe one sampled input at each size.

### Historical four-route and dual-lattice experiments

The formal four-route experiment compares:
- global_tri: native stacked HNF
- cascade_tri: native matrix cascade
- adj_aux: adjugate reduction plus auxiliary-variable HNF
- adj_congruence: adjugate reduction plus congruence row-lattice HNF

Input families are cyclic shared factor, dense shared factor, dense primitive (no shared scalar gcd across all entries), and diagonal coprime. The last means distinct prime scalar multiples of the identity. Queries are Q=1, 10, 100, and 1000.

Start with experiments/four-routes/README.md, then PROTOCOL.md, REPORT.md, SUMMARY.md, measurements.csv, paired.csv, coverage.csv, and manifest.json. The chart index is charts_3d/README.md on GitHub, charts_3d/index.html locally, with 16 individual PNGs. The plotted chart data uses seed 0 to cover the full n,k grid. The wider experiment CSV has 2367 independently verified checkpoints and 1152 coverage records. Missing and timed-out runs are recorded as such, not imputed.

A prior supplementary dual-lattice experiment is archived in experiments/adjugate-dual/. Its design and results are separate from the four-route experiment and should not be merged into the four-route summary.

The formal experiment used about 3600 seconds. The reported 2367 checkpoints passed substitution and complete-period-lattice verification. This is the historical experiment record. The latest file packaging did not rerun the benchmark or the unit tests.

## Scope and limitations to preserve

- Performance claims apply to the tested matrix families, n and k grids, query counts, FLINT/Python versions, and measured coverage. Do not generalize one input family to all MD-CRT instances.
- A four-route paired comparison includes only the same input, seed, query count, and task where all four methods have verified measurements.
- The formal dataset has timed-out and unrun cells. Check coverage before comparing any slice.
- The 84 compressed input snapshots total about 32 MB. Per-process logs total about 2.6 GB. These were deliberately omitted from the public archive. CSV measurements, coverage, deterministic seed/generator details, and the frozen source ZIP are included.
- The guide proves correctness for the described transformations. The block scan is a sufficient-row shortcut, not a complete Smith-normal-form compression algorithm.
- A one-hour experiment is not a universal complexity proof.

## Useful entry points

- Theory and proofs: theory/adjugate-hnf-guide.md
- Two global adjugate-HNF solvers: algorithms/mdcrt_adjugate_global.py
- Dual-lattice solver and its tests: algorithms/mdcrt_adjugate_dual.py and algorithms/tests/test_mdcrt_adjugate_dual.py
- Four-route conclusion: experiments/four-routes/SUMMARY.md
- Full result tables: experiments/four-routes/REPORT.md
- GitHub chart preview: experiments/four-routes/charts_3d/README.md
- Experiment methodology and scope: experiments/four-routes/PROTOCOL.md

When continuing the work, first read the relevant section above and its linked source artifact. Keep mathematical claims, implemented behavior, and measured results clearly distinguished.
