# Quantifying-entanglement-structures
This repository implements polynomial-level observables for detecting entanglement structure in multipartite quantum systems. 
# Nine-Qubit Cluster State: Entanglement Depth and Separability Analysis

Wolfram Mathematica code for evaluating entanglement-depth and \(k\)-separability inequalities for a nine-qubit cluster-type state using measured Pauli correlators and computational-basis populations.

## Overview

The target nine-qubit state is

\[
\begin{aligned}
|CL_9\rangle=\frac{1}{2\sqrt{2}}\big(&|000000000\rangle
+|000000111\rangle+|000111011\rangle+|000111100\rangle\\
&+|111000011\rangle+|111000100\rangle
+|111111000\rangle+|111111111\rangle\big).
\end{aligned}
\]

The set of four computational-basis pairs used by the witness is

\[
\Omega=\{(1,8),(2,7),(3,6),(4,5)\},
\]

where the indices refer to the ordered basis strings in the state above. Each pair consists of bitwise-complementary strings.

The main notebook evaluates:

- **Entanglement-depth bounds:** Eq. S36 for \(d=1,2,3,4\), and Eq. S11 for \(d=5,\ldots,9\).
- **\(k\)-separability bounds:** the corresponding \(k\)-partition construction for \(k=3,\ldots,8\).
- **Experimental four-\(Y\) correlators:** signed sums of the \(\binom{9}{4}=126\) Pauli strings containing four \(Y\) operators and five \(X\) operators.
- **Numerical comparisons:** `LEFT`, `RIGHT`, `Delta =2 LEFT - 2RIGHT`, and `ViolationQ` for each tested bound.

The Mathematica implementations use the weighted swap-set constructions denoted by \(\Gamma\) and \(\Upsilon\), with addition and maximum unions as specified in the associated supplementary definitions (S3–S6). Consult the original paper for formal definitions and the precise interpretation of each inequality.

## Files

| File | Description |
|---|---|
| `CL9_9Q_4Y_S36_S11_k3_to_k8.nb` | Main Mathematica notebook covering depth bounds, separability bounds, experimental imports, and an ideal-state example. |

## Requirements

- Wolfram Mathematica / Wolfram Language with support for `Association`, `Dataset`, and `Import`.
- Local measurement files for experimental evaluation. They are **not included** in this repository unless separately provided.
- Sufficient memory and time for enumerating nine-qubit partitions and weighted exchange sets.

The notebook has not been benchmarked across Mathematica versions; please record the version you use in your analysis.

## Measurement conventions

The implementation uses **1-based Mathematica row indices** and the local Pauli encoding

- `1` = \(X\)
- `2` = \(Y\)
- `3` = \(Z\)

The first measurement row corresponds to the all-\(X\) setting; row **340** corresponds to the all-\(Z\) setting. The four-\(Y\) measurement rows are explicitly stored in the notebook, along with rows for the other even numbers of \(Y\) operators.

Each repetition is expected to contain 341 files named

```text
basis_measurement0.txt
basis_measurement1.txt
...
basis_measurement340.txt
```

Each file is expected to supply 512 measurement outcomes in the first data column. Data are normalized per measurement setting before calculating the parity correlator.

**Indexing is important:** Mathematica row `r` corresponds to `basis_measurement(r - 1).txt`. The output-basis ordering used in the analysis assigns qubit 1 to the least significant bit:

\[
x=1+\sum_{q=1}^{9}(i_q-1)2^{q-1},\qquad i_q\in\{1,2\}.
\]

Check this convention against the acquisition/export software before evaluating experimental data.

## Getting started

1. Download or clone the repository and open `CL9_9Q_4Y_S36_S11_k3_to_k8.nb` in Mathematica.
2. Evaluate the notebook cells in order to define the state, measurement mapping, exchange sets, and comparison functions.
3. Run the included *ideal-state* calculation first. It is a consistency check that does not require local measurement files.
4. For experiment data, edit the `baseDir` and `folders` variables to point to the actual acquisition directories. Ensure different repetitions refer to genuinely different datasets.
5. Run the desired comparison.

Example Mathematica commands:

```mathematica
(* Ideal-state example: h = 2 *)
compareAllDepths[2, idealCorr, idealP]
allKS9[2, idealCorr, idealP]

(* Experimental data: one set of paths for each actual repetition *)
depthResults = allRepetitionsBoth[folders, 2];
separabilityResults = allRepeatedKS9[folders, 2];

depthResults
separabilityResults
```

In these calls, \(h=2\) controls the S36 swap-block selection and the \(k\)-separability construction. The S11 part does not take an \(h\) parameter.

## Interpreting the results

For each evaluated inequality, the notebook reports

\[
\Delta=2\mathrm{LEFT}-2\mathrm{RIGHT}.
\]

A positive `Delta` (`ViolationQ = True`) indicates a violation **under the assumptions and definitions of the corresponding witness bound**. Entanglement depth \(d\) and \(k\)-separability are distinct notions and should not be interpreted interchangeably.

The code also reports intermediate terms such as witness expectation values and swap-set correction terms to facilitate debugging and independent checks.



## Reproducibility notes

- Keep original measurement ordering unchanged; the code relies on its explicit row-to-Pauli mapping.
- Record all experiment repetition folders and acquisition settings.
- For measured results, retain the 512-outcome all-\(Z\) populations alongside the Pauli correlators from the same repetition.
- If you modify the swap-set definitions or normalization factors, rerun the checked example and compare intermediate values before drawing physical conclusions.
- Consider removing machine-specific paths or non-public experimental data before publishing.

## Citation

If you use this code in research, cite the paper containing the witness inequalities (in particular supplementary Eqs. S11 and S36 and Notations S3–S6). Add the full bibliographic citation and DOI here once confirmed.

## License

No license is specified in this README. If you intend to grant others permission to reuse or modify the code, add an explicit `LICENSE` file (for example, MIT, if appropriate for your project).
