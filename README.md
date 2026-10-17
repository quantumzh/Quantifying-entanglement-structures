# Multipartite Entanglement Analysis: 5-Qubit and 9-Qubit States

Wolfram Mathematica notebooks for entanglement-depth analysis of five-qubit and nine-qubit systems, together with separability tests for the five-qubit and nine-qubit cluster-state cases.

## Repository contents

| System | Notebook | Implemented analyses |
|---|---|---|
| 5-qubit | `5qubit/5qubit_merged_S36_S10_k3_k4_separability.nb` | Entanglement depth (S36, S10); 3- and 4-separability |
| 9-qubit GHZ | `9qubit-GHZ/9qubit_S36_S10_corrected_full.nb` | Entanglement depth only (S36, S10),j=2; separability modules is the same as those in 9qubit-Cluster/CL9_9Q_4Y_S36_S11_k3_to_k8.nb |Entanglement depth only (S36, S10),j=2
| 9-qubit GHZ| `GHZ9_4Y_h4_FULL_original_based.nb` |
| 9-qubit Cluster | `9qubit-Cluster/CL9_9Q_4Y_S36_S11_k3_to_k8.nb` | Entanglement depth (S36, S11); 3–8 separability |

```

## 1. Five-qubit analysis

The five-qubit notebook contains entanglement-depth criteria from supplementary Eqs. **S36** and **S10**. It additionally implements **3-separability and 4-separability** using weighted exchange-set constructions.

Example commands for the added separability routines:

```mathematica
paper5KSTenRepeats[1500, 3, 1]
paper5KSTenRepeats[1500, 4, 1]

(* h = 2 *)
paper5KSTenRepeats[1500, 3, 2]
paper5KSTenRepeats[1500, 4, 2]
```

These functions reuse the five-qubit notebook's existing experimental-data handling and exchange-set machinery.

## 2. Nine-qubit GHZ entanglement-depth analysis

**Notebook:** `9qubit_S36_S10_corrected_full.nb` for j=2 inequalities;`GHZ9_4Y_h4_FULL_original_based.nb' for j=4 inequalities;

This notebook implements entanglement-depth bounds for a nine-qubit GHZ-state analysis. It uses supplementary **Eq. S36, Eq. S10 and Eq. S11**, including partition construction, weighted exchange sets, and multiplicity handling. 

The uploaded notebook includes:


- Experimental evolution times and repetition times.
- Measurement-row configuration, data loading, and witness expectation values.
- Depth-bound construction with supplementary Eqs. S36 and S10.
- A unified criterion, summaries, example calls, and checks of partition enumeration and measurement-row layout.

The default local data root is configured with `q9Root`. Change it to your actual experimental directory before running. The notebook uses `q9ExperimentFolder[time]` to name evolution-time folders such as `9Q_GHZ_1500ns`, along with the configured readout subdirectory and `basis_measurement` file prefix.

Relevant Mathematica functions in the notebook include `q9Case`, `q9Summary`, and `q9AllTimes`. Evaluate the initialization cells and consult the notebook's **“6. Examples”** and **“8. Run experiment calculations”** sections for their exact argument conventions before launching a complete dataset run.


## 3. Nine-qubit cluster-state analysis

**Notebook:** `CL9_9Q_4Y_S36_S11_k3_to_k8.nb`

The cluster-type state used in this notebook is

\[
\begin{aligned}
|CL_9\rangle=\frac{1}{2\sqrt{2}}\big(&|000000000\rangle+|000000111\rangle
+|000111011\rangle+|000111100\rangle\\
&+|111000011\rangle+|111000100\rangle
+|111111000\rangle+|111111111\rangle\big).
\end{aligned}
\]

The four complementary pairs are

\[
\Omega=\{(1,8),(2,7),(3,6),(4,5)\},
\]

where the indices enumerate the eight computational-basis terms above.

Its analysis covers:

- **Entanglement-depth criteria:** S36 for \(d=1,\ldots,4\), S11 for \(d=5,\ldots,9\).
- **Separability criteria:** \(k=3,\ldots,8\).
- **Four-\(Y\) witnesses:** 126 measured Pauli strings with four \(Y\) and five \(X\) factors.
- Weighted swap-set operations, computational-basis populations, and per-repetition comparisons.

Example Mathematica commands:

```mathematica
(* Ideal-state example *)
compareAllDepths[2, idealCorr, idealP]
allKS9[2, idealCorr, idealP]

(* Experimental results after configuring the data folders *)
depthResults = allRepetitionsBoth[folders, 2];
separabilityResults = allRepeatedKS9[folders, 2];
```

The cluster notebook uses Pauli codes `1 = X`, `2 = Y`, `3 = Z`; its 1-based measurement row 340 is the all-\(Z\) basis. It expects `basis_measurement0.txt` through `basis_measurement340.txt` for each repetition, with 512 outcomes per measurement setting.

## Interpretation of LEFT and RIGHT

Each implemented witness compares twice of the left-hand side and right-hand side of its relevant supplementary inequality:

\[
\Delta = \mathrm{LEFT} - \mathrm{RIGHT}.
\]

When \(\Delta>0\), the corresponding bound is violated **subject to the assumptions of that specific criterion and the correctness of its implementation**. Entanglement depth and \(k\)-separability are different physical classifications.

The S10/S11/S36 branches and the \(k\)-separability extensions should not be treated as interchangeable just because their output tables share `LEFT`, `RIGHT`, and `Delta` labels.

## Running and reproducibility

1. Open the relevant `.nb` file in Wolfram Mathematica.
2. Configure its experiment root, time points, repetition directories, and file layout.
3. Evaluate notebook cells in order.
4. Run a verified reference example before processing a complete experimental dataset.
5. Record `LEFT`, `RIGHT`, the difference, and the setting parameters for reproducibility.

All experiment paths in the notebooks refer to local files; data are not included automatically. Avoid publishing confidential data and machine-specific file paths unintentionally.



## References and license

The notebooks refer to supplementary **Eqs. S10, S11, S36 and Notations S3–S6**. Please add the source paper's complete bibliographic reference and DOI.

No redistribution license is asserted here. Add a `LICENSE` file if you want to grant reuse permissions.
