<div align="right">
  <a href="./README.md"><img src="https://img.shields.io/badge/English-Active-0A66C2?style=flat-square" alt="English"></a>
  <a href="./README.pl.md"><img src="https://img.shields.io/badge/Polski-Przelacz_jezyk-dc143c?style=flat-square" alt="Polski"></a>
</div>

# Tabular Simplex Solver with Big-M Method in Python
### Operations Research & Mathematical Optimization

**Author:** Filip Gorczyca  
**Affiliation:** Economic Analytics, Krakow University of Economics (UEK)

---

## Project Overview
A from-scratch implementation of the **Tabular Simplex Algorithm** coupled with the **Big-M penalty method** in Python (`NumPy`), developed for Operations Research coursework. Rather than relying on black-box optimization solvers (`scipy.optimize` or `PuLP`), this project builds the complete matrix transformation and Gauss-Jordan pivoting engine to solve general Linear Programming (LP) problems with mixed inequality and equality constraints ($\le, \ge, =$).

---

## Mathematical Formulation & Algorithm Architecture

### 1. General Linear Programming Problem
The solver accepts any linear objective function (maximization or minimization, where $\min Z = -\max(-Z)$ ) subject to $m$ linear constraints and $n$ decision variables:

$$\max (\text{or } \min) \quad Z = \mathbf{c}^T \mathbf{x}$$

$$\text{subject to:} \quad \mathbf{A}\mathbf{x} \; (\le, \ge, =) \; \mathbf{b}, \quad \mathbf{x} \ge \mathbf{0}$$

### 2. Canonical Augmentation (`macierz` function)
To construct an initial basic feasible solution (identity submatrix $\mathbf{I}_m$), each constraint $i = 1, \dots, m$ is automatically augmented based on its relational operator:
* **$\le$ constraint:** Adds $+1$ **slack variable** $s_i \ge 0$ with objective coefficient $c_{s_i} = 0$ (enters initial basis).
* **$\ge$ constraint:** Subtracts $-1$ **surplus variable** $e_i \ge 0$ ($c_{e_i} = 0$) and adds $+1$ **artificial variable** $a_i \ge 0$ with a penalty coefficient $c_{a_i} = -M$ ($M = 1000$) in the objective function (enters initial basis).
* **$=$ constraint:** Adds $+1$ **artificial variable** $a_i \ge 0$ with penalty $c_{a_i} = -M$ (enters initial basis).

### 3. Iterative Simplex Tableau (`simpleks` function)
At each iteration, given the augmented tableau $\mathbf{T} = [\mathbf{A}_{\text{ext}} \mid \mathbf{b}]$ and current basis indices $B$:

1. **Shadow Prices & Opportunity Costs ($Z_j$):**

$$Z_j = \sum_{i=1}^{m} c_{B_i} T_{i,j}$$

2. **Optimality Criterion (Net Evaluation Row $C_j - Z_j$):**

$$\Delta_j = C_j - Z_j \le 0 \quad \forall j = 1, \dots, N_{\text{ext}}$$

If all $\Delta_j \le 0$, the current basis is optimal.

3. **Entering Variable Selection (Dantzig's Rule):**

$$k = \arg\max_{j} (C_j - Z_j)$$

4. **Leaving Variable Selection (Minimum Ratio Test):**

$$r = \arg\min_{i \,:\, T_{i,k} > 0} \left( \frac{T_{i, \text{rhs}}}{T_{i,k}} \right)$$

5. **Gauss-Jordan Row Operations:** Normalizes the pivot row $T_r \leftarrow T_r / T_{r,k}$ and eliminates column $k$ across all other rows $i \neq r$:

$$T_i \leftarrow T_i - T_{i,k} \cdot T_r$$

---

## Demonstration Case Study

The script solves the following mixed-constraint LP problem out of the box:

**Objective Function:**

$$\max Z = 3x_1 + 5x_2 + 2x_3$$

**Subject to Constraints:**
* $x_1 + x_2 + x_3 \le 7$
* $2x_1 + x_2 - x_3 \ge 2$
* $x_1 - x_2 + x_3 = 3$
* $x_1, x_2, x_3 \ge 0$

### Augmented Matrix & Execution Output
* **Decision variables:** $x_1, x_2, x_3$
* **Slack / Surplus variables:** $x_4$ (slack for constraint 1), $x_5$ (surplus for constraint 2)
* **Artificial variables ($M = 1000$):** $x_6$ (constraint 2), $x_7$ (constraint 3)
* **Initial Basis:** $B_0 = [x_4, x_6, x_7]$

```text
Macierz ograniczen po transformacji (A_ext):
[[ 1.  1.  1.  1.  0.  0.  0.]
 [ 2.  1. -1.  0. -1.  1.  0.]
 [ 1. -1.  1.  0.  0.  0.  1.]]
Wektor wspolczynnikow funkcji celu z karami Big-M (C_ext):
[    3.     5.     2.     0.     0. -1000. -1000.]
Baza startowa: ['x4', 'x6', 'x7']

Przebieg iteracji (wiersz wskaznikow optymalnosci Cj - Zj):
Cj - Zj: [3003.0, 5.0, 2.0, 0.0, -1000.0, 0.0, 0.0]
Cj - Zj: [0.0, -1496.5, 1503.5, 0.0, 501.5, -1501.5, 0.0]
Cj - Zj: [0.0, 7.0, 0.0, 0.0, 0.33, -1000.33, -1002.33]
Cj - Zj: [0.0, 0.0, 0.0, -3.5, 0.33, -1000.33, -998.83]
Cj - Zj: [0.0, 0.0, -1.0, -4.0, 0.0, -1000.0, -999.0]

Wszystkie wskazniki Cj - Zj <= 0 -> osiagnieto rozwiazanie optymalne.

ROZWIĄZANIE OPTYMALNE:
  x1 = 5.00
  x2 = 2.00
  x3 = 0.00
  x4 = 0.00
  x5 = 10.00
  x6 = 0.00
  x7 = 0.00

Optymalna wartosc funkcji celu (Z) = 25.00
```

### Verification of Optimal Solution
* **Optimal decision vector:** $x_1^* = 5, \; x_2^* = 2, \; x_3^* = 0$ (with surplus $x_5^* = 10$ and both artificial variables $x_6^* = x_7^* = 0$ driven out of the basis).
* **Objective value:** $Z^* = 3(5) + 5(2) + 2(0) = 25.00$.

---

## Quick Start
```bash
pip install numpy
python simplex_solver.py
```
