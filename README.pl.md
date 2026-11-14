<div align="right">
  <a href="./README.md"><img src="https://img.shields.io/badge/English-Switch_language-0A66C2?style=flat-square" alt="English"></a>
  <a href="./README.pl.md"><img src="https://img.shields.io/badge/Polski-Aktywny-dc143c?style=flat-square" alt="Polski"></a>
</div>

# Algorytm Simpleks (Metoda Dużej M) w języku Python
### Badania Operacyjne i Optymalizacja Matematyczna

**Autor:** Filip Gorczyca  
**Afiliacja:** Analityka Gospodarcza, Uniwersytet Ekonomiczny w Krakowie (UEK)

---

## Opis projektu
Autorska implementacja **tablicowego algorytmu Simpleks** z wykorzystaniem **metody dużej M (Big-M)** w języku Python (`NumPy`), zrealizowana w ramach przedmiotu **Badania Operacyjne**. Zamiast korzystać z gotowych solverów optymalizacyjnych (`scipy.optimize` czy `PuLP`), skrypt implementuje od podstaw pełną transformację macierzową do postaci kanonicznej oraz iteracyjną eliminację Gaussa-Jordana dla zadań programowania liniowego (PL) z ograniczeniami mieszanymi ($\le, \ge, =$).

---

## Sformułowanie matematyczne i architektura algorytmu

### 1. Ogólna postać zadania Programowania Liniowego
Solver obsługuje zadania maksymalizacji lub minimalizacji (sprowadzane przez $\min Z = -\max(-Z)$ ) liniowej funkcji celu przy $m$ warunkach ograniczających i $n$ zmiennych decyzyjnych:

$$\max (\text{lub } \min) \quad Z = \mathbf{c}^T \mathbf{x}$$

$$\text{przy ograniczeniach:} \quad \mathbf{A}\mathbf{x} \; (\le, \ge, =) \; \mathbf{b}, \quad \mathbf{x} \ge \mathbf{0}$$

### 2. Sprowadzenie do postaci kanonicznej (funkcja `macierz`)
W celu utworzenia macierzy jednostkowej $\mathbf{I}_m$ (początkowego rozwiązania bazowego), każde ograniczenie $i = 1, \dots, m$ jest automatycznie rozszerzane w zależności od typu relacji:
* **Ograniczenie typu $\le$:** Dodaje $+1$ **zmienną luzu (dopełniającą)** $s_i \ge 0$ ze współczynnikiem $c_{s_i} = 0$ w funkcji celu (wchodzi do bazy startowej).
* **Ograniczenie typu $\ge$:** Odejmuje $-1$ **zmienną nadmiarową** $e_i \ge 0$ ($c_{e_i} = 0$) oraz dodaje $+1$ **zmienną sztuczną** $a_i \ge 0$ obciążoną w funkcji celu dużą karą $c_{a_i} = -M$ ($M = 1000$), która wchodzi do bazy startowej.
* **Ograniczenie typu $=$:** Dodaje $+1$ **zmienną sztuczną** $a_i \ge 0$ z karą $c_{a_i} = -M$ w funkcji celu (wchodzi do bazy startowej).

### 3. Iteracyjna tablica simpleksowa (funkcja `simpleks`)
W każdej iteracji dla rozszerzonej tablicy $\mathbf{T} = [\mathbf{A}_{\text{ext}} \mid \mathbf{b}]$ oraz bieżącego wektora indeksów bazowych $B$:

1. **Wyliczenie wiersza $Z_j$:**

$$Z_j = \sum_{i=1}^{m} c_{B_i} T_{i,j}$$

2. **Kryterium optymalności (wiersz wskaźników optymalności $C_j - Z_j$):**

$$\Delta_j = C_j - Z_j \le 0 \quad \forall j = 1, \dots, N_{\text{ext}}$$

Jeżeli wszystkie różnice $\Delta_j \le 0$, bieżące rozwiązanie bazowe jest optymalne.

3. **Wybór zmiennej wchodzącej do bazy (kolumna kluczowa):**

$$k = \arg\max_{j} (C_j - Z_j)$$

4. **Wybór zmiennej wychodzącej z bazy (wiersz kluczowy – test ilorazowy):**

$$r = \arg\min_{i \,:\, T_{i,k} > 0} \left( \frac{T_{i, \text{rhs}}}{T_{i,k}} \right)$$

5. **Eliminacja Gaussa-Jordana:** Podzielenie wiersza kluczowego przez element centralny ($T_r \leftarrow T_r / T_{r,k}$) i wyzerowanie pozostałych elementów kolumny kluczowej $k$ dla $i \neq r$:

$$T_i \leftarrow T_i - T_{i,k} \cdot T_r$$

---

## Przykład obliczeniowy

Wbudowany w skrypt przykład rozwiązuje następujące zadanie programowania liniowego z ograniczeniami mieszanymi:

**Funkcja celu:**

$$\max Z = 3x_1 + 5x_2 + 2x_3$$

**Warunki ograniczające:**
* $x_1 + x_2 + x_3 \le 7$
* $2x_1 + x_2 - x_3 \ge 2$
* $x_1 - x_2 + x_3 = 3$
* $x_1, x_2, x_3 \ge 0$

### Przebieg obliczeń i wynik z konsoli
* **Zmienne decyzyjne:** $x_1, x_2, x_3$
* **Zmienna luzu i nadmiarowa:** $x_4$ (luz w 1. równaniu), $x_5$ (nadmiar w 2. równaniu)
* **Zmienne sztuczne ($M = 1000$):** $x_6$ (2. równanie), $x_7$ (3. równanie)
* **Baza startowa:** $B_0 = [x_4, x_6, x_7]$

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

### Interpretacja rozwiązania
* **Optymalny wektor decyzyjny:** $x_1^* = 5, \; x_2^* = 2, \; x_3^* = 0$ (nadwyżka w drugim ograniczeniu wynosi $x_5^* = 10$, a obie zmienne sztuczne zostały poprawnie wyrugowane z bazy: $x_6^* = x_7^* = 0$).
* **Maksymalna wartość funkcji celu:** $Z^* = 3 \cdot 5 + 5 \cdot 2 + 2 \cdot 0 = 25,00$.

---

## Uruchomienie skryptu
```bash
pip install numpy
python simplex_solver.py
```
