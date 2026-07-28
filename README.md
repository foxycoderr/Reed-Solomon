# Reed–Solomon Codec (`fullRS.py`)

A from-scratch Reed–Solomon encoder and decoder over the prime field $\mathbb{F}_{251}$, implemented in Python with NumPy and SymPy.

---

## 1. Preface: What is Reed–Solomon?

Reed–Solomon (RS) codes are a classical family of **error-correcting codes**. They take a message of $k$ symbols and expand it into a longer codeword of $n$ symbols ($n > k$) so that, after transmission through a noisy channel, a decoder can recover the original message even if some symbols were corrupted.

### The idea in one sentence

Treat the message as a polynomial, systematically force the transmitted polynomial to have known roots, then use evaluations at those roots (the *syndromes*) to locate and remove errors.

### Parameters

| Symbol | Meaning |
|--------|---------|
| $n$ | Codeword length (number of transmitted symbols) |
| $k$ | Message length (number of information symbols) |
| $n - k$ | Number of parity / redundancy symbols |
| $t = \lfloor (n-k)/2 \rfloor$ | Maximum number of symbol errors the code can correct |

This project works over **$\mathbb{F}_{251}$** (integers modulo the prime 251). Every coefficient, evaluation, and matrix operation is done mod 251. A fixed field element

$$
\alpha = 6 \in \mathbb{F}_{251}
$$

is used as a primitive generator for building the code (the multiplicative order of 6 mod 251 is 250, so distinct powers $\alpha^i$ can label positions in blocks shorter than 250).

### Encoding (high level)

1. Interpret the message $m_0, m_1, \ldots, m_{k-1}$ as a polynomial

$$
m(x) = m_0 + m_1 x + \cdots + m_{k-1} x^{k-1}.
$$

2. Build the **generator polynomial**

$$
g(x) = \prod_{j=1}^{n-k} \bigl(x - \alpha^j\bigr).
$$

3. Form the codeword polynomial

$$
c(x) = m(x)\, g(x).
$$

By construction, $c(\alpha^j) = 0$ for $j = 1, \ldots, n-k$.

This is a **non-systematic** RS construction: the message is not visible as a contiguous block of the codeword; recovering it later requires dividing the corrected codeword by $g(x)$.

### Channel / errors

The receiver sees

$$
r(x) = c(x) + e(x),
$$

where the error polynomial $e(x)$ is sparse: if there are $\nu \le t$ errors at positions $i_1, \ldots, i_\nu$ with magnitudes $Y_1, \ldots, Y_\nu$, then

$$
e(x) = \sum_{r=1}^{\nu} Y_r\, x^{i_r}.
$$

The **error locators** are $X_r = \alpha^{i_r}$.

### Decoding (high level)

This codebase follows a **Peterson–Gorenstein–Zierler**-style decoder:

1. **Syndromes.** Evaluate the received word at the known roots:

$$
S_j = r(\alpha^j), \quad j = 1, \ldots, n-k.
$$

If all syndromes are zero, there are no errors (within the code’s design).

2. **Error-locator polynomial (ELP).** Solve a linear system built from the syndromes for

$$
\Lambda(x) = 1 + \Lambda_1 x + \cdots + \Lambda_\nu x^\nu = \prod_{r=1}^{\nu} (1 - X_r x).
$$

The roots of $\Lambda$ are the **inverses** $X_r^{-1}$.

3. **Chien search.** Test every field element to find those roots, then invert them to recover the locators $X_r$.

4. **Error magnitudes.** Solve another linear system (Forney / Vandermonde style) for the $Y_r$.

5. **Correct.** Subtract the reconstructed $e(x)$ from $r(x)$ to recover $c(x)$.

6. **Unencode.** Divide $c(x)$ by $g(x)$ to recover $m(x)$.

---

## 2. How to use these scripts

### Requirements

- Python 3
- NumPy
- SymPy

### Run the built-in demo

`fullRS.py` ends with an example call:

```python
run(20, 10, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 3)
```

That means: block length $n = 20$, message length $k = 10$, message symbols `[1..10]`, inject $e = 3$ errors, then encode → corrupt → decode → recover the message.

```bash
python fullRS.py
```

You will see verbose prints for syndromes, ELP search, Gaussian elimination, and success / failure checks.

### Call the API yourself

Import and use `run`, or call the stages separately:

```python
from fullRS import encode, transmit, decode, find_message, run

# One-shot pipeline (encode, corrupt, decode, verify)
run(n=20, k=10, message=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], e=3)

# Or step by step:
codeword = encode(20, 10, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
received_poly, err_locs, err_vals = transmit(codeword, e=3)
corrected_coeffs = decode(20, 10, received_poly)
message = find_message(corrected_coeffs, 20, 10)
```

### Constraints to keep in mind

- Message symbols must be integers in $\{0, 1, \ldots, 250\}$.
- Message length must equal $k$.
- Number of injected errors $e$ should satisfy $e \le t = \lfloor (n-k)/2 \rfloor$ for reliable correction.
- Prefer $n < 250$ so that position labels $\alpha^0, \ldots, \alpha^{n-1}$ stay distinct.

### Changing the demo

Edit the final line of `fullRS.py`, for example:

```python
run(16, 8, [3, 1, 4, 1, 5, 9, 2, 6], 2)
```

Or comment that line out and drive the functions from another script / REPL.

---

## 3. Complete function documentation

Globals used throughout:

```python
alpha = 6    # primitive element for generator roots / locators
mod = 251    # field characteristic (prime)
```

NumPy is imported as `isaac_newton` (matrices, `floor`, modular reduction helpers). Polynomials use `numpy.polynomial.polynomial.Polynomial` (via `from numpy import polynomial`). SymPy `Matrix` is used for modular matrix inversion.

The pipeline orchestrated by `run` is:

```text
encode → transmit → decode → find_message
                 ↘ (inside decode)
                   find_syndromes → find_elp → solve_elp
                   → find_error_polynomial_coefficients_Yr
                   → generate_error_polynomial_coefficients
                   → subtract errors
```

---

### `polymul_mod(poly1, poly2, modulus=251)`

**Role:** Multiply two NumPy `Polynomial` objects with **integer coefficients reduced modulo `modulus`**, avoiding float overflow / precision issues from naïve NumPy polynomial multiplication.

**How it works:**

1. Round each coefficient to an int and reduce mod `modulus`.
2. Compute the Cauchy product (convolution) with running `% modulus`.
3. Return a new `Polynomial` whose coefficients are those modular results.

**Used by:** `encode`, `find_message` (building $g(x)$), and anywhere a modular product of field polynomials is required.

---

### `reduce_poly_mod_251(poly)`

**Role:** Take a `Polynomial` and return a new one whose coefficients are integers in $\{0, \ldots, 250\}$.

**How it works:** Read `.coef`, cast to int, reduce each entry mod 251, rebuild the polynomial.

**Used by:** `find_message`, after polynomial subtraction / multiplication steps during polynomial division over $\mathbb{F}_{251}$.

---

### `encode(n, k, message)`

**Role:** Encode a length-$k$ message into a length-$\approx n$ codeword polynomial $c(x) = m(x)\, g(x)$.

**Inputs:**

| Arg | Meaning |
|-----|---------|
| `n` | Designed block length |
| `k` | Message length |
| `message` | List of $k$ integers in $0..250$ |

**Steps:**

1. Validate `len(message) == k` and that every symbol is an int in range; otherwise print an error and `exit(1)`.
2. Build $m(x)$ from `message` (coefficient of $x^i$ is `message[i]`).
3. Build the generator

$$
g(x) = (x - \alpha)(x - \alpha^2)\cdots(x - \alpha^{n-k})
$$

by starting from $[-\alpha \bmod 251,\; 1]$ (i.e. $x - \alpha$) and successively multiplying by $(x - \alpha^i)$ via `polymul_mod` for $i = 2, \ldots, n-k$.

4. Compute $c(x) = g(x)\, m(x)$ with `polymul_mod`.
5. Print the clean codeword and return it as a `Polynomial`.

**Output:** Codeword polynomial (coefficients are the transmitted symbols in positional order $c_0, c_1, \ldots$).

**Note:** This is non-systematic encoding. The returned polynomial degree is typically $n-1$ when leading coefficients are nonzero.

---

### `transmit(message, e)`

**Role:** Simulate a noisy channel: inject (up to) $e$ symbol errors into a codeword.

**Inputs:**

| Arg | Meaning |
|-----|---------|
| `message` | Codeword (`Polynomial` or coefficient-iterable) |
| `e` | Number of error-injection attempts |

**Steps:**

1. Convert the codeword to a mutable coefficient list.
2. For each of $e$ trials:
   - Pick a random index in `[0, len(message)-1]`.
   - Draw a random magnitude in `1..250` and record it.
   - Record the index as an error location.
   - Update that coefficient (currently: after adding the magnitude, the coefficient is set to `index % 251`, i.e. the stored “corrupted” symbol is the position index itself).
3. Wrap the list back into a `Polynomial`.

**Output:** A list of three items:

```text
[received_polynomial, error_locations, error_values]
```

- `error_locations` / `error_values` are the *attempted* injections (useful for debugging). The same index may be chosen more than once.

**Used by:** `run`.

---

### `eval_mod_poly(coeffs, x, p=251)`

**Role:** Evaluate a polynomial given by coefficient list `coeffs` at point $x$, entirely in modular arithmetic (Horner’s method).

**Formula:** For $f(x) = c_0 + c_1 x + \cdots + c_d x^d$,

$$
f(x) \equiv (((c_d)\, x + c_{d-1})\, x + \cdots)\, x + c_0 \pmod{p}.
$$

**Why it exists:** NumPy’s `Polynomial.__call__` uses floating-point arithmetic. Modular Horner stays exact in $\mathbb{F}_{251}$ for syndromes and for Chien search over all field elements.

**Used by:** `find_syndromes`, `solve_elp`.

---

### `find_syndromes(n, k, r_polynomial)`

**Role:** Compute the $n-k$ syndromes of the received word.

**Inputs:** Block parameters $n$, $k$, and received polynomial $r_polynomial$.

**Steps:**

1. Extract integer coefficients from the received polynomial.  
   *(Implementation detail: it reads `r_polynomial.coef[0]` as an iterable. That matches `decode` / `run`, where the received object is a `Polynomial` wrapped again by `Polynomial(...)`, so `.coef[0]` is itself the inner polynomial / coefficient source.)*
2. For each $j = 1, \ldots, n-k$:

$$
S_j = r(\alpha^j) \bmod 251,
$$

computed with `eval_mod_poly`.

3. Print each syndrome and return the list $[S_1, S_2, \ldots, S_{n-k}]$ (0-based Python list: index `0` holds $S_1$).

**Interpretation:** If $r = c + e$ and $c(\alpha^j) = 0$, then $S_j = e(\alpha^j)$. Nonzero syndromes indicate errors.

---

### `find_elp(syndromes, t)`

**Role:** Find the **error-locator polynomial** $\Lambda(x)$ using the Peterson–Gorenstein–Zierler approach: try candidate error counts $p = t, t-1, \ldots$ until a consistent invertible system is found.

**Inputs:**

| Arg | Meaning |
|-----|---------|
| `syndromes` | List $[S_1,\ldots,S_{n-k}]$ from `find_syndromes` |
| `t` | Designed correction capability $\lfloor(n-k)/2\rfloor$ |

**For each candidate $p$:**

1. Build the $p \times p$ syndrome matrix

$$
M =
\begin{bmatrix}
S_1 & S_2 & \cdots & S_p \\
S_2 & S_3 & \cdots & S_{p+1} \\
\vdots & & \ddots & \vdots \\
S_p & S_{p+1} & \cdots & S_{2p-1}
\end{bmatrix}.
$$

2. Try to invert $M$ modulo 251 with SymPy (`inv_mod`). If non-invertible, try a smaller $p$.
3. Form the right-hand side

$$
\mathbf{b} = \bigl(-S_{p+1},\; -S_{p+2},\; \ldots,\; -S_{2p}\bigr)^\top.
$$

4. Solve $M \mathbf{v} = \mathbf{b}$ for $\mathbf{v}$.
5. **Consistency check:** verify $M \mathbf{v}$ matches $\mathbf{b}$ (printed as “Hooray!…”). If the matrix was invertible but the check fails, a smaller true error count is assumed and the loop continues.
6. Because of the column ordering of $M$ relative to the classical key equation

$$
S_{j+p} + \Lambda_1 S_{j+p-1} + \cdots + \Lambda_p S_j = 0,
$$

the raw solution vector is $(\Lambda_p,\ldots,\Lambda_1)$. The code **reverses** it, then prepends $1$, yielding

$$
\Lambda(x) = 1 + \Lambda_1 x + \cdots + \Lambda_p x^p.
$$

**Output:** NumPy `Polynomial` for $\Lambda(x)$.

**Used by:** `decode`.

---

### `solve_elp(error_correcting_polynomial)`

**Role:** Find the roots of the ELP over $\mathbb{F}_{251}$ (Chien-style exhaustive search), then convert each root $X_r^{-1}$ into the error locator $X_r$.

**Steps:**

1. Read the ELP coefficients as integers from `error_correcting_polynomial.coef`.
2. For every field element $i = 0, 1, \ldots, 250$, evaluate $\Lambda(i)$ with `eval_mod_poly` (modular Horner).
3. If the value is $0$, append $i$ as a root $X_r^{-1}$.
4. Invert each root with Fermat’s little theorem via `pow(a, 251-2, 251)` (vectorized), obtaining

$$
X_r = (X_r^{-1})^{-1} = \alpha^{i_r}.
$$

**Output:** List of integer locators `Xr`.

**Used by:** `decode`.

---

### `find_error_polynomial_coefficients_Yr(syndromes, Xr, p)`

**Role:** Given the error locators $X_r$, solve for the error **magnitudes** $Y_r$.

**Mathematical system:** For $j = 1, \ldots, p$ (here $p =$ number of locators),

$$
\sum_{r=1}^{p} Y_r\, X_r^{\,j} = S_j.
$$

**Steps:**

1. Build an augmented matrix with rows

$$
\bigl[X_1^j,\; X_2^j,\; \ldots,\; X_p^j \;\big|\; S_j\bigr], \quad j = 1,\ldots,p.
$$

2. Run **Gaussian elimination modulo 251**:
   - Scale the pivot row by the modular inverse of the pivot.
   - Clear the pivot column in all other rows by adding a suitable multiple of the pivot row (using `251 - mirror` so addition implements subtraction mod 251).
3. Read the last column of the reduced system as $Y_1, \ldots, Y_p$.

**Output:** List `Yr` of error magnitudes (same order as `Xr`).

**Used by:** `decode` (with `p = len(Xr)`).

---

### `generate_error_polynomial_coefficients(Xr, Yr, n)`

**Role:** Turn locator / magnitude pairs into an explicit length-$n$ coefficient list for $e(x)$.

**Steps:**

1. For each locator $X$, find the position index $i$ such that $\alpha^i \equiv X \pmod{251}$ (search $i = 0..249$).
2. Allocate a zero list of length $n$.
3. Set `error_polynomial_coefficients[i_r] = Y_r` for each pair.

**Output:** Coefficients of

$$
e(x) = \sum_r Y_r\, x^{i_r}.
$$

**Used by:** `decode`.

---

### `decode(n, k, received)`

**Role:** Full decoder from received word to **corrected codeword coefficients** (not yet the original message).

**Inputs:** $n$, $k$, and `received` (typically the `Polynomial` returned as the first element of `transmit`).

**Pipeline:**

1. Wrap `received` as `r_polynomial = Polynomial(received)`.
2. `syndromes = find_syndromes(n, k, r_polynomial)`.
3. $t = \lfloor (n-k)/2 \rfloor$.
4. $\Lambda =$ `find_elp(syndromes, t)`.
5. `Xr = solve_elp(Λ)`.
6. `Yr = find_error_polynomial_coefficients_Yr(syndromes, Xr, len(Xr))`.
7. `e_coeffs = generate_error_polynomial_coefficients(Xr, Yr, n)`.
8. Correct symbolwise:

$$
\hat{c}_i = (r_i - e_i) \bmod 251.
$$

**Output:** List of $n$ integers — the recovered codeword coefficients $\hat{c}(x)$.

**Note:** Prints a warning if `len(received) != n` but does not abort on that check alone.

---

### `find_message(decoded_coeff, n, k)`

**Role:** Recover the original message polynomial $m(x)$ from a corrected codeword by dividing by the generator $g(x)$ in $\mathbb{F}_{251}$.

**Why needed:** Encoding used $c(x) = m(x)\, g(x)$ (non-systematic), so

$$
m(x) = c(x) / g(x)
$$

when $c$ is a valid codeword.

**Steps:**

1. Rebuild the same generator $g(x)$ as in `encode`.
2. Perform polynomial division of the decoded coefficient list by $g(x)$ over $\mathbb{F}_{251}$, working from high degree downward:
   - Leading coefficient of $g$ is 1, so each quotient coefficient equals the current leading coefficient of the dividend.
   - Form the corresponding monomial times $g$, subtract mod 251 (`reduce_poly_mod_251`), and repeat for $n-k$ steps (enough to peel off the degree contributed by $g$).
3. Collect quotient coefficients into `message_coeffs`.

**Output:** List of $k$ message symbols (ideally matching the original message).

**Used by:** `run`, after a successful codeword correction.

---

### `run(n, k, message, e)`

**Role:** End-to-end demonstration / integration test of the whole codec.

**Steps:**

1. `codeword = encode(n, k, message)` — print coefficients.
2. `transmitted = transmit(codeword, e)` — print received poly, true error locations, and attempted error values.
3. `decoded = decode(n, k, transmitted[0])` — corrected codeword coeffs.
4. Compare `decoded` to the original codeword coefficients:
   - match → print success;
   - mismatch → print failure and `exit(1)`.
5. `message_dec = find_message(decoded, n, k)`.
6. Compare to the original `message`:
   - match → print success and `exit(0)`;
   - mismatch → print bug warning.

**Default example** (bottom of `fullRS.py`):

```python
run(20, 10, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 3)
```

---

## Appendix: Symbol map (code ↔ math)

| Code name | Mathematical meaning |
|-----------|----------------------|
| `alpha` | Primitive element $\alpha = 6$ |
| `mod` | Field size / modulus $251$ |
| `message` / `m(x)` | Information polynomial |
| `generator_polynomial` / `g(x)` | Generator $\prod (x-\alpha^j)$ |
| `codeword` / `c(x)` | Encoded word $m\cdot g$ |
| `syndromes[j-1]` | $S_j = r(\alpha^j)$ |
| `error_correcting_polynomial` | ELP $\Lambda(x)$ |
| `Xr` | Error locators $X_r = \alpha^{i_r}$ |
| `Yr` | Error magnitudes $Y_r$ |
| `ir` | Error positions $i_r$ |
| `t` | Designed correction radius $\lfloor(n-k)/2\rfloor$ |
| `p` | Assumed / detected number of errors $\nu$ |
