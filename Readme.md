# Week 3 – Z-Transform (MATLAB)

## How to Run
1. Open MATLAB.
2. Run the script **A_tasks_Z_transform.m**.
3. All results and plots are saved automatically in the `figures/` folder.

---

## A1 – Finite Sequences → Polynomials
- (i) X(z) = 1 + 2z⁻¹ + 5z⁻²  
- (ii) X(z) = 3z⁻¹ + 4z⁻³  
**ROC:** Entire z-plane except possibly z=0 or ∞ (no poles).

---

## A2 – Infinite Sequences & ROC
| Sequence | X(z) | ROC |
|-----------|------|-----|
| aⁿu[n], a=0.6 | 1/(1−0.6z⁻¹) | |z|>0.6 |
| (−0.8)ⁿu[n] | 1/(1+0.8z⁻¹) | |z|>0.8 |
| −(0.9)ⁿu[−n−1] | −z/(z−0.9) | |z|<0.9 |

ROC is outside poles for right‑sided and inside poles for left‑sided sequences.

---

## A3 – Linearity & Shifting
- Z{2x₁[n]−3x₂[n]} = 2/(1−0.5z⁻¹) − 3/(1+0.5z⁻¹)
- Z{x₁[n−3]} = z⁻³/(1−0.5z⁻¹)

---

## A4 – Inverse Z‑Transform
| X(z) | x[n] |
|------|------|
| 1/(1−0.7z⁻¹) | (0.7)ⁿu[n] |
| (1−0.5z⁻¹)/(1−0.8z⁻¹) | (0.8)ⁿu[n] − 0.5(0.8)ⁿ⁻¹u[n−1] |

---

## A5 – H(z), Poles/Zeros & Frequency Response
H(z) = (1−2.4z⁻¹+2.88z⁻²)/(1−0.8z⁻¹+0.64z⁻²)  
- Poles inside unit circle → stable.  
- Band‑pass‑like shape near ω≈π/3.  
Plots saved as **A5_zplane.png** and **A5_freq_response.png**.

---

### Files
- `A_tasks_Z_transform.m`
- `figures/` (plots + text summaries)

