import Mathlib

/-!
# Parity lemma for walks on bipartite graphs (p22, Sec. "parity trap")

If every nonzero entry of a transition (or adjacency) matrix `P` joins vertices of different
colour under a two-colouring `c`, then `(P ^ t) i j = 0` unless `c j = xor (c i) (t odd)`.
Consequence: every diagonal entry of an odd power is zero, so the return probability
`q(t) = tr(P^t)/N` of a non-lazy walk on a bipartite graph vanishes at odd `t`.
-/

open Matrix

variable {n : Type*} [Fintype n] [DecidableEq n] {R : Type*} [Semiring R]

/-- Entries of `P ^ t` vanish between vertices whose colours do not match the parity of `t`. -/
theorem bipartite_pow_apply_eq_zero (c : n → Bool) (P : Matrix n n R)
    (hP : ∀ i j, c i = c j → P i j = 0) :
    ∀ (t : ℕ) (i j : n), c j ≠ xor (c i) (decide (t % 2 = 1)) → (P ^ t) i j = 0 := by
  intro t
  induction t with
  | zero =>
    intro i j h
    have hij : i ≠ j := by
      rintro rfl; simp at h
    simp [Matrix.one_apply_ne hij]
  | succ t ih =>
    intro i j h
    rw [pow_succ, Matrix.mul_apply]
    apply Finset.sum_eq_zero
    intro k _
    by_cases hk : c k = xor (c i) (decide (t % 2 = 1))
    · have hpar : decide ((t + 1) % 2 = 1) = !decide (t % 2 = 1) := by
        by_cases ht : t % 2 = 1
        · have : (t + 1) % 2 ≠ 1 := by omega
          simp [ht, this]
        · have : (t + 1) % 2 = 1 := by omega
          simp [ht, this]
      rw [hpar] at h
      have hkj : c k = c j := by
        rw [hk]; revert h
        cases c i <;> cases c j <;> cases decide (t % 2 = 1) <;> simp
      rw [hP k j hkj, mul_zero]
    · rw [ih i k hk, zero_mul]

/-- Every diagonal entry of an odd power of a bipartite transition matrix is zero. -/
theorem bipartite_odd_pow_diag (c : n → Bool) (P : Matrix n n R)
    (hP : ∀ i j, c i = c j → P i j = 0) (t : ℕ) (ht : t % 2 = 1) (i : n) :
    (P ^ t) i i = 0 :=
  bipartite_pow_apply_eq_zero c P hP t i i (by simp [ht])

/-- Hence the trace (sum of return probabilities over start nodes) of an odd power is zero. -/
theorem bipartite_odd_pow_trace (c : n → Bool) (P : Matrix n n R)
    (hP : ∀ i j, c i = c j → P i j = 0) (t : ℕ) (ht : t % 2 = 1) :
    (P ^ t).trace = 0 :=
  Finset.sum_eq_zero fun i _ => bipartite_odd_pow_diag c P hP t ht i
