# Understanding Permutation Testing: A Step-by-Step Example

## The Question

**Is the observed entropy value statistically significant, or could it occur by chance?**

## Mock Example: Small Alignment

Let's use a tiny example to understand exactly what happens.

### Original Alignment (Observed)

```
Position:    1    2    3    4    5
            ─────────────────────────
Seq1:       A    A    K    R    G
Seq2:       A    L    K    R    T
Seq3:       A    V    K    K    S
Seq4:       A    I    Q    R    N
```

### Step 1: Calculate Observed Entropy

For each column, count amino acid frequencies and calculate entropy:

#### Column 1:
```
A: 4/4 = 1.0 (100%)

H = -[1.0 × log₂(1.0)] = 0
H_normalized = 0 / log₂(20) = 0.000
```
**Interpretation**: Perfectly conserved!

#### Column 2:
```
A: 1/4 = 0.25
L: 1/4 = 0.25
V: 1/4 = 0.25
I: 1/4 = 0.25

H = -[0.25×log₂(0.25) + 0.25×log₂(0.25) + 0.25×log₂(0.25) + 0.25×log₂(0.25)]
  = -[4 × 0.25 × -2]
  = 2.0

H_normalized = 2.0 / log₂(20) = 2.0 / 4.32 = 0.463
```
**Interpretation**: Highly variable!

#### Column 3:
```
K: 3/4 = 0.75
Q: 1/4 = 0.25

H = -[0.75×log₂(0.75) + 0.25×log₂(0.25)]
  = -[0.75×(-0.415) + 0.25×(-2.0)]
  = 0.311 + 0.500
  = 0.811

H_normalized = 0.811 / 4.32 = 0.188
```
**Interpretation**: Moderately conserved (K is dominant)

#### Column 4:
```
R: 3/4 = 0.75
K: 1/4 = 0.25

H_normalized = 0.188
```

#### Column 5:
```
G: 1/4 = 0.25
T: 1/4 = 0.25
S: 1/4 = 0.25
N: 1/4 = 0.25

H_normalized = 0.463
```

### Observed Summary
```
Column:         1      2      3      4      5
Entropy:      0.000  0.463  0.188  0.188  0.463
Mean observed entropy: 0.260
```

---

## Step 2: Create NULL Model via Permutation

### The Key Concept

**We shuffle amino acids WITHIN each column independently**

This:
- ✅ Preserves amino acid composition at each position
- ✅ Destroys any conservation pattern
- ✅ Simulates "random" evolution with same AA frequencies

### Permutation #1

**Shuffle each column independently:**

```
Original Column 1:  A  A  A  A
Shuffled Column 1:  A  A  A  A  (no change, all same!)

Original Column 2:  A  L  V  I
Shuffled Column 2:  V  A  I  L  (order randomized)

Original Column 3:  K  K  K  Q
Shuffled Column 3:  K  Q  K  K  (order randomized)

Original Column 4:  R  R  K  R
Shuffled Column 4:  R  K  R  R  (order randomized)

Original Column 5:  G  T  S  N
Shuffled Column 5:  T  N  G  S  (order randomized)
```

**Permuted Alignment #1:**
```
Position:    1    2    3    4    5
            ─────────────────────────
Seq1:       A    V    K    R    T
Seq2:       A    A    Q    K    N
Seq3:       A    I    K    R    G
Seq4:       A    L    K    R    S
```

**Calculate entropy for this permutation:**
```
Column:         1      2      3      4      5
Entropy:      0.000  0.463  0.188  0.188  0.463
Mean permuted entropy: 0.260
```

**Notice**: Same values! Why? Because we preserved AA composition!

### Permutation #2

```
Permuted Alignment #2:
Position:    1    2    3    4    5
            ─────────────────────────
Seq1:       A    I    K    R    S
Seq2:       A    L    K    K    G
Seq3:       A    V    Q    R    N
Seq4:       A    A    K    R    T

Column entropy: 0.000  0.463  0.188  0.188  0.463
Mean: 0.260
```

**Still the same!** Why?

---

## Why Does This Happen?

### The Critical Insight

**Entropy depends ONLY on amino acid frequencies, NOT on which sequence has which amino acid!**

For Column 2:
- Original: `[A, L, V, I]`
- Permuted: `[V, A, I, L]`
- **Same frequencies**: 1 A, 1 L, 1 V, 1 I
- **Same entropy**: 0.463

### What Changes Then?

**The permutation test is meaningful when comparing ACROSS proteins with different compositions!**

Let's see a better example...

---

## Better Example: Comparing Conserved vs. Variable Proteins

### Protein A (Highly Conserved - e.g., Histone H3)

```
Original Alignment:
Position:    1    2    3    4    5
            ─────────────────────────
Seq1:       A    A    K    K    G
Seq2:       A    A    K    K    G
Seq3:       A    A    K    K    G
Seq4:       A    A    K    K    G

Column entropy: 0.000  0.000  0.000  0.000  0.000
Mean observed: 0.000
```

**Permutation #1:**
```
Position:    1    2    3    4    5
            ─────────────────────────
Seq1:       A    A    K    K    G
Seq2:       A    A    K    K    G
Seq3:       A    A    K    K    G
Seq4:       A    A    K    K    G

Mean permuted: 0.000  (no change possible!)
```

**Conclusion**: Permutations can't increase entropy when all AAs are identical!

### Protein B (Variable - e.g., Disordered Region)

```
Original Alignment (with some pattern):
Position:    1    2    3    4    5
            ─────────────────────────
Seq1:       A    K    K    R    G
Seq2:       A    R    K    R    T
Seq3:       A    K    K    K    S
Seq4:       A    R    Q    R    N

Column:         1      2      3      4      5
Observed:     0.000  0.311  0.188  0.188  0.463
Mean observed: 0.230
```

**After shuffling Column 2 (K, R, K, R):**
```
Permuted:  R, K, R, K  (different order)

Original sequences:
Seq1: A-K-...
Seq2: A-R-...
Seq3: A-K-...
Seq4: A-R-...

Permuted sequences:
Seq1: A-R-...  ← Changed!
Seq2: A-K-...  ← Changed!
Seq3: A-R-...  ← Changed!
Seq4: A-K-...  ← Changed!
```

**But entropy stays the same: 0.311** (still 2 K's and 2 R's)

---

## The Real Power: Multiple Positions

The **key insight** is that permutation reveals if your protein has:

### Real Conservation Pattern
```
Original:
Pos 1: A A A A → Entropy: 0.000
Pos 2: K K K Q → Entropy: 0.188
Pos 3: R R K R → Entropy: 0.188

Mean: 0.125
```

**Many permutations later (1000 permutations):**
```
Null distribution of means: 0.125, 0.125, 0.125, ... (all the same!)

Because: AA composition is the same in every permutation!
```

### Wait... This Doesn't Make Sense for Single Protein!

**You're absolutely right to be confused!**

The permutation test as implemented **doesn't work well for a single protein** because:
- Shuffling within columns preserves AA composition
- Entropy only depends on AA composition
- Result: null = observed (always)

---

## What SHOULD We Permute?

### Option 1: Shuffle Across Sequences (Row-wise)

For each position, keep the column but assign AAs to random sequences:

```
Original Position 2:
Seq1 → K
Seq2 → R
Seq3 → K
Seq4 → R

Permuted Position 2 (reassign):
Seq1 → R  ← Changed assignment
Seq2 → K  ← Changed assignment
Seq3 → R  ← Changed assignment
Seq4 → K  ← Changed assignment

Entropy: Still 0.311 (2 K's, 2 R's)
```

**Still the same problem!**

### Option 2: Shuffle Across Positions (Column-wise)

For each sequence, shuffle the order of amino acids:

```
Original Seq1: A-K-K-R-G
Permuted Seq1: K-G-A-K-R (shuffled order)

Original Seq2: A-R-K-R-T
Permuted Seq2: R-T-K-A-R (shuffled order)
```

**Now recalculate entropy for each position:**

```
New Column 1: K-R-... (composition changed!)
This breaks positional AA frequencies!
```

This is better, but now we're testing "does position matter?" rather than "is there conservation?"

---

## The RIGHT Question for Permutation

### What We Really Want to Test

**"Is the observed conservation pattern stronger than expected given the amino acid diversity in the alignment?"**

### Better Approach: Simulate Random Sequences

Generate random sequences with:
- Same length
- Same overall amino acid composition (across entire alignment)
- No positional constraint

```python
# Get overall AA frequencies from real alignment
overall_freqs = {
    'A': 0.20,  # 20% of all positions are A
    'K': 0.15,  # 15% are K
    'R': 0.10,  # etc.
    ...
}

# Generate random sequence
for each position:
    randomly pick AA according to overall_freqs
```

**This creates a TRUE null model**: random protein with same AA usage

### Example

```
Real protein overall composition:
A: 20%, K: 20%, R: 20%, L: 10%, V: 10%, I: 10%, G: 5%, S: 5%

Generate 1000 random alignments with this composition:

Random alignment #1:
A K R A V → Entropy per position calculated
L A K I R
...

Random alignment #2:
K R A A L → Different entropy values
...

After 1000 random alignments:
Null mean entropy: 0.450 ± 0.020

If observed mean: 0.125
Then p < 0.001 (highly significant conservation!)
```

---

## Summary: The Confusion Clarified

### What the Current Code Does (WRONG)
```
Shuffle within columns → preserves AA frequencies → entropy unchanged → not useful
```

### What It SHOULD Do (RIGHT)
```
Generate random sequences with same overall AA composition → different positional frequencies → meaningful null distribution
```

## Fixed Implementation Coming Next!

The permutation test needs to be redesigned to:
1. Calculate overall AA frequencies from the alignment
2. Generate random sequences matching these frequencies
3. Calculate entropy on random alignments
4. Compare observed vs. null

This will give us a **meaningful** test of whether conservation is real or just due to limited AA diversity.

---

## Visual Summary

```
CURRENT (BROKEN):
Original:  A K R G    Entropy: [0.00, 0.31, 0.18, 0.46]
           A R K T
           A K K S
           A R R N

Shuffle:   A R K S    Entropy: [0.00, 0.31, 0.18, 0.46]  ← SAME!
           A K R N
           A R K G
           A K R T

Result: Null = Observed (always) ❌


FIXED (CORRECT):
Original:  A K R G    Entropy: [0.00, 0.31, 0.18, 0.46]
           A R K T    Mean: 0.24
           A K K S
           A R R N

Random:    R K A S    Entropy: [0.46, 0.31, 0.46, 0.46]  ← DIFFERENT!
           A G R K    Mean: 0.42
           K A R A
           R R K G

After 1000 random alignments:
Null mean: 0.40 ± 0.05
Observed: 0.24
p-value: 0.001 *** (significant conservation!) ✅
```
