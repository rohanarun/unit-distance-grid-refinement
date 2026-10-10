# Tweet drafts (one per PR)

## PR #2 — scaling symmetry + collinear anchors (Oct 10, 2026)

**Main post (ELI25):**

OpenAI's unit-distance proof ends by checking a 20×5 grid of points (100 edges); by noticing the equations don't care if you rescale every velocity by the same number, we got that down to 8×4 = 32, and to 12 when everything sits on one line. Conditional on their lemmas, no new exponent, sympy-checked: github.com/rohanarun/unit-distance-grid-refinement

**Alt text for the chart:**
Step chart of "edges in the final algebraic witness" over time. A black dot marks OpenAI's preprint on Sep 23 at 100 edges (20×5 grid). On Oct 9 the general case drops to 36 and the common-line case to 18. On Oct 10 three branches end at 32 (general, 8×4), 24 (one locus is a line, 6×4) and 12 (both loci the same line, 4×3).

**Thread reply (why it stops here):**
Why not go further? The symmetry group is exactly 5-dimensional (2 slides, 1 spin, 1 shift, 1 scale), so 3j−5 is the floor for normalizing anchor data. Four anchors are forced in the general case: three equations fix three unknowns, the fourth supplies the contradiction. Getting below 32 needs a genuinely new idea, not another normalization.

---

## PR #1 — common-line branch (Oct 9, 2026)

The final step of OpenAI's unit-distance proof needs a grid of points. When both loci are the same line, two equations plus one sign flip are enough: 6×3 = 18 edges instead of 100. Conditional on the same lemmas; no new exponent.
