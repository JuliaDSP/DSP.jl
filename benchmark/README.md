# Remez rewrite validation

Run from a checkout containing the original revision:

```sh
julia --startup-file=no benchmark/remez.jl
julia --startup-file=no benchmark/remez_compatibility.jl
julia --project -e 'using Pkg; Pkg.test(; julia_args=["--check-bounds=yes"])'
```

The first two scripts use only standard libraries. Their optional first argument
selects a baseline git revision; the default is
`7c798756cca39251da14bb5d38146feb2cdd1717`. They load the original source and the
working-tree source in separate modules in the same Julia process. They neither
change the checkout nor install benchmark dependencies.

The timing script checks exact coefficient equality for each case, warms up both
implementations, and reports the median of 31 samples of five calls. Sample order
alternates between old/new and new/old. Allocated bytes include the returned
coefficient vector and are measured after compilation. Compilation time and
process resident memory are not measured. Timing results depend on the CPU and
Julia version; they are not assertions in the package tests.

The differential script compares coefficients bit for bit, warning messages, and
exception types/messages on 600 seeded random designs. It varies tap count,
symmetry, number of bands, weights, grid density, partial/full frequency coverage,
and iteration limits. Both scripts require the baseline revision in git history;
the ordinary package tests do not.

## Measured results

Julia 1.12.7, `znver3`, one Julia thread. Microseconds per call:

| Design | Original | Rewrite | Original bytes | Rewrite bytes |
| --- | ---: | ---: | ---: | ---: |
| 15-tap low-pass | 13.76 | 11.45 | 5,152 | 3,776 |
| 151-tap low-pass | 863.71 | 620.49 | 45,656 | 34,328 |
| 152-tap weighted low-pass | 1,297.99 | 922.68 | 45,592 | 34,328 |
| 51-tap high-pass | 59.44 | 42.47 | 16,080 | 11,976 |
| 180-tap band-pass | 2,430.06 | 1,771.71 | 53,064 | 39,848 |
| 20-tap Hilbert | 10.71 | 9.20 | 6,304 | 4,736 |
| 21-tap Hilbert | 9.35 | 8.53 | 6,304 | 4,736 |
| 200-tap differentiator | 954.59 | 715.02 | 59,808 | 45,272 |
| 201-tap differentiator | 1,001.03 | 770.01 | 55,712 | 42,200 |
| 63-tap partial-band design | 78.45 | 56.03 | 15,432 | 11,400 |
| 63-tap low-pass, density 64 | 270.58 | 186.17 | 65,512 | 49,032 |
| 511-tap low-pass | 5,852.15 | 4,957.51 | 151,888 | 114,768 |
| 201-tap inverse sinc | 1,240.34 | 929.87 | 59,024 | 44,808 |
| Three-argument weighted low-pass | 1,311.01 | 934.69 | 45,944 | 34,680 |
| Three-argument Hilbert, 20 taps | 11.42 | 9.76 | 6,672 | 5,104 |
| Three-argument Hilbert, 21 taps | 9.75 | 8.51 | 6,672 | 5,104 |
| Three-argument differentiator, 200 taps | 964.16 | 729.10 | 60,400 | 45,816 |
| Three-argument differentiator, 201 taps | 1,041.57 | 782.92 | 56,304 | 42,744 |

These cases take approximately 9–31% less time and allocate 24–27% fewer bytes.
The rewrite fills the frequency grid directly and reuses interpolation workspace
for the partial-band basis conversion. Its barycentric evaluator groups four
independent divisions while preserving accumulation order with fused operations.
The exchange step keeps the original search order, endpoint margin, and stopping
rules, including the returned approximation when `maxiter` is reached.
