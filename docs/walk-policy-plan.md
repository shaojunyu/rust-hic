# Planned walk policy support

Status: deferred. This document records future work; the additional policies
below are not implemented. The default BWA-MEM3, BWA-MEM2, and BWA-MEM paths
currently use the built-in `5unique`-style parser. Both local and Slurm runs
use this parser. Chromap produces pairs directly and does not use it.

The goal is to support the walk policies exposed by runHiC while retaining
`5unique` as the default. The reference for future comparisons is
[pairtools 1.1.3 parse](https://github.com/open2c/pairtools/blob/v1.1.3/pairtools/lib/parse.py),
followed by runHiC's existing pair-type selection. Compatibility should be
described precisely rather than assumed from the policy name alone.

| Policy | Intended behavior for walks | Status |
| --- | --- | --- |
| `5unique` | Select the 5′-most unique alignment on each side | Implemented, with documented edge-case differences |
| `3unique` | Select the 3′-most unique alignment on each side | Planned |
| `5any` / `3any` | Select the 5′-most / 3′-most alignment, regardless of uniqueness | Planned |
| `mask` | Mask unrescuable walks | Planned |
| `all` | Report contacts from complex walks, accounting for overlaps within a molecule | Planned |

For the single-output policies, attempt single-ligation rescue before applying
the policy to an unrescuable walk. Continue to retain only `UU`, `UR`, and
`RU` pairs in the BWA pipeline, matching runHiC's current filtering for
pairtools 1.1.3. Consequently, masked `WW` pairs and non-unique pairs chosen
by an `any` policy will not appear in the final selected pairs. Full unfiltered
pairtools output would be a separate feature.

The `all` policy is more involved: one input read pair can yield multiple
contacts, and pairtools handles overlaps between the two read ends within a
molecule. It does not simply emit all combinations of alignments. The current
parser returns at most one pair and stores only the reported 5′ coordinate.
Supporting `all` requires retaining the additional alignment coordinates,
implementing the reference overlap and contact-orientation rules, emitting
multiple selected pairs, and checking downstream reporting. In particular,
selected contact counts may exceed the number of input read pairs, so the
report must not present that ratio as a fraction of reads. Input read counts
and `expected_pairs` must continue to count FASTQ pairs.

When this work is resumed:

- Expose one validated policy in both `map` and `fastq.sam` configuration;
  keep the default and existing output unchanged when the option is omitted.
- Apply policies to all BWA paths, including structured BWA-MEM3 alignments,
  external BWA SAM input, local runs, and Slurm runs. Reject a non-default
  policy for Chromap rather than silently ignoring it.
- Compare each policy with pairtools 1.1.3 followed by the same pair-type
  filter, using targeted and randomized inputs. Cover rescue, mapped and
  non-unique segments, gaps, both strands, and overlapping read ends; for
  `all`, also cover multiple contacts per read pair and output orientation.
- Preserve the previously accepted fixes for leading clipping, internal gap
  detection, and rescue when the first chimeric alignment is unusable. Record
  intentional differences from pairtools when checking compatibility.
- Verify report counts, local and Slurm workflows, optional BAM retention,
  configuration examples, and unchanged default results.
