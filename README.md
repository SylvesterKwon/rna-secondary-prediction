# RNA Secondary structure prediction
This program calculates the structure of RNA's secondary structure that has maximum matching score and satisfy specific conditions.

Written in C++

## Conditions
1. Watson-Crick pair and wobble pair is both allowed. If you want to set an arbitrary matching score, modify the argument parameters.
1. Hairpin is allowed.
1. Pseudoknot is allowed.

## Restrictions
1. Section of the stem of the two Hairpin constituting the Pseudoknot is not allowed to overlap.
1. Consider only when Hairpin and Pseudoknot are connected in series. This means that the structure surrounding the existing structure (nested structure) is not considered. (If you want this feature, you can modify it)

## Performance
- TIME COMPLEXITY: $O(N^5)$
- SPACE COMPLEXITY: $O(N^4)$

## See also...
- AJOU SOFTCON (poster & videos): https://softcon.ajou.ac.kr/works/works.asp?uid=518&category=M
