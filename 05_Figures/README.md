# 05_figures proposal

This is a minimal, cleaner mirror of the current `04_Figures` workflow.

## Design intent
- Preserve the current statistical logic and visual identity.
- Centralize styling so palettes, font sizes, label sizes, and axis title formatting stay uniform.
- Keep each panel script self-contained at the script level:
  - it loads the figure setup if needed
  - it loads the shared style helper
  - it returns a named plot object for the composite
- Minimize scientific churn by reusing the existing `04_Figures` panel logic where possible.

## What this package is
This is a **refactor-ready bridge**:
- **F1** is intentionally very close to the current implementation.
- **F2** and **F3** are restructured to match your requested layout more directly.
- The complex panel internals are kept minimal by wrapping current production scripts rather than rewriting the biology/statistics.

## Important note on F2/F3 panel E
Your requested middle-right panel is "RRHO2 + ORA bars". In the current `04_Figures`, the nearest stable production element is the existing pathway summary panel. To keep this refactor minimal and safe, the provided implementation uses:
- the existing RRHO2 map, and
- the existing pathway summary plot as a bridge panel.

That means the **layout contract is correct now**, while the **exact ORA-bars rendering can be swapped in later** without changing the composite structure.
