# MEMORY.md

Session-by-session log of decisions and progress. Read at the start of every session.

---

## 2026-08-28 — .claude/ knowledge base set up

Created `.claude/CLAUDE.md` (working rules, architecture, key params), modelled on the more
extensive setup in `scQC-flow/.claude/` but kept minimal for a smaller pipeline. Carried over:
the core rules, the memory/decision-log convention, and the "who Alex is" section. Left out
`README.md`, `PIPELINE_DESIGN.md`, `ERRORS.md`, `FEATURE_IDEAS.md` for now — add them if/when the
project grows enough to need them.

Decided: no Claude/Anthropic co-authorship on commits in this repo (Co-Authored-By or
Claude-Session trailers) — Alex is sole author on every commit, logged in `CLAUDE.md` under
"Commits".
