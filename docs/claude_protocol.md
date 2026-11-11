# Claude PI Protocol

## Session setup

- URL: https://claude.ai/new
- Model: Claude Opus 5.5, effort **High** (fall back to the highest available if the
  UI changes).
- Language: English for research discussion.

## Prompt contract

Every round must contain, in order:

1. One-paragraph state delta since the PI's last decision (new evidence only).
2. Concrete artifacts: metrics tables, failure modes, code/commit links.
3. The specific decision requested.
4. Constraints (compute, time, data).

Do **not** re-explain the whole project each round.

## Response capture

- Save the full PI response to `docs/claude/transcripts/YYYY-MM-DD_HHMM_<topic>.md`.
- Extract the decision into `docs/decision_log.md`.
- Extract actionable items into `state/next_actions.md`.

## Rate limits and interruption

1. On limit, record reset time shown by the UI in `state/project_state.json`.
2. Continue queued local/server work that does not need the PI.
3. Resume only after the displayed reset time; never bypass limits.

