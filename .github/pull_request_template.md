## Description

Summarize your changes and the motivation behind them.

Closes #(issue number).

### Type of change (delete unused ones)

- Bug fix
- New feature
- Refactor
- Documentation
- Other (describe)

## Testing

How did you test your changes?

## Checklist

__Check these like this `[x]` to indicate which of the below applies.__

- [ ] I added or updated tests for new behavior
- [ ] I updated documentation if user-facing behavior changed

See the [developer guide](https://mflowcode.github.io/documentation/contributing.html) for full coding standards.

<details>
<summary><strong>GPU changes</strong> (expand if you modified <code>src/simulation/</code>)</summary>

- [ ] GPU results match CPU results
- [ ] Tested on NVIDIA GPU or AMD GPU

</details>

## AI code reviews

Reviews are not retriggered automatically. To request a review, comment on the PR:
- `@claude full review` — Claude full review (also triggers on PR open/reopen/ready)
- Or add label `claude-full-review` — Claude full review via label
