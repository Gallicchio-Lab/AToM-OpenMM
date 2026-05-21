# Contributing

Contributions to AToM-OpenMM are welcome. Useful contributions include bug reports, fixes, new examples, clearer documentation, tests, and workflow improvements for ABFE and RBFE calculations.

## Before You Start

- Check the [issue tracker](https://github.com/Gallicchio-Lab/AToM-OpenMM/issues) for related reports or ongoing work.
- Open an issue for larger changes so the approach can be discussed before implementation.
- Keep changes focused and describe the scientific or workflow motivation clearly.

## Development Setup

Use an isolated conda or mamba environment with the runtime dependencies described in [Installation](installation.md), then install the repository from source:

```bash
git clone https://github.com/Gallicchio-Lab/AToM-OpenMM.git
cd AToM-OpenMM
pip install .
```

For changes to package code, run the relevant tests or a small representative example before opening a pull request:

```bash
pytest
```

Some tests and examples may require OpenMM-compatible hardware and the full runtime environment.

## Pull Requests

- Create a feature branch from the current upstream branch.
- Include a concise summary of the change and any tests or examples you ran.
- Avoid mixing unrelated refactors with behavior changes.
- Update examples or user-facing documentation when behavior, inputs, or outputs change.

## Research Software Note

AToM-OpenMM is research software under active development. Please include enough detail in issues and pull requests for maintainers to reproduce and review the scientific workflow affected by the change.
