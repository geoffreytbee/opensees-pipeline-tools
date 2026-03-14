# Contributing

Contributions are welcome. This project values engineering rigor and clear documentation over rapid feature development.

## Development Setup

```bash
git clone https://github.com/geoffreybee/opensees-pipeline-tools.git
cd opensees-pipeline-tools
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev,docs]"
```

## Code Standards

- **Formatting:** black (line length 88)
- **Linting:** ruff (rules: E, F, I, W)
- **Docstrings:** Google style on all public functions
- **Type hints:** on all function signatures
- **Units:** every numerical constant must have its unit in a comment
- **Tests:** pytest; add tests for any new functionality

## Workflow

1. Create a feature branch from `main`
2. Write code with tests and documentation
3. Run `black .` and `ruff check .` before committing
4. Run `pytest` to verify all tests pass
5. Open a pull request with a clear description of changes

## Commit Messages

Use conventional commits:

- `feat:` new feature
- `fix:` bug fix
- `docs:` documentation changes
- `test:` adding or modifying tests
- `refactor:` code restructuring without behavior change
- `chore:` build, CI, dependency updates

## Engineering Documentation

When adding new analysis capabilities:

- Document the theory with equation references in `docs/theory/`
- List assumptions in `docs/assumptions.md`
- Add validation benchmarks with quantitative error metrics
- Reference specific sections of ASCE ALA 2001, PRCI, or other standards
