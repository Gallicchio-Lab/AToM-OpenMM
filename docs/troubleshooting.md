# Troubleshooting

## MkDocs Command Not Found

Create and activate the docs conda environment:

```bash
conda env create -f docs/environment.yml
conda activate atom-openmm-docs
```

## Strict Build Fails

Run the build from the repository root:

```bash
mkdocs build --strict
```

Common issues include broken Markdown links, missing pages listed in `mkdocs.yml`, or invalid YAML indentation.

## Images Do Not Appear

Store image files under `docs/assets/images/` and reference them from Markdown with relative links from the page that uses them.
