# Quickstart

After creating and activating the docs environment, preview the site locally from the repository root:

```bash
mkdocs serve
```

MkDocs will print a local URL, usually `http://127.0.0.1:8000/`. Open that URL in a browser and edit Markdown files under `docs/`; the preview rebuilds automatically as files change.

For a production-style check, run:

```bash
mkdocs build --strict
```

The generated site is written to `site/`. That directory is ignored by git and should not be committed.
