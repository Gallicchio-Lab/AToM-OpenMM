# Contributing

Documentation changes use the same pull request workflow as code changes.

## Docs Workflow

1. Create the docs environment:

    ```bash
    conda env create -f docs/environment.yml
    ```

2. Activate it:

    ```bash
    conda activate atom-openmm-docs
    ```

3. Preview the site:

    ```bash
    mkdocs serve
    ```

4. Edit Markdown files under `docs/`.

5. Store images under `docs/assets/images/` and reference them with relative Markdown links.

6. Before opening a pull request, run:

    ```bash
    mkdocs build --strict
    ```

7. Open a pull request. The docs workflow will check the strict MkDocs build automatically.
