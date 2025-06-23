# Continuous Integration for COBRA Tutorials

This document explains in detail how the CI/CD system works for automatically converting `.mlx` tutorial files into `.html`, `.pdf`, and `.m` formats, and then syncing them into the `gh-pages` branch of the `cobratoolbox` repository. This guide assumes zero background in CI/CD and is intended for future developers or contributors.

---

## ✨ Objective

Whenever a `.mlx` tutorial file is pushed to the `master` branch of the `COBRA.tutorials` repository, this CI workflow:

1. Converts the `.mlx` file into `.html`, `.pdf`, and `.m` formats using MATLAB.
2. Saves the `.html` file to the `gh-pages` branch of the `cobratoolbox` repo (not COBRA.tutorials).
3. Updates the index page to include or update the tutorial.
4. Pushes updates to both repositories appropriately.

---

## 📄 Workflow Breakdown

### 1. **Workflow Trigger**

```yaml
on:
  push:
    branches: [ master ]
    paths:
      - '**.mlx'
```

This means the workflow runs whenever any `.mlx` file is changed or added to the `master` branch.

### 2. **Checkout Source Repo**

```yaml
- name: Checkout Source Repo
  uses: actions/checkout@v2
  with:
    repository: '${{ github.repository_owner }}/COBRA.tutorials'
    token: ${{ secrets.GITHUB_TOKEN }}
    fetch-depth: 0
```

This checks out the `COBRA.tutorials` repo to the self-hosted runner, with full history (needed for Git diff).

### 3. **Get Repository Owner**

```yaml
- name: Get the repository's owner name
  run: |
    echo "REPO_OWNER=${{ github.repository_owner }}" >> $GITHUB_ENV
```

Stores the repo owner into an environment variable. Mostly for flexibility.

### 4. **Clone Destination Repository (gh-pages)**

```yaml
- name: Clone the destination repository
  run: |
    rm -rf cobratoolbox
    git clone --depth 1 --branch gh-pages https://x-access-token:${{ secrets.DEST_REPO_TOKEN }}@github.com/opencobra/cobratoolbox.git
```

This checks out the `gh-pages` branch of the `cobratoolbox` repo. This is where tutorials will be published as `.html`.

### 5. **Setup Python**

```yaml
- name: Set up Python
  uses: actions/setup-python@v2
  with:
    python-version: '3.x'
```

Python is needed to run the `extract_info.py` script which updates the `index.html`.

### 6. **Install Dependencies**

```yaml
- name: Install Dependencies
  run: |
    python -m pip install --upgrade pip
    pip install beautifulsoup4
```

Installs Python dependencies used for parsing and modifying HTML.

### 7. **Determine Changed `.mlx` Files**

```yaml
last_sync_commit=$(git log --grep="created .pdf, .mlx and .m files" -n 1 --pretty=format:%H)
```

This finds the most recent commit that pushed `.pdf`, `.mlx`, and `.m` files. It marks the last known successful sync.

If that commit is found, we compare the current HEAD to that commit:

```bash
changed_files=$(git diff --name-only "$last_sync_commit"..HEAD | grep '\.mlx' || true)
```

If no such commit exists (first-time run), all `.mlx` files are considered.

### 8. **Convert Files and Update Index**

For each changed `.mlx` file:

* Absolute path is resolved.
* Converted to `.html`, `.pdf`, `.m` using MATLAB:

```bash
xvfb-run /usr/local/MATLAB/R2024a/bin/matlab -batch "export('$ABSOLUTE_FILE_PATH', '$OUTPUT_FILE_PATH', 'Format', '<html/pdf/m>')"
```

* HTML output goes directly into `cobratoolbox/stable/tutorials/...`
* Then:

```bash
python stable/extract_info.py "$HTML_RELATIVE_PATH"
```

This Python script updates `index.html` to reflect the new tutorial or updates the entry.

### 9. **Push Changes to Both Repos**

```bash
cd cobratoolbox
# Commit to gh-pages
cd ..
rm -rf cobratoolbox
# Commit to COBRA.tutorials (for .m/.pdf/.mlx)
```

The first commit updates the `gh-pages` (for web), the second commits changes in `COBRA.tutorials` itself (mostly PDF/M files).

---

## 🪧 Setup Requirements

To make this CI work, you need:

* A self-hosted GitHub runner with MATLAB and `xvfb` installed.
* The token `DEST_REPO_TOKEN` as a secret in the GitHub repo settings for `COBRA.tutorials`, allowing write access to `cobratoolbox`.
* The `stable/extract_info.py` script in the `gh-pages` branch of `cobratoolbox`, with `HOLDER_TEMPLATE.html` and `index.html` present and editable.
* MATLAB must be licensed and accessible as `/usr/local/MATLAB/R2024a/bin/matlab` on the runner.

---

## ⚡ Recovery If CI Breaks

* **Index not updating?**

  * Check logs, especially the output of `extract_info.py`.
  * Confirm `index.html` and `HOLDER_TEMPLATE.html` exist and are correct.

* **Nothing runs?**

  * Check if the `.mlx` file was committed to `master`.
  * Confirm `.mlx` is inside a valid path and not excluded.

* **MATLAB export failing?**

  * Test manually with `xvfb-run matlab -batch ...`.
  * Ensure all file paths are valid and writable.

* **Workflow is skipping files?**

  * Ensure the commit message `created .pdf, .mlx and .m files` exists in history.
  * Push a `.mlx` change and monitor the CI.

---

## 🔗 Related Repositories

* [COBRA.tutorials](https://github.com/opencobra/COBRA.tutorials): Source `.mlx` tutorials.
* [cobratoolbox](https://github.com/opencobra/cobratoolbox): Published `.html` tutorials on `gh-pages` branch.

---

For any further questions, contact the COBRA Toolbox maintainers or check issues in the repositories.

Happy modelling!
