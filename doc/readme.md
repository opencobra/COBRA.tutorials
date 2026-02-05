# Continuous Integration for COBRA Tutorials (For future developers)

This document explains in detail how the CI/CD system works for automatically converting `.mlx` tutorial files into `.html` and `.m` formats, and then syncing them into the `gh-pages` branch of the `cobratoolbox` repository. This guide assumes zero background in CI/CD and is intended for future developers or contributors.

---

## ✨ Objective

Whenever a `.mlx` tutorial file is pushed to the `master` branch of the `COBRA.tutorials` repository, the CI workflows:

1. **Update Tutorial Workflow**: Converts the `.mlx` file into `.html` and `.m` formats using MATLAB.
2. Saves the `.html` file to the `gh-pages` branch of the `cobratoolbox` repo (not COBRA.tutorials).
3. Saves the `.m` file back to the `COBRA.tutorials` repo.
4. Updates the index page to include or update the tutorial.
5. **Remove Deleted Tutorials Workflow**: Detects when `.mlx` files are deleted and removes corresponding `.html` files from the website.

---

## 📄 Workflow 1: Update Tutorial

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
  uses: actions/checkout@v4
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
    echo "REPO_OWNER=${{ github.repository_owner }}" >> "$GITHUB_ENV"
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
  uses: actions/setup-python@v5
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

### 7. **Locate MATLAB**
```yaml
- name: Locate MATLAB
  run: |
    if command -v matlab >/dev/null 2>&1; then
      ML_BIN="$(command -v matlab)"
    else
      ML_BIN="$(ls -1d /usr/local/MATLAB/*/bin/matlab 2>/dev/null | sort -V | tail -n1 || true)"
    fi
    echo "ML_BIN=${ML_BIN}" >> "$GITHUB_ENV"
```

This dynamically locates the MATLAB installation on the runner, ensuring compatibility across different MATLAB versions.

### 8. **Determine Changed `.mlx` Files**
```yaml
last_sync_commit=$(git log --grep="created .mlx and .m files" -n 1 --pretty=format:%H)
```

This finds the most recent commit that pushed `.mlx` and `.m` files. It marks the last known successful sync.

If that commit is found, we compare the current HEAD to that commit:
```bash
changed_files=$(git diff --name-only "$last_sync_commit"..HEAD | grep '\.mlx' | grep -v '^deprecated/' || true)
```

If no such commit exists (first-time run), all `.mlx` files are considered:
```bash
changed_files=$(git ls-files '*.mlx')
```

### 9. **Convert Files**

For each changed `.mlx` file:

* Absolute path is resolved.
* Converted to `.html` and `.m` using MATLAB:
```bash
xvfb-run -a "$ML_BIN" -batch "try, export('$ABSOLUTE_FILE_PATH', '$HTML_FILE_PATH', 'Format', 'html'); catch e, disp(getReport(e,'extended')); exit(1); end"
xvfb-run -a "$ML_BIN" -batch "try, export('$ABSOLUTE_FILE_PATH', '$M_FILE_PATH', 'Format', 'm'); catch e, disp(getReport(e,'extended')); exit(1); end"
```

* HTML output goes directly into `cobratoolbox/stable/tutorials/...`
* `.m` file is saved in the source `COBRA.tutorials` repository

**Note**: PDF generation has been removed to keep the repository lighter.

### 10. **Update Index**

After converting all files, the workflow rebuilds ALL tutorial wrapper pages:
```bash
for raw_html in $(find stable/tutorials -mindepth 2 -type f -name "*.html"); do
  python stable/extract_info.py "$raw_html"
done
```

This Python script:
- Wraps each raw MATLAB-exported HTML file with the `HOLDER_TEMPLATE.html` template
- Updates `index.html` to reflect the new tutorial or updates the existing entry
- Creates download links for `.mlx` and `.m` files on GitHub

### 11. **Push Changes to Both Repos**
```bash
cd cobratoolbox
git add .
git commit -m "Sync tutorial files"
git push [...] gh-pages

cd ..
rm -rf cobratoolbox
git add .
git commit -m "created .mlx and .m files"
git push origin master
```

The first commit updates the `gh-pages` branch (for the website), the second commits `.m` files back to `COBRA.tutorials`.

---

## 📄 Workflow 2: Remove Deleted Tutorials

This workflow handles the case when a `.mlx` tutorial file is deleted from the repository.

### 1. **Workflow Trigger**

Same as Update Tutorial workflow - triggers on any `.mlx` file changes, including deletions.

### 2. **Detect Deleted Files**
```bash
last_sync_commit=$(git log --grep="removed .html and .m files" -n 1 --pretty=format:%H)
```

Looks for the last deletion sync commit. If none found, falls back to looking for creation sync commits.
```bash
deleted_files=$(git diff --name-only --diff-filter=D "$last_sync_commit"..HEAD | grep '\.mlx' || true)
```

The `--diff-filter=D` flag specifically identifies deleted files.

### 3. **Remove Files from Website**

For each deleted `.mlx` file:
```bash
HTML_FILE_PATH="stable/tutorials/$(dirname "$file")/$HTML_FILE_NAME"
git rm "$HTML_FILE_PATH" || rm -f "$HTML_FILE_PATH"
```

Removes the corresponding HTML file from the website.

### 4. **Update Index**
```bash
python stable/remove_from_index.py "$HTML_FILE_PATH"
```

This script removes the tutorial entry from `index.html` and deletes the tutorial wrapper page.

### 5. **Push Changes**
```bash
git add -A
git commit -m "Remove deleted tutorials from source repo"
git push [...] gh-pages
```

Commits and pushes the deletions to the `gh-pages` branch.

**Note**: This workflow does NOT commit anything back to `COBRA.tutorials` since the files are already deleted there.

---

## 🪧 Setup Requirements

To make this CI work, you need:

* A self-hosted GitHub runner with MATLAB and `xvfb` installed.
* The token `DEST_REPO_TOKEN` as a secret in the GitHub repo settings for `COBRA.tutorials`, allowing write access to `cobratoolbox`.
* The following scripts in the `gh-pages` branch of `cobratoolbox`:
  - `stable/extract_info.py` - wraps tutorials and updates index
  - `stable/remove_from_index.py` - removes deleted tutorials from index
  - `stable/tutorials/HOLDER_TEMPLATE.html` - template for wrapping tutorials
  - `stable/tutorials/index.html` - main tutorial index page
* MATLAB must be licensed and accessible on the runner (the workflow auto-detects the MATLAB binary).

---

## ⚡ Recovery If CI Breaks

* **Index not updating?**

  * Check logs, especially the output of `extract_info.py`.
  * Confirm `index.html` and `HOLDER_TEMPLATE.html` exist and are correct.
  * Verify that BeautifulSoup4 is installed properly.

* **Nothing runs?**

  * Check if the `.mlx` file was committed to `master`.
  * Confirm `.mlx` is inside a valid path and not excluded.
  * Look at the GitHub Actions logs for any errors.

* **MATLAB export failing?**

  * Test manually with `xvfb-run matlab -batch ...`.
  * Ensure all file paths are valid and writable.
  * Check MATLAB license status on the runner.
  * Verify MATLAB version supports the export formats.

* **Workflow is skipping files?**

  * Ensure the commit message `created .mlx and .m files` exists in history.
  * Push a `.mlx` change and monitor the CI logs.
  * Check if files are in the `deprecated/` folder (which is excluded).

* **Deleted files not being removed?**

  * Verify that `remove_from_index.py` exists in the `gh-pages` branch.
  * Check that the deletion workflow ran (both workflows trigger on `.mlx` changes).
  * Look for the "removed .html and .m files" commit message.

* **First run after PDF removal processes all files?**

  * This is expected behavior - the workflow won't find the old commit message pattern.
  * After the first run, subsequent runs will only process changed files.
  * This is actually beneficial as it regenerates all tutorials without PDFs.

---

## 🔗 Related Repositories

* [COBRA.tutorials](https://github.com/opencobra/COBRA.tutorials): Source `.mlx` tutorials and generated `.m` files.
* [cobratoolbox](https://github.com/opencobra/cobratoolbox): Published `.html` tutorials on `gh-pages` branch.

---

## 📝 File Format Summary

| Format | Generated? | Location | Purpose |
|--------|-----------|----------|---------|
| `.mlx` | No (source) | COBRA.tutorials/master | Original MATLAB Live Script tutorial |
| `.html` | Yes | cobratoolbox/gh-pages | Web-viewable tutorial for documentation site |
| `.m` | Yes | COBRA.tutorials/master | Plain MATLAB script version |
| `.pdf` | No (removed) | N/A | Previously generated, now removed to reduce repo size |

---

For any further questions, contact the COBRA Toolbox maintainers or check issues in the repositories.

Happy modelling!
