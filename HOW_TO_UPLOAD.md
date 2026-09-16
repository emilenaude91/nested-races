# Putting this on GitHub — two ways

## A. In the browser (no git needed, ~10 minutes)
1. github.com → "+" → **New repository**. Name: `nested-races`. Public. Do **not** tick "Add a README" (one is included here). Create.
2. On the empty-repo page click **uploading an existing file**. Drag the *contents* of this folder in (README.md, CITATION.cff, LICENSE, CHANGELOG.md, .gitignore, and the `paper/` and `replication/` folders — drag the folders themselves; GitHub keeps the structure). Commit message: `v1.0 — first public version`. **Commit changes**.
3. Open README.md in the repo and replace `emilenaude91` (two places, plus one in CITATION.cff) with your GitHub username. Commit.
4. **Releases → Create a new release**: tag `v1.0`, title `Version 1.0 (15 September 2026)`, paste the CHANGELOG entry, attach `paper/naude_2026_nested_races_v1.0.pdf`. Publish release. A release is the citable, frozen snapshot; it is also what Zenodo archives when the GitHub–Zenodo link works for you later (one DOI per release, automatically).

## B. With git (if installed)
```bash
cd nested-races
git init && git add . && git commit -m "v1.0 — first public version"
git branch -M main
git remote add origin https://github.com/emilenaude91/nested-races.git
git push -u origin main
git tag v1.0 && git push origin v1.0
```
Then create the release from the tag in the browser as in A.4.

## Afterwards
- The paper's public URL will be `https://github.com/emilenaude91/nested-races/blob/main/paper/naude_2026_nested_races_v1.0.pdf`; use it in SSRN/MPRA abstracts until a DOI exists.
- Google Scholar does not index GitHub; SSRN, MPRA or arXiv still do that job. GitHub is the canonical, versioned home; the others are discovery channels.
- When Zenodo works: Zenodo → GitHub integration → flip the switch for `nested-races` → publish a release → DOI minted automatically and shown on the release page.
