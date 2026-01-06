Run example
python build_tree_from_sites_dir.py \
  --sites_dir=sites_for_tree \
  --pattern="*.sites" \
  --outprefix=ALL_SITES \
  --n_jobs=8 \
  --recursive



## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
