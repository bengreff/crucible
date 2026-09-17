# Pre-pivot archive

`pre-pivot-2026-09-17/` preserves the repository contents from before the human-facing sandbox reset. Its documents and agent instructions are historical and have no authority over the active project.

The source revision was `67b1f237499a9c2094dd2ca23eb8a0922a2d80c3` on branch `s13-gpu-residency`. All 236 tracked files and four previously untracked planning documents were preserved byte-for-byte with their file modes. The former root `.github/` and `.claude/` were moved too, so their workflows and local settings are not active root configuration. The root Git history itself remains in place.

[The manifest](pre-pivot-2026-09-17/ARCHIVE_MANIFEST.json) records each of those 240 files' original path, SHA-256 digest, permissions, and tracked status. Paths are relative to the archived root. It also lists the moved top-level entries. The archived source tree itself is the source snapshot; the archive README and manifest document the move.

Ignored build caches, generated runs, local settings, and virtual environments were also moved locally, without deletion. They remain ignored and are **not included in the commit**. A Git checkout alone will not recover those local artifacts. Data or runs that become needed evidence should be deliberately curated, with provenance, rather than added wholesale. Local settings may contain machine-specific or private information and should not be published.

The old environment was not revalidated during this documentation/archive operation. Virtual environments and cached binaries can contain absolute paths; recreate the environment before attempting execution from its new location. Historical build instructions, machine dependencies, and external GPU references may require adjustment. The prior Git revision remains available if an original-layout checkout is needed.

Archive integrity verifies preservation, not physical correctness. Existing certificates and session reports retain their original qualifications and dates. See [rebuild notes](../REBUILD_NOTES.md) before interpreting or reusing them.

To check the archived snapshot without running the old simulator, from the repository root:

```sh
python3 - <<'PY'
import hashlib, json, os, stat
from pathlib import Path
root = Path('archive/pre-pivot-2026-09-17')
manifest = json.loads((root / 'ARCHIVE_MANIFEST.json').read_text())
for item in manifest['files']:
    path = root / item['path']
    data = os.readlink(path).encode() if path.is_symlink() else path.read_bytes()
    assert hashlib.sha256(data).hexdigest() == item['sha256'], item['path']
    assert stat.S_IMODE(path.lstat().st_mode) == item['mode'], item['path']
print(f"Verified {len(manifest['files'])} archived files")
PY
```
