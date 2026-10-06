import json
import re
from datetime import datetime
from pathlib import Path

# =========================================================
# Versioned output folders
# =========================================================
#
# Every full run of an analysis writes to a new folder
# <gold>/<prefix>_<N> and records what it is in:
#   * <folder>/run_manifest.json  (description + all parameters)
#   * <gold>/<prefix>_versions.md (one row per version)

REGISTRY_HEADER = (
    '| Version | Date | Description | Results |\n'
    '|---|---|---|---|\n'
)


def _version_numbers(gold_path, prefix):

    pattern = re.compile(rf'^{re.escape(prefix)}_(\d+)$')

    return [
        int(match.group(1))
        for path in Path(gold_path).iterdir()
        if path.is_dir() and (match := pattern.match(path.name))
    ]


def next_version_dir(gold_path, prefix='objective_criteria'):
    """
    Create and return the next free folder <gold_path>/<prefix>_<N+1>.
    """

    numbers = _version_numbers(gold_path, prefix)
    version = max(numbers, default=0) + 1

    out_dir = Path(gold_path) / f'{prefix}_{version}'
    out_dir.mkdir(parents=True, exist_ok=False)

    return out_dir


def _jsonable(value):

    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}

    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]

    if isinstance(value, Path):
        return str(value)

    if hasattr(value, 'item'):   # numpy scalars
        return value.item()

    return value


def write_manifest(out_dir, description, **params):
    """
    Write run_manifest.json (version, date, description and parameters).
    """

    manifest = {
        'version': Path(out_dir).name,
        'date': datetime.now().isoformat(timespec='seconds'),
        'description': description,
        **_jsonable(params),
    }

    with open(Path(out_dir) / 'run_manifest.json', 'w', encoding='utf-8') as f:
        json.dump(manifest, f, indent=2)

    return manifest


def append_registry(gold_path, version, description, results, prefix='objective_criteria'):
    """
    Append one row to <gold_path>/<prefix>_versions.md.
    """

    registry = Path(gold_path) / f'{prefix}_versions.md'

    if not registry.exists():
        registry.write_text(f'# {prefix} versions\n\n' + REGISTRY_HEADER, encoding='utf-8')

    row = f"| {version} | {datetime.now():%Y-%m-%d} | {description} | {results} |\n"

    with open(registry, 'a', encoding='utf-8') as f:
        f.write(row)

    return registry
