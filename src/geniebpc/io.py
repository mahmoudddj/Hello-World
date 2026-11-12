"""Loading of GENIE BPC tables from a secure, environment-configured directory. Never writes patient-level data."""
from __future__ import annotations
import os
from pathlib import Path
import pandas as pd
from . import config as C


def data_dir(path: str | None = None) -> Path:
    p = path or os.environ.get(C.DATA_ENV)
    if not p:
        raise SystemExit(f"No data directory: pass --data-dir or set {C.DATA_ENV}. "
                         "GENIE BPC data are access-controlled and are not distributed with this repository.")
    p = Path(p)
    if not p.is_dir():
        raise SystemExit(f"Data directory does not exist: {p}")
    return p


def _find(d: Path, globs: list[str]) -> Path | None:
    for g in globs:
        hits = sorted(x for x in d.rglob(g) if "non_index" not in x.name)
        if hits:
            return hits[0]
    return None


def _read(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix == ".txt" else ","
    return pd.read_csv(path, sep=sep, comment="#", low_memory=False)


def load_tables(path: str | None = None) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """Returns (tables, schema_report). schema_report has one row per required column with found/missing."""
    d = data_dir(path)
    tables, rows = {}, []
    for name, globs in C.TABLE_GLOBS.items():
        f = _find(d, globs)
        if f is None:
            rows.append(dict(table=name, file="", column="<table>", present=False, optional=name in C.OPTIONAL_TABLES))
            continue
        tables[name] = _read(f)
        rows.append(dict(table=name, file=f.name, column="<table>", present=True, optional=name in C.OPTIONAL_TABLES))
        for col in C.REQUIRED_COLUMNS.get(name, []):
            rows.append(dict(table=name, file=f.name, column=col, present=col in tables[name].columns, optional=False))
    return tables, pd.DataFrame(rows)


def schema_ok(report: pd.DataFrame) -> bool:
    return bool(report.loc[~report.optional, "present"].all())
