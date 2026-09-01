#!/usr/bin/env python3
"""Convert all rotor--stator numerical data to Tecplot ASCII DAT files.

The converter is non-destructive: source CSV, JSON, and NPZ files remain
unchanged.  The output tree mirrors ``rotor_stator/data`` below
``rotor_stator/tecplot/data``.  Profile and modal tables containing a
``branch`` column are written as one Tecplot ZONE per branch.  Text fields are
preserved as zone AUXDATA because classic Tecplot variables are numeric.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
ROTOR_ROOT = SCRIPT_DIR.parent
DEFAULT_DATA = ROTOR_ROOT / "data"
DEFAULT_REFERENCE = ROTOR_ROOT / "reference" / "baseflow_Res1000.npz"
DEFAULT_OUTPUT = ROTOR_ROOT / "tecplot"


@dataclass
class ConversionResult:
    source: str
    output: str
    source_type: str
    rows: int
    zones: int
    numeric_variables: list[str]
    text_aux_fields: list[str]


def clean_name(value: str) -> str:
    value = value.strip()
    value = re.sub(r"\s+", "_", value)
    value = re.sub(r"[^0-9A-Za-z_]+", "_", value)
    value = value.strip("_")
    if not value:
        return "value"
    if value[0].isdigit():
        value = "v_" + value
    return value


def quote(value: Any) -> str:
    return str(value).replace("\\", "/").replace('"', "\\\"")


def parse_number(value: Any) -> float:
    if value is None:
        return float("nan")
    if isinstance(value, bool):
        return 1.0 if value else 0.0
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value)
    text = str(value).strip()
    if not text:
        return float("nan")
    if text.lower() == "true":
        return 1.0
    if text.lower() == "false":
        return 0.0
    return float(text)


def is_numeric_column(rows: list[dict[str, Any]], column: str) -> bool:
    saw_value = False
    for row in rows:
        value = row.get(column)
        if value is None or str(value).strip() == "":
            continue
        saw_value = True
        try:
            parse_number(value)
        except (TypeError, ValueError):
            return False
    return saw_value


def format_number(value: Any) -> str:
    number = parse_number(value)
    if math.isnan(number):
        return "NaN"
    if math.isinf(number):
        return "Inf" if number > 0.0 else "-Inf"
    return f"{number:.16g}"


def normalized_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            return [], []
        original = list(reader.fieldnames)
        headers = [clean_name(name) for name in original]
        if len(set(headers)) != len(headers):
            raise ValueError(f"Duplicate normalized columns in {path}")
        rows = []
        for raw in reader:
            row = {}
            for old, new in zip(original, headers):
                value = raw.get(old, "")
                row[new] = value.strip() if isinstance(value, str) else value
            rows.append(row)
    return headers, rows


def flatten_json(path: Path) -> tuple[list[str], list[dict[str, Any]]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(payload, list):
        if not payload:
            return ["record_count"], [{"record_count": 0}]
        rows = [dict(item) for item in payload]
        headers = list(dict.fromkeys(key for row in rows for key in row))
        return [clean_name(key) for key in headers], [
            {clean_name(key): row.get(key, "") for key in headers} for row in rows
        ]

    if not isinstance(payload, dict):
        return ["value"], [{"value": payload}]

    # temporal_convergence.json: branch -> collocation order -> metrics.
    if payload and all(isinstance(value, dict) for value in payload.values()):
        nested_rows: list[dict[str, Any]] = []
        nested = True
        for outer_key, outer_value in payload.items():
            if not outer_value or not all(
                isinstance(value, dict) for value in outer_value.values()
            ):
                nested = False
                break
            for inner_key, metrics in outer_value.items():
                row = {
                    "branch": outer_key,
                    "collocation_points": inner_key,
                    **metrics,
                }
                nested_rows.append(row)
        if nested:
            headers = list(
                dict.fromkeys(key for row in nested_rows for key in row)
            )
            return headers, nested_rows

    headers = [clean_name(key) for key in payload]
    return headers, [{clean_name(key): value for key, value in payload.items()}]


def zone_groups(
    rows: list[dict[str, Any]],
    numeric_columns: list[str],
    text_columns: list[str],
) -> list[tuple[str, list[dict[str, Any]], dict[str, str]]]:
    if not rows:
        return [("empty", [{"record_count": 0}], {})]

    group_columns = list(text_columns)
    if "branch" in numeric_columns:
        group_columns.append("branch")

    if not group_columns:
        return [("data", rows, {})]

    grouped: dict[tuple[str, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        key = tuple(str(row.get(column, "")) for column in group_columns)
        grouped[key].append(row)

    zones = []
    for index, (key, members) in enumerate(grouped.items(), start=1):
        metadata = {
            clean_name(column): value
            for column, value in zip(group_columns, key)
            if column in text_columns
        }
        parts = []
        for column, value in zip(group_columns, key):
            if value:
                parts.append(f"{clean_name(column)}={value}")
        title = ", ".join(parts) if parts else f"zone_{index}"
        zones.append((title[:180], members, metadata))
    return zones


def write_tecplot(
    output: Path,
    title: str,
    source_relative: str,
    headers: list[str],
    rows: list[dict[str, Any]],
) -> tuple[int, list[str], list[str]]:
    numeric_columns = [
        column for column in headers if is_numeric_column(rows, column)
    ]
    text_columns = [column for column in headers if column not in numeric_columns]
    if not numeric_columns:
        numeric_columns = ["record_count"]
        rows = [{**row, "record_count": 1} for row in rows] or [
            {"record_count": 0}
        ]

    zones = zone_groups(rows, numeric_columns, text_columns)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(f'TITLE = "{quote(title)}"\n')
        stream.write(
            "VARIABLES = "
            + " ".join(f'"{quote(column)}"' for column in numeric_columns)
            + "\n"
        )
        stream.write(f'DATASETAUXDATA Source_File = "{quote(source_relative)}"\n')
        stream.write('DATASETAUXDATA Converter = "convert_all_to_tecplot.py"\n')
        for zone_title, members, metadata in zones:
            stream.write(
                f'ZONE T="{quote(zone_title)}", I={len(members)}, '
                "DATAPACKING=POINT\n"
            )
            for key, value in metadata.items():
                stream.write(f'AUXDATA {key} = "{quote(value)}"\n')
            for row in members:
                stream.write(
                    " ".join(format_number(row.get(column)) for column in numeric_columns)
                    + "\n"
                )
    return len(zones), numeric_columns, text_columns


def convert_table(
    source: Path,
    source_base: Path,
    output: Path,
) -> ConversionResult:
    if source.suffix.lower() == ".csv":
        headers, rows = normalized_csv(source)
        source_type = "csv"
    elif source.suffix.lower() == ".json":
        headers, rows = flatten_json(source)
        source_type = "json"
    else:
        raise ValueError(f"Unsupported tabular source: {source}")

    relative = source.relative_to(source_base).as_posix()
    zones, variables, text_fields = write_tecplot(
        output,
        title=source.stem,
        source_relative=relative,
        headers=headers,
        rows=rows,
    )
    return ConversionResult(
        source=relative,
        output=output.as_posix(),
        source_type=source_type,
        rows=len(rows),
        zones=zones,
        numeric_variables=variables,
        text_aux_fields=text_fields,
    )


def convert_npz(
    source: Path,
    source_base: Path,
    output: Path,
) -> ConversionResult:
    with np.load(source) as archive:
        keys = list(archive.files)
        if "x" in keys:
            keys = ["x", *[key for key in keys if key != "x"]]
        arrays = {key: np.asarray(archive[key]).reshape(-1) for key in keys}
    lengths = {array.size for array in arrays.values()}
    if len(lengths) != 1:
        raise ValueError(f"NPZ arrays have inconsistent lengths: {source}")
    row_count = lengths.pop()
    rows = [
        {key: arrays[key][index] for key in keys}
        for index in range(row_count)
    ]
    relative = source.relative_to(source_base).as_posix()
    zones, variables, text_fields = write_tecplot(
        output,
        title=source.stem,
        source_relative=relative,
        headers=keys,
        rows=rows,
    )
    return ConversionResult(
        source=relative,
        output=output.as_posix(),
        source_type="npz",
        rows=row_count,
        zones=zones,
        numeric_variables=variables,
        text_aux_fields=text_fields,
    )


def relative_output(output_root: Path, path: Path) -> str:
    return path.relative_to(output_root).as_posix()


def write_index(output_root: Path, results: list[ConversionResult]) -> None:
    manifest = []
    for result in results:
        manifest.append(
            {
                **result.__dict__,
                "output": relative_output(output_root, Path(result.output)),
            }
        )
    (output_root / "conversion_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    lines = [
        "# Rotor--stator Tecplot data",
        "",
        "本目录由 `convert_all_to_tecplot.py` 从原始 CSV、JSON 和 NPZ 数据只读生成。",
        "所有文件采用 Tecplot ASCII、`DATAPACKING=POINT`。包含 `branch` 的剖面或谱数据按 branch 分为多个 `ZONE`；文本字段写入 zone `AUXDATA`。",
        "",
        "| DAT 文件 | 来源 | 行数 | ZONE 数 | 数值变量数 |",
        "|---|---|---:|---:|---:|",
    ]
    for result in results:
        output = relative_output(output_root, Path(result.output))
        lines.append(
            f"| `{output}` | `{result.source}` | {result.rows} | "
            f"{result.zones} | {len(result.numeric_variables)} |"
        )
    lines.extend(
        [
            "",
            "## 直接计算的基本流剖面",
            "",
            "`data/baseflow_profiles/principal_branch_Tw1.00_1.16_step0.04.dat` 由 `rotor_stator/scripts/generate_principal_baseflow_tecplot.py` 直接计算生成（若已运行该脚本），包含等温相连主支在 $T_w=1.00,1.04,1.08,1.12,1.16$ 的五个 2001 点 ZONE。",
            "",
            "## Tecplot 导入",
            "",
            "在 Tecplot 360 中选择 `File > Load Data File(s)`，数据类型选择 `Tecplot Data Loader`，再选择一个或多个 `.dat` 文件。",
            "缺失数值写成 `NaN`。空的折返点 JSON 会转换成一个 `record_count=0` 的占位 zone。",
            "",
        ]
    )
    (output_root / "README.md").write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--reference-npz", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data_root = args.data_root.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    results: list[ConversionResult] = []

    sources: Iterable[Path] = sorted(
        [*data_root.rglob("*.csv"), *data_root.rglob("*.json")]
    )
    for source in sources:
        destination = (
            output_root / "data" / source.relative_to(data_root)
        ).with_suffix(".dat")
        results.append(convert_table(source, data_root, destination))

    reference = args.reference_npz.resolve()
    reference_output = output_root / "reference" / f"{reference.stem}.dat"
    results.append(
        convert_npz(reference, reference.parent, reference_output)
    )
    write_index(output_root, results)

    total_rows = sum(result.rows for result in results)
    total_zones = sum(result.zones for result in results)
    print(
        f"Converted {len(results)} files to {output_root}: "
        f"{total_rows} rows, {total_zones} zones"
    )


if __name__ == "__main__":
    main()
