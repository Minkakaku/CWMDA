from __future__ import annotations

import argparse
import itertools
import re
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

pd = None


def require_pandas():
    global pd
    if pd is None:
        try:
            import pandas as pandas_module
        except ModuleNotFoundError as exc:
            raise ModuleNotFoundError(
                "sample_discovery.py requires pandas for discover/suggest/assign commands."
            ) from exc
        pd = pandas_module
    return pd


OUTS_DIR_NAMES = {"outs"}
MATRIX_DIR_NAMES = {"filtered_feature_bc_matrix"}
MATRIX_FILE_NAMES = {"filtered_feature_bc_matrix.h5"}
MANIFEST_COLUMNS = ["sample", "group", "sample_dir", "outs_dir", "matrix_dir", "matrix_type"]


def _as_path(path_value: str) -> Path:
    return Path(path_value).expanduser().resolve()


def _write_table(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)


def read_table(path_value: str) -> pd.DataFrame:
    pd = require_pandas()
    path = _as_path(path_value)
    suffix = path.suffix.lower()
    sep = "\t" if suffix in {".tsv", ".txt", ".xls"} else ","
    return pd.read_csv(path, sep=sep)


def discover_samples(data_root: str) -> pd.DataFrame:
    pd = require_pandas()
    root = _as_path(data_root)
    if not root.exists():
        raise FileNotFoundError("Data root not found: {}".format(root))

    records = []
    for path in sorted(root.rglob("*")):
        if not (path.is_dir() or path.is_file()):
            continue

        path_name = path.name.lower()
        if path_name not in MATRIX_DIR_NAMES and path_name not in MATRIX_FILE_NAMES:
            continue

        outs_dir = path.parent
        if outs_dir.name.lower() not in OUTS_DIR_NAMES:
            continue

        sample_dir = outs_dir.parent
        sample_name = sample_dir.name
        records.append(
            {
                "sample": sample_name,
                "group": sample_name,
                "sample_dir": str(sample_dir),
                "outs_dir": str(outs_dir),
                "matrix_dir": str(path),
                "matrix_type": "h5" if path.is_file() else "mtx_dir",
            }
        )

    if not records:
        outs_like_dirs = []
        for candidate in sorted(root.rglob("*")):
            if candidate.is_dir() and "out" in candidate.name.lower():
                outs_like_dirs.append(str(candidate))
            if len(outs_like_dirs) >= 10:
                break

        hint = ""
        if outs_like_dirs:
            hint = " Nearby outs-like directories: {}".format("; ".join(outs_like_dirs))
        raise ValueError(
            "No sample/outs/filtered_feature_bc_matrix or filtered_feature_bc_matrix.h5 paths were found under {}.{}".format(root, hint)
        )

    manifest = pd.DataFrame(records).drop_duplicates(subset=["matrix_dir"]).sort_values("sample").reset_index(drop=True)

    duplicate_samples = manifest["sample"][manifest["sample"].duplicated()].unique()
    if len(duplicate_samples) > 0:
        duplicate_text = ", ".join(map(str, duplicate_samples))
        raise ValueError(
            "Duplicated sample names were discovered automatically: {}. Rename the sample folders or maintain sample_sheet manually.".format(duplicate_text)
        )

    return manifest[MANIFEST_COLUMNS]


def tokenize_sample_name(sample_name: str) -> List[str]:
    normalized = sample_name.strip()
    normalized = re.sub(r"([a-z])([A-Z])", r"\1_\2", normalized)
    normalized = re.sub(r"([A-Za-z])(\d)", r"\1_\2", normalized)
    normalized = re.sub(r"(\d)([A-Za-z])", r"\1_\2", normalized)
    tokens = [token for token in re.split(r"[_\-\.\s]+", normalized) if token]
    return tokens or [sample_name]


def _format_positions(positions: Sequence[int]) -> str:
    return ",".join(str(position) for position in positions)


def _group_sizes(mapping: Mapping[str, List[str]]) -> str:
    return ",".join(str(len(samples)) for samples in sorted(mapping.values(), key=lambda items: (-len(items), items[0])))


def enumerate_group_suggestions(sample_names: Sequence[str]) -> pd.DataFrame:
    pd = require_pandas()
    sample_names = list(sample_names)
    if len(sample_names) < 2:
        return pd.DataFrame(
            columns=[
                "suggestion_id",
                "positions",
                "group_count",
                "redundancy_saved",
                "balance_score",
                "groups",
                "group_sizes",
            ]
        )

    tokenized = {sample_name: tokenize_sample_name(sample_name) for sample_name in sample_names}
    max_token_count = max(len(tokens) for tokens in tokenized.values())
    suggestions = []
    seen_groupings = set()
    suggestion_index = 1

    searchable_token_count = min(max_token_count, 8)
    for position_count in range(1, searchable_token_count + 1):
        for positions in itertools.combinations(range(max_token_count), position_count):
            if any(any(position >= len(tokens) for position in positions) for tokens in tokenized.values()):
                continue

            group_map: Dict[str, List[str]] = {}
            for sample_name, tokens in tokenized.items():
                group_name = "_".join(tokens[position] for position in positions)
                group_map.setdefault(group_name, []).append(sample_name)

            group_count = len(group_map)
            if group_count <= 1 or group_count >= len(sample_names):
                continue

            grouping_signature = tuple(sorted((group, tuple(sorted(samples))) for group, samples in group_map.items()))
            if grouping_signature in seen_groupings:
                continue
            seen_groupings.add(grouping_signature)

            group_sizes = [len(samples) for samples in group_map.values()]
            balance_score = max(group_sizes) - min(group_sizes)
            suggestions.append(
                {
                    "suggestion_id": "S{:02d}".format(suggestion_index),
                    "positions": _format_positions(positions),
                    "group_count": group_count,
                    "redundancy_saved": len(sample_names) - group_count,
                    "balance_score": balance_score,
                    "groups": "; ".join(
                        "{}:{}".format(group_name, ",".join(sorted(samples)))
                        for group_name, samples in sorted(group_map.items())
                    ),
                    "group_sizes": _group_sizes(group_map),
                }
            )
            suggestion_index += 1

    suggestion_df = pd.DataFrame(suggestions)
    if suggestion_df.empty:
        return pd.DataFrame(
            columns=[
                "suggestion_id",
                "positions",
                "group_count",
                "redundancy_saved",
                "balance_score",
                "groups",
                "group_sizes",
            ]
        )

    return suggestion_df.sort_values(
        ["group_count", "balance_score", "redundancy_saved", "positions", "suggestion_id"],
        ascending=[True, True, False, True, True],
    ).reset_index(drop=True)


def choose_group_assignment(
    manifest: pd.DataFrame,
    group_count: int,
) -> Tuple[pd.DataFrame, Optional[pd.Series], pd.DataFrame]:
    suggestions = enumerate_group_suggestions(manifest["sample"].tolist())
    sample_count = len(manifest)

    if group_count < 1 or group_count > sample_count:
        raise ValueError("group_count must be between 1 and {}.".format(sample_count))

    sample_sheet = manifest.copy()
    if group_count == 1:
        sample_sheet["group"] = "group1"
        return sample_sheet, None, suggestions

    if group_count == sample_count:
        sample_sheet["group"] = sample_sheet["sample"]
        return sample_sheet, None, suggestions

    matched = suggestions[suggestions["group_count"] == group_count].copy()
    if matched.empty:
        available = sorted(set(suggestions["group_count"].tolist()))
        available_text = ", ".join(map(str, available)) if available else "no available suggestions"
        raise ValueError(
            "No exact auto-grouping rule can collapse samples into {} groups. Available group counts: {}.".format(group_count, available_text)
        )

    matched = matched.sort_values(
        ["balance_score", "redundancy_saved", "positions", "suggestion_id"],
        ascending=[True, False, True, True],
    )
    chosen = matched.iloc[0]
    positions = [int(position) for position in str(chosen["positions"]).split(",") if position != ""]

    sample_sheet["group"] = sample_sheet["sample"].map(
        lambda sample_name: "_".join(tokenize_sample_name(sample_name)[position] for position in positions)
    )
    return sample_sheet, chosen, suggestions


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Discover 10x sample directories and generate exact sample-name grouping suggestions."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    discover_parser = subparsers.add_parser("discover", help="Scan the data root and build a sample manifest.")
    discover_parser.add_argument("--data-root", required=True, help="Sample data root.")
    discover_parser.add_argument(
        "--output-manifest",
        required=True,
        help="Output manifest TSV path.",
    )

    suggest_parser = subparsers.add_parser("suggest", help="Generate auto-grouping suggestions from sample names.")
    suggest_parser.add_argument("--manifest", required=True, help="Manifest generated by discover.")
    suggest_parser.add_argument(
        "--output-report",
        required=True,
        help="Output suggestion report TSV path.",
    )

    assign_parser = subparsers.add_parser("assign", help="Generate sample_sheet automatically for a target group count.")
    assign_parser.add_argument("--manifest", required=True, help="Manifest generated by discover.")
    assign_parser.add_argument("--group-count", required=True, type=int, help="Target number of groups.")
    assign_parser.add_argument(
        "--output-sample-sheet",
        required=True,
        help="Output sample_sheet TSV path.",
    )
    assign_parser.add_argument(
        "--output-report",
        default=None,
        help="Optional: also write the full suggestion report TSV.",
    )

    return parser


def main() -> None:
    parser = build_argument_parser()
    args = parser.parse_args()

    if args.command == "discover":
        manifest = discover_samples(args.data_root)
        output_manifest = _as_path(args.output_manifest)
        _write_table(manifest, output_manifest)
        print("Discovered {} samples. Manifest written to: {}".format(len(manifest), output_manifest))
        return

    if args.command == "suggest":
        manifest = read_table(args.manifest)
        suggestions = enumerate_group_suggestions(manifest["sample"].tolist())
        output_report = _as_path(args.output_report)
        _write_table(suggestions, output_report)
        if suggestions.empty:
            print("No auto-grouping rules were found. Empty report written to: {}".format(output_report))
        else:
            print("Generated {} grouping suggestions. Report written to: {}".format(len(suggestions), output_report))
            print("Available group counts: {}".format(", ".join(map(str, sorted(set(suggestions["group_count"].tolist()))))))
        return

    if args.command == "assign":
        manifest = read_table(args.manifest)
        sample_sheet, chosen, suggestions = choose_group_assignment(manifest, args.group_count)
        output_sample_sheet = _as_path(args.output_sample_sheet)
        _write_table(sample_sheet[MANIFEST_COLUMNS], output_sample_sheet)
        print("sample_sheet written to: {}".format(output_sample_sheet))

        if args.output_report:
            output_report = _as_path(args.output_report)
            _write_table(suggestions, output_report)
            print("Full suggestion report written to: {}".format(output_report))

        if chosen is None:
            print("The requested grouping uses the default direct assignment rule.")
        else:
            print(
                "Selected suggestion {} using token positions [{}], collapsed into {} groups.".format(
                    chosen["suggestion_id"],
                    chosen["positions"],
                    chosen["group_count"],
                )
            )
        return


if __name__ == "__main__":
    main()
