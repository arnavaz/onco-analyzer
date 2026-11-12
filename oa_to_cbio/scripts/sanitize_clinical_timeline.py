import argparse
import csv
import os
import tempfile
from pathlib import Path


def sanitize_file(filepath, apply, strip_hash_rows, deduplicate):
    seen_rows = set()
    header_seen = False
    duplicate_count = 0
    hash_rows_remaining = 4 if strip_hash_rows else 0
    hash_rows_stripped = 0
    temporary_path = None

    try:
        with filepath.open("r", encoding="utf-8", newline="") as source:
            if apply:
                temporary_file = tempfile.NamedTemporaryFile(
                    mode="w",
                    encoding="utf-8",
                    newline="",
                    dir=filepath.parent,
                    delete=False,
                )
                temporary_path = Path(temporary_file.name)
            else:
                temporary_file = None

            try:
                for line in source:
                    content = line.rstrip("\r\n")

                    if hash_rows_remaining:
                        if not content.startswith("#"):
                            raise ValueError(
                                "Expected four leading hash-prefixed metadata rows."
                            )
                        hash_rows_remaining -= 1
                        hash_rows_stripped += 1
                        continue

                    if not content or content.startswith("#"):
                        if temporary_file:
                            temporary_file.write(line)
                        continue

                    row = next(csv.reader([content], delimiter="\t", quotechar='"'))
                    if not header_seen:
                        header_seen = True
                    elif deduplicate and tuple(row) in seen_rows:
                        duplicate_count += 1
                        continue
                    elif deduplicate:
                        seen_rows.add(tuple(row))

                    if temporary_file:
                        temporary_file.write(line)
            finally:
                if temporary_file:
                    temporary_file.close()

        if hash_rows_remaining:
            raise ValueError("Expected four leading hash-prefixed metadata rows.")

        if apply and temporary_path:
            if duplicate_count or hash_rows_stripped:
                os.chmod(temporary_path, filepath.stat().st_mode)
                os.replace(temporary_path, filepath)
            else:
                temporary_path.unlink()
            temporary_path = None

        actions = []
        if strip_hash_rows:
            metadata_verb = "Stripped" if apply else "Would strip"
            actions.append(f"{metadata_verb} {4 - hash_rows_remaining} leading hash rows")
        if deduplicate:
            duplicate_verb = "Removed" if apply else "Found"
            actions.append(f"{duplicate_verb} {duplicate_count} duplicate rows")
        print(f"{'; '.join(actions)} in {filepath}.")
    finally:
        if temporary_path and temporary_path.exists():
            temporary_path.unlink()


def main():
    parser = argparse.ArgumentParser(
        description="Clean timeline TSV files or all timelines in an import directory."
    )
    parser.add_argument("path", type=Path, help="Timeline file or import directory")
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Rewrite files; otherwise only report what would change.",
    )
    parser.add_argument(
        "--strip-hash-rows",
        action="store_true",
        help="Remove four leading hash-prefixed metadata rows from a single file.",
    )
    parser.add_argument(
        "--deduplicate",
        action="store_true",
        help="Remove exact duplicate data rows from a single file.",
    )
    args = parser.parse_args()

    if args.path.is_dir():
        timeline_files = sorted(args.path.glob("data_clinical_timeline_*.txt"))
        strip_hash_rows = True
        if not timeline_files:
            print(f"No clinical timeline files found in {args.path}.")
            return
    elif args.path.is_file():
        timeline_files = [args.path]
        strip_hash_rows = args.strip_hash_rows
    else:
        parser.error(f"Path not found or not a file/directory: {args.path}")

    for filepath in timeline_files:
        if args.path.is_dir():
            deduplicate = filepath.name == "data_clinical_timeline_specimen.txt"
        else:
            deduplicate = args.deduplicate or not args.strip_hash_rows
        sanitize_file(filepath, args.apply, strip_hash_rows, deduplicate)


if __name__ == "__main__":
    main()