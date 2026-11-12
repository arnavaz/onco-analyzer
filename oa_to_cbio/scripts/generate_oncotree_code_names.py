import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API_URL = "https://oncotree.mskcc.org/api/tumorTypes"
VERSION = "oncotree_latest_stable"
OUTPUT_PATH = Path(__file__).with_name("oncotree_code_names.json")


def main():
    url = f"{API_URL}?{urlencode({'version': VERSION})}"
    request = Request(url, headers={"User-Agent": "OncoTree-code-map-generator"})

    with urlopen(request, timeout=30) as response:
        tumor_types = json.load(response)

    if not isinstance(tumor_types, list):
        raise RuntimeError("Unexpected OncoTree API response; expected a list.")

    code_to_name = {
        tumor_type["code"]: tumor_type["name"]
        for tumor_type in tumor_types
        if tumor_type.get("code") and tumor_type.get("name")
    }
    if not code_to_name:
        raise RuntimeError("OncoTree API returned no tumor type code/name pairs.")

    OUTPUT_PATH.write_text(
        json.dumps(code_to_name, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {len(code_to_name)} OncoTree code/name pairs to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()