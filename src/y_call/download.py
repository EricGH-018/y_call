import argparse
import pathlib
import urllib.request

def download_database(database: str):
    urls = {
        "YFull": "https://github.com/EricGH-018/y_call/releases/download/v1.0.0/YFull.zip",
        "OY": "https://github.com/EricGH-018/y_call/releases/download/v1.0.0/OY.zip",
        "translation": "https://github.com/EricGH-018/y_call/releases/download/v1.0.0/YF-translations.csv",
        "ages": "https://github.com/EricGH-018/y_call/releases/download/v1.0.0/ages.csv"
    }

    if database not in urls:
        raise ValueError(f"Unknown database '{database}'. Available: {list(urls.keys())}")

    url = urls[database]

    out_dir = pathlib.Path.home() / "data" / "input"
    out_dir.mkdir(parents=True, exist_ok=True)

    dest = out_dir / f"{database}.zip"

    print("\n--- RUNNING ANALYSIS TO DOWNLOAD REFERENCE DATASETS ---\n")

    print(f"Downloading {database} dataset...")
    urllib.request.urlretrieve(url, dest)
    print(f"Saved to {dest}")

def main():
    parser = argparse.ArgumentParser(description="Download reference datasets for y_call")
    parser.add_argument("--database", required=True, help="Name of the dataset to download (YFull, OY, translation or ages)")
    args = parser.parse_args()

    download_database(args.database)

