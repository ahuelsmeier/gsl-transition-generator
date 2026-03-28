import requests
import os
from datetime import datetime
from openpyxl import Workbook, load_workbook

# ===== CONFIG =====
OWNER = "ahuelsmeier"
REPO = "gsl-transition-generator"
# TOKEN = "ghp_your_token"
TOKEN = os.getenv("GITHUB_TOKEN")

if not TOKEN:
    raise ValueError("Missing GITHUB_TOKEN environment variable")

EXCEL_FILE = "github_traffic_log.xlsx"

BASE_URL = f"https://api.github.com/repos/{OWNER}/{REPO}/traffic"

headers = {
    "Authorization": f"Bearer {TOKEN}",
    "Accept": "application/vnd.github+json"
}

ZENODO_RECORD_ID = "18901586"

# ===== FETCH FUNCTIONS =====

def fetch_zenodo_stats():
    """Fetches views and downloads from Zenodo."""
    url = f"https://zenodo.org/api/records/{ZENODO_RECORD_ID}"
    r = requests.get(url) # Zenodo stats are public; no token required for GET
    r.raise_for_status()
    data = r.json()

    stats = data.get("stats", {})
    return {
        "zenodo_views": stats.get("views", 0),
        "zenodo_downloads": stats.get("downloads", 0),
        "zenodo_unique_views": stats.get("unique_views", 0),
        "zenodo_unique_downloads": stats.get("unique_downloads", 0)
    }

def fetch(endpoint):
    url = f"{BASE_URL}/{endpoint}"
    r = requests.get(url, headers=headers)
    r.raise_for_status()
    return r.json()

def get_today_summary():
    views = fetch("views")
    clones = fetch("clones")

    today = datetime.utcnow().strftime("%Y-%m-%d")

    return {
        "date": today,
        "views_total": views.get("count", 0),
        "views_unique": views.get("uniques", 0),
        "clones_total": clones.get("count", 0),
        "clones_unique": clones.get("uniques", 0),
    }

# ===== EXCEL HANDLING =====
def init_excel():
    wb = Workbook()
    ws = wb.active
    ws.title = "Traffic"

    headers = [
        "date",
        "views_total",
        "views_unique",
        "clones_total",
        "clones_unique"
    ]
    ws.append(headers)
    wb.save(EXCEL_FILE)

def append_to_excel(data):
    if not os.path.exists(EXCEL_FILE):
        init_excel()

    wb = load_workbook(EXCEL_FILE)
    ws = wb["Traffic"]

    # Avoid duplicate entries
    dates = [row[0].value for row in ws.iter_rows(min_row=2)]
    if data["date"] in dates:
        print(f"Entry for {data['date']} already exists. Skipping.")
        return

    ws.append([
        data["date"],
        data["views_total"],
        data["views_unique"],
        data["clones_total"],
        data["clones_unique"]
    ])

    wb.save(EXCEL_FILE)
    print(f"Logged data for {data['date']}")

# ===== MAIN =====
if __name__ == "__main__":
    data = get_today_summary()
    append_to_excel(data)
