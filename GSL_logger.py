import requests
import os
from datetime import datetime
from openpyxl import Workbook, load_workbook

# ===== CONFIG =====
OWNER = "ahuelsmeier"
REPO = "gsl-transition-generator"
TOKEN = os.getenv("GITHUB_TOKEN")

# REPLACE THIS with your actual Zenodo Record ID (the number at the end of the URL)
ZENODO_RECORD_ID = "YOUR_ZENODO_ID_HERE" 

if not TOKEN:
    raise ValueError("Missing GITHUB_TOKEN environment variable")

EXCEL_FILE = "github_traffic_log.xlsx"

GITHUB_BASE_URL = f"https://api.github.com/repos/{OWNER}/{REPO}/traffic"
ZENODO_BASE_URL = f"https://zenodo.org/api/records/{ZENODO_RECORD_ID}"

headers_gh = {
    "Authorization": f"Bearer {TOKEN}",
    "Accept": "application/vnd.github+json"
}

# ===== FETCH FUNCTIONS =====
def fetch_github(endpoint):
    url = f"{GITHUB_BASE_URL}/{endpoint}"
    r = requests.get(url, headers=headers_gh)
    r.raise_for_status()
    return r.json()

def fetch_zenodo():
    """Fetches cumulative stats from Zenodo."""
    r = requests.get(ZENODO_BASE_URL)
    r.raise_for_status()
    data = r.json()
    return data.get("stats", {})

def get_today_summary():
    views = fetch_github("views")
    clones = fetch_github("clones")
    zenodo = fetch_zenodo()

    today = datetime.utcnow().strftime("%Y-%m-%d")

    return {
        "date": today,
        "views_total": views.get("count", 0),
        "views_unique": views.get("uniques", 0),
        "clones_total": clones.get("count", 0),
        "clones_unique": clones.get("uniques", 0),
        "zenodo_views": zenodo.get("views", 0),
        "zenodo_downloads": zenodo.get("downloads", 0)
    }

# ===== EXCEL HANDLING =====
def init_excel():
    wb = Workbook()
    ws = wb.active
    ws.title = "Traffic"

    # Updated headers to include Zenodo
    headers = [
        "date",
        "views_total",
        "views_unique",
        "clones_total",
        "clones_unique",
        "zenodo_views",
        "zenodo_downloads"
    ]
    ws.append(headers)
    wb.save(EXCEL_FILE)

def append_to_excel(data):
    if not os.path.exists(EXCEL_FILE):
        init_excel()

    wb = load_workbook(EXCEL_FILE)
    ws = wb["Traffic"]

    # Avoid duplicate entries for the same day
    dates = [row[0].value for row in ws.iter_rows(min_row=2)]
    if data["date"] in dates:
        print(f"Entry for {data['date']} already exists. Skipping.")
        return

    ws.append([
        data["date"],
        data["views_total"],
        data["views_unique"],
        data["clones_total"],
        data["clones_unique"],
        data["zenodo_views"],
        data["zenodo_downloads"]
    ])

    wb.save(EXCEL_FILE)
    print(f"Logged GitHub and Zenodo data for {data['date']}")

# ===== MAIN =====
if __name__ == "__main__":
    data = get_today_summary()
    append_to_excel(data)
