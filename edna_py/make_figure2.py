"""
Figure 2: screenshots of the Python edition of the application with ND-LEEF 2021 data loaded
(early open phase, pea gravel stream, steelhead trout), taken with headless Chromium (Playwright)
against a running instance, then composed into one six-panel figure.

Usage: python3 make_figure2.py <app url> <results dir> <figures dir>
       (start the app first: python -m shiny run --port 7799 app.py)
"""
import sys, os, time
from playwright.sync_api import sync_playwright
from PIL import Image, ImageDraw, ImageFont

url, res, out = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out, exist_ok=True)
samp = os.path.abspath(f"{res}/app_inputs/LEEF2021_sampling_perSize.csv")
TREAT, TARGET, F04 = "LB_PG", "STLHD", "0.306"                 # first-station 0.4 um fraction (reference-station design)
SMIM = {"Us": "0.1006", "Ds": "0.0303", "lambdas": "0.085", "betas": "0.845", "rss": "2.1e-5"}   # early open, pea gravel

with sync_playwright() as p:
    b = p.chromium.launch()
    pg = b.new_page(viewport={"width": 1400, "height": 1000}, device_scale_factor=2)
    pg.goto(url); pg.wait_for_selector("text=eDNA analysis tool"); time.sleep(1)
    pg.set_input_files("#file1", samp); time.sleep(2)
    pg.select_option("#treat_var", TREAT); pg.select_option("#target_var", TARGET); time.sleep(1.5)
    pg.fill("#fr00", F04); pg.keyboard.press("Tab"); time.sleep(1)

    def shot(name, full=False, h=None):
        path = f"{out}/{name}.png"
        if full:
            pg.screenshot(path=path, full_page=True)
        else:
            pg.locator("main").first.screenshot(path=path)
        im = Image.open(path)
        if h:
            im = im.crop((0, 0, im.width, min(im.height, h))); im.save(path)
        print("wrote", path); return path

    shots = {}
    pg.get_by_role("tab", name="Sampling Data").click(); time.sleep(1.5); shots["a"] = shot("fig2a_sampling", full=True, h=1180)
    pg.get_by_role("tab", name="Source Data").click(); time.sleep(1.2); shots["b"] = shot("fig2b_source", h=520)
    pg.get_by_role("tab", name="Prediction").click(); time.sleep(1); pg.click("#submit"); time.sleep(9); shots["c"] = shot("fig2c_prediction", h=1700)
    pg.get_by_role("tab", name="Sensitivity").click(); time.sleep(1.5); shots["d"] = shot("fig2d_sensitivity", h=1150)
    pg.get_by_role("tab", name="LOD").click(); time.sleep(1.5); shots["e"] = shot("fig2e_lod", h=1150)
    pg.get_by_role("tab", name="SMIM").click(); time.sleep(1)
    for k, v in SMIM.items():
        pg.fill(f"#{k}", v); pg.keyboard.press("Tab")
    pg.click("body"); time.sleep(2.5); shots["f"] = shot("fig2f_smim")
    b.close()

# ---- composite
ims = {k: Image.open(v) for k, v in shots.items()}
tw = 1400
scaled = {k: im.resize((tw, int(im.height * tw / im.width)), Image.LANCZOS) for k, im in ims.items()}
order = [["a", "b"], ["c", "d"], ["e", "f"]]
gap = 36; band = 70
H = sum(max(scaled[k].height for k in r) for r in order) + (gap + band) * len(order) + gap
canvas = Image.new("RGB", (2 * tw + 3 * gap, H), "white"); dr = ImageDraw.Draw(canvas)
try:
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 46)
except Exception:
    font = ImageFont.load_default()
y = gap
for r in order:
    x = gap
    for k in r:
        dr.text((x, y), k, fill="black", font=font); canvas.paste(scaled[k], (x, y + band)); x += tw + gap
    y += max(scaled[k].height for k in r) + gap + band
dest = os.path.join(os.path.dirname(os.path.abspath(out.rstrip("/"))), "fig2_app_screens.png")
canvas.save(dest); print("wrote", dest)
