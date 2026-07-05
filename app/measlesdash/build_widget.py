#!/usr/bin/env python3
"""Inline the dashboard data into the widget template -> a self-contained page.

widget.html carries a __DASH_DATA__ placeholder in a <script type=application/json>
block; this replaces it with data/dashboard.json so the result needs no network
(the static-baked dexetera pattern). Output: app/measlesdash/index.html.

Usage: python3 app/measlesdash/build_widget.py
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    template = open(os.path.join(HERE, "widget.html")).read()
    data = open(os.path.join(HERE, "data", "dashboard.json")).read()
    # Guard against a stray </script> prematurely closing the JSON block.
    data = data.replace("</", "<\\/")
    html = template.replace("__DASH_DATA__", data)
    out = os.path.join(HERE, "index.html")
    open(out, "w").write(html)
    print(f"wrote {out}: {os.path.getsize(out)} bytes")


if __name__ == "__main__":
    sys.exit(main())
