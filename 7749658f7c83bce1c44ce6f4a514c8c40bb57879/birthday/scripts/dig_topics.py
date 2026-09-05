# -*- coding: utf-8 -*-
"""在最近两年聊天记录中挖掘出题素材片段。"""
import io
import re
import sys
from datetime import datetime, timedelta

LINE_RE = re.compile(
    r"^\[(\d{4})-(\d{2})-(\d{2}) (\d{2}):(\d{2}):(\d{2})\]\s*(?:\[([^\]]+)\]\s*)?([^:：]+?)\s*[:：]\s*(.*)$"
)

def parse(path):
    records = []
    with io.open(path, "r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            m = LINE_RE.match(raw.rstrip("\r\n"))
            if not m:
                continue
            y, mo, d, h, mi, s, kind, name, text = m.groups()
            records.append((
                datetime(int(y), int(mo), int(d), int(h), int(mi), int(s)),
                kind or "",
                name.strip(),
                text.strip(),
            ))
    return records

def is_she(name):
    return "宝" in name or "莎" in name


def main():
    src = sys.argv[1]
    queries = sys.argv[2:]
    records = parse(src)
    end_time = records[-1][0]
    start_time = end_time - timedelta(days=730)
    window = [r for r in records if r[0] > start_time]

    out = io.open(sys.argv[-1] + ".out.txt", "w", encoding="utf-8") if len(sys.argv) > 3 else sys.stdout

    for q in queries:
        out.write("\n===== " + q + " =====\n")
        count = 0
        for dt, kind, name, text in window:
            if q in text:
                who = "她" if is_she(name) else "我"
                out.write("%s [%s] %s: %s\n" % (
                    dt.strftime("%Y-%m-%d %H:%M"), who, kind or "文本", text[:150]))
                count += 1
                if count >= 40:
                    break
        if count == 0:
            out.write("(无结果)\n")
    out.flush()


if __name__ == "__main__":
    main()
