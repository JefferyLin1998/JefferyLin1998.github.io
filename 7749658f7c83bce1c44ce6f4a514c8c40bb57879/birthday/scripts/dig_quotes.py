# -*- coding: utf-8 -*-
"""挖掘两人的代表性语录。

用法:
  python dig_quotes.py <聊天记录.txt>
"""
import io
import re
import sys
from collections import Counter

LINE_RE = re.compile(
    r"^\[(\d{4})-(\d{2})-(\d{2}) (\d{2}):(\d{2}):(\d{2})\]\s*(?:\[([^\]]+)\]\s*)?([^:：]+?)\s*[:：]\s*(.*)$"
)


def is_she(name):
    return "宝" in name or "莎" in name


def main():
    src = sys.argv[1]
    she_quotes = Counter()
    me_quotes = Counter()
    # 常用语录特征
    she_patterns = {
        "安排规划": ["我要", "我打算", "我计划", "我来", "我先"],
        "分享欲": ["你看", "给你看", "你看这个", "跟你说", "我跟你说"],
        "砍价省钱": ["便宜", "划算", "省了", "减了", "券"],
        "牵挂关心": ["你记得", "别忘了", "你要", "多穿", "早点"],
    }
    me_patterns = {
        "宠溺回应": ["好的老婆", "听你的", "都行", "你说了算", "可以呀"],
        "行动派": ["我来", "我买", "我去", "我订", "我点"],
        "简单直接": ["好的", "行", "嗯嗯", "OK", "ok"],
    }

    with io.open(src, "r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            m = LINE_RE.match(raw.rstrip("\r\n"))
            if not m:
                continue
            y, mo, d, h, mi, s, kind, name, text = m.groups()
            text = text.strip()
            if not text or (kind or ""):
                continue
            she = is_she(name)
            pats = she_patterns if she else me_patterns
            target = she_quotes if she else me_quotes
            for group, words in pats.items():
                for w in words:
                    if w in text and 8 <= len(text) <= 60:
                        target[group + " | " + text[:50]] += 1
                        break

    for label, quotes in (("Sasa", she_quotes), ("Lin", me_quotes)):
        print("\n===== %s 代表语录 =====" % label)
        seen = set()
        count = 0
        for q, n in quotes.most_common(300):
            key = q.split(" | ", 1)[1][:20]
            if key in seen:
                continue
            seen.add(key)
            print("(%d次) %s" % (n, q))
            count += 1
            if count >= 25:
                break


if __name__ == "__main__":
    main()
