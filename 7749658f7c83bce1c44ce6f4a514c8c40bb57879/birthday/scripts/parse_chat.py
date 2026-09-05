# -*- coding: utf-8 -*-
"""解析微信聊天记录导出文件，提取统计信息和出题素材。

用法:
  python parse_chat.py <聊天记录.txt> <输出目录>
"""
import io
import json
import os
import re
import sys
from collections import Counter, defaultdict

LINE_RE = re.compile(
    r"^\[(\d{4})-(\d{2})-(\d{2}) (\d{2}):(\d{2}):(\d{2})\]\s*(?:\[([^\]]+)\]\s*)?([^:：]+?)\s*[:：]\s*(.*)$"
)


def parse(path):
    records = []
    with io.open(path, "r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n").rstrip("\r")
            m = LINE_RE.match(line)
            if not m:
                continue
            y, mo, d, h, mi, s, kind, name, text = m.groups()
            records.append({
                "dt": (int(y), int(mo), int(d), int(h), int(mi), int(s)),
                "kind": kind or "",
                "name": name.strip(),
                "text": text.strip(),
            })
    return records


def is_she(name):
    return "宝" in name or "莎" in name


def main():
    src = sys.argv[1]
    out_dir = sys.argv[2]
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    records = parse(src)
    print("total parsed:", len(records))

    by_year = Counter()
    by_month = defaultdict(Counter)     # year -> month -> count
    by_hour = Counter()
    talker = Counter()
    emoji_counter = Counter()
    sticker_her = Counter()
    sticker_me = Counter()
    day_active = defaultdict(set)       # year -> set of (y,m,d)
    char_count = Counter()              # "she"/"me" -> chars
    msg_per_day_top = Counter()         # (y,m,d) -> count

    # 最近两年 = 2024-09 ~ 2026-09 (导出时间 2026-08-08，取 2025 全年 + 2024 部分)
    # 报告年度区间直接按 2025-01-01 ~ 2025-12-31 + 2026-01-01 ~ 2026-08-08 统计

    samples = defaultdict(list)         # keyword -> list of (dt, name, text)
    keywords = [
        "青海湖", "毕业旅行", "男朋友", "女朋友", "浏阳之光", "领证", "求婚",
        "老公", "老婆", "爱你", "喜欢你", "想你", "第一次", "生日快乐",
        "西安", "三亚", "上海", "戒指", "婚纱", "婚礼",
    ]

    emoji_re = re.compile(
        "[\U0001F300-\U0001FAFF\u2600-\u27BF\uFE0F\u200D\u2B50\u2764]"
    )

    for r in records:
        y, mo, d, h, mi, s = r["dt"]
        she = is_she(r["name"])
        who = "she" if she else "me"
        by_year[y] += 1
        by_month[y][mo] += 1
        by_hour[h] += 1
        talker[who] += 1
        char_count[who] += len(r["text"])
        day_active[y].add((y, mo, d))
        msg_per_day_top[(y, mo, d)] += 1

        text = r["text"]
        for ch in emoji_re.findall(text):
            emoji_counter[ch] += 1
        if r["kind"] == "表情" or re.match(r"^\[[^\]]{1,6}\]$", text):
            (sticker_her if she else sticker_me)[text] += 1

        for kw in keywords:
            if kw in text:
                if len(samples[kw]) < 60:
                    samples[kw].append({
                        "dt": "%04d-%02d-%02d %02d:%02d" % (y, mo, d, h, mi),
                        "who": who,
                        "text": text[:120],
                    })

    # 各年统计
    years = sorted(by_year)
    stats = {
        "total": len(records),
        "years": {
            str(y): {
                "messages": by_year[y],
                "activeDays": len(day_active[y]),
                "months": {str(m): by_month[y][m] for m in range(1, 13)},
            } for y in years
        },
        "talker": dict(talker),
        "chars": dict(char_count),
        "byHour": {str(h): by_hour[h] for h in range(24)},
        "topEmoji": emoji_counter.most_common(30),
        "topStickersHer": sticker_her.most_common(15),
        "topStickersMe": sticker_me.most_common(15),
        "busiestDays": [
            {"date": "%04d-%02d-%02d" % k, "count": v}
            for k, v in msg_per_day_top.most_common(15)
        ],
    }

    def dump(name, obj):
        p = os.path.join(out_dir, name)
        with io.open(p, "w", encoding="utf-8") as f:
            if isinstance(obj, str):
                f.write(obj)
            else:
                json.dump(obj, f, ensure_ascii=False, indent=1)
        print("wrote", p)

    dump("chat_stats.json", stats)
    dump("chat_samples.json", dict(samples))

    # 控制台摘要
    print("\n== 年度消息数 ==")
    for y in years:
        print(y, by_year[y], "activeDays:", len(day_active[y]))
    print("\n== 说话比例 ==", dict(talker))
    print("\n== 最忙的15天 ==")
    for k, v in msg_per_day_top.most_common(15):
        print("%04d-%02d-%02d" % k, v)
    print("\n== 高频emoji前15 ==")
    for ch, n in emoji_counter.most_common(15):
        print(repr(ch), n)
    print("\n== 她的常用表情(文字型)前10 ==")
    for t, n in sticker_her.most_common(10):
        print(t, n)


if __name__ == "__main__":
    main()
