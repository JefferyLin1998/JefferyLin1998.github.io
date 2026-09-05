# -*- coding: utf-8 -*-
"""生成最近两年（滚动 730 天）聊天年度报告统计数据。

用法:
  python report_stats.py <聊天记录.txt> <输出.json>
"""
import io
import json
import re
import sys
from collections import Counter, defaultdict
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


PHRASES = {
    "爱你": ["爱你"],
    "想你": ["想你"],
    "老婆": ["老婆"],
    "老公": ["老公"],
    "宝": ["宝"],
    "哈哈哈": ["哈哈哈"],
    "嗯嗯": ["嗯嗯"],
    "好的": ["好的"],
    "晚安": ["晚安"],
    "早": ["早"],
    "抱抱": ["抱抱"],
    "亲亲": ["亲亲"],
    "啦": ["啦"],
    "嘿嘿": ["嘿嘿"],
    "呜呜": ["呜呜"],
    "吃饭": ["吃饭"],
    "奶茶": ["奶茶"],
    "蛋糕": ["蛋糕"],
    "红包": ["红包"],
    " fireworks": [],
}

EMOJI_RE = re.compile("[\U0001F300-\U0001FAFF\u2600-\u27BF\u2764\u2B50]")


def main():
    src, out = sys.argv[1], sys.argv[2]
    records = parse(src)
    if not records:
        print("no records")
        return

    end_time = records[-1][0]
    start_time = end_time - timedelta(days=730)
    window = [r for r in records if r[0] > start_time]
    print("window:", start_time, "->", end_time, "msgs:", len(window))

    days_span = (end_time.date() - start_time.date()).days

    per_day = Counter()
    per_day_she = Counter()
    per_day_me = Counter()
    hour_all = Counter()
    hour_she = Counter()
    hour_me = Counter()
    kind_counter = Counter()
    phrase_her = Counter()
    phrase_me = Counter()
    emoji_her = Counter()
    emoji_me = Counter()
    voice_her = voice_me = 0
    img_her = img_me = 0
    call_count = 0
    longest_msg_her = ("", None)
    longest_msg_me = ("", None)
    sticker_her = Counter()
    sticker_me = Counter()

    # 每日聊天时段
    night_days = set()
    dawn_days = set()
    # 回复速度：相邻两条不同人消息间隔
    reply_gaps_me = []   # 她说完 -> 我回
    reply_gaps_her = []  # 我说完 -> 她回
    last_by = None
    last_time = None
    last_date = None
    # 每天谁先说话
    first_of_day = Counter()

    # 关键词计数
    kw_counter = Counter()
    keywords = ["宝宝", "宝贝", "老婆", "老公", "爱你", "想你", "哈哈哈", "晚安",
                "早安", "吃饭", "奶茶", "蛋糕", "红包", "旅游", "婚礼", "婚纱",
                "领证", "戒指", "生日", "西安", "三亚", "浏阳", "上海"]

    for dt, kind, name, text in window:
        she = is_she(name)
        who = "she" if she else "me"
        day = dt.date().isoformat()

        per_day[day] += 1
        (per_day_she if she else per_day_me)[day] += 1
        hour_all[dt.hour] += 1
        (hour_she if she else hour_me)[dt.hour] += 1
        kind_counter[kind or "text"] += 1

        if last_date != day:
            first_of_day[who] += 1
            last_date = day

        if dt.hour >= 23 or dt.hour < 1:
            night_days.add(day)
        if 0 <= dt.hour < 6:
            dawn_days.add(day)

        if kind == "语音" or text == "[语音]":
            if she:
                voice_her += 1
            else:
                voice_me += 1
        if kind == "图片" or text == "[图片]":
            if she:
                img_her += 1
            else:
                img_me += 1
        if kind == "网络电话" or text == "[网络电话]":
            call_count += 1

        m = re.match(r"^\[([^\]]{1,8})\]$", text)
        if kind == "表情" or m:
            (sticker_her if she else sticker_me)[text] += 1

        for ch in EMOJI_RE.findall(text):
            (emoji_her if she else emoji_me)[ch] += 1

        for kw in keywords:
            if kw in text:
                kw_counter[kw] += 1

        for label, keys in PHRASES.items():
            if not keys:
                continue
            hit = any(k in text for k in keys)
            if hit:
                (phrase_her if she else phrase_me)[label] += 1

        if she and len(text) > len(longest_msg_her[0]):
            longest_msg_her = (text, dt.strftime("%Y-%m-%d"))
        if not she and len(text) > len(longest_msg_me[0]):
            longest_msg_me = (text, dt.strftime("%Y-%m-%d"))

        # 回复间隔
        if last_by and last_by != who and (dt - last_time).total_seconds() <= 3600:
            gap = (dt - last_time).total_seconds()
            if who == "me":
                reply_gaps_me.append(gap)
            else:
                reply_gaps_her.append(gap)
        last_by = who
        last_time = dt

    def avg(lst):
        if not lst:
            return 0
        return round(sum(lst) / len(lst))

    def median(lst):
        if not lst:
            return 0
        s = sorted(lst)
        return round(s[len(s) // 2])

    top_day = per_day.most_common(1)[0]
    top_day_msg = None
    # 找最忙那天的一句话
    for dt, kind, name, text in window:
        if dt.date().isoformat() == top_day[0] and len(text) > 4 and kind == "":
            top_day_msg = {"dt": dt.strftime("%m-%d %H:%M"), "who": "she" if is_she(name) else "me", "text": text[:60]}
            break

    night_count_she = sum(v for h, v in hour_she.items() if h >= 23 or h < 1)
    night_count_me = sum(v for h, v in hour_me.items() if h >= 23 or h < 1)

    report = {
        "generatedAt": end_time.strftime("%Y-%m-%d %H:%M"),
        "range": {
            "start": start_time.date().isoformat(),
            "end": end_time.date().isoformat(),
            "days": days_span,
            "label": "最近两年",
        },
        "total": len(window),
        "activeDays": len(per_day),
        "avgPerDay": round(len(window) / max(1, len(per_day))),
        "talker": {"she": sum(per_day_she.values()), "me": sum(per_day_me.values())},
        "chars": {
            "she": sum(len(t) for dt, k, n, t in window if is_she(n)),
            "me": sum(len(t) for dt, k, n, t in window if not is_she(n)),
        },
        "busiestDay": {"date": top_day[0], "count": top_day[1], "msg": top_day_msg},
        "peakHour": max(hour_all, key=hour_all.get),
        "hourly": {str(h): hour_all[h] for h in range(24)},
        "hourlyShe": {str(h): hour_she[h] for h in range(24)},
        "hourlyMe": {str(h): hour_me[h] for h in range(24)},
        "nightOwl": {
            "lateMsgsShe": night_count_she,
            "lateMsgsMe": night_count_me,
            "lateDays": len(night_days | dawn_days),
        },
        "firstOfDay": dict(first_of_day),
        "replySpeed": {
            "meAvgSec": avg(reply_gaps_me),
            "herAvgSec": avg(reply_gaps_her),
            "meMedianSec": median(reply_gaps_me),
            "herMedianSec": median(reply_gaps_her),
            "meCount": len(reply_gaps_me),
            "herCount": len(reply_gaps_her),
        },
        "media": {
            "voiceShe": voice_her, "voiceMe": voice_me,
            "imgShe": img_her, "imgMe": img_me,
            "calls": call_count,
        },
        "phrasesHer": phrase_her.most_common(12),
        "phrasesMe": phrase_me.most_common(12),
        "emojiHer": [(c, n) for c, n in emoji_her.most_common(20) if len(c) > 1 or not c.isascii()],
        "emojiMe": [(c, n) for c, n in emoji_me.most_common(20) if len(c) > 1 or not c.isascii()],
        "stickersHer": sticker_her.most_common(10),
        "stickersMe": sticker_me.most_common(10),
        "keywords": kw_counter.most_common(20),
        "longestHer": {"text": longest_msg_her[0][:200], "date": longest_msg_her[1]},
        "longestMe": {"text": longest_msg_me[0][:200], "date": longest_msg_me[1]},
        "monthly": {},
    }

    # 月度统计
    per_month = Counter()
    for day, cnt in per_day.items():
        per_month[day[:7]] += cnt
    report["monthly"] = {k: v for k, v in sorted(per_month.items())}

    with io.open(out, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=1)
    print("wrote", out)


if __name__ == "__main__":
    main()
