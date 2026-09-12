# -*- coding: utf-8 -*-
"""分析两人的聊天风格与兴趣，生成人物画像素材。

用法:
  python profile_stats.py <聊天记录.txt> <输出.json>
"""
import io
import json
import re
import sys
from collections import Counter

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
            records.append((kind or "", name.strip(), text.strip()))
    return records


def is_she(name):
    return "宝" in name or "莎" in name


# 兴趣话题分组
TOPIC_GROUPS = {
    "美食干饭": ["吃饭", "吃啥", "好吃", "晚饭", "午饭", "早饭", "外卖", "干饭", "食堂", "火锅", "烧烤", "面条", "米饭", "饿了"],
    "奶茶咖啡": ["奶茶", "喜茶", "霸王茶姬", "星巴克", "咖啡", "拿铁", "美式", "抹茶", "啵啵"],
    "购物剁手": ["买", "购物", "淘宝", "拼多多", "京东", "下单", "退货", "快递", "到货", "种草", "拔草"],
    "薅羊毛": ["券", "红包", "优惠", "打折", "买一送一", "免费", "白嫖", "薅", "立减", "积分"],
    "旅游出行": ["旅游", "玩", "景点", "机票", "酒店", "火车", "高铁", "攻略", "打卡", "旅行"],
    "工作吐槽": ["上班", "下班", "加班", "老板", "同事", "开会", "项目", "领导", "实习", "辞职", "简历"],
    "游戏动漫": ["游戏", "王者", "原神", " Switch", "switch", "动漫", "番", "漫画", "追番", "steam", "Steam"],
    "影视综艺": ["电影", "剧", "综艺", "看片", "追剧", "电视剧", "纪录片", "Netflix", "视频"],
    "拍照修图": ["照片", "拍照", "修图", "自拍", "相机", "p图", "P图"],
    "养生健康": ["睡觉", "熬夜", "养生", "体检", "感冒", "发烧", "头疼", "肚子", "拉肚子", "减肥", "运动", "散步"],
    "家人朋友": ["妈妈", "爸", "家里", "婆婆", "奶奶", "外婆", "红燕", "萌雪", "朱雅文", "龙哥"],
    "音乐": ["歌", "听歌", "唱歌", "网易云", "歌单", "演唱会"],
    "数码科技": ["电脑", "手机", "苹果", "iPhone", " iPad", "显卡", "CPU", "键盘", "耳机", "大疆"],
    "理财搞钱": ["股票", "基金", "理财", "赚钱", "存钱", "省钱", "预算", "花呗", "还贷", "房价"],
    "婚礼家庭": ["婚礼", "婚纱", "结婚", "领证", "戒指", "婚", "宝宝", "怀孕", "生娃"],
    "穿搭美妆": ["衣服", "裙子", "口红", "化妆", "护肤品", "香水", "穿搭", "包包", "鞋子", "高跟鞋"],
}

# 性格特征探测词
TRAIT_SIGNALS = {
    "哈哈大笑型": ["哈哈", "哈哈哈", "哈哈哈哈"],
    "颜文字可爱型": ["~", "～", "嘻嘻", "嘿嘿", "啦啦", "呀"],
    "撒娇粘人型": ["抱抱", "么么", "亲亲", "想你", "想你啦", "嘛", "哼"],
    "操心唠叨型": ["记得", "注意", "小心", "别忘", "早点睡", "多喝", "穿暖"],
    "严谨分析型": ["分析", "查了", "看了下", "研究了", "对比", "总结", "计划"],
    "情绪感叹型": ["啊啊啊", "天哪", "我的天", "救命", "绝了", "牛", "卧槽", "离谱"],
    "务实直球型": ["好的", "行", "嗯", "ok", "OK", "可以", "没问题"],
    "好奇提问型": ["为什么", "怎么", "？", "啥"],
    "甜蜜表达型": ["爱你", "喜欢", "宝", "老婆", "老公", "亲爱"],
}

EMOJI_RE = re.compile("[\U0001F300-\U0001FAFF\U00002600-\U000027BF]")


def main():
    src, out = sys.argv[1], sys.argv[2]
    records = parse(src)
    print("total:", len(records))

    stats = {"she": {"topics": Counter(), "traits": Counter(), "emoji": Counter(),
                     "sticker": Counter(), "chars": 0, "msgs": 0, "long_msgs": 0,
                     "question": 0, "exclaim": 0, "sticker_msgs": 0},
             "me": {"topics": Counter(), "traits": Counter(), "emoji": Counter(),
                    "sticker": Counter(), "chars": 0, "msgs": 0, "long_msgs": 0,
                    "question": 0, "exclaim": 0, "sticker_msgs": 0}}

    for kind, name, text in records:
        who = "she" if is_she(name) else "me"
        s = stats[who]
        s["msgs"] += 1
        s["chars"] += len(text)

        if len(text) >= 50:
            s["long_msgs"] += 1
        if "？" in text or "?" in text:
            s["question"] += 1
        if "！" in text or "!" in text:
            s["exclaim"] += 1

        for ch in EMOJI_RE.findall(text):
            s["emoji"][ch] += 1

        if kind == "表情" or re.match(r"^\[[^\]]{1,6}\]$", text):
            s["sticker_msgs"] += 1
            s["sticker"][text] += 1

        for group, words in TOPIC_GROUPS.items():
            if any(w in text for w in words):
                s["topics"][group] += 1

        for trait, words in TRAIT_SIGNALS.items():
            hit = 0
            for w in words:
                if w in text:
                    hit += 1
            if hit:
                s["traits"][trait] += 1

    result = {}
    for who in ("she", "me"):
        s = stats[who]
        result[who] = {
            "msgs": s["msgs"],
            "avgLen": round(s["chars"] / max(1, s["msgs"]), 1),
            "longMsgs": s["long_msgs"],
            "longRatio": round(s["long_msgs"] / max(1, s["msgs"]) * 100, 1),
            "questionRatio": round(s["question"] / max(1, s["msgs"]) * 100, 1),
            "exclaimRatio": round(s["exclaim"] / max(1, s["msgs"]) * 100, 1),
            "stickerMsgs": s["sticker_msgs"],
            "topics": s["topics"].most_common(10),
            "traits": s["traits"].most_common(10),
            "emoji": s["emoji"].most_common(10),
            "stickers": s["sticker"].most_common(8),
        }

    with io.open(out, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=1)
    print("wrote", out)

    for who, label in (("she", "Sasa"), ("me", "Lin")):
        r = result[who]
        print("\n===== %s =====" % label)
        print("消息数:", r["msgs"], " 平均长度:", r["avgLen"])
        print("长消息占比: %s%%  问句占比: %s%%  感叹占比: %s%%" % (r["longRatio"], r["questionRatio"], r["exclaimRatio"]))
        print("话题TOP:", r["topics"][:6])
        print("特质TOP:", r["traits"][:6])
        print("emoji:", r["emoji"][:6])
        print("表情:", r["stickers"][:5])


if __name__ == "__main__":
    main()
