# -*- coding: utf-8 -*-
"""Organize Yunnan photos: exact/near dedupe + people-focused sorting."""
from __future__ import annotations

import hashlib
import json
import os
import shutil
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageOps, UnidentifiedImageError

ROOT = Path(r"D:\GitLocal\JefferyLin1998.github.io\7749658f7c83bce1c44ce6f4a514c8c40bb57879\birthday\pictures\云南")
KEEP_DIR = ROOT / "_整理" / "人物重点"
LANDSCAPE_DIR = ROOT / "_整理" / "风景或其他"
DUP_DIR = ROOT / "_整理" / "重复备份"
NEAR_DUP_DIR = ROOT / "_整理" / "近似重复备份"
REPORT = ROOT / "_整理" / "整理报告.json"

IMG_EXTS = {".jpg", ".jpeg", ".png", ".heic", ".webp", ".bmp"}


def iter_images(root: Path):
    for p in root.rglob("*"):
        if not p.is_file():
            continue
        if "_整理" in p.parts:
            continue
        if p.suffix.lower() in IMG_EXTS:
            yield p


def file_md5(path: Path, chunk=1024 * 1024) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def average_hash(path: Path, size: int = 16) -> str | None:
    try:
        with Image.open(path) as im:
            im = ImageOps.exif_transpose(im)
            im = im.convert("L").resize((size, size), Image.Resampling.BILINEAR)
            pixels = list(im.getdata())
    except (UnidentifiedImageError, OSError, ValueError):
        return None
    avg = sum(pixels) / len(pixels)
    bits = "".join("1" if px >= avg else "0" for px in pixels)
    # pack to hex
    return f"{int(bits, 2):0{size * size // 4}x}"


def hamming(a: str, b: str) -> int:
    if a is None or b is None or len(a) != len(b):
        return 999
    x = int(a, 16) ^ int(b, 16)
    return x.bit_count()


def detect_faces(path: Path, cascade) -> tuple[int, float]:
    """Return (face_count, largest_face_area_ratio)."""
    import cv2
    import numpy as np

    try:
        with Image.open(path) as im:
            im = ImageOps.exif_transpose(im)
            im = im.convert("RGB")
            # downscale for speed
            max_side = 960
            w, h = im.size
            scale = min(1.0, max_side / max(w, h))
            if scale < 1:
                im = im.resize((int(w * scale), int(h * scale)), Image.Resampling.BILINEAR)
            arr = np.array(im)[:, :, ::-1]  # RGB->BGR
    except Exception:
        return 0, 0.0

    gray = cv2.cvtColor(arr, cv2.COLOR_BGR2GRAY)
    faces = cascade.detectMultiScale(gray, scaleFactor=1.1, minNeighbors=5, minSize=(40, 40))
    if len(faces) == 0:
        return 0, 0.0
    img_area = arr.shape[0] * arr.shape[1]
    largest = max(int(fw) * int(fh) for (_, _, fw, fh) in faces)
    return len(faces), largest / img_area


def safe_move(src: Path, dest_dir: Path) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / src.name
    if dest.exists():
        stem, suf = src.stem, src.suffix
        i = 1
        while True:
            cand = dest_dir / f"{stem}_{i}{suf}"
            if not cand.exists():
                dest = cand
                break
            i += 1
    shutil.move(str(src), str(dest))
    return dest


def main():
    files = list(iter_images(ROOT))
    print(f"found {len(files)} images")

    # 1) exact duplicates by md5
    md5_map: dict[str, list[Path]] = defaultdict(list)
    for i, p in enumerate(files, 1):
        if i % 100 == 0:
            print(f"hashing {i}/{len(files)}")
        try:
            md5_map[file_md5(p)].append(p)
        except OSError as e:
            print("skip", p.name, e)

    exact_keep = []
    exact_dups = []
    for paths in md5_map.values():
        paths = sorted(paths, key=lambda x: (-x.stat().st_size, x.name))
        exact_keep.append(paths[0])
        exact_dups.extend(paths[1:])

    print(f"exact unique={len(exact_keep)}, exact dups={len(exact_dups)}")

    # 2) near-duplicates by average hash among unique files
    hashes = {}
    for i, p in enumerate(exact_keep, 1):
        if i % 100 == 0:
            print(f"ahash {i}/{len(exact_keep)}")
        hashes[p] = average_hash(p)

    near_groups = []
    used = set()
    items = [(p, h) for p, h in hashes.items() if h]
    for i, (p, h) in enumerate(items):
        if p in used:
            continue
        group = [p]
        used.add(p)
        for q, hq in items[i + 1 :]:
            if q in used:
                continue
            if hamming(h, hq) <= 5:  # very similar
                group.append(q)
                used.add(q)
        if len(group) > 1:
            near_groups.append(group)

    near_keep = []
    near_dups = []
    for group in near_groups:
        group = sorted(group, key=lambda x: (-x.stat().st_size, x.name))
        near_keep.append(group[0])
        near_dups.extend(group[1:])

    # remaining uniques not in any near group
    near_dup_set = set(near_dups)
    final_candidates = [p for p in exact_keep if p not in near_dup_set]
    print(f"after near-dedupe candidates={len(final_candidates)}, near dups={len(near_dups)}")

    # 3) face detection if available
    people = []
    others = []
    face_stats = {}
    try:
        import cv2

        cascade_path = cv2.data.haarcascades + "haarcascade_frontalface_default.xml"
        cascade = cv2.CascadeClassifier(cascade_path)
        has_cv = True
        print("opencv face detection enabled")
    except Exception as e:
        has_cv = False
        cascade = None
        print("no opencv, skip face detection:", e)

    if has_cv:
        for i, p in enumerate(final_candidates, 1):
            if i % 50 == 0:
                print(f"faces {i}/{len(final_candidates)}")
            n, ratio = detect_faces(p, cascade)
            face_stats[str(p.name)] = {"faces": n, "largest_ratio": round(ratio, 4)}
            # people-focused: at least 1 face and largest face not tiny
            if n >= 1 and ratio >= 0.01:
                people.append((p, n, ratio))
            else:
                others.append(p)
        people.sort(key=lambda x: (-x[2], -x[1], x[0].name))
    else:
        others = list(final_candidates)

    # 4) move files
    moved = {"exact_dups": [], "near_dups": [], "people": [], "others": []}

    for p in exact_dups:
        dest = safe_move(p, DUP_DIR)
        moved["exact_dups"].append(dest.name)

    for p in near_dups:
        if p.exists():
            dest = safe_move(p, NEAR_DUP_DIR)
            moved["near_dups"].append(dest.name)

    for p, n, ratio in people:
        if p.exists():
            dest = safe_move(p, KEEP_DIR)
            moved["people"].append({"file": dest.name, "faces": n, "largest_ratio": round(ratio, 4)})

    for p in others:
        if p.exists():
            dest = safe_move(p, LANDSCAPE_DIR)
            moved["others"].append(dest.name)

    report = {
        "total_input": len(files),
        "exact_unique": len(exact_keep),
        "exact_duplicates_moved": len(moved["exact_dups"]),
        "near_duplicates_moved": len(moved["near_dups"]),
        "people_focused": len(moved["people"]),
        "landscape_or_other": len(moved["others"]),
        "people_top": moved["people"][:40],
        "folders": {
            "人物重点": str(KEEP_DIR),
            "风景或其他": str(LANDSCAPE_DIR),
            "重复备份": str(DUP_DIR),
            "近似重复备份": str(NEAR_DUP_DIR),
        },
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({k: report[k] for k in report if k != "people_top"}, ensure_ascii=False, indent=2))
    print("report:", REPORT)


if __name__ == "__main__":
    main()
