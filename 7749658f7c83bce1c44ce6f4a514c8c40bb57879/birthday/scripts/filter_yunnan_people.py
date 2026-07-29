# -*- coding: utf-8 -*-
"""Face-detect Yunnan photos already in 风景或其他, move people-focused out."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import cv2
import numpy as np
from PIL import Image, ImageOps, UnidentifiedImageError

ROOT = Path(r"D:\GitLocal\JefferyLin1998.github.io\7749658f7c83bce1c44ce6f4a514c8c40bb57879\birthday\pictures\云南")
SRC = ROOT / "_整理" / "风景或其他"
KEEP = ROOT / "_整理" / "人物重点"
REPORT = ROOT / "_整理" / "人物筛选报告.json"

IMG_EXTS = {".jpg", ".jpeg", ".png", ".heic", ".webp", ".bmp"}


def detect_faces(path: Path, cascade) -> tuple[int, float]:
    try:
        with Image.open(path) as im:
            im = ImageOps.exif_transpose(im)
            im = im.convert("RGB")
            max_side = 960
            w, h = im.size
            scale = min(1.0, max_side / max(w, h))
            if scale < 1:
                im = im.resize((int(w * scale), int(h * scale)), Image.Resampling.BILINEAR)
            arr = np.array(im)[:, :, ::-1]
    except (UnidentifiedImageError, OSError, ValueError):
        return 0, 0.0

    gray = cv2.cvtColor(arr, cv2.COLOR_BGR2GRAY)
    faces = cascade.detectMultiScale(gray, scaleFactor=1.1, minNeighbors=5, minSize=(40, 40))
    if len(faces) == 0:
        return 0, 0.0
    img_area = arr.shape[0] * arr.shape[1]
    largest = max(int(fw) * int(fh) for (_, _, fw, fh) in faces)
    return int(len(faces)), float(largest / img_area)


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
    files = [p for p in SRC.iterdir() if p.is_file() and p.suffix.lower() in IMG_EXTS]
    print(f"scanning {len(files)} images in {SRC}", flush=True)

    cascade = cv2.CascadeClassifier(cv2.data.haarcascades + "haarcascade_frontalface_default.xml")
    people = []
    for i, p in enumerate(files, 1):
        if i % 50 == 0 or i == 1:
            print(f"faces {i}/{len(files)}", flush=True)
        n, ratio = detect_faces(p, cascade)
        if n >= 1 and ratio >= 0.012:
            people.append((p, n, ratio))

    people.sort(key=lambda x: (-x[2], -x[1], x[0].name))
    moved = []
    for p, n, ratio in people:
        dest = safe_move(p, KEEP)
        moved.append({"file": dest.name, "faces": n, "largest_ratio": round(ratio, 4)})

    report = {
        "scanned": len(files),
        "people_focused": len(moved),
        "remaining_in_landscape": len(files) - len(moved),
        "people": moved,
        "criteria": "faces>=1 and largest_face_area_ratio>=0.012",
    }
    REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({
        "scanned": report["scanned"],
        "people_focused": report["people_focused"],
        "remaining_in_landscape": report["remaining_in_landscape"],
        "top10": moved[:10],
    }, ensure_ascii=False, indent=2), flush=True)
    print("report:", REPORT, flush=True)


if __name__ == "__main__":
    main()
