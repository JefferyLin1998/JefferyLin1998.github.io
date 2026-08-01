(function () {
  "use strict";

  var STORAGE_KEY = "birthday-gift-progress-v9";

  var content = null;
  var stationIndex = 0;
  var quizIndex = 0;
  var judgeStreak = 0;
  var judgeUsedInSession = [];
  var judgeCurrentStmt = null;
  var puzzleSize = 3;
  var puzzleTiles = [];
  var puzzleMoves = 0;
  var puzzleSolved = false;
  var puzzleSelected = -1;
  var touchStartX = 0;
  var completedStations = {};
  var unlockedUpTo = 0;
  var activeStation = null;
  var trainRunTimer = null;
  var comicQueue = [];
  var comicIndex = 0;
  var comicCorrect = 0;
  var comicBusy = false;
  var comicShowIntroOnce = false;

  var el = {};

  function $(id) {
    return document.getElementById(id);
  }

  function initElements() {
    el.cover = $("screen-cover");
    el.journey = $("screen-journey");
    el.stationsContainer = $("rail-stations");
    el.dots = $("journey-dots");
    el.train = $("rail-train");
    el.lineFill = $("rail-line-fill");
    el.railMap = $("rail-map");
    el.railMapWrap = $("rail-map-wrap");
    el.focusIcon = $("focus-icon");
    el.focusCity = $("focus-city");
    el.focusTagline = $("focus-tagline");
    el.focusMemory = $("focus-memory");
    el.focusStationLabel = $("focus-station-label");
    el.journeyTitle = $("journey-title");
    el.journeySubtitle = $("journey-subtitle");
    el.journeyTrainNo = $("journey-train-no");
    el.btnNext = $("btn-next");
    el.btnPrev = $("btn-prev");
    el.btnEnter = $("btn-enter-station");
    el.journeyHint = $("journey-hint");
    el.herNameCover = $("her-name-cover");
    el.judgeChat = $("judge-chat");
    el.judgeProgress = $("judge-progress");
    el.judgeFeedback = $("judge-feedback");
    el.judgeIntro = $("judge-intro");
    el.judgeTitle = $("judge-title");
    el.questionText = $("question-text");
    el.options = $("options");
    el.hintMsg = $("hint-msg");
    el.quizProgress = $("quiz-progress");
    el.letterBody = $("letter-body");
    el.daysCount = $("days-count");
    el.finaleTitle = $("finale-title");
    el.finaleSubtitle = $("finale-subtitle");
    el.cake = $("finale-cake");
    el.memoryBody = $("memory-body");
    el.memoryIcon = $("memory-icon");
    el.memoryCity = $("memory-city");
    el.memoryTagline = $("memory-tagline");
    el.puzzleBoard = $("puzzle-board");
    el.puzzleStatus = $("puzzle-status");
    el.puzzleTitle = $("puzzle-title");
    el.puzzleIntro = $("puzzle-intro");
    el.puzzleMsg = $("puzzle-msg");
    el.puzzlePeek = $("puzzle-peek");
    el.puzzlePeekImg = $("puzzle-peek-img");
    el.proposalScene = $("proposal-scene");
    el.proposalStage = $("proposal-stage");
    el.silStage = $("sil-stage");
    el.silSceneMedia = $("sil-scene-media");
    el.silSky = $("sil-sky");
    el.silCelestial = $("sil-celestial");
    el.silStageBadge = $("silhouette-stage");
    el.ringStepPanel = $("ring-step-panel");
    el.ringStepLabel = $("ring-step-label");
    el.comicIntro = $("comic-intro");
    el.comicImage = $("comic-image");
    el.comicOptions = $("comic-options");
    el.comicFeedback = $("comic-feedback");
    el.comicProgress = $("comic-progress");
    el.comicGuide = $("sx-guide");
    el.comicGuideName = $("comic-guide-name");
    el.comicGuideAvatar = $("comic-guide-avatar");
  }

  function loadContent() {
    return fetch("./data/content.json")
      .then(function (res) {
        if (!res.ok) throw new Error("load failed");
        return res.json();
      })
      .then(function (data) {
        content = data;
        applyContent();
        buildRailMap();
      });
  }

  function applyContent() {
    el.herNameCover.textContent = content.herName;
    el.journeyTitle.textContent = content.journey.title;
    el.journeySubtitle.textContent = content.journey.subtitle;
    el.journeyTrainNo.textContent = content.journey.trainNo || "G1314";
    el.finaleTitle.textContent = content.finale.title;
    el.finaleSubtitle.textContent = content.finale.subtitle;

    if (content.togetherDate) {
      var start = new Date(content.togetherDate + "T00:00:00");
      var now = new Date();
      var days = Math.floor((now - start) / (1000 * 60 * 60 * 24));
      if (days >= 0) {
        el.daysCount.textContent = days;
        $("days-together").classList.remove("hidden");
      }
    }
  }

  function saveProgress() {
    try {
      localStorage.setItem(
        STORAGE_KEY,
        JSON.stringify({
          seenCover: true,
          stationIndex: stationIndex,
          quizIndex: quizIndex,
          judgeStreak: judgeStreak,
          completed: completedStations,
          unlockedUpTo: unlockedUpTo,
        })
      );
    } catch (e) {}
  }

  function loadProgress() {
    try {
      var raw = localStorage.getItem(STORAGE_KEY);
      if (raw) return JSON.parse(raw);
    } catch (e) {}
    return null;
  }

  function getStationOrdinalLabel(index) {
    var station = content.stations[index];
    var cn = ["一", "二", "三", "四", "五", "六", "七", "八", "九", "十"];
    if (station && station.id === "tobecontinued") {
      return "终点站 · 未完待续";
    }
    if (index === content.stations.length - 1) {
      return "终点站";
    }
    if (index < cn.length) {
      return "第" + cn[index] + "站";
    }
    return "第 " + (index + 1) + " 站";
  }

  function scrollStationIntoView(index, smooth) {
    var nodes = el.stationsContainer.querySelectorAll(".rail-station");
    var wrap = el.railMapWrap;
    if (!nodes[index] || !wrap) return;
    var node = nodes[index];
    var wrapRect = wrap.getBoundingClientRect();
    var nodeRect = node.getBoundingClientRect();
    var target =
      wrap.scrollLeft +
      (nodeRect.left + nodeRect.width / 2) -
      (wrapRect.left + wrapRect.width / 2);
    wrap.scrollTo({ left: target, behavior: smooth ? "smooth" : "auto" });
  }

  function markCompleted(stationId) {
    if (completedStations[stationId]) return false;
    completedStations[stationId] = true;

    var idx = findStationIndexById(stationId);
    if (idx > -1) {
      var count = content.stations.length;
      if (idx + 1 > unlockedUpTo) {
        unlockedUpTo = Math.min(idx + 1, count - 1);
      }
    }

    saveProgress();
    updateStationStates();
    return true;
  }

  function findStationIndexById(id) {
    for (var i = 0; i < content.stations.length; i++) {
      if (content.stations[i].id === id) return i;
    }
    return -1;
  }

  function isStationUnlocked(index) {
    return index <= unlockedUpTo;
  }

  function advanceToNextStation() {
    var count = content.stations.length;
    if (stationIndex < unlockedUpTo) {
      goToStation(stationIndex + 1, true);
      return;
    }
    if (stationIndex + 1 < count) {
      unlockedUpTo = Math.min(stationIndex + 1, count - 1);
      saveProgress();
      updateStationStates();
      goToStation(stationIndex + 1, true);
    }
  }

  function showCover() {
    el.cover.classList.add("active");
    el.journey.classList.remove("active");
  }

  function showJourney() {
    el.cover.classList.remove("active");
    el.journey.classList.add("active");
    goToStation(stationIndex, false);
    saveProgress();
    requestAnimationFrame(function () {
      scrollStationIntoView(stationIndex, false);
    });
  }

  function buildRailMap() {
    var stations = content.stations;
    var minWidth = Math.max(340, stations.length * 64 + 56);
    el.railMap.style.minWidth = minWidth + "px";

    el.stationsContainer.innerHTML = "";
    el.dots.innerHTML = "";

    stations.forEach(function (station, i) {
      var node = document.createElement("button");
      node.type = "button";
      node.className = "rail-station";
      node.dataset.index = String(i);
      node.innerHTML =
        '<span class="station-city">' +
        station.city +
        "</span>" +
        '<span class="station-dot-wrap">' +
        '<span class="station-dot"></span>' +
        '<span class="station-check hidden">✓</span>' +
        '<span class="station-lock">🔒</span>' +
        "</span>" +
        '<span class="station-icon">' +
        station.icon +
        "</span>";

      node.addEventListener("click", function () {
        if (i > unlockedUpTo) {
          flashLockedHint();
          return;
        }
        goToStation(i, true);
      });
      el.stationsContainer.appendChild(node);

      var dot = document.createElement("button");
      dot.type = "button";
      dot.className = "journey-dot";
      dot.setAttribute("aria-label", station.city);
      dot.addEventListener("click", function () {
        if (i > unlockedUpTo) {
          flashLockedHint();
          return;
        }
        goToStation(i, true);
      });
      el.dots.appendChild(dot);
    });

    updateStationStates();
    requestAnimationFrame(function () {
      positionTrain(false);
    });
  }

  function goToStation(index, animate) {
    var count = content.stations.length;
    var target = Math.max(0, Math.min(index, count - 1));
    if (target > unlockedUpTo) {
      target = unlockedUpTo;
    }
    stationIndex = target;
    updateStationFocus();
    positionTrain(animate);
    saveProgress();
    if (animate) {
      scrollStationIntoView(stationIndex, true);
    }
  }

  function positionTrain(animate) {
    var nodes = el.stationsContainer.querySelectorAll(".rail-station");
    var activeNode = nodes[stationIndex];
    if (!activeNode || !el.train) return;

    var mapRect = el.railMap.getBoundingClientRect();
    var dotEl = activeNode.querySelector(".station-dot");
    var dotRect = dotEl.getBoundingClientRect();
    var count = content.stations.length;
    var trainW = el.train.offsetWidth || 88;

    var left = dotRect.left - mapRect.left + dotRect.width / 2 - trainW / 2;
    var trackEl = el.railMap.querySelector(".hsr-route-track");
    var trackRect = trackEl ? trackEl.getBoundingClientRect() : mapRect;
    var top = trackRect.top - mapRect.top + trackRect.height / 2 - 28;

    el.train.classList.toggle("no-transition", !animate);
    el.train.style.left = left + "px";
    el.train.style.top = top + "px";
    el.train.style.transform = "";

    if (animate) {
      el.train.classList.add("is-running");
      el.railMap.classList.add("is-train-running");
      if (trainRunTimer) clearTimeout(trainRunTimer);
      trainRunTimer = setTimeout(function () {
        el.train.classList.remove("is-running");
        el.railMap.classList.remove("is-train-running");
      }, 750);
    }

    if (!animate) {
      requestAnimationFrame(function () {
        el.train.classList.remove("no-transition");
      });
    }

    var fillPct = count <= 1 ? 100 : (stationIndex / (count - 1)) * 100;
    if (el.lineFill) {
      el.lineFill.style.width = fillPct + "%";
    }
  }

  function updateStationFocus() {
    var station = content.stations[stationIndex];

    el.focusStationLabel.textContent = getStationOrdinalLabel(stationIndex);
    el.focusIcon.textContent = station.icon;
    el.focusCity.textContent = station.city;
    el.focusTagline.textContent = station.tagline;
    el.focusMemory.textContent = station.memory || "";

    var dots = el.dots.querySelectorAll(".journey-dot");
    dots.forEach(function (dot, i) {
      dot.classList.toggle("active", i === stationIndex);
    });

    var nodes = el.stationsContainer.querySelectorAll(".rail-station");
    nodes.forEach(function (node, i) {
      node.classList.toggle("is-active", i === stationIndex);
    });

    updateJourneyHint();
    updateNavButtons();
  }

  function updateJourneyHint() {
    if (!el.journeyHint) return;
    var count = content.stations.length;
    var station = content.stations[stationIndex];
    var isDone = !!completedStations[station.id];
    var isLast = stationIndex === count - 1;

    if (el.btnEnter) {
      if (isLast && isDone) {
        el.btnEnter.textContent = "再来一次 🎂";
      } else if (isDone) {
        el.btnEnter.textContent = "再回忆一下";
      } else if (isLast) {
        el.btnEnter.textContent = "进入终点站";
      } else {
        el.btnEnter.textContent = "进入这一站";
      }
    }

    if (isLast && isDone) {
      el.journeyHint.textContent = "🎉 全程通关！感谢你和我的这段旅程";
      return;
    }

    if (isDone) {
      el.journeyHint.textContent = "这一站已打卡 ✓ 进入下一站继续旅程";
    } else {
      el.journeyHint.textContent = "完成当前挑战，解锁下一站";
    }
  }

  function updateNavButtons() {
    if (!el.btnPrev || !el.btnNext) return;
    var count = content.stations.length;
    el.btnPrev.classList.toggle("is-disabled", stationIndex <= 0);
    el.btnNext.classList.toggle("is-disabled", stationIndex >= unlockedUpTo);
    el.btnNext.classList.toggle("hidden", stationIndex >= count - 1);
  }

  var lockedHintTimer = null;
  function flashLockedHint() {
    if (!el.journeyHint) return;
    var original = el.journeyHint.textContent;
    el.journeyHint.textContent = "🔒 下一站还没解锁，先完成当前这一站吧～";
    el.journeyHint.classList.add("locked-flash");
    if (lockedHintTimer) clearTimeout(lockedHintTimer);
    lockedHintTimer = setTimeout(function () {
      el.journeyHint.classList.remove("locked-flash");
      updateJourneyHint();
    }, 2000);
  }

  function updateStationStates() {
    if (!el.stationsContainer) return;
    content.stations.forEach(function (station, i) {
      var node = el.stationsContainer.querySelectorAll(".rail-station")[i];
      if (!node) return;
      var done = !!completedStations[station.id];
      node.classList.toggle("is-done", done);
      node.classList.toggle("is-locked", i > unlockedUpTo);

      var check = node.querySelector(".station-check");
      if (check) check.classList.toggle("hidden", !done);

      var lock = node.querySelector(".station-lock");
      if (lock) lock.classList.toggle("hidden", i <= unlockedUpTo);

      var iconSpan = node.querySelector(".station-icon");
      if (iconSpan) iconSpan.classList.toggle("hidden", i > unlockedUpTo);
    });
  }

  function enterCurrentStation() {
    if (stationIndex > unlockedUpTo) {
      flashLockedHint();
      return;
    }
    var station = content.stations[stationIndex];
    activeStation = station;
    openOverlay(station.type, station);
  }

  function openOverlay(type, station) {
    closeAllOverlays();
    if (station) activeStation = station;
    var overlay = $("overlay-" + type);
    if (!overlay) return;
    overlay.hidden = false;
    document.body.classList.add("overlay-open");

    if (type === "judge") {
      startJudgeGame();
    }

    if (type === "puzzle" && station) {
      var puzzleData = getPuzzleData();
      $("puzzle-title").textContent =
        puzzleData.title || station.icon + " " + station.city + "站 · 拼图";
      startPuzzle();
    }

    if (type === "proposal" && station) {
      startProposal(station);
    }

    if (type === "silhouette" && station) {
      startSilhouette(station);
    }

    if (type === "ringdesign" && station) {
      startRingDesign(station);
    }

    if (type === "comicmatch" && station) {
      startComicMatch(station);
    }

    if (type === "memory" && station) {
      $("overlay-memory-city").textContent = station.city;
      el.memoryIcon.textContent = station.icon;
      el.memoryCity.textContent = station.city;
      el.memoryTagline.textContent = station.tagline;
      el.memoryBody.textContent = station.memory || "";
    }

    if (type === "quiz") {
      if (station) {
        $("quiz-station-title").textContent = station.icon + " " + station.city + "站 · 回忆问答";
      }
      renderQuestion();
    }

    if (type === "letter") {
      if (station) {
        var letterBadge = $("overlay-letter-city");
        if (letterBadge) letterBadge.textContent = station.icon + " " + station.city;
      }
      startLetter();
    }

    if (type === "finale") {
      if (station) {
        var finaleBadge = $("overlay-finale-city");
        if (finaleBadge) finaleBadge.textContent = station.icon + " " + station.city;
      }
      spawnConfetti();
    }
  }

  function closeOverlay(type) {
    var overlay = $("overlay-" + type);
    if (overlay) overlay.hidden = true;
    activeStation = null;
    if (type === "ringdesign") {
      var resEl = $("ring-result-content");
      if (resEl) resEl.remove();
    }
    if (!document.querySelector(".project-overlay:not([hidden])")) {
      document.body.classList.remove("overlay-open");
    }
    requestAnimationFrame(function () {
      positionTrain(false);
    });
  }

  function closeAllOverlays() {
    document.querySelectorAll(".project-overlay").forEach(function (o) {
      o.hidden = true;
    });
    document.body.classList.remove("overlay-open");
    activeStation = null;
  }


  function getJudgeData() {
    return content.liuyangJudge || { questionBank: [], characters: [], passCount: 6 };
  }

  function getJudgeBank() {
    var data = getJudgeData();
    return data.questionBank || data.statements || [];
  }

  function getJudgePassCount() {
    var data = getJudgeData();
    return data.passCount || 6;
  }

  function pickRandomStatement() {
    var bank = getJudgeBank();
    if (!bank.length) return null;

    var available = [];
    for (var i = 0; i < bank.length; i++) {
      if (judgeUsedInSession.indexOf(i) === -1) available.push(i);
    }
    if (!available.length) {
      judgeUsedInSession = [];
      for (var j = 0; j < bank.length; j++) available.push(j);
    }

    var pick = available[Math.floor(Math.random() * available.length)];
    judgeUsedInSession.push(pick);
    return bank[pick];
  }

  function getStatementSide(stmt) {
    if (stmt && typeof stmt.speaker === "number") {
      return stmt.speaker === 0 ? "left" : "right";
    }
    return judgeStreak % 2 === 0 ? "left" : "right";
  }

  function appendJudgeBubble(stmt, side) {
    var data = getJudgeData();
    if (!stmt) return;
    var char = data.characters[stmt.speaker] || { name: "娃", avatar: "" };
    var isRight = side === "right";

    var row = document.createElement("div");
    row.className = "chat-row " + (isRight ? "chat-row-right" : "chat-row-left");

    var avatar = document.createElement("img");
    avatar.className = "chat-avatar";
    avatar.src = char.avatar;
    avatar.alt = char.name;
    avatar.loading = "lazy";

    var body = document.createElement("div");
    body.className = "chat-body";

    var name = document.createElement("span");
    name.className = "chat-name";
    name.textContent = char.name;

    var bubble = document.createElement("div");
    bubble.className = "chat-bubble";
    bubble.textContent = stmt.text;

    body.appendChild(name);
    body.appendChild(bubble);

    if (isRight) {
      row.appendChild(body);
      row.appendChild(avatar);
    } else {
      row.appendChild(avatar);
      row.appendChild(body);
    }

    el.judgeChat.appendChild(row);
  }

  function appendJudgeResult(isCorrect, customText) {
    var row = document.createElement("div");
    row.className = "chat-row chat-row-system";
    var bubble = document.createElement("div");
    bubble.className = "chat-system-msg" + (isCorrect ? "" : " chat-system-msg-fail");
    bubble.textContent =
      customText || (isCorrect ? "✓ 判断正确" : "✗ 答错了，从头再来");
    row.appendChild(bubble);
    el.judgeChat.appendChild(row);
  }

  function scrollJudgeChatToBottom() {
    if (!el.judgeChat) return;
    requestAnimationFrame(function () {
      el.judgeChat.scrollTop = el.judgeChat.scrollHeight;
    });
  }

  function updateJudgeProgress() {
    var passCount = getJudgePassCount();
    if (el.judgeProgress) {
      el.judgeProgress.textContent = "连对 " + judgeStreak + " / " + passCount;
    }
  }

  function setJudgeActionsEnabled(enabled) {
    $("btn-judge-true").disabled = !enabled;
    $("btn-judge-false").disabled = !enabled;
  }

  function showCurrentJudgeQuestion() {
    judgeCurrentStmt = pickRandomStatement();
    if (!judgeCurrentStmt) {
      el.judgeFeedback.textContent = "题库为空，请先配置题目～";
      setJudgeActionsEnabled(false);
      return;
    }
    appendJudgeBubble(judgeCurrentStmt, getStatementSide(judgeCurrentStmt));
    updateJudgeProgress();
    scrollJudgeChatToBottom();
  }

  function restartJudgeRound() {
    judgeStreak = 0;
    judgeUsedInSession = [];
    judgeCurrentStmt = null;
    saveProgress();
    el.judgeChat.innerHTML = "";
    el.judgeFeedback.textContent = "答错了，从头再来！加油～";
    showCurrentJudgeQuestion();
    setJudgeActionsEnabled(true);
  }

  function startJudgeGame() {
    var data = getJudgeData();
    if (el.judgeTitle) {
      el.judgeTitle.textContent = data.title || "浏阳站 · 真假大对决";
    }
    if (el.judgeIntro) {
      el.judgeIntro.textContent = data.intro || "";
    }
    el.judgeFeedback.textContent = "";
    el.judgeChat.innerHTML = "";
    setJudgeActionsEnabled(true);
    showCurrentJudgeQuestion();
  }

  function onJudgeAnswer(userSaysTrue) {
    var stmt = judgeCurrentStmt;
    if (!stmt) return;

    setJudgeActionsEnabled(false);
    var correct = userSaysTrue === stmt.isTrue;

    if (correct) {
      judgeStreak += 1;
      saveProgress();
      el.judgeFeedback.textContent = "答对啦！❤️";
      appendJudgeResult(true);
      scrollJudgeChatToBottom();
      updateJudgeProgress();

      setTimeout(function () {
        el.judgeFeedback.textContent = "";
        var passCount = getJudgePassCount();

        if (judgeStreak >= passCount) {
          var station = activeStation || content.stations.find(function (s) {
            return s.type === "judge";
          });
          if (station) markCompleted(station.id);
          judgeStreak = 0;
          judgeUsedInSession = [];
          judgeCurrentStmt = null;
          saveProgress();
          closeOverlay("judge");
          advanceToNextStation();
          return;
        }

        showCurrentJudgeQuestion();
        setJudgeActionsEnabled(true);
      }, 700);
    } else {
      el.judgeFeedback.textContent = stmt.wrongHint || "不对哦～";
      appendJudgeResult(false);
      scrollJudgeChatToBottom();
      setTimeout(restartJudgeRound, 1400);
    }
  }

  /* ── 北海站：交换拼图 ── */
  function getPuzzleData() {
    return content.beihaiPuzzle || { image: "", size: 3, title: "", intro: "", success: "" };
  }

  function isPuzzleSolved() {
    for (var i = 0; i < puzzleTiles.length; i++) {
      if (puzzleTiles[i] !== i) return false;
    }
    return true;
  }

  function puzzleSwap(i, j) {
    var tmp = puzzleTiles[i];
    puzzleTiles[i] = puzzleTiles[j];
    puzzleTiles[j] = tmp;
  }

  function shufflePuzzle() {
    var n = puzzleSize;
    var total = n * n;
    puzzleTiles = [];
    for (var i = 0; i < total; i++) puzzleTiles.push(i);
    puzzleSolved = false;
    puzzleMoves = 0;
    puzzleSelected = -1;

    do {
      for (var i = total - 1; i > 0; i--) {
        var j = Math.floor(Math.random() * (i + 1));
        puzzleSwap(i, j);
      }
    } while (isPuzzleSolved());
  }

  function renderPuzzle() {
    var data = getPuzzleData();
    var n = puzzleSize;
    var total = n * n;
    el.puzzleBoard.innerHTML = "";
    el.puzzleBoard.style.setProperty("--puzzle-n", n);

    for (var i = 0; i < total; i++) {
      var val = puzzleTiles[i];
      var cell = document.createElement("button");
      cell.type = "button";
      cell.className = "puzzle-cell";
      cell.dataset.pos = String(i);

      var tr = Math.floor(val / n);
      var tc = val % n;
      var posPct = n > 1 ? 100 / (n - 1) : 0;
      cell.style.backgroundImage = "url('" + data.image + "')";
      cell.style.backgroundSize = n * 100 + "% " + n * 100 + "%";
      cell.style.backgroundPosition = tc * posPct + "% " + tr * posPct + "%";

      var badge = document.createElement("span");
      badge.className = "puzzle-num";
      badge.textContent = val + 1;
      cell.appendChild(badge);

      if (i === puzzleSelected) {
        cell.classList.add("puzzle-cell-selected");
      }
      if (puzzleSolved) {
        cell.classList.add("puzzle-cell-done");
      }

      (function (idx) {
        cell.addEventListener("click", function () {
          onPuzzleTileClick(idx);
        });
      })(i);

      var r = Math.floor(i / n);
      var c = i % n;
      cell.style.left = (c / n) * 100 + "%";
      cell.style.top = (r / n) * 100 + "%";
      cell.style.width = 100 / n + "%";
      cell.style.height = 100 / n + "%";
      el.puzzleBoard.appendChild(cell);
    }

    if (el.puzzleStatus) {
      el.puzzleStatus.textContent = puzzleSolved
        ? puzzleMoves + " 次交换 · 通关"
        : puzzleMoves + " 次交换";
    }
  }

  function onPuzzleTileClick(index) {
    if (puzzleSolved) return;

    if (puzzleSelected === -1) {
      puzzleSelected = index;
      renderPuzzle();
      if (el.puzzleMsg) el.puzzleMsg.textContent = "再点另一块，就能交换啦～";
      return;
    }

    if (puzzleSelected === index) {
      puzzleSelected = -1;
      renderPuzzle();
      if (el.puzzleMsg) el.puzzleMsg.textContent = "";
      return;
    }

    puzzleSwap(puzzleSelected, index);
    puzzleMoves += 1;
    puzzleSelected = -1;
    renderPuzzle();

    if (isPuzzleSolved()) {
      puzzleSolved = true;
      var data = getPuzzleData();
      if (el.puzzleMsg) el.puzzleMsg.textContent = data.success || "拼好啦！❤️";
      var doneBtn = $("btn-puzzle-done");
      if (doneBtn) doneBtn.classList.remove("hidden");
      renderPuzzle();
    } else {
      if (el.puzzleMsg) el.puzzleMsg.textContent = "";
    }
  }

  function startPuzzle() {
    var data = getPuzzleData();
    puzzleSize = data.size || 3;
    if (el.puzzleIntro) el.puzzleIntro.textContent = data.intro || "";
    if (el.puzzleMsg) el.puzzleMsg.textContent = "";
    if (el.puzzlePeekImg) el.puzzlePeekImg.src = data.image;
    if (el.puzzlePeek) el.puzzlePeek.classList.add("hidden");
    var doneBtn = $("btn-puzzle-done");
    if (doneBtn) doneBtn.classList.add("hidden");
    shufflePuzzle();
    renderPuzzle();
  }

  /* ── 上海站：求婚仪式 ── */
  var proposalPhase = "intro";
  var proposalNoCount = 0;
  var proposalMatchLit = false;

  function getProposalData() {
    return content.shanghaiProposal || {};
  }

  function showProposalPhase(name) {
    proposalPhase = name;
    var phases = el.proposalScene.querySelectorAll(".proposal-phase");
    phases.forEach(function (p) {
      p.classList.toggle("active", p.id === "proposal-phase-" + name);
    });
    if (el.proposalStage) {
      var labels = {
        intro: "那一刻",
        match: "点亮",
        memory: "回忆",
        ask: "求婚",
        accept: "我愿意"
      };
      el.proposalStage.textContent = labels[name] || name;
    }
  }

  function daysSince(dateStr) {
    if (!dateStr) return 0;
    var start = new Date(dateStr + "T00:00:00");
    if (isNaN(start.getTime())) return 0;
    var now = new Date();
    var d = Math.floor((now - start) / 86400000);
    return d >= 0 ? d : 0;
  }

  function buildProposalStars() {
    var box = el.proposalScene.querySelector(".proposal-stars");
    if (!box || box.childElementCount) return;
    for (var i = 0; i < 30; i++) {
      var s = document.createElement("span");
      s.className = "proposal-star";
      s.textContent = "✦";
      s.style.left = Math.random() * 100 + "%";
      s.style.top = Math.random() * 60 + "%";
      s.style.fontSize = 6 + Math.random() * 10 + "px";
      s.style.animationDelay = Math.random() * 3 + "s";
      box.appendChild(s);
    }
  }

  function startProposal(station) {
    var data = getProposalData();
    buildProposalStars();

    $("proposal-intro-text").textContent = data.phases && data.phases.intro
      ? data.phases.intro
      : "有些瞬间，一辈子都不会忘。";
    $("proposal-match-tip").textContent = data.phases && data.phases.matchHint
      ? data.phases.matchHint
      : "用手指划过火柴，点燃蜡烛 →";
    $("proposal-memory-text").textContent = data.phases && data.phases.memory
      ? data.phases.memory
      : "";
    $("proposal-photo-caption").textContent = data.phases && data.phases.memoryCaption
      ? data.phases.memoryCaption
      : "";
    $("proposal-photo").src = data.image || "./pictures/求婚.jpg";
    $("proposal-question").textContent = data.phases && data.phases.ask
      ? data.phases.ask
      : "你愿意嫁给我吗？";
    $("btn-proposal-yes").textContent = data.phases && data.phases.yesText
      ? data.phases.yesText
      : "我愿意 ❤️";
    $("proposal-accept-text").textContent = data.phases && data.phases.accept
      ? data.phases.accept
      : "";

    proposalMatchLit = false;
    proposalNoCount = 0;
    bindMatchGesture();
    showProposalPhase("intro");
  }

  function gotoMatchPhase() {
    showProposalPhase("match");
  }

  function gotoMemoryPhase() {
    showProposalPhase("memory");
    var photo = $("proposal-photo");
    if (photo) {
      photo.classList.remove("is-revealed");
      void photo.offsetWidth;
      photo.classList.add("is-revealed");
    }
  }

  function gotoAskPhase() {
    showProposalPhase("ask");
    proposalNoCount = 0;
    positionProposalNoButton();
  }

  function gotoAcceptPhase() {
    showProposalPhase("accept");
    var data = getProposalData();
    var days = daysSince(data.proposalDate);
    var label = (data.phases && data.phases.daysLabel) || "从那一天到今天，我们已经携手走过";
    $("proposal-days").textContent = label + " " + days + " 天 ❤️";
    spawnProposalFireworks();
  }

  function bindMatchGesture() {
    var match = $("match");
    var candle = $("candle");
    if (!match || !candle || match.dataset.bound === "1") return;
    match.dataset.bound = "1";

    var lit = false;
    function light() {
      if (lit || proposalPhase !== "match") return;
      lit = true;
      proposalMatchLit = true;
      match.classList.add("is-struck");
      candle.classList.add("is-lit");
      setTimeout(gotoMemoryPhase, 900);
    }

    function bindDrag(target) {
      function handler(e) {
        e.preventDefault();
        light();
      }
      target.addEventListener("click", handler);
      target.addEventListener("touchstart", handler, { passive: false });
    }
    bindDrag(match);
    bindDrag(candle);
  }

  function positionProposalNoButton() {
    var no = $("btn-proposal-no");
    if (!no) return;
    no.style.position = "";
    no.style.left = "";
    no.style.top = "";
    no.style.transform = "";
  }

  function dodgeProposalNo() {
    var no = $("btn-proposal-no");
    if (!no) return;
    proposalNoCount += 1;
    var data = getProposalData();
    var teases = (data.phases && data.phases.noTease) || ["再想想嘛～"];
    no.textContent = teases[(proposalNoCount - 1) % teases.length];

    var sceneRect = el.proposalScene.getBoundingClientRect();
    var btnW = no.offsetWidth || 90;
    var btnH = no.offsetHeight || 44;
    var padding = 12;
    var maxX = sceneRect.width - btnW - padding;
    var maxY = sceneRect.height - btnH - padding;
    var x = padding + Math.random() * Math.max(0, maxX - padding);
    var y = padding + Math.random() * Math.max(0, maxY - padding);
    no.style.position = "absolute";
    no.style.left = x + "px";
    no.style.top = y + "px";
  }

  function spawnProposalFireworks() {
    var colors = ["#ffd54f", "#ffb300", "#ff8a65", "#f06292", "#81c784", "#fff"];
    var scene = el.proposalScene;
    for (var i = 0; i < 8; i++) {
      (function (k) {
        setTimeout(function () {
          var burst = document.createElement("div");
          burst.className = "proposal-firework";
          burst.style.left = 15 + Math.random() * 70 + "%";
          burst.style.top = 10 + Math.random() * 40 + "%";
          var color = colors[k % colors.length];
          for (var p = 0; p < 12; p++) {
            var particle = document.createElement("span");
            particle.className = "proposal-spark";
            var ang = (Math.PI * 2 * p) / 12;
            var dist = 30 + Math.random() * 20;
            particle.style.setProperty("--dx", Math.cos(ang) * dist + "px");
            particle.style.setProperty("--dy", Math.sin(ang) * dist + "px");
            particle.style.background = color;
            particle.style.boxShadow = "0 0 6px " + color;
            burst.appendChild(particle);
          }
          scene.appendChild(burst);
          setTimeout(function () {
            burst.remove();
          }, 1200);
        }, k * 250);
      })(i);
    }
  }

  /* ── 南京站：剪影叙事 ── */
  var silIndex = 0;

  function getSilData() {
    return content.nanjingSilhouette || { scenes: [] };
  }

  function setSilSky(index, total) {
    if (!el.silSky || !el.silCelestial) return;
    var ratio = total <= 1 ? 0 : index / (total - 1);
    var topPct = 8 + ratio * 55;
    var hue1 = 32 - ratio * 32;
    var hue2 = 210 + ratio * 30;
    var skyTop = "hsl(" + hue1 + ", 85%, " + (62 - ratio * 30) + "%)";
    var skyMid = "hsl(" + (hue1 - 8) + ", 70%, " + (48 - ratio * 28) + "%)";
    var skyBottom = "hsl(" + hue2 + ", 55%, " + (22 - ratio * 12) + "%)";
    el.silSky.style.background =
      "linear-gradient(180deg, " + skyTop + " 0%, " + skyMid + " 45%, " + skyBottom + " 100%)";

    if (ratio < 0.55) {
      el.silCelestial.className = "sil-celestial sil-celestial-sun";
      el.silCelestial.style.background =
        "radial-gradient(circle, #fff5c2 0%, #ffd54f 45%, rgba(255,138,101,0.5) 75%, transparent 100%)";
    } else {
      el.silCelestial.className = "sil-celestial sil-celestial-moon";
      el.silCelestial.style.background =
        "radial-gradient(circle, #fff 0%, #e0e7ff 50%, rgba(160,180,255,0.4) 80%, transparent 100%)";
    }
    el.silCelestial.style.top = topPct + "%";
    el.silCelestial.style.right = 14 + ratio * 30 + "%";
  }

  var SIL_ART = {
    military: function () {
      return ''
        + '<div class="sil-art sil-art-military">'
        +   '<div class="sil-ground"></div>'
        +   '<div class="sil-sun-disk"></div>'
        +   '<div class="sil-figure-group">'
        +     '<div class="sil-figure sil-figure-m1"></div>'
        +     '<div class="sil-figure sil-figure-m2"></div>'
        +     '<div class="sil-figure sil-figure-m3"></div>'
        +     '<div class="sil-figure sil-figure-m4"></div>'
        +     '<div class="sil-figure sil-figure-m5"></div>'
        +   '</div>'
        + '</div>';
    },
    dorm: function () {
      return ''
        + '<div class="sil-art sil-art-dorm">'
        +   '<div class="sil-room">'
        +     '<div class="sil-bunk sil-bunk-l"></div>'
        +     '<div class="sil-bunk sil-bunk-r"></div>'
        +     '<div class="sil-table"></div>'
        +     '<div class="sil-cake"></div>'
        +     '<div class="sil-figure sil-figure-sit sil-figure-sit1"></div>'
        +     '<div class="sil-figure sil-figure-sit sil-figure-sit2"></div>'
        +     '<div class="sil-figure sil-figure-sit sil-figure-sit3"></div>'
        +     '<div class="sil-flag"></div>'
        +   '</div>'
        + '</div>';
    },
    club: function () {
      return ''
        + '<div class="sil-art sil-art-club">'
        +   '<div class="sil-floor"></div>'
        +   '<div class="sil-desk"></div>'
        +   '<div class="sil-laptop"></div>'
        +   '<div class="sil-figure sil-figure-type sil-figure-type1"></div>'
        +   '<div class="sil-figure sil-figure-type sil-figure-type2"></div>'
        +   '<div class="sil-bubble"></div>'
        + '</div>';
    },
    library: function () {
      return ''
        + '<div class="sil-art sil-art-library">'
        +   '<div class="sil-bookshelf"></div>'
        +   '<div class="sil-lamp"></div>'
        +   '<div class="sil-lamp-light"></div>'
        +   '<div class="sil-figure sil-figure-sit sil-figure-study"></div>'
        +   '<div class="sil-formula">∑ ∫ lim</div>'
        + '</div>';
    },
    graduation: function () {
      return ''
        + '<div class="sil-art sil-art-graduation">'
        +   '<div class="sil-ground"></div>'
        +   '<div class="sil-figure sil-figure-grad sil-figure-grad1"></div>'
        +   '<div class="sil-figure sil-figure-grad sil-figure-grad2"></div>'
        +   '<div class="sil-figure sil-figure-grad sil-figure-grad3"></div>'
        +   '<div class="sil-figure sil-figure-grad sil-figure-grad4"></div>'
        +   '<div class="sil-cap"></div>'
        + '</div>';
    }
  };

  function renderSilScene(index) {
    var data = getSilData();
    var scene = data.scenes[index];
    if (!scene) return;
    var total = data.scenes.length;

    setSilSky(index, total);

    $("sil-year").textContent = scene.year;
    $("sil-scene-title").textContent = scene.title;
    $("sil-scene-location").textContent = "📍 " + scene.location;

    if (scene.photo) {
      el.silSceneMedia.innerHTML =
        '<div class="sil-photo-wrap"><img class="sil-photo" src="' + scene.photo +
        '" alt="' + scene.title + '" loading="lazy" /><div class="sil-photo-polaroid-tape"></div></div>';
      el.silSceneMedia.className = "sil-scene-media sil-scene-media-photo";
    } else {
      var art = (SIL_ART[scene.scene] || SIL_ART.military)();
      el.silSceneMedia.innerHTML = art;
      el.silSceneMedia.className = "sil-scene-media sil-scene-media-art";
      el.silSceneMedia.dataset.scene = scene.scene || "";
    }

    var narration = $("sil-narration");
    narration.textContent = "";
    el.silSceneMedia.style.opacity = "0";
    void el.silSceneMedia.offsetWidth;
    el.silSceneMedia.style.opacity = "1";

    if (el.silStageBadge) {
      el.silStageBadge.textContent = "第 " + (index + 1) + " 幕 / " + total;
    }

    typeText(narration, scene.narration, 28);

    var prog = $("sil-progress");
    prog.innerHTML = "";
    for (var i = 0; i < total; i++) {
      var dot = document.createElement("span");
      dot.className = "sil-dot" + (i === index ? " active" : "") + (i < index ? " done" : "");
      prog.appendChild(dot);
    }

    $("btn-sil-prev").disabled = index === 0;
    $("btn-sil-next").textContent = index === total - 1 ? "走完这四年" : "下一幕";
  }

  function typeText(node, text, speed) {
    if (!node) return;
    node.textContent = "";
    var i = 0;
    function step() {
      if (i <= text.length) {
        node.textContent = text.slice(0, i);
        i += 1;
        setTimeout(step, speed);
      }
    }
    step();
  }

  function startSilhouette(station) {
    var data = getSilData();
    silIndex = 0;
    $("sil-title").textContent = data.title || "我的大学四年";
    $("sil-subtitle").textContent = data.subtitle || "";
    $("sil-intro-text").textContent = data.intro || "";

    var avatar = $("sil-avatar");
    var ph = avatar.querySelector(".sil-avatar-placeholder");
    if (data.avatar) {
      var img = document.createElement("img");
      img.src = data.avatar;
      img.alt = "我";
      img.onerror = function () {
        img.style.display = "none";
        if (ph) ph.style.display = "block";
      };
      avatar.appendChild(img);
      if (ph) ph.style.display = "none";
    }

    ["sil-intro", "sil-scene-phase", "sil-outro"].forEach(function (id) {
      $(id).classList.toggle("active", id === "sil-intro");
    });

    buildSilStars();
  }

  function buildSilStars() {
    var box = el.silStage.querySelector(".sil-stars-bg");
    if (!box || box.childElementCount) return;
    for (var i = 0; i < 40; i++) {
      var s = document.createElement("span");
      s.className = "sil-star";
      s.textContent = "✦";
      s.style.left = Math.random() * 100 + "%";
      s.style.top = Math.random() * 70 + "%";
      s.style.fontSize = 5 + Math.random() * 9 + "px";
      s.style.animationDelay = Math.random() * 3 + "s";
      box.appendChild(s);
    }
  }

  function silNext() {
    var data = getSilData();
    if (silIndex < data.scenes.length - 1) {
      silIndex += 1;
      renderSilScene(silIndex);
    } else {
      $("sil-outro-text").textContent = data.outro || "";
      ["sil-intro", "sil-scene-phase", "sil-outro"].forEach(function (id) {
        $(id).classList.toggle("active", id === "sil-outro");
      });
    }
  }

  function silPrev() {
    if (silIndex > 0) {
      silIndex -= 1;
      renderSilScene(silIndex);
    }
  }

  /* ── 三亚站：选钻小课堂 ── */
  var ringParamIndex = 0;
  var ringPhase = "param";
  var ringChoices = {};

  function getRingData() {
    return content.sanyaRing || { paramSteps: [], realRing: {} };
  }

  function getRingParamSteps() {
    return getRingData().paramSteps || [];
  }

  function getRingOption(stepId, optionId) {
    var step = getRingParamSteps().find(function (s) {
      return s.id === stepId;
    });
    if (!step) return null;
    return (step.options || []).find(function (o) {
      return o.id === optionId;
    }) || null;
  }

  function escapeHtml(str) {
    return String(str || "")
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

  function getRingPreviewState() {
    return {
      carat: ringChoices.carat || "50",
      color: ringChoices.color || "gh",
      clarity: ringChoices.clarity || "vs",
      cut: ringChoices.cut || "3ex",
      metal: ringChoices.metal || "platinum",
      style: ringChoices.style || "solitaire",
      locked: {
        carat: !!ringChoices.carat,
        color: !!ringChoices.color,
        clarity: !!ringChoices.clarity,
        cut: !!ringChoices.cut,
        metal: !!ringChoices.metal,
        style: !!ringChoices.style
      }
    };
  }

  function ringMetalPalette(metalId) {
    if (metalId === "rosegold") {
      return { band: "#d4a08c", light: "#f6ddd0", dark: "#8f5c4c", mid: "#e0b8a8" };
    }
    if (metalId === "white18k") {
      return { band: "#e5e8ec", light: "#ffffff", dark: "#8992a0", mid: "#cfd4da" };
    }
    return { band: "#c8d0d8", light: "#f8fafc", dark: "#6d7986", mid: "#aeb8c2" };
  }

  function ringDiamondPalette(colorId) {
    if (colorId === "df") {
      return {
        body: "#e9f2fb",
        glow: "#ffffff",
        edge: "#88a4c2",
        facetLight: "#ffffff",
        facetDark: "#cfe1f2",
        tableCore: "#ffffff"
      };
    }
    if (colorId === "ij") {
      return {
        body: "#f5ead1",
        glow: "#fffaf0",
        edge: "#bd9a63",
        facetLight: "#fffaf0",
        facetDark: "#ecd6a8",
        tableCore: "#fffdf6"
      };
    }
    return {
      body: "#edf0f3",
      glow: "#ffffff",
      edge: "#8f9ca9",
      facetLight: "#ffffff",
      facetDark: "#dbe1e7",
      tableCore: "#ffffff"
    };
  }

  function ringCaratRadius(caratId) {
    if (caratId === "30") return 11;
    if (caratId === "100") return 20;
    return 15;
  }

  function ringCutSparkle(cutId) {
    if (cutId === "good") return { rays: 4, opacity: 0.28, facetBoost: 0.5 };
    if (cutId === "vg") return { rays: 6, opacity: 0.5, facetBoost: 0.78 };
    return { rays: 8, opacity: 0.82, facetBoost: 1 };
  }

  function ringPolarPt(cx, cy, r, deg) {
    var rad = (deg * Math.PI) / 180;
    return [
      (cx + r * Math.cos(rad)).toFixed(1),
      (cy + r * Math.sin(rad)).toFixed(1)
    ];
  }

  function ringClarityMarks(clarityId, cx, cy, r) {
    if (clarityId === "vvs") return "";
    if (clarityId === "vs") {
      return (
        '<circle cx="' +
        (cx + r * 0.22) +
        '" cy="' +
        (cy - r * 0.08) +
        '" r="0.8" fill="rgba(70,80,90,0.2)"/>'
      );
    }
    return (
      '<circle cx="' +
      (cx + r * 0.26) +
      '" cy="' +
      (cy - r * 0.1) +
      '" r="1.1" fill="rgba(60,70,80,0.3)"/>' +
      '<circle cx="' +
      (cx - r * 0.2) +
      '" cy="' +
      (cy + r * 0.24) +
      '" r="0.7" fill="rgba(60,70,80,0.24)"/>'
    );
  }

  function ringSparkleRays(cx, cy, r, cut) {
    var html = "";
    var n = cut.rays;
    for (var i = 0; i < n; i++) {
      var a = (Math.PI * 2 * i) / n - Math.PI / 2;
      var x1 = cx + Math.cos(a) * (r + 3);
      var y1 = cy + Math.sin(a) * (r + 3);
      var x2 = cx + Math.cos(a) * (r + 10 + (i % 2) * 3);
      var y2 = cy + Math.sin(a) * (r + 10 + (i % 2) * 3);
      html +=
        '<line x1="' +
        x1.toFixed(1) +
        '" y1="' +
        y1.toFixed(1) +
        '" x2="' +
        x2.toFixed(1) +
        '" y2="' +
        y2.toFixed(1) +
        '" stroke="rgba(255,255,255,' +
        cut.opacity +
        ')" stroke-width="1.2" stroke-linecap="round"/>';
    }
    return html;
  }

  /**
   * Round brilliant top-view: octagon table + 8 kite/star facets radiating
   * to the girdle, plus specular flashes tuned by cut quality.
   */
  function ringFacetedGem(cx, cy, r, gem, cut, clarityId, uid) {
    var tableR = r * 0.44;
    var angles = [];
    var i;
    for (i = 0; i < 8; i++) angles.push(-90 + i * 45);

    function tablePt(idx) {
      return ringPolarPt(cx, cy, tableR, angles[((idx % 8) + 8) % 8]);
    }
    function girdlePt(idx) {
      return ringPolarPt(cx, cy, r, angles[((idx % 8) + 8) % 8] + 22.5);
    }

    var html =
      '<circle cx="' +
      cx +
      '" cy="' +
      cy +
      '" r="' +
      r +
      '" fill="url(#' +
      uid +
      'gem)" stroke="' +
      gem.edge +
      '" stroke-width="1"/>';

    for (i = 0; i < 8; i++) {
      var p1 = tablePt(i);
      var p2 = tablePt(i + 1);
      var p3 = girdlePt(i);
      var shade = i % 2 === 0 ? gem.facetLight : gem.facetDark;
      html +=
        '<polygon points="' +
        p1.join(",") +
        " " +
        p2.join(",") +
        " " +
        p3.join(",") +
        '" fill="' +
        shade +
        '" stroke="' +
        gem.edge +
        '" stroke-width="0.4" opacity="0.92"/>';
    }

    var tablePts = [];
    for (i = 0; i < 8; i++) tablePts.push(tablePt(i).join(","));
    html +=
      '<polygon points="' +
      tablePts.join(" ") +
      '" fill="url(#' +
      uid +
      'table)" stroke="' +
      gem.edge +
      '" stroke-width="0.5"/>';

    var flashA = [tablePt(0), tablePt(1), girdlePt(0)]
      .map(function (p) {
        return p.join(",");
      })
      .join(" ");
    var flashB = [tablePt(4), tablePt(5), girdlePt(4)]
      .map(function (p) {
        return p.join(",");
      })
      .join(" ");
    html +=
      '<polygon points="' +
      flashA +
      '" fill="#ffffff" opacity="' +
      (0.55 * cut.facetBoost).toFixed(2) +
      '"/>';
    html +=
      '<polygon points="' +
      flashB +
      '" fill="#ffffff" opacity="' +
      (0.3 * cut.facetBoost).toFixed(2) +
      '"/>';

    var glint = ringPolarPt(cx, cy, tableR * 0.4, -55);
    html +=
      '<circle cx="' +
      glint[0] +
      '" cy="' +
      glint[1] +
      '" r="' +
      Math.max(1, r * 0.1).toFixed(1) +
      '" fill="#ffffff" opacity="' +
      cut.opacity +
      '"/>';

    html += ringClarityMarks(clarityId, cx, cy, r);
    html += ringSparkleRays(cx, cy, r, cut);
    return html;
  }

  /** Small round stones (pave etc.) — gradient body + a single catch-light dot. */
  function ringSimpleGem(cx, cy, r, gem, uid) {
    var hl = ringPolarPt(cx, cy, r * 0.32, -50);
    return (
      '<circle cx="' +
      cx.toFixed(1) +
      '" cy="' +
      cy.toFixed(1) +
      '" r="' +
      r.toFixed(1) +
      '" fill="url(#' +
      uid +
      'gem)" stroke="' +
      gem.edge +
      '" stroke-width="0.5"/>' +
      '<circle cx="' +
      hl[0] +
      '" cy="' +
      hl[1] +
      '" r="' +
      Math.max(0.6, r * 0.32).toFixed(1) +
      '" fill="#ffffff" opacity="0.85"/>'
    );
  }

  function ringSvgMarkup(state) {
    var s = state || getRingPreviewState();
    var metal = ringMetalPalette(s.metal);
    var gem = ringDiamondPalette(s.color);
    var r = ringCaratRadius(s.carat);
    var cut = ringCutSparkle(s.cut);
    var style = s.style || "solitaire";
    var uid = "rg" + Math.floor(Math.random() * 1e6);

    var defs =
      "<defs>" +
      '<linearGradient id="' +
      uid +
      'band" x1="0" y1="0" x2="0" y2="1">' +
      '<stop offset="0%" stop-color="' +
      metal.light +
      '"/>' +
      '<stop offset="45%" stop-color="' +
      metal.band +
      '"/>' +
      '<stop offset="100%" stop-color="' +
      metal.dark +
      '"/>' +
      "</linearGradient>" +
      '<radialGradient id="' +
      uid +
      'gem" cx="35%" cy="30%" r="72%">' +
      '<stop offset="0%" stop-color="#fff"/>' +
      '<stop offset="38%" stop-color="' +
      gem.glow +
      '"/>' +
      '<stop offset="100%" stop-color="' +
      gem.body +
      '"/>' +
      "</radialGradient>" +
      '<radialGradient id="' +
      uid +
      'table" cx="50%" cy="35%" r="75%">' +
      '<stop offset="0%" stop-color="#ffffff"/>' +
      '<stop offset="55%" stop-color="' +
      gem.tableCore +
      '"/>' +
      '<stop offset="100%" stop-color="' +
      gem.facetLight +
      '"/>' +
      "</radialGradient>" +
      '<filter id="' +
      uid +
      'soft"><feGaussianBlur stdDeviation="2"/></filter>' +
      "</defs>";

    var band =
      '<ellipse cx="120" cy="118" rx="78" ry="26" fill="none" stroke="url(#' +
      uid +
      'band)" stroke-width="14"/>' +
      '<ellipse cx="120" cy="118" rx="74" ry="22.5" fill="none" stroke="' +
      metal.dark +
      '" stroke-width="1" opacity="0.4"/>' +
      '<path d="M 66 115 Q 120 90 174 115" fill="none" stroke="#ffffff" stroke-width="2.2" opacity="0.5" stroke-linecap="round"/>';

    var body = "";

    if (style === "pave") {
      var stones = "";
      for (var i = 0; i < 11; i++) {
        var t = i / 10;
        var px = 48 + t * 144;
        var py = 104 + Math.abs(t - 0.5) * 28;
        var pr = 3.4 + (1 - Math.abs(t - 0.5) * 1.4) * (r / 16);
        stones +=
          '<circle cx="' +
          px.toFixed(1) +
          '" cy="' +
          py.toFixed(1) +
          '" r="' +
          (pr + 1.1).toFixed(1) +
          '" fill="none" stroke="' +
          metal.dark +
          '" stroke-width="1.2" opacity="0.7"/>' +
          ringSimpleGem(px, py, pr, gem, uid);
      }
      body =
        '<path d="M42 112 Q120 78 198 112" fill="none" stroke="url(#' +
        uid +
        'band)" stroke-width="13" stroke-linecap="round"/>' +
        '<path d="M48 112 Q120 86 192 112" fill="none" stroke="' +
        metal.light +
        '" stroke-width="2" opacity="0.5"/>' +
        stones +
        ringSparkleRays(120, 88, r * 0.55, cut);
    } else if (style === "geometric") {
      var accentR = Math.max(4.5, r * 0.4);
      body =
        '<rect x="40" y="104" width="160" height="16" rx="3" fill="url(#' +
        uid +
        'band)"/>' +
        '<path d="M96 104 L120 78 L144 104" fill="none" stroke="url(#' +
        uid +
        'band)" stroke-width="10" stroke-linejoin="round"/>' +
        '<path d="M104 104 L120 86 L136 104" fill="none" stroke="' +
        metal.light +
        '" stroke-width="2" opacity="0.65"/>' +
        '<circle cx="120" cy="72" r="' +
        (accentR + 2.2).toFixed(1) +
        '" fill="' +
        gem.glow +
        '" opacity="0.3" filter="url(#' +
        uid +
        'soft)"/>' +
        ringSimpleGem(120, 72, accentR, gem, uid) +
        ringSparkleRays(120, 72, accentR + 1, cut);
    } else {
      var prongTips = "";
      var prongLines = "";
      for (var p = 0; p < 6; p++) {
        var ang = (Math.PI * 2 * p) / 6 - Math.PI / 2;
        var px0 = 120 + Math.cos(ang) * (r * 0.5);
        var py0 = 78 + Math.sin(ang) * (r * 0.5);
        var px1 = 120 + Math.cos(ang) * (r + 4.5);
        var py1 = 78 + Math.sin(ang) * (r + 4.5);
        prongLines +=
          '<line x1="' +
          px0.toFixed(1) +
          '" y1="' +
          py0.toFixed(1) +
          '" x2="' +
          px1.toFixed(1) +
          '" y2="' +
          py1.toFixed(1) +
          '" stroke="' +
          metal.dark +
          '" stroke-width="2.4" stroke-linecap="round"/>';
        prongTips +=
          '<circle cx="' +
          px1.toFixed(1) +
          '" cy="' +
          py1.toFixed(1) +
          '" r="1.1" fill="' +
          metal.light +
          '"/>';
      }
      body =
        band +
        '<circle cx="120" cy="78" r="' +
        (r + 3) +
        '" fill="' +
        gem.glow +
        '" opacity="0.4" filter="url(#' +
        uid +
        'soft)"/>' +
        prongLines +
        ringFacetedGem(120, 78, r, gem, cut, s.clarity, uid) +
        prongTips;
    }

    return (
      '<svg viewBox="0 0 240 170" xmlns="http://www.w3.org/2000/svg">' +
      defs +
      '<ellipse cx="120" cy="148" rx="72" ry="9" fill="rgba(0,0,0,0.45)"/>' +
      '<ellipse cx="120" cy="146" rx="48" ry="5" fill="rgba(0,0,0,0.35)"/>' +
      body +
      "</svg>"
    );
  }

  function ringPreviewCaptionText(state) {
    var s = state || getRingPreviewState();
    var parts = [];
    if (s.locked.style) {
      var st = getRingOption("style", s.style);
      if (st) parts.push(st.name);
    }
    if (s.locked.metal) {
      var mt = getRingOption("metal", s.metal);
      if (mt) parts.push(mt.name);
    }
    if (s.locked.carat) {
      var ct = getRingOption("carat", s.carat);
      if (ct) parts.push(ct.name);
    }
    if (parts.length) return "预览 · " + parts.join(" · ");
    return "选择参数，预览会跟着变";
  }

  function formatRingBudget(amount) {
    var n = Math.round(Number(amount) || 0);
    return "¥" + n.toLocaleString("zh-CN");
  }

  function formatRingBudgetDelta(amount) {
    var n = Math.round(Number(amount) || 0);
    if (n <= 0) return "约 +¥0";
    return "约 +" + formatRingBudget(n);
  }

  function getRingBudgetSummary() {
    var steps = getRingParamSteps();
    var total = 0;
    var filled = 0;
    steps.forEach(function (step) {
      var opt = getRingOption(step.id, ringChoices[step.id]);
      if (!opt) return;
      filled += 1;
      total += Number(opt.budget) || 0;
    });
    return { total: total, filled: filled, count: steps.length };
  }

  function updateRingBudget() {
    var data = getRingData();
    var summary = getRingBudgetSummary();
    var labelEl = $("ring-budget-label");
    var valueEl = $("ring-budget-value");
    var hintEl = $("ring-budget-hint");
    var box = $("ring-budget");
    if (labelEl) labelEl.textContent = data.budgetLabel || "预估预算";
    if (!valueEl || !hintEl) return;

    if (summary.filled === 0) {
      valueEl.textContent = "—";
      hintEl.textContent = data.budgetEmpty || "选参数后同步估算";
      if (box) box.classList.remove("is-active");
      return;
    }

    valueEl.textContent = formatRingBudget(summary.total);
    if (summary.filled < summary.count) {
      hintEl.textContent = (
        data.budgetPartialHint ||
        "已计入 {n}/{total} 项，未选项暂不计入"
      )
        .replace("{n}", String(summary.filled))
        .replace("{total}", String(summary.count));
    } else {
      hintEl.textContent = data.budgetHint || "示意价 · 非门店报价";
    }
    if (box) box.classList.add("is-active");
  }

  function updateRingPreview() {
    var box = $("ring-preview-svg");
    var caption = $("ring-preview-caption");
    var state = getRingPreviewState();
    if (box) {
      box.innerHTML = ringSvgMarkup(state);
      box.classList.toggle("is-dim", !ringChoices.style && !ringChoices.metal);
    }
    if (caption && ringPhase !== "done") {
      caption.textContent = ringPreviewCaptionText(state);
    }
    updateRingBudget();
  }

  function updateRingSheet() {
    var data = getRingData();
    var steps = getRingParamSteps();
    var sheet = $("ring-sheet");
    var progress = $("ring-progress");
    var caption = $("ring-preview-caption");
    var titleEl = $("ring-sheet-title");
    if (titleEl) titleEl.textContent = data.sheetTitle || "你的配置单";
    if (!sheet) return;

    sheet.innerHTML = "";
    var filled = 0;
    steps.forEach(function (step, i) {
      var opt = getRingOption(step.id, ringChoices[step.id]);
      var li = document.createElement("li");
      li.className =
        "ring-sheet-item" +
        (opt ? " is-filled" : "") +
        (ringPhase === "param" && i === ringParamIndex ? " is-current" : "");
      li.innerHTML =
        '<span class="ring-sheet-key">' +
        escapeHtml(step.short || step.title) +
        "</span>" +
        '<span class="ring-sheet-val">' +
        escapeHtml(opt ? opt.name : data.sheetEmpty || "待选") +
        "</span>";
      sheet.appendChild(li);
      if (opt) filled += 1;
    });

    if (progress) {
      progress.innerHTML = "";
      steps.forEach(function (_, i) {
        var dot = document.createElement("span");
        dot.className =
          "ring-progress-dot" +
          (i < filled || (ringPhase === "done" && i < steps.length)
            ? " is-done"
            : "") +
          (ringPhase === "param" && i === ringParamIndex ? " is-current" : "");
        progress.appendChild(dot);
      });
    }

    if (caption && ringPhase === "done") {
      caption.textContent = "配置单已齐 · 这是按你参数生成的预览";
    }

    updateRingPreview();
  }

  function buildRingParamOptions() {
    var steps = getRingParamSteps();
    var step = steps[ringParamIndex];
    var box = $("ring-param-opts");
    var tipEl = $("ring-tip");
    var hintEl = $("ring-param-hint");
    if (!box || !step) return;

    if (hintEl) hintEl.textContent = step.hint || "";
    box.innerHTML = "";

    var selectedId = ringChoices[step.id] || null;
    (step.options || []).forEach(function (opt) {
      var btn = document.createElement("button");
      btn.type = "button";
      btn.className =
        "ring-choice" + (selectedId === opt.id ? " is-selected" : "");
      btn.innerHTML =
        '<span class="ring-choice-mark" aria-hidden="true"></span>' +
        '<span class="ring-choice-text"><strong>' +
        escapeHtml(opt.name) +
        "</strong><span>" +
        escapeHtml(opt.desc || "") +
        '</span><em class="ring-choice-budget">' +
        formatRingBudgetDelta(opt.budget) +
        "</em></span>";
      btn.addEventListener("click", function () {
        ringChoices[step.id] = opt.id;
        buildRingParamOptions();
        updateRingSheet();
      });
      box.appendChild(btn);
    });

    if (tipEl) {
      var selected = getRingOption(step.id, selectedId);
      if (selected && selected.tip) {
        tipEl.textContent = "小课堂 · " + selected.tip;
        tipEl.classList.remove("hidden");
      } else {
        tipEl.textContent = "";
        tipEl.classList.add("hidden");
      }
    }

    updateRingPreview();
  }

  function buildRingSummary() {
    var data = getRingData();
    var steps = getRingParamSteps();
    var list = $("ring-summary-list");
    var titleEl = $("ring-summary-title");
    if (titleEl) titleEl.textContent = data.summaryTitle || "配置完成";
    if (!list) return;
    list.innerHTML = "";
    steps.forEach(function (step) {
      var opt = getRingOption(step.id, ringChoices[step.id]);
      var li = document.createElement("li");
      li.innerHTML =
        "<strong>" +
        escapeHtml(step.short || step.title) +
        "</strong><span>" +
        escapeHtml(opt ? opt.name : "—") +
        "</span>";
      list.appendChild(li);
    });
  }

  function showRingPhase(phase) {
    ringPhase = phase;
    var data = getRingData();
    var steps = getRingParamSteps();

    if (el.ringStepPanel) {
      el.ringStepPanel.querySelectorAll(".ring-step").forEach(function (s) {
        s.classList.toggle("active", s.getAttribute("data-step") === phase);
      });
    }

    var introEl = $("ring-intro");
    if (introEl) {
      introEl.style.display =
        phase === "param" && ringParamIndex === 0 ? "" : "none";
    }

    if (phase === "param") {
      var step = steps[ringParamIndex];
      if (el.ringStepLabel && step) {
        el.ringStepLabel.textContent =
          ringParamIndex + 1 + " / " + steps.length + " · " + step.title;
      }
      buildRingParamOptions();
    } else {
      if (el.ringStepLabel) {
        el.ringStepLabel.textContent = "完成 · 配置单";
      }
      buildRingSummary();
      var teases = data.teases || ["不错。"];
      $("ring-tease").textContent =
        teases[Math.floor(Math.random() * teases.length)];
      $("ring-finalize-btn").textContent = data.finalCta || "揭晓我们的对戒";
    }

    updateRingSheet();
    $("ring-prev-btn").disabled = phase === "param" && ringParamIndex <= 0;
    $("ring-next-btn").style.display = phase === "done" ? "none" : "";
  }

  function ringCanAdvance() {
    if (ringPhase !== "param") return true;
    var step = getRingParamSteps()[ringParamIndex];
    return !!(step && ringChoices[step.id]);
  }

  function ringNext() {
    var steps = getRingParamSteps();
    if (!ringCanAdvance()) {
      if (el.ringStepLabel) {
        var step = steps[ringParamIndex];
        el.ringStepLabel.textContent =
          "请先选择" + ((step && step.title) || "一项");
      }
      return;
    }
    if (ringParamIndex < steps.length - 1) {
      ringParamIndex += 1;
      showRingPhase("param");
      return;
    }
    showRingPhase("done");
  }

  function ringPrev() {
    if (ringPhase === "done") {
      ringParamIndex = Math.max(0, getRingParamSteps().length - 1);
      showRingPhase("param");
      return;
    }
    if (ringParamIndex > 0) {
      ringParamIndex -= 1;
      showRingPhase("param");
    }
  }

  function finalizeRing() {
    var data = getRingData();
    showRingResult(data);
  }

  function showRingResult(data) {
    var overlay = $("overlay-ringdesign");
    overlay.hidden = false;
    document.body.classList.add("overlay-open");
    activeStation =
      activeStation ||
      content.stations.find(function (s) {
        return s.type === "ringdesign";
      });

    var real = data.realRing || {};
    var hers = real.hers || {};
    var his = real.his || {};
    var steps = getRingParamSteps();
    var configHtml = steps
      .map(function (step) {
        var opt = getRingOption(step.id, ringChoices[step.id]);
        return (
          "<li><strong>" +
          escapeHtml(step.short || step.title) +
          "</strong><span>" +
          escapeHtml(opt ? opt.name : "—") +
          "</span></li>"
        );
      })
      .join("");

    var html = "";
    html += '<div class="ring-result-overlay" id="ring-result-content">';
    html += '  <div class="ring-result-card">';
    html +=
      '    <h2 class="ring-result-title">' +
      escapeHtml(data.compareTitle || "对比") +
      "</h2>";
    html += '    <div class="ring-result-compare">';
    html += '      <div class="ring-result-col ring-result-yours">';
    html +=
      '        <p class="ring-result-label">' +
      escapeHtml(data.yourConfigLabel || "你刚才选的") +
      "</p>";
    html +=
      '        <div class="ring-result-preview">' +
      ringSvgMarkup(getRingPreviewState()) +
      "</div>";
    html +=
      '        <p class="ring-result-budget">' +
      escapeHtml(data.budgetLabel || "预估预算") +
      " · " +
      formatRingBudget(getRingBudgetSummary().total) +
      "</p>";
    html += '        <ul class="ring-result-sheet">' + configHtml + "</ul>";
    html += "      </div>";
    html += '      <div class="ring-result-vs">VS</div>';
    html += '      <div class="ring-result-col ring-result-real">';
    html +=
      '        <p class="ring-result-label">' +
      escapeHtml(real.label || "真实对戒") +
      "</p>";
    html += '        <div class="ring-result-pair">';
    html += '          <figure class="ring-result-photo">';
    if (hers.image) {
      html +=
        '            <img src="' +
        escapeHtml(hers.image) +
        '" alt="' +
        escapeHtml(hers.label || "女戒") +
        '" loading="lazy" />';
    }
    html +=
      "            <figcaption>" +
      escapeHtml((hers.label || "女戒") + " · " + (hers.name || "")) +
      "</figcaption>";
    html += "          </figure>";
    html += '          <figure class="ring-result-photo">';
    if (his.image) {
      html +=
        '            <img src="' +
        escapeHtml(his.image) +
        '" alt="' +
        escapeHtml(his.label || "男戒") +
        '" loading="lazy" />';
    }
    html +=
      "            <figcaption>" +
      escapeHtml((his.label || "男戒") + " · " + (his.name || "")) +
      "</figcaption>";
    html += "          </figure>";
    html += "        </div>";
    html +=
      '        <p class="ring-result-brand">' +
      escapeHtml(
        [real.brand, real.metal].filter(Boolean).join(" · ") || ""
      ) +
      "</p>";
    html +=
      '        <p class="ring-result-desc">' +
      escapeHtml(real.desc || "") +
      "</p>";
    html += "      </div>";
    html += "    </div>";
    html +=
      '    <p class="ring-result-vow">' +
      escapeHtml(data.finalVow || "") +
      "</p>";
    html +=
      '    <button type="button" class="btn-primary ring-confirm-btn" id="ring-result-done">打卡这一站 ✓</button>';
    html += "  </div></div>";

    var old = $("ring-result-content");
    if (old) old.remove();
    $("ring-design-wrap").insertAdjacentHTML("beforebegin", html);

    $("ring-result-done").addEventListener("click", function () {
      var resEl = $("ring-result-content");
      if (resEl) resEl.remove();
      if (activeStation) markCompleted(activeStation.id);
      closeOverlay("ringdesign");
      advanceToNextStation();
    });
  }

  function startRingDesign(station) {
    var data = getRingData();
    ringParamIndex = 0;
    ringPhase = "param";
    ringChoices = {};

    $("ring-intro").textContent = data.intro || "";

    var old = $("ring-result-content");
    if (old) old.remove();

    showRingPhase("param");
  }

  /* ── 云南站：看图猜地点 ── */
  function getComicMatchData() {
    return content.yunnanMatch || { title: "", intro: "", passCount: 6, rounds: [] };
  }

  function shuffleArray(arr) {
    var a = arr.slice();
    for (var i = a.length - 1; i > 0; i--) {
      var j = Math.floor(Math.random() * (i + 1));
      var t = a[i];
      a[i] = a[j];
      a[j] = t;
    }
    return a;
  }

  function updateComicProgress() {
    var data = getComicMatchData();
    var need = data.passCount || 6;
    if (el.comicProgress) {
      el.comicProgress.textContent = comicCorrect + " / " + need;
    }
  }

  function setComicGuideLine(text, mood) {
    if (el.comicIntro) el.comicIntro.textContent = text || "";
    if (!el.comicGuide) return;
    el.comicGuide.classList.remove("is-talk", "is-wrong");
    if (mood === "correct") el.comicGuide.classList.add("is-talk");
    if (mood === "wrong") el.comicGuide.classList.add("is-wrong");
  }

  function getComicGuideLines() {
    var data = getComicMatchData();
    return data.guideLines || {};
  }

  function loadComicPanel(src, altText) {
    if (!el.comicImage) return;
    el.comicImage.classList.add("is-loading");
    el.comicImage.onload = function () {
      el.comicImage.classList.remove("is-loading");
    };
    el.comicImage.onerror = function () {
      el.comicImage.classList.remove("is-loading");
      el.comicImage.alt = "图片加载失败";
    };
    el.comicImage.alt = altText || "旅途照片";
    el.comicImage.decoding = "async";
    el.comicImage.loading = "eager";
    if (el.comicImage.getAttribute("src") === src && el.comicImage.complete) {
      el.comicImage.classList.remove("is-loading");
      return;
    }
    el.comicImage.src = src;
  }

  function restartComicMatch(message) {
    var data = getComicMatchData();
    var lines = getComicGuideLines();
    comicBusy = true;
    comicCorrect = 0;
    comicIndex = 0;
    comicQueue = shuffleArray(data.rounds || []);
    updateComicProgress();
    setComicGuideLine(message || lines.wrong || "答错了，从头再来～", "wrong");
    setTimeout(function () {
      renderComicRound();
    }, 1100);
  }

  function startComicMatch(station) {
    var data = getComicMatchData();
    comicBusy = false;
    comicCorrect = 0;
    comicIndex = 0;
    comicShowIntroOnce = true;
    comicQueue = shuffleArray(data.rounds || []);

    if (el.comicGuideName) {
      el.comicGuideName.textContent = data.guideName || "铜脸向导";
    }
    if (el.comicGuideAvatar && data.guideAvatar) {
      el.comicGuideAvatar.src = data.guideAvatar;
    }
    updateComicProgress();

    (data.rounds || []).forEach(function (round) {
      if (!round || !round.image) return;
      var pre = new Image();
      pre.decoding = "async";
      pre.src = round.image;
    });

    renderComicRound();
  }

  function renderComicRound() {
    var data = getComicMatchData();
    var lines = getComicGuideLines();
    var round = comicQueue[comicIndex];
    if (!round) {
      comicQueue = shuffleArray(data.rounds || []);
      comicIndex = 0;
      round = comicQueue[0];
    }
    if (!round) return;

    comicBusy = false;
    if (el.comicFeedback) el.comicFeedback.textContent = "";

    if (comicShowIntroOnce) {
      comicShowIntroOnce = false;
      setComicGuideLine(data.intro || lines.start || "", "idle");
    } else {
      setComicGuideLine(
        lines.start || "瞪大眼睛看仔细——这张，汝可知何处？",
        "idle"
      );
    }
    loadComicPanel(round.image, round.comicTitle || "旅途照片");

    var options = shuffleArray(round.options || []);
    el.comicOptions.innerHTML = "";
    options.forEach(function (opt, i) {
      var btn = document.createElement("button");
      btn.type = "button";
      btn.className = "option-btn photo-option";
      btn.innerHTML =
        '<span class="photo-opt-mark">' +
        String.fromCharCode(65 + i) +
        '</span><span class="photo-opt-text">' +
        opt +
        "</span>";
      btn.addEventListener("click", function () {
        onComicOption(opt, btn, round);
      });
      el.comicOptions.appendChild(btn);
    });
  }

  function onComicOption(choice, btn, round) {
    if (comicBusy) return;
    comicBusy = true;

    var buttons = el.comicOptions.querySelectorAll(".option-btn");
    buttons.forEach(function (b) {
      b.classList.remove("correct", "wrong");
      b.disabled = true;
    });

    var data = getComicMatchData();
    var lines = getComicGuideLines();
    var need = data.passCount || 6;

    if (choice === round.correct) {
      btn.classList.add("correct");
      comicCorrect += 1;
      updateComicProgress();
      setComicGuideLine(lines.correct || "答对啦！", "correct");

      setTimeout(function () {
        if (comicCorrect >= need) {
          setComicGuideLine(lines.pass || "过关！", "correct");
          var station = activeStation;
          setTimeout(function () {
            if (station && station.followUp === "puzzle") {
              openOverlay("puzzle", station);
              return;
            }
            if (station) markCompleted(station.id);
            closeOverlay("comicmatch");
            advanceToNextStation();
          }, 700);
          return;
        }
        comicIndex += 1;
        renderComicRound();
      }, 550);
    } else {
      btn.classList.add("wrong");
      var hint = round.wrongHint
        ? round.wrongHint + " " + (lines.wrong || "进度清零，重新开始～")
        : lines.wrong || "答错了，进度清零，重新开始～";
      restartComicMatch(hint);
    }
  }

  function renderQuestion() {
    var q = content.questions[quizIndex];
    el.questionText.textContent = q.text;
    el.hintMsg.textContent = "";
    el.quizProgress.textContent =
      "第 " + (quizIndex + 1) + " / " + content.questions.length + " 题";

    el.options.innerHTML = "";
    q.options.forEach(function (opt, i) {
      var btn = document.createElement("button");
      btn.type = "button";
      btn.className = "option-btn";
      btn.textContent = opt;
      btn.addEventListener("click", function () {
        onOptionClick(i, btn);
      });
      el.options.appendChild(btn);
    });
  }

  function onOptionClick(index, btn) {
    var q = content.questions[quizIndex];
    var buttons = el.options.querySelectorAll(".option-btn");
    buttons.forEach(function (b) {
      b.classList.remove("correct", "wrong");
    });

    if (index === q.correctIndex) {
      btn.classList.add("correct");
      el.hintMsg.textContent = "答对啦！❤️";
      setTimeout(function () {
        quizIndex += 1;
        saveProgress();
        if (quizIndex >= content.questions.length) {
          var station = activeStation || content.stations.find(function (s) {
            return s.type === "quiz";
          });
          if (station) markCompleted(station.id);
          quizIndex = 0;
          saveProgress();
          closeOverlay("quiz");
          advanceToNextStation();
        } else {
          renderQuestion();
        }
      }, 600);
    } else {
      btn.classList.add("wrong");
      el.hintMsg.textContent = q.wrongHint || "再试一次吧～";
    }
  }

  function startLetter() {
    var letter = content.letter;
    el.letterBody.innerHTML = "";
    $("letter-done-btn").classList.add("hidden");

    var greeting = document.createElement("p");
    greeting.className = "letter-greeting";
    greeting.textContent = letter.greeting;
    el.letterBody.appendChild(greeting);

    var fullText = letter.paragraphs.join("\n\n");
    var paraEl = document.createElement("p");
    el.letterBody.appendChild(paraEl);

    var sigWrap = document.createElement("div");
    sigWrap.className = "letter-signature hidden";
    sigWrap.innerHTML =
      "<p>" +
      letter.signature +
      "</p><p>" +
      (letter.signatureDate || "") +
      "</p>";
    el.letterBody.appendChild(sigWrap);

    var i = 0;
    function type() {
      if (i <= fullText.length) {
        paraEl.textContent = fullText.slice(0, i);
        i += 1;
        setTimeout(type, 40);
      } else {
        sigWrap.classList.remove("hidden");
        $("letter-done-btn").classList.remove("hidden");
      }
    }
    type();
  }

  function spawnConfetti() {
    var colors = ["#4caf50", "#81c784", "#a5d6a7", "#ffd54f", "#c8e6c9"];
    for (var n = 0; n < 40; n++) {
      (function (j) {
        setTimeout(function () {
          var piece = document.createElement("div");
          piece.className = "confetti-piece";
          piece.style.left = Math.random() * 100 + "vw";
          piece.style.top = "-10px";
          piece.style.background = colors[j % colors.length];
          piece.style.animationDuration = 2 + Math.random() * 2 + "s";
          document.body.appendChild(piece);
          setTimeout(function () {
            piece.remove();
          }, 4000);
        }, j * 30);
      })(n);
    }
  }

  function initFloats() {
    var container = $("hearts-bg");
    var floats = ["🚄", "✨", "🍃", "🌿", "💚"];
    for (var i = 0; i < 10; i++) {
      var span = document.createElement("span");
      span.className = "heart-float";
      span.textContent = floats[i % floats.length];
      span.style.left = Math.random() * 100 + "%";
      span.style.animationDuration = 8 + Math.random() * 8 + "s";
      span.style.animationDelay = Math.random() * 5 + "s";
      container.appendChild(span);
    }
  }

  function initStars() {
    var container = $("stars-bg");
    for (var i = 0; i < 12; i++) {
      var star = document.createElement("span");
      star.className = "star-twinkle";
      star.textContent = "✦";
      star.style.left = Math.random() * 100 + "%";
      star.style.top = Math.random() * 40 + "%";
      container.appendChild(star);
    }
  }

  var resizeTimer = null;

  function bindSwipe() {
    var wrap = el.railMapWrap;
    wrap.addEventListener(
      "touchstart",
      function (e) {
        touchStartX = e.changedTouches[0].screenX;
      },
      { passive: true }
    );
    wrap.addEventListener(
      "touchend",
      function (e) {
        var diff = e.changedTouches[0].screenX - touchStartX;
        if (Math.abs(diff) > 40) {
          var next = stationIndex + (diff < 0 ? 1 : -1);
          if (next > unlockedUpTo) {
            flashLockedHint();
            return;
          }
          if (next < 0) return;
          goToStation(next, true);
        }
      },
      { passive: true }
    );
  }

  function bindEvents() {
    $("btn-start").addEventListener("click", showJourney);
    $("btn-enter-station").addEventListener("click", enterCurrentStation);
    $("btn-prev").addEventListener("click", function () {
      if (stationIndex <= 0) return;
      goToStation(stationIndex - 1, true);
    });
    $("btn-next").addEventListener("click", function () {
      if (stationIndex >= unlockedUpTo) {
        flashLockedHint();
        return;
      }
      goToStation(stationIndex + 1, true);
    });
    $("btn-judge-true").addEventListener("click", function () {
      onJudgeAnswer(true);
    });
    $("btn-judge-false").addEventListener("click", function () {
      onJudgeAnswer(false);
    });

    $("btn-puzzle-shuffle").addEventListener("click", function () {
      if (el.puzzleMsg) el.puzzleMsg.textContent = "";
      var doneBtn = $("btn-puzzle-done");
      if (doneBtn) doneBtn.classList.add("hidden");
      shufflePuzzle();
      renderPuzzle();
    });

    $("btn-puzzle-peek").addEventListener("click", function () {
      if (!el.puzzlePeek) return;
      el.puzzlePeek.classList.toggle("hidden");
    });

    $("btn-puzzle-done").addEventListener("click", function () {
      if (activeStation) markCompleted(activeStation.id);
      closeOverlay("puzzle");
      advanceToNextStation();
    });

    $("btn-proposal-start").addEventListener("click", gotoMatchPhase);
    $("btn-proposal-ask").addEventListener("click", gotoAskPhase);
    $("btn-proposal-yes").addEventListener("click", gotoAcceptPhase);
    $("btn-proposal-no").addEventListener("click", dodgeProposalNo);
    $("btn-proposal-done").addEventListener("click", function () {
      if (activeStation) markCompleted(activeStation.id);
      closeOverlay("proposal");
      advanceToNextStation();
    });

    $("btn-sil-start").addEventListener("click", function () {
      var data = getSilData();
      silIndex = 0;
      ["sil-intro", "sil-scene-phase", "sil-outro"].forEach(function (id) {
        $(id).classList.toggle("active", id === "sil-scene-phase");
      });
      renderSilScene(0);
    });
    $("btn-sil-next").addEventListener("click", silNext);
    $("btn-sil-prev").addEventListener("click", silPrev);
    $("btn-sil-done").addEventListener("click", function () {
      if (activeStation) markCompleted(activeStation.id);
      closeOverlay("silhouette");
      advanceToNextStation();
    });

    $("ring-prev-btn").addEventListener("click", ringPrev);
    $("ring-next-btn").addEventListener("click", ringNext);
    $("ring-finalize-btn").addEventListener("click", finalizeRing);

    $("letter-done-btn").addEventListener("click", function () {
      var station = content.stations.find(function (s) {
        return s.type === "letter";
      });
      if (station) markCompleted(station.id);
      closeOverlay("letter");
      advanceToNextStation();
    });

    $("memory-done-btn").addEventListener("click", function () {
      if (activeStation) markCompleted(activeStation.id);
      closeOverlay("memory");
      advanceToNextStation();
    });

    el.cake.addEventListener("click", function () {
      el.cake.classList.add("tapped");
      spawnConfetti();
      var station = content.stations.find(function (s) {
        return s.type === "finale";
      });
      if (station) markCompleted(station.id);
      var count = content.stations.length;
      if (stationIndex === count - 1) {
        unlockedUpTo = count - 1;
        saveProgress();
        updateStationStates();
      }
      setTimeout(function () {
        el.cake.classList.remove("tapped");
      }, 600);
    });

    document.querySelectorAll(".btn-back").forEach(function (btn) {
      btn.addEventListener("click", function () {
        closeOverlay(btn.dataset.close);
      });
    });

    window.addEventListener("resize", function () {
      if (!el.journey.classList.contains("active")) return;
      if (resizeTimer) clearTimeout(resizeTimer);
      resizeTimer = setTimeout(function () {
        positionTrain(false);
      }, 150);
    });
  }

  function restoreOrStart() {
    var saved = loadProgress();
    if (saved) {
      stationIndex = saved.stationIndex || 0;
      quizIndex = saved.quizIndex || 0;
      judgeStreak = saved.judgeStreak || 0;
      completedStations = saved.completed || {};
      unlockedUpTo = typeof saved.unlockedUpTo === "number" ? saved.unlockedUpTo : 0;
      if (saved.seenCover) {
        showJourney();
      } else {
        showCover();
      }
    } else {
      showCover();
    }
    updateStationStates();
  }

  function showLoadError() {
    document.querySelector(".app").innerHTML =
      '<div class="card"><h2>加载失败</h2><p class="subtitle">请检查网络，或通过本地服务器 / GitHub Pages 访问。</p></div>';
  }

  initElements();
  initFloats();
  initStars();
  bindEvents();
  bindSwipe();

  loadContent()
    .then(restoreOrStart)
    .catch(showLoadError);
})();
