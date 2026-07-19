(function () {
  "use strict";

  var STORAGE_KEY = "birthday-gift-progress-v7";

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
    var minWidth = Math.max(320, stations.length * 58 + 48);
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
    var trainW = el.train.offsetWidth || 72;

    var left = dotRect.left - mapRect.left + dotRect.width / 2 - trainW / 2;
    var trackEl = el.railMap.querySelector(".hsr-route-track");
    var trackRect = trackEl ? trackEl.getBoundingClientRect() : mapRect;
    var top = trackRect.top - mapRect.top + trackRect.height / 2 - 18;

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
    var overlay = $("overlay-" + type);
    if (!overlay) return;
    overlay.hidden = false;
    document.body.classList.add("overlay-open");

    if (type === "judge") {
      startJudgeGame();
    }

    if (type === "puzzle" && station) {
      $("puzzle-title").textContent = station.icon + " " + station.city + "站 · 拼图";
      startPuzzle();
    }

    if (type === "proposal" && station) {
      startProposal(station);
    }

    if (type === "silhouette" && station) {
      startSilhouette(station);
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
