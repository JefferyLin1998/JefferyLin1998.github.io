(function () {
  "use strict";

  var STORAGE_KEY = "birthday-gift-progress-v7";

  var content = null;
  var stationIndex = 0;
  var quizIndex = 0;
  var judgeStreak = 0;
  var judgeUsedInSession = [];
  var judgeCurrentStmt = null;
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
