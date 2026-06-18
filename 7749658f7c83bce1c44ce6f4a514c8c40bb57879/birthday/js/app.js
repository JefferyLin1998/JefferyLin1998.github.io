(function () {
  "use strict";

  var STORAGE_KEY = "birthday-gift-progress-v4";

  var content = null;
  var stationIndex = 0;
  var quizIndex = 0;
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
    el.secretHint = $("secret-hint");
    el.secretInput = $("secret-input");
    el.secretError = $("secret-error");
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
    el.secretHint.textContent = content.secretHint;
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
  }

  function buildRailMap() {
    var stations = content.stations;
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

    if (animate) {
      el.train.classList.add("is-running");
      if (trainRunTimer) clearTimeout(trainRunTimer);
      trainRunTimer = setTimeout(function () {
        el.train.classList.remove("is-running");
      }, 750);
    }

    if (!animate) {
      requestAnimationFrame(function () {
        el.train.classList.remove("no-transition");
      });
    }

    var fillPct = count <= 1 ? 100 : (stationIndex / (count - 1)) * 100;
    el.lineFill.style.width = fillPct + "%";
  }

  function updateStationFocus() {
    var station = content.stations[stationIndex];
    var ordinals = ["第一站", "第二站", "第三站", "第四站", "第五站", "第六站"];

    el.focusStationLabel.textContent = ordinals[stationIndex] || "第 " + (stationIndex + 1) + " 站";
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

    if (type === "secret" && station) {
      $("overlay-secret-city").textContent = station.city;
      $("secret-station-title").textContent = station.city + "站";
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
      startLetter();
    }

    if (type === "finale") {
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
      positionTrain(true);
    });
  }

  function closeAllOverlays() {
    document.querySelectorAll(".project-overlay").forEach(function (o) {
      o.hidden = true;
    });
    document.body.classList.remove("overlay-open");
    activeStation = null;
  }

  function normalizeAnswer(s) {
    return String(s || "")
      .trim()
      .toLowerCase();
  }

  function checkSecret() {
    var val = normalizeAnswer(el.secretInput.value);
    var ans = normalizeAnswer(content.secretAnswer);
    var station = activeStation || content.stations[0];
    if (val === ans) {
      el.secretError.textContent = "答对啦！❤️ 高铁即将发车…";
      markCompleted(station.id);
      setTimeout(function () {
        closeOverlay("secret");
        el.secretInput.value = "";
        el.secretError.textContent = "";
        advanceToNextStation();
      }, 900);
    } else {
      el.secretError.textContent = "不对哦，再想想～";
      el.secretInput.value = "";
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
    for (var i = 0; i < 20; i++) {
      var star = document.createElement("span");
      star.className = "star-twinkle";
      star.textContent = "✦";
      star.style.left = Math.random() * 100 + "%";
      star.style.top = Math.random() * 40 + "%";
      star.style.animationDelay = Math.random() * 3 + "s";
      container.appendChild(star);
    }
  }

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
    $("btn-secret").addEventListener("click", checkSecret);
    el.secretInput.addEventListener("keydown", function (e) {
      if (e.key === "Enter") checkSecret();
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
      if (el.journey.classList.contains("active")) {
        positionTrain(false);
      }
    });
  }

  function restoreOrStart() {
    var saved = loadProgress();
    if (saved) {
      stationIndex = saved.stationIndex || 0;
      quizIndex = saved.quizIndex || 0;
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
