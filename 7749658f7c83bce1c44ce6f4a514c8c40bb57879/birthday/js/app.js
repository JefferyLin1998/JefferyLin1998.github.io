(function () {
  "use strict";

  var SCREENS = ["cover", "secret", "quiz", "letter", "finale"];
  var STORAGE_KEY = "birthday-gift-progress";

  var content = null;
  var currentScreen = 0;
  var quizIndex = 0;
  var selectedOption = -1;

  var el = {
    screens: {},
    progressFill: null,
    stepLabel: null,
    herNameCover: null,
    secretInput: null,
    secretHint: null,
    secretError: null,
    questionText: null,
    options: null,
    hintMsg: null,
    nextQuizBtn: null,
    letterBody: null,
    daysCount: null,
    finaleTitle: null,
    finaleSubtitle: null,
    cake: null,
  };

  function $(id) {
    return document.getElementById(id);
  }

  function initElements() {
    SCREENS.forEach(function (name) {
      el.screens[name] = $("screen-" + name);
    });
    el.progressFill = $("progress-fill");
    el.stepLabel = $("step-label");
    el.herNameCover = $("her-name-cover");
    el.secretInput = $("secret-input");
    el.secretHint = $("secret-hint");
    el.secretError = $("secret-error");
    el.questionText = $("question-text");
    el.options = $("options");
    el.hintMsg = $("hint-msg");
    el.nextQuizBtn = $("next-quiz-btn");
    el.letterBody = $("letter-body");
    el.daysCount = $("days-count");
    el.finaleTitle = $("finale-title");
    el.finaleSubtitle = $("finale-subtitle");
    el.cake = $("finale-cake");
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
      });
  }

  function applyContent() {
    el.herNameCover.textContent = content.herName;
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

  function saveProgress(screenIndex, quiz) {
    try {
      localStorage.setItem(
        STORAGE_KEY,
        JSON.stringify({ screen: screenIndex, quiz: quiz || 0 })
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

  function updateProgress() {
    var total = 3 + content.questions.length;
    var done = 0;
    if (currentScreen > 0) done += 1;
    if (currentScreen > 1) done += 1;
    if (currentScreen >= 2) done += quizIndex + (currentScreen > 2 ? content.questions.length - quizIndex : 0);
    if (currentScreen > 2) done = 3 + content.questions.length;

    var pct = Math.min(100, Math.round((done / total) * 100));
    if (el.progressFill) el.progressFill.style.width = pct + "%";

    var labels = ["", "暗号", "回忆问答", "情书", "祝福"];
    if (currentScreen >= 2 && currentScreen < 3) {
      el.stepLabel.textContent =
        "第 " + (quizIndex + 1) + " / " + content.questions.length + " 题";
    } else if (labels[currentScreen]) {
      el.stepLabel.textContent = labels[currentScreen];
    }
  }

  function showScreen(index) {
    currentScreen = index;
    SCREENS.forEach(function (name, i) {
      el.screens[name].classList.toggle("active", i === index);
    });
    $("progress-wrap").classList.toggle("hidden", index === 0 || index === 4);
    updateProgress();
    saveProgress(index, quizIndex);
  }

  function normalizeAnswer(s) {
    return String(s || "")
      .trim()
      .toLowerCase();
  }

  function checkSecret() {
    var val = normalizeAnswer(el.secretInput.value);
    var ans = normalizeAnswer(content.secretAnswer);
    if (val === ans) {
      el.secretError.textContent = "";
      quizIndex = 0;
      renderQuestion();
      showScreen(2);
    } else {
      el.secretError.textContent = "不对哦，再想想～";
      el.secretInput.value = "";
    }
  }

  function renderQuestion() {
    var q = content.questions[quizIndex];
    el.questionText.textContent = q.text;
    el.hintMsg.textContent = "";
    selectedOption = -1;
    el.nextQuizBtn.disabled = true;
    el.nextQuizBtn.classList.add("hidden");

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
    updateProgress();
  }

  function onOptionClick(index, btn) {
    var q = content.questions[quizIndex];
    var buttons = el.options.querySelectorAll(".option-btn");
    buttons.forEach(function (b) {
      b.classList.remove("selected", "correct", "wrong");
    });

    if (index === q.correctIndex) {
      btn.classList.add("correct");
      el.hintMsg.textContent = "答对啦！❤️";
      selectedOption = index;
      setTimeout(function () {
        quizIndex += 1;
        if (quizIndex >= content.questions.length) {
          startLetter();
          showScreen(3);
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
    var speed = 40;

    function type() {
      if (i <= fullText.length) {
        paraEl.textContent = fullText.slice(0, i);
        i += 1;
        setTimeout(type, speed);
      } else {
        sigWrap.classList.remove("hidden");
        $("letter-next-btn").classList.remove("hidden");
      }
    }
    type();
  }

  function spawnConfetti() {
    var colors = ["#e91e63", "#f48fb1", "#fce4ec", "#ffeb3b", "#ce93d8"];
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

  function initHearts() {
    var container = $("hearts-bg");
    var hearts = ["💕", "💖", "💗", "🌸", "✨"];
    for (var i = 0; i < 12; i++) {
      var span = document.createElement("span");
      span.className = "heart-float";
      span.textContent = hearts[i % hearts.length];
      span.style.left = Math.random() * 100 + "%";
      span.style.animationDuration = 8 + Math.random() * 8 + "s";
      span.style.animationDelay = Math.random() * 5 + "s";
      container.appendChild(span);
    }
  }

  function bindEvents() {
    $("btn-start").addEventListener("click", function () {
      showScreen(1);
    });

    $("btn-secret").addEventListener("click", checkSecret);

    el.secretInput.addEventListener("keydown", function (e) {
      if (e.key === "Enter") checkSecret();
    });

    $("letter-next-btn").addEventListener("click", function () {
      showScreen(4);
      spawnConfetti();
    });

    el.cake.addEventListener("click", function () {
      el.cake.classList.add("tapped");
      spawnConfetti();
      setTimeout(function () {
        el.cake.classList.remove("tapped");
      }, 600);
    });
  }

  function restoreOrStart() {
    var saved = loadProgress();
    if (saved && saved.screen > 0) {
      quizIndex = saved.quiz || 0;
      if (saved.screen === 2) {
        renderQuestion();
      }
      if (saved.screen === 3) {
        startLetter();
      }
      showScreen(Math.min(saved.screen, SCREENS.length - 1));
    } else {
      showScreen(0);
    }
  }

  function showLoadError() {
    document.querySelector(".app").innerHTML =
      '<div class="card"><h2>加载失败</h2><p class="subtitle">请检查网络，或通过 GitHub Pages 访问本页面。</p></div>';
  }

  initElements();
  initHearts();
  bindEvents();

  loadContent()
    .then(restoreOrStart)
    .catch(showLoadError);
})();
