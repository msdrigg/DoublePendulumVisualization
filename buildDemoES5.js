"use strict";

function _instanceof(left, right) { if (right != null && typeof Symbol !== "undefined" && right[Symbol.hasInstance]) { return !!right[Symbol.hasInstance](left); } else { return left instanceof right; } }

function _classCallCheck(instance, Constructor) { if (!_instanceof(instance, Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var DoublePendulumDemo = /*#__PURE__*/function () {
  function DoublePendulumDemo(constantsDict, animationLock) {
    _classCallCheck(this, DoublePendulumDemo);

    _defineProperty(this, "stepSize", 0.001);

    _defineProperty(this, "repeats", 32);

    _defineProperty(this, "randomness", .4);

    _defineProperty(this, "started", false);

    _defineProperty(this, "currentCoords", void 0);

    _defineProperty(this, "currentTime", void 0);

    _defineProperty(this, "currentCartCoords", void 0);

    _defineProperty(this, "width", void 0);

    _defineProperty(this, "height", void 0);

    _defineProperty(this, "origin", void 0);

    _defineProperty(this, "dots1", void 0);

    _defineProperty(this, "lines1", void 0);

    _defineProperty(this, "dots2", void 0);

    _defineProperty(this, "lines2", void 0);

    _defineProperty(this, "trails", void 0);

    _defineProperty(this, "continueLooping", false);

    _defineProperty(this, "restartable", false);

    _defineProperty(this, "colorscaleVariable", 0);

    _defineProperty(this, "currentStep", 0);

    _defineProperty(this, "lineResizer", 0.3);

    _defineProperty(this, "dotResizer", 0.1);

    _defineProperty(this, "lineGenerator", d3.line().curve(d3.curveNatural));

    this.isReset = true;
    this.lock = animationLock;
    this.constants = [constantsDict.l1, constantsDict.m1, constantsDict.l2, constantsDict.m2];
    this.phi1Init = constantsDict.phi1Init;
    this.phi2Init = constantsDict.phi2Init;
    this.pendulumNumber = constantsDict.pendulumNumber;
    this.deviation = constantsDict.deviation;
    this.location = constantsDict.location;
    this.trailUpdateInterval = constantsDict.trailUpdateInterval;
    this.trailLength = constantsDict.trailLength;
    this.useTrails = constantsDict.trails;
    this.trailPaths = new Array(this.pendulumNumber);
    this.trailPathsNumeric = new Array(this.pendulumNumber);
    this.trailCounter = 0;
    this.windowSize = 100 / (2 * 1.05 * (this.l1_ + this.l2_));
    this.dotsResizer = 0.1 / this.windowSize;
    this.caption = constantsDict.caption;

    if ('colorScale' in constantsDict) {
      this.colorScale = constants.colorScale;
    } else {
      this.colorScale = d3.interpolateInferno;
    }

    this.currentTime = 0;
    this.initialCoords = this.getInitialCoords();
    this.lastCoords = this.initialCoords.clone();
    this.initialCoordsCart = convertToCoordinates(this.initialCoords, this.constants);
    this.lastTrailUpdateCoords = this.initialCoordsCart.slice(2);
    this.svg = d3.select(this.location).append("svg").attr("viewBox", "0 0 100 100");
    var explanation = "Hover to play. Click to restart";

    if (constantsDict.usingMobile) {
      explanation = "Tap to play. Double tap to restart";
    }

    var caption = d3.select(this.location).append("p").attr("class", "demoCaption").html(this.caption);

    if (constantsDict.explain) {
      caption.append("div").attr("class", "demoSubcaption").html(explanation);
    }
  }

  _createClass(DoublePendulumDemo, [{
    key: "start",
    value: function start() {
      if (this.lock.request()) {
        return;
      }

      if (!this.continueLooping) {
        this.isReset = false;
        this.continueLooping = true;
        var demoContext = this;
        requestAnimationFrame(function () {
          return demoContext.updateDisplay();
        });
      }
    }
  }, {
    key: "stop",
    value: function stop() {
      if (this.continueLooping) this.continueLooping = false;
    }
  }, {
    key: "restart",
    value: function restart() {
      this.stop();

      for (var i = 0; i < this.pendulumNumber; i++) {
        this.trailPaths[i] = "M " + (this.initialCoordsCart[i][2] * this.windowSize).toPrecision(5) + " " + (this.initialCoordsCart[i][3] * this.windowSize).toPrecision(5) + " ";
      }

      this.currentCoords = this.initialCoords;
      this.lastCoords = this.currentCoords.clone();
      this.currentCartCoords = this.initialCoordsCart;
      this.drawPendulumsAt(this.initialCoordsCart);
      this.isReset = true;
    }
  }, {
    key: "getInitialCoords",
    value: function getInitialCoords() {
      return initPendulums(this.phi1Init + (Math.random() - 0.5) * this.randomness, this.phi2Init + (Math.random() - 0.5) * this.randomness, this.deviation, this.pendulumNumber);
    }
  }, {
    key: "init",
    value: function init() {
      var _this = this;

      for (var i = 0; i < this.pendulumNumber; i++) {
        this.trailPaths[i] = "M " + (this.initialCoordsCart[i][2] * this.windowSize).toPrecision(5) + " " + (this.initialCoordsCart[i][3] * this.windowSize).toPrecision(5) + " ";
      }

      this.origin = this.svg.append('g').attr("transform", "translate(" + 50 + "," + 50 + ")").attr("class", "origin");
      this.trails = this.origin.selectAll(".trails").data(this.initialCoordsCart).enter().append('path').attr("fill", "transparent").attr("class", "trails").attr("stroke-width", this.lineResizer).attr("d", function (d) {
        return "M " + (d[2] * _this.windowSize).toPrecision(5) + " " + (d[3] * _this.windowSize).toPrecision(5) + " ";
      }).attr("stroke", this.colorScaleFunc.bind(this));
      var pendulums = this.origin.selectAll('.pendulums').data(this.initialCoordsCart).enter().append('g').attr("class", "pendulums");
      this.lines1 = pendulums.append('line').attr("class", "pendulumLine1").attr("x1", 0).attr("y1", 0).attr("stroke-width", 2 * this.lineResizer).attr("x2", function (d) {
        return d[0] * _this.windowSize;
      }).attr("y2", function (d) {
        return d[1] * _this.windowSize;
      }).attr("stroke", this.colorScaleFunc.bind(this));
      this.lines2 = pendulums.append('line').attr("class", "pendulumLine2").attr("stroke-width", 2 * this.lineResizer).attr("x1", function (d) {
        return d[0] * _this.windowSize;
      }).attr("y1", function (d) {
        return d[1] * _this.windowSize;
      }).attr("x2", function (d) {
        return d[2] * _this.windowSize;
      }).attr("y2", function (d) {
        return d[3] * _this.windowSize;
      }).attr("stroke", this.colorScaleFunc.bind(this));
      this.dots1 = pendulums.append('circle').attr("class", "dots1").attr("fill", this.colorScaleFunc.bind(this)).attr("cx", function (d) {
        return d[0] * _this.windowSize;
      }).attr("cy", function (d) {
        return d[1] * _this.windowSize;
      }).attr("r", this.constants[1] * this.dotResizer * this.windowSize);
      this.dots2 = pendulums.append('circle').attr("class", "dots2").attr("fill", this.colorScaleFunc.bind(this)).attr("cx", function (d) {
        return d[2] * _this.windowSize;
      }).attr("cy", function (d) {
        return d[3] * _this.windowSize;
      }).attr("r", this.constants[3] * this.dotResizer * this.windowSize);
      this.origin.append("circle").attr("class", "origin-mark").attr("r", Math.max(this.constants[2], this.constants[1]) * 1.5 * this.dotResizer * this.windowSize).style("fill", this.colorScaleFunc.bind(this)(null, (this.pendulumNumber - 1) / (this.pendulumNumber + 5)));
      this.currentCoords = this.initialCoords.clone();
    }
  }, {
    key: "colorScaleFunc",
    value: function colorScaleFunc(d, i) {
      if (this.pendulumNumber == 1) {
        if (d != null) return this.colorScale(.95);
        return d3.rgb(this.colorScale(.95)).hex();
      }

      if (d == null) return d3.rgb(this.colorScale(.95)).hex();
      return this.colorScale(i / (this.pendulumNumber + 5) + 2 / this.pendulumNumber);
    }
  }, {
    key: "updateTrailDataWith",
    value: function updateTrailDataWith(newCoords) {
      this.trailCounter++;
      if (this.trailLength == 0) return;

      if (this.trailCounter > this.trailLength) {
        for (var i = 0; i < this.pendulumNumber; i++) {
          this.trailPaths[i] = "M" + this.trailPaths[i].slice(this.trailPaths[i].indexOf("L") + 1) + "L " + (newCoords[i][2] * this.windowSize).toPrecision(5) + " " + (newCoords[i][3] * this.windowSize).toPrecision(5) + " ";
        }
      } else {
        for (var _i = 0; _i < this.pendulumNumber; _i++) {
          this.trailPaths[_i] += "L " + (newCoords[_i][2] * this.windowSize).toPrecision(5) + " " + (newCoords[_i][3] * this.windowSize).toPrecision(5) + " ";
        }
      }

      this.lastCartCoords = this.newCoords;
    }
  }, {
    key: "drawPendulumsAt",
    value: function drawPendulumsAt(newCoords) {
      this.svg.selectAll("*").remove();
      this.initialCoords = this.getInitialCoords();
      this.initialCoordsCart = convertToCoordinates(this.initialCoords, this.constants);
      this.trailCounter = 0;
      this.init();
    }
  }, {
    key: "updateDisplay",
    value: function updateDisplay() {
      var _this2 = this;

      this.lastCoords = this.currentCoords.clone();

      for (var i = 0; i < this.repeats; i++) {
        this.currentCoords = RK4LA(derivativeLA, this.stepSize, this.currentTime, this.currentCoords, this.constants);
      }

      this.currentTime += this.stepSize * this.repeats;
      this.currentCartCoords = convertToCoordinates(this.currentCoords, this.constants);
      this.dots1.data(this.currentCartCoords).attr("cx", function (d) {
        return d[0] * _this2.windowSize;
      }).attr("cy", function (d) {
        return d[1] * _this2.windowSize;
      });
      this.dots2.data(this.currentCartCoords).attr("cx", function (d) {
        return d[2] * _this2.windowSize;
      }).attr("cy", function (d) {
        return d[3] * _this2.windowSize;
      });
      this.lines1.data(this.currentCartCoords).attr("x2", function (d) {
        return d[0] * _this2.windowSize;
      }).attr("y2", function (d) {
        return d[1] * _this2.windowSize;
      });
      this.lines2.data(this.currentCartCoords).attr("x1", function (d) {
        return d[0] * _this2.windowSize;
      }).attr("y1", function (d) {
        return d[1] * _this2.windowSize;
      }).attr("x2", function (d) {
        return d[2] * _this2.windowSize;
      }).attr("y2", function (d) {
        return d[3] * _this2.windowSize;
      });

      if (this.useTrails && !(this.currentStep % this.trailUpdateInterval)) {
        this.lastTrailUpdateCoords = this.currentCartCoords.slice(2);
        this.updateTrailDataWith(this.currentCartCoords);
        this.trails.data(this.trailPaths).attr("d", function (d) {
          return d;
        });
      }

      this.currentStep += 1;

      if (this.continueLooping) {
        requestAnimationFrame(function () {
          return _this2.updateDisplay();
        });
      } else {
        this.lock.release();
      }
    }
  }, {
    key: "l1_",
    get: function get() {
      return this.constants[0];
    },
    set: function set(newValue) {
      this.constants[0] = newValue;
      this.windowSize = 100 / (2 * 1.05 * (this.l1_ + this.l2_));
      this.dotsResizer = 0.1 / this.windowSize;
      this.restart.bind(this)();
    }
  }, {
    key: "m1_",
    get: function get() {
      return this.constants[1];
    },
    set: function set(newValue) {
      this.constants[1] = newValue;
      this.restart.bind(this)();
    }
  }, {
    key: "l2_",
    get: function get() {
      return this.constants[2];
    },
    set: function set(newValue) {
      this.constants[2] = newValue;
      this.windowSize = 100 / (2 * 1.05 * (this.l1_ + this.l2_));
      this.dotsResizer = 0.1 / this.windowSize;
      this.restart.bind(this)();
    }
  }, {
    key: "m2_",
    get: function get() {
      return this.constants[3];
    },
    set: function set(newValue) {
      this.constants[3] = newValue;
      this.restart.bind(this)();
    }
  }, {
    key: "pendulumNumber_",
    set: function set(newValue) {
      this.pendulumNumber = newValue;
      this.initialCoords = this.getInitialCoords();
      this.currentCoords = this.initialCoords.clone();
      this.currentCartCoords = convertToCoordinates(this.currentCoords, this.constants);
      this.initialCoordsCart = convertToCoordinates(this.initialCoords, this.constants); // this.restart.bind(this)();
    },
    get: function get() {
      return this.pendulumNumber;
    }
  }, {
    key: "phi1Init_",
    set: function set(newValue) {
      this.phi1Init = newValue; // this.restart.bind(this)();
    },
    get: function get() {
      return this.phi1Init;
    }
  }, {
    key: "phi2Init_",
    set: function set(newValue) {
      this.phi2Init = newValue; // this.restart.bind(this)();
    },
    get: function get() {
      return this.phi2Init;
    }
  }, {
    key: "trailLength_",
    set: function set(newValue) {
      this.trailLength = newValue; // this.restart.bind(this)();
    },
    get: function get() {
      return this.trailLength;
    }
  }]);

  return DoublePendulumDemo;
}();