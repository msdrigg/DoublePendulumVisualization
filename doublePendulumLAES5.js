"use strict";

var g = 9.81;

Float64Array.prototype.add = function (a) {
  for (var i = 0; i < this.length; i++) {
    this[i] += a[i];
  }

  return this;
};

Float64Array.prototype.mul = function (a) {
  for (var i = 0; i < this.length; i++) {
    this[i] *= a[i];
  }

  return this;
};

Float64Array.prototype.div = function (a) {
  for (var i = 0; i < this.length; i++) {
    this[i] /= a[i];
  }

  return this;
};

Float64Array.prototype.sub = function (a) {
  for (var i = 0; i < this.length; i++) {
    this[i] -= a[i];
  }

  return this;
};

Float64Array.prototype.scalarAdd = function (a) {
  for (var i = 0; i < this.length; i++) {
    this[i] += a;
  }

  return this;
};

Float64Array.prototype.scalarMul = function (a) {
  for (var i = 0; i < this.length; i++) {
    this[i] *= a;
  }

  return this;
};

Float64Array.prototype.scalarDiv = function (a) {
  for (var i = 0; i < this.length; i++) {
    this[i] /= a;
  }

  return this;
};

Float64Array.prototype.scalarSub = function (a) {
  for (var i = 0; i < this.length; i++) {
    this[i] -= a;
  }

  return this;
};

Float64Array.prototype.square = function () {
  for (var i = 0; i < this.length; i++) {
    this[i] = this[i] * this[i];
  }

  return this;
};

Float64Array.prototype.cos = function () {
  for (var i = 0; i < this.length; i++) {
    this[i] = Math.cos(this[i]);
  }

  return this;
};

Float64Array.prototype.sin = function () {
  for (var i = 0; i < this.length; i++) {
    this[i] = Math.sin(this[i]);
  }

  return this;
};

Float64Array.prototype.concat = function (a) {
  var destination = new Float64Array(this.length + a.length);
  destination.set(this);
  destination.set(a, this.length);
  return destination;
};

function concat(a, b, c, d) {
  var destination = new Float64Array(a.length + b.length + c.length + d.length);
  destination.set(a);
  destination.set(b, a.length);
  destination.set(c, a.length + b.length);
  destination.set(d, a.length + b.length + c.length);
  return destination;
}

Float64Array.prototype.clone = function () {
  var destination = new Float64Array(this.length);
  destination.set(this);
  return destination;
};

Float64Array.prototype.print = function () {
  console.log(this);
};

Array.prototype.transpose = function () {
  var newArray = new Array(this[0].length);

  for (var i = 0; i < this[0].length; i++) {
    var rowArray = new Float64Array(this.length);

    for (var j = 0; j < this.length; j++) {
      rowArray[j] = this[j][i];
    }

    ;
    newArray[i] = rowArray;
  }

  ;
  return newArray;
};

function linspace(start, stop, number) {
  var arr = new Float64Array(number);
  var step;
  if (number > 1) step = (stop - start) / (number - 1);else step = 0;

  for (var i = 0; i < number; i++) {
    arr[i] = start + step * i;
  }

  return arr;
}

function initPendulums(phi1, phi2, deviation, number) {
  // Phi1 Array contains all the values of phi1 (no deviation)
  // Phi2 Array contains all the values of phi2 but deviated such that the 
  //  	average value is phi2, and the difference between any 2 phi2's is 
  //  	deviation
  var phi1Array = linspace(phi1, phi1, number);
  var phi2Array;

  if (number > 1) {
    phi2Array = linspace(phi2 - deviation / 2 * (number - 1), phi2 + deviation / 2 * (number - 1), number);
  } else {
    phi2Array = linspace(phi2 - deviation / 2, phi2 + deviation / 2, number);
  }

  var p1Array = linspace(0, 0, number);
  var p2Array = p1Array.clone();
  var output = concat(phi1Array, p1Array, phi2Array, p2Array);
  return output;
}

function RK4LA(f, h, t, p, constants) {
  var k1 = f(t, p.clone(), constants).scalarMul(h);
  var k2 = f(t + h / 2, p.clone().add(k1.clone().scalarMul(0.5)), constants).scalarMul(h);
  var k3 = f(t + h / 2, p.clone().add(k2.clone().scalarMul(0.5)), constants).scalarMul(h);
  var k4 = f(t + h, p.clone().add(k3), constants).scalarMul(h);
  return p.clone().add(k1.clone().add(k2).add(k3).add(k4).scalarMul(1 / 6));
}

function derivativeLA(t, p, constants) {
  // p: phi1, p1, phi2, p2
  // constants: l1, m1, l2, m2
  var l1 = constants[0],
      m1 = constants[1],
      l2 = constants[2],
      m2 = constants[3];
  var vectorLength = p.length / 4;
  var phi1 = p.slice(0, vectorLength),
      p1 = p.slice(vectorLength, vectorLength * 2),
      phi2 = p.slice(vectorLength * 2, vectorLength * 3),
      p2 = p.slice(vectorLength * 3, vectorLength * 4);
  var cosdif = phi1.clone().sub(phi2).cos();
  var sindif = phi1.clone().sub(phi2).sin();
  var divisor = sindif.clone().square().scalarMul(m2).scalarAdd(m1);
  var h1 = p1.clone().mul(p2).mul(sindif).div(divisor.clone().scalarMul(l1 * l2));
  var h2 = p1.clone().square().scalarMul(m2 * l2 * l2).add(p2.clone().square().scalarMul(l1 * l1 * (m1 + m2))).sub(cosdif.clone().mul(p1).mul(p2).scalarMul(2 * m2 * l2 * l2)).div(divisor.clone().square().scalarMul(2 * l1 * l1 * l2 * l2));
  var dphi1 = p1.clone().scalarMul(l2).sub(p2.clone().scalarMul(m1).mul(cosdif)).div(divisor.clone().scalarMul(l1 * l1 * l2));
  var dphi2 = p2.clone().scalarMul(l1 * (m1 + m2)).sub(p1.clone().mul(cosdif).scalarMul(m2 * l2)).div(divisor.clone().scalarMul(m2 * l1 * l2 * l2));
  var dp1 = h2.clone().mul(sindif).mul(cosdif).scalarMul(2).sub(h1).sub(phi1.clone().sin().scalarMul(g * l1 * (m1 + m2)));
  var dp2 = h1.clone().sub(h2.clone().mul(sindif).mul(cosdif).scalarMul(2)).sub(phi2.clone().sin().scalarMul(g * l2 * m2));
  return concat(dphi1, dp1, dphi2, dp2);
}

function convertToCoordinates(p, constants) {
  var l1 = constants[0],
      m1 = constants[1],
      l2 = constants[2],
      m2 = constants[3];
  var vectorLength = p.length / 4;
  var phi1 = p.slice(0, vectorLength),
      p1 = p.slice(vectorLength, vectorLength * 2),
      phi2 = p.slice(vectorLength * 2, vectorLength * 3),
      p2 = p.slice(vectorLength * 3, vectorLength * 4);
  var x1 = phi1.clone().sin().scalarMul(l1),
      y1 = phi1.clone().cos().scalarMul(l1);
  var x2 = x1.clone().add(phi2.clone().sin().scalarMul(l2)),
      y2 = y1.clone().add(phi2.clone().cos().scalarMul(l2));
  return [x1, y1, x2, y2].transpose();
}