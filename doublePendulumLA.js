const g = 9.81

Float64Array.prototype.add = function (a) {
	for (let i = 0; i < this.length; i++){
		this[i] += a[i];
	}
	return this;
};

Float64Array.prototype.mul = function (a) {
	for (let i = 0; i < this.length; i++){
		this[i] *= a[i];
	}
	return this;
}

Float64Array.prototype.div = function (a) {
	for (let i = 0; i < this.length; i++){
		this[i] /= a[i];
	}
	return this;
}

Float64Array.prototype.sub = function (a) {
	for (let i = 0; i < this.length; i++){
		this[i] -= a[i];
	}
	return this;
}

Float64Array.prototype.scalarAdd = function (a) {
	for (let i = 0; i < this.length; i++){
		this[i] += a;
	}
	return this;
};

Float64Array.prototype.scalarMul = function (a) {
	for (let i = 0; i < this.length; i++){
		this[i] *= a;
	}
	return this;
}

Float64Array.prototype.scalarDiv = function (a) {
	for (let i = 0; i < this.length; i++){
		this[i] /= a;
	}
	return this;
}

Float64Array.prototype.scalarSub = function (a) {
	for (let i = 0; i < this.length; i++){
		this[i] -= a;
	}
	return this;
}

Float64Array.prototype.square = function () {
	for (let i = 0; i < this.length; i++){
		this[i] = this[i]*this[i];
	}
	return this;
}

Float64Array.prototype.cos = function () {
	for (let i = 0; i < this.length; i++){
		this[i] = Math.cos(this[i]);
	}
	return this;
}

Float64Array.prototype.sin = function () {
	for (let i = 0; i < this.length; i++){
		this[i] = Math.sin(this[i]);
	}
	return this;
}

Float64Array.prototype.concat = function (a) {
	let destination = new Float64Array(this.length + a.length);
	destination.set(this);
	destination.set(a, this.length);
	return destination;
}

function concat (a, b, c, d) {
	let destination = new Float64Array(a.length + b.length + c.length + d.length);
	destination.set(a);
	destination.set(b, a.length);
	destination.set(c, a.length + b.length);
	destination.set(d, a.length + b.length + c.length);
	return destination;
}

Float64Array.prototype.clone = function () {
	let destination = new Float64Array(this.length);
	destination.set(this);
	return destination;
}

Float64Array.prototype.print = function() {
	console.log(this);
}

Array.prototype.transpose = function() {
    let newArray = new Array(this[0].length);
    for(var i = 0; i < this[0].length; i++){
		let rowArray = new Float64Array(this.length);
        for(var j = 0; j < this.length; j++){
            rowArray[j] = this[j][i];
        };
		newArray[i] = rowArray;
    };
    return newArray;
}

function linspace (start, stop, number) {
	let arr = new Float64Array(number);
	let step;
	if (number > 1) step = (stop - start) / (number - 1);
	else step = 0;
	for (let i = 0; i < number; i++) {
	arr[i] = start + (step * i);
	}
	return arr;
}

function initPendulums(phi1, phi2, deviation, number){
	// Phi1 Array contains all the values of phi1 (no deviation)
	// Phi2 Array contains all the values of phi2 but deviated such that the 
	//  	average value is phi2, and the difference between any 2 phi2's is 
	//  	deviation
	const phi1Array = linspace(phi1, phi1, number);
	let phi2Array;
	if (number > 1) {
		phi2Array = linspace(phi2 - deviation/2*(number - 1), 
						    phi2 + deviation/2*(number - 1), 
							number);
	}
	else {
		phi2Array = linspace(phi2 - deviation/2, 
						    phi2 + deviation/2, 
							number);
	}
	const p1Array = linspace(0, 0, number);
	const p2Array = p1Array.clone();
	const output = concat(phi1Array, p1Array, phi2Array, p2Array);
	return output;
}

function RK4LA (f, h, t, p, constants) {
	const k1 = f(t, p.clone(), constants).scalarMul(h);
	const k2 = f(t + h/2, p.clone().add(k1.clone().scalarMul(0.5)), constants).scalarMul(h);
	const k3 = f(t + h/2, p.clone().add(k2.clone().scalarMul(0.5)), constants).scalarMul(h);
	const k4 = f(t + h, p.clone().add(k3), constants).scalarMul(h);
	return p.clone().add(k1.clone().add(k2).add(k3).add(k4).scalarMul(1/6));
}

function derivativeLA (t, p, constants) {
	// p: phi1, p1, phi2, p2
	// constants: l1, m1, l2, m2
	const l1 = constants[0], m1 = constants[1],
		l2 = constants[2], m2 = constants[3];
	const vectorLength = p.length/4
	const phi1 = p.slice(0, vectorLength), p1 = p.slice(vectorLength, vectorLength*2),
		phi2 = p.slice(vectorLength*2, vectorLength*3), p2 = p.slice(vectorLength*3, vectorLength*4);
	
	const cosdif = phi1.clone().sub(phi2).cos();
	const sindif = phi1.clone().sub(phi2).sin();
	const divisor = sindif.clone().square().scalarMul(m2).scalarAdd(m1);
	
	const h1 = p1.clone().mul(p2).mul(sindif)
		.div(divisor.clone().scalarMul(l1*l2));
	const h2 = p1.clone().square().scalarMul(m2*l2*l2)
		.add(p2.clone().square().scalarMul(l1*l1*(m1 + m2)))
		.sub(cosdif.clone().mul(p1).mul(p2).scalarMul(2*m2*l2*l2))
		.div(divisor.clone().square().scalarMul(2*l1*l1*l2*l2));
	
	const dphi1 = 
		p1.clone().scalarMul(l2).sub(p2.clone().scalarMul(m1).mul(cosdif))
		.div(divisor.clone().scalarMul(l1*l1*l2));
	const dphi2 = p2.clone().scalarMul(l1*(m1 + m2)).sub(p1.clone().mul(cosdif).scalarMul(m2*l2))
		.div(divisor.clone().scalarMul(m2*l1*l2*l2));
	
	const dp1 = h2.clone().mul(sindif).mul(cosdif).scalarMul(2).sub(h1)
		.sub(phi1.clone().sin().scalarMul(g*l1*(m1 + m2)));
	const dp2 = h1.clone().sub(h2.clone().mul(sindif).mul(cosdif).scalarMul(2))
		.sub(phi2.clone().sin().scalarMul(g*l2*m2));
	
	return concat(dphi1, dp1, dphi2, dp2);
}

function convertToCoordinates(p, constants) {
	const l1 = constants[0], m1 = constants[1],
		l2 = constants[2], m2 = constants[3];
	const vectorLength = p.length/4;
	const phi1 = p.slice(0, vectorLength), p1 = p.slice(vectorLength, vectorLength*2),
		phi2 = p.slice(vectorLength*2, vectorLength*3), p2 = p.slice(vectorLength*3, vectorLength*4);
	
	const x1 = phi1.clone().sin().scalarMul(l1),
		y1 = phi1.clone().cos().scalarMul(l1);
	
	const x2 = x1.clone().add(phi2.clone().sin().scalarMul(l2)),
		y2 = y1.clone().add(phi2.clone().cos().scalarMul(l2));
	return [x1, y1, x2, y2].transpose();
}
