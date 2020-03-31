const g = 9.81

Float64Array.prototype.add = function (a) {
	let destination = new Float64Array(this.length);
	destination.set(this);
	for (let i = 0; i < this.length; i++){
		destination[i] += a[i];
	}
	return destination;
};

Float64Array.prototype.mul = function (a) {
	let destination = new Float64Array(this.length);
	destination.set(this);
	for (let i = 0; i < this.length; i++){
		destination[i] *= a[i];
	}
	return destination;
}

Float64Array.prototype.div = function (a) {
	let destination = new Float64Array(this.length);
	destination.set(this);
	for (let i = 0; i < this.length; i++){
		destination[i] /= a[i];
	}
	return destination;
}

Float64Array.prototype.sub = function (a) {
	let destination = new Float64Array(this.length);
	destination.set(this);
	for (let i = 0; i < this.length; i++){
		destination[i] -= a[i];
	}
	return destination;
}

Float64Array.prototype.scalarAdd = function (a) {
	let destination = new Float64Array(this.length);
	destination.set(this);
	for (let i = 0; i < this.length; i++){
		destination[i] += a;
	}
	return destination;
};

Float64Array.prototype.scalarMul = function (a) {
	let destination = new Float64Array(this.length);
	destination.set(this);
	for (let i = 0; i < this.length; i++){
		destination[i] *= a;
	}
	return destination;
}

Float64Array.prototype.scalarDiv = function (a) {
	let destination = new Float64Array(this.length);
	destination.set(this);
	for (let i = 0; i < this.length; i++){
		destination[i] /= a;
	}
	return destination;
}

Float64Array.prototype.scalarSub = function (a) {
	let destination = new Float64Array(this.length);
	destination.set(this);
	for (let i = 0; i < this.length; i++){
		destination[i] -= a;
	}
	return destination;
}

Float64Array.prototype.square = function () {
	let destination = new Float64Array(this.length);
	destination.set(this);
	for (let i = 0; i < this.length; i++){
		destination[i] = destination[i]*destination[i];
	}
	return destination;
}

Float64Array.prototype.cos = function () {
	let destination = new Float64Array(this.length);
	destination.set(this);
	for (let i = 0; i < this.length; i++){
		destination[i] = Math.cos(destination[i]);
	}
	return destination;
}

Float64Array.prototype.sin = function () {
	let destination = new Float64Array(this.length);
	destination.set(this);
	for (let i = 0; i < this.length; i++){
		destination[i] = Math.sin(destination[i]);
	}
	return destination;
}

Float64Array.prototype.concat = function (a) {
	let destination = new Float64Array(this.length + a.length);
	destination.set(this);
	destination.set(a, this.length);
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
	let step = (stop - start) / (number - 1);
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
	const phi2Array = linspace(phi2 - deviation/2*(number - 1), 
						    phi2 + deviation/2*(number - 1), 
							number);
	const p1Array = linspace(0, 0, number);
	const p2Array = p1Array.clone();
	const output = phi1Array.concat(p1Array).concat(phi2Array).concat(p2Array);
	return output;
}

function RK4LA (f, h, t, p, constants) {
	const k1 = f(t, p, constants).scalarMul(h);
	const k2 = f(t + h/2, p.add(k1.scalarMul(0.5)), constants).scalarMul(h);
	const k3 = f(t + h/2, p.add(k2.scalarMul(0.5)), constants).scalarMul(h);
	const k4 = f(t + h, p.add(k3), constants).scalarMul(h);
	return p.add(k1.add(k2).add(k3).add(k4).scalarMul(1/6));
}

function derivativeLA (t, p, constants) {
	// p: phi1, p1, phi2, p2
	// constants: l1, m1, l2, m2 (all tf scalars)
	const l1 = constants[0], m1 = constants[1],
		l2 = constants[2], m2 = constants[3];
	const vectorLength = p.length/4
	const phi1 = p.slice(0, vectorLength), p1 = p.slice(vectorLength, vectorLength*2)
		phi2 = p.slice(vectorLength*2, vectorLength*3), p2 = p.slice(vectorLength*3, vectorLength*4);
	
	const cosdif = phi1.sub(phi2).cos();
	const sindif = phi1.sub(phi2).sin();
	const divisor = sindif.square().scalarMul(m2).scalarAdd(m1);
	
	const h1 = p1.mul(p2).mul(sindif)
		.div(divisor.scalarMul(l1*l2));
	const h2 = p1.square().scalarMul(m2*l2*l2)
		.add(p2.square().scalarMul(l1*l1*(m1 + m2)))
		.sub(cosdif.mul(p1).mul(p2).scalarMul(2*m2*l2*l2))
		.div(divisor.square().scalarMul(2*l1*l1*l2*l2));
	
	const dphi1 = 
		p1.scalarMul(l2).sub(p2.scalarMul(m1).mul(cosdif))
		.div(divisor.scalarMul(l1*l1*l2));
	const dphi2 = p2.scalarMul(l1*(m1 + m2)).sub(p1.mul(cosdif).scalarMul(m2*l2))
		.div(divisor.scalarMul(m2*l1*l2*l2));
	
	const dp1 = h2.mul(sindif).mul(cosdif).scalarMul(2).sub(h1)
		.sub(phi1.sin().scalarMul(g*l1*(m1 + m2)));
	const dp2 = h1.sub(h2.mul(sindif).mul(cosdif).scalarMul(2))
		.sub(phi2.sin().scalarMul(g*l2*m2));
	
	return dphi1.concat(dp1).concat(dphi2).concat(dp2);
}

function convertToCoordinates(p, constants) {
	const l1 = constants[0], m1 = constants[1],
		l2 = constants[2], m2 = constants[3];
	const vectorLength = p.length/4
	const phi1 = p.slice(0, vectorLength), p1 = p.slice(vectorLength, vectorLength*2)
		phi2 = p.slice(vectorLength*2, vectorLength*3), p2 = p.slice(vectorLength*3, vectorLength*4);
	
	const x1 = phi1.sin().scalarMul(l1),
		y1 = phi1.cos().scalarMul(l1);
	
	const x2 = x1.add(phi2.sin().scalarMul(l2)),
		y2 = y1.add(phi2.cos().scalarMul(l2));
	return [x1, y1, x2, y2].transpose();
}
