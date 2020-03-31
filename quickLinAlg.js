Float64Array.prototype.add = (a, b) => {
	let destination = new Float64Array(a.length);
	destination.set(a);
	for (let i = 0; i < len(a); i++){
		destination[i] += b[i];
	}
	return destination;
};

Float64Array.prototype.mul = (a, b) => {
	let destination = new Float64Array(a.length);
	destination.set(a);
	for (let i = 0; i < len(a); i++){
		destination[i] *= b[i];
	}
	return destination;
}

Float64Array.prototype.div = (a, b) => {
	let destination = new Float64Array(a.length);
	destination.set(a);
	for (let i = 0; i < len(a); i++){
		destination[i] /= b[i];
	}
	return destination;
}

Float64Array.prototype.sub = (a, b) => {
	let destination = new Float64Array(a.length);
	destination.set(a);
	for (let i = 0; i < len(a); i++){
		destination[i] -= b[i];
	}
	return destination;
}

Float64Array.prototype.scalarAdd = (a, b) => {
	let destination = new Float64Array(a.length);
	destination.set(a);
	for (let i = 0; i < len(a); i++){
		destination[i] += b;
	}
	return destination;
};

Float64Array.prototype.scalarMul = (a, b) => {
	let destination = new Float64Array(a.length);
	destination.set(a);
	for (let i = 0; i < len(a); i++){
		destination[i] *= b;
	}
	return destination;
}

Float64Array.prototype.scalarDiv = (a, b) => {
	let destination = new Float64Array(a.length);
	destination.set(a);
	for (let i = 0; i < len(a); i++){
		destination[i] /= b;
	}
	return destination;
}

Float64Array.prototype.scalarSub = (a, b) => {
	let destination = new Float64Array(a.length);
	destination.set(a);
	for (let i = 0; i < len(a); i++){
		destination[i] -= b;
	}
	return destination;
}

Float64Array.prototype.square = (a) => {
	let destination = new Float64Array(a.length);
	destination.set(a);
	for (let i = 0; i < len(a); i++){
		destination[i] = destination[i]*destination[i];
	}
	return destination;
}

Float64Array.prototype.cos = (a) => {
	let destination = new Float64Array(a.length);
	destination.set(a);
	for (let i = 0; i < len(a); i++){
		destination[i] = Math.cos(destination[i]);
	}
	return destination;
}

Float64Array.prototype.sin = (a) => {
	let destination = new Float64Array(a.length);
	destination.set(a);
	for (let i = 0; i < len(a); i++){
		destination[i] = Math.sin(destination[i]);
	}
	return destination;
}

Float64Array.prototype.clone = (a) => {
	let destination = new Float64Array(a.length);
	destination.set(a);
	return destination;
}

Float64Array.prototype.concat = (a, b) => {
	let destination = new Float64Array(a.length + b.length);
	destination.set(a);
	destination.set(b, a.length);
	return destination;
}

// Float64Array.prototype.repeat = (times, axis=0) => {
	// let destination = new Float64Array(a.length + b.length);
	// destination.set(a);
	// destination.set(b, a.length);
	// return destination;
// }

function linspace (start, stop, number) {
	var arr = new Float64Array(number);
	var step = (stopValue - startValue) / (number - 1);
	for (var i = 0; i < number; i++) {
	arr[i] = startValue + (step * i);
	}
	return arr;
}