class DoublePendulumDemo {
	stepSize = 0.001;
	repeats = 32;
	randomness = .4;
	started = false;
	
	currentCoords;
	currentTime;
	currentCartCoords;
	width;
	height;

	origin;
	dots1;
	lines1;
	dots2;
	lines2;
	trails;
	
	continueLooping = false;
	restartable = false;
	colorscaleVariable = 0;
	currentStep = 0;
	lineResizer = 0.3;
	dotResizer = 0.1;
	
	lineGenerator = d3.line().curve(d3.curveNatural);
	
	constructor(constantsDict, animationLock) {
		this.isReset = true;
		this.lock = animationLock;
		this.constants = [constantsDict.l1, 
			constantsDict.m1,
			constantsDict.l2,
			constantsDict.m2];
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
		this.windowSize = 100/(2*1.05*(this.l1_ + this.l2_));
		this.dotsResizer = 0.1/this.windowSize;
		this.caption = constantsDict.caption;
		if ('colorScale' in constantsDict) {
			this.colorScale = constants.colorScale;
		}
		else {
			this.colorScale = d3.interpolateInferno;
		}
		this.currentTime = 0;
		this.initialCoords = this.getInitialCoords();
		this.lastCoords = this.initialCoords.clone();
		this.initialCoordsCart = convertToCoordinates(this.initialCoords, this.constants);

		this.lastTrailUpdateCoords = this.initialCoordsCart.slice(2);
		this.svg = d3.select(this.location)
			.append("svg")
			.attr("viewBox", "0 0 100 100");
		let explanation = "Hover to play. Click to restart";
		if (constantsDict.usingMobile) {
			explanation = "Tap to play. Double tap to restart";
		}
		let caption = d3.select(this.location)
			.append("p")
			.attr("class", "demoCaption")
			.html(this.caption);
		if (constantsDict.explain){
			caption.append("div")
				.attr("class", "demoSubcaption")
				.html(explanation);
		}
	}
	
	get l1_() {
		return this.constants[0];
	}
	get m1_() {
		return this.constants[1];
	}
	get l2_() {
		return this.constants[2];
	}
	get m2_() {
		return this.constants[3];
	}
	set l1_(newValue) {
		this.constants[0] = newValue;
		this.windowSize = 100/(2*1.05*(this.l1_ + this.l2_));
		this.dotsResizer = 0.1/this.windowSize;
		this.restart.bind(this)();
	}
	set m1_(newValue) {
		this.constants[1] = newValue;
		this.restart.bind(this)();
	}
	set l2_(newValue) {
		this.constants[2] = newValue;
		this.windowSize = 100/(2*1.05*(this.l1_ + this.l2_));
		this.dotsResizer = 0.1/this.windowSize;
		this.restart.bind(this)();
	}
	set m2_(newValue) {
		this.constants[3] = newValue;
		this.restart.bind(this)();
	}
	set pendulumNumber_(newValue) {
		this.pendulumNumber = newValue;
		this.initialCoords = this.getInitialCoords();
		this.currentCoords = this.initialCoords.clone();
		this.currentCartCoords = convertToCoordinates(this.currentCoords, this.constants);
		this.initialCoordsCart = convertToCoordinates(this.initialCoords, this.constants);
		// this.restart.bind(this)();
	}
	set phi1Init_(newValue) {
		this.phi1Init = newValue;
		// this.restart.bind(this)();
	}
	set phi2Init_(newValue) {
		this.phi2Init = newValue;
		// this.restart.bind(this)();
	}
	set trailLength_ (newValue) {
		this.trailLength = newValue;
		// this.restart.bind(this)();
	}
	get pendulumNumber_() {
		return this.pendulumNumber;
	}
	get phi1Init_() {
		return this.phi1Init;
	}
	get phi2Init_() {
		return this.phi2Init;
	}
	get trailLength_() {
		return this.trailLength;
	}
	
	
	start() {
		if (this.lock.request()) {
			return;
		}
		if (! this.continueLooping ) {
			this.isReset = false;
			this.continueLooping = true;
			const demoContext = this;
			requestAnimationFrame(()=>demoContext.updateDisplay());
		}
	}

	stop() {
		if (this.continueLooping) this.continueLooping = false;
	}
	
	restart() {
		this.stop();
		for (let i = 0; i < this.pendulumNumber; i++) {
			this.trailPaths[i] = "M " + (this.initialCoordsCart[i][2]*this.windowSize).toPrecision(5) + " " + (this.initialCoordsCart[i][3]*this.windowSize).toPrecision(5) + " ";
		}
		this.currentCoords = this.initialCoords;
		this.lastCoords = this.currentCoords.clone();
		this.currentCartCoords = this.initialCoordsCart;
		this.drawPendulumsAt(this.initialCoordsCart);
		this.isReset = true;
	}

	getInitialCoords() {
		return initPendulums(this.phi1Init + (Math.random() - 0.5)*this.randomness, 
			this.phi2Init + (Math.random() - 0.5)*this.randomness, 
			this.deviation, 
			this.pendulumNumber);
	}

	init() {
		for (let i = 0; i < this.pendulumNumber; i++) {
			this.trailPaths[i] = "M " + (this.initialCoordsCart[i][2]*this.windowSize).toPrecision(5) + " " + (this.initialCoordsCart[i][3]*this.windowSize).toPrecision(5) + " ";
		}
		this.origin = this.svg.append('g')
			.attr("transform", "translate(" + 50 + "," + 50 + ")")
			.attr("class", "origin");

		this.trails = this.origin.selectAll(".trails")
			.data(this.initialCoordsCart).enter()
			.append('path')
			.attr("fill", "transparent")
			.attr("class", "trails")
			.attr("stroke-width", this.lineResizer)
			.attr("d", (d)=>  { return "M " + (d[2]*this.windowSize).toPrecision(5) + " " + (d[3]*this.windowSize).toPrecision(5) + " ";} )
			.attr("stroke", this.colorScaleFunc.bind(this));

		let pendulums = this.origin.selectAll('.pendulums')
			.data(this.initialCoordsCart).enter()
			.append('g').attr("class", "pendulums");
		
			
		this.lines1 = pendulums
			.append('line')
			.attr("class", "pendulumLine1")
			.attr("x1", 0)
			.attr("y1", 0)
			.attr("stroke-width", 2*this.lineResizer)
			.attr("x2", (d)=> { return d[0]*this.windowSize;})
			.attr("y2", (d)=>  { return d[1]*this.windowSize;})
			.attr("stroke", this.colorScaleFunc.bind(this));
			
		this.lines2 = pendulums
			.append('line')
			.attr("class", "pendulumLine2")
			.attr("stroke-width", 2*this.lineResizer)
			.attr("x1", (d)=>  { return d[0]*this.windowSize;})
			.attr("y1", (d)=>  { return d[1]*this.windowSize;})
			.attr("x2", (d)=>  { return d[2]*this.windowSize;})
			.attr("y2", (d)=>  { return d[3]*this.windowSize;})
			.attr("stroke", this.colorScaleFunc.bind(this));

		this.dots1 = pendulums
			.append('circle')
			.attr("class", "dots1")
			.attr("fill", this.colorScaleFunc.bind(this))
			.attr("cx", (d)=> { return d[0]*this.windowSize;})
			.attr("cy", (d)=> { return d[1]*this.windowSize;})
			.attr("r", this.constants[1]*this.dotResizer*this.windowSize);
			
		this.dots2 = pendulums
			.append('circle')
			.attr("class", "dots2")
			.attr("fill", this.colorScaleFunc.bind(this))
			.attr("cx", (d)=>  { return d[2]*this.windowSize;})
			.attr("cy", (d)=>  { return d[3]*this.windowSize;})
			.attr("r", this.constants[3]*this.dotResizer*this.windowSize);

		this.origin.append("circle")
			.attr("class", "origin-mark")
			.attr("r", Math.max(this.constants[2], this.constants[1])*1.5*this.dotResizer*this.windowSize)
			.style("fill", this.colorScaleFunc.bind(this)(null, (this.pendulumNumber - 1)/(this.pendulumNumber + 5)));
				
		this.currentCoords = this.initialCoords.clone();
	}
	
	colorScaleFunc (d, i) {
		if (this.pendulumNumber == 1) {
			if (d != null) return this.colorScale(.95);
			return d3.rgb(this.colorScale(.95)).hex();
		}
		if (d == null) return d3.rgb(this.colorScale(.95)).hex()
		return this.colorScale(i/(this.pendulumNumber + 5) + 2/this.pendulumNumber);
	}

	updateTrailDataWith(newCoords) {
		this.trailCounter++;
		if (this.trailLength == 0) return;
		if (this.trailCounter > this.trailLength) {
			for ( let i = 0; i < this.pendulumNumber; i++) {
				this.trailPaths[i] = "M" + this.trailPaths[i].slice(this.trailPaths[i].indexOf("L") + 1) + "L " + (newCoords[i][2]*this.windowSize).toPrecision(5) + " " + (newCoords[i][3]*this.windowSize).toPrecision(5) + " ";
			}
		}
		else {
			for ( let i = 0; i < this.pendulumNumber; i++) {
				this.trailPaths[i] += "L " + (newCoords[i][2]*this.windowSize).toPrecision(5) + " " + (newCoords[i][3]*this.windowSize).toPrecision(5) + " ";
			}
		}
		this.lastCartCoords = this.newCoords;
	}

	drawPendulumsAt(newCoords) {
		this.svg.selectAll("*").remove();
		this.initialCoords = this.getInitialCoords();
		this.initialCoordsCart = convertToCoordinates(this.initialCoords, this.constants);
		this.trailCounter = 0;
		this.init();
	}


	updateDisplay() {
		this.lastCoords = this.currentCoords.clone();
		for (let i = 0; i < this.repeats; i++) {
			this.currentCoords = RK4LA(derivativeLA, this.stepSize, this.currentTime, this.currentCoords, this.constants);
		}
		this.currentTime += this.stepSize*this.repeats;
		this.currentCartCoords = convertToCoordinates(this.currentCoords, this.constants);

		this.dots1.data(this.currentCartCoords)
			.attr("cx", (d) => { return d[0]*this.windowSize;})
			.attr("cy", (d) => { return d[1]*this.windowSize;});
		this.dots2.data(this.currentCartCoords)
			.attr("cx", (d) => { return d[2]*this.windowSize;})
			.attr("cy", (d) => { return d[3]*this.windowSize;});
			
		this.lines1.data(this.currentCartCoords)
			.attr("x2", (d) => { return d[0]*this.windowSize;})
			.attr("y2", (d) => { return d[1]*this.windowSize;});
		
		this.lines2.data(this.currentCartCoords)
			.attr("x1", (d) => { return d[0]*this.windowSize;})
			.attr("y1", (d) => { return d[1]*this.windowSize;})
			.attr("x2", (d) => { return d[2]*this.windowSize;})
			.attr("y2", (d) => { return d[3]*this.windowSize;});
		
		if (this.useTrails&&!(this.currentStep%this.trailUpdateInterval)) {
			this.lastTrailUpdateCoords = this.currentCartCoords.slice(2);
			this.updateTrailDataWith(this.currentCartCoords);
			this.trails.data(this.trailPaths)
				.attr("d", (d)=> { return d;} );
		}
		this.currentStep += 1;
		if (this.continueLooping) {
			requestAnimationFrame(()=>this.updateDisplay());
		}
		else {
			this.lock.release();
		}
	}
}