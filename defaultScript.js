class AnimationLock {
	locked;
	constructor(){ 
		this.locked = false;
	}
	request() {
		if (this.locked) {
			return false;
		}
		this.locked = true;
	}
	release() {
		this.locked = false;
	}
}

function mobileCheck() {
  var check = false;
  (function(a){if(/(android|bb\d+|meego).+mobile|avantgo|bada\/|blackberry|blazer|compal|elaine|fennec|hiptop|iemobile|ip(hone|od)|iris|kindle|lge |maemo|midp|mmp|mobile.+firefox|netfront|opera m(ob|in)i|palm( os)?|phone|p(ixi|re)\/|plucker|pocket|psp|series(4|6)0|symbian|treo|up\.(browser|link)|vodafone|wap|windows ce|xda|xiino|android|ipad|playbook|silk/i.test(a)||/1207|6310|6590|3gso|4thp|50[1-6]i|770s|802s|a wa|abac|ac(er|oo|s\-)|ai(ko|rn)|al(av|ca|co)|amoi|an(ex|ny|yw)|aptu|ar(ch|go)|as(te|us)|attw|au(di|\-m|r |s )|avan|be(ck|ll|nq)|bi(lb|rd)|bl(ac|az)|br(e|v)w|bumb|bw\-(n|u)|c55\/|capi|ccwa|cdm\-|cell|chtm|cldc|cmd\-|co(mp|nd)|craw|da(it|ll|ng)|dbte|dc\-s|devi|dica|dmob|do(c|p)o|ds(12|\-d)|el(49|ai)|em(l2|ul)|er(ic|k0)|esl8|ez([4-7]0|os|wa|ze)|fetc|fly(\-|_)|g1 u|g560|gene|gf\-5|g\-mo|go(\.w|od)|gr(ad|un)|haie|hcit|hd\-(m|p|t)|hei\-|hi(pt|ta)|hp( i|ip)|hs\-c|ht(c(\-| |_|a|g|p|s|t)|tp)|hu(aw|tc)|i\-(20|go|ma)|i230|iac( |\-|\/)|ibro|idea|ig01|ikom|im1k|inno|ipaq|iris|ja(t|v)a|jbro|jemu|jigs|kddi|keji|kgt( |\/)|klon|kpt |kwc\-|kyo(c|k)|le(no|xi)|lg( g|\/(k|l|u)|50|54|\-[a-w])|libw|lynx|m1\-w|m3ga|m50\/|ma(te|ui|xo)|mc(01|21|ca)|m\-cr|me(rc|ri)|mi(o8|oa|ts)|mmef|mo(01|02|bi|de|do|t(\-| |o|v)|zz)|mt(50|p1|v )|mwbp|mywa|n10[0-2]|n20[2-3]|n30(0|2)|n50(0|2|5)|n7(0(0|1)|10)|ne((c|m)\-|on|tf|wf|wg|wt)|nok(6|i)|nzph|o2im|op(ti|wv)|oran|owg1|p800|pan(a|d|t)|pdxg|pg(13|\-([1-8]|c))|phil|pire|pl(ay|uc)|pn\-2|po(ck|rt|se)|prox|psio|pt\-g|qa\-a|qc(07|12|21|32|60|\-[2-7]|i\-)|qtek|r380|r600|raks|rim9|ro(ve|zo)|s55\/|sa(ge|ma|mm|ms|ny|va)|sc(01|h\-|oo|p\-)|sdk\/|se(c(\-|0|1)|47|mc|nd|ri)|sgh\-|shar|sie(\-|m)|sk\-0|sl(45|id)|sm(al|ar|b3|it|t5)|so(ft|ny)|sp(01|h\-|v\-|v )|sy(01|mb)|t2(18|50)|t6(00|10|18)|ta(gt|lk)|tcl\-|tdg\-|tel(i|m)|tim\-|t\-mo|to(pl|sh)|ts(70|m\-|m3|m5)|tx\-9|up(\.b|g1|si)|utst|v400|v750|veri|vi(rg|te)|vk(40|5[0-3]|\-v)|vm40|voda|vulc|vx(52|53|60|61|70|80|81|83|85|98)|w3c(\-| )|webc|whit|wi(g |nc|nw)|wmlb|wonu|x700|yas\-|your|zeto|zte\-/i.test(a.substr(0,4))) check = true;})(navigator.userAgent||navigator.vendor||window.opera);
  return check;
};

const usingMobile = mobileCheck();
let goodColorScales = [d3.interpolateOranges, d3.interpolateInferno, d3.interpolatePlasma];

let constants1 = {
	"l1": 8, "m1": 5, "l2": 6, "m2": 5, 
	"phi1Init": 0.5, "phi2Init": 0.5, 
	"pendulumNumber": 1, "deviation": 0, 
	"trailUpdateInterval": 3, "trailLength": 60,
	"location": "#demoSmall", "trails": false, 
	"colorscale": goodColorScales[1], 
	"explain": true, "usingMobile": usingMobile,
	"caption": "A double pendulum with a small oscillation"
};

let constants2 = JSON.parse(JSON.stringify(constants1));
constants2.location = "#demoVanilla";
constants2.pendulumNumber = 60;
constants2.explain = false;
constants2.phi1Init = 1.5;
constants2.phi2Init = 2.7;
constants2.deviation = 0.0005;
constants2.colorscale = goodColorScales[1];
constants2.caption = "60 pendulums 0.03 degrees apart initially";

let constants3 = JSON.parse(JSON.stringify(constants2));
constants3.trails = true;
constants3.location = "#demoTrails";
constants3.caption = "50 pendulums with their paths traced out";
	
let constants4 = JSON.parse(JSON.stringify(constants3));
constants4.location = "#demoPlayable";
constants4.trailLength = 0;
constants4.caption = "Change the parameters below to adjust the setup<br>" + 
					 "Note: more pendulums and longer paths will cause the demo to run slower.";

let lock = new AnimationLock();
let demo1 = new DoublePendulumDemo(constants1, lock);
let demo2 = new DoublePendulumDemo(constants2, lock);
let demo3 = new DoublePendulumDemo(constants3, lock);
let demo4 = new DoublePendulumDemo(constants4, lock);


let demos = [demo1, demo2, demo3, demo4]
demos.forEach((demo)=>{
	demo.init();
	const demoLocation = document.getElementById(demo.location.slice(1));
	if (!usingMobile) {
		demoLocation.onmouseenter = demo.start.bind(demo);
		demoLocation.onmouseleave = demo.stop.bind(demo);
		demoLocation.onclick = demo.restart.bind(demo);
	}
	else {
		demoLocation.onclick = function (){
			demos.forEach ((demo) => {
				if (demo.continueLooping) {
					demo.stop.bind(demo)();
				}
			});
			demo.start.bind(demo)();
		};
		demoLocation.ondblclick = demo.restart.bind(demo);
	}
});

// Deviation: Slider 0, 1
// Colorscale: Button (iterating through choices
class Slider {
	constructor(name, parameter, startValue, stopValue, isInteger, initialValue) {
		this.name = name;
		this.startValue = startValue;
		this.stopValue = stopValue;
		this.isInteger = isInteger;
		this.parameter = parameter;
		this.initialValue = initialValue;
	}
	
	publish (containerID, number, total) {
		let container = document.getElementById(containerID);
		let fullDiv = document.createElement('div');
		fullDiv.className = "sliderWrapper";
		let textHolder = document.createElement('div');
		textHolder.innerHTML = (this.name + ": " + '<span style="font-weight: bold">' +  this.initialValue + "</span>");
		let slider = document.createElement("input");
		slider.setAttribute("type", "range");
		slider.setAttribute("min", this.startValue);
		slider.setAttribute("max", this.stopValue);
		slider.setAttribute("value", this.initialValue);
		
		slider.className = "slider";
		slider.id = "input" + number;
		
		if (this.isInteger) {slider.setAttribute("step", 1);}
		else {slider.setAttribute("step", 0.1);}
		
		slider.oninput = (d)=>{
			textHolder.innerHTML = this.name + ": " + '<span style="font-weight: bold">' +  d.target.value + "</span>";
			demo4[this.parameter + "_"] = d.target.valueAsNumber;
		};
		slider.onmouseup = demo4.restart.bind(demo4);
		fullDiv.appendChild(slider);
		fullDiv.appendChild(textHolder);
		container.append(fullDiv);
	}
}

let slidersConsts = [["Pendulum Number", "pendulumNumber", 1, 200, true, 50], 
	["Trail Length", "trailLength", 0, 200, true, 0], 
	["Rod 1 Length", "l1", 3, 10, false, 5],
	["Rod 2 Length", "l2", 3, 10, false, 5],
	["Mass 1", "m1", 1, 10, false, 5],
	["Mass 2", "m2", 1, 10, false, 5],
	["Phi 1 Initial", "phi1Init", 0, 2*Math.PI, false, 1.5],
	["Phi 2 Initial", "phi2Init", 0, 2*Math.PI, false, 2.5]];

slidersConsts.forEach((item, i)=>{
	let slider = new Slider(...item);
	if ( i%2 == 1 ){
		slider.publish("rightCol", i, slidersConsts.length);
	}
	else {
		slider.publish("leftCol", i, slidersConsts.length);
	}
});