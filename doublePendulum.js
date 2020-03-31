const tf2 = tf.scalar(2, 'float32'),
	tfd6 = tf.scalar(1/6, 'float32'),
	tfd2 = tf.scalar(1/2, 'float32');

const g = tf.scalar(9.81, 'float32');

function initPendulums(phi1, phi2, deviation, number){
	// Phi1 Array contains all the values of phi1 (no deviation)
	// Phi2 Array contains all the values of phi2 but deviated such that the 
	//  	average value is phi2, and the difference between any 2 phi2's is 
	//  	deviation
	phi1Array = tf.linspace(phi1, phi1, number);
	phi2Array = tf.linspace(phi2 - deviation/2*(number - 1), 
						    phi2 + deviaiton/2*(number - 1), 
							number);
	p1Array = tf.linspace(2, 4, number);
	p2Array = p1Array.clone();
	return tf.stack([phi1Array, p1Array, phi2Array, p2Array], axis=1);
}

function RK4TF (f, h, t, p, constants) {
	const k1 = h.mul(f(t, p, constants));
	const k2 = h.mul(f(t.add(h.mul(tfd2)), p.add(k1.mul(tfd2)), constants));
	const k3 = h.mul(f(t.add(h.mul(tfd2)), p.add(k2.mul(tfd2)), constants));
	const k4 = h.mul(f(t.add(h), p.add(k3), constants));
	return p.add(tfd6.mul(k1.add(k2).add(k3).add(k4)));
}

function derivativeTF (t, p, constants) {
	// p: phi1, p1, phi2, p2
	// constants: l1, m1, l2, m2 (all tf scalars)
	const l1 = constants[0], m1 = constants[1],
		l2 = constants[2], m2 = constants[3];
	
	const positionUnstacked = tf.unstack(p, 1);
	const phi1 = positionUnstacked[0], phi2 = positionUnstacked[2],
		p1 = positionUnstacked[1], p2 = positionUnstacked[3];
	
	const cosdif = phi1.sub(phi2).cos();
	const sindif = phi1.sub(phi2).sin();
	const divisor = m1.add(m2.mul(sindif.square()));
	
	const h1 = p1.mul(p2).mul(sindif)
		.div(l1.mul(l2).mul(divisor));
	const h2 = m2.mul(l2.square()).mul(p1.square())
		.add(m1.add(m2).mul(l1.square()).mul(p2.square()))
		.sub(m2.mul(l1).mul(l2).mul(p1).mul(p2).mul(tf2).mul(cosdif))
		.div(tf2.mul(l1.square()).mul(l2.square()).mul(divisor.square()));
	
	const dphi1 = 
		l2.mul(p1).sub(l1.mul(p2).mul(cosdif))
		.div(l1.square().mul(l2).mul(divisor));
	const dphi2 = m1.add(m2).mul(l2).mul(p2).sub(m2.mul(l2).mul(p1).mul(cosdif))
		.div(m2.mul(l1).mul(l2.square()).mul(divisor));
	
	const dp1 = h2.mul(sindif).mul(cosdif).mul(tf2).sub(h1)
		.sub(m1.add(m2).mul(l1).mul(phi1.sin()).mul(g));
	const dp2 = h1.sub(h2.mul(sindif))
		.sub(m2.mul(l2).mul(phi2.sin()).mul(g));
	
	return tf.stack([dphi1, dp1, dphi2, dp2], 1);
}
