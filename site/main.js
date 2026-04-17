// Port of experiments/wireless_k3_opt2.py simulator for the A_3 algorithm,
// plus a canvas renderer.

const TWO_PI = 2 * Math.PI;
const Y_OPT = 1.2158578321292429;
const ROBOT_COLORS = ["#2a9d8f", "#e76f51", "#264653"];
const STAR_COLOR = "#c0392b";

// --- Phase construction ---------------------------------------------------

function simplePhases(a, L, s) {
  return { approachStart: [0], a: [a], L: [L], s: [s] };
}

function twoPhases(a, La, sa, b, Lb, sb) {
  return {
    approachStart: [0, 2 + La],
    a: [a, b],
    L: [La, Lb],
    s: [sa, sb],
  };
}

function a3Robots(y) {
  const L = Math.PI - y / 2;
  const aB = TWO_PI - y;
  const aRedeploy = Math.PI - y / 2;
  return [
    simplePhases(0, L, +1),                          // r_1
    simplePhases(aB, L, -1),                         // r_2
    twoPhases(aB, y, +1, aRedeploy, 0, +1),          // r_3
  ];
}

// --- Position / hit-time --------------------------------------------------

function position(phases, t) {
  const K = phases.approachStart.length;
  for (let k = 0; k < K; k++) {
    const ts = phases.approachStart[k];
    const tb = ts + 1;
    const te = tb + phases.L[k];
    const ak = phases.a[k];
    const sk = phases.s[k];
    const Lk = phases.L[k];

    if (t >= ts && t <= tb) {
      const frac = t - ts;
      return [frac * Math.cos(ak), frac * Math.sin(ak)];
    }
    if (t > tb && t <= te) {
      const u = t - tb;
      const ang = ak + sk * u;
      return [Math.cos(ang), Math.sin(ang)];
    }
    if (k + 1 < K && t > te && t <= te + 1) {
      // transit boundary endpoint -> origin (linear)
      const endAng = ak + sk * Lk;
      const ex = Math.cos(endAng);
      const ey = Math.sin(endAng);
      const frac = t - te;
      return [(1 - frac) * ex, (1 - frac) * ey];
    }
  }
  // past last phase -> parked at last scan endpoint
  const last = K - 1;
  const endAng = phases.a[last] + phases.s[last] * phases.L[last];
  return [Math.cos(endAng), Math.sin(endAng)];
}

function hitTime(phases, theta) {
  let best = Infinity;
  for (let k = 0; k < phases.approachStart.length; k++) {
    const ak = phases.a[k];
    const sk = phases.s[k];
    const Lk = phases.L[k];
    let delta = ((theta - ak) * sk) % TWO_PI;
    if (delta < 0) delta += TWO_PI;
    if (delta <= Lk + 1e-12) {
      const cand = phases.approachStart[k] + 1 + delta;
      if (cand < best) best = cand;
    }
  }
  return best;
}

function evacForTheta(robots, theta) {
  const hits = robots.map((r) => hitTime(r, theta));
  const tFound = Math.min(...hits);
  const finder = hits.indexOf(tFound);
  if (!isFinite(tFound)) {
    return { evac: Infinity, tFound: Infinity, finder: -1, bottleneck: -1 };
  }
  const exitX = Math.cos(theta);
  const exitY = Math.sin(theta);
  let maxChord = 0;
  let bot = 0;
  robots.forEach((r, i) => {
    const [x, y] = position(r, tFound);
    const d = Math.hypot(x - exitX, y - exitY);
    if (d > maxChord) {
      maxChord = d;
      bot = i;
    }
  });
  return { evac: tFound + maxChord, tFound, finder, bottleneck: bot };
}

function worstCase(robots, n = 2001) {
  let worstE = -Infinity;
  let worstT = 0;
  for (let i = 0; i < n; i++) {
    const theta = (i / n) * TWO_PI;
    const { evac } = evacForTheta(robots, theta);
    if (evac > worstE) {
      worstE = evac;
      worstT = theta;
    }
  }
  return { worstEvac: worstE, worstTheta: worstT };
}

// --- Canvas rendering -----------------------------------------------------

const canvas = document.getElementById("viz");
const ctx = canvas.getContext("2d");
const W = canvas.width;
const H = canvas.height;
const cx = W / 2;
const cy = H / 2;
const R = Math.min(W, H) / 2 - 30;

function toScreen([x, y]) {
  return [cx + x * R, cy - y * R];
}

function drawDisk() {
  ctx.strokeStyle = "#333";
  ctx.lineWidth = 1.2;
  ctx.beginPath();
  ctx.arc(cx, cy, R, 0, TWO_PI);
  ctx.stroke();

  // origin marker
  ctx.fillStyle = "#555";
  ctx.beginPath();
  ctx.arc(cx, cy, 3, 0, TWO_PI);
  ctx.fill();
}

function drawTrajectory(phases, color, tMax = 5) {
  ctx.strokeStyle = color;
  ctx.lineWidth = 1.5;
  ctx.globalAlpha = 0.5;
  ctx.beginPath();
  const dt = 0.01;
  const first = toScreen(position(phases, 0));
  ctx.moveTo(first[0], first[1]);
  for (let t = dt; t <= tMax; t += dt) {
    const sc = toScreen(position(phases, t));
    ctx.lineTo(sc[0], sc[1]);
  }
  ctx.stroke();
  ctx.globalAlpha = 1;
}

function drawRobotDot(phases, t, color) {
  const [x, y] = position(phases, t);
  const [sx, sy] = toScreen([x, y]);
  ctx.fillStyle = color;
  ctx.beginPath();
  ctx.arc(sx, sy, 7, 0, TWO_PI);
  ctx.fill();
  ctx.strokeStyle = "white";
  ctx.lineWidth = 1.5;
  ctx.stroke();
}

function drawExitStar(theta) {
  const [sx, sy] = toScreen([Math.cos(theta), Math.sin(theta)]);
  ctx.fillStyle = STAR_COLOR;
  ctx.beginPath();
  const outer = 10;
  const inner = 4;
  for (let i = 0; i < 10; i++) {
    const a = -Math.PI / 2 + (i * Math.PI) / 5;
    const r = i % 2 === 0 ? outer : inner;
    const x = sx + r * Math.cos(a);
    const y = sy + r * Math.sin(a);
    if (i === 0) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
  }
  ctx.closePath();
  ctx.fill();
}

function drawChord(fromXY, toXY, color) {
  const [fx, fy] = toScreen(fromXY);
  const [tx, ty] = toScreen(toXY);
  ctx.strokeStyle = color;
  ctx.setLineDash([5, 4]);
  ctx.lineWidth = 1.8;
  ctx.beginPath();
  ctx.moveTo(fx, fy);
  ctx.lineTo(tx, ty);
  ctx.stroke();
  ctx.setLineDash([]);
}

// --- State + UI wiring -----------------------------------------------------

const yEl = document.getElementById("y-slider");
const thetaEl = document.getElementById("theta-slider");
const yValEl = document.getElementById("y-val");
const thetaValEl = document.getElementById("theta-val");
const evacNowEl = document.getElementById("evac-now");
const evacWorstEl = document.getElementById("evac-worst");
const tFoundEl = document.getElementById("t-found");
const bottleneckEl = document.getElementById("bottleneck");

let y = Y_OPT;
let theta = 1.486472;

function render() {
  const robots = a3Robots(y);
  const { evac, tFound, bottleneck } = evacForTheta(robots, theta);
  const { worstEvac } = worstCase(robots);

  yValEl.textContent = y.toFixed(4);
  thetaValEl.textContent = theta.toFixed(4);
  evacNowEl.textContent = evac.toFixed(6);
  evacWorstEl.textContent = worstEvac.toFixed(6);
  tFoundEl.textContent = tFound.toFixed(4);
  bottleneckEl.textContent = "r_" + (bottleneck + 1);

  ctx.clearRect(0, 0, W, H);
  drawDisk();
  robots.forEach((r, i) => drawTrajectory(r, ROBOT_COLORS[i]));
  robots.forEach((r, i) => drawRobotDot(r, tFound, ROBOT_COLORS[i]));
  drawExitStar(theta);

  // Chord from the bottleneck robot to the exit.
  const bottleneckPos = position(robots[bottleneck], tFound);
  const exitXY = [Math.cos(theta), Math.sin(theta)];
  drawChord(bottleneckPos, exitXY, ROBOT_COLORS[bottleneck]);
}

yEl.addEventListener("input", (e) => {
  y = +e.target.value;
  render();
});
thetaEl.addEventListener("input", (e) => {
  theta = +e.target.value;
  render();
});
document.getElementById("jump-worst").addEventListener("click", () => {
  const robots = a3Robots(y);
  const { worstTheta } = worstCase(robots);
  theta = worstTheta;
  thetaEl.value = worstTheta;
  render();
});
document.getElementById("reset-y").addEventListener("click", () => {
  y = Y_OPT;
  yEl.value = Y_OPT;
  render();
});

// Click on the boundary/disk to set theta.
canvas.addEventListener("click", (e) => {
  const rect = canvas.getBoundingClientRect();
  const scaleX = canvas.width / rect.width;
  const scaleY = canvas.height / rect.height;
  const x = (e.clientX - rect.left) * scaleX - cx;
  const yPix = -((e.clientY - rect.top) * scaleY - cy);
  let a = Math.atan2(yPix, x);
  if (a < 0) a += TWO_PI;
  theta = a;
  thetaEl.value = theta.toFixed(3);
  render();
});

render();
