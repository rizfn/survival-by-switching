// ----------------------------
// Environment switching
// ----------------------------
function Kenv(t, freq, startInA) {
  if (freq <= 0) return startInA ? 1 : 2;
  const T = 1 / freq;
  const phase = (t % T) / T;
  return phase < 0.5 ? (startInA ? 1 : 2) : (startInA ? 2 : 1);
}

// ----------------------------
// Ricker consumption F_E(x) = x e^{-alpha_E x}
// ----------------------------
function FE(x, alphaE) {
  return x * Math.exp(-alphaE * x);
}

// ----------------------------
// RK4 integrator
// ----------------------------
function rk4Step(x, y, t, dt, params) {
  const f = (x, y, t) => {
    const envA = Kenv(t, params.freq, params.startInA) === 1;

    const K = envA ? params.k1 : params.k2;
    const alphaE = envA ? params.alpha1 : params.alpha2;
    const gamma = 1;

    const Fx = FE(x, alphaE);

    const dx = params.r * x * (1 - x / K) - Fx * y;
    const dy = y * (Fx - gamma);

    return [dx, dy];
  };

  const [stage1x, stage1y] = f(x, y, t);
  const [stage2x, stage2y] = f(x + 0.5*dt*stage1x, y + 0.5*dt*stage1y, t + 0.5*dt);
  const [stage3x, stage3y] = f(x + 0.5*dt*stage2x, y + 0.5*dt*stage2y, t + 0.5*dt);
  const [stage4x, stage4y] = f(x + dt*stage3x, y + dt*stage3y, t + dt);

  const xn = x + (dt/6)*(stage1x + 2*stage2x + 2*stage3x + stage4x);
  const yn = y + (dt/6)*(stage1y + 2*stage2y + 2*stage3y + stage4y);

  return [xn, yn];
}

// ----------------------------
// Simulation
// ----------------------------
function runSimulation() {
  const params = {
    freq: parseFloat(document.getElementById("freq").value),
    r: parseFloat(document.getElementById("r").value),
    alpha1: parseFloat(document.getElementById("alpha1").value),
    alpha2: parseFloat(document.getElementById("alpha2").value),
    k1: parseFloat(document.getElementById("k1").value),
    k2: parseFloat(document.getElementById("k2").value),
    startInA: document.getElementById("startInA").checked
  };

  let x = parseFloat(document.getElementById("x0").value);
  let y = parseFloat(document.getElementById("y0").value);

  const Ttotal = parseFloat(document.getElementById("Ttotal").value);
  const dt = 0.01;

  const ts = [], xs = [], ys = [];

  let t = 0;
  while (t < Ttotal) {
    ts.push(t);
    xs.push(x);
    ys.push(y);

    [x, y] = rk4Step(x, y, t, dt, params);
    t += dt;
  }

  drawTimeSeries(ts, xs, ys);
  drawPhase(xs, ys);
}

// ----------------------------
// Plotting
// ----------------------------
function drawTimeSeries(ts, xs, ys) {
  const svg = d3.select("#timeseries");
  svg.selectAll("*").remove();

  const width = +svg.attr("width");
  const height = +svg.attr("height");
  const margin = { top: 20, right: 20, bottom: 30, left: 40 };
  const innerW = width - margin.left - margin.right;
  const innerH = height - margin.top - margin.bottom;

  const g = svg.append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);

  const xScale = d3.scaleLinear().domain(d3.extent(ts)).range([0, innerW]);
  const yScale = d3.scaleLinear()
    .domain([Math.min(d3.min(xs), d3.min(ys), 0),
             Math.max(d3.max(xs), d3.max(ys))])
    .nice()
    .range([innerH, 0]);

  g.append("g").attr("transform", `translate(0,${innerH})`).call(d3.axisBottom(xScale));
  g.append("g").call(d3.axisLeft(yScale));

  const line = (arr, color) =>
    g.append("path")
      .datum(arr)
      .attr("fill", "none")
      .attr("stroke", color)
      .attr("stroke-width", 1.5)
      .attr("d", d3.line()
        .x((d,i)=>xScale(ts[i]))
        .y((d,i)=>yScale(arr[i])));

  line(xs, "steelblue");
  line(ys, "tomato");
}

function drawPhase(xs, ys) {
  const svg = d3.select("#phase");
  svg.selectAll("*").remove();

  const width = +svg.attr("width");
  const height = +svg.attr("height");
  const margin = { top: 20, right: 20, bottom: 30, left: 40 };
  const innerW = width - margin.left - margin.right;
  const innerH = height - margin.top - margin.bottom;

  const g = svg.append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);

  const xScale = d3.scaleLinear().domain(d3.extent(xs)).nice().range([0, innerW]);
  const yScale = d3.scaleLinear().domain(d3.extent(ys)).nice().range([innerH, 0]);

  g.append("g").attr("transform", `translate(0,${innerH})`).call(d3.axisBottom(xScale));
  g.append("g").call(d3.axisLeft(yScale));

  g.append("path")
    .datum(xs)
    .attr("fill", "none")
    .attr("stroke", "gray")
    .attr("stroke-width", 1.2)
    .attr("d", d3.line()
      .x((d,i)=>xScale(xs[i]))
      .y((d,i)=>yScale(ys[i])));
}

// ----------------------------
// Events
// ----------------------------
["freq","r","alpha1","alpha2","k1","k2","x0","y0","Ttotal","startInA"]
  .forEach(id => document.getElementById(id).addEventListener("input", runSimulation));

document.getElementById("runBtn").addEventListener("click", runSimulation);

window.addEventListener("load", runSimulation);
