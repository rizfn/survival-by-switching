// ----------------------------
// Plotting
// ----------------------------
function drawPhasePlane(svgId, r, K, gamma, options) {
  const svg = d3.select(svgId);
  svg.selectAll("*").remove();

  const markerId = `arrow-${svgId.substring(1)}`;
  svg.append("defs").append("marker")
      .attr("id", markerId)
      .attr("viewBox", "0 -5 10 10")
      .attr("refX", 8)
      .attr("refY", 0)
      .attr("markerWidth", 5)
      .attr("markerHeight", 5)
      .attr("orient", "auto")
    .append("path")
      .attr("d", "M0,-4L10,0L0,4")
      .attr("fill", "#aaa");

  const width = +svg.attr("width");
  const height = +svg.attr("height");
  const margin = { top: 20, right: 20, bottom: 40, left: 50 };
  const innerW = width - margin.left - margin.right;
  const innerH = height - margin.top - margin.bottom;

  const xMin = -0.1;
  const xMax = 2.1;
  const yMin = -0.1;
  const yMax = 1.1;

  const xNullData = [];
    for (let x = xMin; x <= xMax + 0.05; x += 0.05) {
      xNullData.push([x, r * (1 - x / K)]);
  }

  const xScale = d3.scaleLinear().domain([xMin, xMax]).range([0, innerW]);
  const yScale = d3.scaleLinear().domain([yMin, yMax]).range([innerH, 0]);

  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

  // Draw vector field
  const gridSpacing = 25;
  for (let px = 0; px <= innerW; px += gridSpacing) {
      for (let py = 0; py <= innerH; py += gridSpacing) {
          const x = xScale.invert(px);
          const y = yScale.invert(py);

          const dx = r * x * (1 - x / K) - x * y;
          const dy = y * (x - gamma);

          const denom = Math.sqrt(dx*dx + dy*dy);
          if (denom > 1e-6) {
              const nx = dx / denom;
              const ny = dy / denom;
              // Length proportional to strength (magnitude), capped to prevent massive overlap
              const len = options?.propLength ? Math.min(denom * 20, 35) : 16;
              
              // Arrow line
              g.append("line")
                  .attr("x1", px)
                  .attr("y1", py)
                  .attr("x2", px + nx * len)
                  .attr("y2", py - ny * len)
                  .attr("stroke", "#ccc")
                  .attr("stroke-width", 1.2)
                  .attr("marker-end", `url(#${markerId})`);
          }
      }
  }

  // Axes
  g.append("g").attr("transform", `translate(0,${innerH})`).call(d3.axisBottom(xScale));
  g.append("g").call(d3.axisLeft(yScale));

  // x-nullcline
  const lineGen = d3.line()
      .x(d => xScale(d[0]))
      .y(d => yScale(d[1]));

  g.append("path")
      .datum(xNullData)
      .attr("fill", "none")
      .attr("stroke", "blue")
      .attr("stroke-width", 2.5)
      .attr("d", lineGen);
      
  g.append("line") // Total extinction line for X
      .attr("x1", xScale(0))
      .attr("x2", xScale(0))
      .attr("y1", yScale(yMin))
      .attr("y2", yScale(yMax))
      .attr("stroke", "blue")
      .attr("stroke-width", 2.5);

  // y-nullcline
  g.append("line") // y=0 extinction line
      .attr("x1", xScale(xMin))
      .attr("x2", xScale(xMax))
      .attr("y1", yScale(0))
      .attr("y2", yScale(0))
      .attr("stroke", "red")
      .attr("stroke-width", 2.5);

    g.append("line")
      .attr("x1", xScale(gamma))
      .attr("x2", xScale(gamma))
      .attr("y1", yScale(yMin))
      .attr("y2", yScale(yMax))
      .attr("stroke", "red")
      .attr("stroke-width", 2.5);

  if (options?.showSim) {
      const dt = 0.05;
      const Ttotal = 100; // Run for enough time to see trajectory
      
      let x = options.x0;
      let y = options.y0;
      const pathData = [[x, y]];
      
      for (let t = 0; t < Ttotal; t += dt) {
          const fx = (X, Y) => r * X * (1 - X / K) - X * Y;
          const fy = (X, Y) => Y * (X - gamma);

          const stage1x = fx(x, y);
          const stage1y = fy(x, y);
          
          const stage2x = fx(x + 0.5 * dt * stage1x, y + 0.5 * dt * stage1y);
          const stage2y = fy(x + 0.5 * dt * stage1x, y + 0.5 * dt * stage1y);
          
          const stage3x = fx(x + 0.5 * dt * stage2x, y + 0.5 * dt * stage2y);
          const stage3y = fy(x + 0.5 * dt * stage2x, y + 0.5 * dt * stage2y);
          
          const stage4x = fx(x + dt * stage3x, y + dt * stage3y);
          const stage4y = fy(x + dt * stage3x, y + dt * stage3y);

          x += (dt / 6) * (stage1x + 2 * stage2x + 2 * stage3x + stage4x);
          y += (dt / 6) * (stage1y + 2 * stage2y + 2 * stage3y + stage4y);
          
          // Stop drawing if it leaves arbitrary bounds completely
          if (x < xMin-1 || x > xMax+1 || y < yMin-1 || y > yMax+10) break;
          
          pathData.push([x, y]);
      }
      
      const trajGen = d3.line()
          .x(d => xScale(d[0]))
          .y(d => yScale(d[1]));
          
      g.append("path")
          .datum(pathData)
          .attr("fill", "none")
          .attr("stroke", "purple")
          .attr("stroke-width", 2.5)
          .attr("d", trajGen);
          
      g.append("circle")
          .attr("cx", xScale(options.x0))
          .attr("cy", yScale(options.y0))
          .attr("r", 5)
          .attr("fill", "purple");
  }
}

// ----------------------------
// Update loop
// ----------------------------
function updatePlots() {
  const options = {
    r: parseFloat(document.getElementById("r").value),
    k1: parseFloat(document.getElementById("k1").value),
    k2: parseFloat(document.getElementById("k2").value),
    gamma2: parseFloat(document.getElementById("gamma2").value),
    x0: parseFloat(document.getElementById("x0").value),
    y0: parseFloat(document.getElementById("y0").value),
    propLength: document.getElementById("propLength").checked,
    showSim: document.getElementById("showSim").checked
  };

  drawPhasePlane("#plot1", options.r, options.k1, 1, options);
  drawPhasePlane("#plot2", options.r, options.k2, options.gamma2, options);
}

["r", "k1", "k2", "gamma2", "x0", "y0", "propLength", "showSim"].forEach(id => {
  document.getElementById(id).addEventListener("input", updatePlots);
});

window.addEventListener("load", updatePlots);
