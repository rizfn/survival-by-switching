const canvas = document.getElementById('canvas');
const ctx = canvas.getContext('2d');
const envText = document.getElementById('env-text');
const timeControls = document.getElementById('time-controls');
const btnTMinus = document.getElementById('btn-t-minus');
const btnTPlus = document.getElementById('btn-t-plus');

let width, height;
let currentEnv = 0;
let activeEnv = 0;
let lastSwitchTime = Date.now();
let switchIntervalSec = 1.0; // Seconds between oscillation
let isAverageMode = false;
const NUM_PARTICLES = 4000;
let particles = [];

// x from 0.0 to 2.5, y from -0.2 to 2.0, but negative y will not be rendered.
const xMin = 0.0, xMax = 3.0;
const yMin = -0.1, yMax = 1.0;

function resize() {
  width = canvas.width = window.innerWidth;
  height = canvas.height = window.innerHeight;
  initParticles();
  positionHtmlLabels();
}
window.addEventListener('resize', resize);

function toScreen(x, y) {
  return {
    px: ((x - xMin) / (xMax - xMin)) * width,
    py: height - ((y - yMin) / (yMax - yMin)) * height
  };
}

class Particle {
  constructor() {
    this.reset(true);
  }
  
  reset(randomizeAge = false) {
    this.x = xMin + Math.random() * (xMax - xMin);
    this.y = Math.random() * yMax;

    this.age = randomizeAge ? Math.random() * 200 : 0;
    this.maxAge = 150 + Math.random() * 200; // Live enough to trace full paths
    this.dx = 0;
    this.dy = 0;
  }
  
  update(dt) {
    let dx = 0, dy = 0;

    if (isAverageMode) {
      const dx1 = 1.0 * this.x * (1.0 - this.x / 1.0) - this.x * this.y;
      const dy1 = this.y * (this.x - 1.0);
      const dx2 = 4.0 * this.x * (1.0 - this.x / 2.0) - this.x * this.y;
      const dy2 = this.y * (this.x - 2.2);
      dx = (dx1 + dx2) / 2.0;
      dy = (dy1 + dy2) / 2.0;
    } else if (activeEnv === 0) {
      // Env 1: K=1, r=1, gamma=1
      dx = 1.0 * this.x * (1.0 - this.x / 1.0) - this.x * this.y;
      dy = this.y * (this.x - 1.0);
      
    } else if (activeEnv === 1) {
      // Env 2: K=2.0, r=4.0, gamma=2.2
      dx = 4.0 * this.x * (1.0 - this.x / 2.0) - this.x * this.y;
      dy = this.y * (this.x - 2.2);
      
    }

    // Use natural velocities without normalization
    this.dx = dx;
    this.dy = dy;

    // Scale real-world dt (in seconds) by ~0.5 to keep speeds reasonable
    const simDt = dt * 0.5;
    this.x += dx * simDt;
    this.y += dy * simDt;
    
    // Original age was ++ per frame (~60/sec)
    this.age += dt * 60;

    // Reset if off-screen, out of life, or dipped below y=0
    if (this.x < xMin || this.x > xMax || this.y < 0 || this.y > yMax || this.age > this.maxAge) {
      this.reset();
    }
  }
  
  draw() {
    if (this.y < 0) return; // Do not render if inadvertently negative (should not happen mathematically if y>0 initially)

    const pos = toScreen(this.x, this.y);
    
    // Color based on environment
    if (isAverageMode) {
        ctx.fillStyle = 'rgba(255, 255, 255, 0.8)'; // White
    } else if (activeEnv === 0) {
        ctx.fillStyle = 'rgba(59, 130, 246, 0.8)'; // Blue
    } else if (activeEnv === 1) {
        ctx.fillStyle = 'rgba(239, 68, 68, 0.8)'; // Red
    }
    
    ctx.beginPath();
    ctx.arc(pos.px, pos.py, 1.4, 0, Math.PI * 2);
    ctx.fill();
  }
}

function positionHtmlLabels() {
    const p1 = toScreen(1.0, 0);
    const lbl1 = document.getElementById('label-1');
    if (lbl1) {
        lbl1.style.left = p1.px + 'px';
        lbl1.style.top = (p1.py + 15) + 'px';
        lbl1.innerText = '(1, 0)';
    }

    const p2 = toScreen(2.0, 0);
    const lbl2 = document.getElementById('label-2');
    if (lbl2) {
        lbl2.style.left = p2.px + 'px';
        lbl2.style.top = (p2.py + 15) + 'px';
        lbl2.innerText = '(2, 0)';
    }
}

function initParticles() {
  particles = [];
  for (let i = 0; i < NUM_PARTICLES; i++) {
    particles.push(new Particle());
  }
}

function drawFixedPoint(x, y, color, isSink) {
    const p = toScreen(x, y);
    
    if (isSink) {
        const time = Date.now() / 1000;
        const throb = 1 + 0.15 * Math.sin(time * 5); // throbs smoothly with reduced amplitude
        
        ctx.beginPath();
        ctx.arc(p.px, p.py, 12 * throb, 0, Math.PI * 2);
        ctx.fillStyle = color;
        ctx.globalAlpha = 0.5;
        ctx.fill();
        
        ctx.globalAlpha = 1.0;
        ctx.beginPath();
        ctx.arc(p.px, p.py, 6, 0, Math.PI * 2);
        ctx.fillStyle = color;
        ctx.fill();
    } else {
        // Source is now minimal, no glow
        ctx.globalAlpha = 1.0;
        ctx.beginPath();
        ctx.arc(p.px, p.py, 3, 0, Math.PI * 2);
        ctx.fillStyle = '#aaa';
        ctx.fill();
    }
}

function drawGrid() {
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.05)';
    ctx.lineWidth = 1;
    ctx.beginPath();
    for(let x = xMin; x <= xMax; x += 0.25) {
        const top = toScreen(x, yMax);
        const bot = toScreen(x, 0); // Don't vertically draw grid beyond y=0
        ctx.moveTo(top.px, top.py);
        ctx.lineTo(bot.px, bot.py);
    }
    for(let y = 0.0; y <= yMax; y += 0.25) {
        const left = toScreen(xMin, y);
        const right = toScreen(xMax, y);
        ctx.moveTo(left.px, left.py);
        ctx.lineTo(right.px, right.py);
    }
    ctx.stroke();

    // Emphasize y=0 (X-axis base)
    const left0 = toScreen(xMin, 0);
    const right0 = toScreen(xMax, 0);
    ctx.strokeStyle = '#fff';
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(left0.px, left0.py);
    ctx.lineTo(right0.px, right0.py);
    ctx.stroke();

    // Emphasize x=0 (Y-axis base)
    const top0 = toScreen(0, yMax);
    const bot0 = toScreen(0, 0);
    ctx.beginPath();
    ctx.moveTo(top0.px, top0.py);
    ctx.lineTo(bot0.px, bot0.py);
    ctx.stroke();
}

let lastFrameTime = performance.now();

function animate() {
  const now = performance.now();
  let dt = (now - lastFrameTime) / 1000.0;
  if (dt > 0.1) dt = 0.1; // Cap to avoid huge jumps across idle tabs
  lastFrameTime = now;

  if (currentEnv === 2) {
    if (!isAverageMode && Date.now() - lastSwitchTime >= switchIntervalSec * 1000) {
      activeEnv = activeEnv === 0 ? 1 : 0;
      lastSwitchTime = Date.now();
    }
  } else {
    activeEnv = currentEnv;
    lastSwitchTime = Date.now();
  }

  ctx.fillStyle = 'rgba(0, 0, 0, 0.08)';
  ctx.fillRect(0, 0, width, height);

  drawGrid();

  for (let i = 0; i < particles.length; i++) {
    particles[i].update(dt);
    particles[i].draw();
  }

  // Draw points above Trails
  if (isAverageMode) {
      drawFixedPoint(0.0, 0, '#aaa', false); // Source at origin
      drawFixedPoint(1.6, 0.1, '#ffffff', true); // White Sink
  } else if (activeEnv === 0) {
      drawFixedPoint(0.0, 0, '#aaa', false); // Source at origin
      drawFixedPoint(1.0, 0, '#3b82f6', true); // Blue Sink
  } else if (activeEnv === 1) {
      drawFixedPoint(0.0, 0, '#aaa', false); // Source at origin
      drawFixedPoint(2.0, 0, '#ef4444', true); // Red Sink
  }

  positionHtmlLabels(); // Keep them locked nicely

  requestAnimationFrame(animate);
}

function updateEnvText() {
    isAverageMode = (currentEnv === 2 && switchIntervalSec === 0);
    if (currentEnv === 2) {
        if (switchIntervalSec === 0) {
            envText.innerText = `Infinite Switching (Averaged Env)`;
        } else {
            envText.innerText = `Switching Environment every ${switchIntervalSec.toFixed(1)}s`;
        }
        if(timeControls) timeControls.style.display = 'block';
    } else if (currentEnv === 0) {
        envText.innerText = `Environment 1 (K=1, r=1, γ=1)`;
        if(timeControls) timeControls.style.display = 'none';
    } else if (currentEnv === 1) {
        envText.innerText = `Environment 2 (K=2.0, r=4.0, γ=2.2)`;
        if(timeControls) timeControls.style.display = 'none';
    }
}

function toggleEnvironment() {
    currentEnv = (currentEnv + 1) % 3;
    if (currentEnv === 2) {
        activeEnv = 0;
        lastSwitchTime = Date.now();
    }
    updateEnvText();
    
    // Clear fully on transition
    ctx.fillStyle = '#000';
    ctx.fillRect(0, 0, width, height);
}

document.addEventListener('keydown', (e) => {
  if (e.code === 'Space') {
    e.preventDefault();
    toggleEnvironment();
  } else if (e.code === 'ArrowUp' || e.key === 'w' || e.key === 'W') {
    e.preventDefault();
    if (currentEnv === 2) {
        switchIntervalSec = Math.min(10.0, switchIntervalSec + 0.1);
        updateEnvText();
    }
  } else if (e.code === 'ArrowDown' || e.key === 's' || e.key === 'S') {
    e.preventDefault();
    if (currentEnv === 2) {
        switchIntervalSec = Math.max(0.0, switchIntervalSec - 0.1);
        if (switchIntervalSec < 0.05) switchIntervalSec = 0;
        updateEnvText();
    }
  }
});

const toggleBtn = document.getElementById('toggle-btn');
if (toggleBtn) {
    toggleBtn.addEventListener('click', (e) => {
        e.preventDefault();
        toggleEnvironment();
        toggleBtn.blur(); // Remove focus to prevent spacebar double-trigger
    });
}
if (btnTMinus) {
    btnTMinus.addEventListener('click', (e) => {
        e.preventDefault();
        if (currentEnv === 2) {
            switchIntervalSec = Math.max(0.0, switchIntervalSec - 0.1);
            if (switchIntervalSec < 0.05) switchIntervalSec = 0;
            updateEnvText();
        }
        btnTMinus.blur();
    });
}
if (btnTPlus) {
    btnTPlus.addEventListener('click', (e) => {
        e.preventDefault();
        if (currentEnv === 2) {
            switchIntervalSec = Math.min(20.0, switchIntervalSec + 0.1);
            updateEnvText();
        }
        btnTPlus.blur();
    });
}

// Start
resize();
animate();
