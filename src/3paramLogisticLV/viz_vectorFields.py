import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.gridspec import GridSpec
from pathlib import Path


def env_params():
	# default parameters (can be changed when calling)
	params = {
		'A': {'r': 1.0, 'K': 1.0, 'gamma': 1.0},
		'B': {'r': 4.0, 'K': 2.0, 'gamma': 2.2},
	}
	return params


def lotka_prey_growth(x, y, r, K):
	return r * x * (1 - x / K) - x * y


def lotka_predator_growth(x, y, gamma):
	return y * (x - gamma)


def vector_field(x, y, env):
	r = env['r']
	K = env['K']
	gamma = env['gamma']
	dx = lotka_prey_growth(x, y, r, K)
	dy = lotka_predator_growth(x, y, gamma)
	return dx, dy


def averaged_params(envA, envB):
	return {
		'r': 0.5 * (envA['r'] + envB['r']),
		'gamma': 0.5 * (envA['gamma'] + envB['gamma'])
	}


def averaged_fixed_points(envA, envB):
	# The averaged system keeps both carrying capacities in the prey-nullcline.
	rA = envA['r']
	KA = envA['K']
	gammaA = envA['gamma']
	rB = envB['r']
	KB = envB['K']
	gammaB = envB['gamma']

	fixed = [(0.0, 0.0)]

	# Prey-only equilibrium: solve the averaged prey equation at y = 0.
	prey_denominator = (rA / KA) + (rB / KB)
	if prey_denominator > 0:
		xstar = (rA + rB) / prey_denominator
		fixed.append((xstar, 0.0))

	# Coexistence equilibrium: x is pinned by the averaged predator equation.
	xstar = 0.5 * (gammaA + gammaB)
	ystar = 0.5 * (
		rA * (1 - xstar / KA) +
		rB * (1 - xstar / KB)
	)
	if ystar > 0:
		fixed.append((xstar, ystar))

	return fixed



def jacobian(env, x, y):
	r = env['r']
	K = env['K']
	gamma = env['gamma']
	return np.array([
		[r * (1 - 2 * x / K) - y, -x],
		[y, x - gamma],
	])


def averaged_jacobian(envA, envB, x, y):
	return 0.5 * (jacobian(envA, x, y) + jacobian(envB, x, y))


def is_stable_equilibrium(jac, tol=1e-8):
    eigenvalues = np.linalg.eigvals(jac)
    return np.all(np.real(eigenvalues) <= tol)


def averaged_vector_field(x, y, envA, envB):
	# average the vector fields (arithmetic average)
	dx1, dy1 = vector_field(x, y, envA)
	dx2, dy2 = vector_field(x, y, envB)
	return 0.5 * (dx1 + dx2), 0.5 * (dy1 + dy2)


def find_fixed_points(env, x_range=(0.01, 3.0), tol=1e-8):
	# For Lotka-type system, trivial fixed points are (0,0) and (K,0), and coexistence when x=gamma
	r = env['r']
	K = env['K']
	gamma = env['gamma']
	fixed = []
	fixed.append((0.0, 0.0))
	fixed.append((K, 0.0))
	# interior fixed point candidate: x = gamma, solve for y from dx/dt = 0
	xstar = gamma
	if xstar > 0 and K is not None:
		# compute y* from prey equation: 0 = r x*(1-x*/K) - x* y* => y* = r (1 - x*/K)
		ystar = r * (1 - xstar / K)
		if ystar > 0:
			fixed.append((xstar, ystar))
	return fixed


def plot_vector_fields(params=None, outdir=None):
	if params is None:
		params = env_params()

	# increase default font sizes for readability in SVG
	plt.rcParams.update({
		'font.size': 18,
		'axes.titlesize': 22,
		'axes.labelsize': 20,
		'xtick.labelsize': 18,
		'ytick.labelsize': 18,
		'legend.fontsize': 18,
	})

	envA = params['A']
	envB = params['B']
	envAvg = averaged_params(envA, envB)

	# grid for plotting
	x = np.linspace(0.01, max(envA['K'], envB['K']) * 1.05, 30)
	y = np.linspace(0.0, max(envA['K'], envB['K']) * 0.4, 30)
	X, Y = np.meshgrid(x, y)

	# compute vector fields
	DX_A, DY_A = vector_field(X, Y, envA)
	DX_B, DY_B = vector_field(X, Y, envB)
	DX_avg, DY_avg = averaged_vector_field(X, Y, envA, envB)

	# magnitudes for unified colorbar
	M_A = np.hypot(DX_A, DY_A)
	M_B = np.hypot(DX_B, DY_B)
	M_avg = np.hypot(DX_avg, DY_avg)
	global_vmin = float(np.nanmin([np.nanmin(M_A), np.nanmin(M_B), np.nanmin(M_avg)]))
	global_vmax = float(np.nanmax([np.nanmax(M_A), np.nanmax(M_B), np.nanmax(M_avg)]))
	global_norm = plt.Normalize(vmin=global_vmin, vmax=global_vmax)

	# plotting style using GridSpec: three equal plot columns + narrow colorbar column
	fig = plt.figure(figsize=(20, 6))
	fig.patch.set_facecolor('none')
	gs = GridSpec(nrows=1, ncols=4, width_ratios=[1, 1, 1, 0.08], wspace=0.25)
	ax0 = fig.add_subplot(gs[0, 0])
	ax1 = fig.add_subplot(gs[0, 1])
	ax2 = fig.add_subplot(gs[0, 2])
	cax = fig.add_subplot(gs[0, 3])
	axes = [ax0, ax1, ax2]

	fields = [ (DX_A, DY_A, envA, 'Environment A'),
			   (DX_B, DY_B, envB, 'Environment B'),
			   (DX_avg, DY_avg, envAvg, 'Averaged Field') ]

	# compute small negative padding so axes extend slightly left/down of 0
	x_lower = min(0.0, x.min()) - 0.05 * (x.max() - x.min())
	x_upper = x.max()
	y_lower = min(0.0, y.min()) - 0.05 * (y.max() - y.min())
	y_upper = y.max()

	mags = [M_A, M_B, M_avg]
	for ax, (DX, DY, env, title), M in zip(axes, fields, mags):
		# normalize for quiver
		DXn = DX / (M + 1e-8)
		DYn = DY / (M + 1e-8)
		cmap = cm.inferno
		q = ax.quiver(X, Y, DXn, DYn, M, cmap=cmap, norm=global_norm, pivot='mid', scale=30, alpha=0.95)
		ax.set_title(title)
		ax.set_xlabel('Prey x')
		ax.set_ylabel('Predator y')
		ax.set_xlim(x_lower, x_upper)
		ax.set_ylim(y_lower, y_upper)
		ax.grid(True, linestyle='--', alpha=0.4)
		ax.set_axisbelow(True)
		ax.patch.set_facecolor('none')

		# add colorbar for this axis
		# (single shared colorbar will be added after the loop)

		# highlight fixed points for the appropriate field
		if title == 'Averaged Field':
			fixed = averaged_fixed_points(envA, envB)
		else:
			fixed = find_fixed_points(env)

		for (fx, fy) in fixed:
			if fx is None or fy is None:
				continue
			if title == 'Averaged Field':
				stable = is_stable_equilibrium(averaged_jacobian(envA, envB, fx, fy))
			else:
				stable = is_stable_equilibrium(jacobian(env, fx, fy))
			if not stable:
				continue
			ax.plot(fx, fy, 'o', color='red', ms=11, mec='k', mew=1.2)
			ax.text(fx, fy, f'  ({fx:.2f},{fy:.2f})', color='red', fontsize=16, va='bottom')

	# add a single shared colorbar in the dedicated colorbar axis
	cmap = cm.inferno
	sm = cm.ScalarMappable(cmap=cmap, norm=global_norm)
	sm.set_array([])
	fig.colorbar(sm, cax=cax)
	cax.tick_params(labelsize=16)
	cax.yaxis.labelpad = 6
	cax.set_ylabel('vector magnitude', fontsize=16)

	# adjust spacing explicitly (avoids tight_layout warning with GridSpec + colorbar axis)
	fig.subplots_adjust(left=0.06, right=0.98, top=0.95, bottom=0.10, wspace=0.25)

	# save
	script_dir = Path(__file__).parent
	outdir = Path(outdir) if outdir is not None else script_dir / 'plots/vector_fields'
	outdir.mkdir(parents=True, exist_ok=True)
	outpath = outdir / 'vector_fields.svg'
	plt.savefig(outpath, dpi=300, format='svg', bbox_inches='tight', facecolor='none', transparent=True)
	plt.savefig(outpath.with_suffix('.pdf'), dpi=300, format='pdf', bbox_inches='tight', facecolor='none', transparent=True)
	plt.close(fig)
	print(f"Saved vector-field figure to: {outpath}")


if __name__ == '__main__':
	plot_vector_fields()
