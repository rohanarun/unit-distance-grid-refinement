"""Render progress/progress.png from progress/progress.json.

Usage: python -B progress/make_chart.py   (needs matplotlib)
"""
from datetime import date, datetime, timedelta
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE / 'progress.json').read_text())

SURFACE = '#fcfcfb'
INK = '#0b0b0b'
INK2 = '#52514e'
MUTED = '#898781'
GRID = '#e1e0d9'
AXIS = '#c3c2b7'
SERIES = ['#2a78d6', '#eb6834', '#1baf7a']  # validated categorical slots 1-3

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica Neue', 'Helvetica', 'Arial', 'DejaVu Sans'],
    'axes.edgecolor': AXIS, 'axes.labelcolor': INK2,
    'xtick.color': MUTED, 'ytick.color': MUTED,
    'text.color': INK,
})


def at(iso, hours=0):
    return datetime.combine(date.fromisoformat(iso), datetime.min.time()) + timedelta(hours=12 + hours)


fig, ax = plt.subplots(figsize=(11, 6.2), dpi=200)
fig.patch.set_facecolor(SURFACE)
ax.set_facecolor(SURFACE)

base = DATA['baseline']
series = DATA['series']
first_change = min(p['date'] for s in series for p in s['points'][1:])
last_change = max(p['date'] for s in series for p in s['points'])
end = at(last_change, 48)

# Shared history: every branch inherits the source's grid until the first refinement.
ax.plot([at(base['date']), at(first_change)], [base['edges']] * 2, color=INK2, lw=2,
        solid_capstyle='round', zorder=2)

# Overlapping vertical drops are offset by a few hours so each branch stays visible.
jitter = {0: -5, 1: 0, 2: 5}
for k, s in enumerate(series):
    pts = s['points']
    xs = [at(first_change, jitter[k])]
    ys = [base['edges']]
    for p in pts[1:]:
        xs.append(at(p['date'], jitter[k]))
        ys.append(p['edges'])
    xs.append(end)
    ys.append(pts[-1]['edges'])
    ax.step(xs, ys, where='post', color=SERIES[k], lw=2, solid_capstyle='round', zorder=3)
    for p in pts[1:]:
        ax.plot(at(p['date'], jitter[k]), p['edges'], 'o', ms=8, color=SERIES[k],
                markeredgecolor=SURFACE, markeredgewidth=2, zorder=4)
    final = pts[-1]
    ax.annotate(f"{s['name']}: {final['grid']} = {final['edges']}", xy=(end, final['edges']),
                xytext=(6, 0), textcoords='offset points', va='center', ha='left',
                fontsize=10, color=INK)

# Intermediate labels, de-duplicated when two branches share a value.
seen = set()
for k, s in enumerate(series):
    for p in s['points'][1:-1]:
        key = (p['date'], p['edges'])
        if key in seen:
            continue
        seen.add(key)
        ax.annotate(f"{p['grid']} = {p['edges']}", xy=(at(p['date']), p['edges']),
                    xytext=(-9, 0), textcoords='offset points', ha='right', va='center',
                    fontsize=9, color=INK2)

ax.plot(at(base['date']), base['edges'], 'o', ms=9, color=INK, markeredgecolor=SURFACE,
        markeredgewidth=2, zorder=5)
ax.annotate(f"{base['label']} preprint, {date.fromisoformat(base['date']):%b %d}: "
            f"{base['grid']} grid = {base['edges']} edges",
            xy=(at(base['date']), base['edges']), xytext=(10, -3), textcoords='offset points',
            va='top', fontsize=10, color=INK, fontweight='bold')

ax.set_ylim(0, 110)
ax.set_xlim(at(base['date'], -24), end + timedelta(days=10))
ax.set_ylabel('edges in the final algebraic witness', color=INK2)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=3))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
ax.yaxis.grid(True, color=GRID, lw=0.8)
ax.set_axisbelow(True)
for side in ('top', 'right'):
    ax.spines[side].set_visible(False)
ax.tick_params(length=0)

handles = [Line2D([], [], color=INK2, lw=2, label='shared history (source grid)')]
handles += [Line2D([], [], color=SERIES[k], lw=2, label=s['name']) for k, s in enumerate(series)]
ax.legend(handles=handles, loc='upper right', frameon=False, fontsize=9, labelcolor=INK2)

fig.suptitle('Unit distances: shrinking the final generic grid in the OpenAI argument',
             x=0.06, ha='left', fontsize=14, fontweight='bold', color=INK)
ax.set_title('Left points x right anchors needed for the final sign-change contradiction. '
             'Conditional on the same source lemmas; no new exponent.',
             loc='left', fontsize=10, color=INK2, pad=8)
fig.text(0.06, 0.015, 'github.com/rohanarun/unit-distance-grid-refinement  ·  '
         'baseline: OpenAI, "A power saving for planar unit distances" (Sep 23, 2026)',
         fontsize=8.5, color=MUTED)
fig.subplots_adjust(left=0.08, right=0.77, top=0.86, bottom=0.12)
out = HERE / 'progress.png'
fig.savefig(out, facecolor=SURFACE)
print(f'wrote {out}')
