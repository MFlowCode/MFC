// Shared by the dashboard pages (index.html, benchmarks.html, ci.html, tests.html): where the data is, how it
// is read, and the few things more than one page works out from it.
"use strict";
const MFCCI = (() => {
    const BASE = "https://raw.githubusercontent.com/MFlowCode/MFC-ci-data/main/";
    const REPO = "https://github.com/MFlowCode/MFC";
    const lines = text => text.split("\n").filter(Boolean).map(l => JSON.parse(l));
    const get = path => fetch(BASE + path).then(r => r.ok ? r.text() : "").catch(() => "");
    const med = a => a.length ? [...a].sort((x, y) => x - y)[a.length >> 1] : null;
    const esc = t => String(t ?? "").replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]);
    const day = t => new Date(t).toLocaleDateString(undefined, { month: "short", day: "numeric" });
    const PR_OF = r => r.pr ?? (/\(#(\d+)\)\s*$/.exec(r.subject || r.title || "") || [])[1];
    const commitLink = (sha, pr) => pr ? `<a href="${REPO}/pull/${pr}" target="_blank" rel="noopener">#${pr}</a>`
                                       : `<a href="${REPO}/commit/${sha}" target="_blank" rel="noopener"><code>${sha.slice(0, 7)}</code></a>`;

    // ---- lanes: a machine, a compiler and an offload model (and on Phoenix, for benchmarks, a GPU model) ----
    // The lanes are whatever the records name (ci-record.py: "<cluster>-<device>-<offload>"). The clusters
    // known here get their names and their order; any other shows as its record spells it.
    const CLUSTERS = { frontier: ["Frontier", "CCE"], frontier_amd: ["Frontier", "AMD"], phoenix: ["Phoenix", "NVIDIA"] };
    const OFFLOAD = { acc: "OpenACC", omp: "OpenMP", none: "CPU" }, MODULES = { frontier: "f", frontier_amd: "famd", phoenix: "p" };
    function describe(key) {
        const [cluster, device, offload] = key.split("-"), words = cluster.split("_").map(w => w[0].toUpperCase() + w.slice(1));
        const [machine, compiler] = CLUSTERS[cluster] || [words[0], (words[1] || "").toUpperCase()];
        const how = OFFLOAD[offload] || offload;
        return { key, cluster, device, machine, compiler, offload: how, label: [machine === "Phoenix" && compiler === "NVIDIA" ? "" : machine, compiler, how].filter(Boolean).join(" · ") };
    }
    const rank = l => [Object.keys(CLUSTERS).indexOf(l.cluster) + 1 || 99, l.device === "gpu" ? 0 : 1, l.cluster, l.offload];
    const CI_LANES = [];                                         // filled, in place, as records load
    function seeLanes(keys) {
        for (const key of new Set(keys)) if (!CI_LANES.some(l => l.key === key)) CI_LANES.push(describe(key));
        CI_LANES.sort((a, b) => { const x = rank(a), y = rank(b); for (let i = 0; i < x.length; i++) if (x[i] !== y[i]) return x[i] < y[i] ? -1 : 1; return 0; });
    }
    seeLanes(["frontier-gpu-acc", "frontier-gpu-omp", "frontier_amd-gpu-omp", "phoenix-gpu-acc", "phoenix-gpu-omp"]);
    const LANE_KEY = `<details class="key"><summary>Lane names</summary><p><strong>Frontier</strong>: OLCF, AMD GPUs. <strong>Phoenix</strong>: Georgia Tech, NVIDIA GPUs.
        <strong>CCE</strong>: Cray compiler. <strong>AMD</strong>: AMD flang. <strong>NVIDIA</strong>: NVHPC compiler, on Phoenix. <strong>OpenACC</strong>, <strong>OpenMP</strong>: GPU offload models. <strong>CPU</strong>: no offload.</p></details>`;
    const FRONTIER = { "frontier-gpu-acc": "CCE · OpenACC", "frontier-gpu-omp": "CCE · OpenMP", "frontier_amd-gpu-omp": "AMD · OpenMP" };
    const GPUS = ["h200", "h100", "a100", "l40s", "v100"];
    const laneKey = r => `${r.cluster}-${r.device}-${r.interface}${r.hardware ? `-${r.hardware}` : ""}`;
    function laneOf(r) {
        const key = laneKey(r);
        if (r.cluster.startsWith("frontier")) return { key, view: "Frontier", label: FRONTIER[key] || key, order: Object.keys(FRONTIER).indexOf(key) };
        const gpu = r.hardware || "unknown GPU", at = GPUS.indexOf(gpu);
        return { key, view: `NVIDIA · ${r.interface === "acc" ? "OpenACC" : "OpenMP"}`, label: gpu.toUpperCase().replace("UNKNOWN GPU", "GPU not recorded"), order: at < 0 ? GPUS.length : at };
    }
    const benchLanes = records => { const seen = new Map(); for (const r of records) if (!seen.has(r.lane)) seen.set(r.lane, laneOf(r)); return [...seen.values()]; };

    // ---- loading ----
    const cache = {};
    const months = (n = 4) => [...Array(n).keys()].map(i => { const d = new Date(); d.setUTCDate(1); d.setUTCMonth(d.getUTCMonth() - i); return `ci/${d.toISOString().slice(0, 7)}.jsonl`; });
    const load = {
        bench: () => cache.bench ??= get("benchmarks.jsonl").then(t => lines(t).map(r => ({ ...r, lane: laneKey(r), t: Date.parse(r.date) })).sort((a, b) => a.t - b.t)),
        runs:  () => cache.runs ??= Promise.all(months().map(get)).then(ms => { const rs = ms.flatMap(lines).map(r => ({ ...r, t: Date.parse(r.date) })).sort((a, b) => a.t - b.t);
            seeLanes(rs.slice(-400).flatMap(r => r.jobs.map(j => j.l))); return rs; }),
        // What each lane loads today: its lines of toolchain/modules on master.
        modules: () => cache.modules ??= fetch("https://raw.githubusercontent.com/MFlowCode/MFC/master/toolchain/modules").then(r => r.ok ? r.text() : "").catch(() => "").then(text => l => {
            const slug = MODULES[l.cluster], found = [];
            if (slug) for (const line of text.split("\n")) { const m = line.match(new RegExp(`^${slug}-(all|${l.device})\\s+(.*)`)); if (m) found.push(...m[2].trim().split(/\s+/).filter(w => !w.includes("="))); }
            return found; }),
        usage: () => cache.usage ??= get("usage.jsonl").then(lines),
        builds: () => cache.builds ??= get("builds.jsonl").then(lines),
        scaling: () => cache.scaling ??= get("scaling.jsonl").then(lines, () => []).then(rs => rs.map(r => ({ ...r, t: Date.parse(r.date) })).sort((a, b) => a.t - b.t)),
    };

    // ---- what is running now: asked of GitHub when the page loads (the history holds finished runs only) ----
    // Three unauthenticated requests, kept for two minutes; GitHub allows a visitor 60 an hour. If it
    // cannot be reached, a snapshot beside the data (live.json) is used when there is one.
    const API = "https://api.github.com/repos/MFlowCode/MFC";
    const JOB = /(Frontier|Phoenix) \((CCE|AMD|NVHPC)\) \((gpu|cpu)-(acc|omp|none)(?: \[(\d+\/\d+)\])?\)/;
    const CLUSTER = { "Frontier CCE": "frontier", "Frontier AMD": "frontier_amd", "Phoenix NVHPC": "phoenix" };
    const liveRun = r => ({ id: r.id, workflow: (r.path || "").endsWith("bench.yml") ? "bench" : "test", name: r.name, event: r.event, branch: r.head_branch, sha: r.head_sha,
                            title: r.display_title, status: r.status, t: Date.parse(r.run_started_at || r.created_at) });
    const liveJob = j => { const m = JOB.exec(j.name || ""); if (!m) return null;
        return { l: `${CLUSTER[`${m[1]} ${m[2]}`]}-${m[3]}-${m[4]}`, k: j.name.startsWith("Bench") ? "bench" : j.name.startsWith("Case Opt") ? "case-opt" : "test", s: m[5] || null,
                 status: j.status, c: j.conclusion, created: Date.parse(j.created_at), started: j.started_at ? Date.parse(j.started_at) : null, id: j.id }; };
    const isMasterTest = r => r.workflow === "test" && r.branch === "master" && ["push", "schedule", "workflow_dispatch"].includes(r.event);
    function live() {
        if (cache.live) return cache.live;
        try { const kept = JSON.parse(sessionStorage.getItem("mfcci-live") || "null"); if (kept && Date.now() - kept.at < 120e3) return cache.live = Promise.resolve(kept); } catch (e) { /* no storage */ }
        const ask = url => fetch(url, { headers: { Accept: "application/vnd.github+json" } }).then(r => { if (!r.ok) throw new Error(r.status); return r.json(); });
        return cache.live = (async () => {
            try {
                const runs = [...(await ask(`${API}/actions/runs?status=in_progress&per_page=50`)).workflow_runs, ...(await ask(`${API}/actions/runs?status=queued&per_page=50`)).workflow_runs].map(liveRun);
                const m = runs.filter(isMasterTest).sort((a, b) => b.t - a.t)[0];
                const jobs = m ? (await ask(`${API}/actions/runs/${m.id}/jobs?per_page=100`)).jobs.map(liveJob).filter(Boolean) : [];
                const out = { at: Date.now(), runs, master: m ? { run: m, jobs } : null };
                try { sessionStorage.setItem("mfcci-live", JSON.stringify(out)); } catch (e) { /* no storage */ }
                return out;
            } catch (e) {
                const snap = await get("live.json");
                try { return snap ? JSON.parse(snap) : null; } catch (e2) { return null; }
            }
        })();
    }
    // A lane's state in the master run now in progress: null when there is none or the lane is not in it.
    function liveLane(now, lane, kind) {
        const js = now?.master ? now.master.jobs.filter(j => j.l === lane && j.k === kind) : []; if (!js.length) return null;
        const state = js.some(j => j.status === "in_progress") ? "running" : js.some(j => j.status !== "completed") ? "queued"
                    : js.some(j => j.c === "failure") ? "failed" : js.every(j => j.c === "success") ? "passed" : "no result";
        return { state, run: now.master.run, jobs: js };
    }
    const clock = t => new Date(t).toLocaleTimeString(undefined, { hour: "numeric", minute: "2-digit" });
    const ago = t => { const m = Math.max(0, Math.round((Date.now() - t) / 60e3)); return m < 90 ? `${m} min` : `${(m / 60).toFixed(1)} h`; };
    // "12 runs in progress, 3 queued, oldest queued 2.1 h" for a status line; "" when nothing is known.
    function liveSummary(now) {
        if (!now) return null;
        // A run queued for days is waiting on a person (a fork's first run needs approval), not on a runner.
        const fresh = now.runs.filter(r => now.at - r.t < 3 * 864e5), running = fresh.filter(r => r.status === "in_progress"), queued = fresh.filter(r => r.status !== "in_progress");
        return { text: `${running.length} running, ${queued.length} queued` + (queued.length ? ` (oldest ${ago(Math.min(...queued.map(r => r.t)))})` : ""), at: now.at };
    }

    // ---- state that follows you between pages and into a link: ?range=30&view=... ----
    const param = (name, value) => {
        const u = new URL(location.href);
        if (value === undefined) { let v = u.searchParams.get(name); if (v === null) try { v = localStorage.getItem(`mfcci-${name}`); } catch (e) { /* no storage */ } return v; }
        if (value === null || value === "") u.searchParams.delete(name); else u.searchParams.set(name, value);
        try { history.replaceState(null, "", u); } catch (e) { /* sandboxed */ }
        try { if (name === "range") localStorage.setItem("mfcci-range", value); } catch (e) { /* no storage */ }
        return value;
    };
    const withParams = href => { const r = param("range"), [page, hash] = href.split("#"); return (r ? `${page}${page.includes("?") ? "&" : "?"}range=${encodeURIComponent(r)}` : page) + (hash ? `#${hash}` : ""); };

    // ---- the header and footer every page wears ----
    const NAV = [["index.html", "Status"], ["benchmarks.html", "Performance"], ["ci.html", "CI health"], ["tests.html", "Test health"], ["../index.html", "MFC home"], ["../documentation/index.html", "Documentation"]];
    function chrome(active) {
        const head = document.getElementById("site-header"), foot = document.getElementById("site-footer");
        if (head) head.innerHTML = `<div class="wrap"><a class="brand" href="index.html"><img src="../res/banner.png" alt="MFC"><span>Dashboard</span></a>
            <nav>${NAV.map(([href, text]) => `<a href="${/^(index|benchmarks|ci|tests)\.html$/.test(href) ? withParams(href) : href}"${href === active ? ' aria-current="page"' : ""}>${text}</a>`).join("")}</nav></div>`;
        if (foot) foot.innerHTML = `<div class="wrap"><strong>&copy; ${new Date().getFullYear()} Spencer H. Bryngelson and Tim Colonius</strong> &middot; <a href="https://comp-physics.group/">Computational Physics Group @ Georgia Tech</a></div>`;
    }

    // ---- newest record, and a visible warning when it is old ----
    function updated(el, t, staleDays) {
        if (!el || !t) return;
        const days = Math.floor((Date.now() - t) / 864e5);
        el.className = days >= staleDays ? "updated stale" : "updated";
        el.textContent = `Updated ${new Date(t).toLocaleString(undefined, { dateStyle: "medium", timeStyle: "short" })}` + (days >= staleDays ? ` · ${days} days old` : "");
    }

    // ---- sortable tables: click a header ----
    function sortable(table) {
        table.addEventListener("click", e => {
            const th = e.target.closest("th"); if (!th || !table.contains(th) || th.closest("tr") !== table.rows[0]) return;
            const col = [...th.parentNode.children].indexOf(th), dir = th.dataset.dir === "asc" ? -1 : 1;
            const num = s => { const m = /-?\d[\d,]*\.?\d*/.exec(s.replace(/[▲▼✓✕]/g, "")); return m ? parseFloat(m[0].replace(/,/g, "")) * (/^\s*[▼]|−/.test(s) && !/-/.test(m[0]) ? 1 : 1) : NaN; };
            const rows = [...table.rows].slice(1).filter(r => !r.classList.contains("detail") && r.cells.length > col);
            rows.sort((a, b) => { const x = a.cells[col].textContent.trim(), y = b.cells[col].textContent.trim(), nx = num(x), ny = num(y);
                return dir * (isNaN(nx) || isNaN(ny) ? x.localeCompare(y) : nx - ny); });
            for (const o of th.parentNode.children) delete o.dataset.dir; th.dataset.dir = dir === 1 ? "asc" : "desc";
            for (const r of rows) table.tBodies[0].appendChild(r);
        });
    }

    // How a job failed, in a few words, by the cause the recorder gave it.
    const CAUSE = { tests: "tests failed", build: "build failed", cluster: "cluster down", queue: "never left the queue", timeout: "SLURM time limit", node: "bad node",
                    checkout: "git checkout failed", upload: "artifact upload failed", runner: "runner stopped", other: "cause not recognized", "not read": "log not read" };
    // One line per failed job of a cell: the cause, then the log's own line or the failed tests.
    const howFailed = jobs => jobs.filter(j => j.c === "failure").map(j => {
        const cause = CAUSE[j.why || "not read"] || j.why, f = j.f || [];
        return { j, cause, detail: f.length ? f.slice(0, 3).map(t => `${t[1]} (${t[2]}${t[3] ? `: ${t[3]}` : ""})`).join("; ") + (f.length > 3 ? `; and ${f.length - 3} more` : "") : j.e || "" };
    });

    // ---- scaling: the rules of .github/scripts/scaling-check.py ----
    const SCALING = { points: 5, grind: 5, weakFloor: 80, runs: 4 };
    // Per size of the latest run: efficiency against its smallest size, what is usual (the median of
    // the last four runs before it on the same lane), and whether it is outside the tolerances.
    function scalingStatus(records) {
        const latest = records.at(-1); if (!latest) return null;
        const measure = r => ["weak", "strong"].flatMap(kind => { const rs = r.runs.filter(x => x.kind === kind).sort((a, b) => a.nodes - b.nodes), base = rs[0];
            return rs.map(x => ({ ...x, base: x === base, eff: base.grind * base.ranks / (x.grind * x.ranks) * 100 })); });
        const past = records.slice(0, -1).filter(r => r.cluster === latest.cluster && r.interface === latest.interface).slice(-SCALING.runs).map(measure);
        const rows = measure(latest).map(x => {
            const before = past.map(m => m.find(y => y.kind === x.kind && y.nodes === x.nodes)).filter(Boolean), usual = med(before.map(y => x.base ? y.grind : y.eff));
            const bad = x.base ? usual !== null && (x.grind / usual - 1) * 100 > SCALING.grind
                               : usual !== null ? usual - x.eff > SCALING.points : x.kind === "weak" && x.eff < SCALING.weakFloor;
            return { ...x, usual, bad };
        });
        return { latest, rows, bad: rows.filter(r => r.bad), measure };
    }

    // ---- master: each run's result per lane, and each lane's standing now ----
    const WORST = ["success", "cancelled", "timed_out", "failure"];
    // Master Test Suite runs with, per row (lane and kind), the worst result among the row's latest jobs.
    function masterRuns(runs) {
        return runs.filter(r => r.master && r.workflow === "test").map(r => {
            const last = Math.max(...r.jobs.map(j => j.a)), cells = {};
            for (const j of r.jobs.filter(j => j.a === last || !r.jobs.some(o => o.a > j.a && o.l === j.l && o.k === j.k && o.s === j.s))) {
                const key = `${j.l}|${j.k}`, c = cells[key] ||= { c: "success", jobs: [] };
                c.jobs.push(j); if (WORST.indexOf(j.c) > WORST.indexOf(c.c)) c.c = j.c;
            }
            const all = Object.values(cells).map(c => c.c);
            return { ...r, cells, worst: all.every(c => c === "cancelled") ? "cancelled" : all.includes("failure") ? "failure" : all.includes("timed_out") ? "timed_out" : all.includes("cancelled") ? "cancelled" : "success" };
        });
    }
    function laneStatus(runs) {
        const ms = masterRuns(runs), recent = ms.filter(r => r.t >= Date.now() - 7 * 864e5), rows = [];     // a lane is current if master ran it this week
        for (const l of CI_LANES) for (const kind of ["test", "case-opt"]) {
            const k = `${l.key}|${kind}`; if (!recent.some(r => r.cells[k])) continue;
            const had = ms.filter(r => r.cells[k] && r.cells[k].c !== "cancelled"), last = had.at(-1), passed = [...had].reverse().find(r => r.cells[k].c === "success");
            const firstBad = last && last.cells[k].c !== "success" ? had[passed ? had.indexOf(passed) + 1 : 0] : null;
            rows.push({ key: k, lane: l, kind, last, passed, firstBad, ok: !!last && last.cells[k].c === "success", badRuns: firstBad ? had.length - had.indexOf(firstBad) : 0 });
        }
        return { rows, allGreen: [...ms].reverse().find(r => r.worst === "success") };
    }

    // ---- tests: cases that failed or passed only on a retry, per lane ----
    function testStats(runs, days, scope) {
        const from = Date.now() - days * 864e5, js = [], tests = new Map(), byName = new Map();
        for (const r of runs) if (r.t >= from && (scope === "all" || r.master)) for (const j of r.jobs) if (j.k === "test" && j.n !== undefined) js.push({ r, j });
        const at = (uuid, name) => {
            const key = uuid || byName.get(name) || `name:${name}`;
            if (uuid && !name.endsWith("...")) byName.set(name, uuid);
            const t = tests.get(key) || { key, uuid, name, lanes: new Map(), cls: new Map(), where: [] };
            if (name.length > t.name.length || t.name.endsWith("...") && !name.endsWith("...")) t.name = name;
            tests.set(key, t); return t;
        };
        for (const kind of ["f", "r"]) for (const { r, j } of js) for (const [uuid, name, cls] of j[kind] || []) {
            const t = at(uuid, name), c = t.lanes.get(j.l) || { f: 0, r: 0 };
            c[kind]++; t.lanes.set(j.l, c); t.cls.set(cls, (t.cls.get(cls) || 0) + 1); t.where.push({ r, j, kind });
        }
        const rows = [...tests.values()].map(t => ({ ...t, f: [...t.lanes.values()].reduce((a, c) => a + c.f, 0), r: [...t.lanes.values()].reduce((a, c) => a + c.r, 0) }));
        return { jobs: js, tests: rows.sort((a, b) => b.f - a.f || b.r - a.r) };
    }

    // ---- node-hours on Frontier over the last days, from recorded SLURM jobs and, before those, accounting ----
    function frontierHours(runs, usage, days) {
        const from = new Date(Date.now() - days * 864e5).toISOString().slice(0, 10); let live = 0, first = "9999";
        for (const r of runs) for (const j of r.jobs) if (j.l.startsWith("frontier") && (j.u || []).length) { const d = r.date.slice(0, 10); if (d < first) first = d;
            if (d >= from) for (const [, nodes, elapsed] of j.u) live += nodes * elapsed / 3600; }
        return live + usage.filter(u => u.cluster === "frontier" && u.date >= from && u.date < first).reduce((a, u) => a + u.node_hours, 0);
    }

    // A lasting change in level, per lane and case: at each benchmarked commit, the median of the (up to) six
    // commits before against the six from it on, at least four on each side. It counts when it is at least
    // 3% and five times the scatter on either side, and it is placed at the commit where that contrast peaks.
    // Runs across a change to the case or to the modules are never compared; a run over 3x its neighbours
    // and a run whose commit is only estimated are left out. One case moving alone has to move 10%;
    // otherwise two cases on the lane must.
    // The suspects are the commits after the lane's last benchmark at the old level, up to its first at the
    // new one. Lanes that changed over overlapping spans are one event, and its suspects are the commits
    // every one of those lanes allows.
    function shifts(records, lanes) {
        const med = a => { const b = [...a].sort((x, y) => x - y); return b[b.length >> 1]; };
        const scatter = a => { const m = med(a); return 1.4826 * med(a.map(v => Math.abs(v - m))) / m; };
        const W = 6, MIN = 4, spans = [];
        // master's commits in order: each benchmarked commit lists the ones it covers
        const order = [], at = new Map();
        for (const r of records) for (const c of r.covers?.length ? r.covers : [[r.sha, r.subject || "", PR_OF(r)]]) if (!at.has(c[0])) { at.set(c[0], order.length); order.push(c); }
        const slugs = [...new Set(records.flatMap(r => Object.keys(r.cases)))];
        for (const lane of lanes) {
            // A lane with few runs (a Phoenix GPU model CI rarely lands on) cannot carry a judgement.
            if (records.filter(r => r.lane === lane.key && !r.approx).length < 60) continue;
            const hits = [];
            for (const slug of slugs) {
                const runs = records.filter(r => r.lane === lane.key && r.cases[slug] && !r.approx && at.has(r.sha)).sort((x, y) => at.get(x.sha) - at.get(y.sha)), by = new Map();
                runs.forEach((r, i) => {
                    const near = runs.slice(Math.max(0, i - 7), i + 8).map(q => q.cases[slug].grind);
                    if (r.cases[slug].grind > 3 * med(near)) return;                       // a stray run
                    const c = by.get(r.sha) || { r, seg: `${r.cases[slug].def || ""}|${r.env || ""}`, vs: [] };
                    c.vs.push(r.cases[slug].grind); by.set(r.sha, c);
                });
                const cs = [...by.values()].map(c => ({ ...c, v: med(c.vs) })), score = cs.map(() => null);
                for (let i = 1; i < cs.length; i++) {
                    const before = [], after = [];
                    for (let k = i - 1; k >= 0 && before.length < W && cs[k].seg === cs[i].seg; k--) before.push(cs[k].v);
                    for (let k = i; k < cs.length && after.length < W && cs[k].seg === cs[i].seg; k++) after.push(cs[k].v);
                    if (before.length < MIN || after.length < MIN) continue;
                    const pct = med(after) / med(before) - 1, noise = Math.max(scatter(before), scatter(after), 0.002);
                    if (Math.abs(pct) >= 0.03 && Math.abs(pct) >= 5 * noise) score[i] = pct;
                }
                for (let i = 1; i < cs.length; i++) {
                    if (score[i] === null) continue;
                    // the step sits where the contrast is largest among neighbouring commits flagged the same way
                    let peak = true;
                    for (let k = Math.max(1, i - 3); k <= Math.min(cs.length - 1, i + 3); k++)
                        if (k !== i && score[k] !== null && Math.sign(score[k]) === Math.sign(score[i]) && (Math.abs(score[k]) > Math.abs(score[i]) || Math.abs(score[k]) === Math.abs(score[i]) && k < i)) peak = false;
                    if (peak) hits.push({ lane, slug, pct: score[i] * 100, lo: at.get(cs[i - 1].r.sha) + 1, hi: at.get(cs[i].r.sha), r: cs[i].r });
                }
            }
            // per lane: cases that stepped over overlapping spans, in the same direction, are one change
            hits.sort((x, y) => x.lo - y.lo);
            for (const h of hits) {
                const s = spans.find(s => s.lane === lane && Math.sign(s.hits[0].pct) === Math.sign(h.pct) && h.lo <= s.hi && h.hi >= s.lo);
                if (s) { s.lo = Math.max(s.lo, h.lo); s.hi = Math.min(s.hi, h.hi); s.hits.push(h); } else spans.push({ lane, lo: h.lo, hi: h.hi, hits: [h] });
            }
        }
        // Kept: two cases together, or one moving 10% on a lane with a single hardware (a Phoenix lane is
        // one GPU model, with fewer runs, so it needs two); and a span short enough to point at something.
        const solid = spans.filter(s => (s.hits.length >= 2 || Math.abs(s.hits[0].pct) >= 10 && s.lane.view === "Frontier") && s.hi - s.lo < 8).sort((x, y) => x.lo - y.lo), events = [];
        for (const s of solid) {
            const e = events.find(e => Math.sign(e.spans[0].hits[0].pct) === Math.sign(s.hits[0].pct) && s.lo <= e.hi && s.hi >= e.lo);
            if (e) { e.lo = Math.max(e.lo, s.lo); e.hi = Math.min(e.hi, s.hi); e.spans.push(s); } else events.push({ lo: s.lo, hi: s.hi, spans: [s] });
        }
        const first = new Map(); for (const r of records) if (!first.has(r.sha)) first.set(r.sha, r);
        return events.map(e => ({ ...e, suspects: order.slice(e.lo, e.hi + 1), r: first.get(order[e.hi][0]) || e.spans[0].hits[0].r })).sort((a, b) => b.hi - a.hi);
    }
    return { BASE, REPO, CI_LANES, LANE_KEY, lines, get, med, esc, day, PR_OF, commitLink, laneKey, laneOf, benchLanes, load, param, withParams, chrome, SCALING, scalingStatus, updated, sortable, CAUSE, howFailed, live, liveLane, liveSummary, clock, ago,
             masterRuns, laneStatus, testStats, frontierHours, shifts, CAP: 2688 };
})();
