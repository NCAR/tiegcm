#!/usr/bin/env python
"""Job-script building blocks for any machine and scheduler, configured by
config/machines.yaml and config/batch.yaml.

Identical copies live in tiegcm/tiegcmrun and kaiju-private/scripts; keep them in sync.
"""
import getpass, os, re, math, shlex, socket, subprocess, yaml


def account_value(value):
    """Expand $USER / ${USER} in a machines.yaml account_default to the login name."""
    if not isinstance(value, str) or "USER" not in value:
        return value
    return re.sub(r"\$\{USER\}|\$USER\b", getpass.getuser(), value)


# NCAR project-code group names, e.g. nhao0001 or p28100036.
PROJECT_GROUP = re.compile(r"^(?:[a-z]{4}\d{4}|[a-z]\d{8})$", re.IGNORECASE)
_PROJECTS = []


def user_projects():
    """The user's project codes from `id -Gn`, upper-case (cached; [] if none)."""
    if not _PROJECTS:
        try:
            out = subprocess.run(["id", "-Gn"], capture_output=True, text=True, timeout=10).stdout
        except (OSError, subprocess.SubprocessError):
            out = ""
        _PROJECTS.append(list(dict.fromkeys(g.upper() for g in out.split() if PROJECT_GROUP.match(g))))
    return list(_PROJECTS[0])


def dot_relative(path):
    """Prefix a bare relative name with './' so bash does not search CDPATH or PATH for it."""
    path = os.fspath(path)
    if not path:
        raise ValueError("dot_relative: empty path")
    if os.path.isabs(path) or path in (".", "..") or path.startswith(("./", "../")):
        return path
    return "./" + path


def run_relative(target, base_dir, run_root, dot=True):
    """Spell `target` for a reader whose cwd is `base_dir`: relative if inside `run_root`,
    absolute otherwise or when a `..` through a symlinked directory would reach another file."""
    t = os.path.abspath(os.fspath(target))
    b = os.path.abspath(os.fspath(base_dir))
    r = os.path.abspath(os.fspath(run_root))
    if os.path.commonpath([r, t]) != r:
        return t
    rel = os.path.relpath(t, b)
    if os.path.realpath(os.path.join(b, rel)) != os.path.realpath(t):
        return t
    return dot_relative(rel) if dot else rel


def spell(path, base_dir, run_root=None, dot=True):
    """run_relative() that passes through None, '' and values containing `$`."""
    if path is None:
        return None
    path = os.fspath(path)
    if not path or "$" in path:
        return path
    return run_relative(path, base_dir, base_dir if run_root is None else run_root, dot=dot)


def cdpath_reset(shell="bash"):
    """The line that unsets CDPATH (cdpath in tcsh/csh) for a generated script."""
    return "unset cdpath" if os.path.basename(shell) in ("tcsh", "csh") else "unset CDPATH"


def _dq(expr):
    return '"' + expr.replace('"', '\\"') + '"'


class JobGen:
    dot_relative = staticmethod(dot_relative)
    run_relative = staticmethod(run_relative)
    spell = staticmethod(spell)

    def __init__(self, config_dir, machine=None, hostname=None, node_type=None):
        with open(os.path.join(config_dir, "machines.yaml")) as f:
            self.machines = yaml.safe_load(f)["machines"]
        with open(os.path.join(config_dir, "batch.yaml")) as f:
            self.batch_systems = yaml.safe_load(f)["batch_systems"]
        self.name = self._resolve(machine, hostname or socket.gethostname())
        self.m = self.machines[self.name]
        self.b = self.batch_systems[self.m["scheduler"]]
        self.nt, self.node_type_name = None, None
        nts = self.m.get("node_types")
        if nts:
            ntn = node_type or self.m.get("node_type_default") or next(iter(nts))
            if ntn not in nts:
                raise ValueError(f"Unknown node_type {ntn!r} for {self.name}; known: {sorted(nts)}")
            self.nt, self.node_type_name = nts[ntn], ntn

    def _cores(self):
        return (self.nt or {}).get("cores_per_node", self.m["cores_per_node"])

    def _modules(self, stack=None):
        # node_type modules win, then the stack's, then the machine's
        nt = self.nt or {}
        if "modules" in nt:
            return nt["modules"]
        stacks = self.m.get("stacks") or {}
        if stacks:
            if stack is None:
                if len(stacks) == 1:
                    stack = next(iter(stacks))
                else:
                    raise ValueError(
                        f"machine {self.name!r} defines stacks {sorted(stacks)}; pass stack=")
            if stack not in stacks:
                raise ValueError(
                    f"machine {self.name!r} has no {stack!r} stack (known: {sorted(stacks)})")
            return stacks[stack].get("modules", [])
        return self.m.get("modules", [])

    def _module_use(self, stack=None):
        stacks = self.m.get("stacks") or {}
        if stack is None and len(stacks) == 1:
            stack = next(iter(stacks))
        st = stacks.get(stack) or {}
        return (list(self.m.get("module_use", [])) + list(st.get("module_use", []))
                + list((self.nt or {}).get("module_use", [])))

    def _node_model(self):
        return (self.nt or {}).get("model")

    def _resolve(self, override, hostname):
        if override:
            if override not in self.machines:
                raise ValueError(f"Unknown machine {override!r}; known: {sorted(self.machines)}")
            return override
        # a catch-all '.*' machine is tried last
        ordered = sorted(self.machines, key=lambda n: self.machines[n]["nodename_regex"] == ".*")
        for n in ordered:
            if re.search(self.machines[n]["nodename_regex"], hostname):
                return n
        raise ValueError(f"No machine matches hostname {hostname!r}; pass machine= explicitly.")

    def sched(self, key, **fmt):
        """Scheduler environment token, e.g. sched('jobid') -> '$PBS_JOBID' ('' if undefined)."""
        v = (self.b.get("sched_env") or {}).get(key, "")
        return v.format(**fmt) if (v and fmt) else (v or "")

    def shebang(self):
        return f"#!/bin/{self.m['shell']}"

    def directive_lines(self, values):
        """Directive lines for {logical_name: value}; None, '' and False are omitted."""
        pre = self.b["directive_prefix"]
        out = []
        for name, tmpl in self.b["directives"].items():
            v = values.get(name)
            if tmpl is None or v is False or v in (None, "", "None"):
                continue
            out.append(f"{pre} {tmpl.format(v=v)}")
        return out

    def _derive(self, group):
        threads = int(group.get("threads", 1))
        tpn = group.get("tasks_per_node") or (self._cores() // threads)
        nodes = math.ceil(int(group["mpi_tasks"]) / tpn)
        return nodes, tpn, threads

    def queue_option(self, queue, key, default=None):
        """machines.<m>.queues.<queue>.<key>, or `default`."""
        qs = self.m.get("queues")
        if not isinstance(qs, dict) or queue in (None, ""):
            return default
        return (qs.get(str(queue)) or {}).get(key, default)

    def resource_lines(self, task_groups, walltime, queue=None):
        """The resource directive line(s). Queues with a mem_per_chunk (shared queues such as
        derecho develop) also get a per-chunk memory request."""
        r = self.b["resource"]
        if not r.get("chunk"):
            return []
        model = self._node_model()
        node_sel = r["node_selector"].format(v=model) if (model and r.get("node_selector")) else ""
        mem = self.queue_option(queue, "mem_per_chunk")
        mem_part = r["mem_part"].format(v=mem) if (mem and r.get("mem_part")) else ""
        chunks = []
        for g in task_groups:
            nodes, tpn, threads = self._derive(g)
            tp = r["threads_part"].format(threads=threads) if threads > 1 else ""
            chunks.append(r["chunk"].format(
                nodes=nodes, cpus_per_node=self._cores(), tasks_per_node=tpn,
                threads=threads, threads_part=tp, node_sel=node_sel) + mem_part)
        joined = r["join"].join(chunks)
        line = r["directive"].format(chunks=joined, walltime=walltime)
        pre = self.b["directive_prefix"]
        out = [f"{pre} {ln}" if not ln.startswith(pre) else ln for ln in line.split("\n")]
        if model and r.get("node_selector_directive"):
            out.append(f"{pre} {r['node_selector_directive'].format(v=model)}")
        return out

    def cdpath_reset(self):
        return cdpath_reset(self.m["shell"])

    def env_lines(self, env):
        if self.m["shell"] == "tcsh":
            return [f"setenv {k} {v}" for k, v in env.items()]
        return [f'export {k}="{v}"' for k, v in env.items()]

    def module_lines(self, stack=None):
        """module purge/use/load lines; `stack` picks a coupled stack's set (e.g. "grt")."""
        out = []
        if self.m.get("module_purge"):
            out.append(self.m["module_purge"])
        for d in self._module_use(stack):
            out.append(f"module use -a {d}")
        for mod in self._modules(stack):
            out.append(f"module load {mod}")
        return out

    def stack_option(self, stack, key, default=None):
        """machines.<m>.stacks.<stack>.<key>, or `default`."""
        return ((self.m.get("stacks") or {}).get(stack) or {}).get(key, default)

    def total_tasks(self, task_groups):
        return sum(int(g["mpi_tasks"]) for g in task_groups)

    def run_command(self, task_groups, exe, args, log=None):
        cmd = self.m["mpirun"].format(total_tasks=self.total_tasks(task_groups))
        cmd = f"{cmd} {exe} {args}".rstrip()
        if log:
            cmd += f" >&! {log}" if self.m["shell"] == "tcsh" else f" > {log} 2>&1"
        return cmd

    def submit_command(self, script, depends_on=None, subdir="."):
        """Submit command for `script`, run from `subdir` so the job's submit directory is its
        own; waits (afterok) on `depends_on` job ids or shell expressions."""
        sub = self.b["submit"]
        dep = ""
        if depends_on and self.b.get("dependency"):
            dep = self.b["dependency"].format(jobids=":".join(depends_on)) + " "
        cmd = f"{sub} {dep}{shlex.quote(dot_relative(script))}"
        if os.fspath(subdir) in ("", "."):
            return cmd
        return f"cd {shlex.quote(dot_relative(subdir))} && {cmd}"

    def submit_chain_script(self, chains, header=None, jobids=None):
        """Text of a fail-fast bash script that submits every job in order.

        chains: [(comment, var, [(subdir, script), ...], first_deps), ...]; each job waits on the
        previous one in its chain, the first on first_deps. jobids: file collecting script=id lines.
        """
        chained = bool(self.b.get("dependency"))
        dep = ""
        if chained:
            dep = "${2:+" + self.b["dependency"].format(jobids="$2") + "} "
        body = f'id=$(cd "$1" && {self.b["submit"]} {dep}"$3") || exit 1; '
        if chained:
            body += '[ -n "$id" ] || { echo "no job id for $3" >&2; exit 1; }; '
        if jobids:
            body += f'echo "${{3#./}}=$id" >> {shlex.quote(jobids)}; '
        L = ["#!/bin/bash"] + ([f"# {header}"] if header else []) + [
            "set -euo pipefail", cdpath_reset("bash"), 'cd "$(dirname "$0")"']
        if jobids:
            L.append(f": > {shlex.quote(jobids)}")
        L.append(f'submit() {{ local id; {body}echo "$id"; }}')
        for comment, var, jobs, first_deps in chains:
            if comment:
                L.append(f"# {comment}")
            for i, (subdir, script) in enumerate(jobs):
                deps = ":".join(first_deps or []) if i == 0 else f"${var}"
                sd = dot_relative("." if subdir in (None, "") else os.fspath(subdir))
                L.append(f"{var}=$(submit {shlex.quote(sd)} {_dq(deps)} "
                         f"{shlex.quote(dot_relative(os.fspath(script)))})")
                if not jobids:
                    L.append(f"echo ${var}")
        if jobids:
            L.append(f"cat {shlex.quote(jobids)}")
        return "\n".join(L) + "\n"
