"""Locate the FoldX executable wherever it happens to live.

FoldX has no installer. You unzip it, and the folder is called whatever the archive was
called - `foldx_windows`, `foldx5`, `FoldX5_1`, `foldx_20251231`, or nothing in particular if
you renamed it. Nothing registers itself, nothing goes on PATH. So the only reliable way to
find it is to look.

Matching is deliberately loose: any file whose name contains "foldx" and looks executable,
in any folder. The distribution's YASARA plugin is excluded - it is an executable containing
"foldx" but is not the command-line binary.

Search order, cheapest first:

    1. an explicit path the caller passed
    2. $FOLDX
    3. PATH
    4. likely folders - drive roots, Program Files, home, Desktop, Downloads, Documents,
       and any *foldx* folder inside them
    5. a bounded walk of the home directory
    6. every fixed drive, only when asked (deep=True); minutes, not seconds

Used by step 3 and step 4, and by find_foldx.py at the repo root.
"""
import os
import subprocess
import sys

# Names that contain "foldx" but are not the command-line binary.
NOT_FOLDX = ("yasara", "plugin", "readme", "manual", "licence", "license")

# Directories never worth descending into: huge, or system-owned. Matched on the WHOLE
# folder name, not as substrings - "windows" as a substring would reject foldx_windows,
# which is exactly what the CRG archive unzips to.
SKIP = frozenset(("windows", "winsxs", "$recycle.bin", "system volume information",
                  "node_modules", ".git", "__pycache__", "site-packages", "temp", "tmp",
                  "cache", ".venv", "venv", "anaconda3", "miniconda3", "programdata",
                  "$windows.~bt", "recovery", "perflogs"))


def _is_exe(path):
    """Does this look like a runnable FoldX binary?"""
    name = os.path.basename(path).lower()
    if "foldx" not in name:
        return False
    if any(bad in name for bad in NOT_FOLDX):
        return False
    if not os.path.isfile(path):
        return False
    if sys.platform == "win32":
        return name.endswith(".exe")
    # POSIX: either the executable bit, or no extension at all
    return os.access(path, os.X_OK) and "." not in name.replace("foldx", "")


def _scan(d, out, depth=0, max_depth=2):
    """Collect FoldX-looking files in d, descending at most max_depth levels."""
    try:
        entries = list(os.scandir(d))
    except (OSError, PermissionError):
        return
    for e in entries:
        try:
            if e.is_file() and _is_exe(e.path):
                out.append(e.path)
            elif e.is_dir() and depth < max_depth:
                low = e.name.lower()
                named = "foldx" in low
                # a FoldX-named folder is always worth entering, whatever else it is called
                if not named and (low.startswith(".") or low in SKIP):
                    continue
                # descend freely near the top, but deeper only through *foldx* folders
                if depth == 0 or named:
                    _scan(e.path, out, depth + 1, max_depth)
        except (OSError, PermissionError):
            continue


def _drives():
    if sys.platform != "win32":
        return ["/opt", "/usr/local", "/usr/local/bin"]
    out = []
    for letter in "CDEFGHIJKLMNOPQRSTUVWXYZ":
        root = "%s:/" % letter
        if os.path.isdir(root):
            out.append(root)
    return out


def likely_dirs():
    """Folders worth checking before resorting to a full scan."""
    home = os.path.expanduser("~")
    out = []
    for d in _drives():
        out.append(d)
    for base in (home, os.path.join(home, "Desktop"), os.path.join(home, "Downloads"),
                 os.path.join(home, "Documents"), "C:/Program Files",
                 "C:/Program Files (x86)", "/opt", "/usr/local/bin",
                 os.path.join(home, "bin")):
        if os.path.isdir(base):
            out.append(base)
    return [d for d in out if os.path.isdir(d)]


def companions(exe):
    """Support files FoldX needs, which must sit beside the executable.

    FoldX 4 reads rotabase.txt at startup and dies without it. FoldX 5 has that data
    compiled in but still reads molecules/ for some commands. Reporting which are present
    turns 'RepairPDB produced nothing' into something diagnosable.
    """
    d = os.path.dirname(os.path.abspath(exe))
    found = {}
    for name in ("rotabase.txt", "molecules", "yasaraPlugin.py"):
        p = os.path.join(d, name)
        found[name] = os.path.exists(p)
    return found


def verify(exe, timeout=30):
    """Run it. Returns (ok, the version line).

    The flag is `-version`, one dash. FoldX parses `--version` as an option that takes a
    value and answers "the required argument for option '--version' is missing" - it still
    prints the banner, so it looks like it worked, which is how the wrong flag survives.
    """
    try:
        p = subprocess.Popen([exe, "-version"], stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, cwd=os.path.dirname(
                                 os.path.abspath(exe)) or None)
        try:
            out = p.communicate(timeout=timeout)[0]
        except Exception:  # noqa: BLE001 - old pythons lack the timeout kwarg
            p.kill()
            return False, "timed out"
    except OSError as e:
        return False, str(e)
    text = (out or b"").decode("utf-8", "replace")
    # the banner is a block of asterisks; the useful line is the one naming the version
    line = ""
    for ln in text.splitlines():
        if "foldx" in ln.lower():
            line = ln.strip().strip("*").strip()
            break
    if not line:
        rest = text.strip().splitlines()
        line = rest[0].strip() if rest else ""
    # "No pdbs for the run found" follows -version and is expected; we asked for no PDB
    return ("foldx" in text.lower()), line


def rank(paths):
    """Best candidate first: has companions, then newest."""
    def key(p):
        c = companions(p)
        return (-(c["rotabase.txt"] or c["molecules"]), -os.path.getmtime(p))
    return sorted(set(paths), key=key)


def find(explicit=None, deep=False, log=None):
    """Return (path or None, list of all candidates found).

    Honours the same precedence the pipeline documents, so a path passed on the command
    line or set in $FOLDX always wins over anything discovered on disk.
    """
    say = log or (lambda *a: None)

    if explicit:
        if os.path.isfile(explicit):
            return explicit, [explicit]
        say("the path given does not exist: %s" % explicit)

    env = os.environ.get("FOLDX")
    if env and os.path.isfile(env):
        say("found via $FOLDX")
        return env, [env]
    if env:
        say("$FOLDX is set but points at nothing: %s" % env)

    # PATH
    for d in os.environ.get("PATH", "").split(os.pathsep):
        if not d:
            continue
        try:
            for n in os.listdir(d):
                p = os.path.join(d, n)
                if _is_exe(p):
                    say("found on PATH")
                    return p, [p]
        except (OSError, PermissionError):
            pass

    hits = []
    say("searching likely folders ...")
    for d in likely_dirs():
        _scan(d, hits, max_depth=2)

    if not hits:
        say("searching your home folder ...")
        _scan(os.path.expanduser("~"), hits, max_depth=4)

    if not hits and deep:
        say("searching every drive - this takes a few minutes ...")
        for d in _drives():
            _scan(d, hits, max_depth=6)

    hits = rank(hits)
    return (hits[0] if hits else None), hits


def resolve(explicit=None, required=True):
    """Convenience for the pipeline steps: find it or exit with advice."""
    exe, _ = find(explicit)
    if exe:
        return exe
    if not required:
        return None
    sys.exit("FoldX not found.\n"
             "  Run:  python find_foldx.py\n"
             "  It searches your disk, checks the executable runs, and prints the one\n"
             "  line to set $FOLDX. Or pass --foldx <path to the executable>.")
