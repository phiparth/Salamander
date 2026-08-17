#!/usr/bin/env python3
"""Find FoldX on this machine, check it runs, and print the one line that points
Salamander at it.

FoldX ships as a zip with no installer. Wherever you unzipped it - `foldx_windows` on your
Desktop, `foldx5` in Downloads, a renamed folder on drive D - this finds it. You do not have
to move anything.

    python find_foldx.py                    look for it and report
    python find_foldx.py --deep             also scan every drive (slower, thorough)
    python find_foldx.py --all              list every copy found, not just the best
    python find_foldx.py --copy-to C:/FoldX put a copy somewhere predictable
    python find_foldx.py --move-to C:/FoldX move it there instead of copying

Nothing is moved, copied or changed unless you pass --copy-to or --move-to.
"""
import argparse
import os
import shutil
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
from strip2 import foldx as fx  # noqa: E402

# what a working FoldX folder needs to bring along
BUNDLE = ("rotabase.txt", "molecules", "yasaraPlugin.py")


def human(p):
    try:
        return "%.1f MB" % (os.path.getsize(p) / (1024.0 * 1024.0))
    except OSError:
        return "?"


def report(exe):
    print("  executable : %s" % exe)
    print("  size       : %s" % human(exe))
    print("  folder     : %s" % os.path.dirname(os.path.abspath(exe)))
    c = fx.companions(exe)
    for name in ("rotabase.txt", "molecules"):
        print("  %-11s: %s" % (name, "present" if c[name] else "not in this folder"))
    if not c["rotabase.txt"] and not c["molecules"]:
        print("  WARNING: neither rotabase.txt nor molecules/ sits beside the executable.")
        print("           FoldX 4 needs rotabase.txt and will exit immediately without it.")
        print("           If you unzipped only the .exe, unzip the rest into this folder.")


def transfer(exe, dest, move):
    """Copy or move the executable and its support files into dest."""
    src_dir = os.path.dirname(os.path.abspath(exe))
    dest = os.path.abspath(dest)
    if os.path.abspath(src_dir) == dest:
        print("already in %s - nothing to do" % dest)
        return os.path.join(dest, os.path.basename(exe))

    if os.path.exists(dest) and os.listdir(dest):
        print("WARNING: %s already exists and is not empty." % dest)
        print("         Files with the same names will be overwritten.")
        ans = raw_input("continue? [y/N] ") if sys.version_info[0] == 2 \
            else input("continue? [y/N] ")
        if ans.strip().lower() not in ("y", "yes"):
            sys.exit("cancelled - nothing was changed")

    if not os.path.isdir(dest):
        os.makedirs(dest)

    op = shutil.move if move else shutil.copy2
    verb = "moved" if move else "copied"

    target = os.path.join(dest, os.path.basename(exe))
    op(exe, target)
    print("%s %s" % (verb, os.path.basename(exe)))

    for name in BUNDLE:
        s = os.path.join(src_dir, name)
        if not os.path.exists(s):
            continue
        d = os.path.join(dest, name)
        if os.path.isdir(s):
            if os.path.isdir(d):
                shutil.rmtree(d)
            (shutil.move if move else shutil.copytree)(s, d)
        else:
            op(s, d)
        print("%s %s" % (verb, name))

    print("\nFoldX is now in %s" % dest)
    return target


def env_lines(exe):
    p = os.path.abspath(exe)
    print("\nPoint Salamander at it. Run the ONE line matching your shell - you do not")
    print("need more than one:\n")
    if sys.platform == "win32":
        print("  Command Prompt, permanent (reopen the terminal afterwards):")
        print('    setx FOLDX "%s"' % p)
        print("\n  PowerShell, this session only:")
        print('    $env:FOLDX = "%s"' % p.replace("\\", "/"))
    else:
        print("  bash / zsh, this session only:")
        print("    export FOLDX=%s" % p)
        print("\n  to make it permanent, add that line to ~/.bashrc or ~/.zshrc")
    print("\nOr skip the variable entirely and pass it per command:")
    print("    --foldx %s" % p)


def main():
    p = argparse.ArgumentParser(description="find the FoldX executable")
    p.add_argument("--deep", action="store_true",
                   help="scan every fixed drive; use if the quick search fails")
    p.add_argument("--all", action="store_true", help="list every copy found")
    p.add_argument("--copy-to", metavar="DIR",
                   help="copy the executable and its support files to DIR")
    p.add_argument("--move-to", metavar="DIR", help="move them to DIR instead")
    p.add_argument("--path", help="check this path instead of searching")
    p.add_argument("--no-verify", action="store_true", help="skip running the executable")
    a = p.parse_args()

    if a.copy_to and a.move_to:
        sys.exit("pass either --copy-to or --move-to, not both")

    log = lambda *m: print(*m, flush=True)
    print("Searching for FoldX ...")
    exe, hits = fx.find(a.path, deep=a.deep, log=log)

    if not exe:
        print("\nFoldX not found.")
        if not a.deep:
            print("\nThe quick search covers drive roots, Program Files, your home folder,")
            print("Desktop, Downloads and Documents. If you unzipped it somewhere else, try:")
            print("    python find_foldx.py --deep")
        else:
            print("\nNothing matching 'foldx' anywhere on your fixed drives. Either it is not")
            print("unzipped yet, or it is on a network/removable drive. Point at it directly:")
            print("    python find_foldx.py --path D:/wherever/foldx_windows/foldx.exe")
        print("\nFoldX is a free academic download, but you must register first:")
        print("    https://foldxsuite.crg.eu/")
        return 1

    print("\nFound %d candidate%s.\n" % (len(hits), "" if len(hits) == 1 else "s"))
    if a.all and len(hits) > 1:
        for i, h in enumerate(hits, 1):
            print("[%d] %s  (%s)" % (i, h, human(h)))
        print("\nBest match - has its support files, or is newest:")
    report(exe)

    if not a.no_verify:
        print("\nRunning it ...")
        ok, line = fx.verify(exe)
        if ok:
            print("  OK: %s" % (line or "it starts"))
        else:
            print("  FAILED: %s" % line)
            print("  The file was found but will not run. Usual causes: it is not the")
            print("  command-line binary, the licence has expired, or on macOS/Linux it")
            print("  lacks the executable bit (chmod +x).")

    if a.copy_to or a.move_to:
        print("")
        exe = transfer(exe, a.copy_to or a.move_to, bool(a.move_to))

    env_lines(exe)
    return 0


if __name__ == "__main__":
    sys.exit(main())
