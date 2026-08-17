#!/usr/bin/env python3
"""Report everything needed to diagnose a failed install, in one command.

    python check_env.py

Uses only the standard library, so it runs even when nothing else is installed - which is
the case it exists for. Nothing is downloaded and nothing is modified; it only looks.

Each check prints OK or a FAIL line naming the fix. Exit code is the number of failures.
"""
import os
import sys

OK, BAD = "  OK   ", "  FAIL "
fails = []


def fail(msg, fix):
    fails.append(msg)
    print(BAD + msg)
    for line in fix.splitlines():
        print("         -> " + line)


def ok(msg):
    print(OK + msg)


print("=" * 78)
print("Salamander install check")
print("=" * 78)

# ---- 1. the interpreter -------------------------------------------------------------
print("\n[1] Python")
print("      executable : %s" % sys.executable)
print("      version    : %s" % sys.version.split()[0])
print("      working dir: %s" % os.getcwd())

if sys.version_info < (3, 9):
    fail("Python %d.%d is too old; Salamander needs 3.9 or newer."
         % sys.version_info[:2],
         "On Windows the bare name `python` often points at some other program's\n"
         "bundled Python (MGLTools, ArcGIS, Cygwin). Use `py -3` instead of `python`,\n"
         "or activate the virtual environment from section 2 of the README.")
else:
    ok("Python %s is new enough." % sys.version.split()[0])

in_venv = sys.prefix != getattr(sys, "base_prefix", sys.prefix)
if in_venv:
    ok("Running inside a virtual environment (%s)." % sys.prefix)
else:
    print("      note   : not inside a virtual environment. That is allowed, but the")
    print("               README's `.venv` keeps these packages away from your system.")

# ---- 2. are we in the right folder? -------------------------------------------------
print("\n[2] Repository layout")
here = os.path.dirname(os.path.abspath(__file__))
if os.path.abspath(os.getcwd()) != here:
    fail("You are not in the Salamander folder.",
         "Every command in the README uses paths relative to the repo root.\n"
         "Run this first:  cd %s" % here)
else:
    ok("Current directory is the repo root.")

for rel in ("pipeline/step1_conserved_sites.py", "src/strip2/core.py", "requirements.txt"):
    if os.path.exists(os.path.join(here, rel)):
        ok("found %s" % rel)
    else:
        fail("missing %s" % rel,
             "The download is incomplete. Re-clone the repository.")

# ---- 3. dependencies ----------------------------------------------------------------
print("\n[3] Python packages")
need = [("torch", "PyTorch"), ("transformers", "Transformers"), ("numpy", "NumPy"),
        ("pandas", "pandas"), ("sentencepiece", "SentencePiece")]
versions = {}
for mod, label in need:
    try:
        m = __import__(mod)
        versions[mod] = getattr(m, "__version__", "?")
        ok("%-14s %s" % (mod, versions[mod]))
    except ImportError:
        fail("%s (%s) is not installed in THIS interpreter." % (mod, label),
             "Activate the virtual environment, then:  pip install -r requirements.txt\n"
             "If the prompt does not show (.venv), activation did not take effect.")

# transformers 5 dropped support for the .bin checkpoints ProtT5 ships as
tv = versions.get("transformers")
if tv:
    try:
        major = int(tv.split(".")[0])
    except ValueError:
        major = 0
    if major >= 5:
        fail("transformers %s is too new." % tv,
             "ProtT5 is distributed as pytorch_model.bin, and transformers 5 removed\n"
             "support for .bin checkpoints. Pin it back:\n"
             '  pip install "transformers>=4.30,<5"')
    else:
        ok("transformers major version %d supports ProtT5's .bin weights." % major)

# transformers >= 4.56 refuses to torch.load a .bin unless torch is >= 2.6
tvv, tor = versions.get("transformers"), versions.get("torch")
if tvv and tor:
    def tup(s):
        out = []
        for p in s.split("+")[0].split("."):
            try:
                out.append(int(p))
            except ValueError:
                break
        return tuple(out)
    if tup(tvv) >= (4, 56) and tup(tor) < (2, 6):
        fail("torch %s is too old for transformers %s to load a .bin model." % (tor, tvv),
             "transformers 4.56+ refuses .bin checkpoints on torch below 2.6.\n"
             "Either:  pip install --upgrade torch\n"
             'Or:      pip install "transformers>=4.30,<4.56"')
    else:
        ok("torch %s and transformers %s are a workable pair." % (tor, tvv))

# ---- 4. the ProtT5 weights ----------------------------------------------------------
print("\n[4] ProtT5 weights")
plm = os.environ.get("PLM")
print("      PLM env var: %s" % (plm if plm else "(not set)"))

cand = [p for p in (plm, os.path.join(here, "models", "prot_t5"), "models/prot_t5") if p]
found = None
for p in cand:
    if os.path.isdir(p):
        found = p
        break

if not found:
    fail("No ProtT5 folder found.",
         "Looked in: %s\n"
         "Either download it (README section 3, Option B), or drop the --plm flag\n"
         "entirely and let step 2 pull it from Hugging Face automatically."
         % ", ".join(cand))
else:
    print("      folder     : %s" % os.path.abspath(found))
    want = {"config.json": 0.0005, "pytorch_model.bin": 1000.0, "spiece.model": 0.2,
            "tokenizer_config.json": 0.0}
    have = {}
    for f in os.listdir(found):
        fp = os.path.join(found, f)
        if os.path.isfile(fp):
            have[f] = os.path.getsize(fp) / (1024.0 * 1024.0)

    weights = [f for f in have if f.endswith((".bin", ".safetensors"))]
    for f in sorted(want):
        if f in have:
            ok("%-24s %8.1f MB" % (f, have[f]))
        elif f == "pytorch_model.bin" and weights:
            ok("weights present as %s" % ", ".join(weights))
        else:
            fail("%s is missing from %s" % (f, found),
                 "The download did not finish. Delete the folder and re-run the\n"
                 "snapshot_download command from README section 3.")

    incomplete = [f for f in have if f.endswith((".incomplete", ".lock"))]
    if incomplete:
        fail("Found %d partial download file(s), e.g. %s"
             % (len(incomplete), incomplete[0]),
             "The download was interrupted. Delete the folder and download again.")

    for w in weights:
        if have[w] < 500:
            fail("%s is only %.1f MB - far too small." % (w, have[w]),
                 "A real ProtT5 encoder checkpoint is roughly 2.2-2.5 GB. This is a\n"
                 "truncated download or a git-lfs pointer file. Download it again.")

    # the actual load, only if transformers imported
    if "transformers" in versions and not fails:
        print("\n      loading it (this takes 30-60 s and needs ~3 GB of RAM) ...")
        try:
            from transformers import T5EncoderModel
            T5EncoderModel.from_pretrained(found)
            ok("ProtT5 loaded successfully.")
        except Exception as e:  # noqa: BLE001 - the message is the whole point
            fail("ProtT5 failed to load: %s: %s" % (type(e).__name__, e),
                 "Send this line to whoever maintains the repo - it names the cause.")
    elif "transformers" in versions:
        print("\n      skipping the load test until the failures above are fixed.")

# ---- summary -----------------------------------------------------------------------
print("\n" + "=" * 78)
if not fails:
    print("All checks passed. Continue with README section 4 (FoldX).")
else:
    print("%d check(s) failed:" % len(fails))
    for f in fails:
        print("  - %s" % f)
print("=" * 78)
sys.exit(len(fails))
