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


# Windows NTSTATUS codes that arrive as process exit codes. A hard crash produces no Python
# traceback at all, so the exit code is the only evidence of what happened.
CRASH = {
    3221225477: ("0xC0000005 ACCESS_VIOLATION",
                 "A segfault inside a compiled library. Usually a broken install or a\n"
                 "conflicting DLL of the same name found earlier on PATH - Anaconda,\n"
                 "MATLAB, Intel tools and GPU drivers all ship libiomp5md.dll / mkl DLLs."),
    3221225781: ("0xC0000135 DLL_NOT_FOUND",
                 "A dependent DLL is missing. Install the Visual C++ Redistributable."),
    3221225595: ("0xC000007B INVALID_IMAGE_FORMAT",
                 "A 32-bit DLL was loaded into 64-bit Python, or the reverse."),
    3221225725: ("0xC00000FD STACK_OVERFLOW", "Runaway recursion inside the library."),
    1073741795: ("0xC000001D ILLEGAL_INSTRUCTION",
                 "The CPU does not support an instruction the build requires, normally\n"
                 "AVX2. Install an older build that does not assume it."),
    3221225501: ("0xC000001D ILLEGAL_INSTRUCTION",
                 "The CPU does not support an instruction the build requires, normally AVX2."),
}


def probe(label, code, timeout=600):
    """Run a risky snippet in a SEPARATE interpreter and report how it ended.

    Importing torch or loading a 2.4 GB checkpoint can kill the process outright rather
    than raise. Doing it in a subprocess means this script survives to report the exit
    code, which is the only clue such a crash leaves behind.
    """
    import subprocess
    print("      %s ..." % label)
    try:
        r = subprocess.run([sys.executable, "-c", code], stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, timeout=timeout)
    except Exception as e:  # noqa: BLE001
        fail("%s could not be run: %s" % (label, e), "Unexpected; report this.")
        return False
    out = (r.stdout or b"").decode("utf-8", "replace").strip()
    rc = r.returncode
    if rc == 0:
        ok("%s: %s" % (label, out.splitlines()[-1] if out else "OK"))
        return True
    if rc in CRASH:
        name, why = CRASH[rc]
        fail("%s CRASHED - exit code %d (%s)." % (label, rc, name), why)
    else:
        tail = out.splitlines()[-4:] if out else ["(no output)"]
        fail("%s failed with exit code %d." % (label, rc), "\n".join(tail))
    return False


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

# FAMSA and PyTorch are both SIMD-compiled; on a CPU without AVX2 they do not raise, they
# abort the process. Reported always, because a silent exit gives nothing else to go on.
if sys.platform == "win32":
    try:
        import ctypes
        _f = ctypes.windll.kernel32.IsProcessorFeaturePresent
        _avx, _avx2 = bool(_f(39)), bool(_f(40))
        print("      cpu        : AVX=%s AVX2=%s" % (_avx, _avx2))
        if not _avx2:
            fail("This CPU does not support AVX2.",
                 "PyTorch wheels and pyfamsa are built assuming AVX2. Without it they\n"
                 "do not raise an error - the process dies with no traceback at all,\n"
                 "which is what a step 1 run ending at 'using N orthologs' looks like.\n"
                 "Workarounds:\n"
                 "  step 1:  add  --aligner biopython\n"
                 "  torch :  pip install \"torch==2.4.1\""
                 " --index-url https://download.pytorch.org/whl/cpu")
    except Exception:  # noqa: BLE001 - diagnostic only
        pass

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
# pyfamsa and Bio are compiled extensions, like torch - they fail together when the
# Visual C++ runtime is missing, and step 1 needs one of them to build its alignment
need = [("torch", "PyTorch"), ("transformers", "Transformers"), ("numpy", "NumPy"),
        ("pandas", "pandas"), ("sentencepiece", "SentencePiece"),
        ("pyfamsa", "FAMSA aligner, used by step 1"),
        ("Bio", "Biopython, step 1's fallback aligner")]
versions = {}
dll_trouble = []
for mod, label in need:
    try:
        m = __import__(mod)
        versions[mod] = getattr(m, "__version__", "?")
        ok("%-14s %s" % (mod, versions[mod]))
    except ImportError as e:
        # a genuinely absent package, or an installed one whose native library will not load
        msg = str(e)
        if "DLL" in msg or "dynamic link" in msg:
            dll_trouble.append(mod)
            fail("%s is installed but its native libraries will not load: %s" % (mod, msg),
                 "This is a system library problem, not a Python one. See section [3b].")
        else:
            fail("%s (%s) is not installed in THIS interpreter." % (mod, label),
                 "Activate the virtual environment, then:  pip install -r requirements.txt\n"
                 "If the prompt does not show (.venv), activation did not take effect.")
    except OSError as e:
        # WinError 1114: "A dynamic link library (DLL) initialization routine failed"
        dll_trouble.append(mod)
        fail("%s is installed but failed to initialise: %s" % (mod, e),
             "This is a system library problem, not a Python one. See section [3b].")
    except Exception as e:  # noqa: BLE001 - report anything rather than traceback
        fail("%s failed to import: %s: %s" % (mod, type(e).__name__, e),
             "Unexpected. Send this whole report to whoever maintains the repo.")

# ---- 3b. why a native library would not load (Windows) ----------------------------
if dll_trouble:
    print("\n[3b] Native library diagnosis (because %s failed above)"
          % ", ".join(dll_trouble))
    import platform
    bits = platform.architecture()[0]
    print("      python build : %s on %s" % (bits, platform.machine()))
    if bits != "64bit":
        fail("This Python is %s. PyTorch is 64-bit only." % bits,
             "Install 64-bit Python 3 and rebuild the virtual environment.")

    if sys.platform == "win32":
        import ctypes

        # torch links against the Visual C++ runtime, which Windows does not ship
        # and the wheel does not bundle
        missing_rt = []
        for dll in ("msvcp140.dll", "vcruntime140.dll", "vcruntime140_1.dll"):
            try:
                ctypes.WinDLL(dll)
            except OSError:
                missing_rt.append(dll)
        if missing_rt:
            fail("The Visual C++ runtime is missing: %s" % ", ".join(missing_rt),
                 "THIS IS ALMOST CERTAINLY THE CAUSE. Install the Microsoft Visual C++\n"
                 "Redistributable (x64), then open a NEW terminal and retry:\n"
                 "  https://aka.ms/vs/17/release/vc_redist.x64.exe")
        else:
            ok("Visual C++ runtime present (msvcp140, vcruntime140, vcruntime140_1).")

        # modern torch wheels are built assuming AVX2
        try:
            ipfp = ctypes.windll.kernel32.IsProcessorFeaturePresent
            avx = bool(ipfp(39))   # PF_AVX_INSTRUCTIONS_AVAILABLE
            avx2 = bool(ipfp(40))  # PF_AVX2_INSTRUCTIONS_AVAILABLE
            print("      cpu features : AVX=%s AVX2=%s" % (avx, avx2))
            if not avx2:
                fail("This CPU does not support AVX2.",
                     "Current PyTorch wheels assume AVX2 and abort while loading c10.dll.\n"
                     "Install an older build that does not:\n"
                     '  pip install "torch==2.4.1" --index-url https://download.pytorch.org/whl/cpu')
            else:
                ok("CPU supports AVX2, which current PyTorch wheels require.")
        except Exception:  # noqa: BLE001 - diagnostic only
            print("      cpu features : could not be queried")

        print("      If both checks above are OK, the remaining causes are a partial")
        print("      download (pip uninstall torch, then reinstall) or antivirus blocking")
        print("      the DLLs.")

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

# ---- 3c. does torch actually compute? -----------------------------------------------
# Importing torch and using it are different tests. A CPU or DLL problem often survives the
# import and kills the process at the first real matrix multiply, which is where ProtT5's
# load lands. Separating the two says whether the problem is torch or the checkpoint.
if "torch" in versions:
    print("\n[3c] PyTorch runtime")
    probe("torch imports and multiplies a 512x512 matrix",
          "import torch\n"
          "x = torch.randn(512, 512)\n"
          "y = (x @ x).sum().item()\n"
          "print('matmul ok, threads=%d' % torch.get_num_threads())",
          timeout=300)

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

    # exact published size of Rostlab/prot_t5_xl_half_uniref50-enc pytorch_model.bin.
    # A truncated download is the usual cause of an access violation (exit code
    # 3221225477 / 0xC0000005) while loading the model: torch reads past the end of the
    # file and the process dies with no Python error at all.
    EXPECT = 2416373051
    SHA256 = "7f51ba885541c7dc569d46b796af57cc7a2ba7945107dced4f19d1b5ec091157"
    for w in weights:
        n = int(round(have[w] * 1024.0 * 1024.0))
        exact = os.path.getsize(os.path.join(found, w))
        if w == "pytorch_model.bin" and exact != EXPECT:
            fail("%s is %d bytes; the published file is %d (%+d)."
                 % (w, exact, EXPECT, exact - EXPECT),
                 "The download did not complete, so the file is truncated. Delete the\n"
                 "folder and download it again:\n"
                 "  rmdir /s /q models\\prot_t5      (Windows)\n"
                 "  rm -rf models/prot_t5            (macOS/Linux)\n"
                 "then re-run the snapshot_download command from README section 3.\n"
                 "Verify afterwards - this must print %s:\n"
                 "  certutil -hashfile models\\prot_t5\\pytorch_model.bin SHA256"
                 % SHA256[:16])
        elif exact < 500 * 1024 * 1024:
            fail("%s is only %.1f MB - far too small." % (w, have[w]),
                 "A real ProtT5 encoder checkpoint is 2.25 GB. This is a truncated\n"
                 "download or a git-lfs pointer file. Download it again.")

    # the actual load, only if transformers imported
    if "transformers" in versions and not fails:
        print("")
        probe("ProtT5 loads (30-60 s, needs ~3 GB of RAM)",
              "from transformers import T5EncoderModel\n"
              "T5EncoderModel.from_pretrained(r'''%s''')\n"
              "print('loaded')" % found)
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
