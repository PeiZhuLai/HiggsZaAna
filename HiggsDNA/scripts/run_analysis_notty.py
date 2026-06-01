"""
Thin launcher for HiggsDNA run_analysis that avoids the os.get_terminal_size()
crash when stdout is not a TTY (background / pipe), without needing a PTY
(`script`), which was injecting EINTR into dill.dump during the merge step.

We monkeypatch os.get_terminal_size to a fixed (200, 50) and then exec the
real run_analysis main with argv passed through.
"""
import os
import sys


def _fixed_terminal_size(*args, **kwargs):
    return os.terminal_size((200, 50))


os.get_terminal_size = _fixed_terminal_size

# Make run_analysis importable as a module path
sys.argv[0] = "scripts/run_analysis.py"
runpy_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else ".",
)

# Just exec the script in-process
script = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts/run_analysis.py"
with open(script) as f:
    code = compile(f.read(), script, "exec")
exec(code, {"__name__": "__main__", "__file__": script})
