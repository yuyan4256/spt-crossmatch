#!/usr/bin/env python
"""spt — one terminal entry point for this repo's scripts.

Usage:
  spt <command> [args...]     run a command; args go to the script unchanged
  spt <command> --help        that command's own help
  spt list                    every command with a one-line description
  spt run <path> [args...]    any other script under scripts/, by relative path
                              e.g. spt run null_resample/rescore_floors.py

Runs the script with the same python this file was launched with, so alias it
to the anaconda interpreter in ~/.zshrc:

  alias spt="$HOME/anaconda3/bin/python $HOME/Desktop/SPT-crossmatch-clean/scripts/spt.py"

Scripts that rewrite outputs/v4_crossmatch_table.csv in place (add_milliquas,
add_decaps, update_sigma_wan2025, ...) are deliberately NOT given short names —
reach them through `spt run` so the table is never touched by a typo.
"""
import ast
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

COMMANDS = {
    'plot':        'plot_source_panels.py',
    'plot-csc2':   'plot_csc2_cutout_unwise.py',
    'precache':    'precache_fits.py',
    'build-table': 'build_v4_crossmatch_table.py',
}


def _first_doc_line(path):
    try:
        with open(path) as fh:
            doc = ast.get_docstring(ast.parse(fh.read()))
        return (doc or '').strip().splitlines()[0]
    except Exception:
        return ''


def usage():
    print('spt — entry point for SPT-crossmatch-clean/scripts\n')
    print('commands:')
    for name, fn in COMMANDS.items():
        print(f'  {name:12s} {_first_doc_line(os.path.join(HERE, fn))}')
    print(f'  {"run":12s} any other script by path, e.g. '
          'spt run null_resample/rescore_floors.py')
    print(f'  {"list":12s} this list')
    print('\n`spt <command> --help` shows that command\'s own options.')


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ('list', 'help', '-h', '--help'):
        usage()
        return
    cmd, rest = sys.argv[1], sys.argv[2:]
    if cmd == 'run':
        if not rest:
            raise SystemExit('usage: spt run <path-under-scripts/> [args...]')
        script, rest = os.path.join(HERE, rest[0]), rest[1:]
    elif cmd in COMMANDS:
        script = os.path.join(HERE, COMMANDS[cmd])
    else:
        usage()
        raise SystemExit(f"\nunknown command {cmd!r}")
    if not os.path.isfile(script):
        raise SystemExit(f'no such script: {script}')
    # Replace this process with the script, same interpreter, args untouched.
    os.execv(sys.executable, [sys.executable, script] + rest)


if __name__ == '__main__':
    main()
