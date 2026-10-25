"""Small producer-neutral command line for prepared visualization scenes."""
import argparse
import json
from pathlib import Path
import sys
import time

from .visualization import RenderSettings, SceneSequence, render_sequence


def parser():
    root = argparse.ArgumentParser(prog='topwater', description='Topwater workflows')
    commands = root.add_subparsers(dest='command', required=True)
    render = commands.add_parser('render', help='Render a prepared top-down hull/water animation')
    render.add_argument('scene', type=Path, help='Producer-neutral prepared SceneSequence .npz')
    render.add_argument('output', type=Path, help='Complete destination MP4 file path; fixed 800 x 800 pixels')
    render.add_argument('--settings', type=Path, help='JSON object of RenderSettings fields')
    return root


def main(argv=None):
    started = time.monotonic()
    args = parser().parse_args(argv)
    try:
        settings = RenderSettings(**({} if args.settings is None else json.loads(args.settings.read_text(encoding='utf-8-sig'))))
        loading = time.monotonic()
        scene = SceneSequence.load(args.scene)
        input_load_seconds = time.monotonic() - loading
        result = render_sequence(scene, args.output, settings)
        result['input_load_validation_seconds'] = input_load_seconds
        result['cli_main_wall_seconds'] = time.monotonic() - started
        result['cli_timing_note'] = 'main wall includes input and rendering; interpreter/package imports precede main and require an external process wall timer'
    except (ValueError, TypeError, ImportError, FileExistsError, OSError) as error:
        print(f'topwater render: {error}', file=sys.stderr)
        return 2
    print(json.dumps(result))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
