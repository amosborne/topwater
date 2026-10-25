"""Small producer-neutral command line for prepared visualization scenes."""
import argparse
import json
from pathlib import Path
import sys

from .visualization import RenderSettings, SceneSequence, render_sequence


def parser():
    root = argparse.ArgumentParser(prog='topwater', description='Topwater workflows')
    commands = root.add_subparsers(dest='command', required=True)
    render = commands.add_parser('render', help='Render prepared hull/water meshes using PyVista')
    render.add_argument('scene', type=Path, help='Producer-neutral prepared SceneSequence .npz')
    render.add_argument('output', type=Path, help='Complete destination MP4 file path; fixed 800 x 800 pixels')
    render.add_argument('--settings', type=Path, help='JSON object of RenderSettings fields')
    return root


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        settings = RenderSettings(**({} if args.settings is None else json.loads(args.settings.read_text(encoding='utf-8-sig'))))
        scene = SceneSequence.load(args.scene)
        result = render_sequence(scene, args.output, settings)
    except (ValueError, TypeError, ImportError, FileExistsError, OSError) as error:
        print(f'topwater render: {error}', file=sys.stderr)
        return 2
    print(json.dumps(result))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
