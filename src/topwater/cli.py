"""Directory-based generation, simulation and rendering workflows."""
import argparse
from pathlib import Path
import sys


def parser():
    root = argparse.ArgumentParser(prog='topwater', description='Topwater workflows')
    commands = root.add_subparsers(dest='command', required=True)
    generate = commands.add_parser('generate-sea',
        help='Generate explicit waves; docs/cli/generate-sea.md',
        epilog='Documentation: docs/cli/generate-sea.md')
    generate.add_argument('directory', nargs='?', type=Path, default=Path('.'),
                          help='Example directory (default: current directory)')
    generate.add_argument('--record', required=True, help='Packaged observed sea-state ID')
    generate.add_argument('--seed', type=int, help='Random seed; omitted seed is drawn and printed')
    solve = commands.add_parser('solve', help='Solve an example; docs/cli/solve.md',
                               epilog='Documentation: docs/cli/solve.md')
    solve.add_argument('directory', nargs='?', type=Path, default=Path('.'),
                       help='Example directory (default: current directory)')
    solve.add_argument('duration', type=float, help='Positive simulation duration in seconds')
    solve.add_argument('--solver', choices=('native', 'wecsim'), default='native',
                       help='Backend (default: native; currently unavailable)')
    render = commands.add_parser('render', help='Render result.npz; docs/cli/render.md',
                                epilog='Documentation: docs/cli/render.md')
    render.add_argument('directory', nargs='?', type=Path, default=Path('.'),
                        help='Example directory (default: current directory)')
    return root


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        if args.command == 'generate-sea':
            from . import sea_states
            from .workflow_data import Sea
            waves = sea_states.realize(args.record, seed=args.seed)
            sea = Sea(waves.amplitudes, waves.angular_frequencies, waves.directions, waves.phases)
            args.directory.mkdir(parents=True, exist_ok=True)
            destination = args.directory / 'sea.npz'
            sea.save(destination, overwrite=True)
            print(f'Wrote {destination}: {len(waves.amplitudes)} components; '
                  f'record={args.record}; seed={waves.seed}')
        elif args.command == 'solve':
            import math
            if not math.isfinite(args.duration) or args.duration <= 0:
                raise ValueError('duration must be positive and finite')
            if args.solver == 'native':
                raise RuntimeError('native coupled dynamics are not implemented; '
                                   'use --solver=wecsim with its documented dependencies')
            from .workflow_data import Boat, Sea
            from .wecsim_backend import solve
            boat = Boat.load(args.directory / 'boat.npz')
            sea = Sea.load(args.directory / 'sea.npz')
            result = solve(boat, sea, args.duration)
            destination = args.directory / 'result.npz'
            result.save(destination, overwrite=True)
            print(f'Wrote {destination}; existing video.mp4 was left unchanged')
        else:
            from .workflow_data import Result
            from .result_render import render_result
            result = Result.load(args.directory / 'result.npz')
            destination = args.directory / 'video.mp4'
            render_result(result, destination)
            print(f'Wrote {destination}')
    except KeyboardInterrupt:
        print(f'topwater {args.command}: interrupted; previous outputs preserved', file=sys.stderr)
        return 130
    except (ValueError, TypeError, ImportError, OSError, RuntimeError, ArithmeticError) as error:
        print(f'topwater {args.command}: {error}', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
