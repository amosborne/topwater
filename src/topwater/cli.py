"""Directory-based generation, simulation and rendering workflows."""
import argparse
from pathlib import Path
import sys


def parser():
    root = argparse.ArgumentParser(prog='topwater', description='Topwater workflows')
    commands = root.add_subparsers(dest='command', required=True)
    boat = commands.add_parser('boat', help='Generate a standard boat; docs/workflow.md',
                               epilog='Documentation: docs/workflow.md')
    boat.add_argument('directory', nargs='?', type=Path, default=Path('.'),
                      help='Example directory (default: current directory)')
    boat.add_argument('--length', required=True, type=float, help='Hull length in metres')
    boat.add_argument('--beam', required=True, type=float, help='Maximum hull beam in metres')
    boat.add_argument('--deadrise', required=True, type=float, help='Transom deadrise in degrees')
    boat.add_argument('--mass', required=True, type=float, help='Total operating mass in kilograms')
    generate = commands.add_parser('sea',
        help='Generate explicit waves; docs/workflow.md',
        epilog='Documentation: docs/workflow.md')
    generate.add_argument('directory', nargs='?', type=Path, default=Path('.'),
                          help='Example directory (default: current directory)')
    generate.add_argument('--record', required=True, help='Packaged observed sea-state ID')
    generate.add_argument('--seed', type=int, help='Random seed; omitted seed is drawn and recorded in sea.log')
    solve = commands.add_parser('solve', help='Solve an example; docs/workflow.md',
                               epilog='Documentation: docs/workflow.md')
    solve.add_argument('directory', nargs='?', type=Path, default=Path('.'),
                       help='Example directory (default: current directory)')
    from topwater.solve.settings import DEFAULTS
    from topwater.solve.settings import DEFAULT_DWELL_S
    solve.add_argument('--dwell', type=float, default=DEFAULT_DWELL_S,
                       help='Positive full-strength dwell in seconds; total includes rest/ramp on both sides (default: 100)')
    solve.add_argument('--rest', type=float, default=10.,
                       help='Quiet target interval on each side in seconds (default: 10)')
    solve.add_argument('--ramp', type=float, default=5.,
                       help='Quintic ramp on each side in seconds (default: 5)')
    solve.add_argument('--solver', choices=('native', 'wecsim'), default='native',
                       help='Backend (default: native; wecsim requires MATLAB)')
    solve.add_argument('--device', default='auto',
                       help='Native reconstruction: auto, cpu, gpu, or GPU name (default: auto)')
    controls = solve.add_argument_group('calculation and output controls')
    for flag, field, kind, units in (
        ('time-step', 'time_step_s', float, 'Motion timestep, seconds'),
        ('packet-tolerance', 'packet_tolerance_m', float, 'Maximum packet envelope fitting error, metres')):
        default = getattr(DEFAULTS, field)
        controls.add_argument('--'+flag, dest=field, type=kind,
            default=default,
            help=f'{units} (default: {default if default is not None else "auto"})')

    render = commands.add_parser('render', help='Render result.npz; docs/workflow.md',
                                epilog='Documentation: docs/workflow.md')
    render.add_argument('directory', nargs='?', type=Path, default=Path('.'),
                        help='Example directory (default: current directory)')
    return root



def _solve_timing(args):
    """Validate CLI timing before loading data or entering a solver backend."""
    import math
    from topwater.solve.settings import DEFAULTS
    time_step = getattr(args, "time_step_s", DEFAULTS.time_step_s)
    for name, value, positive in (('dwell', args.dwell, True), ('ramp', args.ramp, True),
                                   ('rest', args.rest, False)):
        if not math.isfinite(value) or (value <= 0 if positive else value < 0):
            raise ValueError(f'{name} must be {"positive" if positive else "nonnegative"} and finite')
        steps = value/time_step
        if not math.isfinite(steps) or not math.isclose(value, round(steps)*time_step,
                rel_tol=0., abs_tol=8*sys.float_info.epsilon*max(abs(value), time_step)):
            raise ValueError(f'{name} must align to the selected {time_step:g}s timestep')
    total = 2*args.rest+2*args.ramp+args.dwell
    if not math.isfinite(total):
        raise ValueError('total packet duration must be finite')
    return total, total/2, time_step

def main(argv=None):
    args = parser().parse_args(argv)
    try:
        from topwater.command_log import CommandLog
        names = {'boat': 'boat.log', 'sea': 'sea.log',
                 'solve': 'result.log', 'render': 'video.log'}
        details = ' '.join(f'{name}={value}' for name, value in vars(args).items()
                           if name != 'command')
        with CommandLog(args.directory / names[args.command], args.command, details) as log:
            if args.command == 'boat':
                from topwater.boat.generation import generate_boat
                log.event('INFO', 'STAGE generate convex hull, integrate mass and validate flotation')
                log.event('INFO', 'MODEL offshore-warped-v/2 control-point convex envelope; '
                    'structure_kg=0.1*length_ft^3+70; hull/deck=70/30; '
                    'deck_depth_fraction=0.6; payload_split=50/50; '
                    'payload_xyz_over_length=(0.05,0,0.03),(0.50,0,0.20); '
                    'x_from_transom; z_from_lowest_keel_datum')
                boat = generate_boat(args.length, args.beam, args.deadrise, args.mass)
                log.event('INFO', f'MASS cg_m={boat.cg_m.tolist()} '
                          f'inertia_cg_kg_m2={boat.inertia_cg_kg_m2.tolist()}')
                destination = args.directory / 'boat.npz'
                log.event('INFO', 'STAGE publish validated boat archive')
                boat.save(destination, overwrite=True)
                log.published(destination)
                print(f'Wrote {destination}: {len(boat.vertices_m)} vertices, '
                      f'{len(boat.faces)} physical panels; mass={boat.mass_kg:g} kg')
            elif args.command == 'sea':
                from topwater.sea import states as sea_states
                from topwater.workflow_data import Sea
                log.event('INFO', 'STAGE generate observed sea realization')
                waves = sea_states.realize(args.record, seed=args.seed)
                sea = Sea(waves.amplitudes, waves.angular_frequencies, waves.directions, waves.phases)
                destination = args.directory / 'sea.npz'
                sea.save(destination, overwrite=True)
                log.published(destination)
                log.event('INFO', f'SEA record={args.record} seed={waves.seed}')
                print(f'Wrote {destination}: {len(waves.amplitudes)} components; record={args.record}')
            elif args.command == 'solve':
                log.event('INFO', 'STAGE validate solve inputs and backend')
                from dataclasses import asdict, fields
                from topwater.solve.settings import DEFAULTS
                from topwater.solve.settings import SolveOptions
                calculation = SolveOptions(**{field.name: getattr(args, field.name, getattr(DEFAULTS, field.name)) for field in fields(SolveOptions)})
                total, midpoint, timestep = _solve_timing(args)
                log.event('INFO', 'OPTIONS '+str(asdict(calculation)))
                log.event('INFO', f'PACKET dwell_s={args.dwell:g} rest_each_s={args.rest:g} '
                    f'ramp_each_s={args.ramp:g} total_s={total:g} '
                    f'dwell_midpoint_s={midpoint:g} time_step_s={timestep:g}')
                log.event('INFO', 'SEA phases_rad specify target carrier phases DIRECTLY at dwell midpoint; '
                    'phase zero means a target carrier crest at the solver wave reference point; '
                    'finite Fourier truncation may depart from that target')
                from topwater.workflow_data import Boat
                from topwater.workflow_data import Sea
                from topwater.solve import simulate as solve
                from topwater.solve.outputs import save_solve_outputs
                boat = Boat.load(args.directory / 'boat.npz')
                sea = Sea.load(args.directory / 'sea.npz')
                log.event('INFO', f'STAGE solve motion and total surface solver={args.solver}')
                options = dict(options=calculation, solver=args.solver)
                if args.solver == 'native':
                    options.update(device=args.device)
                elif args.device != 'auto':
                    raise ValueError('--device requires --solver=native')
                result = solve(boat, sea, args.dwell, rest_s=args.rest, ramp_s=args.ramp, **options)
                log.event('INFO', 'STAGE prepare and publish result archive and motion/sea-spectrum plot')
                destination, plot = save_solve_outputs(result, args.directory, sea=sea)
                log.published(destination)
                log.published(plot)
                print(f'Wrote {destination} and {plot}; existing video.mp4 was left unchanged')
            else:
                from topwater.workflow_data import Result
                from topwater.render.workflow import render_result
                from topwater.render.scene import RenderSettings
                log.event('INFO', 'STAGE load self-contained result')
                result = Result.load(args.directory / 'result.npz')
                destination = args.directory / 'video.mp4'
                settings = RenderSettings()
                log.event('INFO', 'STAGE render and verify complete video')
                timings = render_result(result, destination, settings=settings)
                log.published(destination)
                log.event('INFO',
                    f'VIDEO frames={timings["frames"]} fps={timings["fps"]:g} '
                    f'size={timings["size"][0]}x{timings["size"][1]} '
                    f'source_end_s={timings["source_end_s"]:.3f} '
                    f'encoded_duration_s={timings["video_duration_s"]:.3f} verified=true')
                timing_keys = ('preparation_seconds', 'setup_seconds', 'geometry_seconds',
                               'draw_seconds', 'encoder_start_seconds', 'encoder_send_seconds',
                               'encoder_close_seconds', 'verification_seconds', 'elapsed_s')
                log.event('INFO', 'TIMING render ' + ' '.join(
                    f'{key.removesuffix("_seconds")}_s={timings[key]:.3f}'
                    if key != 'elapsed_s' else f'elapsed_s={timings[key]:.3f}'
                    for key in timing_keys))
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
