"""Stage progress and bounded Capytaine warning reporting."""
from functools import wraps
import logging
import sys
from time import perf_counter


class _StageProgress:
    """Wall-time progress at existing work boundaries; never a physical clock."""
    def __init__(self, name, total, unit):
        self.name, self.total, self.unit = name, int(total), unit
        self.started = self.last_report = perf_counter()
        print(f'STAGE {name} started total={self.total} unit={unit}', file=sys.stderr)

    def update(self, completed):
        now = perf_counter()
        if now-self.last_report >= 5.:
            print(f'PROGRESS {self.name} completed={completed}/{self.total} '
                  f'unit={self.unit} elapsed_s={now-self.started:.1f}', file=sys.stderr)
            self.last_report = now

    def finish(self):
        print(f'STAGE {self.name} complete completed={self.total}/{self.total} '
              f'unit={self.unit} elapsed_s={perf_counter()-self.started:.1f}', file=sys.stderr)

class _CapytaineWarningFilter(logging.Filter):
    """Keep the first physical warning of each known family, then summarize."""
    def __init__(self):
        super().__init__()
        self.counts = {}

    def filter(self, record):
        text = record.getMessage()
        family = next((prefix for prefix in ('Mesh resolution', 'Irregular frequencies', 'Water depth')
                       if text.startswith(prefix)), None)
        if family is None:
            return True
        self.counts[family] = self.counts.get(family, 0)+1
        return self.counts[family] == 1

def _summarize_capytaine_warnings(function):
    @wraps(function)
    def call(*args, **kwargs):
        logger = logging.getLogger('capytaine.bem.problems_checks')
        summary = _CapytaineWarningFilter()
        logger.addFilter(summary)
        try:
            return function(*args, **kwargs)
        finally:
            logger.removeFilter(summary)
            for family, count in summary.counts.items():
                if count > 1:
                    print(f'Capytaine: {family} warning applies to {count} solves in '
                          f'{function.__name__}; first warning retained. These limitations remain unqualified.',
                          file=sys.stderr)
    return call
