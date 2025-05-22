from signal import signal, SIGINT, SIGTERM
import sys
import warnings


class TerminationGuard:

    def __enter__(self):
        self.terminate = False
        self.sigint = signal(SIGINT, self)
        self.sigterm = signal(SIGTERM, self)

    def __call__(self, *_):
        self.terminate = True
        warnings.warn("The process termination is delayed until a critical operation is completed")

    def __exit__(self, *_):
        signal(SIGINT, self.sigint)
        signal(SIGTERM, self.sigterm)
        if self.terminate:
            sys.exit(0)
