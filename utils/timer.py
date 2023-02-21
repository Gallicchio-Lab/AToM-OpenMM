from contextlib import contextmanager
import time

@contextmanager
def Timer(logger, message):
    logger(f"Started: {message}")
    start = time.monotonic()
    yield
    duration = time.monotonic() - start
    logger(f"Finished: {message} (duration: {duration} s)")
