from contextlib import contextmanager
import time

# class Timer2:

#     def __init__(self, logger, message):
#         self.logger = logger
#         self.message = message

#     def __enter__(self):
#         self.start = time.monotonic()
#         self.logger(f"Started: {self.message}")
#         return None

#     def __exit__(self, exc_type, exc_value, exc_tb):
#         duration = time.monotonic() - self.start
#         self.logger(f"Finished: {self.message} (duration: {duration} s)")

@contextmanager
def Timer(logger, message):
    logger(f"Started: {message}")
    start = time.monotonic()
    yield
    duration = time.monotonic() - start
    logger(f"Finished: {message} (duration: {duration} s)")
