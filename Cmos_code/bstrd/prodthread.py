import signal

from threading import Thread, Event
from .prodqueue import ProdQueue, Full, Locked


def prepend_signal(sig, func):
    handler = signal.getsignal(sig)

    def wrapper(*args, **kwargs):
        func(*args, **kwargs)
        handler(*args, **kwargs)

    signal.signal(sig, wrapper)



class ProdThread:
    """
    Upon call of start(),
    the provided func will be executed with the argument
    running = threading.Event()
    in a separate thread.

    The result is expected to be iterable
    and the yielded values will be filled into a queue.

    If the queue is full, the values will be dropped.

    The oldest entry can be retrieved/removed via get().

    The iterator should obey the state of the running Event.

    Calling stop() clears the running Event and joins the thread.
    """

    def __init__(self, func, maxsize=0):
        self.func = func

        self.thread = None
        self.queue = ProdQueue(maxsize=maxsize)
        self.get = self.queue.get

        self.running = Event()

        prepend_signal(signal.SIGINT, self.stop)


    def __repr__(self):
        tn = type(self).__name__
        running = "running" if self.running.is_set() else "stopped"
        return f"{tn}: {running}"


    def target(self):
        self.running.set()
        gen = self.func(self.running)
        for data in gen:
            try:
                self.queue.put_nowait(data)
            except (Full, Locked):
                pass
        self.running.clear()


    def start(self):
        if not self.thread:
            self.thread = thread = Thread(target=self.target)
            thread.start()

    def stop(self, *args): # signal() gives some extra args
        self.running.clear()
        if self.thread:
            self.thread.join()
            self.thread = None



