from queue import Queue, Full 


class ProdQueue(Queue):

    def put_nowait(self, *args, **kwargs):
        if self.mutex.locked():
            raise Locked

        super().put_nowait(*args, **kwargs)


    def clear(self):
        with self.mutex:
            self.queue.clear()



class Locked(Exception):
    "Exception raised by ProdQueue.put_nowait() if mutex is locked."
    pass



