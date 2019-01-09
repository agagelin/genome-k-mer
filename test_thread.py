from threading import Thread
from time import sleep

class TestThread(Thread):

    def run(self):
        self.what = False
        for i in range(3):
            sleep(1)
        self.what = True

N = 4
thds = list()
for i in range(N):
    thd = TestThread()
    thd.start()
    thds.append(thd)

for thd in thds:
    thd.join()
    print(thd.what)
