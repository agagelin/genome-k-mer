from multiprocessing import Process, Manager
import time
import random
import sys

def task_(e):
    print("start")
    time.sleep(2)
    print("end")
    return e

def proc_task(results, i, e):
    results[i] = task_(e)

N = 2
manager = Manager()
results = manager.list([None]*N)
print(results)
jobs = list()
for i in range(N):
    job = Process(target=proc_task, args=(results, i, i))
    job.start()
    jobs.append(job)

for j in jobs:
    j.join()
    print(results)
print(results)

print("Hand made")
proc_task(results, 0, "what")
print(results)
