import os
import subprocess
import time
import json
import shutil

nThreads = 8

threads = []
logfiles = []

cwd = os.getcwd()

shutil.rmtree("./log", ignore_errors=True)
os.makedirs("./log/", exist_ok=True)

for t in range(nThreads):
    threadId = "th" + str(t)
    os.makedirs("./log/" + threadId, exist_ok=True)
    os.makedirs("./log/" + threadId + "/coords", exist_ok=True)

    fLog = open("./log/" + threadId + "/log.txt", "w")
    thread = subprocess.Popen(["./mycell", threadId, "config.json"], cwd=cwd, stdout=fLog)
    #stdout, stderr = p.communicate()
    threads.append(thread)
    logfiles.append(fLog)

for t in range(nThreads):
    threads[t].wait()

for t in range(nThreads):
    logfiles[t].close()
