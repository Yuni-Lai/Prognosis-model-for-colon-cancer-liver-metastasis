# -*- coding: utf-8 -*-
from concurrent.futures import ProcessPoolExecutor
import concurrent.futures
import time

start = time.perf_counter()


def do_something(seconds):
    print(f'Sleeping {seconds} second(s)...')
    time.sleep(seconds)
    return f'Done Sleeping...{seconds}'
#do_something(5)


secs = [5,4,3,2,1]
if __name__ == "__main__":
    #results = executor.map(do_something, secs)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        secs = [5, 4, 3, 2, 1]
        results = executor.map(do_something, secs)
    for result in results:
        print(result)
    
    finish = time.perf_counter()
    
    print(f'Finished in {round(finish-start, 2)} second(s)')