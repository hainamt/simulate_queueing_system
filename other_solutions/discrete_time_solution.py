import numpy as np

lambda_arrival = 4
lambda_service = 100
num_events = 5000

interarrival_times = [0.0] + np.random.exponential(scale=1/lambda_arrival, size=num_events - 1).tolist()
work_times = np.random.exponential(scale=1/lambda_service, size=num_events).tolist()

arrival_times = np.cumsum(interarrival_times).tolist()

print(f"Arrival Times: {arrival_times}")

served_times = []
finish_times = []

busy_until = -1.0
for arrival, work_time in zip(arrival_times, work_times):
    served_time = max(arrival, busy_until)
    busy_until = served_time + work_time
    served_times.append(served_time)
    finish_times.append(busy_until)

print(f"Time service begins: {served_times}")
print(f"Time finish works: {finish_times}")
print(f"Busy until: {busy_until}")

current_time = 0
num_customer_serving = []
num_customer_waiting = []
num_customer_sys = []
time_step=0.1
num_samples = int(busy_until/time_step)
print(f"Sample {num_samples} times")

for _ in range(num_samples):
    serving = sum(served_time <= current_time < served_time + work for (served_time, work) in zip(served_times, work_times))
    waiting = sum(arrival <= current_time < served_time for (arrival, served_time) in zip(arrival_times, served_times))
    
    num_customer_serving.append(serving)
    num_customer_waiting.append(waiting)
    num_customer_sys.append(serving + waiting)

    current_time += time_step