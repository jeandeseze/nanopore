import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


plt.figure()
plt.plot(
    [100,75,50,25,0][::-1],
    [100*6624/9017,100*6330/17303,100*3210/16527,100*989/8459,100*1041/19290][::-1],
    'bo'
)
plt.ylabel("percentage of A in sequence")
plt.xlabel("capping percentage")
plt.savefig("test.png")


plt.figure()
plt.plot(
    [100,75,50,25,0][::-1],
    [0.936, 0.954,0.968, 0.973,0.980][::-1],
    'bo'
)
plt.ylabel("percentage of G in sequence")
plt.xlabel("capping percentage")
plt.savefig("test1.png")

