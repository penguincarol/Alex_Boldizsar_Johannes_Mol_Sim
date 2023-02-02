import pandas
import math
import matplotlib.pyplot as plt

# read the excel file (everything on "text" setting)
dataframe = pandas.read_excel('TestFile.xlsx')
print()

# how many different Particle amounts there are
amountOfmeasures = len(dataframe.index)
# the ParticleAmounts themselves
particles = []
# the time data
firstmeasure = []
secondmeasure = []

# loads the particle amounts in a list
for i in range(amountOfmeasures):
    particles += [dataframe["Particles"].loc[dataframe.index[i]]]

# loads the first measurements in a list and computes the mean of 1-3 measurements
for i in range(amountOfmeasures):
    added = 0
    counter = 0
    for j in range(1, 4):
        probe = dataframe['1.{}'.format(j)].loc[dataframe.index[i]]
        if not math.isnan(probe):
            counter += 1
            added += probe
    if added != 0:
        firstmeasure += [added / counter]

# loads the second measurements in a list and computes the mean of 1-3 measurements
for i in range(amountOfmeasures):
    added = 0
    counter = 0
    for j in range(1, 4):
        probe = dataframe['2.{}'.format(j)].loc[dataframe.index[i]]
        if not math.isnan(probe):
            counter += 1
            added += probe
    if added != 0:
        secondmeasure += [added / counter]

# plotting
plt.plot(particles, firstmeasure, color='green', linestyle='dashed', linewidth=1.5,
         marker='o', markerfacecolor='blue', markersize=7)
plt.plot(particles, secondmeasure, color='blue', linestyle='dashed', linewidth=1.5,
         marker='x', markerfacecolor='red', markersize=7)

# setting x and y axis range
# plt.ylim(1, 8)
# plt.xlim(1, 8)

# naming the x axis
plt.xlabel('Time')
# naming the y axis
plt.ylabel('y - axis')

# giving a title to the graph
plt.title('<Title>')
plt.legend(["first", "second"], loc="lower left")

# function to show the plot
plt.show()

# combining the list for the table + time unit
combilist = []
for i in range(len(firstmeasure)):
    combilist += [[str(firstmeasure[i]) + " ms"] + [str(secondmeasure[i]) + " ms"]]
# the colums/rows
val1 = ["Benchmark 1", "Benchmark 2"]
val2 = ["{} Particles".format(particles[i]) for i in range(5)]
val3 = combilist

# the table
fig, ax = plt.subplots()
ax.set_axis_off()
table = ax.table(
    cellText=val3,
    rowLabels=val2,
    colLabels=val1,
    rowColours=["lightblue"] * 10,
    colColours=["lightblue"] * 10,
    cellLoc='center',
    loc='upper left')
# the title for the table
ax.set_title('<Title>',
             fontweight="bold")

plt.show()
