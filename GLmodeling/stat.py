import plotly.offline as py
import plotly.graph_objs as go


# Create random data with numpy
Y = []

file = open('EnVStime.txt', 'r')
# file.readline()
# file.readline()
X = [int(i) for i in range(1, 7000)]
i = 0
for line in file:
    row = line.split()
    Y.append(row[0])
file.close()
# Create a trace
trace = go.Scattergl(
x = X,
y = Y,
mode = 'lines+markers',
marker = dict(size = 3)

)

data = [trace]

py.plot(data, filename='EnVStime')
