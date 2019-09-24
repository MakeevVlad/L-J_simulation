import plotly.offline as py
import plotly.graph_objs as go

# Create random data with numpy
X = []
Y = []

file = open('data.txt', 'r')
# file.readline()
# file.readline()
for line in file:
    row = line.split()
    X.append(row[0])
    Y.append(row[1])
file.close()
# Create a trace
trace = go.Scattergl(
x = X,
y = Y,
mode = 'markers',
marker = dict(size = 3)

)

data = [trace]

py.plot(data, filename='plot')
