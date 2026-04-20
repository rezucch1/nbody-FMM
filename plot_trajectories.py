#! /bin/python

import sys
INPUT_FILE = sys.argv[1] if len(sys.argv[1:]) \
   else "output.csv"

import matplotlib.pyplot as plt
import csv

labels = []
data = {}

with open(INPUT_FILE, newline='' ) as csvfile:
  csv_reader = csv.reader(csvfile, skipinitialspace=True)
  labels = csv_reader.__next__()
  for l in labels:
    data[l] = []
  for row in csv_reader:
    for idx in range(len(row)):
      data[labels[idx]].append(float(row[idx]))

for idx in range((len(labels) - 1) // 2):
  plt.plot(data[f'x{idx}'], data[f'y{idx}'], label=f'{idx}')

for idx in range((len(labels) - 1) // 2):
  plt.scatter(data[f'x{idx}'][-1], data[f'y{idx}'][-1])

plt.axis('equal')

plt.legend(loc='center right')
plt.show() 