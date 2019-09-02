import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import math
import time
import pdb

def plot_env(R,H):
	rect1 = mpatches.Rectangle(xy =(0.0,0.0) ,width =20 , height = H, alpha = 0.5)
	rect2 = mpatches.Rectangle(xy=(0,H), width = 20, height = R, alpha = 0.5)
	wedge = mpatches.Wedge(center = (20,H), r = R,  theta1 = 0 , theta2= 90, alpha=0.5)
	rect3 = mpatches.Rectangle(xy = (20,0), width = R, height = H, alpha =0.5)
	rect4 = mpatches.Rectangle(xy = (20+R, 0), width = 50, height = R, alpha = 0.5)
	ax.add_patch(wedge)
	ax.add_patch(rect1)
	ax.add_patch(rect2)
	ax.add_patch(rect3)
	ax.add_patch(rect4)

def plot_morphology(data, L, F, R,i):
	d = data[0]*100
	a = data[1]*100
	alpha = data[2]
	beta = data[3]
	theta = data[4]
	t = data[5]*100
	#print(t)

	S2 = (20+d, a)
	S1 = (S2[0]-L * math.cos(theta), S2[1]+L*math.sin(theta))
	S0 = (S1[0] - F * math.cos(theta+alpha), S1[1] + F*math.sin(theta+alpha))
	S3 = (S2[0]+F*math.cos(beta-theta), S2[1]+F*math.sin(beta-theta))
	#x, y = np.array([[S0[0],S1[0],S2[0]], [S0[1], S1[1], S2[1]]])
	
	x, y = np.array([[S0[0],S1[0],S2[0], S3[0]], [S0[1], S1[1], S2[1], S3[1]]])
	line = mlines.Line2D(x, y, lw=2., color = "grey", alpha=1)
	CircleS0 = mpatches.Circle(S0, 1, color = "orange",ec="none")
	CircleS1 = mpatches.Circle(S1, 1, color = "orange", ec="none")
	CircleS2 = mpatches.Circle(S2, 1, color = "orange", ec="none")
	CircleS3 = mpatches.Circle(S3, 1, color = "orange", ec="none")

	if (t>0):
		T = (S2[0]-t * math.cos(theta), S2[1]+t*math.sin(theta))
		CircleT = mpatches.Circle(T, 1, color = "red",  ec="none")
		ax.add_patch(CircleT)

	ax.add_patch(CircleS0)
	ax.add_patch(CircleS1)
	ax.add_patch(CircleS2)
	ax.add_patch(CircleS3)

	ax.add_line(line)
	plt.axis('scaled')
	plt.axis('equal')
	plt.savefig("./Plot/"+str(i)+".pdf",format = 'pdf')

file = open("PairData.txt","r")
Data = eval(file.readline())

print(len(Data))

i = 0
for data in Data:
	fig = plt.figure()
	ax = fig.add_subplot(111)
        if i == 47:
            print(data)
	plot_env(3.5, 9.5)
	plot_morphology(data , 14.5, 13.5, 3.5,i)
	i+=1
