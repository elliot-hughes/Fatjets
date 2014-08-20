from ROOT import *
import math
import sys

#msg_test = 1150
#msq_test = 150
msg_test = 850
msq_test = 250
R_test = 1.2
trim_R_test = 0.7
trim_ptf_test = 0.1

do_true = False
do_prune = False
do_trim = True
save_plots = True

Rs = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
Rs = [0.8, 1.2, 1.6, 2.0,]#  2.4, 2.8]
#Rs = [1.0, 1.1, 1.2, 1.3]
#Rs = [1.2, 1.4, 1.6]
#Rs = [1.2, 1.4]
Rs = [1.2]

points = [
	(850, 250),
	(1150, 150),
]
points = [
	(850, 250)
]

colors = [kBlack, kRed, kBlue, kSpring+2, kViolet, kGray+2, kOrange-6, kCyan+2]

cut1_dr = 0.1
cut2_dr_low = 2.75
cut2_dr_high = 3.75

gROOT.SetBatch()
#gStyle.SetOptStat("eroum")
gStyle.SetOptStat(0)
c1 = TCanvas("c1", "c1", 1000, 500)
c1.SetCanvasSize(500, 500)

h1 = {#	name			x-axis label		bins	xmin	xmax	W		S		E		N
	"njets":	["", "Number Of Fatjets",		60,		0,		60,		0.7,	0.65,	0.85,	0.85, "r"],
	"m0":		["", "m^{0} (GeV)",				50,		0,		1000,	0.7,	0.65,	0.85,	0.85, "r"],
	"m0_cut1":	["", "m^{0} (GeV)",				50,		0,		1000,	0.7,	0.65,	0.85,	0.85, "r"],
	"m0_cut2":	["", "m^{0} (GeV)",				50,		0,		1000,	0.7,	0.65,	0.85,	0.85, "r"],
	"m1":		["", "m^{1} (GeV)",				50,		0,		1000,	0.7,	0.65,	0.85,	0.85, "r"],
	"pt0":		["", "p_{T} (GeV)",				50,		0,		1000,	0.7,	0.65,	0.85,	0.85, "r"],
	"dr_fj":	["", "#Delta R^{01}",			50,		0,		5,		0.7,	0.65,	0.85,	0.85, "r"],
	"dr_sq":	["", "#Delta R Between Squarks",50,		0,		5,		0.7,	0.65,	0.85,	0.85, "p"],
	"dr_fj0":	["", "#Delta R^{0}_{min}",		50,		0,		5,		0.7,	0.65,	0.85,	0.85, "r"],
	"dm":		["", "\"Mass Drop\"",			40,		0,		4,		0.7,	0.65,	0.85,	0.85, "r"],
}
h2 = {#	name			x-axis label		bins	xmin	xmax	y-axis label		bins		ymin		ymax
	"y-phi":	["", "y",					50,		-5,		5,		"#phi",				50,			0,		6.4],
	"pt-dr":	["", "p_{T}^{0}",			50,		0,		1000,	"#Delta R^{0}",		50,			0,		3],
	"pt-m":		["", "p_{T}^{0}",			50,		0,		1000,	"m^{0}",			50,			0,		700]
}
th1 = {}
th2 = {}
for key, values in h1.iteritems():
	if (values[9] == "r"):
		for R in Rs:
			th1[key, "fat_gen", R] = TH1F("{0}_{1}_R{2}".format(key, "fat_gen", R), values[0], values[2], values[3], values[4])
	elif (values[9] == "p"):
		for point in points:
			th1[key, "fat_gen", point] = TH1F("{0}_{1}_msg{2}_msq{3}".format(key, "fat_gen", point[0], point[1]), values[0], values[2], values[3], values[4])
for key, values in h2.iteritems():
	th2[key, "fat_gen"] = TH2F("{0}_{1}".format(key, "fat_gen"), values[0], values[2], values[3], values[4], values[6], values[7], values[8])

output_string = ""
if do_true:
	output_string += "_true"
if do_prune:
	output_string += "_pruned"
if do_trim:
	output_string += "_trimmed_{0}_{1}".format(int(trim_R_test*10), int(trim_ptf_test*10))
print output_string
for R in Rs:
	print R
	tf_in = TFile("fatjets_sgtosq_msg{0}_msq{1}_R{2}{3}.root".format(msg_test, msq_test, int(R*10), output_string))	# Open the TFile
	tt_in = {}

	tt_in["fat_gen"] = TTree()		# Make an empty TTree to fill
	tf_in.GetObject("{0}/{1};1".format("analyzer", "fat_gen"), tt_in["fat_gen"])		# Get the TTree
	n_event = -1
	for event in tt_in["fat_gen"]:
		n_event += 1
		if (n_event % 500 == 0):
			print "Event {0}".format(n_event)
		y0 = event.y[0]
		y1 = event.y[1]
		phi0 = event.phi[0]
		phi1 = event.phi[1]
		dphi = abs(phi0 - phi1)
		if (dphi > 2*math.pi):
			dphi -= 2*math.pi
		dr_fj = ( (y0 - y1)**2.0 + (dphi)**2.0 )**(0.5)
		m0 = event.m[0]
		m1 = event.m[1]
		m2 = event.m[2]
		pt0 = event.pt[0]
		dr0_min = event.dr_fj[0]
		dm = (m1-m2)/(m0-m1)
		th1["njets", "fat_gen", R].Fill(len(event.m))
		th1["m0", "fat_gen", R].Fill(m0)
#		th1["m1", "fat_gen", R].Fill(event.m[1])
		th1["pt0", "fat_gen", R].Fill(event.pt[0])
		th1["dr_fj", "fat_gen", R].Fill(dr_fj)
		th1["dr_fj0", "fat_gen", R].Fill(dr0_min)
		th1["dm", "fat_gen", R].Fill(dm)
		if (dr0_min <= cut1_dr):
			th1["m0_cut1", "fat_gen", R].Fill(m0)
		if (dr_fj <= cut2_dr_high and dr_fj >= cut2_dr_low):
			th1["m0_cut2", "fat_gen", R].Fill(m0)
		if (n_event == 0 and R == min(Rs)):
			for fj in range(len(event.m)):
				th2["y-phi", "fat_gen"].Fill(event.y[fj], event.phi[fj])
		if (R_test == 1.2):
			th2["pt-dr", "fat_gen"].Fill(pt0, dr0_min)
			th2["pt-m", "fat_gen"].Fill(pt0, m0)

for point in points:
	msg = point[0]
	msq = point[1]
	print point
	tf_in = TFile("fatjets_sgtosq_msg{0}_msq{1}_R{2}{3}.root".format(msg, msq, int(R_test*10), output_string))	# Open the TFile
	tt_in = {}

	tt_in["fat_gen"] = TTree()		# Make an empty TTree to fill
	tf_in.GetObject("{0}/{1};1".format("analyzer", "fat_gen"), tt_in["fat_gen"])		# Get the TTree
	n_event = -1
	counter = 0
	for event in tt_in["fat_gen"]:
		n_event += 1
		if (event.m and len(event.m) > 1):
			if (event.m[0] < event.m[1]):
				counter += 1
		th1["dr_sq", "fat_gen", point].Fill(event.dr_sq[0])
	print counter

# PRINT H1s:
gPad.SetBottomMargin(0.15)
gPad.SetLeftMargin(0.2)
n_h1 = -1
for key, values in h1.iteritems():
	n_h1 += 1
	c1.Clear()
	tl = TLegend(values[5], values[6], values[7], values[8])	# W S E N
#	tl.SetHeader("The Legend Title")
	tl.SetFillColor(kWhite)
	maxy = 0
	if (values[9] == "r"):
		for R in Rs:
			y = th1[key, "fat_gen", R].GetMaximum()
			if (y > maxy):
				maxy = y
		i = -1
		for R in Rs:
			i += 1
			th1[key, "fat_gen", R].SetLineColor(colors[i])
			th1[key, "fat_gen", R].SetLineWidth(2)
			if (i == 0):
				th1[key, "fat_gen", R].SetMaximum(1.1*maxy)
				th1[key, "fat_gen", R].GetXaxis().SetTitle(values[1])
				th1[key, "fat_gen", R].GetXaxis().CenterTitle(1)
				th1[key, "fat_gen", R].GetXaxis().SetLabelSize(0.035)
				th1[key, "fat_gen", R].GetXaxis().SetTitleSize(0.05)
			#	th1[key, "fat_gen", R].GetXaxis().SetTitleOffset(1.8)
				th1[key, "fat_gen", R].GetYaxis().SetTitle("Events")
				th1[key, "fat_gen", R].GetYaxis().CenterTitle(1)
				th1[key, "fat_gen", R].GetYaxis().SetTitleSize(0.05)
				th1[key, "fat_gen", R].GetYaxis().SetTitleOffset(1.8)
				th1[key, "fat_gen", R].Draw()
			else:
				th1[key, "fat_gen", R].Draw("same")
			tl.AddEntry(th1[key, "fat_gen", R], "R = {0}".format(R), "l")
		tl.Draw()
		if save_plots:
			c1.SaveAs("plots/h1_{0}_{1}_msg{2}_msq{3}{4}.pdf".format(key, "fat_gen", msg_test, msq_test, output_string))
			c1.SaveAs("plots/h1_{0}_{1}_msg{2}_msq{3}{4}.svg".format(key, "fat_gen", msg_test, msq_test, output_string))
	elif (values[9] == "p"):
		for point in points:
			y = th1[key, "fat_gen", point].GetMaximum()
			if (y > maxy):
				maxy = y
		i = -1
		for point in points:
			i += 1
			th1[key, "fat_gen", point].SetLineColor(colors[i])
			th1[key, "fat_gen", point].SetLineWidth(2)
			if (i == 0):
				th1[key, "fat_gen", point].SetMaximum(1.1*maxy)
				th1[key, "fat_gen", point].GetXaxis().SetTitle(values[1])
				th1[key, "fat_gen", point].GetXaxis().CenterTitle(1)
				th1[key, "fat_gen", point].GetXaxis().SetLabelSize(0.035)
				th1[key, "fat_gen", point].GetXaxis().SetTitleSize(0.05)
			#	th1[key, "fat_gen", point].GetXaxis().SetTitleOffset(1.8)
				th1[key, "fat_gen", point].GetYaxis().SetTitle("Events")
				th1[key, "fat_gen", point].GetYaxis().CenterTitle(1)
				th1[key, "fat_gen", point].GetYaxis().SetTitleSize(0.05)
				th1[key, "fat_gen", point].GetYaxis().SetTitleOffset(1.8)
				th1[key, "fat_gen", point].Draw()
			else:
				th1[key, "fat_gen", point].Draw("same")
			tl.AddEntry(th1[key, "fat_gen", point], "{0}, {1}".format(point[0], point[1]), "l")
		tl.Draw()
		if save_plots:
			c1.SaveAs("plots/h1_{0}_{1}_R{2}{3}.pdf".format(key, "fat_gen", int(10*R_test), output_string))
			c1.SaveAs("plots/h1_{0}_{1}_R{2}{3}.svg".format(key, "fat_gen", int(10*R_test), output_string))
gPad.SetBottomMargin(0.15)
gPad.SetLeftMargin(0.1)
gPad.SetRightMargin(0.2)
for key, values in h2.iteritems():
	c1.Clear()
	tl.Clear()
	th2[key, "fat_gen"].GetXaxis().SetTitle(values[1])
	th2[key, "fat_gen"].GetXaxis().CenterTitle(1)
	th2[key, "fat_gen"].GetXaxis().SetLabelSize(0.035)
	th2[key, "fat_gen"].GetXaxis().SetTitleSize(0.05)
#	th2[key, "fat_gen"].GetXaxis().SetTitleOffset(1.8)
	th2[key, "fat_gen"].GetYaxis().SetTitle(values[5])
	th2[key, "fat_gen"].GetYaxis().CenterTitle(1)
	th2[key, "fat_gen"].GetYaxis().SetTitleSize(0.05)
	th2[key, "fat_gen"].GetYaxis().SetTitleOffset(0.7)
	th2[key, "fat_gen"].GetZaxis().SetTitle("Fatjets")
	th2[key, "fat_gen"].GetZaxis().CenterTitle(1)
	th2[key, "fat_gen"].GetZaxis().SetTitleSize(0.05)
	th2[key, "fat_gen"].GetZaxis().SetTitleOffset(1.0)
	th2[key, "fat_gen"].Draw("colz")
#	tl.AddEntry(th1[key, "fat_gen"], "{0}".format("fat_gen"), "l")
#	tl.Draw()
	if save_plots:
		c1.SaveAs("plots/h2_{0}_{1}_msg{2}_msq{3}_R{4}{5}.pdf".format(key, "fat_gen", msg_test, msq_test, int(R_test*10), output_string))
		c1.SaveAs("plots/h2_{0}_{1}_msg{2}_msq{3}_R{4}{5}.svg".format(key, "fat_gen", msg_test, msq_test, int(R_test*10), output_string))
