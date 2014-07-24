from ROOT import *
import math
import sys

msg = 850
msq = 250

Rs = [0.8, 1.0, 1.2, 1.4]

colors = [kBlack, kRed, kBlue, kSpring+2, kViolet]

gROOT.SetBatch()
#gStyle.SetOptStat("eroum")
gStyle.SetOptStat(0)
c1 = TCanvas("c1", "c1", 1000, 500)
c1.SetCanvasSize(500, 500)

h1 = {#	name			x-axis label		bins	xmin	xmax	W		S		E		N
	"njets":	["", "Number Of Fatjets",		50,		10,		60,		0.7,	0.65,	0.85,	0.85],
	"m0":		["", "m^{0} (GeV)",				50,		0,		1000,	0.7,	0.65,	0.85,	0.85],
	"m1":		["", "m^{1} (GeV)",				50,		0,		1000,	0.7,	0.65,	0.85,	0.85],
}
h2 = {#	name			x-axis label		bins	xmin	xmax	y-axis label		bins		ymin		ymax
	"y-phi":	["", "y",						50,		-5,		5,		"#phi",				50,			0,		6.4]
}
th1 = {}
th2 = {}
for R in Rs:
	for key, values in h1.iteritems():
		th1[key, "fat_gen", R] = TH1F("{0}_{1}_R{2}".format(key, "fat_gen", R), values[0], values[2], values[3], values[4])
for key, values in h2.iteritems():
	th2[key, "fat_gen"] = TH2F("{0}_{1}".format(key, "fat_gen"), values[0], values[2], values[3], values[4], values[6], values[7], values[8])

for R in Rs:
	tf_in = TFile("fatjets_sgtosq_msg{0}_msq{1}_R{2}.root".format(msg, msq, int(R*10)))	# Open the TFile
	tt_in = {}

	tt_in["fat_gen"] = TTree()		# Make an empty TTree to fill
	tf_in.GetObject("{0}/{1};1".format("analyzer", "fat_gen"), tt_in["fat_gen"])		# Get the TTree
	n_event = -1
	for event in tt_in["fat_gen"]:
#		print n_event
		n_event += 1
		th1["njets", "fat_gen", R].Fill(len(event.m))
		th1["m0", "fat_gen", R].Fill(event.m[0])
		th1["m1", "fat_gen", R].Fill(event.m[1])
		if (n_event == 0 and R == min(Rs)):
			for fj in range(len(event.m)):
				th2["y-phi", "fat_gen"].Fill(event.y[fj], event.phi[fj])

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
	c1.SaveAs("plots/h1_{0}_{1}.pdf".format(key, "fat_gen"))
	c1.SaveAs("plots/h1_{0}_{1}.svg".format(key, "fat_gen"))
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
	c1.SaveAs("plots/h2_{0}_{1}.pdf".format(key, "fat_gen"))
	c1.SaveAs("plots/h2_{0}_{1}.svg".format(key, "fat_gen"))
