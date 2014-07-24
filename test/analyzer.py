from ROOT import *
import math
import sys

colors = [kBlack, kRed, kBlue, kSpring+2, kViolet]

studies = ["ca", "kt", "ktr"]

#gStyle.SetOptStat("eroum")
gStyle.SetOptStat(0)
c1 = TCanvas("c1", "c1", 1000, 500)
c1.SetCanvasSize(500, 500)

h1_fat = {#	name			x-axis label		bins	xmin	xmax	W		S		E		N
	"njets":	["", "Number Of Fatjets",		50,		0,		50,		0.7,	0.65,	0.85,	0.85],
	"pt0":		["", "p_{T}^{0} (GeV)",			50,		0,		1000,	0.7,	0.65,	0.85,	0.85],
	"pt1":		["", "p_{T}^{1} (GeV)",			50,		0,		1000,	0.7,	0.65,	0.85,	0.85],
	"m0":		["", "m^{0} (GeV)",				50,		0,		1000,	0.7,	0.65,	0.85,	0.85],
	"m1":		["", "m^{1} (GeV)",				50,		0,		1000,	0.7,	0.65,	0.85,	0.85],
}
h2_fat = {#	name			x-axis label		bins	xmin	xmax	y-axis label		bins		ymin		ymax
	"y-phi":	["", "y",						50,		-5,		5,		"#phi",				50,			-3.2,		3.2]
}
h1_studies = {#	name			x-axis label		bins	xmin	xmax	W		S		E		N
	"n":		["", "Number Of Subjets",			21,		-0.5,	20.5,	0.7,	0.65,	0.85,	0.85],
	"dpt0":		["", "#Delta p_{T}^{0} (GeV)",		50,		0,		1000,	0.7,	0.65,	0.85,	0.85],
	"dpt1":		["", "#Delta p_{T}^{1} (GeV)",		50,		0,		1000,	0.7,	0.65,	0.85,	0.85],
	"pt00":		["", "p_{T}^{00} (GeV)",			50,		0,		1000,	0.7,	0.65,	0.85,	0.85],
	"m00":		["", "m^{00} (GeV)",				50,		0,		500,	0.7,	0.65,	0.85,	0.85],
	"m01":		["", "m^{01} (GeV)",				50,		0,		500,	0.7,	0.65,	0.85,	0.85],
}

th1 = {}
th2 = {}
for key, values in h1_fat.iteritems():
	th1[key, "fat"] = TH1F("{0}_{1}".format(key, "fat"), values[0], values[2], values[3], values[4])
for key, values in h2_fat.iteritems():
	th2[key, "fat"] = TH2F("{0}_{1}".format(key, "fat"), values[0], values[2], values[3], values[4], values[6], values[7], values[8])
for key, values in h1_studies.iteritems():
	for study in studies:
		th1[key, study] = TH1F("{0}_{1}".format(key, study), values[0], values[2], values[3], values[4])

tf_in = TFile("ntuple_fatjet.root")	# Open the TFile
tt_in = {}

tt_in["fat"] = TTree()		# Make an empty TTree to fill
tf_in.GetObject("{0}/{1};1".format("analyzer", "ttree_fat"), tt_in["fat"])		# Get the TTree
print "============ fat ================"
n_event = -1
for event in tt_in["fat"]:
	n_event += 1
	th1["njets", "fat"].Fill(event.njets)
	th1["pt0", "fat"].Fill(event.pt[0])
	th1["pt1", "fat"].Fill(event.pt[1])
	th1["m0", "fat"].Fill(event.m[0])
	th1["m1", "fat"].Fill(event.m[1])
	if (n_event == 0):
		for fj in range(len(event.y)):
			th2["y-phi", "fat"].Fill(event.y[fj], event.phi[fj])

for study in studies:
	print "============ {0} =================".format(study)
	tt_in[study] = TTree()		# Make an empty TTree to fill
	tf_in.GetObject("{0}/{1};1".format("analyzer", "ttree_" + study), tt_in[study])		# Get the TTree
	n_event = -1
	for event in tt_in[study]:
		n_event += 1
#		print event.dpt0
#		print event.pt0_sub[0]
#		print "   " + str(event.pt0_sub[1])
#		print event.dpt1
#		print event.pt1_sub[0]
#		print "   " + str(event.pt1_sub[1])
#		print "------------------"
		th1["n", study].Fill(len(event.pt0_sub))
		dpt0 = event.pt0_sub[0]-event.pt0_sub[-1]
		dpt1 = event.pt1_sub[0]-event.pt1_sub[-1]
		if (dpt0 > 0):
			th1["dpt0", study].Fill(dpt0)
		if (dpt1 > 0):
			th1["dpt1", study].Fill(dpt1)
		th1["pt00", study].Fill(event.pt0_sub[0])
		th1["m00", study].Fill(event.m0_sub[0])
		if (len(event.m0_sub) > 1):
			th1["m01", study].Fill(event.m0_sub[1])


gPad.SetBottomMargin(0.15)
gPad.SetLeftMargin(0.2)
tl = TLegend(values[5], values[6], values[7], values[8])	# W S E N
#tl.SetHeader("The Legend Title")
tl.SetFillColor(kWhite)
for key, values in h1_fat.iteritems():
	c1.Clear()
	tl.Clear()
	th1[key, "fat"].SetLineColor(colors[0])
	th1[key, "fat"].SetLineWidth(2)
	th1[key, "fat"].GetXaxis().SetTitle(values[1])
	th1[key, "fat"].GetXaxis().CenterTitle(1)
	th1[key, "fat"].GetXaxis().SetLabelSize(0.035)
	th1[key, "fat"].GetXaxis().SetTitleSize(0.05)
#	th1[key, "fat"].GetXaxis().SetTitleOffset(1.8)
	th1[key, "fat"].GetYaxis().SetTitle("Events")
	th1[key, "fat"].GetYaxis().CenterTitle(1)
	th1[key, "fat"].GetYaxis().SetTitleSize(0.05)
	th1[key, "fat"].GetYaxis().SetTitleOffset(1.8)
	th1[key, "fat"].Draw()
#	tl.AddEntry(th1[key, "fat"], "{0}".format("fat"), "l")
#	tl.Draw()
	c1.SaveAs("plots_pruned/h1_{0}_{1}.pdf".format(key, "fat"))
#	c1.SaveAs("plots_pruned/h1_{0}_{1}.svg".format(key, "fat"))
gPad.SetBottomMargin(0.15)
gPad.SetLeftMargin(0.1)
gPad.SetRightMargin(0.2)
for key, values in h2_fat.iteritems():
	c1.Clear()
	tl.Clear()
	th2[key, "fat"].GetXaxis().SetTitle(values[1])
	th2[key, "fat"].GetXaxis().CenterTitle(1)
	th2[key, "fat"].GetXaxis().SetLabelSize(0.035)
	th2[key, "fat"].GetXaxis().SetTitleSize(0.05)
#	th2[key, "fat"].GetXaxis().SetTitleOffset(1.8)
	th2[key, "fat"].GetYaxis().SetTitle(values[5])
	th2[key, "fat"].GetYaxis().CenterTitle(1)
	th2[key, "fat"].GetYaxis().SetTitleSize(0.05)
	th2[key, "fat"].GetYaxis().SetTitleOffset(0.7)
	th2[key, "fat"].GetZaxis().SetTitle("Fatjets")
	th2[key, "fat"].GetZaxis().CenterTitle(1)
	th2[key, "fat"].GetZaxis().SetTitleSize(0.05)
	th2[key, "fat"].GetZaxis().SetTitleOffset(1.0)
	th2[key, "fat"].Draw("colz")
#	tl.AddEntry(th1[key, "fat"], "{0}".format("fat"), "l")
#	tl.Draw()
	c1.SaveAs("plots_pruned/h2_{0}_{1}.pdf".format(key, "fat"))
#	c1.SaveAs("plots_pruned/h2_{0}_{1}.svg".format(key, "fat"))
gPad.SetBottomMargin(0.15)
gPad.SetLeftMargin(0.2)
gPad.SetRightMargin(0.1)
for key, values in h1_studies.iteritems():
	c1.Clear()
	tl.Clear()
	maxy = 0
	for study in studies:
		y = th1[key, study].GetMaximum()
		if (y > maxy):
			maxy = y
	i = -1
	for study in studies:
		i += 1
		th1[key, study].SetLineColor(colors[i])
		th1[key, study].SetLineWidth(2)
		if i == 0:
			th1[key, study].SetMaximum(1.1*maxy)
			th1[key, study].GetXaxis().SetTitle(values[1])
			th1[key, study].GetXaxis().CenterTitle(1)
			th1[key, study].GetXaxis().SetLabelSize(0.035)
			th1[key, study].GetXaxis().SetTitleSize(0.05)
#			th1[key, study].GetXaxis().SetTitleOffset(1.8)
			th1[key, study].GetYaxis().SetTitle("Events")
			th1[key, study].GetYaxis().CenterTitle(1)
			th1[key, study].GetYaxis().SetTitleSize(0.05)
			th1[key, study].GetYaxis().SetTitleOffset(1.8)
			th1[key, study].Draw()
		else:
			th1[key, study].Draw("same")
		tl.AddEntry(th1[key, study], "{0}".format(study), "l")
		tl.Draw()
		c1.SaveAs("plots_pruned/h1_{0}_{1}.pdf".format(key, "studies"))
#		c1.SaveAs("plots_pruned/h1_{0}_{1}.svg".format(key, "studies"))
