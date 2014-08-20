from ROOT import *
import math
import sys

msg_test = 850
msq_test = 250
#msg_test = 850
#msq_test = 250
R_test = 1.2
trim_R_test = 0.7
trim_ptf_test = 0.1
#trim_test = (trim_R_test, trim_ptf_test)

do_gen = True
do_pf = True
do_true = False
do_prune = False
do_trim = False
save_plots = True

Rs = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
Rs = [0.8, 1.2, 1.6, 2.0,]#  2.4, 2.8]
#Rs = [1.0, 1.1, 1.2, 1.3]
#Rs = [1.2, 1.4, 1.6]
#Rs = [1.2, 1.4]
Rs = [0.8, 1.0, 1.2, 1.6]
Rs = [1.2]

points = [
	(850, 250),
#	(1150, 150),
]

trims = [
	(0.5, 0.1),
	(0.7, 0.1),
#	(0.8, 0.1),
	(1.0, 0.1),
]
#trims = [
#	(0.8, 0.1),
#	(0.8, 0.2),
##	(0.8, 0.3),
#]


fj_types = []
if (do_gen): fj_types.append("gen")
if (do_pf): fj_types.append("pf")
if len(fj_types) == 0: print "Error: What type of fatjet do you want to look at?"
#print fj_types

colors = [kBlack, kRed, kBlue, kSpring+2, kViolet, kGray+2, kOrange-6, kCyan+2]
if ( len(Rs) > len(colors) or len(points) > len(colors) ): print "Error: Add more colors!"

cut1_dr = 0.1
cut2_dr_low = 2.75
cut2_dr_high = 3.75

gROOT.SetBatch()
#gStyle.SetOptStat("eroum")
gStyle.SetOptStat(0)
c1 = TCanvas("c1", "c1", 1000, 500)
c1.SetCanvasSize(500, 500)

h1 = {#	name			x-axis label		bins	xmin	xmax	W		S		E		N
	"njets":	["", "Number Of Fatjets",		60,		0,		60,		0.7,	0.65,	0.85,	0.85, "r", "genpf"],
	"m0":		["", "m^{0} (GeV)",				40,		0,		1000,	0.7,	0.65,	0.85,	0.85, "rt", "genpf"],
#	"m0_cut1":	["", "m^{0} (GeV)",				50,		0,		1000,	0.7,	0.65,	0.85,	0.85, "r", "genpf"],
	"m0_cut2":	["", "m^{0} (GeV)",				50,		0,		1000,	0.7,	0.65,	0.85,	0.85, "r", "genpf"],
	"m1":		["", "m^{1} (GeV)",				50,		0,		1000,	0.7,	0.65,	0.85,	0.85, "r", "genpf"],
	"pt0":		["", "p_{T} (GeV)",				50,		0,		1000,	0.7,	0.65,	0.85,	0.85, "r", "genpf"],
	"dr_fj":	["", "#Delta R^{01}",			50,		0,		5,		0.7,	0.65,	0.85,	0.85, "r", "genpf"],
	"dr_sq":	["", "#Delta R Between Squarks",50,		0,		5,		0.7,	0.65,	0.85,	0.85, "p", "gen"],
	"dr_fj0":	["", "#Delta R^{0}_{min}",		50,		0,		5,		0.7,	0.65,	0.85,	0.85, "r", "gen"],
	"dm":		["", "\"Mass Drop\"",			40,		0,		4,		0.7,	0.65,	0.85,	0.85, "r", "genpf"],
}
h2 = {#	name			x-axis label		bins	xmin	xmax	y-axis label		bins		ymin		ymax
	"y-phi":	["", "y",					50,		-5,		5,		"#phi",				50,			0,		6.4],
	"pt-dr":	["", "p_{T}^{0}",			50,		0,		1000,	"#Delta R^{0}",		50,			0,		3],
	"pt-m":		["", "p_{T}^{0}",			50,		0,		1000,	"m^{0}",			50,			0,		700]
}
th1 = {}
th2 = {}
for fj_type in fj_types:
	for key, values in h1.iteritems():
		if fj_type in values[10]:
			if "r" in values[9]:
				for R in Rs:
					th1[key, "fat_{0}".format(fj_type), R] = TH1F("{0}_fat{1}_R{2}".format(key, fj_type, R), values[0], values[2], values[3], values[4])
			if "p" in values[9]:
				for point in points:
					th1[key, "fat_{0}".format(fj_type), point] = TH1F("{0}_fat{1}_msg{2}_msq{3}".format(key, fj_type, point[0], point[1]), values[0], values[2], values[3], values[4])
			if "t" in values[9]:
				for trim in trims:
					th1[key, "fat_{0}".format(fj_type), trim] = TH1F("{0}_fat{1}_tR{2}_ptf{3}".format(key, fj_type, trim[0], trim[1]), values[0], values[2], values[3], values[4])
	for key, values in h2.iteritems():
			th2[key, "fat_{0}".format(fj_type)] = TH2F("{0}_fat{1}".format(key, fj_type), values[0], values[2], values[3], values[4], values[6], values[7], values[8])

output_string = ""
if do_true:
	output_string += "_true"
if do_prune:
	output_string += "_pruned"
if do_trim:
	output_string += "_trimmed_{0}_{1}".format(int(trim_R_test*10), int(trim_ptf_test*10))
print "The output string is \"{0}\".".format(output_string)

for fj_type in fj_types:
	print "Making {0} fatjet R plots:".format(fj_type)
	for R in Rs:
		print ">> Looking at distance parameter {0}.".format(R)
		tf_in = TFile("fatjets_sgtosq_msg{0}_msq{1}_R{2}{3}.root".format(msg_test, msq_test, int(R*10), output_string))	# Open the TFile
		tt_in = {}
		tt_in["fat_{0}".format(fj_type)] = TTree()		# Make an empty TTree to fill
		tf_in.GetObject("{0}/{1};1".format("analyzer", "fat_{0}".format(fj_type)), tt_in["fat_{0}".format(fj_type)])		# Get the TTree
		n_event = -1
		for event in tt_in["fat_{0}".format(fj_type)]:
			n_event += 1
			if (n_event % 500 == 0):
				print ">> Processing Event {0}".format(n_event)
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
			if fj_type == "gen": dr0_min = event.dr_fj[0]
			dm = (m1-m2)/(m0-m1)
			th1["njets", "fat_{0}".format(fj_type), R].Fill(len(event.m))
			th1["m0", "fat_{0}".format(fj_type), R].Fill(m0)
	#		th1["m1", "fat_{0}".format(fj_type), R].Fill(event.m[1])
			th1["pt0", "fat_{0}".format(fj_type), R].Fill(event.pt[0])
			th1["dr_fj", "fat_{0}".format(fj_type), R].Fill(dr_fj)
			if fj_type == "gen": th1["dr_fj0", "fat_{0}".format(fj_type), R].Fill(dr0_min)
			th1["dm", "fat_{0}".format(fj_type), R].Fill(dm)
#			if (dr0_min <= cut1_dr):
#				th1["m0_cut1", "fat_{0}".format(fj_type), R].Fill(m0)
			if (dr_fj <= cut2_dr_high and dr_fj >= cut2_dr_low):
				th1["m0_cut2", "fat_{0}".format(fj_type), R].Fill(m0)
			if (n_event == 0 and R == min(Rs)):
				for fj in range(len(event.m)):
					th2["y-phi", "fat_{0}".format(fj_type)].Fill(event.y[fj], event.phi[fj])
			if (R == R_test):
				th2["pt-dr", "fat_{0}".format(fj_type)].Fill(pt0, dr0_min)
				th2["pt-m", "fat_{0}".format(fj_type)].Fill(pt0, m0)

for fj_type in fj_types:
	print "Making {0} fatjet mass-point plots:".format(fj_type)
	for point in points:
		msg = point[0]
		msq = point[1]
		print ">> Looking at guino mass of {0} and squark mass of {1}.".format(msg, msq)
		tf_in = TFile("fatjets_sgtosq_msg{0}_msq{1}_R{2}{3}.root".format(msg, msq, int(R_test*10), output_string))	# Open the TFile
		tt_in = {}
		tt_in["fat_{0}".format(fj_type)] = TTree()		# Make an empty TTree to fill
		tf_in.GetObject("{0}/{1};1".format("analyzer", "fat_{0}".format(fj_type)), tt_in["fat_{0}".format(fj_type)])		# Get the TTree
		n_event = -1
		for event in tt_in["fat_{0}".format(fj_type)]:
			n_event += 1
			if (n_event % 500 == 0):
				print ">> Processing Event {0}".format(n_event)
			if fj_type == "gen": th1["dr_sq", "fat_{0}".format(fj_type), point].Fill(event.dr_sq[0])

if (do_trim and not do_prune and not do_true):		# Otherwise the output_string gets jacked up.
	for fj_type in fj_types:
		print "Making {0} fatjet trim plots:".format(fj_type)
		for trim in trims:
			tR = trim[0]
			ptf = trim[1]
			output_string_temp = "_trimmed_{0}_{1}".format(int(tR*10), int(ptf*10))
			print ">> Looking at trim raidius of {0} and pt fraction of {1}.".format(tR, ptf)
			tf_in = TFile("fatjets_sgtosq_msg{0}_msq{1}_R{2}{3}.root".format(msg_test, msq_test, int(R_test*10), output_string_temp))	# Open the TFile
			tt_in = {}
			tt_in["fat_{0}".format(fj_type)] = TTree()		# Make an empty TTree to fill
			tf_in.GetObject("{0}/{1};1".format("analyzer", "fat_{0}".format(fj_type)), tt_in["fat_{0}".format(fj_type)])		# Get the TTree
			n_event = -1
			for event in tt_in["fat_{0}".format(fj_type)]:
				n_event += 1
				if (n_event % 500 == 0):
					print ">> Processing Event {0}".format(n_event)
				th1["m0", "fat_{0}".format(fj_type), trim].Fill(event.m[0])

# PRINT H1s:
gPad.SetBottomMargin(0.15)
gPad.SetLeftMargin(0.2)
for fj_type in fj_types:
	n_h1 = -1
	for key, values in h1.iteritems():
		if fj_type in values[10]:
			n_h1 += 1
			c1.Clear()
			maxy = 0
			if "r" in values[9]:
				tl = TLegend(values[5], values[6], values[7], values[8])	# W S E N
			#	tl.SetHeader("The Legend Title")
				tl.SetFillColor(kWhite)
				for R in Rs:
					y = th1[key, "fat_{0}".format(fj_type), R].GetMaximum()
					if (y > maxy):
						maxy = y
				i = -1
				for R in Rs:
					i += 1
					th1[key, "fat_{0}".format(fj_type), R].SetLineColor(colors[i])
					th1[key, "fat_{0}".format(fj_type), R].SetLineWidth(2)
					if (i == 0):
						th1[key, "fat_{0}".format(fj_type), R].SetMaximum(1.1*maxy)
						th1[key, "fat_{0}".format(fj_type), R].GetXaxis().SetTitle(values[1])
						th1[key, "fat_{0}".format(fj_type), R].GetXaxis().CenterTitle(1)
						th1[key, "fat_{0}".format(fj_type), R].GetXaxis().SetLabelSize(0.035)
						th1[key, "fat_{0}".format(fj_type), R].GetXaxis().SetTitleSize(0.05)
					#	th1[key, "fat_{0}".format(fj_type), R].GetXaxis().SetTitleOffset(1.8)
						th1[key, "fat_{0}".format(fj_type), R].GetYaxis().SetTitle("Events")
						th1[key, "fat_{0}".format(fj_type), R].GetYaxis().CenterTitle(1)
						th1[key, "fat_{0}".format(fj_type), R].GetYaxis().SetTitleSize(0.05)
						th1[key, "fat_{0}".format(fj_type), R].GetYaxis().SetTitleOffset(1.8)
						th1[key, "fat_{0}".format(fj_type), R].Draw()
					else:
						th1[key, "fat_{0}".format(fj_type), R].Draw("same")
					tl.AddEntry(th1[key, "fat_{0}".format(fj_type), R], "R = {0}".format(R), "l")
				tl.Draw()
				if save_plots:
					c1.SaveAs("plots/h1_{0}_{1}_msg{2}_msq{3}{4}.pdf".format(key, "fat_{0}".format(fj_type), msg_test, msq_test, output_string))
					c1.SaveAs("plots/h1_{0}_{1}_msg{2}_msq{3}{4}.svg".format(key, "fat_{0}".format(fj_type), msg_test, msq_test, output_string))
			if "p" in values[9]:
				tl = TLegend(values[5], values[6], values[7], values[8])	# W S E N
			#	tl.SetHeader("The Legend Title")
				tl.SetFillColor(kWhite)
				for point in points:
					y = th1[key, "fat_{0}".format(fj_type), point].GetMaximum()
					if (y > maxy):
						maxy = y
				i = -1
				for point in points:
					i += 1
					th1[key, "fat_{0}".format(fj_type), point].SetLineColor(colors[i])
					th1[key, "fat_{0}".format(fj_type), point].SetLineWidth(2)
					if (i == 0):
						th1[key, "fat_{0}".format(fj_type), point].SetMaximum(1.1*maxy)
						th1[key, "fat_{0}".format(fj_type), point].GetXaxis().SetTitle(values[1])
						th1[key, "fat_{0}".format(fj_type), point].GetXaxis().CenterTitle(1)
						th1[key, "fat_{0}".format(fj_type), point].GetXaxis().SetLabelSize(0.035)
						th1[key, "fat_{0}".format(fj_type), point].GetXaxis().SetTitleSize(0.05)
					#	th1[key, "fat_{0}".format(fj_type), point].GetXaxis().SetTitleOffset(1.8)
						th1[key, "fat_{0}".format(fj_type), point].GetYaxis().SetTitle("Events")
						th1[key, "fat_{0}".format(fj_type), point].GetYaxis().CenterTitle(1)
						th1[key, "fat_{0}".format(fj_type), point].GetYaxis().SetTitleSize(0.05)
						th1[key, "fat_{0}".format(fj_type), point].GetYaxis().SetTitleOffset(1.8)
						th1[key, "fat_{0}".format(fj_type), point].Draw()
					else:
						th1[key, "fat_{0}".format(fj_type), point].Draw("same")
					tl.AddEntry(th1[key, "fat_{0}".format(fj_type), point], "{0}, {1}".format(point[0], point[1]), "l")
				tl.Draw()
				if save_plots:
					c1.SaveAs("plots/h1_{0}_{1}_R{2}{3}.pdf".format(key, "fat_{0}".format(fj_type), int(10*R_test), output_string))
					c1.SaveAs("plots/h1_{0}_{1}_R{2}{3}.svg".format(key, "fat_{0}".format(fj_type), int(10*R_test), output_string))
			if ("t" in values[9]) and do_trim:
				tl = TLegend(values[5], values[6], values[7], values[8])	# W S E N
			#	tl.SetHeader("The Legend Title")
				tl.SetFillColor(kWhite)
				for trim in trims:
					y = th1[key, "fat_{0}".format(fj_type), trim].GetMaximum()
					if (y > maxy):
						maxy = y
				i = -1
				for trim in trims:
					i += 1
					th1[key, "fat_{0}".format(fj_type), trim].SetLineColor(colors[i])
					th1[key, "fat_{0}".format(fj_type), trim].SetLineWidth(2)
					if (i == 0):
						th1[key, "fat_{0}".format(fj_type), trim].SetMaximum(1.1*maxy)
						th1[key, "fat_{0}".format(fj_type), trim].GetXaxis().SetTitle(values[1])
						th1[key, "fat_{0}".format(fj_type), trim].GetXaxis().CenterTitle(1)
						th1[key, "fat_{0}".format(fj_type), trim].GetXaxis().SetLabelSize(0.035)
						th1[key, "fat_{0}".format(fj_type), trim].GetXaxis().SetTitleSize(0.05)
					#	th1[key, "fat_{0}".format(fj_type), trim].GetXaxis().SetTitleOffset(1.8)
						th1[key, "fat_{0}".format(fj_type), trim].GetYaxis().SetTitle("Events")
						th1[key, "fat_{0}".format(fj_type), trim].GetYaxis().CenterTitle(1)
						th1[key, "fat_{0}".format(fj_type), trim].GetYaxis().SetTitleSize(0.05)
						th1[key, "fat_{0}".format(fj_type), trim].GetYaxis().SetTitleOffset(1.8)
						th1[key, "fat_{0}".format(fj_type), trim].Draw()
					else:
						th1[key, "fat_{0}".format(fj_type), trim].Draw("same")
					tl.AddEntry(th1[key, "fat_{0}".format(fj_type), trim], "{0}, {1}".format(trim[0], trim[1]), "l")
				tl.Draw()
				if save_plots:
					c1.SaveAs("plots/h1_{0}_{1}_msg{2}_msq{3}_R{4}_trimmed.pdf".format(key, "fat_{0}".format(fj_type), msg_test, msq_test, int(10*R_test)))
					c1.SaveAs("plots/h1_{0}_{1}_msg{2}_msq{3}_R{4}_trimmed.svg".format(key, "fat_{0}".format(fj_type), msg_test, msq_test, int(10*R_test)))
gPad.SetBottomMargin(0.15)
gPad.SetLeftMargin(0.1)
gPad.SetRightMargin(0.2)
for fj_type in fj_types:
	for key, values in h2.iteritems():
		c1.Clear()
		tl.Clear()
		th2[key, "fat_{0}".format(fj_type)].GetXaxis().SetTitle(values[1])
		th2[key, "fat_{0}".format(fj_type)].GetXaxis().CenterTitle(1)
		th2[key, "fat_{0}".format(fj_type)].GetXaxis().SetLabelSize(0.035)
		th2[key, "fat_{0}".format(fj_type)].GetXaxis().SetTitleSize(0.05)
	#	th2[key, "fat_{0}".format(fj_type)].GetXaxis().SetTitleOffset(1.8)
		th2[key, "fat_{0}".format(fj_type)].GetYaxis().SetTitle(values[5])
		th2[key, "fat_{0}".format(fj_type)].GetYaxis().CenterTitle(1)
		th2[key, "fat_{0}".format(fj_type)].GetYaxis().SetTitleSize(0.05)
		th2[key, "fat_{0}".format(fj_type)].GetYaxis().SetTitleOffset(0.7)
		th2[key, "fat_{0}".format(fj_type)].GetZaxis().SetTitle("Fatjets")
		th2[key, "fat_{0}".format(fj_type)].GetZaxis().CenterTitle(1)
		th2[key, "fat_{0}".format(fj_type)].GetZaxis().SetTitleSize(0.05)
		th2[key, "fat_{0}".format(fj_type)].GetZaxis().SetTitleOffset(1.0)
		th2[key, "fat_{0}".format(fj_type)].Draw("colz")
	#	tl.AddEntry(th1[key, "fat_{0}".format(fj_type)], "{0}".format("fat_{0}".format(fj_type)), "l")
	#	tl.Draw()
		if save_plots:
			c1.SaveAs("plots/h2_{0}_{1}_msg{2}_msq{3}_R{4}{5}.pdf".format(key, "fat_{0}".format(fj_type), msg_test, msq_test, int(R_test*10), output_string))
			c1.SaveAs("plots/h2_{0}_{1}_msg{2}_msq{3}_R{4}{5}.svg".format(key, "fat_{0}".format(fj_type), msg_test, msq_test, int(R_test*10), output_string))
if do_trim:
	c1.Clear()
	th1["m0", "fat_pf", trims[1]].Draw()
	th1["m0", "fat_pf", 1.2].Draw("same")
	if save_plots:
		c1.SaveAs("plots/h1_test.pdf")
		c1.SaveAs("plots/h2_test.svg")
